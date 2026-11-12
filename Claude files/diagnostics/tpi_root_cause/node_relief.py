# WHERE A CONSTRAINT-RELIEF TPI WOULD BUY (2026-11-06, READ-ONLY): AT EVERY NODE OF THE APPENDIX-B 12% BOX WHERE THE D
# CONSTRAINT BINDS WITHOUT THE TPI, THE ONE-OFF PURCHASE (RULES HELD FIXED) THAT SWITCHES IT OFF, AGAINST THE D BANK'S BOOK.
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

if __name__ == "__main__":
    from crisis_shock import setup12
    from one_off_purchase import read_forced
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    from global_projection.solver_recursive.state_grid import STATE_NAMES, IS, IBDD, IM, default_prob
    cal, ss, sproc, rules = setup12()
    g = rules.grid
    B = cal["B_gov_D_ss"]
    rows = []
    for i in range(g.n):
        S = g.points[i]
        x = np.array([rules.vals[k][0][i] for k in SOLVE])
        o = point_residuals(S, 0, x, rules, cal, ss, sproc, n_gh=rules.n_gh or 5)[1]
        u = g.points_unit[i]
        lab = ",".join(f"{STATE_NAMES[j]}{u[j]:+.0f}" for j in range(g.d) if j != IS and abs(u[j]) > 1e-9) or "axis"
        if o["mu_D"] <= 1e-6:
            rows.append((i, lab, float(default_prob(S[IS])), 0.0, 0.0, o["b_D_D_new"] / B, S[IM] / B))
            continue
        f = lambda m: read_forced(rules, cal, ss, sproc, S, m)["mu_D"] - 1e-6
        hi = next((m for m in (0.1, 0.2, 0.3, 0.45, 0.6, 0.8) if f(m) < 0), None)
        m = brentq(f, 0.0, hi, xtol=1e-4) if hi else np.inf
        rows.append((i, lab, float(default_prob(S[IS])), float(o["mu_D"]), m, o["b_D_D_new"] / B, S[IM] / B))
    bind = [r for r in rows if r[3] > 0]
    print(f"D constraint binds at {len(bind)}/{g.n} nodes (no TPI, no-default regime)")
    print(f"{'node':>4} {'where':>18} {'p^d':>6} {'mu_D':>6} {'buy to relieve':>15} {'bank book':>10} {'M':>7}")
    for i, lab, pd, mu, m, b, M in sorted(bind, key=lambda r: -r[4]):
        print(f"{i:4d} {lab:>18} {100*pd:5.2f}% {mu:6.3f} {100*m:14.1f}% {100*b:9.1f}% {100*M:6.1f}%")
    ms = np.array([r[4] for r in bind])
    print(f"\npurchase needed: median {100*np.median(ms):.1f}%, max {100*np.max(ms):.1f}% of the debt; "
          f"{int(np.sum(ms > 0.28))} nodes need more than the current book band (28%)")
