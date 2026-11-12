# WHICH RISK-BOX TOP KEEPS THE TPI'S PURCHASES SMALL: THE 23 TOP-SLICE NODE CONFIGURATIONS OF THE 6.2% GRID, READ AT A
# LOWER p^d ON THE SOLVED NO-TPI RULES (EXACT NODE SOLVES), THEN THE ONE-OFF PURCHASE THAT HOLDS A 200 bp CAP THERE.
import pickle
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import root

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))
from rung_anatomy import economy
from global_projection.solver_recursive.decision_rules import SOLVE
from global_projection.solver_recursive.recursive_experiment import _sov_bp
from global_projection.solver_recursive.state_grid import IS
from global_projection.solver_recursive import point_map as PM

ROOT = HERE.parents[2]
cal, ss, sproc = economy()
rules = pickle.load(open(ROOT / "results" / "GLOBAL_cache" / "risk_c9717a6268e5.pkl", "rb"))
g, B, ngh = rules.grid, cal["B_gov_D_ss"], rules.n_gh or 5
top = [i for i in range(g.n) if g.points[i][IS] > g.hi[IS] - 1e-9]
cap = float(sys.argv[1]) if len(sys.argv) > 1 else 200.0


def node(S, tpi):
    c = dict(cal, tpi_on=tpi, tpi_cap_bp=cap)
    x = np.array([float(rules.eval(k, 0, S[None, :])[0]) for k in SOLVE])
    if not tpi:
        x[-1] = 0.0
    f = lambda z: PM.point_residuals(S, 0, z, rules, c, ss, sproc, n_gh=ngh)[0]
    sol = root(f, x, method="hybr", tol=1e-13)
    r, o = PM.point_residuals(S, 0, sol.x, rules, c, ss, sproc, n_gh=ngh)
    return o, np.max(np.abs(r))


print(f"cap {cap:.0f} bp; {len(top)} top-slice configurations; M box edge = 27.9% of B")
print("  p^d    no-TPI spread: min / centre / max      buy%B at the cap: max (nodes buying)   worst resid")
for pd in (0.015, 0.0175, 0.02, 0.0225, 0.025, 0.03):
    sp, buy, res = [], [], []
    for i in top:
        S = g.points[i].copy(); S[IS] = np.log(pd / (1 - pd))
        o0, r0 = node(S, False)
        o1, r1 = node(S, True)
        sp.append(_sov_bp(o0, cal)); buy.append(100 * o1["m_cb"] / B); res.append(max(r0, r1))
    c = sp[[k for k, i in enumerate(top) if np.allclose(g.to_unit(g.points[i:i + 1])[0][np.arange(12) != IS], 0)][0]]
    print(f"  {100 * pd:4.2f}%  {min(sp):7.1f} / {c:7.1f} / {max(sp):7.1f}      {max(buy):6.1f} ({sum(b > 1e-3 for b in buy):2d})"
          f"                      {max(res):.0e}")
