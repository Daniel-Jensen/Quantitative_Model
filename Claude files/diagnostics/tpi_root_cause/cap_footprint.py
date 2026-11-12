# THE 200 bp CAP'S FOOTPRINT ON A SOLVED NO-TPI RULE SET: EACH NODE'S SPREAD AND THE ONE-OFF PURCHASE THAT HOLDS THE CAP
# THERE (NODE ALONE, CONTINUATION FIXED), AGAINST THE BOX'S EUROSYSTEM-BOOK EDGE. Usage: cap_footprint.py RULES.pkl [cap]
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
from global_projection.solver_recursive.state_grid import IS, IM, STATE_NAMES, default_prob
from global_projection.solver_recursive import point_map as PM

rules = pickle.load(open(sys.argv[1], "rb"))
cap = float(sys.argv[2]) if len(sys.argv) > 2 else 200.0
cal, ss, sproc = economy()
g, B, ngh = rules.grid, cal["B_gov_D_ss"], rules.n_gh or 5
M_edge = np.max(np.abs(g.points[:, IM]))
u = g.to_unit(g.points)
rows = []
for i in range(g.n):
    S = g.points[i]
    x = np.array([rules.vals[k][0][i] for k in SOLVE])
    c0 = dict(cal, tpi_on=False)
    o0 = PM.point_residuals(S, 0, x, rules, c0, ss, sproc, n_gh=ngh)[1]
    sp = _sov_bp(o0, cal)
    m = Mp = np.nan
    if sp > cap:
        c1 = dict(cal, tpi_on=True, tpi_cap_bp=cap)
        f = lambda z: PM.point_residuals(S, 0, z, rules, c1, ss, sproc, n_gh=ngh)[0]
        sol = root(f, x, method="hybr", tol=1e-13)
        r, o1 = PM.point_residuals(S, 0, sol.x, rules, c1, ss, sproc, n_gh=ngh)
        m, Mp = 100 * o1["m_cb"] / B, o1["M_cb_new"] / M_edge
        if np.max(np.abs(r)) > 1e-8:
            m = -m  # flagged: the node solve did not converge
    edge = ",".join(f"{STATE_NAMES[j]}{'+' if u[i, j] > 0 else '-'}" for j in range(12) if abs(u[i, j]) > 0.99 and j != IS)
    rows.append((sp, i, 100 * default_prob(S[IS]), m, Mp, edge or "centre"))
print(f"cap {cap:.0f} bp; box top p^d {100 * default_prob(g.hi[IS]):.2f}%; M edge {100 * M_edge / B:.1f}% of B")
print(f"nodes above the cap: {sum(r[0] > cap for r in rows)}; highest spread {max(r[0] for r in rows):.1f} bp")
for sp, i, pdq, m, Mp, edge in sorted(rows, reverse=True)[:28]:
    tag = "" if np.isnan(m) else (f"buys {m:5.1f}% of B, M'/edge {Mp:5.2f}" if m >= 0 else f"NODE SOLVE FAILED ({-m:.1f}%)")
    print(f"  node {i:2d}  p^d {pdq:5.2f}%  spread {sp:6.1f}  {edge:10s} {tag}")
