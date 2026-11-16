# THE EUROSYSTEM-BOOK EDGES OF A RULE SET: THE D BANK'S BOOK AT THE M+ / M- NODES (A SHORT BOOK IS A CORNER IN THE NO-TPI
# MAP), AND ONE NODE'S CAP SOLVE AT TWO SMOOTHINGS (SOLVER OR RULES?). Usage: RULES.pkl [node]
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
from global_projection.solver_recursive.state_grid import IS, IM, IBDD, default_prob
from global_projection.solver_recursive import point_map as PM

rules = pickle.load(open(sys.argv[1], "rb"))
node = int(sys.argv[2]) if len(sys.argv) > 2 else 8
cal, ss, sproc = economy()
g, B, ngh = rules.grid, cal["B_gov_D_ss"], rules.n_gh or 5
u = g.to_unit(g.points)
for i in np.where(np.abs(u[:, IM]) > 0.99)[0]:
    x = np.array([rules.vals[k][0][i] for k in SOLVE])
    o = PM.point_residuals(g.points[i], 0, x, rules, dict(cal, tpi_on=False), ss, sproc, n_gh=ngh)[1]
    print(f"node {i:2d} M{'+' if u[i, IM] > 0 else '-'} p^d {100 * default_prob(g.points[i, IS]):5.2f}%: M_lag {100 * g.points[i, IM] / B:6.1f}% of B, "
          f"b_DD_lag {100 * g.points[i, IBDD] / B:6.1f}%, b_DD' {100 * o['b_D_D_new'] / B:6.1f}%, spread {_sov_bp(o, cal):6.1f}")
S = g.points[node]
x0 = np.array([rules.vals[k][0][node] for k in SOLVE])
for eps in (1e-2, 1e-3, 1e-4):
    c = dict(cal, tpi_on=True, tpi_cap_bp=200.0, tpi_eps=eps)
    f = lambda z: PM.point_residuals(S, 0, z, rules, c, ss, sproc, n_gh=ngh)[0]
    sol = root(f, x0, method="hybr", tol=1e-13)
    r, o = PM.point_residuals(S, 0, sol.x, rules, c, ss, sproc, n_gh=ngh)
    print(f"node {node} cap 200, smoothing {eps:.0e}: buys {100 * o['m_cb'] / B:5.1f}% of B, spread {_sov_bp(o, cal):6.1f}, "
          f"b_DD' {100 * o['b_D_D_new'] / B:5.1f}%, max|r| {np.max(np.abs(r)):.0e} (eq {np.argmax(np.abs(r)) + 1})")
