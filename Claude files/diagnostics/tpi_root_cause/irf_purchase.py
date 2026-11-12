# THE ONE-OFF PURCHASE THAT HOLDS A CAP AT THE IRF'S IMPACT STATE (THE NO-TPI REST POINT WITH p^d AT THE SHOCK), RULES
# HELD: HOW MUCH OF THE DEBT EACH bp OF THE BREACH COSTS WHERE THE ILLUSTRATION HAPPENS. Usage: RULES.pkl [shock]
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
from global_projection.solver_recursive.recursive_experiment import _sov_bp, stochastic_rest_point
from global_projection.solver_recursive.state_grid import IS
from global_projection.solver_recursive import point_map as PM

rules = pickle.load(open(sys.argv[1], "rb"))
shock = float(sys.argv[2]) if len(sys.argv) > 2 else 0.02
cal, ss, sproc = economy()
B, ngh = cal["B_gov_D_ss"], rules.n_gh or 5
S = stochastic_rest_point(rules, dict(cal, tpi_on=False), ss, sproc, verbose=False)
S[IS] = np.log(shock / (1 - shock))
x = np.array([float(rules.eval(k, 0, S[None, :])[0]) for k in SOLVE]); x[-1] = 0.0
print(f"IRF impact state, p^d {100 * shock:.2f}%:   cap -> spread, bought % of D debt, D-bank book % of D debt, n_D vs no-TPI")
n0 = None
for cap in (1e5, 220, 215, 210, 205, 200):
    c = dict(cal, tpi_on=cap < 1e4, tpi_cap_bp=cap)
    f = lambda z: PM.point_residuals(S, 0, z, rules, c, ss, sproc, n_gh=ngh)[0]
    sol = root(f, x, method="hybr", tol=1e-13)
    r, o = PM.point_residuals(S, 0, sol.x, rules, c, ss, sproc, n_gh=ngh)
    n0 = o["n_D"] if n0 is None else n0
    print(f"   {'none' if cap > 1e4 else f'{cap:.0f}':>5}   {_sov_bp(o, cal):6.1f} bp   {100 * o['m_cb'] / B:5.1f}%   "
          f"{100 * o['b_D_D_new'] / B:5.1f}%   {100 * (o['n_D'] / n0 - 1):+5.2f}%   resid {np.max(np.abs(r)):.0e}")
    if np.max(np.abs(r)) < 1e-8:
        x = sol.x
