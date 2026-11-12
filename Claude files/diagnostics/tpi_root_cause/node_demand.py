# ONE NODE'S PURCHASE AS THE CAP FALLS, RULES HELD AT THE 602.6 bp RUNG (ONE-OFF: THE CONTINUATION IS NOT RE-SOLVED).
# Each cap is solved at the node alone (its 14 unknowns, regime 0) with scipy root; it shows the shape of the private
# demand the TPI leans on, the cap where next period's Eurosystem book leaves the M box, and the corner b_DD' = 0.
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
from global_projection.solver_recursive.state_grid import default_prob, IS
from global_projection.solver_recursive import point_map as PM

ROOT = HERE.parents[2]
cal, ss, sproc = economy()
sv = pickle.load(open(ROOT / "results" / "GLOBAL_tpi_walk" / "tpi_rung_00060260.pkl", "rb"))
rules = sv["rules"]
B = cal["B_gov_D_ss"]
M_half = 0.32 * (B - cal["b_D_F_ss"])
nodes = [int(a) for a in sys.argv[1:]] or [42, 58, 49]
for i in nodes:
    S = rules.grid.points[i]
    x = np.array([rules.vals[k][0][i] for k in SOLVE])
    print(f"\nnode {i}: p^d {100 * default_prob(S[IS]):.2f}%   (M box edge = {100 * M_half / B:.1f}% of B)")
    print("   cap_bp   sov_bp   buy%B   M'/M_half   b_DD'%B  b_DF'%B    foc_D      mu_D   max|r|")
    for cap in (650, 600, 550, 500, 450, 400, 350, 300, 250, 200):
        cal["tpi_cap_bp"] = float(cap)
        f = lambda z: PM.point_residuals(S, 0, z, rules, cal, ss, sproc, n_gh=rules.n_gh or 5)[0]
        sol = root(f, x, method="hybr", tol=1e-13)
        r, o = PM.point_residuals(S, 0, sol.x, rules, cal, ss, sproc, n_gh=rules.n_gh or 5)
        print(f"   {cap:6.0f} {_sov_bp(o, cal):8.1f} {100 * o['m_cb'] / B:7.2f} {o['M_cb_new'] / M_half:10.2f}"
              f" {100 * o['b_D_D_new'] / B:9.2f} {100 * o['b_D_F_new'] / B:8.2f} {o['foc_D']:9.2e} {o['mu_D']:8.4f}"
              f"  {np.max(np.abs(r)):.0e}")
        if np.max(np.abs(r)) < 1e-8:
            x = sol.x
