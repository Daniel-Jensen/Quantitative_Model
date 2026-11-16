# HOW WIDE CAN THE EUROSYSTEM-BOOK BAND GO: AT THE M+ EDGE STATE OF A CANDIDATE BAND (BOOK AT band x THE D BANKS' SS BOOK,
# THE BANK HOLDING THE REST, THE CLAIM ON THE SHEAR'S CENTRE LINE), THE NO-TPI PERIOD MAP SOLVED AT THE NODE ON THE CLEAN
# 0.65-BOX RULES: THE D BANK'S NEXT BOOK b_DD' IN EACH REGIME (NEGATIVE = THE NO-TPI CORNER THAT BROKE THE 1.15 BOX).
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
from global_projection.solver_recursive.state_grid import IS, IM, IBDD, IO
from global_projection.solver_recursive import point_map as PM

rules = pickle.load(open(HERE.parents[2] / "results" / "GLOBAL_cache" / "risk_80883b6daaf8.pkl", "rb"))
cal, ss, sproc = economy()
c_off = dict(cal, tpi_on=False)
g, B, ngh = rules.grid, cal["B_gov_D_ss"], rules.n_gh or 5
c_bDD = cal["B_gov_D_ss"] - cal["b_D_F_ss"]
centre = g.points[int(np.argmin(np.abs(g.to_unit(g.points)).sum(axis=1)))]
print("band  M edge %B   s at   regime 0: b_DD' %B   regime 1: b_DD' %B")
for band in (0.65, 0.8, 0.9, 0.95, 1.0, 1.15):
    for pd in (0.001, 0.02):
        S = centre.copy()
        S[IM] = band * c_bDD
        S[IBDD] = centre[IBDD] - S[IM]
        S[IO] = centre[IO] + (1.0 + cal["r_dep_D_target"]) * ss["Q_bD_ss"] * S[IM]
        S[IS] = np.log(pd / (1 - pd))
        out = []
        for d in (0, 1):
            x = np.array([float(rules.eval(k, d, S[None, :])[0]) for k in SOLVE]); x[-1] = 0.0
            f = lambda z: PM.point_residuals(S, d, z, rules, c_off, ss, sproc, n_gh=ngh)[0][:-1]
            sol = root(lambda z: f(np.append(z, 0.0)), x[:-1], method="hybr", tol=1e-12)
            o = PM.point_residuals(S, d, np.append(sol.x, 0.0), rules, c_off, ss, sproc, n_gh=ngh)[1]
            raw = o["Bp_D"] - o["b_D_F_new"] - o["M_cb_new"]
            out.append(f"{100 * raw / B:7.1f} (resid {np.max(np.abs(sol.fun)):.0e})")
        print(f"{band:4.2f}  {100 * band * c_bDD / B:8.1f}   {100 * pd:4.1f}%   {out[0]:>22}   {out[1]:>22}")
