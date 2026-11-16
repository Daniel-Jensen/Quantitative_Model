# WHERE THE SMOOTHING WALK BITES: THE 200 bp CHECKPOINT (SOLVED AT ITS SAVED SMOOTHING) EVALUATED AT A SMALLER ONE, NODE BY
# NODE: WHICH NODES AND EQUATIONS MOVE, AND THE STATE OF THEIR TWO COMPLEMENTARITY PAIRS. Usage: [eps] (default 5e-3)
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))
from rung_anatomy import economy
from global_projection.solver_recursive.decision_rules import SOLVE
from global_projection.solver_recursive.recursive_experiment import _sov_bp
from global_projection.solver_recursive.state_grid import IS, IM, STATE_NAMES, default_prob
from global_projection.solver_recursive import point_map as PM

EQ = ("cap_D", "cap_F", "lab_D", "lab_F", "euler_D", "uip", "goods_D", "bondD_D", "bondD_F", "euler_F", "dep_clear",
      "bondF_F", "bondF_D", "tpi")
sv = pickle.load(open(HERE.parents[2] / "results" / "GLOBAL_tpi_walk" / "tpi_rung_00020000.pkl", "rb"))
rules, eps0 = sv["rules"], sv["eps"]
eps1 = float(sys.argv[1]) if len(sys.argv) > 1 else 5e-3
cal, ss, sproc = economy()
cal["tpi_cap_bp"] = 200.0
g, B, ngh = rules.grid, cal["B_gov_D_ss"], rules.n_gh or 5
u = g.to_unit(g.points)
print(f"checkpoint cap {sv['cap']} bp at smoothing {eps0}; re-evaluated at {eps1}")
rows = []
for i in range(g.n):
    x = np.array([rules.vals[k][0][i] for k in SOLVE])
    r0, o = PM.point_residuals(g.points[i], 0, x, rules, dict(cal, tpi_eps=eps0), ss, sproc, n_gh=ngh)
    r1, _ = PM.point_residuals(g.points[i], 0, x, rules, dict(cal, tpi_eps=eps1), ss, sproc, n_gh=ngh)
    j = int(np.argmax(np.abs(r1)))
    edge = ",".join(f"{STATE_NAMES[k]}{'+' if u[i, k] > 0 else '-'}" for k in range(12) if abs(u[i, k]) > 0.99 and k != IS)
    rows.append((np.max(np.abs(r1)), i, EQ[j], 100 * default_prob(g.points[i, IS]), edge or "centre", _sov_bp(o, cal),
                 100 * o["m_cb"] / B, (o["Q_bD"] - o["Q_floor"]) / ss["Q_bD_ss"], 100 * o["b_D_D_new"] / B, o["foc_D"],
                 np.max(np.abs(r0))))
print(" node  p^d%  edge        spread  buy%B   price gap   b_DD'%B    foc_D      |r| at 5e-3 (eq)   |r| at saved")
for r in sorted(rows, reverse=True)[:12]:
    res, i, eq, pdq, edge, sp, m, gap, bdd, foc, r0 = r
    print(f"  {i:2d}  {pdq:5.2f}  {edge:10s} {sp:6.1f}  {m:6.2f}  {gap:+.2e}  {bdd:7.2f}  {foc:+.2e}   {res:.1e} ({eq})   {r0:.0e}")
