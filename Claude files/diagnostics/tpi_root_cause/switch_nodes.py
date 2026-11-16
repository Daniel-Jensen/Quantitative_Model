# WHICH NODES SWITCH NEXT AT A SAVED RUNG, AND WHERE THEY SIT IN THE BOX (EDGE CODES: STATE+/- IN UNIT COORDINATES).
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
from global_projection.solver_recursive.state_grid import default_prob, IS, STATE_NAMES
from global_projection.solver_recursive import point_map as PM

ROOT = HERE.parents[2]
cal, ss, sproc = economy()
for name in sys.argv[1:] or ["tpi_rung_00059539.pkl"]:
    sv = pickle.load(open(ROOT / "results" / "GLOBAL_tpi_walk" / name, "rb"))
    rules, cap = sv["rules"], sv["cap"]
    cal["tpi_cap_bp"] = cap
    g = rules.grid
    u = g.to_unit(g.points)
    print(f"\nrung {cap:.2f} bp: nodes by no-TPI-style spread (top 12), with their box position")
    rows = []
    for i in range(g.n):
        x = np.array([rules.vals[k][0][i] for k in SOLVE])
        _, o = PM.point_residuals(g.points[i], 0, x, rules, cal, ss, sproc, n_gh=rules.n_gh or 5)
        edge = ",".join(f"{STATE_NAMES[j]}{'+' if u[i, j] > 0 else '-'}" for j in range(u.shape[1])
                        if abs(u[i, j]) > 0.99 and j != IS)
        rows.append((_sov_bp(o, cal), i, 100 * default_prob(g.points[i][IS]), 100 * o["m_cb"] / cal["B_gov_D_ss"],
                     edge or "centre"))
    for sp, i, pdq, m, edge in sorted(rows, reverse=True)[:14]:
        print(f"  node {i:2d}  p^d {pdq:5.2f}%  spread {sp:7.1f} bp  buying {m:6.2f}% of B   edge: {edge}")
    n_top = sum(1 for r in rows if r[2] > 6.0)
    n_47 = sum(1 for r in rows if 4.0 < r[2] < 6.0)
    print(f"  nodes at p^d 6.2%: {n_top}; at 4.7%: {n_47}")
