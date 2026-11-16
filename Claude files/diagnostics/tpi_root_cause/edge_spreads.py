# THE NO-TPI SPREAD AT EACH TOP-SLICE EDGE CONFIGURATION, READ AT ONE p^d ON THE SOLVED 6.2%-BOX RULES, AND AT THE
# STOCHASTIC REST POINT WITH THE SAME p^d (THE IRF'S IMPACT STATE).
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))
from rung_anatomy import economy
from box_top_scan import node, rules, g, top, cal, ss, sproc
from global_projection.solver_recursive.recursive_experiment import _sov_bp, stochastic_rest_point
from global_projection.solver_recursive.state_grid import IS, STATE_NAMES

pd = float(sys.argv[2]) if len(sys.argv) > 2 else 0.02
u = g.to_unit(g.points)
rows = []
for i in top:
    S = g.points[i].copy(); S[IS] = np.log(pd / (1 - pd))
    o, r = node(S, False)
    edge = ",".join(f"{STATE_NAMES[j]}{'+' if u[i, j] > 0 else '-'}" for j in range(12) if abs(u[i, j]) > 0.99 and j != IS)
    rows.append((_sov_bp(o, cal), edge or "centre"))
print(f"p^d {100 * pd:.2f}%: " + "  ".join(f"{e} {s:.0f}" for s, e in sorted(rows, reverse=True)))
S = stochastic_rest_point(rules, dict(cal, tpi_on=False), ss, sproc, verbose=False)
S[IS] = np.log(pd / (1 - pd))
o, r = node(S, False)
print(f"rest point + shock: {_sov_bp(o, cal):.0f} bp (resid {r:.0e})")
