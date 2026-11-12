# ONE-OFF PURCHASE AGAINST THE CAP AT p^d 2%, FOR THE P_D+ EDGE AND THE CENTRE CONFIGURATIONS (6.2%-BOX NO-TPI RULES):
# HOW MANY % OF THE DEBT EACH bp OF COMPRESSION COSTS AT THIS RISK LEVEL.
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.argv = sys.argv[:1]
import box_top_scan as T
from global_projection.solver_recursive.recursive_experiment import _sov_bp
from global_projection.solver_recursive.state_grid import IS, IP_D

u = T.g.to_unit(T.g.points)
pick = {"centre": [i for i in T.top if np.allclose(np.delete(u[i], IS), 0)][0],
        "P_D+": [i for i in T.top if u[i, IP_D] > 0.99][0]}
for name, i in pick.items():
    S = T.g.points[i].copy(); S[IS] = np.log(0.02 / 0.98)
    print(f"\n{name} at p^d 2%:  cap  ->  buy%B  spread  resid")
    for cap in (300, 280, 260, 240, 220, 200):
        T.cap = float(cap)
        o, r = T.node(S, True)
        print(f"      {cap:4d}      {100 * o['m_cb'] / T.B:6.1f}  {_sov_bp(o, T.cal):6.1f}  {r:.0e}")
