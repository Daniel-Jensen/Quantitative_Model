# WHAT THE CONTINUATION SEES WHEN NEXT PERIOD'S EUROSYSTEM BOOK PASSES THE M EDGE: THE TRUE (b_DD', M') AGAINST THE
# STATE THE RULES ARE READ AT (grid.clip, IN THE SHEARED BOX COORDINATES). Node 42, rules of the 602.6 bp rung.
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))
from rung_anatomy import economy
from global_projection.solver_recursive.state_grid import IBDD, IM

ROOT = HERE.parents[2]
cal, ss, sproc = economy()
rules = pickle.load(open(ROOT / "results" / "GLOBAL_tpi_walk" / "tpi_rung_00060260.pkl", "rb"))["rules"]
g, B = rules.grid, cal["B_gov_D_ss"]
S = g.points[42].copy()
tot = 0.8604 * B  # the D bank's book plus the Eurosystem's, held fixed: purchases move bonds from one to the other
print("  bought%B   true b_DD'%B  true M'%B  |  read b_DD'%B  read M'%B")
for buy in (0.0, 0.10, 0.20, 0.30, 0.40, 0.60, 0.86):
    Sn = S.copy(); Sn[IM] = buy * B; Sn[IBDD] = tot - buy * B
    Sr = g.clip(Sn[None, :])[0]
    print(f"  {100 * buy:7.0f}  {100 * Sn[IBDD] / B:12.1f} {100 * Sn[IM] / B:10.1f}  |  {100 * Sr[IBDD] / B:12.1f}"
          f" {100 * Sr[IM] / B:10.1f}")
