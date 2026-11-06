# THE TPI AGAINST NO TPI ALONG THE RISK AXIS (OTHER STATES AT THE BOX CENTRE, NO EUROSYSTEM BOOK): THE FITTED SPREAD AT THE
# s NODES AND BETWEEN THEM, TO SEE WHETHER THE CAP'S EFFECT AT THE TOP LEAKS DOWN THE s POLYNOMIAL.
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from global_projection.solver_recursive.state_grid import IS, default_prob

off = pickle.load(open(ROOT / "results" / "GLOBAL" / "data" / "rules_risk.pkl", "rb"))
on = pickle.load(open(ROOT / "results" / "GLOBAL" / "data" / "rules_tpi.pkl", "rb"))
dD, dF = 0.056, 0.039
sp = lambda QD, QF: 4e4 * (dD * (1 / QD - 1) - dF * (1 / QF - 1))
g = on.grid
u = g.to_unit(g.points)
centre = g.points[int(np.argmin(np.abs(u).sum(axis=1)))]
s_nodes = sorted(set(np.round(g.points[np.all(np.abs(np.delete(u, IS, axis=1)) < 1e-9, axis=1), IS], 10)))
print("   s        p^d%   node?   spread no TPI   spread TPI   diff (bp)")
for s in np.linspace(g.points[:, IS].min(), g.points[:, IS].max(), 17):
    S = centre.copy(); S[IS] = s
    a = sp(off.eval("Q_bD", 0, S[None, :])[0], off.eval("Q_bF", 0, S[None, :])[0])
    b = sp(on.eval("Q_bD", 0, S[None, :])[0], on.eval("Q_bF", 0, S[None, :])[0])
    node = "node" if np.min(np.abs(np.array(s_nodes) - s)) < 1e-6 else ""
    print(f"{s:8.3f}  {100 * default_prob(s):6.3f}  {node:5s}   {a:12.1f}   {b:10.1f}   {b - a:+8.1f}")
for s in s_nodes:
    S = centre.copy(); S[IS] = s
    a = sp(off.eval("Q_bD", 0, S[None, :])[0], off.eval("Q_bF", 0, S[None, :])[0])
    b = sp(on.eval("Q_bD", 0, S[None, :])[0], on.eval("Q_bF", 0, S[None, :])[0])
    print(f"node s {s:8.3f} p^d {100 * default_prob(s):6.3f}%: no TPI {a:6.1f}  TPI {b:6.1f}  diff {b - a:+6.1f}")
