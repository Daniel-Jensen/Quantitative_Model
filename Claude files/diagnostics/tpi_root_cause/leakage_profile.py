# WHERE IN s THE TPI RULES DIFFER FROM THE NO-TPI RULES (2026-11-06, READ-ONLY): THE CAP BINDS ONLY AT THE TOP s NODE,
# SO ANY GAP FAR BELOW IT IS THE INTERPOLANT CARRYING THE TOP NODE'S PURCHASES ACROSS THE RANGE.
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

RULES = ("Q_bD", "alpha_D", "C_D", "N_D", "Kp_D", "p", "A_D")

if __name__ == "__main__":
    from post_impact import setup
    from global_projection.solver_recursive.state_grid import IS, IM, default_prob
    from global_projection.solver_recursive.decision_rules import SOLVE
    cal, ss, sproc, r_off, r_on, c_off, c_on, run = setup()
    rests = np.load(HERE / "runs" / "exact_rest_points.npy", allow_pickle=True).item()
    g = r_on.grid
    S0 = rests["off"].copy()
    nodes_s = np.unique(np.round(g.points[:, IS], 6))

    print("THE s-AXIS NODES (other states at the box centre): TPI node value / no-TPI node value - 1, in %")
    on_axis = [i for i in range(g.n) if np.all(np.abs(np.delete(g.points_unit[i], IS)) < 1e-12)]
    for i in sorted(on_axis, key=lambda i: g.points[i, IS]):
        s = g.points[i, IS]
        xcb = r_on.vals["x_cb"][0][i]
        print(f"  s {s:7.3f} (pd {100*float(default_prob(s)):5.2f}%) x_cb {xcb:+.4f}  " +
              "  ".join(f"{k} {100*(r_on.vals[k][0][i]/r_off.vals[k][0][i]-1):+.4f}" for k in RULES))

    print("\nTHE FITTED RULES ALONG s AT THE NO-TPI REST POINT'S OTHER STATES: TPI / no-TPI - 1, in %")
    for s in np.linspace(g.lo[IS], g.hi[IS], 27):
        S = S0.copy(); S[IS] = s
        mark = " <- node" if np.min(np.abs(nodes_s - s)) < 1e-3 else ""
        print(f"  s {s:7.3f} (pd {100*float(default_prob(s)):5.2f}%) " +
              "  ".join(f"{k} {100*(float(r_on.eval(k, 0, np.atleast_2d(S))[0])/float(r_off.eval(k, 0, np.atleast_2d(S))[0])-1):+.4f}"
                        for k in RULES) + mark)

    print("\nTHE CROSS NODES s = top x each state = +-1: TPI / no-TPI - 1 for Q_bD and alpha_D, and the purchase x_cb")
    from global_projection.solver_recursive.state_grid import STATE_NAMES
    for i in range(g.n):
        u = g.points_unit[i]
        if abs(u[IS] - 1.0) < 1e-12:
            j = [k for k in range(g.d) if k != IS and abs(u[k]) > 1e-12]
            lab = f"{STATE_NAMES[j[0]]}={u[j[0]]:+.0f}" if j else "axis"
            print(f"  {lab:10s} x_cb {r_on.vals['x_cb'][0][i]:+.4f}  Q_bD {100*(r_on.vals['Q_bD'][0][i]/r_off.vals['Q_bD'][0][i]-1):+.3f}"
                  f"  alpha_D {100*(r_on.vals['alpha_D'][0][i]/r_off.vals['alpha_D'][0][i]-1):+.3f}"
                  f"  C_D {100*(r_on.vals['C_D'][0][i]/r_off.vals['C_D'][0][i]-1):+.3f}")
