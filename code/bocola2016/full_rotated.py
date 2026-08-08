# DECORRELATED FULL SOLVE: PCA-ROTATED 6-STATE GRID TO CUT INFEASIBLE CORNERS
# AND RESOLVE THE BARELY-BINDING CONSTRAINT (mu ~ 0.001, below the axis-aligned
# solve's resolution). Rotates the three endogenous states (K,B,P) on the
# ergodic covariance; keeps physical Dz, g, s (the s box is held to the solved
# range). Warm-starts from the saved axis-aligned solution.
import time

import numpy as np

from calibration import solve_bgp
from smolyak import SmolyakGrid
from rotated_grid import rotated_from_path
from expectations_full import gauss_hermite_3d, default_prob
from time_iteration_full import (time_iterate_full, point_block_full,
                                 initial_coef_full_from_saved,
                                 initial_coef_full)
from simulate_full import simulate_full

OUT = "/Users/Huawei/Quantitative_Model/code/bocola2016/output/full_solution.npz"


def _seed_rotated(rgrid, coef_saved, old_grid):
    # WARM START THE ROTATED GRID FROM THE SAVED AXIS-ALIGNED RULES (BOTH REGIMES).
    pts = old_grid.clip(rgrid.points)
    rv = {k: (old_grid.eval(coef_saved[k][0], pts),
              old_grid.eval(coef_saved[k][1], pts))
          for k in ("C", "R", "QB", "av")}
    coef = {k: (rgrid.fit(rv[k][0]), rgrid.fit(rv[k][1])) for k in rv}
    return coef, {k: (rv[k][0].copy(), rv[k][1].copy()) for k in rv}


def main():
    cal = solve_bgp(delta_mode="standard")
    d = np.load(OUT)
    a_grid = SmolyakGrid(d["lo"], d["hi"], mu=2)
    coef_saved = {k: (d[f"coef_{k}_0"], d[f"coef_{k}_1"])
                  for k in ("C", "R", "QB", "av")}
    quad = gauss_hermite_3d(3)

    # ergodic cloud of the axis-aligned solution -> PCA-rotated grid.
    # s box covers the ergodic s spread (incl. high-risk tail) + a margin, so
    # the persistent s process is not clipped; capped at s=-1.5 (pi^d~0.18).
    path = simulate_full(dict(grid=a_grid, coef=coef_saved), cal, T=6000, burn=1000)
    s_span = path[:, 5].max() - path[:, 5].min()
    exog_lo = np.array([path[:, 3].mean() - 3.5 * path[:, 3].std(),
                        path[:, 4].mean() - 3.5 * path[:, 4].std(),
                        path[:, 5].min() - 0.1 * s_span])
    exog_hi = np.array([path[:, 3].mean() + 3.5 * path[:, 3].std(),
                        path[:, 4].mean() + 3.5 * path[:, 4].std(),
                        min(path[:, 5].max() + 0.1 * s_span, -1.5)])
    rgrid = rotated_from_path(path, mu=2, k_endog=4.0, exog_lo=exog_lo,
                              exog_hi=exog_hi)

    init = _seed_rotated(rgrid, coef_saved, a_grid)
    t0 = time.perf_counter()
    sol = time_iterate_full(rgrid, cal, quad=quad, damp=0.25, tol=1e-6,
                            maxit=120, verbose=True, init=init)
    sol["grid"] = rgrid
    bad = (sol["mu"] >= 0.10).mean()
    print(f"\nrotated full solve {time.perf_counter()-t0:.0f}s  "
          f"conv={sol['conv']:.2e}  grid_bad_d0={bad:.0%}")

    base = np.array([cal["K_bg"], cal["B_bg"], cal["P_bg"], cal["gamma"],
                     cal["g_star"], cal["s_star"]])
    print(f"\n  mu resolution (BGP mu ~ {cal['mu_bg']:.4f}):")
    print(f"  {'s':>6} {'pi^d':>7} {'mu':>9} {'lev':>7} {'Q_b':>8}")
    for s in [-9.06, -7.06, -5.06, -4.06, -3.66]:
        st = base.copy(); st[5] = s
        C = rgrid.eval(sol["coef"]["C"][0], st)[0]
        R = rgrid.eval(sol["coef"]["R"][0], st)[0]
        QB = rgrid.eval(sol["coef"]["QB"][0], st)[0]
        _, mu, _, sm = point_block_full([C, R, QB], st, 0, rgrid, sol["coef"],
                                        cal, quad)
        print(f"  {s:+6.2f} {default_prob(s):7.4f} {mu:9.5f} "
              f"{sm['assets']/sm['N']:7.3f} {QB:8.5f}")

    np.savez("/Users/Huawei/Quantitative_Model/code/bocola2016/output/"
             "full_rotated_solution.npz",
             mean=rgrid.mean, V=rgrid.V, lo=rgrid._g.lo, hi=rgrid._g.hi,
             **{f"coef_{k}_{dd}": sol["coef"][k][dd]
                for k in ("C", "R", "QB", "av") for dd in (0, 1)})
    print("\nsaved rotated full solution.")


if __name__ == "__main__":
    main()
