# CONTINUE THE FULL SOLVE FROM THE SAVED SOLUTION TO TIGHTEN CONVERGENCE.
# The default IRF's flat real block traces to mu reading 0 (the barely-binding
# constraint is below the under-converged solve's resolution). Warm-starting
# and iterating much further should resolve mu (~0.001 at the BGP, rising with
# risk) so the liquidity channel can activate.
import time

import numpy as np

from calibration import solve_bgp
from smolyak import SmolyakGrid
from expectations_full import gauss_hermite_3d, default_prob
from time_iteration_full import (time_iterate_full,
                                 initial_coef_full_from_saved, point_block_full)

OUT = "/Users/Huawei/Quantitative_Model/code/bocola2016/output/full_solution.npz"


def main():
    cal = solve_bgp(delta_mode="standard")
    d = np.load(OUT)
    grid = SmolyakGrid(d["lo"], d["hi"], mu=2)
    coef_saved = {k: (d[f"coef_{k}_0"], d[f"coef_{k}_1"])
                  for k in ("C", "R", "QB", "av")}
    init = initial_coef_full_from_saved(grid, coef_saved)
    quad = gauss_hermite_3d(3)

    t0 = time.perf_counter()
    sol = time_iterate_full(grid, cal, quad=quad, damp=0.2, tol=1e-6,
                            maxit=400, verbose=True, init=init)
    sol["grid"] = grid
    print(f"\ncontinued solve {time.perf_counter()-t0:.0f}s  conv={sol['conv']:.2e}")

    # does mu now resolve? check the d=0 slice through the BGP as s rises
    base = np.array([cal["K_bg"], cal["B_bg"], cal["P_bg"], cal["gamma"],
                     cal["g_star"], cal["s_star"]])
    print(f"\n  mu resolution (BGP mu should be ~{cal['mu_bg']:.4f}):")
    print(f"  {'s':>6} {'pi^d':>7} {'mu':>9} {'lev':>7} {'Q_b':>8}")
    for s in [-9.06, -7.06, -5.06, -4.06, -3.66]:
        st = base.copy(); st[5] = s
        C = grid.eval(sol["coef"]["C"][0], st)[0]
        R = grid.eval(sol["coef"]["R"][0], st)[0]
        QB = grid.eval(sol["coef"]["QB"][0], st)[0]
        _, mu, _, sm = point_block_full([C, R, QB], st, 0, grid, sol["coef"],
                                        cal, quad)
        print(f"  {s:+6.2f} {default_prob(s):7.4f} {mu:9.5f} "
              f"{sm['assets']/sm['N']:7.3f} {QB:8.5f}")

    np.savez(OUT, lo=grid.lo, hi=grid.hi,
             **{f"coef_{k}_{dd}": sol["coef"][k][dd]
                for k in ("C", "R", "QB", "av") for dd in (0, 1)})
    print("\nre-saved tightened full solution.")


if __name__ == "__main__":
    main()
