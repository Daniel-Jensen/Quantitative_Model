# FULL-MODEL DRIVER: 6-STATE, TWO-REGIME SOLVE WITH PRICED SOVEREIGN DEFAULT.
# Seeds both regimes from the restricted perturbation (flat in s), solves the
# two-branch fixed point, and checks the pass-through signature: a higher risk
# factor s (higher priced default probability) lowers the bond price Q_B and,
# through bank net worth, tightens the constraint.
import time

import numpy as np

from calibration import solve_bgp, bgp_report
from smolyak import SmolyakGrid
from perturbation import solve_perturbation
from expectations_full import gauss_hermite_3d, default_prob
from time_iteration_full import time_iterate_full, point_block_full


def full_box(cal, k=0.15, b=0.25, p=0.28, s_lo=-11.0, s_hi=-3.0):
    # 6-DIM BOX. s spans low-risk (pi^d~2e-5) to Bocola's stress (pi^d~5%);
    # the extreme s->0 (pi^d=0.5) region is excluded (rarely visited, and it
    # would force the near-insolvency corner into the grid).
    sd = lambda sig, rho: sig / np.sqrt(1 - rho ** 2)
    K, B, P = cal["K_bg"], cal["B_bg"], cal["P_bg"]
    lo = np.array([(1 - k) * K, (1 - b) * B, (1 - p) * P,
                   cal["gamma"] - 3.5 * sd(cal["sig_z"], cal["rho_z"]),
                   np.log(cal["g_star"]) - 3.5 * sd(cal["sig_g"], cal["rho_g"]),
                   s_lo])
    hi = np.array([(1 + k) * K, (1 + b) * B, (1 + p) * P,
                   cal["gamma"] + 3.5 * sd(cal["sig_z"], cal["rho_z"]),
                   np.log(cal["g_star"]) + 3.5 * sd(cal["sig_g"], cal["rho_g"]),
                   s_hi])
    lo[4], hi[4] = np.exp(lo[4]), np.exp(hi[4])
    return lo, hi


def solve_full(cal, mu=2, quad_n=3, damp=0.3, tol=1e-6, maxit=300, verbose=True):
    # PERTURBATION-SEEDED TWO-REGIME SOLVE WITH DEFAULT ON.
    pert = solve_perturbation(cal)
    lo, hi = full_box(cal)
    grid = SmolyakGrid(lo, hi, mu=mu)
    quad = gauss_hermite_3d(quad_n)
    sol = time_iterate_full(grid, cal, pert, quad=quad, damp=damp, tol=tol,
                            maxit=maxit, verbose=verbose)
    sol["pert"] = pert
    return sol


def pass_through_check(sol, cal):
    # ALONG A LOW-RISK -> HIGH-RISK SLICE (d=0), Q_B SHOULD FALL AS s RISES.
    grid, coef = sol["grid"], sol["coef"]
    s_grid = np.linspace(cal["s_star"] - 3.0, cal["s_star"] + 3.0, 7)
    base = np.array([cal["K_bg"], cal["B_bg"], cal["P_bg"], cal["gamma"],
                     cal["g_star"], cal["s_star"]])
    print("  s        pi^d      Q_b        C")
    QBs = []
    for s in s_grid:
        st = base.copy(); st[5] = s
        QB = grid.eval(coef["QB"][0], st)[0]; C = grid.eval(coef["C"][0], st)[0]
        QBs.append(QB)
        print(f"  {s:+.2f}   {default_prob(s):.4f}   {QB:.5f}   {C:.5f}")
    print(f"  => Q_b monotonically falling in s: "
          f"{all(np.diff(QBs) < 1e-6)}  (Bocola pass-through)")


def main():
    cal = solve_bgp(delta_mode="standard")
    print(bgp_report(cal)); print(f"  pi^d(s*) = {default_prob(cal['s_star']):.5f}\n")
    t0 = time.perf_counter()
    sol = solve_full(cal, damp=0.25, maxit=150, verbose=True)
    print(f"\nfull solve {time.perf_counter()-t0:.0f}s  conv={sol['conv']:.2e}  "
          f"grid_bad_d0={(sol['mu']>=0.1).mean():.0%}")
    print("\npass-through (d=0 slice through the BGP, varying s):")
    pass_through_check(sol, cal)
    np.savez("/Users/Huawei/Quantitative_Model/code/bocola2016/output/"
             "full_solution.npz", lo=sol["grid"].lo, hi=sol["grid"].hi,
             **{f"coef_{k}_{d}": sol["coef"][k][d]
                for k in ("C", "R", "QB", "av") for d in (0, 1)})
    print("\nsaved full solution to output/full_solution.npz")


if __name__ == "__main__":
    main()
