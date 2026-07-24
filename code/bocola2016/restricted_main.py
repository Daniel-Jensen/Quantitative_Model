# RESTRICTED-MODEL DRIVER: PERTURBATION SEED -> GLOBAL SOLVE -> ACCURACY/IRF.
# The no-sovereign-risk model (drop s, d) -- Bocola's own first estimation step
# and the validation benchmark for the full model. Perturbation-seeded global
# solve; accuracy reported on-grid, off-grid, and (what matters) along the
# ergodic path; global-vs-perturbation IRF cross-check.
import time

import numpy as np

from calibration import solve_bgp, bgp_report
from smolyak import SmolyakGrid
from perturbation import solve_perturbation, irf as pirf, NX
from time_iteration import (time_iterate, initial_coef_perturbation,
                            initial_coef_from_coef, point_block)
from simulate import simulate
from rotated_grid import rotated_from_path
from diagnostics import accuracy_report, print_accuracy


def build_box(cal, k=0.15, b=0.25, p=0.28, wide=1.0):
    # STATE BOX AROUND THE BGP. The capital band is WIDE on purpose: capital is
    # a near-unit-root state (perturbation eig 0.98) whose ergodic wandering
    # must stay strictly interior, else box-edge clipping distorts the sim.
    # Corners are infeasible but provably never visited (see interior check).
    sd_z = cal["sig_z"] / np.sqrt(1 - cal["rho_z"] ** 2)
    sd_g = cal["sig_g"] / np.sqrt(1 - cal["rho_g"] ** 2)
    K, B, P = cal["K_bg"], cal["B_bg"], cal["P_bg"]
    lo = np.array([(1 - wide * k) * K, (1 - wide * b) * B, (1 - wide * p) * P,
                   cal["gamma"] - 3.5 * sd_z, np.log(cal["g_star"]) - 3.5 * sd_g])
    hi = np.array([(1 + wide * k) * K, (1 + wide * b) * B, (1 + wide * p) * P,
                   cal["gamma"] + 3.5 * sd_z, np.log(cal["g_star"]) + 3.5 * sd_g])
    lo[4], hi[4] = np.exp(lo[4]), np.exp(hi[4])
    return lo, hi


def solve_on_box(cal, pert, lo, hi, mu=2, damp=0.3, tol=1e-6, maxit=250,
                 verbose=False):
    # ONE PERTURBATION-SEEDED GLOBAL SOLVE ON A GIVEN BOX.
    grid = SmolyakGrid(lo, hi, mu=mu)
    init = initial_coef_perturbation(grid, cal, pert)
    sol = time_iterate(grid, cal, damp=damp, tol=tol, maxit=maxit,
                       verbose=verbose, init=init)
    sol["grid"] = grid
    sol["pert"] = pert
    return sol


def solve_restricted(cal, mu=2, damp=0.3, tol=1e-6, maxit=250, verbose=True):
    # PERTURBATION-SEEDED GLOBAL SOLVE ON A BGP-CENTERED BOX (SINGLE PASS).
    pert = solve_perturbation(cal)
    lo, hi = build_box(cal)
    return solve_on_box(cal, pert, lo, hi, mu, damp, tol, maxit, verbose)


def solve_restricted_adapted(cal, mu=2, rounds=4, margin=0.08, verbose=True):
    # BOX ADAPTATION: TIGHTEN THE BOX TO THE ERGODIC CLOUD, RE-SOLVE, REPEAT.
    # Round 1 seeds from the perturbation (correct ergodic set); later rounds
    # WARM-START from the previous round's rules (fast). Tightens directly to
    # the ergodic [0.3%,99.7%] quantiles + margin, dropping the never-visited
    # corners a rectangular box overcovers.
    pert = solve_perturbation(cal)
    lo, hi = build_box(cal)
    prev = None
    for r in range(1, rounds + 1):
        grid = SmolyakGrid(lo, hi, mu=mu)
        if prev is None:
            init = initial_coef_perturbation(grid, cal, pert)
            maxit = 250
        else:
            init = initial_coef_from_coef(grid, prev["coef"], prev["grid"])
            maxit = 120
        sol = time_iterate(grid, cal, damp=0.3, tol=1e-6, maxit=maxit,
                           verbose=False, init=init)
        sol["grid"] = grid; sol["pert"] = pert
        path = simulate(grid, sol["coef"], cal, T=8000, burn=1500)
        q_lo = np.quantile(path, 0.003, axis=0)
        q_hi = np.quantile(path, 0.997, axis=0)
        span = q_hi - q_lo
        lo_new, hi_new = q_lo - margin * span, q_hi + margin * span
        interior = (path.min(0) >= lo_new).all() and (path.max(0) <= hi_new).all()
        bad = (sol["mu"] >= 0.10).mean()
        if verbose:
            print(f"  box round {r}: conv={sol['conv']:.1e}  grid_bad={bad:.0%}  "
                  f"interior={interior}  "
                  f"P[{lo[2]:.2f},{hi[2]:.2f}]->[{lo_new[2]:.2f},{hi_new[2]:.2f}]",
                  flush=True)
        prev = sol
        if np.allclose(lo_new, lo, rtol=0.03) and np.allclose(hi_new, hi, rtol=0.03):
            break
        lo, hi = lo_new, hi_new
    # final solve on the converged box, warm-started
    grid = SmolyakGrid(lo, hi, mu=mu)
    init = initial_coef_from_coef(grid, prev["coef"], prev["grid"])
    sol = time_iterate(grid, cal, damp=0.3, tol=1e-6, maxit=200,
                       verbose=False, init=init)
    sol["grid"] = grid; sol["pert"] = pert
    return sol


def ergodic_mu(sol, cal, path, n=2000):
    # CONSTRAINT MULTIPLIER ALONG A SAMPLE OF THE ERGODIC PATH.
    grid, coef, quad = sol["grid"], sol["coef"], sol["quad"]
    idx = np.linspace(0, len(path) - 1, n).astype(int)
    mus = np.empty(n)
    for i, k in enumerate(idx):
        s = path[k]
        C = grid.eval(coef["C"], s)[0]; R = grid.eval(coef["R"], s)[0]
        QB = grid.eval(coef["QB"], s)[0]
        _, mus[i], _, _ = point_block([C, R, QB], s, grid, coef, cal, quad)
    return mus


def solve_restricted_rotated(cal, mu=2, rot_rounds=2, verbose=True):
    # DECORRELATED SOLVE: axis-aligned seed -> PCA-rotated grid on the ergodic
    # cloud -> re-solve, aligning the box with the (correlated) ergodic set so
    # the never-visited anti-correlated corners are excluded, not overcovered.
    pert = solve_perturbation(cal)
    lo, hi = build_box(cal)
    grid = SmolyakGrid(lo, hi, mu=mu)
    init = initial_coef_perturbation(grid, cal, pert)
    sol = time_iterate(grid, cal, damp=0.3, tol=1e-6, maxit=250, verbose=False,
                       init=init)
    sol["grid"] = grid; sol["pert"] = pert
    if verbose:
        print(f"  seed (axis-aligned): conv={sol['conv']:.1e}  "
              f"grid_bad={(sol['mu']>=0.1).mean():.0%}", flush=True)

    for r in range(1, rot_rounds + 1):
        path = simulate(sol["grid"], sol["coef"], cal, T=8000, burn=1500)
        rgrid = rotated_from_path(path, mu=mu, k_endog=4.0)
        init = initial_coef_from_coef(rgrid, sol["coef"], sol["grid"])
        # run to FULL convergence on the rotated grid (an under-converged
        # rotated solve simulates worse than the axis-aligned seed)
        rsol = time_iterate(rgrid, cal, damp=0.3, tol=1e-6, maxit=400,
                            verbose=False, init=init)
        rsol["grid"] = rgrid; rsol["pert"] = pert
        sol = rsol
        if verbose:
            print(f"  rotation round {r}: conv={sol['conv']:.1e}  "
                  f"grid_bad={(sol['mu']>=0.1).mean():.0%}", flush=True)
    return sol


def interior_fraction(grid, path):
    # FRACTION OF THE ERGODIC PATH STRICTLY INSIDE THE BOX (per dimension).
    # If ~1.0 the never-visited corners' inaccuracy is provably irrelevant.
    lo, hi = grid.lo, grid.hi
    inside = np.all((path > lo) & (path < hi), axis=1)
    return inside.mean()


def main():
    cal = solve_bgp(delta_mode="standard")
    print(bgp_report(cal)); print()

    # Wide axis-aligned box: best ERGODIC accuracy (path stays interior, so the
    # infeasible corners are never visited). Decorrelation (solve_restricted_-
    # rotated) gives machine-precision ON-grid accuracy and 3x fewer bad corners
    # but slightly clips the near-unit-root capital -> worse ergodic behavior;
    # kept available but not the primary. See docs/bocola2016_replication.md.
    t0 = time.perf_counter()
    sol = solve_restricted(cal, damp=0.3, tol=1e-6, maxit=350, verbose=True)
    print(f"\nsolve {time.perf_counter()-t0:.0f}s  conv={sol['conv']:.2e}  "
          f"grid_bad={(sol['mu']>=0.1).mean():.0%}")

    path = simulate(sol["grid"], sol["coef"], cal, T=12000, burn=2000)
    rep = accuracy_report(sol["grid"], sol["coef"], cal, sol["quad"],
                          n_off=1500, path=path)
    print(); print_accuracy(rep)

    frac = interior_fraction(sol["grid"], path)
    print(f"\nergodic path strictly interior to the box: {frac:.2%}  "
          f"(=> corner inaccuracy is never visited)")
    mus = ergodic_mu(sol, cal, path)
    print(f"ergodic mu: mean {mus.mean():.5f}  median {np.median(mus):.5f}  "
          f"p99 {np.quantile(mus,0.99):.5f}  binding>1e-6 {(mus>1e-6).mean():.0%}")
    print("  (Bocola: mu near zero most of the time, spikes in stress)")

    np.savez("/Users/Huawei/Quantitative_Model/code/bocola2016/output/"
             "restricted_solution.npz",
             lo=sol["grid"].lo, hi=sol["grid"].hi,
             **{f"coef_{k}": v for k, v in sol["coef"].items()})
    print("\nsaved restricted solution to output/restricted_solution.npz")


if __name__ == "__main__":
    main()
