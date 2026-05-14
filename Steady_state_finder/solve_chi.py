from scipy.optimize import least_squares
from types import SimpleNamespace
import numpy as np


def solve_chi(model, calibration, A_target, B_target,
              xtol=1e-3, tol=1e-4, verbose=True):

    beta = 0.99
    sol = None
    chi0_init = calibration['chi0']   # fixed starting point — reset for every beta
    chi1_init = calibration['chi1']

    while beta >= 0.8:
        calibration['beta'] = beta

        if verbose:
            print(f"\nTrying beta = {beta:.3f}")

        n_success = [0]  # count feasible evaluations inside res

        def res(x):
            chi0, chi1 = x           # beta is NOT in x — held fixed by the while loop
            calibration['chi0'] = chi0
            calibration['chi1'] = chi1

            try:
                ss = model.steady_state(calibration)
                A, B = ss['A'], ss['B']
                n_success[0] += 1
                if verbose:
                    print(f"  chi0={chi0:.4f}, chi1={chi1:.4f}, beta={beta:.4f} | A={A:.4f}, B={B:.4f}")
            except Exception:
                if verbose:
                    print(f"  chi0={chi0:.4f}, chi1={chi1:.4f}, beta={beta:.4f} | FAILED")
                return [1e6, 1e6]

            return [A - A_target, B - B_target]

        x0 = [chi0_init, chi1_init]
        n_success[0] = 0

        r0 = np.array(res(x0))
        if np.linalg.norm(r0) < tol:
            if verbose:
                print(f"  Initial guess already satisfies tolerance — skipping optimizer")
            sol = SimpleNamespace(x=np.array(x0), fun=r0)
            break

        try:
            sol = least_squares(res, x0, xtol=xtol)
        except Exception:
            if verbose:
                print(f"  Solver crashed for beta={beta:.4f}, decreasing beta by 0.001")
            beta -= 0.001
            continue

        # If every evaluation failed, the adjustment cost surface is entirely
        # infeasible at this beta — decrease immediately without checking norm.
        if n_success[0] == 0:
            if verbose:
                print(f"  No feasible adjustment costs at beta={beta:.4f}, decreasing beta by 0.001")
            beta -= 0.001
            continue

        if np.linalg.norm(sol.fun) < tol:
            break

        beta -= 0.001

    # assign final values
    chi0, chi1 = sol.x
    calibration['chi0'] = chi0
    calibration['chi1'] = chi1
    calibration['beta'] = beta

    ss = model.steady_state(calibration)

    if verbose:
        print("\nFinal:")
        print(f"chi0 = {chi0:.6f}")
        print(f"chi1 = {chi1:.6f}")
        print(f"beta = {beta:.6f}")
        print(f"A = {ss['A']:.6f} (target {A_target:.6f})")
        print(f"B = {ss['B']:.6f} (target {B_target:.6f})")

    return calibration, ss