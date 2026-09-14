# HOUSEHOLD BLOCK: EGM FOR A ONE-ASSET INCOMPLETE-MARKETS PROBLEM WITH GHH UTILITY.
import numpy as np


def make_asset_grid(cal, country="D"):
    # NON-UNIFORM ASSET GRID, DENSER NEAR THE BORROWING CONSTRAINT.
    n_a   = cal[f"n_a_{country}"]
    a_min = cal[f"a_min_{country}"]
    a_max = cal[f"a_max_{country}"]
    curve = cal[f"a_curve_{country}"]
    return a_min + (a_max - a_min) * np.linspace(0, 1, n_a) ** curve


def egm_step(c_next, a_grid, Pi, r_today, r_next, y_e, beta, sigma, a_min,
             vN_today=0.0, vN_next=0.0):
    # ONE BACKWARD EGM STEP: GIVEN c_{t+1}, SOLVE FOR c_t AND THE SAVINGS POLICY.
    n_a, n_e = c_next.shape

    x_next  = np.maximum(c_next - vN_next, 1e-11)              # GHH composite tomorrow
    Eu_next = (x_next ** (-sigma)) @ Pi.T
    x_endo  = (beta * (1 + r_next) * Eu_next) ** (-1 / sigma)  # Euler inversion
    c_endo  = x_endo + vN_today

    m_endo = c_endo + a_grid[:, None]                          # endogenous cash-on-hand
    a_endo = (m_endo - y_e[None, :]) / (1 + r_today)

    c_today     = np.empty((n_a, n_e))
    a_pol_today = np.empty((n_a, n_e))

    for e in range(n_e):
        c_today[:, e]     = np.interp(a_grid, a_endo[:, e], c_endo[:, e])
        a_pol_today[:, e] = (1 + r_today) * a_grid + y_e[e] - c_today[:, e]

        constrained = a_grid < a_endo[0, e]                    # borrowing constraint binds
        a_pol_today[constrained, e] = a_min
        c_today[constrained, e] = (1 + r_today) * a_grid[constrained] + y_e[e] - a_min

    return c_today, np.maximum(a_pol_today, a_min)


def solve_steady_state_household(a_grid, Pi, r_ss, y_e, beta, sigma, a_min, tol,
                                 maxiter=10_000, vN_ss=0.0):
    # STEADY-STATE HOUSEHOLD POLICIES BY EGM FIXED-POINT ITERATION.
    c = np.maximum((1 + r_ss) * a_grid[:, None] + y_e[None, :] - a_grid[:, None], 1e-11)

    for _ in range(maxiter):
        c_new, a_pol = egm_step(c, a_grid, Pi, r_ss, r_ss, y_e, beta, sigma, a_min,
                                vN_today=vN_ss, vN_next=vN_ss)
        diff = np.max(np.abs(c_new - c))
        c = c_new
        if diff < tol:
            break
    else:
        raise RuntimeError(f"Household EGM did not converge (diff={diff:.2e})")
    return c, a_pol


def solve_backward_transition(a_grid, Pi, r_path, y_path, c_ss, beta, sigma, a_min,
                              vN_path=None, use_fast=True):
    # BACKWARD INDUCTION OVER THE TRANSITION (NUMBA KERNEL OR NUMPY), TERMINAL c_ss.
    from blocks import fast_kernels

    T = y_path.shape[0]
    n_a, n_e = a_grid.shape[0], y_path.shape[1]

    if vN_path is None:
        vN_path = np.zeros(T)

    if use_fast and fast_kernels.HAVE_NUMBA:
        return fast_kernels.hh_backward(
            np.ascontiguousarray(a_grid), np.ascontiguousarray(Pi.T),
            np.ascontiguousarray(r_path, dtype=float),
            np.ascontiguousarray(y_path),
            np.ascontiguousarray(c_ss), float(beta), float(sigma),
            float(a_min), np.ascontiguousarray(vN_path, dtype=float))

    c_path     = np.empty((T, n_a, n_e))
    a_pol_path = np.empty((T, n_a, n_e))

    c_next  = c_ss
    vN_next = vN_path[-1] if len(vN_path) > 0 else 0.0   # terminal period at SS
    for t in range(T - 1, -1, -1):
        c_t, a_pol_t = egm_step(c_next, a_grid, Pi, r_path[t], r_path[t + 1],
                                y_path[t], beta, sigma, a_min,
                                vN_today=vN_path[t], vN_next=vN_next)
        c_path[t]     = c_t
        a_pol_path[t] = a_pol_t
        c_next  = c_t
        vN_next = vN_path[t]

    return c_path, a_pol_path
