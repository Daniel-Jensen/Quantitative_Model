"""Household block: endogenous grid method (EGM) for a one-asset,
incomplete-markets consumption-savings problem with GHH utility.

GHH composite: x = c − v(N) where v(N) = chi·N^(1+1/frisch)/(1+1/frisch).
Utility: u(x) = x^(1−sigma)/(1−sigma).
Labour supply (static FOC): chi·N^(1/frisch) = w/P_CES (income-effect-free).

For the individual household's deposit problem, aggregate N is taken as
given (competitive labour market with idiosyncratic productivity e). The
EGM Euler equation operates on the composite x, and vN (the aggregate
labour disutility evaluated at the period's N) is passed in from outside.
Setting vN=0 everywhere recovers the standard CRRA-separable case.

State: (a, e). Choice: a' (savings), with c = (1+r)·a + y(e) − a'.
Borrowing constraint: a' ≥ a_min.
"""
import numpy as np


def make_asset_grid(cal, country="D"):
    "Non-uniform asset grid; more points near the borrowing constraint."
    n_a    = cal[f"n_a_{country}"]
    a_min  = cal[f"a_min_{country}"]
    a_max  = cal[f"a_max_{country}"]
    curve  = cal[f"a_curve_{country}"]
    return a_min + (a_max - a_min) * np.linspace(0, 1, n_a) ** curve


def egm_step(c_next, a_grid, Pi, r_today, r_next, y_e, beta, sigma, a_min,
             vN_today=0.0, vN_next=0.0):
    """One backward EGM step.

    vN_today, vN_next: aggregate labour disutility v(N_t) and v(N_{t+1}).
    With GHH: Euler operates on composite x = c − vN. Pass vN=0 for
    standard CRRA-separable behaviour.
    """
    n_a, n_e = c_next.shape

    # GHH composite tomorrow: x_{t+1} = c_{t+1} − v(N_{t+1})
    x_next = np.maximum(c_next - vN_next, 1e-11)

    # Expected marginal utility of composite tomorrow
    Eu_next = (x_next ** (-sigma)) @ Pi.T         # (n_a, n_e)

    # EGM inversion: x_t^{−sigma} = beta·(1+r_next)·E[x_{t+1}^{−sigma}]
    x_endo = (beta * (1 + r_next) * Eu_next) ** (-1 / sigma)  # (n_a, n_e)

    # Recover actual consumption: c_t = x_t + v(N_t)
    c_endo = x_endo + vN_today

    # Endogenous cash-on-hand: m = c_t + a' (implies a today)
    m_endo     = c_endo + a_grid[:, None]
    a_endo     = (m_endo - y_e[None, :]) / (1 + r_today)

    c_today    = np.empty((n_a, n_e))
    a_pol_today = np.empty((n_a, n_e))

    for e in range(n_e):
        c_today[:, e]    = np.interp(a_grid, a_endo[:, e], c_endo[:, e])
        a_pol_today[:, e] = (1 + r_today) * a_grid + y_e[e] - c_today[:, e]

        constrained = a_grid < a_endo[0, e]
        a_pol_today[constrained, e] = a_min
        c_today[constrained, e] = (1 + r_today) * a_grid[constrained] + y_e[e] - a_min

    a_pol_today = np.maximum(a_pol_today, a_min)
    return c_today, a_pol_today


def solve_steady_state_household(a_grid, Pi, r_ss, y_e, beta, sigma, a_min, tol,
                                  maxiter=10_000, vN_ss=0.0):
    """Solve the household problem in steady state via EGM iteration.

    vN_ss: aggregate labour disutility at the steady state (scalar).
    """
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
                               vN_path=None):
    """Backward induction over a finite horizon, terminal condition c_ss.

    r_path: length T+1 (r_path[0]=r_dep_ss pre-shock; r_path[t] for t≥1
            is the rate locked in at t−1, realized at t).
    vN_path: (T,) array of aggregate labour disutility per period, or None
             (treated as zeros, i.e. standard CRRA-separable).
    """
    T = y_path.shape[0]
    n_a, n_e = a_grid.shape[0], y_path.shape[1]

    if vN_path is None:
        vN_path = np.zeros(T)

    c_path    = np.empty((T, n_a, n_e))
    a_pol_path = np.empty((T, n_a, n_e))

    c_next = c_ss
    # Terminal-period GHH disutility: period T is permanently at SS, so use vN_path[-1].
    # Using 0.0 here would mis-state the x composite and make households under-save.
    vN_next = vN_path[-1] if len(vN_path) > 0 else 0.0
    for t in range(T - 1, -1, -1):
        r_today = r_path[t]
        r_next  = r_path[t + 1]
        vN_today = vN_path[t]
        c_t, a_pol_t = egm_step(c_next, a_grid, Pi, r_today, r_next, y_path[t],
                                  beta, sigma, a_min,
                                  vN_today=vN_today, vN_next=vN_next)
        c_path[t]    = c_t
        a_pol_path[t] = a_pol_t
        c_next = c_t
        vN_next = vN_today

    return c_path, a_pol_path
