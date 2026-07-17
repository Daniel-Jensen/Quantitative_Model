# **Capital block: Cobb-Douglas capital demand + Jermann (1998) adjustment cost.**
import numpy as np


def gamma_params(cal, country="D"):
    # **Jermann adjustment-cost coefficients (γ0, γ1) pinned to δ, ξ.**
    delta = cal[f"delta_{country}"]
    ksi   = cal[f"ksi_{country}"]
    gamma0 = delta ** ksi / (1 - ksi)
    gamma1 = -delta * ksi / (1 - ksi)
    return gamma0, gamma1


def capital_demand(rk_ss, mc_ss, cal, country="D"):
    # **Steady-state capital stock inverted from the Cobb-Douglas FOC.**
    alpha = cal[f"alpha_{country}"]
    delta = cal[f"delta_{country}"]
    Z_ss  = cal[f"Z_ss_{country}"]
    return (mc_ss * alpha * Z_ss / (rk_ss + delta)) ** (1 / (1 - alpha))


def solve_capital_path(Kap_path, Kap_ss, Q_ss, mpk_path, cal, country="D",
                       quality0=1.0, Kap_lag_path=None):
    # **Investment, capital price Q, and realized return rk along the path.**
    # Timing (Bocola eq. 6): Kap_path[t] is bought/priced at t and produces at
    # t+1; Kap_lag_path[t] is the stock carried INTO t — it is both the Jermann
    # rebuilding base and the production stock (mpk_path must be computed on
    # it). quality0 < 1: GK capital-quality loss at t=0 (branch only) —
    # destroys a fraction 1-quality0 of the incoming stock and scales the t=0
    # claim return.
    delta  = cal[f"delta_{country}"]
    ksi    = cal[f"ksi_{country}"]
    gamma0, gamma1 = gamma_params(cal, country)

    if Kap_lag_path is None:
        Kap_lag_path = np.concatenate(([quality0 * Kap_ss], Kap_path[:-1]))
    bracket = (Kap_path / Kap_lag_path - (1 - delta) - gamma1) / gamma0

    # negative bracket → NaN powers; raise so the outer solver penalizes the guess
    if np.any(bracket < 0):
        raise ValueError(
            f"[{country}] Jermann inversion: negative bracket "
            "(capital falling faster than adjustment-cost technology allows)"
        )

    iota  = bracket ** (1 / (1 - ksi))
    Q     = 1.0 / (gamma0 * (1 - ksi) * iota ** (-ksi))   # marginal Tobin's Q

    Q_lag = np.concatenate(([Q_ss], Q[:-1]))
    rk    = (mpk_path + (1 - delta) * Q) / Q_lag - 1
    if quality0 != 1.0:
        rk[0] = quality0 * (mpk_path[0] + (1 - delta) * Q[0]) / Q_lag[0] - 1

    I          = iota * Kap_lag_path
    # capital producers' rents (Jermann): value of installed new capital − cost.
    # No mpk reconciliation term: firms rent exactly the bank-held vintage.
    cap_profit = Q * (Kap_path - (1 - delta) * Kap_lag_path) - I

    return dict(iota=iota, Q=Q, rk=rk, I=I, cap_profit=cap_profit)
