# CAPITAL BLOCK: COBB-DOUGLAS CAPITAL DEMAND + JERMANN (1998) ADJUSTMENT COST.
import numpy as np


def gamma_params(cal, country="D"):
    # JERMANN ADJUSTMENT-COST COEFFICIENTS (gamma0, gamma1) PINNED TO delta, ksi.
    delta = cal[f"delta_{country}"]
    ksi   = cal[f"ksi_{country}"]
    gamma0 = delta ** ksi / (1 - ksi)
    gamma1 = -delta * ksi / (1 - ksi)
    return gamma0, gamma1


def capital_demand(rk_ss, mc_ss, cal, country="D"):
    # STEADY-STATE CAPITAL STOCK INVERTED FROM THE COBB-DOUGLAS FOC.
    alpha = cal[f"alpha_{country}"]
    delta = cal[f"delta_{country}"]
    Z_ss  = cal[f"Z_ss_{country}"]
    return (mc_ss * alpha * Z_ss / (rk_ss + delta)) ** (1 / (1 - alpha))


def solve_capital_path(Kap_path, Kap_lag0, Q_lag0, mpk_path, cal, country="D",
                       Kap_lag_path=None):
    # INVESTMENT, CAPITAL PRICE Q, AND REALIZED RETURN rk ALONG THE PATH.
    # Timing (Bocola eq. 6): Kap_path[t] is bought/priced at t and produces at
    # t+1; Kap_lag_path[t] is the stock carried INTO t, which is both the
    # Jermann rebuilding base and the production stock mpk_path was built on.
    delta  = cal[f"delta_{country}"]
    ksi    = cal[f"ksi_{country}"]
    gamma0, gamma1 = gamma_params(cal, country)

    if Kap_lag_path is None:
        Kap_lag_path = np.concatenate(([Kap_lag0], Kap_path[:-1]))
    bracket = (Kap_path / Kap_lag_path - (1 - delta) - gamma1) / gamma0

    # negative bracket -> NaN powers; raise so the outer solver penalizes the guess
    if np.any(bracket < 0):
        raise ValueError(
            f"[{country}] Jermann inversion: negative bracket "
            "(capital falling faster than adjustment-cost technology allows)"
        )

    iota  = bracket ** (1 / (1 - ksi))
    Q     = 1.0 / (gamma0 * (1 - ksi) * iota ** (-ksi))   # marginal Tobin's Q

    Q_lag = np.concatenate(([Q_lag0], Q[:-1]))
    rk    = (mpk_path + (1 - delta) * Q) / Q_lag - 1
    I     = iota * Kap_lag_path
    # capital producers' rents: value of installed new capital - cost. No mpk
    # reconciliation term, because firms rent exactly the bank-held vintage.
    cap_profit = Q * (Kap_path - (1 - delta) * Kap_lag_path) - I

    return dict(iota=iota, Q=Q, rk=rk, I=I, cap_profit=cap_profit)
