"Capital block: Cobb-Douglas capital demand and the Jermann (1998) convex capital-adjustment cost (capital producer)."
import numpy as np


def gamma_params(cal, country="D"):
    # JERMANN (1998) CAPITAL-ADJUSTMENT COST PARAMETERS
    delta = cal[f"delta_{country}"]
    ksi   = cal[f"ksi_{country}"]
    gamma0 = delta ** ksi / (1 - ksi)
    gamma1 = -delta * ksi / (1 - ksi)
    return gamma0, gamma1


def capital_demand(rk_ss, mc_ss, cal, country="D"):
    # INVERTED CAPITAL DEMAND FROM COBB-DOUGLAS

    alpha = cal[f"alpha_{country}"]
    delta = cal[f"delta_{country}"]
    Z_ss  = cal[f"Z_ss_{country}"]
    
    return (mc_ss * alpha * Z_ss / (rk_ss + delta)) ** (1 / (1 - alpha))


def solve_capital_path(Kap_path, Kap_ss, Q_ss, mpk_path, cal, country="D"):
    # CAPITAL ALONG THE TRANSTION 

    delta  = cal[f"delta_{country}"]
    ksi    = cal[f"ksi_{country}"]
    # get the jerman adjustment cost parameters
    gamma0, gamma1 = gamma_params(cal, country)
    
    T = len(Kap_path)

    Kap_lag = np.concatenate(([Kap_ss], Kap_path[:-1]))
    bracket = (Kap_path / Kap_lag - (1 - delta) - gamma1) / gamma0

    
    # compute the capital-adjustment cost and the price of capital Q given firm path
    iota  = bracket ** (1 / (1 - ksi))
    Q     = 1.0 / (gamma0 * (1 - ksi) * iota ** (-ksi))

    Q_lag = np.concatenate(([Q_ss], Q[:-1]))
    rk    = (mpk_path + (1 - delta) * Q) / Q_lag - 1

    I          = iota * Kap_lag
    cap_profit = (Q * (Kap_path - (1 - delta) * Kap_lag) - I
                  + mpk_path * (Kap_path - Kap_lag))

    return dict(iota=iota, Q=Q, rk=rk, I=I, cap_profit=cap_profit)
