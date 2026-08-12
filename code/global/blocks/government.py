# GOVERNMENT BLOCK: HATCHONDO-MARTINEZ PERPETUITY BONDS, BOHN (1998) FISCAL RULE.
# Default risk is EXOGENOUS (Bocola 2016 eqs. 11-12): the priced default
# probability is an input path to the transition solver, never a function of
# the debt stock. Debt still evolves endogenously under the Bohn tax.
import numpy as np


def govt_steady_state(cal, rdep_ss, country):
    # STEADY-STATE GOVERNMENT BLOCK (NO DEFAULT, CONSTANT DEBT STOCK).
    delta_b  = cal[f"delta_b_{country}"]
    B_gov_ss = cal[f"B_gov_{country}_ss"]
    G        = cal[f"G_{country}"]

    Q_B_ss = delta_b / (rdep_ss + delta_b)
    Tax_ss = G + delta_b * B_gov_ss * (1.0 - Q_B_ss)   # G + coupon = Tax + issuance
    return dict(Q_B_ss=Q_B_ss, Tax_ss=Tax_ss, b_gov_ss=B_gov_ss)


def govt_transition(cal, gs, Q_B_path, def_real_path, country, b_gov0=None,
                    b_anchor=None, recap_path=None):
    # FORWARD-INTEGRATE THE DEBT STOCK UNDER THE BOHN TAX AT GIVEN BOND PRICES.
    delta_b       = cal[f"delta_b_{country}"]
    recovery_rate = cal.get(f"recovery_rate_{country}", 1.0)   # F never defaults
    phi_lamb      = cal[f"phi_lamb_{country}"]
    G             = cal[f"G_{country}"]
    b_gov_ss      = gs["b_gov_ss"]
    T             = len(Q_B_path)

    if def_real_path is None:
        def_real_path = np.zeros(T)
    if recap_path is None:
        recap_path = np.zeros(T)   # branch bailout outlays, financed by issuance

    if b_anchor is None:
        b_anchor = b_gov_ss
        Tax_base = gs["Tax_ss"]
    else:
        # branches re-anchor to the post-haircut stock, else the haircut becomes a
        # tax-cut windfall -> default expansionary, wrong-signed risk premium
        Tax_base = G + delta_b * b_anchor * (1.0 - gs["Q_B_ss"])

    b_gov_bop = np.empty(T)   # beginning-of-period stock
    b_gov_eop = np.empty(T)   # end-of-period stock (bank-held over t -> t+1)
    Tax       = np.empty(T)
    coupon    = np.empty(T)
    net_iss   = np.empty(T)

    b = float(b_gov_ss if b_gov0 is None else b_gov0)
    for t in range(T):
        b_gov_bop[t] = b
        surv_t    = 1.0 - def_real_path[t] * (1.0 - recovery_rate)
        # Bohn rule on the SURVIVING stock (pre-haircut would spike taxes ~31% GDP)
        Tax[t]    = Tax_base + phi_lamb * (b * surv_t - b_anchor)
        coupon[t] = delta_b * b * surv_t
        new_bonds = (G + recap_path[t] + coupon[t] - Tax[t]) / Q_B_path[t]
        net_iss[t] = Q_B_path[t] * new_bonds
        b = (1.0 - delta_b) * b * surv_t + new_bonds
        b_gov_eop[t] = b

    return dict(Tax=Tax, coupon=coupon, net_issuance=net_iss,
                b_gov=b_gov_bop, b_gov_eop=b_gov_eop)
