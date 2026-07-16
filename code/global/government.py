# **Government block: HM perpetuity bonds, Bohn (1998) fiscal rule, Cole-Kehoe crisis zones.**
import numpy as np


def hm_bond_price_ss(rdep_ss, delta_b):
    # **Steady-state HM perpetuity price.**
    return delta_b / (rdep_ss + delta_b)


def hm_bond_return_ss(Q_B_ss, delta_b):
    # **Realized SS return on the HM perpetuity (= rdep_ss by no-arbitrage).**
    return (delta_b + (1 - delta_b) * Q_B_ss) / Q_B_ss - 1


def govt_steady_state(cal, rdep_ss, country):
    # **Steady-state government block (no default, constant debt stock).**
    delta_b  = cal[f"delta_b_{country}"]
    B_gov_ss = cal[f"B_gov_{country}_ss"]
    G        = cal[f"G_{country}"]

    Q_B_ss    = hm_bond_price_ss(rdep_ss, delta_b)
    rb_ss     = hm_bond_return_ss(Q_B_ss, delta_b)
    coupon_ss = delta_b * B_gov_ss
    Tax_ss    = G + coupon_ss * (1.0 - Q_B_ss)   # budget: G + coupon = Tax + issuance
    return dict(Q_B_ss=Q_B_ss, rb_ss=rb_ss,
                Tax_ss=Tax_ss, b_gov_ss=B_gov_ss, coupon_ss=coupon_ss)


def ck_default_prob(b_gov, Y_ss, cal, sunspot, country):
    # **Cole-Kehoe zone map: 0 (safe) / sunspot (crisis zone) / 1 (certain default).**
    # sunspot = lenders' priced no-rollover probability in the crisis zone; never realized.
    b_low  = cal[f"b_ck_low_{country}"]
    b_high = cal[f"b_ck_high_{country}"]
    b_y    = b_gov / Y_ss
    return float(np.where(b_y < b_low, 0.0,
                          np.where(b_y >= b_high, 1.0, sunspot)))


def govt_transition(cal, gs, Q_B_path, def_real_path, country, b_gov0=None,
                    b_anchor=None, recap_path=None):
    # **Forward-integrate the debt stock under the Bohn tax at given bond prices.**
    # recap_path: default-branch bailout outlays (extra spending financed by issuance).
    delta_b       = cal[f"delta_b_{country}"]
    recovery_rate = cal[f"recovery_rate_{country}"]
    phi_lamb      = cal[f"phi_lamb_{country}"]
    G             = cal[f"G_{country}"]
    Tax_ss        = gs["Tax_ss"]
    b_gov_ss      = gs["b_gov_ss"]
    T             = len(Q_B_path)

    if def_real_path is None:
        def_real_path = np.zeros(T)
    if recap_path is None:
        recap_path = np.zeros(T)

    if b_anchor is None:
        b_anchor = b_gov_ss
        Tax_base = Tax_ss
    else:
        # branches re-anchor to the post-haircut stock (else the haircut becomes a
        # tax-cut windfall → default expansionary, wrong-signed risk premium)
        Tax_base = G + delta_b * b_anchor * (1.0 - gs["Q_B_ss"])

    b_gov_bop = np.empty(T)   # beginning-of-period stock
    b_gov_eop = np.empty(T)   # end-of-period stock (bank-held over t→t+1)
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
