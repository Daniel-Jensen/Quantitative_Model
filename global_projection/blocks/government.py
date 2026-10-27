# GOVERNMENT BLOCK: HATCHONDO-MARTINEZ PERPETUITY BONDS, BOHN (1998) FISCAL RULE.
# Default risk is EXOGENOUS (Bocola 2016 eqs. 11-12): the priced default
# probability is an input path to the transition solver, never a function of
# the debt stock. Debt still evolves endogenously under the Bohn tax.


def govt_steady_state(cal, rdep_ss, country):
    # STEADY-STATE GOVERNMENT BLOCK (NO DEFAULT, CONSTANT DEBT STOCK).
    delta_b  = cal[f"delta_b_{country}"]
    B_gov_ss = cal[f"B_gov_{country}_ss"]
    G        = cal[f"G_{country}"]

    Q_B_ss = delta_b / (rdep_ss + delta_b)
    Tax_ss = G + delta_b * B_gov_ss * (1.0 - Q_B_ss)   # G + coupon = Tax + issuance
    # THE QUARTERLY TAX RESPONSE gamma_tau. The paper sets the annual Bohn coefficient
    # (surplus and debt both over annual GDP), which is 4*gamma_tau. The Bocola calibration
    # instead targets a debt root: with B' = (1-delta_b)*x + (G + delta_b*x - Tax)/Q and
    # Tax = Tax_ss + gamma*(x - B_ss) on the surviving stock x = B*surv, the debt root is
    # dB'/dx = (1-delta_b) + (delta_b-gamma)/Q_B_ss, inverted here for gamma.
    if f"bohn_{country}" in cal:
        gamma_tau = cal[f"bohn_{country}"] / 4.0
    else:
        gamma_tau = delta_b - (cal[f"debt_root_{country}"] - (1.0 - delta_b)) * Q_B_ss
    return dict(Q_B_ss=Q_B_ss, Tax_ss=Tax_ss, b_gov_ss=B_gov_ss, gamma_tau=gamma_tau)


