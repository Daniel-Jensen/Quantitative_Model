# GOVERNMENT BLOCK: HATCHONDO-MARTINEZ PERPETUITY BONDS, BOHN (1998) FISCAL RULE.
# Default risk is EXOGENOUS (Bocola 2016 eqs. 11-12): the priced default
# probability is an input path to the transition solver, never a function of
# the debt stock. Debt still evolves endogenously under the Bohn tax.


def govt_steady_state(cal, country, Q_B_ss):
    # STEADY-STATE GOVERNMENT BLOCK (NO DEFAULT, CONSTANT DEBT STOCK) AT THE SS BOND PRICE.
    delta_b  = cal[f"delta_b_{country}"]
    B_gov_ss = cal[f"B_gov_{country}_ss"]
    G        = cal[f"G_{country}"]

    Tax_ss = G + delta_b * B_gov_ss * (1.0 - Q_B_ss)   # G + coupon = Tax + issuance
    # the Bohn coefficient is the calibrated primitive (2026-10-26, appendix B): the annual response
    # of the primary surplus to debt, both over annual output, is 4*gamma with quarterly Y = 1. With
    # B' = (1-delta_b)*x + (G + delta_b*x - Tax)/Q and Tax = Tax_ss + gamma*(x - B_ss) on the surviving
    # stock x = B*surv, the debt root it implies is dB'/dx = (1-delta_b) + (delta_b-gamma)/Q_B_ss,
    # reported, not targeted (the root was the primitive until 2026-10-26)
    gamma_tau = cal["bohn_annual"] / 4.0
    return dict(Q_B_ss=Q_B_ss, Tax_ss=Tax_ss, b_gov_ss=B_gov_ss, gamma_tau=gamma_tau)
