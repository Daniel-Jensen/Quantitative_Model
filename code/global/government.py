"Government block: Hatchondo-Martinez (2009) geometric-decay perpetuity bonds, Bohn (1998) fiscal rule, Cole-Kehoe (2000) self-fulfilling crisis zones."


import numpy as np


# ── Steady-state helpers ──────────────────────────────────────────────────────

def hm_bond_price_ss(rdep_ss, delta_b):
    # Price of bond in steady state (it pays out coupon delta)
    return delta_b / (rdep_ss + delta_b)


def hm_bond_return_ss(Q_B_ss, delta_b):
    #Realised return at SS on HM perpetuity (= rdep_ss by no-arbitrage).
    return (delta_b + (1 - delta_b) * Q_B_ss) / Q_B_ss - 1


def govt_steady_state(cal, rdep_ss, country):
    # Steady-state government block. Zero risk of default in SS, and bond stock constant

    delta_b  = cal[f"delta_b_{country}"]
    B_gov_ss = cal[f"B_gov_{country}_ss"]
    G        = cal[f"G_{country}"]

    Q_B_ss    = hm_bond_price_ss(rdep_ss, delta_b)
    rb_ss     = hm_bond_return_ss(Q_B_ss, delta_b)
    coupon_ss = delta_b * B_gov_ss
    
    # Government budget at SS: G + coupon = Tax + issuance_proceeds
    # Issuance = delta_b*B_gov new bonds at price Q_B_ss
    Tax_ss    = G + coupon_ss * (1.0 - Q_B_ss)
    return dict(Q_B_ss=Q_B_ss, rb_ss=rb_ss,
                Tax_ss=Tax_ss, b_gov_ss=B_gov_ss, coupon_ss=coupon_ss)


# ── Cole-Kehoe default probability ───────────────────────────────────────────

def ck_default_prob(b_gov, Y_ss, cal, sunspot, country):
    """Cole-Kehoe (2000) default probability for a single period.

    Three zones based on debt-to-output ratio b_gov/Y_ss:
      Safe zone   (b/Y < b_ck_low):           returns 0.
      Crisis zone (b_ck_low ≤ b/Y < b_ck_high): returns sunspot ∈ [0,1].
      Certain-default (b/Y ≥ b_ck_high):      returns 1.

    `sunspot` is the exogenous probability that lenders coordinate on the
    no-rollover equilibrium, conditional on the crisis zone being active —
    the analogue of Bocola (2016)'s exogenous AR(1) default-risk process s_t
    (his eq. 12), restricted to the CK crisis zone.  In the risk-only
    experiment this probability is PRICED but default is never REALIZED.
    """
    b_low  = cal[f"b_ck_low_{country}"]
    b_high = cal[f"b_ck_high_{country}"]
    b_y    = b_gov / Y_ss
    return float(np.where(b_y < b_low, 0.0,
                          np.where(b_y >= b_high, 1.0, sunspot)))


# ── Transition-path government block ─────────────────────────────────────────

def govt_transition(cal, gs, Q_B_path, def_real_path, country, b_gov0=None,
                    b_anchor=None):
  # CALCULATES TAX ADJUSTMENT TAKING PRICES AS GIVEN (NO FEEDBACK FROM TAXES TO PRICES), AND THE ENTIRE SPENDING PATH

    # unpack parameters
    delta_b       = cal[f"delta_b_{country}"]
    recovery_rate = cal[f"recovery_rate_{country}"]
    phi_lamb      = cal[f"phi_lamb_{country}"]
    G             = cal[f"G_{country}"]
    Tax_ss        = gs["Tax_ss"]
    b_gov_ss      = gs["b_gov_ss"]
    T             = len(Q_B_path)

    # initialize arrays
    if def_real_path is None:
        def_real_path = np.zeros(T)

    if b_anchor is None:
        b_anchor = b_gov_ss
        Tax_base = Tax_ss
    else:
        # Budget-balancing tax at the anchor (stationary at b = b_anchor)
        Tax_base = cal[f"G_{country}"] + delta_b * b_anchor * (1.0 - gs["Q_B_ss"])

    b_gov_bop = np.empty(T)   # stock at beginning of period t
    b_gov_eop = np.empty(T)   # stock at end of period t (held by banks over t→t+1)
    Tax       = np.empty(T)
    coupon    = np.empty(T)
    net_iss   = np.empty(T)

    # steady state value of bonds
    b = float(b_gov_ss if b_gov0 is None else b_gov0)
    for t in range(T):
        b_gov_bop[t] = b
        #tax rule
        Tax[t]    = Tax_base + phi_lamb * (b - b_anchor)

        # Coupon payments calculation
        surv_t    = 1.0 - def_real_path[t] * (1.0 - recovery_rate)
        coupon[t] = delta_b * b * surv_t

        #new bonds
        new_bonds = (G + coupon[t] - Tax[t]) / Q_B_path[t]
        net_iss[t] = Q_B_path[t] * new_bonds
        b = (1.0 - delta_b) * b * surv_t + new_bonds
        b_gov_eop[t] = b

    return dict(Tax=Tax, coupon=coupon, net_issuance=net_iss,
                b_gov=b_gov_bop, b_gov_eop=b_gov_eop)
