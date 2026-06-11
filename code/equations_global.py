from sequence_jacobian import simple


@simple
def trade_balance(p, IM_D, IM_F):
    NX_D    = IM_F - p * IM_D
    NX_F    = IM_D - IM_F / p
    return NX_D, NX_F



@simple
def external_account_D(NX_D, q_b_D, q_b_F, b_F_D, b_D_F,rb_actual_F, rb_actual_D):
    receipts_from_F_bonds = (1 + rb_actual_F) * q_b_F(-1) * b_F_D(-1)
    payments_on_D_bonds   = (1 + rb_actual_D) * q_b_D(-1) * b_D_F(-1)
    nfa_D = q_b_F * b_F_D - q_b_D * b_D_F
    ca_res_D = (NX_D + receipts_from_F_bonds - payments_on_D_bonds- nfa_D)
    return nfa_D, ca_res_D


@simple
def global_goods_mkt(goods_mkt_D, goods_mkt_F, p):
    global_goods_res = goods_mkt_D + p * goods_mkt_F
    return global_goods_res



@simple
def domestic_bond_clearing(b_gov_D, b_gov_F, b_D_F, b_F_D):
    b_D_D = b_gov_D - b_D_F
    b_F_F = b_gov_F - b_F_D
    return b_D_D, b_F_F


@simple
def bond_yield(q_b_D, q_b_F, delta_b_D, delta_b_F):
    # Woodford perpetuity holding-period return: rb = delta_b * (1/q_b - 1)
    # This equals rb_actual in SS and gives the correct annualised yield.
    # The old formula 1/q_b - 1 treated q_b as a zero-coupon price and
    # overstated the yield by a factor of 1/delta_b (~20×).
    rb_D      = delta_b_D * (1.0 / q_b_D - 1.0)
    rb_F      = delta_b_F * (1.0 / q_b_F - 1.0)
    spread_rb = rb_D - rb_F
    return rb_D, rb_F, spread_rb


@simple
def portfolio_level_anchors(b_F_D_anchor, b_D_F_anchor):
    b_F_D_ss = b_F_D_anchor
    b_D_F_ss = b_D_F_anchor
    return b_F_D_ss, b_D_F_ss


@simple
def portfolio_adj_cost(rb_actual_F, rb_actual_D, rdep_D, rdep_F,
                       b_F_D, b_D_F,
                       b_F_D_ss, b_D_F_ss,
                       psi_bF_D, psi_bD_F,
                       excess_return_F_D_ss, excess_return_D_F_ss,
                       tau_mp_D, tau_mp_F, p):
    # Level penalty on face-value bond stocks anchors the external position level,
    # not only its composition relative to net worth.
    # Expected D-good return on F-bonds: (1+rb_F)·p(+1)/p − 1
    rb_F_dg_next = (1 + rb_actual_F(+1)) * p(+1) / p - 1
    b_F_D_res = (rb_F_dg_next - rdep_D(+1)) - excess_return_F_D_ss \
                - psi_bF_D * (b_F_D - b_F_D_ss) \
                - tau_mp_D

    # Expected F-good return on D-bonds: (1+rb_D)·p/p(+1) − 1
    rb_D_fg_next = (1 + rb_actual_D(+1)) * p / p(+1) - 1
    b_D_F_res    = (rb_D_fg_next - rdep_F(+1)) - excess_return_D_F_ss \
                   - psi_bD_F * (b_D_F - b_D_F_ss) \
                   - tau_mp_F

    return b_F_D_res, b_D_F_res


@simple
def divert_portfolio_adj(rb_actual_F, rb_actual_D, rdep_D, rdep_F, p,
                         b_F_D, b_D_F, b_F_D_ss, b_D_F_ss, psi_bF_D, psi_bD_F,
                         excess_return_F_D_ss, excess_return_D_F_ss, tau_mp_D, tau_mp_F,
                         psi_spread_D, psi_spread_F, def_rate_D, def_rate_F):
    # D holds F-bonds (F-good claim -> convert with p); issuer = F
    rb_F_dg_next = (1 + rb_actual_F(+1)) * p(+1) / p - 1
    # IC-theory derived required premium: D-bank IC parameters govern D-bank's FOC on F-bonds
    prem_FD      = excess_return_F_D_ss + psi_spread_D * def_rate_F(+1)
    # T-2 fix: deposit rate for the t->t+1 holding period is locked at t (rdep, not rdep(+1)).
    b_F_D_res    = (rb_F_dg_next - rdep_D) - prem_FD \
                   - psi_bF_D * (b_F_D - b_F_D_ss) - tau_mp_D
    # F holds D-bonds (D-good claim -> convert with p); issuer = D
    rb_D_fg_next = (1 + rb_actual_D(+1)) * p / p(+1) - 1
    # IC-theory derived required premium: F-bank IC parameters govern F-bank's FOC on D-bonds
    prem_DF      = excess_return_D_F_ss + psi_spread_F * def_rate_D(+1)
    b_D_F_res    = (rb_D_fg_next - rdep_F) - prem_DF \
                   - psi_bD_F * (b_D_F - b_D_F_ss) - tau_mp_F
    return b_F_D_res, b_D_F_res
