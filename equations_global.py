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
def bond_yield(q_b_D, q_b_F):
    # Implied yields from bond prices — for interpretation only.
    # q_b is the primary model variable; rb is derived here and NOT used elsewhere.
    rb_D      = 1.0 / q_b_D - 1.0
    rb_F      = 1.0 / q_b_F - 1.0
    spread_rb = rb_D - rb_F
    return rb_D, rb_F, spread_rb


@simple
def portfolio_level_anchors(b_F_D_anchor, b_D_F_anchor):
    """Rename fixed calibration scalars so SJ's block graph tracks them as outputs."""
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
                         lambda_BF_D, lambda_BD_F, psi_lambda_B_D, psi_lambda_B_F,
                         def_rate_D, def_rate_F):
    # D holds F-bonds (F-good claim -> convert with p); issuer = F
    rb_F_dg_next = (1 + rb_actual_F(+1)) * p(+1) / p - 1
    lam_FD_eff   = lambda_BF_D + psi_lambda_B_D * def_rate_F(+1)
    prem_FD      = (lam_FD_eff / lambda_BF_D) * excess_return_F_D_ss
    b_F_D_res    = (rb_F_dg_next - rdep_D(+1)) - prem_FD \
                   - psi_bF_D * (b_F_D - b_F_D_ss) - tau_mp_D
    # F holds D-bonds (D-good claim -> convert with p); issuer = D
    rb_D_fg_next = (1 + rb_actual_D(+1)) * p / p(+1) - 1
    lam_DF_eff   = lambda_BD_F + psi_lambda_B_F * def_rate_D(+1)
    prem_DF      = (lam_DF_eff / lambda_BD_F) * excess_return_D_F_ss
    b_D_F_res    = (rb_D_fg_next - rdep_F(+1)) - prem_DF \
                   - psi_bD_F * (b_D_F - b_D_F_ss) - tau_mp_F
    return b_F_D_res, b_D_F_res
