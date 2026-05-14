from sequence_jacobian import simple


@simple
def trade_balance(p, IM_D, IM_F):
    # IM_D: D's imports of F goods, in F-units.  Cost in D-units = p * IM_D.
    # IM_F: F's imports of D goods, in D-units.  Cost in F-units = IM_F / p.
    NX_D    = IM_F - p * IM_D
    NX_F    = IM_D - IM_F / p
    # tob_res is kept as NX_D for legacy output; it is NOT used as a target.
    # The exchange rate p is now pinned by goods_mkt_D = 0 (see targets_ss / targets_tp).
    tob_res = NX_D
    return NX_D, NX_F, tob_res


@simple
def global_goods_mkt(goods_mkt_D, goods_mkt_F, p):
    # Walras identity: goods_mkt_D + p * goods_mkt_F = 0 exactly when all
    # budget constraints hold in consistent units.  Non-zero residual signals
    # a genuine unit inconsistency (e.g. P_CES ≠ 1 mixing bundle/D-good units).
    global_goods_res = goods_mkt_D + p * goods_mkt_F
    return global_goods_res



@simple
def global_bond_market(B_supply_D, B_supply_F, b_D_D, b_D_F, b_F_F, b_F_D):
    bond_mkt_D_res = B_supply_D - (b_D_D + b_D_F)
    bond_mkt_F_res = B_supply_F - (b_F_F + b_F_D)
    return bond_mkt_D_res, bond_mkt_F_res


@simple
def domestic_bond_clearing(b_gov_D, b_gov_F, b_D_F, b_F_D):
    b_D_D = b_gov_D - b_D_F
    b_F_F = b_gov_F - b_F_D
    return b_D_D, b_F_F



@simple
def portfolio_adj_cost(rb_actual_F, rb_actual_D, rdep_D, rdep_F,
                       b_F_D, b_D_F, n_inter_D, n_inter_F,
                       phi_bF_D_ss, phi_bD_F_ss, psi_bF_D, psi_bD_F,
                       excess_return_F_D_ss, excess_return_D_F_ss):

    phi_bF_D  = b_F_D / n_inter_D          # no p: same currency
    b_F_D_res = (rb_actual_F(+1) - rdep_D(+1)) - excess_return_F_D_ss \
                - psi_bF_D * (phi_bF_D - phi_bF_D_ss)

    phi_bD_F  = b_D_F / n_inter_F          # no 1/p: same currency
    b_D_F_res = (rb_actual_D(+1) - rdep_F(+1)) - excess_return_D_F_ss \
                - psi_bD_F * (phi_bD_F - phi_bD_F_ss)

    return b_F_D_res, b_D_F_res


@simple
def terms_of_trade(p, pi_D, pi_F):
    p_res = p - p(-1) * (1 + pi_F) / (1 + pi_D)
    return p_res