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


@simple
def risk_weight(def_rate_D, def_rate_F, alpha_w_D, alpha_w_F,
                shock_w_sov_D, shock_w_sov_F):
    # Sovereign risk weights respond to lagged realised default rate plus an
    # exogenous "credit rating shock".  Rating agencies act on observed
    # default — the lag mimics this convention and also breaks the otherwise
    # simultaneous w_sov ↔ mp_wedge ↔ rb ↔ b_gov ↔ def_rate loop.
    #
    # At SS def_rate_c = 0 and shock_w_sov_c = 0 → w_sov_c = 1 (neutral).
    # alpha_w_c is the rating sensitivity to default; 0 recovers the
    # exogenous-w_sov regime of PRs #6/#7.  shock_w_sov_c is the new
    # exogenous driver — perturbing it simulates a sovereign credit-rating
    # action (Moody's downgrade, ECB collateral haircut announcement, etc.).
    w_sov_D = 1.0 + alpha_w_D * def_rate_D(-1) + shock_w_sov_D
    w_sov_F = 1.0 + alpha_w_F * def_rate_F(-1) + shock_w_sov_F
    return w_sov_D, w_sov_F


@simple
def macroprudential_wedge(b_gov_D, b_gov_F, n_inter_D, n_inter_F,
                          vartheta_D, vartheta_F, w_sov_D, w_sov_F,
                          phi_sov_D_ss, phi_sov_F_ss):
    # Sovereign-capital-charge wedge.  Each sovereign carries ONE wedge that
    # enters every bond FOC for that sovereign — same number for all holders.
    # phi_sov_c is global banking-sector exposure to sovereign c (consolidated
    # issuance / consolidated bank equity).  By bond market clearing
    #     b_gov_D = b_D_D + b_D_F   and   b_gov_F = b_F_F + b_F_D,
    # so this is identical to (b_c_c + b_c_~c) / (n_D + n_F).  We write it in
    # terms of b_gov_c to avoid a degenerate redundancy: in ha_full the block
    # domestic_bond_clearing pins b_D_D = b_gov_D - b_D_F, so taking both
    # b_D_D and b_D_F as inputs would give sequence_jacobian an exact zero
    # composed Jacobian that it stores as an empty SimpleSparse and later
    # fails to compose with.
    #
    # Functional form: mp_wedge_c = vartheta_c · (w_sov_c · phi_sov_c - phi_sov_c_ss).
    # The bracket is the deviation of RISK-WEIGHTED exposure from its SS
    # anchor — Basel-standard framing.  At SS w_sov_c = 1 and
    # phi_sov_c = phi_sov_c_ss, so the bracket is zero by construction.
    # When w_sov_c = 1 this collapses to vartheta_c · (phi_sov_c - phi_sov_c_ss)
    # — identical to the original mp_wedge formula and bit-equal to the
    # PR #6 / #7 IRFs at any vartheta.  When w_sov_c > 1 (rating downgrade),
    # the wedge bites even at unchanged nominal exposure: the regulator
    # treats the same b_gov as a larger "risk-weighted" position.  This
    # gives w_sov_c a nonzero first-order Jacobian in the linearised system
    # (it had none in the old multiplicative form because the linearisation
    # evaluated at phi_sov = phi_sov_ss made the d/dw_sov derivative zero).
    n_total    = n_inter_D + n_inter_F
    phi_sov_D  = b_gov_D / n_total
    phi_sov_F  = b_gov_F / n_total
    mp_wedge_D = vartheta_D * (w_sov_D * phi_sov_D - phi_sov_D_ss)
    mp_wedge_F = vartheta_F * (w_sov_F * phi_sov_F - phi_sov_F_ss)
    return phi_sov_D, phi_sov_F, mp_wedge_D, mp_wedge_F