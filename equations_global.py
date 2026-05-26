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
                       excess_return_F_D_ss, excess_return_D_F_ss,
                       mp_wedge_F, mp_wedge_D):
    # Cross-border bond Euler equations.
    #
    # mp_wedge_F enters the D-banks' FOC for F-bonds: when the regulator
    # raises the capital charge on F-sovereign exposure, D-banks require a
    # higher excess return to hold F-bonds.  This forces b_F_D down (D-banks
    # reduce F-bond holdings) and rb_F up — the portfolio-rebalancing channel.
    #
    # mp_wedge_D enters the F-banks' FOC for D-bonds symmetrically: when D's
    # risk weight rises, F-banks flee D-sovereign bonds (b_D_F falls), which
    # is the capital-flight / market-fragmentation channel.
    #
    # Without these wedge terms the macroprudential policy only affects the
    # within-country domestic bond FOCs (domestic_bond_foc_D/F) and the
    # cross-border portfolio is entirely unresponsive to the regulatory charge.

    phi_bF_D  = b_F_D / n_inter_D          # no p: same currency
    b_F_D_res = (rb_actual_F(+1) - rdep_D(+1)) - excess_return_F_D_ss \
                - psi_bF_D * (phi_bF_D - phi_bF_D_ss) \
                - mp_wedge_F                # wedge on F-sovereign, borne by D-banks

    phi_bD_F  = b_D_F / n_inter_F          # no 1/p: same currency
    b_D_F_res = (rb_actual_D(+1) - rdep_F(+1)) - excess_return_D_F_ss \
                - psi_bD_F * (phi_bD_F - phi_bD_F_ss) \
                - mp_wedge_D                # wedge on D-sovereign, borne by F-banks

    return b_F_D_res, b_D_F_res


@simple
def terms_of_trade(p, pi_D, pi_F):
    p_res = p - p(-1) * (1 + pi_F) / (1 + pi_D)
    return p_res


@simple
def risk_weight(def_rate_D, def_rate_F, alpha_w_D, alpha_w_F,
                shock_w_sov_D, shock_w_sov_F):
    # Sovereign risk weights respond CONTEMPORANEOUSLY to the realised default
    # rate at time t, plus an exogenous credit-rating shock (shock_w_sov_c).
    #
    # Economic timing (all at period t):
    #   shock_def_D hits → def_rate_D rises → w_sov_D rises → mp_wedge_D rises
    #   → domestic_bond_foc_D forces rb_D up → government faces higher borrowing
    #   cost going forward.  This is the "macroprudential surcharge" channel.
    #
    # No contemporaneous DAG cycle: def_rate_D(t) depends only on
    # b_gov_D(t-1)/Y_D(t-1) (lagged) + shock_def_D(t) (exogenous), so
    # risk_weight → macroprudential_wedge is a valid forward edge in the block
    # DAG with no feedback to def_rate_D at time t.
    #
    # At SS def_rate_c = 0 and shock_w_sov_c = 0 → w_sov_c = 1 (neutral).
    w_sov_D = 1.0 + alpha_w_D * def_rate_D + shock_w_sov_D
    w_sov_F = 1.0 + alpha_w_F * def_rate_F + shock_w_sov_F
    return w_sov_D, w_sov_F


@simple
def macroprudential_wedge(b_gov_D, b_gov_F, n_inter_D, n_inter_F,
                          vartheta_D, vartheta_F, w_sov_D, w_sov_F,
                          phi_sov_D_ss, phi_sov_F_ss):
    n_total    = n_inter_D + n_inter_F
    phi_sov_D  = b_gov_D / n_total
    phi_sov_F  = b_gov_F / n_total
    mp_wedge_D = vartheta_D * (w_sov_D * phi_sov_D - phi_sov_D_ss)
    mp_wedge_F = vartheta_F * (w_sov_F * phi_sov_F - phi_sov_F_ss)
    return phi_sov_D, phi_sov_F, mp_wedge_D, mp_wedge_F