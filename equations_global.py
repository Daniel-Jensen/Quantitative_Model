from sequence_jacobian import simple


@simple
def trade_balance(p, IM_D, IM_F):
    NX_D    = IM_F - p * IM_D
    NX_F    = IM_D - IM_F / p
    tob_res = NX_D
    return NX_D, NX_F, tob_res


@simple
def global_goods_mkt(goods_mkt_D, goods_mkt_F, p):
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
def bond_price(SDF_D, SDF_F, def_rate_D, def_rate_F, recovery_rate_D, recovery_rate_F):
    # Forward-looking asset pricing: bond is a one-period claim to face value 1,
    # reduced by haircut · default rate next period:
    #     q_b_t = E_t[SDF · (1 − haircut · def_rate(+1))]
    # Yield (promised) is then implied by the price:  1 + rb = 1 / q_b.
    # This inverts the prior convention (q_b derived from rb) so that q_b is the
    # primitive equilibrium asset price and rb is the YTM that drops out of it.
    haircut_D = 1.0 - recovery_rate_D
    haircut_F = 1.0 - recovery_rate_F
    q_b_D = SDF_D * (1 - haircut_D * def_rate_D(+1))
    q_b_F = SDF_F * (1 - haircut_F * def_rate_F(+1))
    rb_D  = 1.0 / q_b_D - 1.0
    rb_F  = 1.0 / q_b_F - 1.0
    return q_b_D, q_b_F, rb_D, rb_F


@simple
def portfolio_adj_cost(rb_actual_F, rb_actual_D, rdep_D, rdep_F,
                       b_F_D, b_D_F, n_inter_D, n_inter_F,
                       phi_bF_D_ss, phi_bD_F_ss, psi_bF_D, psi_bD_F,
                       excess_return_F_D_ss, excess_return_D_F_ss,
                       tau_mp_D, tau_mp_F):
    # Cross-border bond Euler equations (portfolio-adjustment-cost approach).
    # tau_mp_D is the per-unit macroprudential tax D-banks pay on bond holdings;
    # it enters the D-banks' FOC for F-bonds and raises the required excess return.
    # tau_mp_F enters F-banks' FOC for D-bonds symmetrically.

    phi_bF_D  = b_F_D / n_inter_D
    b_F_D_res = (rb_actual_F(+1) - rdep_D(+1)) - excess_return_F_D_ss \
                - psi_bF_D * (phi_bF_D - phi_bF_D_ss) \
                - tau_mp_D

    phi_bD_F  = b_D_F / n_inter_F
    b_D_F_res = (rb_actual_D(+1) - rdep_F(+1)) - excess_return_D_F_ss \
                - psi_bD_F * (phi_bD_F - phi_bD_F_ss) \
                - tau_mp_F

    return b_F_D_res, b_D_F_res


@simple
def terms_of_trade(p, pi_D, pi_F):
    p_res = p - p(-1) * (1 + pi_F) / (1 + pi_D)
    return p_res