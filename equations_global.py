from sequence_jacobian import simple


@simple
def trade_balance(p, IM_D, IM_F,
                  b_D_F, b_F_D, q_b_D, q_b_F,
                  def_rate_D, def_rate_F, recovery_rate_D, recovery_rate_F):
    # NX_X = X's net exports (in X-units) = X's goods out − X's goods in.
    # Cross-border bond payments are PHYSICAL resource flows that must enter NX:
    #   D→F (D-units): D-gov interest (1−def_D·h)·b_D_F(−1) minus F-bank new lending q_b_D·b_D_F
    #   F→D (F-units): symmetric
    # The augmented NX captures both the goods trade balance and the financial
    # flows; without this, cross-border bond payments leak from Walras.
    haircut_D = 1.0 - recovery_rate_D
    haircut_F = 1.0 - recovery_rate_F
    net_D_out_via_bonds = (1.0 - haircut_D * def_rate_D) * b_D_F(-1) - q_b_D * b_D_F
    net_F_out_via_bonds = (1.0 - haircut_F * def_rate_F) * b_F_D(-1) - q_b_F * b_F_D
    NX_D    = IM_F - p * IM_D + net_D_out_via_bonds - p * net_F_out_via_bonds
    NX_F    = IM_D - IM_F / p + net_F_out_via_bonds - net_D_out_via_bonds / p
    tob_res = NX_D + p * NX_F        # identity check, should be machine-precision 0
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
def bond_yield(q_b_D, q_b_F):
    # Implied yields from bond prices — for interpretation only.
    # q_b is the primary model variable; rb is derived here and NOT used elsewhere.
    rb_D = 1.0 / q_b_D - 1.0
    rb_F = 1.0 / q_b_F - 1.0
    return rb_D, rb_F


@simple
def portfolio_adj_cost(rb_actual_F, rb_actual_D, rdep_D, rdep_F,
                       b_F_D, b_D_F, n_inter_D, n_inter_F,
                       q_b_F, q_b_D,
                       phi_bF_D_ss, phi_bD_F_ss, psi_bF_D, psi_bD_F,
                       excess_return_F_D_ss, excess_return_D_F_ss,
                       tau_mp_D, tau_mp_F):
    # Cross-border bond Euler equations (portfolio-adjustment-cost approach).
    # phi shares use market values (q_b * face_value / net_worth) for consistency
    # with the market-value balance sheet in k_balance_sheet_D/F.
    # psi pins the equilibrium quantity (slope of demand curve).
    # tau_mp / T0 shifts the required excess return (calibrates the price spread).

    phi_bF_D  = q_b_F * b_F_D / n_inter_D
    b_F_D_res = (rb_actual_F(+1) - rdep_D(+1)) - excess_return_F_D_ss \
                - psi_bF_D * (phi_bF_D - phi_bF_D_ss) \
                - tau_mp_D

    phi_bD_F  = q_b_D * b_D_F / n_inter_F
    b_D_F_res = (rb_actual_D(+1) - rdep_F(+1)) - excess_return_D_F_ss \
                - psi_bD_F * (phi_bD_F - phi_bD_F_ss) \
                - tau_mp_F

    return b_F_D_res, b_D_F_res


@simple
def terms_of_trade(p, pi_D, pi_F):
    p_res = p - p(-1) * (1 + pi_F) / (1 + pi_D)
    return p_res