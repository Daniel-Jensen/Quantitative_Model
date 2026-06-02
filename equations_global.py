from sequence_jacobian import simple


@simple
def trade_balance(p, IM_D, IM_F,
                  b_D_F, b_F_D, q_b_D, q_b_F,
                  def_rate_D, def_rate_F, recovery_rate_D, recovery_rate_F):
    # NX_X = X's net exports = X's goods leaving X − X's goods entering X (in X-units).
    # Cross-border bond payments are PHYSICAL resource flows that must enter NX:
    #   D→F (D-units): D-gov pays (1−def_D·h)·b_D_F(-1) at t, F-bank pays q_b_D·b_D_F for new bonds.
    #   F→D (F-units): symmetric.
    haircut_D = 1.0 - recovery_rate_D
    haircut_F = 1.0 - recovery_rate_F
    net_D_out_via_bonds = (1.0 - haircut_D * def_rate_D) * b_D_F(-1) - q_b_D * b_D_F
    net_F_out_via_bonds = (1.0 - haircut_F * def_rate_F) * b_F_D(-1) - q_b_F * b_F_D
    # Augmented bilateral trade balance (goods + cross-border bond flows).
    NX_D    = IM_F - p * IM_D + net_D_out_via_bonds - p * net_F_out_via_bonds
    NX_F    = IM_D - IM_F / p + net_F_out_via_bonds - net_D_out_via_bonds / p
    tob_res = NX_D + p * NX_F        # identity: must be 0
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
def cross_border_returns(rb_actual_D, rb_actual_F, p):
    # D-bonds are face-value in D-units; F-bonds face-value in F-units.
    # p = price of F-good in D-units.  Cross-border holdings must be expressed
    # in the holder's currency:
    #   F-bond payoff in D-units per F-unit invested → multiply by p/p(-1)
    #   D-bond payoff in F-units per D-unit invested → multiply by p(-1)/p
    # These convert REALIZED period-t returns of foreign bonds into the
    # holder-country's accounting unit.
    rb_F_in_D = (1.0 + rb_actual_F) * p      / p(-1) - 1.0
    rb_D_in_F = (1.0 + rb_actual_D) * p(-1)  / p     - 1.0
    return rb_F_in_D, rb_D_in_F


@simple
def portfolio_adj_cost(rb_actual_F, rb_actual_D, rdep_D, rdep_F,
                       b_F_D, b_D_F, n_inter_D, n_inter_F,
                       q_b_F, q_b_D, p,
                       phi_bF_D_ss, phi_bD_F_ss, psi_bF_D, psi_bD_F,
                       excess_return_F_D_ss, excess_return_D_F_ss,
                       tau_mp_D, tau_mp_F):
    # Cross-border bond Euler equations.  Shares and expected returns are
    # expressed in the HOLDING country's accounting unit, requiring p conversion.
    #   D bank holds F-bonds:  market value q_b_F·b_F_D in F-units → ·p in D-units.
    #   F bank holds D-bonds:  market value q_b_D·b_D_F in D-units → /p in F-units.

    phi_bF_D       = q_b_F * b_F_D * p / n_inter_D
    rb_F_in_D_next = (1.0 + rb_actual_F(+1)) * p(+1) / p - 1.0
    b_F_D_res      = (rb_F_in_D_next - rdep_D(+1)) - excess_return_F_D_ss \
                     - psi_bF_D * (phi_bF_D - phi_bF_D_ss) \
                     - tau_mp_D

    phi_bD_F       = q_b_D * b_D_F / (n_inter_F * p)
    rb_D_in_F_next = (1.0 + rb_actual_D(+1)) * p / p(+1) - 1.0
    b_D_F_res      = (rb_D_in_F_next - rdep_F(+1)) - excess_return_D_F_ss \
                     - psi_bD_F * (phi_bD_F - phi_bD_F_ss) \
                     - tau_mp_F

    return b_F_D_res, b_D_F_res


@simple
def terms_of_trade(p, pi_D, pi_F):
    p_res = p - p(-1) * (1 + pi_F) / (1 + pi_D)
    return p_res