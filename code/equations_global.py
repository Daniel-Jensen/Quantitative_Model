from sequence_jacobian import simple


# ── COUNTRY SIZE (2026-08-07) ────────────────────────────────────────────────
# `size_F` is F's size relative to D (Germany/Greece 2010 GDP = 11.697). Every
# F-side variable in the model is PER F CAPITA and stays O(1); D-side variables
# are aggregates (size_D == 1). The four blocks below are the only places the two
# countries meet, so they are the only places the weight appears.
#
# Why this exists: with both countries normalised to Y_ss = 1 the model could not
# match the EBA portfolio-composition moment (DE banks' Greek book / DE bank net
# worth = 0.0075) and the market-structure moment (foreigners hold 12.7% of the
# bank-held Greek stock) at the same time — joint consistency needs n_F/n_D = 8.85
# against the model's 0.761, a gap that is exactly the GDP ratio. Under per-capita
# F variables plus this weight, both hold simultaneously.
#
# Watch the units: b_D_F and IM_F and b_gov_F are PER F CAPITA, so they carry a
# `* size_F` when they meet a D aggregate; b_F_D and IM_D are D aggregates, so
# they carry a `/ size_F` when they meet an F per-capita quantity.


@simple
def trade_balance(p, IM_D, IM_F, size_F):
    # D's exports = what F imports, aggregated over F: size_F * IM_F (D goods).
    # D's imports = IM_D (F goods), valued at p.
    NX_D    = size_F * IM_F - p * IM_D
    # F's exports per F capita = IM_D / size_F (D goods), valued at 1/p in F goods
    # -- but IM_D is already in D goods, so the conversion is the /p on imports.
    NX_F    = IM_D / size_F - IM_F / p
    return NX_D, NX_F



@simple
def external_account_D(NX_D, q_b_D, q_b_F, b_F_D, b_D_F, rb_actual_F, rb_actual_D,
                       size_F):
    # b_F_D is a D aggregate; b_D_F is per F capita and must be aggregated.
    receipts_from_F_bonds = (1 + rb_actual_F) * q_b_F(-1) * b_F_D(-1)
    payments_on_D_bonds   = (1 + rb_actual_D) * q_b_D(-1) * size_F * b_D_F(-1)
    nfa_D = q_b_F * b_F_D - q_b_D * size_F * b_D_F
    ca_res_D = (NX_D + receipts_from_F_bonds - payments_on_D_bonds - nfa_D)
    return nfa_D, ca_res_D


@simple
def global_goods_mkt(goods_mkt_D, goods_mkt_F, p, size_F):
    # goods_mkt_F is per F capita; the union-wide residual needs it aggregated.
    global_goods_res = goods_mkt_D + p * size_F * goods_mkt_F
    return global_goods_res



@simple
def domestic_bond_clearing(b_gov_D, b_gov_F, b_D_F, b_F_D, size_F):
    # D debt (aggregate) is held by D banks plus size_F F-banks' per-capita books.
    b_D_D = b_gov_D - size_F * b_D_F
    # F debt is per F capita; the D aggregate holding is spread over size_F.
    b_F_F = b_gov_F - b_F_D / size_F
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
def terms_of_trade(p, pi_D, pi_F):
    # p = P_F/P_D in euro producer prices. In a monetary union the nominal
    # exchange rate is fixed at 1, so terms-of-trade movement IS the inflation
    # differential. This pins pi_D - pi_F off an unknown that already exists.
    # Zero at SS: p/p(-1) = 1 and pi_D = pi_F = 0.
    tot_res = p / p(-1) - (1.0 + pi_F) / (1.0 + pi_D)
    return tot_res


@simple
def union_inflation(pi_D, pi_F, omega_pi_D):
    # The ECB stabilises union-wide producer-price inflation -- the phi_pi -> inf
    # limit of a Taylor rule, stated as an abstraction and NOT a modelled rule.
    # Financial contracts carry no policy rate, so no Fisher relation is needed
    # to close the nominal side.
    #
    # With terms_of_trade this gives pi_D = -(1 - omega_pi_D)*dlog p. At the
    # capital-key omega_pi_D = 0.071, 93% of any terms-of-trade adjustment is D
    # producer-price deflation and 7% is F inflation -- the 2010-12 internal-
    # devaluation pattern. Do NOT use model GDP weights: the model normalises
    # Y_D_ss ~ Y_F_ss ~ 1, so those would give ~0.5 and split it evenly.
    pi_U         = omega_pi_D * pi_D + (1.0 - omega_pi_D) * pi_F
    union_pi_res = pi_U
    return pi_U, union_pi_res


@simple
def portfolio_level_anchors(b_F_D_anchor, b_D_F_anchor):
    b_F_D_ss = b_F_D_anchor
    b_D_F_ss = b_D_F_anchor
    return b_F_D_ss, b_D_F_ss


@simple
def gk_cross_border_foc(nu_bF_D, nu_K_D, Delta_bF_eff_D, SDF_banker_D, Omega_p1_D,
                        nu_bD_F, nu_K_F, Delta_bD_eff_F, SDF_banker_F, Omega_p1_F,
                        b_F_D, b_D_F, b_F_D_ss, b_D_F_ss, psi_bF_D, psi_bD_F,
                        gk_wedge_F_D_ss, gk_wedge_D_F_ss):
    """GK portfolio optimality on the two CROSS-BORDER sovereign legs.

    Replaces ``divert_portfolio_adj`` (and the never-wired ``portfolio_adj_cost``), which
    required

        (r_cross(+1) - rdep) = excess_return_*_ss
                             + (EL_price_issuer + psi_spread_holder) * def_rate_issuer(+1)
                             - psi * (b - b_ss) - tau_mp.

    THE DOUBLE-COUNT THAT IS NOW GONE. From 2026-08-17 the expected loss also sat inside
    ``intermediation_P1_D/F``, so ``b_F_D`` and ``b_D_F`` had ``EL * def_rate`` netted off
    TWICE — once in the ``nu`` that values the position and again in the block that chose
    it — while the two own-sovereign legs netted it once. Both cross-border legs now read
    the same ``nu``s the own legs do, so the loss enters exactly once on all four legs.
    ``psi_spread_D/F`` and the ``excess_return_*_ss`` anchors are deleted outright.

    THE CONDITION. Same first-order condition as ``gk_bond_foc_D``, plus a quadratic
    portfolio adjustment cost on the cross-border stock:

        nu_bF_D / nu_K_D = Delta_bF_eff_D  +  (adjustment cost).

    Why the cross legs need the cost and the own legs do not: under a binding IC with
    linear payoffs the banker's problem is linear in portfolio shares, so four
    proportionality conditions cannot hold simultaneously against only two bond prices.
    The own-sovereign legs pin ``q_b_D`` and ``q_b_F``; the cross-border legs then pin
    QUANTITIES, with ``psi_bF_D``/``psi_bD_F`` as the frictions that make an interior
    cross-border position optimal at all. This is the standard portfolio-cost device and
    it is NOT a sovereign-spread wedge: it loads on the bond STOCK gap, carries no
    ``def_rate`` term, and is identically zero at the calibrated position.

    UNITS. The residual is divided through by ``SDF_banker * Omega`` so it is stated in
    RETURN units, which is what ``psi_bF_D``/``psi_bD_F`` were calibrated against; the
    division is exact and changes no root (both factors are strictly positive — see
    ``steady_state.assert_gk_well_posed``). Ordering the terms this way also keeps
    ``psi_bF_D`` comparable with its pre-refactor value.

    ``gk_wedge_*_ss`` is the CONSTANT steady-state shadow cost of holding the
    EBA-measured cross-border position — the level at which the adjustment cost sits at
    the calibrated stock. It is a scalar set once in ``_apply_ss_anchors``, it does not
    move with ``def_rate``, and at the preferred calibration it is ~0 because GK
    optimality at a riskless symmetric steady state forces the cross divertability to
    equal the own one. ``steady_state`` prints it and asserts it stays small precisely so
    an inconsistent ``Delta`` cannot hide inside it.
    """
    # D bank holding F bonds. No terms-of-trade conversion: q_b_F is a D-GOOD price
    # (intermediation_P3_D, external_account_D), so nu_bF_D is already in D goods.
    foc_F_D   = (nu_bF_D - Delta_bF_eff_D * nu_K_D) / (SDF_banker_D * Omega_p1_D)
    b_F_D_res = foc_F_D - gk_wedge_F_D_ss - psi_bF_D * (b_F_D - b_F_D_ss)

    # F bank holding D bonds. The p/p(+1) conversion into F goods already happened
    # inside intermediation_P1_F, so nu_bD_F and nu_K_F are both F-good objects here.
    foc_D_F   = (nu_bD_F - Delta_bD_eff_F * nu_K_F) / (SDF_banker_F * Omega_p1_F)
    b_D_F_res = foc_D_F - gk_wedge_D_F_ss - psi_bD_F * (b_D_F - b_D_F_ss)

    return b_F_D_res, b_D_F_res
