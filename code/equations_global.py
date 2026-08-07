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
                         psi_spread_D, psi_spread_F, EL_price_D, EL_price_F, def_rate_D, def_rate_F):
    # D holds F-bonds (F-good claim -> convert with p); issuer = F
    rb_F_dg_next = (1 + rb_actual_F(+1)) * p(+1) / p - 1
    # IC-theory derived required premium: D-bank IC parameters govern D-bank's FOC on F-bonds
    # macro-pru-fix: EL_price_F = fundamental expected-loss loading on F-bonds (issuer=F),
    # independent of psi_lambda_B. See divert_bond_foc_D.
    prem_FD      = excess_return_F_D_ss + (EL_price_F + psi_spread_D) * def_rate_F(+1)
    # T-2 fix: deposit rate for the t->t+1 holding period is locked at t (rdep, not rdep(+1)).
    b_F_D_res    = (rb_F_dg_next - rdep_D) - prem_FD \
                   - psi_bF_D * (b_F_D - b_F_D_ss) - tau_mp_D
    # F holds D-bonds (D-good claim -> convert with p); issuer = D
    rb_D_fg_next = (1 + rb_actual_D(+1)) * p / p(+1) - 1
    # IC-theory derived required premium: F-bank IC parameters govern F-bank's FOC on D-bonds
    # macro-pru-fix: EL_price_D = fundamental expected-loss loading on D-bonds (issuer=D),
    # independent of psi_lambda_B. See divert_bond_foc_D.
    prem_DF      = excess_return_D_F_ss + (EL_price_D + psi_spread_F) * def_rate_D(+1)
    b_D_F_res    = (rb_D_fg_next - rdep_F) - prem_DF \
                   - psi_bD_F * (b_D_F - b_D_F_ss) - tau_mp_F
    return b_F_D_res, b_D_F_res
