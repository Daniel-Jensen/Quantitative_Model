import json
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_EBA_MOMENTS = os.path.join(os.path.dirname(_HERE), "data", "eba_moments.json")


# Which object the model's intermediary represents.
#   "broad" (default) = the whole capital-funding sector. n_inter follows from the
#       MEASURED leverage and the balance sheet, N = (Q*K + sovereign)/theta, so
#       omega_K = 1 and the passive-fund device disappears entirely.
#   "ct1"  = Core Tier 1 of the EBA stress-test sample. Historical; this is the
#       scope that makes omega_K tiny and the accelerator gain ~1/n_inter, i.e.
#       the source of the dynamic instability. Kept for comparison only.
BANK_SCOPE = "broad"


def load_eba_targets(path: str = _EBA_MOMENTS, scope: str | None = None) -> dict:
    """Read the EBA 2011 moment set produced by ``code/eba_calibration.py``.

    Single source of truth: nothing here or in ``steady_state.py`` may carry its
    own copy of these numbers. (The retired ``audit_artifacts/`` harness did
    exactly that and silently tested a different model for weeks.) Regenerate
    with ``python code/eba_calibration.py``.
    """
    key = {"broad": "model_targets_broad", "ct1": "model_targets"}[scope or BANK_SCOPE]
    with open(path) as fh:
        return json.load(fh)[key]


# ─────────────────────────────────────────────────────────────────────────────
# EBA switch. True = the MEASURED EBA 2011 moment set (code/eba_calibration.py ->
# data/eba_moments.json), read at the scope set by BANK_SCOPE above.
# False = the pre-EBA placeholder calibration, kept bit-exact for regression.
#
# LIVE since 2026-07-31, once BANK_SCOPE="broad" resolved the last blocker.
# The three problems found and fixed getting here:
#
# 1. Collateral mapping. ic_delta_calibration closed its Delta back-solve with a
#    hidden ratio = Delta_cross/Delta_own = 2.0, which capped Delta_own <= 0.5
#    against a GK requirement of > ~0.73 under the CT1 scope. Removed; Delta is
#    free and the IC residual is checked directly.
# 2. omega_K as a FIXED SHARE made the passive fund mirror bank deleveraging, so
#    dK/dN = theta/omega_K. New fund_rule=1 (fixed quantity) gives dK/dN = theta
#    with an identical steady state.
# 3. n_inter scope. CT1 of the stress-test sample is not the net worth of the
#    agent intermediating the whole capital stock. Under the broad scope
#    n_inter = (Q*K + sovereign)/theta, omega_K = 1, and the model is STABLE:
#    b_gov_D[499] ~ 1e-9 to 1e-6 with both impact signs correct.
#
# See docs/eba_calibration.md.
EBA_CALIBRATION = True


def get_calibration():
    _eba_json = load_eba_targets()
    # `eba` returns the measured EBA value only when the switch is on;
    # otherwise the pre-EBA fallback passed as the second argument.
    def eba_or(key, pre_eba):
        return _eba_json[key] if EBA_CALIBRATION else pre_eba
    calibration_start = {

        # ── Preferences ───────────────────────────────────────────────────────
        'frisch_D':     0.50,    'frisch_F':     0.50,
        'eis_D':        0.5,     'eis_F':        0.5,

        # ── Rates & Asset Prices ──────────────────────────────────────────────
        # Nominal deposit rate. Deposits are nominal euro contracts; the derived
        # real rates rdep_D/F (ex-ante) and rdep_expost_D/F (realised) come from
        # deposit_rates_D/F. At SS pi = 0, so rdep = i_dep and the SS is
        # unchanged from the real-deposit calibration.
        'i_dep_D':      0.000,   'i_dep_F':      0.000,
        'q_b_D':        0.83,    'q_b_F':        0.83,
        'Q_D':          1.0,     'Q_F':          1.0,

        # ── Production ────────────────────────────────────────────────────────
        'alpha_D':      0.35,    'alpha_F':      0.35,
        'delta_D':      0.025,   'delta_F':      0.025,
        'ksi_D':        0.50,    'ksi_F':        0.50,

        # ── Long-term bonds ───────────────────────────────────────────────────
        # EBA REBUILD (2026-07-31): delta_b is now MEASURED, from the sovereign
        # maturity ladder (EBA worksheet 5, MATURITY_CODE 125..155) repriced at
        # the 31-Dec-2010 market yield. GR banks' GGB book: 5.13y weighted-average
        # residual maturity but only 3.12y MODIFIED DURATION (a 12% discount rate
        # pulls duration far below maturity); DE banks' Bund book: 4.86y / 4.22y.
        # Inverting through the HM perpetuity gives 0.0777 (D) / 0.0568 (F).
        #
        # This retires the long-standing "empirical duration is 7y, port
        # delta_b=0.036/0.038 from bank-cal" item. That target was the average
        # residual maturity of the SOVEREIGN's whole outstanding stock, which is
        # the wrong object: delta_b governs the duration of the book sitting on
        # BANK balance sheets, at the yields those banks actually faced. The
        # F-1 blocker (porting 0.036 needs mv_rule=1 AND phi_lamb=0.60 jointly)
        # therefore does not bind — 0.0777/0.0568 is close to the old 0.10 and
        # runs under the par rule.
        'delta_b_D':    eba_or('delta_b_D', 0.10),   'delta_b_F':    eba_or('delta_b_F', 0.10),

        # ── Aggregate Targets (SS) ────────────────────────────────────────────
        'Y_D':          1.00,    'Y_F':          1.00,
        'Y_ss_D':       1.0,     'Y_ss_F':       1.0,
        'N_D':          1.00,    'N_F':          1.00,
        'w_D':          0.65,    'w_F':          0.65,

        # ── Financial Intermediaries (Gertler-Karadi) ─────────────────────────
        'f_D':          0.12,    'f_F':          0.12,
        'lambda_gk_D':  0.2,     'lambda_gk_F':  0.2,
        'beta_inter_D': 0.9975155088,  'beta_inter_F': 0.9975155088,
        # Divertability of sovereign bonds in the multi-asset IC (higher = worse
        # collateral). Genuine hardcoded structural inputs since the C-1 fix.
        #
        # 0.2/0.4 under BANK_SCOPE="broad". These were briefly raised to 0.85/0.90
        # while the CT1 scope was live: GK well-posedness needs
        #     f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross
        # (steady_state.assert_gk_well_posed), and CT1's phi_own=2.39 violated it
        # at 0.2, sending lambda_gk and Omega NEGATIVE while the solver still
        # converged with machine-zero residuals. That bound was an artifact of the
        # CT1 scope. At the broad-sector phi_own=0.456 the constraint reads
        # f*theta = 0.661 > 0.8*0.456 + 0.6*0.0034 = 0.367 — comfortably satisfied,
        # so the paper's "sovereigns are good collateral" story is restored.
        # Under BANK_SCOPE="ct1" these must be raised above ~0.73.
        #
        # Still UNIDENTIFIED: no EBA counterpart, no moment attached. The
        # feasibility inequality bounds them but does not pin a level.
        'Delta_bD_D':   0.2,     'Delta_bF_F':   0.2,
        'Delta_bF_D':   0.4,     'Delta_bD_F':   0.4,
        'lambda_BD_D':  0.06,    'lambda_BF_F':  0.06,
        'lambda_BF_D':  0.06,    'lambda_BD_F':  0.06,
        # EBA 2011 REBUILD (2026-07-31). Supersedes both the 2026-07-22 EBA build
        # and the 2026-07-30 pre-EBA revert. All values below come from
        # data/eba_moments.json (regenerate: python code/eba_calibration.py);
        # see docs/eba_calibration.md for the identification ledger, including
        # what this moment set does NOT pin down.
        #
        # psi_lambda_B: the one amplification dial, still UNIDENTIFIED by EBA and
        # tuned to the paper's 150bp GR-DE spread target on a 1pp default shock.
        # RETUNED 2026-07-31 for BANK_SCOPE="broad". The mapping is smooth and
        # monotone with no breakdown region — peak annualised spread
        #     5.4 / 23.1 / 58.6 / 111.0 / 142.5 / 157.8 bp
        # at psi_lambda_B = 0 / 1 / 3 / 6 / 8 / 9 — and b_gov_D[499] stays in
        # ~1e-9..1e-4 throughout. 8.5 interpolates to ~150bp.
        # Higher than the historical 0.31 / 1.18 / 3.0 because the broad scope's
        # phi_own = 0.456 is far below the CT1 scope's 2.39, so more of the
        # default loading has to come from the friction. NOTE the old "breakdown
        # above ~1.5-2.0" warning was specific to CT1-thin net worth and does not
        # apply here.
        # Sweep method caveat: psi_spread_D is derived from psi_lambda_B inside
        # _apply_ss_anchors, so a sweep MUST re-solve the SS per point. Patching
        # the flag on an already-solved SS leaves psi_spread stale and inverts the
        # apparent sign of the spread response.
        #
        # RETUNED 2026-08-06 (Task 14, add-nkpc): sticky prices (NKPC blocks) and
        # nominal (non-state-contingent) deposit contracts both raise spread
        # transmission, moving the old 8.5 -> 162.0 bp (was 150.4 bp pre-change,
        # ~8% over target). Re-bisected against the same 150bp GR-DE peak-spread
        # moment on a 1pp default shock, holding everything else fixed:
        #   psi_lambda_B = 8.5  -> peak spread 0.4053 pp -> 162.14 bp
        #   psi_lambda_B = 7.0  -> peak spread 0.3405 pp -> 136.21 bp
        #   psi_lambda_B = 7.8  -> peak spread 0.3729 pp -> 149.16 bp
        #   psi_lambda_B = 7.85 -> peak spread 0.3753 pp -> 150.14 bp  <- adopted
        # b_gov_D[499] stayed in ~1e-5..1e-4 across the whole bracket (no
        # instability); n_inter_D[0] and Y_D[0] both negative throughout
        # (correct doom-loop sign). See docs/STATE.md for the full record.
        'psi_lambda_B_D': 7.85 if EBA_CALIBRATION else 3.0,
        'psi_lambda_B_F': 7.85 if EBA_CALIBRATION else 3.0,
        # Bank net worth = Core Tier 1 / own quarterly nominal GDP.
        # GR 22,778/55,898 = 0.4075; DE 114,317/653,815 = 0.1748.
        'n_inter_D':    eba_or('n_inter_D', 0.75*4),  'n_inter_F':    eba_or('n_inter_F', 0.75*4),
        # GK leverage on the GK-ELIGIBLE book: (corporate ex-CRE + commercial real
        # estate + sovereign) EAD / CT1, own-country. 5.51 (D) / 6.94 (F).
        # NOT CT1/total assets (14.9/32.9) -- theta multiplies only the GK book,
        # and the total-assets version was previously verified not to converge.
        'theta_D':      eba_or('theta_D', 4.0),    'theta_F':      eba_or('theta_F', 4.0),
        # Bank share of the capital stock, MEASURED: (corporate + CRE) EAD / K,
        # with K from the conventional K/Y_annual = 2.7. 0.117 (D) / 0.067 (F).
        # This is no longer the back-solved residual that made an ASSUMED
        # theta=4.0 consistent with a K target; theta and omega_K now come from
        # the same observed balance sheet, and the resulting K is an
        # over-identifying check printed by steady_state.py (expect ~10.8).
        'omega_K_D':    eba_or('omega_K_D', 1.0),  'omega_K_F':    eba_or('omega_K_F', 1.0),
        # Capital-fund behaviour. THIS IS THE FIX for the EBA dynamic instability.
        #   0 = FIXED SHARE  (legacy): fund holds (1-omega_K)*K, so it mechanically
        #       mirrors bank deleveraging and dK/dN = theta/omega_K.
        #   1 = FIXED QUANTITY: fund holds a constant K_fund, bank is the marginal
        #       holder, dK/dN = theta.
        # Identical steady state when K_fund = (1-omega_K)*K_ss; the difference is
        # purely dynamic. Under fixed share, holding K/Y at its conventional target
        # forces omega_K = N(theta-phi)/K, hence dK/dN = theta*K/(N*(theta-phi)),
        # i.e. the accelerator gain is INVERSELY PROPORTIONAL TO BANK NET WORTH.
        # Measured CT1 is 7.4x thinner than the placeholder, which is what made the
        # EBA calibration explosive (b_gov[499] ~ 1e2-1e3). Verified: with the
        # pre-EBA bank block (omega_K=1) the same model is stable at ~1e-8
        # regardless of concentration or Delta. See docs/eba_calibration.md.
        # A fixed share is also the harder assumption to defend: it says non-bank
        # capital holders shrink in lockstep with bank equity, which is the
        # amplification, not an independent behavioural claim.
        'fund_rule_D':  1.0,     'fund_rule_F':  1.0,
        # Fund's fixed capital holding, = (1-omega_K)*K_target with K_target=10.8
        # (K/Y_annual=2.7). Unused when fund_rule=0. At omega_K=1 the fund is empty.
        'K_fund_D':     (1.0 - eba_or('omega_K_D', 1.0)) * 10.8,
        'K_fund_F':     (1.0 - eba_or('omega_K_F', 1.0)) * 10.8,

        # ── Bellman nu risk-discount ───────────────────────────────────────────
        'psi_nu_bD_D':  0.0,     'psi_nu_bD_F':  0.0,
        'psi_nu_bF_D':  0.0,     'psi_nu_bF_F':  0.0,

        # ── Fiscal & Government Debt ──────────────────────────────────────────
        # BANK-HELD sovereign stock / own quarterly GDP: 1.116 (D) / 0.483 (F),
        # i.e. 27.9% / 12.1% of annual own GDP. NOT headline debt/GDP (~150% GR):
        # the non-bank / official / ECB residual is outside the model's bank
        # block by construction. These are start values; steady_state.py
        # overwrites all three from the solved portfolio-share targets.
        'B_supply_D':   eba_or('B_supply_D_qgdp', 0.6*4),  'B_supply_F':   eba_or('B_supply_F_qgdp', 0.6*4),
        'b_gov_D':      eba_or('B_supply_D_qgdp', 0.6*4),  'b_gov_F':      eba_or('B_supply_F_qgdp', 0.6*4),
        'b_gov_ss_D':   eba_or('B_supply_D_qgdp', 0.6*4),  'b_gov_ss_F':   eba_or('B_supply_F_qgdp', 0.6*4),

        # ── Fiscal Rule ───────────────────────────────────────────────────────
        'tau_D':        0.181,   'tau_F':        0.181,
        'lamb_D':       0.85,    'lamb_F':       0.85,
        'lamb_ss_D':    0.85,    'lamb_ss_F':    0.85,
        # ~Bohn (1998) fiscal-feedback magnitude. 0.60, not the pre-EBA 0.15:
        # phi_bD_D=2.39 (measured) amplifies the doom loop roughly 10x relative to
        # the 0.25 placeholder, and 0.60 is what the 2026-07-22 EBA build needed
        # for stationarity. Also clears F-1's near-unit-root zone [0.15,0.18] by a
        # wide margin. Not moment-matched — see the identification ledger.
        'phi_lamb_D':   0.15,    'phi_lamb_F':   0.15,
        # Fiscal-rule debt measure: 0 = par/face value (default), 1 = market value
        # (q_b·b_gov(-1)). mv_gov_ss is recomputed exactly from the solved SS in
        # build_and_solve; these are placeholders (unused when mv_rule=0).
        # Market-value rule REQUIRED under the EBA calibration; par rule for the
        # pre-EBA placeholder. Measured 2026-07-31: at mv_rule=0 with the measured
        # concentration the debt path explodes even with the collateral friction
        # switched off entirely (psi_lambda_B=0 gives b_gov_D[499]=6.4e+03 and a
        # nonsense peak spread), so this is the DEBT/fiscal mode, not amplification.
        # The driver is phi_own=2.39 (a ~10x stronger doom loop than the 0.25
        # placeholder), not duration — which is why the "measured delta_b is close
        # to the old 0.10, so the par rule is fine" reasoning did not carry: it was
        # about duration. mv_rule=1 + phi_lamb=0.60 is the pairing the 2026-07-22
        # EBA build verified stationary, and F-1's hard break at
        # mv_rule=1 + phi_lamb=0.15 is far away.
        'mv_rule_D':    0.0,     'mv_rule_F':    0.0,
        'mv_gov_ss_D':  0.6*4,   'mv_gov_ss_F':  0.6*4,

        # ── Sovereign Default ─────────────────────────────────────────────────
        'shock_def_D':      0.000,  'shock_def_F':      0.0,
        'T_ls_D':           0.000,  'T_ls_F':           0.000,
        'def_rate_D':       0.000,  'def_rate_F':       0.0,
        'def_scale_D':      0.25,   'def_scale_F':      0.25,
        'def_curvature_D':  0.5,    'def_curvature_F':  0.5,
        'def_offset_D':     0.05,   'def_offset_F':     0.05,
        # EL-1 resolution retained through the 2026-07-30 pre-EBA revert (author
        # decision): 0.30 = Greek PSI NPV-recovery framing (Zettelmeyer, Trebesch &
        # Gulati, PIIE WP13-8; 59-65% investor NPV loss). The pre-EBA value was 0.00,
        # i.e. 100% loss-given-default -- counterfactual for Greece, and the *harshest*
        # possible assumption rather than a neutral one.
        # While writeoff_enabled=0 this is live ONLY through EL_price (the realized-
        # haircut terms in bond_return/government_ss/budget_residual are gated by
        # writeoff_enabled). EL_price = (1-rec)*delta_b/q_b: 0.1025 at rec=0.00 ->
        # 0.0717 at rec=0.30. Against psi_spread=0.8385 that moves total default
        # loading by only ~3.3%.
        'recovery_rate_D':  0.30,   'recovery_rate_F':  0.30,
        'zeta_writeoff_D':  0.0,    'zeta_writeoff_F':  0.0,
        'writeoff_enabled_D': 0.0,  'writeoff_enabled_F': 0.0,

        # ── ECB balance sheet (TPI conduit) ───────────────────────────────────
        # Capital-key split of the CB's D-bond programme cash flows between the
        # two treasuries. kappa_cb_F = F share of the two-country renormalised
        # euro-area capital key: Bundesbank 26.1% / Bank of Greece 2.0% of the
        # euro-area key -> 26.1/28.1 ≈ 0.929. Read only by the TPI layer;
        # SS-neutral (cb_buy_ss = 0).
        'kappa_cb_F':       0.929,

        # ── Intermediary Capital Adjustment Cost ──────────────────────────────
        'chi0_D':           0.00,   'chi0_F':           0.00,
        'chi1_D':           0.00,   'chi1_F':           0.00,
        'chi2_D':           2.0,    'chi2_F':           2.0,

        # ── Macroprudential Bond Tax ──────────────────────────────────────────
        'T0_D':             0.000,  'T0_F':             0.000,
        'T1_D':             0.0,    'T1_F':             0.0,

        # ── Trade & Terms of Trade ────────────────────────────────────────────
        'omega':            0.85,
        'epsilon_trade':    1.5,
        'p':                0.50,

        # ── Cross-Border Bond Portfolio ───────────────────────────────────────
        # EBA 2011 cross-holdings / capital: GR banks' Bund book 411/22,778 =
        # 0.018; DE banks' GGB book 7,934/114,317 = 0.069. A 10-30x thinner
        # direct contagion channel than the 0.25 symmetric placeholder.
        # Own-holdings set in steady_state.py from the same moment file.
        'phi_bF_D_ss':  eba_or('phi_bF_D_ss', 0.25),  'phi_bD_F_ss':  eba_or('phi_bD_F_ss', 0.25),
        'psi_bF_D':     0.5,     'psi_bD_F':     0.5,

        # ── Wage Markups ──────────────────────────────────────────────────────
        # Unchanged: wages are flexible. mu_w = 1 is the SS-neutralising device
        # in labor_ss_D/F; there is no wage Phillips curve.
        'mu_w_D':       1.0,     'mu_w_F':       1.0,

        # ── Price Rigidity (Rotemberg) ────────────────────────────────────────
        # mu_p: gross price markup, epsilon_p = 6. FREE TO FIRST ORDER under the
        #   subsidy neutralisation -- the gap (mu_p*mc - 1) linearises to mc_hat
        #   for any mu_p -- so this needs no defending unless live markups are
        #   ever adopted.
        # mc: SS real marginal cost = 1/mu_p. The production subsidy
        #   tau_s = 1 - 1/mu_p makes labour demand collapse to the competitive
        #   w = (1-alpha)Y/N at this value, so the SS is bit-identical to flex.
        # kappa_p: Calvo theta_p = 0.75 at beta = 0.985, slope
        #   (1-theta)(1-beta*theta)/theta = 0.0871. Euro-area IPN median price
        #   duration ~4 quarters (Alvarez et al. 2006; Dhyne et al. 2006).
        #   Agrees with Bi-Foerster-Traum's implied 0.0846 to within 3%.
        # pi: SS producer-price inflation, exactly zero.
        'mu_p_D':       1.20,    'mu_p_F':       1.20,
        'mc_D':    1.0 / 1.20,   'mc_F':    1.0 / 1.20,
        'kappa_p_D':    0.0871,  'kappa_p_F':    0.0871,
        'pi_D':         0.0,     'pi_F':         0.0,

        # omega_pi_D: weight on D in the union producer-price aggregate that the
        # ECB is assumed to stabilise. = 1 - kappa_cb_F, the renormalised
        # two-country capital key (BuBa 26.1 / BoG 2.0 of the euro-area key).
        # DO NOT use model GDP weights: the model normalises Y_D_ss ~ Y_F_ss ~ 1,
        # so they would give ~0.5 and split the terms-of-trade adjustment evenly
        # between Greek deflation and German inflation -- counterfactual for
        # 2010-12. Load-bearing twice over once deposits are nominal, since it
        # scales pi_D and hence the Fisher revaluation on bank balance sheets.
        'omega_pi_D':   0.071,

        # ── Idiosyncratic Income Process (Rouwenhorst) ────────────────────────
        'rho_z_D':  0.90,    'rho_z_F':  0.90,
        'sigma_z_D': 0.3,    'sigma_z_F': 0.3,
        'nZ_D':     15,      'nZ_F':     15,
        'nDep_D':   500,     'nDep_F':   500,
        'Depmax_D': 150,     'Depmax_F': 150,
    }

    # ── Bond Holdings: initial SS guess ──────────────────────────────────────
    _n_D = calibration_start['n_inter_D']
    _n_F = calibration_start['n_inter_F']
    _B_D = calibration_start['B_supply_D']
    _B_F = calibration_start['B_supply_F']

    b_F_D = calibration_start['phi_bF_D_ss'] * _n_D / calibration_start['q_b_F']
    b_D_F = calibration_start['phi_bD_F_ss'] * _n_F / calibration_start['q_b_D']

    calibration_start.update({
        'b_F_D': b_F_D,         'b_D_F': b_D_F,
        'b_D_D': _B_D - b_D_F,  'b_F_F': _B_F - b_F_D,
        'b_F_D_anchor': b_F_D,  'b_D_F_anchor': b_D_F,
        'psi_bD_D': 0.0,        'psi_bF_F': 0.0,
    })

    return calibration_start
