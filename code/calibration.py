import json
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_EBA_MOMENTS = os.path.join(os.path.dirname(_HERE), "data", "eba_moments.json")


def load_eba_targets(path: str = _EBA_MOMENTS) -> dict:
    """Read the EBA 2011 moment set produced by ``code/eba_calibration.py``.

    Single source of truth: nothing here or in ``steady_state.py`` may carry its
    own copy of these numbers. (The retired ``audit_artifacts/`` harness did
    exactly that and silently tested a different model for weeks.) Regenerate
    with ``python code/eba_calibration.py``.
    """
    with open(path) as fh:
        return json.load(fh)["model_targets"]


# ─────────────────────────────────────────────────────────────────────────────
# EBA switch. True = the MEASURED EBA 2011 moment set (code/eba_calibration.py ->
# data/eba_moments.json). False = the pre-EBA placeholder calibration, kept
# bit-exact for regression comparison.
#
# STATUS (2026-07-31): the STEADY STATE is now correct under the measured moments
# — that was the collateral-mapping fix — but the DYNAMICS are still explosive,
# so the default stays False. Two separate problems, one solved:
#
# 1. SOLVED — collateral mapping. The GK block is well-posed only if
#        f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross
#    which needs Delta_own > ~0.73. The inherited Delta_own=0.2 violated it by
#    -1.26 (D) / -1.42 (F), giving lambda_gk_D=-0.087, Omega_D=-0.301 (negative
#    IC multiplier and franchise value) while the solver still converged with
#    machine-zero residuals. What had pinned Delta at 0.2/0.4 was an undocumented
#    convention inside ic_delta_calibration (ratio = Delta_cross/Delta_own = 2.0)
#    which with Delta_cross<=1 capped Delta_own at 0.5. That convention is gone;
#    Delta is now free and the module checks the IC residual directly. At
#    Delta = 0.85/0.90: lambda_gk_D=+0.927, Omega_D=+4.62, K_D=10.80 — and
#    +0.927 is essentially the pre-EBA +0.923, so the amplification block keeps
#    its previous strength with measured concentration.
#
# 2. OPEN — dynamic instability. The measured moments give a financial-accelerator
#    gain theta*phi_own = 13.17, versus 4*0.25 = 1.0 for the placeholder, and the
#    linearised system has an unstable root there: b_gov_D[499] ~ 1e2-1e3 instead
#    of ~1e-5. Diagnosed 2026-07-31 and NOT a fiscal or friction problem:
#      - present at psi_lambda_B = 0 (collateral friction fully off);
#      - flat in phi_lamb up to 25 (peak spread ~1.1e7bp at 0.6, 1.5 AND 25), so
#        it is not the debt/fiscal mode;
#      - mv_rule=1 does not fix it either.
#    `chi1` (intermediation adjustment cost, currently 0) is the strongest lever
#    found: chi1=0.5 cuts the peak spread 1.1e7bp -> 6.0bp and b_gov[499]
#    -2038 -> +70. But no chi1 in [0.2, 5.0] removes the root (b_gov[499] stays
#    70-560). It damps amplitude, not stability.
#
# See docs/eba_calibration.md "GK feasibility" and "Dynamic instability".
EBA_CALIBRATION = False


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
        'rdep_D':       0.000,   'rdep_F':       0.000,
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
        # RAISED 0.2/0.4 -> 0.85/0.90 by the 2026-07-31 EBA rebuild. Not a free
        # choice: the GK block is well-posed only if the banker's franchise value
        # covers the non-divertable part of the sovereign book,
        #     f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross
        # (steady_state.assert_gk_well_posed). At the MEASURED moments
        # (theta=5.51/6.94, phi_own=2.39/2.76, f=0.12) the old 0.2/0.4 violates it
        # by -1.26 (D) / -1.42 (F): lambda_gk and Omega go NEGATIVE while the
        # solver converges with machine-zero residuals. The measured concentration
        # therefore puts a LOWER BOUND on Delta_own of ~0.73 — EBA data partially
        # identifies a parameter with no direct empirical counterpart. The other
        # levers are out of range: f would need > 0.349 (literature 0.03-0.12),
        # theta > 16.03 (measured 5.51).
        #
        # ECONOMICS, and why this is a correction rather than a fudge: measured
        # leverage is only 5.5x on a book that is ~43% sovereign. You cannot also
        # claim sovereigns are excellent collateral — if they were, the bank would
        # lever further and theta=5.5 would not be the binding constraint. High
        # concentration at low leverage *implies* bonds are nearly as divertable as
        # capital. Consistent with 2010-12 Greece (collapsing GGB collateral
        # eligibility, rising ECB haircuts).
        #
        # 0.85 chosen from a Delta sweep (docs/eba_calibration.md): it delivers
        # lambda_gk_D=+0.927, essentially identical to the pre-EBA +0.923, so the
        # amplification block keeps its previous strength while the concentration
        # becomes measured. 0.90/0.95 also works (lambda_gk ~0.49/0.46) but halves
        # it. Delta_cross=0.90 preserves Delta_own < Delta_cross (own sovereign
        # still better collateral than foreign); cross-holdings are ~1% of the
        # book so this barely binds.
        'Delta_bD_D':   0.85 if EBA_CALIBRATION else 0.2,
        'Delta_bF_F':   0.85 if EBA_CALIBRATION else 0.2,
        'Delta_bF_D':   0.90 if EBA_CALIBRATION else 0.4,
        'Delta_bD_F':   0.90 if EBA_CALIBRATION else 0.4,
        'lambda_BD_D':  0.06,    'lambda_BF_F':  0.06,
        'lambda_BF_D':  0.06,    'lambda_BD_F':  0.06,
        # EBA 2011 REBUILD (2026-07-31). Supersedes both the 2026-07-22 EBA build
        # and the 2026-07-30 pre-EBA revert. All values below come from
        # data/eba_moments.json (regenerate: python code/eba_calibration.py);
        # see docs/eba_calibration.md for the identification ledger, including
        # what this moment set does NOT pin down.
        #
        # psi_lambda_B is the one amplification dial and remains UNIDENTIFIED by
        # EBA -- it is tuned to the 150bp GR-DE spread target. What changed is how
        # much work it has to do: the mechanical mark-to-market channel is now
        # measured (phi_own x ladder duration), so psi_lambda_B no longer stands
        # in for a mechanical loss that was ~10x too weak. See the MTM block in
        # data/eba_moments.json.
        'psi_lambda_B_D': 1.0 if EBA_CALIBRATION else 3.0,
        'psi_lambda_B_F': 1.0 if EBA_CALIBRATION else 3.0,
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
        'phi_lamb_D':   0.60 if EBA_CALIBRATION else 0.15,
        'phi_lamb_F':   0.60 if EBA_CALIBRATION else 0.15,
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
        'mv_rule_D':    1.0 if EBA_CALIBRATION else 0.0,
        'mv_rule_F':    1.0 if EBA_CALIBRATION else 0.0,
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
        'mu_w_D':       1.0,     'mu_w_F':       1.0,

        # ── SS Real Variables ─────────────────────────────────────────────────
        'mc_D':         1.0,     'mc_F':         1.0,

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
