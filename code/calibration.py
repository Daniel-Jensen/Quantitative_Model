import numpy as np


def get_calibration():
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
        'delta_b_D':    0.10,    'delta_b_F':    0.10,

        # ── Aggregate Targets (SS) ────────────────────────────────────────────
        'Y_D':          1.00,    'Y_F':          1.00,
        'Y_ss_D':       1.0,     'Y_ss_F':       1.0,
        'N_D':          1.00,    'N_F':          1.00,
        'w_D':          0.65,    'w_F':          0.65,

        # ── Financial Intermediaries (Gertler-Karadi) ─────────────────────────
        'f_D':          0.12,    'f_F':          0.12,
        'lambda_gk_D':  0.2,     'lambda_gk_F':  0.2,
        'beta_inter_D': 0.9975155088,  'beta_inter_F': 0.9975155088,
        'Delta_bD_D':   0.2,     'Delta_bF_F':   0.2,
        'Delta_bF_D':   0.4,     'Delta_bD_F':   0.4,
        'lambda_BD_D':  0.06,    'lambda_BF_F':  0.06,
        'lambda_BF_D':  0.06,    'lambda_BD_F':  0.06,
        # PRE-EBA REVERT (2026-07-30): calibration rolled back to the values in
        # force at abcbb6e, the last commit before the EBA-calibration work
        # (eade414 onward). Structural fixes (C-1 multi-asset lambda_gk, W-1/W-2,
        # T-2, omega_K capital-fund generalisation) are RETAINED -- only parameter
        # values move. Every EBA value is recorded in the adjacent comment so the
        # roll-forward is mechanical. See docs/eba_calibration.md for the EBA map.
        #
        # psi_lambda_B: original default (EBA-calibrated value was 1.1793, tuned to
        # the 150bp GR-DE spread target). NOTE the breakdown warning in
        # docs/eba_calibration.md applies to the EBA calibration's thin net worth;
        # at n_inter=3.0 the pre-EBA model sat comfortably at 3.0.
        'psi_lambda_B_D': 3.0,   'psi_lambda_B_F': 3.0,
        # Bank net worth: 0.75 of annual GDP. (EBA 2011 values were 0.408 / 0.175,
        # = CT1 / quarterly own-GDP; GR 22,778/55,898, DE 114,317/653,815.)
        'n_inter_D':    0.75*4,  'n_inter_F':    0.75*4,
        'theta_D':      4,       'theta_F':      4,
        # Bank capital-intermediation share. omega_K=1 = all capital in banks, the
        # pre-EBA balance sheet: the passive capital fund is empty and div_fund=0,
        # so capital_fund_D/F and the fund terms in smart_steady_D/F are inert.
        # (EBA values were 0.0602 / 0.0190, recomputed in steady_state.py to hold
        # K≈10.8 against EBA-thin net worth.)
        'omega_K_D':    1.0,     'omega_K_F':    1.0,

        # ── Bellman nu risk-discount ───────────────────────────────────────────
        'psi_nu_bD_D':  0.0,     'psi_nu_bD_F':  0.0,
        'psi_nu_bF_D':  0.0,     'psi_nu_bF_F':  0.0,

        # ── Fiscal & Government Debt ──────────────────────────────────────────
        # Government debt at 60% of annual GDP (quarterly units). These are
        # placeholders: steady_state.py overwrites all three from the solved
        # portfolio-share targets. (EBA values were 1.19 / 0.591 — bank-held
        # sovereign only, 27.9% / 12.1% of annual own GDP.)
        'B_supply_D':   0.6*4,   'B_supply_F':   0.6*4,
        'b_gov_D':      0.6*4,   'b_gov_F':      0.6*4,
        'b_gov_ss_D':   0.6*4,   'b_gov_ss_F':   0.6*4,

        # ── Fiscal Rule ───────────────────────────────────────────────────────
        'tau_D':        0.181,   'tau_F':        0.181,
        'lamb_D':       0.85,    'lamb_F':       0.85,
        'lamb_ss_D':    0.85,    'lamb_ss_F':    0.85,
        # phi_lamb raised from 0.02 after T-2 fix: deposit re-dating makes the
        # debt→spread spiral live; phi_lamb < ~0.12 is explosive at current amplification.
        # (EBA calibration raised this to 0.60 because phi_bD_D=2.39 massively
        # amplifies the doom loop; at pre-EBA exposures 0.15 was the committed value.
        # Caveat: F-1 later identified [0.15,0.18] as a near-unit-root zone on the
        # EBA-anchored model — watch the IRF tails here.)
        'phi_lamb_D':   0.15,    'phi_lamb_F':   0.15,
        # Fiscal-rule debt measure: 0 = par/face value (default), 1 = market value
        # (q_b·b_gov(-1)). mv_gov_ss is recomputed exactly from the solved SS in
        # build_and_solve; these are placeholders (unused when mv_rule=0).
        # Par/face-value rule (pre-EBA). mv_rule=1 was tried on 2026-07-30 and
        # REVERTED: at the pre-EBA phi_lamb=0.15 it lands in F-1's near-unit-root
        # zone [0.15,0.18] (which was identified UNDER mv_rule=1) and the model
        # breaks -- n_inter_D[0]=-1554%, Y_D[0]=+0.17% (perverse sign),
        # b_gov_D[499]=1.6e-2. Verified directly: mv_rule=1 needs phi_lamb=0.60
        # (main's EBA-era pairing) to stay healthy (n_inter_D[0]=-5.89%,
        # Y_D[0]=-0.024%, b_gov_D[499]=0.0). mv_rule=1 and phi_lamb=0.15 are NOT
        # a usable pair. Consequence of mv_rule=0: empirical duration
        # (delta_b=0.036/0.038) stays blocked -- F-1 finds the par rule explosive
        # at every phi_lamb in [0.02,0.50] there -- and diagnostics/regimes/
        # assumes the market-value rule in its provenance string.
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
        # Symmetric cross-border share. Own-holdings set in steady_state.py.
        # (EBA 2011 cross-holdings / capital were 0.018 GR-holds-Bund and
        # 0.069 DE-holds-GR — a far thinner direct contagion channel.)
        'phi_bF_D_ss':  0.25,    'phi_bD_F_ss':  0.25,
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
