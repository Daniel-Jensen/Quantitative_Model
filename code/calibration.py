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
        # Collateral-friction doom-loop sensitivity. Amplification ~ psi_lambda_B * phi_own.
        # With EBA phi_bD_D=2.39 (~10x the old 0.25), psi must fall ~10x (0.75/2.39≈0.31)
        # to keep comparable amplification AND keep Delta_eff = Delta + psi*def_rate below 1
        # (Delta_own back-solves to ~0.99; psi=3.0 makes def shocks cross 1 -> sign inversion).
        'psi_lambda_B_D': 0.31,  'psi_lambda_B_F': 0.31,
        # EBA 2011 (31 Dec 2010): CT1 / quarterly own-GDP. GR 22,778/55,898=0.408;
        # DE 114,317/653,815=0.175. (was 0.75*4=3.0 each — overstated bank equity ~7x.)
        'n_inter_D':    0.408,   'n_inter_F':    0.175,
        'theta_D':      4,       'theta_F':      4,
        # Bank capital-intermediation share: banks hold omega_K of the physical
        # capital stock; the residual (1-omega_K)K is held by a passive, deposit-
        # funded capital fund whose spread (rk-rdep) is rebated to households
        # (div_fund). omega_K=1 recovers the original all-capital-in-banks model.
        # Set to preserve K≈10.8 (K/annualY≈2.7) against EBA-thin net worth:
        # omega_K = n·(theta - phi_own - phi_cross)/(Q·K_target). Recomputed exactly
        # in steady_state.py. SS-neutral at omega_K=1 (div_fund=0). See docs/eba_calibration.md.
        'omega_K_D':    0.0602,  'omega_K_F':    0.0190,

        # ── Bellman nu risk-discount ───────────────────────────────────────────
        'psi_nu_bD_D':  0.0,     'psi_nu_bD_F':  0.0,
        'psi_nu_bF_D':  0.0,     'psi_nu_bF_F':  0.0,

        # ── Fiscal & Government Debt ──────────────────────────────────────────
        # Bank-held government debt (decision b1): GR+DE banks held 62.4bn GR / 323bn DE
        # of sovereign (27.9% / 12.1% of annual own GDP). This is the bank-channel debt,
        # NOT headline debt/GDP (~170% GR) — the non-bank/ECB residual enters via SMP
        # (phase 2). B_supply_D = b_D_D + b_D_F ≈ 1.19; B_supply_F ≈ 0.591 (quarterly).
        'B_supply_D':   1.19,    'B_supply_F':   0.591,
        'b_gov_D':      1.19,    'b_gov_F':      0.591,
        'b_gov_ss_D':   1.19,    'b_gov_ss_F':   0.591,

        # ── Fiscal Rule ───────────────────────────────────────────────────────
        'tau_D':        0.181,   'tau_F':        0.181,
        'lamb_D':       0.85,    'lamb_F':       0.85,
        'lamb_ss_D':    0.85,    'lamb_ss_F':    0.85,
        # phi_lamb raised from 0.02 after T-2 fix: deposit re-dating makes the
        # debt→spread spiral live; phi_lamb < ~0.12 is explosive at current amplification.
        # EBA sovereign exposure (phi_bD_D=2.39) massively amplifies the doom loop,
        # so the fiscal feedback must be stronger than the old 0.15 (F-1). Raised to
        # 0.30; market-value rule (mv_rule=1) on to restore stationarity.
        'phi_lamb_D':   0.60,    'phi_lamb_F':   0.60,
        # Fiscal-rule debt measure: 0 = par/face value (default), 1 = market value
        # (q_b·b_gov(-1)). mv_gov_ss is recomputed exactly from the solved SS in
        # build_and_solve; these are placeholders (unused when mv_rule=0).
        # mv_rule=1 (F-1 fix): market-value rule restores stationarity under the
        # high EBA sovereign exposure. writeoff stays OFF (risk-premium framing).
        'mv_rule_D':    1.0,     'mv_rule_F':    1.0,
        'mv_gov_ss_D':  0.6*4,   'mv_gov_ss_F':  0.6*4,

        # ── Sovereign Default ─────────────────────────────────────────────────
        'shock_def_D':      0.000,  'shock_def_F':      0.0,
        'T_ls_D':           0.000,  'T_ls_F':           0.000,
        'def_rate_D':       0.000,  'def_rate_F':       0.0,
        'def_scale_D':      0.25,   'def_scale_F':      0.25,
        'def_curvature_D':  0.5,    'def_curvature_F':  0.5,
        'def_offset_D':     0.05,   'def_offset_F':     0.05,
        'recovery_rate_D':  0.00,   'recovery_rate_F':  0.00,
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
        # EBA 2011 cross-holdings / capital: GR banks' Bund 411/22,778=0.018;
        # DE banks' GR 7,934/114,317=0.069 (was 0.25 each — overstated the direct
        # cross-border channel ~10-30x). Own-holdings set in steady_state.py.
        'phi_bF_D_ss':  0.018,   'phi_bD_F_ss':  0.069,
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
