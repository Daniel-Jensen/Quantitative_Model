"""
Notes:
* Bond denomination: D-bonds are D-good claims (priced by D-bank using rdep_D).
  F-bonds are F-good claims (priced by F-bank using rdep_F).
  Cross-border values: D-bank's F-bond leg in D-goods = p·Q_bF·b_F_D;
                       F-bank's D-bond leg in F-goods = Q_bD·b_D_F/p.
"""


def get_calibration():
    cal = dict(

        # ── HOUSEHOLD PREFERENCES (GHH) ─────────────────────────────────────
        sigma_D=2.0,   sigma_F=2.0,
        frisch_D=0.5,  frisch_F=0.5,
        chi_D=0.5417,  chi_F=0.5417,   # warm start only — OVERWRITTEN in the
                                        # SS solve to pin N_ss = 1 (steady_state.py)

        # ── IDIOSYNCRATIC INCOME PROCESS (Rouwenhorst) ───────────────────────
        n_e_D=2,     n_e_F=2,
        rho_e_D=0.9,   sigma_e_D=0.2,
        rho_e_F=0.9,   sigma_e_F=0.2,

        # ── ASSET GRIDS ──────────────────────────────────────────────────────
        a_min_D=0.0, a_max_D=87.2, n_a_D=250, a_curve_D=2.0,
        a_min_F=0.0, a_max_F=87.2, n_a_F=250, a_curve_F=2.0,

        # ── FIRMS (Cobb-Douglas + full price flexibility) ────────────────────
        epsilon_D=6.0, epsilon_F=6.0,   # demand elasticity → mc = (ε-1)/ε
        Z_ss_D=0.45,   Z_ss_F=0.45,     # warm start only — OVERWRITTEN in the SS
                                          # solve to pin Y_ss = 1 (steady_state.py)

        # ── CAPITAL (Jermann 1998 adjustment cost) ────────────────────────────
        alpha_D=0.35,  alpha_F=0.35,    # capital share
        delta_D=0.025, delta_F=0.025,
        ksi_D=0.50,    ksi_F=0.50,      # adjustment-cost curvature

        # ── FINANCIAL INTERMEDIARY ─────────────────────────────────────────────
        f_D=0.05,              f_F=0.05,
        r_dep_D_target=0.000,   r_dep_F_target=0.000,
        beta_inter_D=0.96,      beta_inter_F=0.96,

        # Gertler-Karadi calibration targets: the single divertability λ and the
        # entrant transfer ω_ent are SOLVED (calibrate_bank_targets in bank.py) so
        # the deterministic SS hits these two moments.  φ = QS/N = total assets/equity.
        leverage_target_D=4.0,          leverage_target_F=4.0,
        credit_spread_target_D=0.005,   credit_spread_target_F=0.005,  # 200 bps/yr ÷ 4 (quarterly)

        # warm start only — OVERWRITTEN by calibrate_bank_targets in steady_state.py
        lambda_K_D=0.22,        lambda_K_F=0.22,
        lambda_bD_D=0.22,       lambda_bD_F=0.22,
        lambda_bF_D=0.22,       lambda_bF_F=0.22,
        omega_ent_D=0.002,      omega_ent_F=0.002,

        entrant_mode="proportional",
        phi_entry=0.0,

        # ── PORTFOLIO ADJUSTMENT COSTS (cross-border bonds) ──────────────────
        # b_D_F_ss / b_F_D_ss ≈ 20% of the respective bond supply: foreign
        # banks' pre-crisis holdings of peripheral debt (union contagion leg).
        psi_bF_D=0.05,          psi_bD_F=0.05,
        b_F_D_ss=0.744,         b_D_F_ss=0.744,
        excess_return_F_D_ss=0.0,               # overwritten after SS solve
        excess_return_D_F_ss=0.0,               # overwritten after SS solve

        # ── GOVERNMENT BONDS ─────────────────
        # delta_b = quarterly amortization rate
        delta_b_D=0.25,        delta_b_F=0.25,
        B_gov_D_ss=3.722,       B_gov_F_ss=3.722,   # 93% of annual GDP (= 0.93×4×Y_ss=1)

        # ── DEFAULT RISK (Cole-Kehoe zones × Bocola pricing) ─────────────────
        # Thresholds are FACE-value debt to quarterly Y_ss ratios.  F is always
        # safe.  SS sits inside the D crisis zone (b/Y_ss ≈ 3.7), so the
        # sunspot is priced immediately; b_ck_high is out of reach (no
        # fundamental default in the risk-only experiment).
        b_ck_low_D=3.00,        b_ck_low_F=99.0,
        b_ck_high_D=6.00,       b_ck_high_F=99.0,
        # Bonds retain 80% of value in default (mild event; feasibility ladder
        # removed — full-scale event solved directly in risk_branch).
        recovery_rate_D=0.80,   recovery_rate_F=0.80,

        # ── DEFAULT STATE (risk-channel branch) ──────────────────────────────
        # Canonical output cost of default (Arellano 2008 tradition): TFP in
        # the post-default branch is Z·(1 − cost·rho^h).  Without it (and
        # with the Bohn windfall removed via re-anchoring) default would be
        # expansionary — debt relief with no pain — and the risk premium
        # would have the wrong sign.  5% on impact, half-life ~7 quarters,
        # is conservative next to the Greek 2012 output collapse.
        def_output_cost_D=0.05, def_output_rho_D=0.90,
        # Pessimistic probability tilt for the risk weighting (EZ-lite dial;
        # 1.0 = off, physical probabilities — the Bocola-faithful baseline).
        chi_tilt=1.0,

        sdf_mode="income",
        kappa_d=2.00,

        # Outer Cole-Kehoe zone-indicator iteration (converges in 1 pass when
        # debt stays inside the crisis zone; damping 1.0 = undamped)
        ck_max_iter=25,
        ck_tol=1e-12,
        ck_damping=1.0,

        # ── FISCAL ────────────────────────────────────────────────────────────
        phi_lamb_D=0.15,        phi_lamb_F=0.15,
        G_D=0.0,                G_F=0.0,

        # ── TRADE / CES BASKET ───────────────────────────────────────────────
        omega_home=0.85,       epsilon_trade=0.5,

        # ── SOLVER SETTINGS ───────────────────────────────────────────────────
        T=500,
        tol_hh=1e-12,
        tol_dist=1e-12,
        tol_mkt=1e-12,

        # ── INITIAL GUESSES FOR STEADY-STATE SOLVER ───────────────────────────
        rk_D_guess=0.0045,      rk_F_guess=0.0045,
        beta_guess_D=0.98,      beta_guess_F=0.98,
    )
    return cal
