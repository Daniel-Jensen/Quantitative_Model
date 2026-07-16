# **All model parameters for the two-country HANK–GK monetary union.**
# Bond denomination: D-bonds are D-good claims (priced with rdep_D), F-bonds
# F-good claims (rdep_F); cross-border legs convert via p (D-goods per F-good).


def get_calibration():
    cal = dict(
        # Household preferences (GHH)
        sigma_D=2.0,   sigma_F=2.0,
        frisch_D=0.5,  frisch_F=0.5,
        chi_D=0.5417,  chi_F=0.5417,   # warm start; SS solve overwrites to pin N_ss=1

        # Idiosyncratic income (Rouwenhorst)
        n_e_D=2,     n_e_F=2,
        rho_e_D=0.9,   sigma_e_D=0.2,
        rho_e_F=0.9,   sigma_e_F=0.2,

        # Asset grids
        a_min_D=0.0, a_max_D=87.2, n_a_D=250, a_curve_D=2.0,
        a_min_F=0.0, a_max_F=87.2, n_a_F=250, a_curve_F=2.0,

        # Firms (Cobb-Douglas, flexible prices)
        epsilon_D=6.0, epsilon_F=6.0,   # demand elasticity → mc = (ε-1)/ε
        Z_ss_D=0.45,   Z_ss_F=0.45,     # warm start; SS solve overwrites to pin Y_ss=1

        # Capital (Jermann 1998 adjustment cost)
        alpha_D=0.35,  alpha_F=0.35,    # capital share
        delta_D=0.025, delta_F=0.025,
        ksi_D=0.50,    ksi_F=0.50,      # adjustment-cost curvature

        # Financial intermediary
        f_D=0.05,              f_F=0.05,
        r_dep_D_target=0.000,   r_dep_F_target=0.000,
        beta_inter_D=0.96,      beta_inter_F=0.96,
        # λ and ω_ent are SOLVED (calibrate_bank_targets) to hit these SS moments
        leverage_target_D=4.0,          leverage_target_F=4.0,
        credit_spread_target_D=0.005,   credit_spread_target_F=0.005,  # 200 bps/yr
        # warm starts; overwritten by calibrate_bank_targets in steady_state.py
        lambda_K_D=0.22,        lambda_K_F=0.22,
        lambda_bD_D=0.22,       lambda_bD_F=0.22,
        lambda_bF_D=0.22,       lambda_bF_F=0.22,
        omega_ent_D=0.002,      omega_ent_F=0.002,
        entrant_mode="proportional",   # "proportional" (unit root) | "anchored"
        phi_entry=0.0,                 # anchored mode only: entrant scale on α^phi

        # Cross-border bond portfolio adjustment costs
        psi_bF_D=0.05,          psi_bD_F=0.05,
        b_F_D_ss=0.744,         b_D_F_ss=0.744,   # ≈20% of each bond supply (contagion leg)
        excess_return_F_D_ss=0.0,               # overwritten after SS solve
        excess_return_D_F_ss=0.0,               # overwritten after SS solve

        # Government bonds. δ_b=0.036 ⇒ HM duration ~7y (long duration → large MTM losses)
        delta_b_D=0.036,        delta_b_F=0.036,
        B_gov_D_ss=3.722,       B_gov_F_ss=3.722,   # 93% of annual GDP

        # Default risk (Cole-Kehoe zones × Bocola pricing); thresholds are b/Y_ss.
        # SS sits in the D crisis zone; F always safe; b_ck_high out of reach.
        b_ck_low_D=3.00,        b_ck_low_F=99.0,
        b_ck_high_D=6.00,       b_ck_high_F=99.0,
        recovery_rate_D=0.45,   recovery_rate_F=0.45,   # 55% haircut (Greek PSI 2012)

        # Default-state branch. Output cost (Arellano 2008) makes the state a
        # recession → correct risk-premium sign; without it default is expansionary.
        def_output_cost_D=0.05, def_output_rho_D=0.90,
        # GK capital-quality loss ξ_K at h=0: stops capital being the branch safe
        # haven (else two-branch pricing drives μ<0 and the channel turns expansionary).
        def_capital_quality_D=0.05,
        # Government recap (HFSF/EFSF): equity injection financed by issuance;
        # makes the full PSI haircut feasible. Always on when > 0.
        recap_share_D=0.5,
        # Branch prices ONE fixed event: branch_haircut_scale of the (1-recovery)
        # haircut (1.0 = full PSI). Ladder = opt-in scale search if infeasible.
        branch_haircut_scale=1.0,
        branch_use_ladder=False,
        chi_tilt=1.0,           # pessimistic probability tilt (1.0 = off, Bocola baseline)
        sdf_mode="income",      # default-state SDF: "income" | "empirical" | "model"
        kappa_d=2.00,           # free SDF loading for sdf_mode="empirical"

        # Working capital (Neumeyer-Perri): firms pre-finance ζ×wage-bill at
        # r_wc=rdep(-1)+λμ/Ω̃. The spread→output channel; ζ=0 nests the model off.
        zeta_wc_D=1.0,          zeta_wc_F=1.0,

        # Cole-Kehoe outer zone-indicator iteration
        ck_max_iter=25,
        ck_tol=1e-12,
        ck_damping=1.0,         # 1.0 = undamped

        # Fiscal
        phi_lamb_D=0.15,        phi_lamb_F=0.15,   # Bohn rule strength
        G_D=0.0,                G_F=0.0,

        # Trade / CES basket
        omega_home=0.85,       epsilon_trade=0.5,

        # Solver settings
        T=200,                 # sunspot horizon (T=100 truncates, T=500 identical at 5× cost)
        tol_hh=1e-12,
        tol_dist=1e-12,
        tol_mkt=1e-12,        # SS stage-1 hybr xtol
        tol_transition=1e-10,  # 7T acceptance; do NOT tighten to 1e-12 (hybr plateaus ~5e-11)
        n_jobs=0,             # FD-Jacobian workers; 0 → os.cpu_count()
        use_numba=True,       # JIT EGM/distribution kernels; numpy fallback otherwise
    )
    return cal
