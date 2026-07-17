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

        # Financial intermediary. f = exit/payout share; kernel weights (1-f) on
        # the franchise value α′ (Bocola's ψ = 1-f = survival). beta_inter proxies
        # the banker's SDF Λ^nd (Bocola: household MRS, Λ_ss = β_hh ≈ 0.99);
        # values ≪ 1/(1+rdep) push α_ss below 1 and mute the franchise channel.
        f_D=0.05,              f_F=0.05,
        r_dep_D_target=0.000,   r_dep_F_target=0.000,
        beta_inter_D=0.99,      beta_inter_F=0.99,
        # λ and ω_ent are SOLVED (calibrate_bank_targets) to hit these SS moments
        leverage_target_D=4.0,          leverage_target_F=4.0,
        credit_spread_target_D=0.005,   credit_spread_target_F=0.005,  # 200 bps/yr
        # warm starts; overwritten by calibrate_bank_targets in steady_state.py
        lambda_K_D=0.22,        lambda_K_F=0.22,
        lambda_bD_D=0.22,       lambda_bD_F=0.22,
        lambda_bF_D=0.22,       lambda_bF_F=0.22,
        omega_ent_D=0.002,      omega_ent_F=0.002,

        # Cross-border bond portfolio adjustment costs
        psi_bF_D=0.05,          psi_bD_F=0.05,
        b_F_D_ss=0.744,         b_D_F_ss=0.744,   # ≈20% of each bond supply (contagion leg)
        excess_return_F_D_ss=0.0,               # overwritten after SS solve
        excess_return_D_F_ss=0.0,               # overwritten after SS solve

        # Government bonds. δ_b=0.036 ⇒ HM duration ~7y (long duration → large MTM losses)
        delta_b_D=0.036,        delta_b_F=0.036,
        B_gov_D_ss=3.722,       B_gov_F_ss=3.722,   # 93% of annual GDP

        # Default risk (Bocola 2016): the PRICED default probability π_t is an
        # exogenous input path to the solver (his s-shock, eqs. 11-12). Only D
        # is risky; F bonds are safe. The feared event is a pure haircut.
        recovery_rate_D=0.45,   # 55% haircut (Greek PSI 2012; Bocola D=0.55)

        # Default-branch DIAGNOSTIC flags — all OFF (=0) at the Bocola-pure
        # baseline; the default-state recession must arise endogenously.
        def_output_cost_D=0.0,  def_output_rho_D=0.90,   # Arellano output cost
        def_capital_quality_D=0.0,   # GK ξ_K capital-quality loss at h=0
        recap_share_D=0.0,           # HFSF-style recap financed by issuance

        # Working capital (Neumeyer-Perri): firms pre-finance ζ×wage-bill at
        # r_wc=rdep(-1)+λμ/Ω̃ (Bocola §V.C open-economy variant, his R_W).
        # The spread→output channel; ζ=0 nests the model off.
        zeta_wc_D=1.0,          zeta_wc_F=1.0,

        # Fiscal
        phi_lamb_D=0.15,        phi_lamb_F=0.15,   # Bohn rule strength
        G_D=0.0,                G_F=0.0,

        # Trade / CES basket
        omega_home=0.85,       epsilon_trade=0.5,

        # Solver settings
        T=200,                 # risk-shock horizon (T=100 truncates, T=500 identical at 5× cost)
        tol_hh=1e-12,
        tol_dist=1e-12,
        tol_mkt=1e-12,        # SS stage-1 hybr xtol
        tol_transition=1e-10,  # 7T acceptance; do NOT tighten to 1e-12 (hybr plateaus ~5e-11)
        n_jobs=0,             # FD-Jacobian workers; 0 → os.cpu_count()
        use_numba=True,       # JIT EGM/distribution kernels; numpy fallback otherwise
    )
    return cal
