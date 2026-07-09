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
        chi_D=1.0,     chi_F=1.0,      # initial guess; overwritten in SS solve to match Nss=1

        # ── IDIOSYNCRATIC INCOME PROCESS (Rouwenhorst) ───────────────────────
        n_e_D=2,     n_e_F=2,
        rho_e_D=0.9,   sigma_e_D=0.2,
        rho_e_F=0.9,   sigma_e_F=0.2,

        # ── ASSET GRIDS ──────────────────────────────────────────────────────
        a_min_D=0.0, a_max_D=87.2, n_a_D=250, a_curve_D=2.0,
        a_min_F=0.0, a_max_F=87.2, n_a_F=250, a_curve_F=2.0,

        # ── FIRMS (Cobb-Douglas + full price flexibility) ────────────────────
        epsilon_D=6.0, epsilon_F=6.0,   # demand elasticity → mc = (ε-1)/ε
        Z_ss_D=0.44804075,  Z_ss_F=0.44804075,  # scaled so Y_ss = 1 (Z = (1/Y_old)^(1-α))

        # ── CAPITAL (Jermann 1998 adjustment cost) ────────────────────────────
        alpha_D=0.35,  alpha_F=0.35,    # capital share
        delta_D=0.025, delta_F=0.025,
        ksi_D=0.50,    ksi_F=0.50,      # adjustment-cost curvature

        # ── FINANCIAL INTERMEDIARY ─────────────────────────────────────────────
        f_D=0.028,              f_F=0.028,
        r_dep_D_target=0.000,   r_dep_F_target=0.000,
        beta_inter_D=0.96,      beta_inter_F=0.96,
        lambda_K_D=0.22,        lambda_K_F=0.22,
        lambda_bD_D=0.22,       lambda_bD_F=0.22,
        lambda_bF_D=0.22,       lambda_bF_F=0.22,
        omega_ent_D=0.002,      omega_ent_F=0.002,

        # ── PORTFOLIO ADJUSTMENT COSTS (cross-border bonds) ──────────────────
        # b_D_F_ss / b_F_D_ss ≈ 20% of the respective bond supply: foreign
        # banks' pre-crisis holdings of peripheral debt (union contagion leg).
        psi_bF_D=0.01,          psi_bD_F=0.01,
        b_F_D_ss=0.744,         b_D_F_ss=0.744,
        excess_return_F_D_ss=0.0,               # overwritten after SS solve
        excess_return_D_F_ss=0.0,               # overwritten after SS solve

        # ── GOVERNMENT BONDS ─────────────────
        # delta_b = quarterly amortization rate
        delta_b_D=0.036,        delta_b_F=0.036,
        B_gov_D_ss=3.722,       B_gov_F_ss=3.722,   # 93% of annual GDP (= 0.93×4×Y_ss=1)

        # ── DEFAULT RISK (Cole-Kehoe zones × Bocola pricing) ─────────────────
        # Thresholds are FACE-value debt to quarterly Y_ss ratios.  F is always
        # safe.  SS sits inside the D crisis zone (b/Y_ss ≈ 3.7), so the
        # sunspot is priced immediately; b_ck_high is out of reach (no
        # fundamental default in the risk-only experiment).
        b_ck_low_D=3.00,        b_ck_low_F=99.0,
        b_ck_high_D=6.00,       b_ck_high_F=99.0,
        # Recovery = 1 − haircut; haircut 0.55 per Greek PSI 2012
        # (Zettelmeyer-Trebesch-Gulati 2013, used by Bocola 2016).
        recovery_rate_D=0.45,   recovery_rate_F=0.45,

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
        # Default-state SDF (risk channel):
        #   sdf_mode="empirical": Λ^d = beta_inter·kappa_d.  DEVIATION from
        #     Bocola (2016), whose Λ is the model-consistent (rep-agent,
        #     log-utility) household SDF.  Two reasons: (i) the comovement
        #     problem makes the model-consistent SDF wrong-signed here (see
        #     sdf_mode="model" below); (ii) in the HA economy there is no
        #     representative family, so "the household SDF" is not a
        #     well-defined object for the banker — Λ^d is an ASSUMPTION: a
        #     banker-specific discount with an empirically disciplined
        #     default-state loading kappa_d.  Calibration target:
        #     risk-channel share of the lending-spread response ≈ 45%.
        #   sdf_mode="model": Λ^d = beta_inter·(x^d/x^nd)^(−σ) from branch
        #     consumption composites.  Currently WRONG-SIGNED: the comovement
        #     problem (deposit-rate collapse) makes branch consumption RISE,
        #     so the household SDF prices default states as good times.  Use
        #     only after the union-deposit-market fix.
        sdf_mode="empirical",
        kappa_d=2.00,

        # Outer Cole-Kehoe zone-indicator iteration (converges in 1 pass when
        # debt stays inside the crisis zone; damping 1.0 = undamped)
        ck_max_iter=25,
        ck_tol=1e-11,
        ck_damping=1.0,

        # ── FISCAL ────────────────────────────────────────────────────────────
        phi_lamb_D=0.15,        phi_lamb_F=0.15,
        G_D=0.0,                G_F=0.0,

        # ── TRADE / CES BASKET ───────────────────────────────────────────────
        omega_home=0.85,       epsilon_trade=0.5,

        # ── SOLVER SETTINGS ───────────────────────────────────────────────────
        T=100,
        tol_hh=1e-11,
        tol_dist=1e-11,
        tol_mkt=1e-11,

        # ── INITIAL GUESSES FOR STEADY-STATE SOLVER ───────────────────────────
        rk_D_guess=0.0045,      rk_F_guess=0.0045,
        beta_guess_D=0.98,      beta_guess_F=0.98,
    )
    return cal
