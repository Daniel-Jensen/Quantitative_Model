# ALL MODEL PARAMETERS FOR THE TWO-COUNTRY HANK-GK MONETARY UNION.
# D-bonds are D-good claims (priced with rdep_D), F-bonds F-good claims
# (rdep_F); cross-border legs convert via p (D-goods per F-good).
# Bocola (2016) Tables 1-2 anchor every parameter with a direct counterpart.
# The deliberate divergences are flagged "vs Bocola" at the parameter itself —
# each one is load-bearing, so read the note before importing his value.


def get_calibration():
    # BUILD THE PARAMETER DICT CONSUMED BY EVERY BLOCK.
    cal = dict(
        # Household preferences. Bocola §II.A.1 uses log utility (NOT Epstein-Zin)
        # => sigma = 1. His nu = 0.5 (Table 1) is the INVERSE Frisch chosen for a
        # Frisch elasticity of 2; `frisch` here is the elasticity itself.
        # GHH (not his separable form) is his own §V.C open-economy fix.
        sigma_D=1.0,   sigma_F=1.0,
        frisch_D=2.0,  frisch_F=2.0,
        chi_D=0.5417,  chi_F=0.5417,   # warm start; SS solve overwrites to pin N_ss=1

        # Idiosyncratic income (Rouwenhorst)
        n_e_D=2,       n_e_F=2,
        rho_e_D=0.9,   sigma_e_D=0.2,
        rho_e_F=0.9,   sigma_e_F=0.2,

        # Asset grids
        a_min_D=0.0, a_max_D=87.2, n_a_D=250, a_curve_D=2.0,
        a_min_F=0.0, a_max_F=87.2, n_a_F=250, a_curve_F=2.0,

        # Firms (Cobb-Douglas, flexible prices)
        epsilon_D=6.0, epsilon_F=6.0,   # demand elasticity -> mc = (eps-1)/eps
        Z_ss_D=0.45,   Z_ss_F=0.45,     # warm start; SS solve overwrites to pin Y_ss=1

        # Capital (Jermann 1998 adjustment cost). ksi = elasticity of Tobin's q
        # wrt I/K (Bocola Table 2 posterior mean 0.42).
        alpha_D=0.30,  alpha_F=0.30,    # capital share (Bocola Table 1)
        delta_D=0.025, delta_F=0.025,
        ksi_D=0.42,    ksi_F=0.42,

        # Financial intermediary. f = exit/payout share, so the Omega kernel puts
        # weight (1-f) on the franchise value (Bocola's psi = 0.96 survival,
        # Table 2 posterior mean).
        f_D=0.04,               f_F=0.04,
        # R^bg = 1.003 quarterly (Bocola Table 1 sample-average risk-free rate).
        r_dep_D_target=0.003,   r_dep_F_target=0.003,
        # beta*R = 1 at the SS under log utility => beta_inter = 1/R^bg. This MUST
        # move with r_dep_target: leaving it at 0.99 collapses the alpha fixed
        # point to a near-tangency and the SS solve fails.
        beta_inter_D=0.997,     beta_inter_F=0.997,
        # lambda and omega_ent are SOLVED (calibrate_bank_targets) to hit these.
        # vs Bocola: leverage 4.0 not his lev^bg = 5.0, and spread 200bp not his
        # mu^bg = 0.001. The two are NOT independently choosable — the franchise
        # fixed point alpha = Omega(1+rdep)/(1-mu) has a FOLD, and at leverage 5
        # with a 200bp spread the calibration would need the UPPER root while
        # _alpha_ss_fixed_point (like value iteration from below) takes the least
        # root, so stage 1 cannot hit the target. Cutting the spread toward his
        # value instead parks the SS bank next to the mu=0 kink and the
        # EXPANSIONARY TFP experiment stops converging (75bp and 150bp both stall;
        # the risk experiment survives because sovereign risk TIGHTENS the IC).
        # Move either only with both experiments re-verified.
        leverage_target_D=4.0,          leverage_target_F=4.0,
        credit_spread_target_D=0.005,   credit_spread_target_F=0.005,  # 200 bps/yr
        # warm starts; overwritten by calibrate_bank_targets in steady_state.py
        lambda_K_D=0.22,        lambda_K_F=0.22,
        lambda_bD_D=0.22,       lambda_bD_F=0.22,
        lambda_bF_D=0.22,       lambda_bF_F=0.22,
        omega_ent_D=0.002,      omega_ent_F=0.002,

        # Cross-border bond portfolio adjustment costs
        psi_bF_D=0.05,          psi_bD_F=0.05,
        b_F_D_ss=0.744,         b_D_F_ss=0.744,   # ~20% of each supply (contagion leg)
        excess_return_F_D_ss=0.0,                 # overwritten after the SS solve
        excess_return_D_F_ss=0.0,                 # overwritten after the SS solve

        # Government bonds. vs Bocola: delta_b = 0.036 (~7y HM duration, the Greek
        # anchor) not his 0.056 (18-month Italian maturity). Long duration is what
        # makes priced risk generate large MTM losses; at delta_b = 0.25 the
        # repricing shrinks ~6x and the risk channel turns EXPANSIONARY.
        delta_b_D=0.036,        delta_b_F=0.036,
        B_gov_D_ss=3.722,       B_gov_F_ss=3.722,   # 93% of annual GDP

        # Default risk (Bocola 2016): the PRICED default probability pi_t is an
        # exogenous input path to the solver (his s-shock), built per experiment
        # in main.py. Only D is risky; the feared event is a pure haircut.
        recovery_rate_D=0.45,   # 55% haircut (Greek PSI 2012; Bocola D = 0.55)

        # Transmission Protection Instrument: Markov-switching CB backstop on the
        # D-sovereign, also built per experiment in main.py (off = no purchases,
        # backstop never doubted). psi_cb_D is the portfolio-balance elasticity
        # turning a Q_bD price gap into a purchase quantity. It plays the same role
        # as psi_bF_D/psi_bD_F but NOT at the same magnitude: those scale
        # return differentials (10s-100s of bps), this scales a PRICE-LEVEL gap
        # (several % of a ~0.75 bond price under stress). At 0.05 purchases reach
        # ~50% of B_gov_D_ss and the Newton solve diverges; at 0.5 they stay ~0.2%
        # and goods_D/goods_F land at the no-TPI tolerance. The single most
        # solver-sensitive parameter — do not lower without re-running
        # tests/test_tpi.py's budget-closure test.
        psi_cb_D=0.5,

        # Working capital (Neumeyer-Perri): firms pre-finance zeta x wage bill at
        # r_wc = rdep(-1) + lambda*mu/Omega. The only spread->output channel;
        # zeta = 0 nests it off exactly.
        zeta_wc_D=1.0,          zeta_wc_F=1.0,

        # Fiscal
        phi_lamb_D=0.15,        phi_lamb_F=0.15,   # Bohn rule strength
        G_D=0.0,                G_F=0.0,

        # Trade / CES basket
        omega_home=0.85,        epsilon_trade=0.5,

        # Solver settings
        T=300,                 # risk-shock horizon (T=100 truncates, T=500 identical)
        tol_hh=1e-12,
        tol_dist=1e-12,
        tol_mkt=1e-12,         # SS stage-1 hybr xtol
        tol_transition=1e-10,  # 7T acceptance; do NOT tighten (hybr plateaus ~5e-11)
        n_jobs=0,              # FD-Jacobian workers; 0 -> os.cpu_count()
        use_numba=True,        # JIT EGM/distribution kernels; numpy fallback otherwise
    )
    return cal
