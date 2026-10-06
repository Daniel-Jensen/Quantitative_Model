# ALL PARAMETERS OF THE GLOBAL MODEL; DEPARTURES FROM BOCOLA (2016) ARE FLAGGED WHERE THEY OCCUR.


def get_calibration():
    # BUILD THE PARAMETER DICT USED BY EVERY BLOCK.
    cal = dict(
        # households: log utility, Frisch elasticity 2 (Bocola), GHH preferences
        sigma_D=1.0,   sigma_F=1.0,
        frisch_D=2.0,  frisch_F=2.0,
        chi_D=0.5417,  chi_F=0.5417,  # overwritten by the SS solve to pin N_ss = 1

        # idiosyncratic income (Rouwenhorst)
        n_e_D=2,       n_e_F=2,
        rho_e_D=0.9,   sigma_e_D=0.2,
        rho_e_F=0.9,   sigma_e_F=0.2,

        # asset grids
        a_min_D=0.0, a_max_D=87.2, n_a_D=250, a_curve_D=2.0,
        a_min_F=0.0, a_max_F=87.2, n_a_F=250, a_curve_F=2.0,

        # firms: Cobb-Douglas, flexible prices
        epsilon_D=6.0, epsilon_F=6.0,  # demand elasticity
        Z_ss_D=0.45,   Z_ss_F=0.45,  # overwritten by the SS solve to pin Y_ss = 1

        # capital: Jermann adjustment costs; ksi = 0.50 is inside Bocola's posterior
        alpha_D=0.30,  alpha_F=0.30,  # capital share (Bocola)
        delta_D=0.025, delta_F=0.025,
        ksi_D=0.50,    ksi_F=0.50,

        # banker exit share: 0.08, vs Bocola's 0.0354 (the bank block's one departure)
        f_D=0.08,               f_F=0.08,
        # risk-free rate 1.003 a quarter (Bocola)
        r_dep_D_target=0.003,   r_dep_F_target=0.003,
        # beta*R = 1 at the SS, so beta_inter must move with r_dep_target
        beta_inter_D=0.997,     beta_inter_F=0.997,
        # leverage 5 (Bocola); spread 100 bp/yr, not his 8, so the constraint binds at rest
        leverage_target_D=5.0,            leverage_target_F=5.0,
        credit_spread_target_D=0.0025,    credit_spread_target_F=0.0025,  # 100 bp/yr
        # warm starts, overwritten by calibrate_bank_targets
        lambda_K_D=0.22,        lambda_K_F=0.22,
        lambda_bD_D=0.22,       lambda_bD_F=0.22,
        lambda_bF_D=0.22,       lambda_bF_F=0.22,
        omega_ent_D=0.002,      omega_ent_F=0.002,

        # Bocola's debt-elastic premium on the cross-border position; keeps the model stationary
        kappa_nfa=0.01,
        # cross-border adjustment costs; they pin each sovereign's split between the banks
        psi_bF_D=2.0,           psi_bD_F=2.0,
        b_F_D_ss=0.196,         b_D_F_ss=0.196,  # ~20% of each stock held abroad
        excess_return_F_D_ss=0.0,  # overwritten after the SS solve
        excess_return_D_F_ss=0.0,  # overwritten after the SS solve

        # bonds mature at 5.6% a quarter (Bocola's pi); long duration drives the repricing
        delta_b_D=0.056,        delta_b_F=0.056,
        # debt sized so sovereigns are 7.6% of D-bank assets (Bocola)
        B_gov_D_ss=0.98,        B_gov_F_ss=0.98,

        # default risk is an exogenous input (Bocola's s-shock); only D is risky
        recovery_rate_D=0.45,  # 55% haircut (Greek PSI; Bocola)
        # Bocola's 0.63 under his physicists' Gauss-Hermite rule, i.e. 0.63/sqrt(2) here
        sigma_s=0.4455,

        # TPI: the D share of the Eurosystem P&L (capital key, as in calibration/ssj.py)
        tpi_key_D=0.071,
        # TPI switch, spread cap (bp/yr over the F bond) and smoothing of its two kinks
        tpi_on=False,
        tpi_cap_bp=200.0,
        tpi_eps=1e-4,
        # the form the TPI rules are read in (the solver roots in the FB form)
        tpi_gz=True,

        # working capital: firms pre-finance the wage bill at r_wc; zeta = 0 switches it off
        zeta_wc_D=1.0,          zeta_wc_F=1.0,

        # the Bohn tax rule, stated as the debt root it delivers
        debt_root_D=0.93,       debt_root_F=0.93,
        G_D=0.0,                G_F=0.0,

        # F is eight times D's size; every variable is per capita of its own country
        size_D=1.0,             size_F=8.0,

        # home bias; omega_home_F is derived below so trade balances at p = 1
        omega_home_D=0.85,      epsilon_trade=0.5,

        # solver settings
        T=300,  # risk-shock horizon
        tol_hh=1e-12,
        tol_dist=1e-12,
        tol_mkt=1e-12,  # SS stage-1 tolerance
        tol_transition=1e-10,  # do not tighten: hybr plateaus ~5e-11
        n_jobs=0,  # Jacobian workers: 0 = every core, 1 = serial
        use_numba=True,  # numba kernels; numpy fallback
    )
    # derived: F's home bias consistent with the country sizes
    cal["omega_home_F"] = 1.0 - (1.0 - cal["omega_home_D"]) * cal["size_D"] / cal["size_F"]
    cal["omega_home"] = cal["omega_home_D"]  # legacy key
    return cal
