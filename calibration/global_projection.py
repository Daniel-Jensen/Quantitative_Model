# ALL PARAMETERS OF THE GLOBAL MODEL; DEPARTURES FROM BOCOLA (2016) ARE FLAGGED WHERE THEY OCCUR.


def get_calibration():
    # BUILD THE PARAMETER DICT USED BY EVERY BLOCK.
    cal = dict(
        # households, GHH (appendix B, 2026-10-26): risk aversion 1.35 and Frisch 0.5 (Smets-Wouters 2003)
        sigma_D=1.35,  sigma_F=1.35,
        frisch_D=0.5,  frisch_F=0.5,
        chi_D=0.5417,  chi_F=0.5417,  # overwritten by the SS solve to pin N_ss = 1

        # idiosyncratic income (Rouwenhorst)
        n_e_D=2,       n_e_F=2,
        rho_e_D=0.9,   sigma_e_D=0.2,
        rho_e_F=0.9,   sigma_e_F=0.2,

        # asset grids
        a_min_D=0.0, a_max_D=87.2, n_a_D=250, a_curve_D=2.0,
        a_min_F=0.0, a_max_F=87.2, n_a_F=250, a_curve_F=2.0,

        # firms: Cobb-Douglas, flexible prices
        epsilon_D=11.0, epsilon_F=11.0,  # demand elasticity: a 10% SS markup
        Z_ss_D=0.45,   Z_ss_F=0.45,  # overwritten by the SS solve to pin Y_ss = 1

        # capital: Jermann adjustment costs; ksi = 0.50 is inside Bocola's posterior
        alpha_D=0.30,  alpha_F=0.30,  # capital share (Bocola)
        delta_D=0.025, delta_F=0.025,
        ksi_D=0.50,    ksi_F=0.50,

        # banker exit share: 0.08, vs Bocola's 0.0354 (the bank block's one departure)
        f_D=0.08,               f_F=0.08,
        # risk-free rate 1.003 a quarter (Bocola); beta_inter = 1/(1 + r) is set below the dict
        r_dep_D_target=0.003,   r_dep_F_target=0.003,
        # THE TARGETS, hit at the STOCHASTIC REST POINT (calibrate_stochastic.py): leverage 8 (GR
        # financial accounts 2002-07: 7.98; not Bocola's 5), the loan spread r_wc - r_dep at
        # 5 bp/yr (appendix B: the 2002Q2-2009Q4 Euribor-Eurepo spread; it was 100 bp until
        # 2026-10-26), and own-sovereign bonds as a share of own-bank assets (EBA 2010: GR 14.4%,
        # DE 6.5%; not Bocola's 7.6%)
        srp_leverage_D=8.0,     srp_leverage_F=8.0,
        srp_spread_D=0.000125,  srp_spread_F=0.000125,
        srp_exposure_D=0.144,   srp_exposure_F=0.065,
        # the deterministic-SS instruments that deliver them (the appendix's closed forms until the
        # rest-point loop has run)
        leverage_target_D=8.0,            leverage_target_F=8.0,
        credit_spread_target_D=0.000125,  credit_spread_target_F=0.000125,
        sov_exposure_D=0.144,             sov_exposure_F=0.065,
        # warm starts, overwritten by calibrate_bank_targets
        lambda_K_D=0.22,        lambda_K_F=0.22,
        lambda_bD_D=0.22,       lambda_bD_F=0.22,
        lambda_bF_D=0.22,       lambda_bF_F=0.22,
        omega_ent_D=0.002,      omega_ent_F=0.002,

        # Bocola's debt-elastic premium on the cross-border position; keeps the model stationary
        kappa_nfa=0.01,
        # cross-border adjustment costs; they pin each sovereign's split between the banks
        psi_bF_D=2.0,           psi_bD_F=2.0,
        # the share of each stock held by the OTHER country's banks (EBA 2010): F banks hold 12.7%
        # of the D stock (German banks' share of Greek and German banks' holdings; appendix B), D
        # banks 0.5% of the F stock (not in the appendix; ~0 in the data, positive so the book stays
        # clear of the solver's floors; kept by decision 2026-10-26)
        foreign_share_D=0.127,  foreign_share_F=0.005,
        # B_gov_D/F_ss, b_D_F_ss, b_F_D_ss and the friction anchors are SOLVED in the SS
        # (steady_state._size_books) from the exposure instruments and these shares
        # no short sales on the cross-border books: each book >= 0 against its FOC <= 0, a
        # Fischer-Burmeister pair, HARD (2026-10-26): the solve ladder works at 1e-4 and polishes
        # the exact pair in at the end (recursive_experiment.HARD_WORK)
        book_eps=0.0,

        # coupon decay (Eurostat irt, appendix B): D 5.6% a quarter, Macaulay duration 4.2 years at
        # the 1.25% SS yield (Greek banks' Greek book, EBA 2010: 4.4 years); F 3.9%
        delta_b_D=0.056,        delta_b_F=0.039,

        # default risk is an exogenous input (Bocola's s-shock); only D is risky
        recovery_rate_D=0.45,  # 55% haircut (Greek PSI; Bocola)
        # the logit-p^d AR(1) (appendix B, risk-process table): persistence 0.9 (the bias-corrected
        # 2001-07 ten-year estimate), innovation s.d. 0.41 (full-sample ten-year, rescaled by the
        # one-year CDS ratio); the probabilists' Gauss-Hermite rule takes it as is
        rho_s=0.9,              sigma_s=0.41,

        # the leverage constraint's KKT switch mu = max(z, 0) is a GENUINE max (2026-10-26 decision):
        # the solve ladder walks the smoothing (z + sqrt(z^2 + eps^2))/2 down from 1e-2 to 1e-4 and
        # polishes the hard max in at the end (recursive_experiment.HARD_WORK); 0 here also gives
        # the appendix's closed-form SS (bank.ic_slack_factor = 0)
        mu_eps=0.0,
        # the funded net-worth floor's max, HARD as well (solved at 1e-3 first)
        floor_eps=0.0,
        # the priced default probability at rest, per quarter (the AR(1) mean of s): 0.10%
        # (appendix B; GR 1y CDS 2008Q3: 0.092%)
        pd_rest=0.001,

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

        # the Bohn tax rule: annual response of the primary surplus to debt, both over annual output.
        # 0.06, the Greek 2001-2019 estimate in appendix B's fiscal table (decision 2026-10-26; the
        # appendix's Mendoza-Ostry 0.04 gave root 0.993 and an explosive rule below q = 0.821, and broke
        # the coarse d0 solve alone); root 0.987 at the SS bond price 0.947, explosive below q = 0.732
        bohn_annual=0.06,
        # government purchases over output (Eurostat nama_10_gdp; appendix B)
        G_D=0.204,              G_F=0.191,

        # F is 11.3 times D's size (DE/GR nominal GDP); every variable is per capita of its own country
        size_D=1.0,             size_F=11.3,

        # home bias (Greek import share 32%, 2012) and trade elasticity 1.5 (Fitzgerald-Haller 2014);
        # omega_home_F is derived below so trade balances at p = 1
        omega_home_D=0.68,      epsilon_trade=1.5,

        # solver settings
        T=300,  # risk-shock horizon
        tol_hh=1e-12,
        tol_dist=1e-12,
        tol_mkt=1e-12,  # SS stage-1 tolerance
        tol_transition=1e-10,  # do not tighten: hybr plateaus ~5e-11
        n_jobs=0,  # Jacobian workers: 0 = every core, 1 = serial
        use_numba=True,  # numba kernels; numpy fallback
    )
    # derived: the banker discount factor with beta_inter*(1 + r) = 1 exactly (appendix B)
    for c in ("D", "F"):
        cal[f"beta_inter_{c}"] = 1.0 / (1.0 + cal[f"r_dep_{c}_target"])
    # derived: F's home bias consistent with the country sizes
    cal["omega_home_F"] = 1.0 - (1.0 - cal["omega_home_D"]) * cal["size_D"] / cal["size_F"]
    cal["omega_home"] = cal["omega_home_D"]  # legacy key
    return cal
