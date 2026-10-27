# ALL PARAMETERS OF THE GLOBAL MODEL, IN TWO CALIBRATIONS.
# "paper" is the viva paper's (sections 2 and 4, appendix B; Overleaf, 2026-10-22). It is
# asymmetric across countries (G, coupon decay, own-sovereign shares), and the steady state as
# written imposes one union deposit market and a symmetric steady state, so it cannot be solved
# until deposits clear nationally (TODO, Code). "bocola" is the Bocola-based calibration behind
# the 6 Oct 2026 reference solution (STATE section 3), kept so that run reproduces.

CALIBRATION = "paper"  # the active calibration: "paper" or "bocola"


def get_calibration(version=None):
    # BUILD THE PARAMETER DICT USED BY EVERY BLOCK.
    version = CALIBRATION if version is None else version
    cal = _bocola()
    if version == "paper":
        for c in ("D", "F"):
            del cal[f"debt_root_{c}"]  # the paper states the fiscal rule by its Bohn coefficient
        cal.update(_paper())
        _size_sovereign_stocks(cal)
    elif version != "bocola":
        raise ValueError(f"unknown calibration {version!r}: 'paper' or 'bocola'")
    # derived: F's home bias consistent with the country sizes
    cal["omega_home_F"] = 1.0 - (1.0 - cal["omega_home_D"]) * cal["size_D"] / cal["size_F"]
    cal["omega_home"] = cal["omega_home_D"]  # legacy key
    return cal


def _paper():
    # THE PAPER'S VALUES WHERE THEY DIFFER FROM THE BOCOLA CALIBRATION (APPENDIX B, TABLE 1).
    return dict(
        # households: Smets-Wouters euro-area risk aversion; Frisch 0.5 (their 1/2.4 = 0.42)
        sigma_D=1.35,  sigma_F=1.35,
        frisch_D=0.5,  frisch_F=0.5,

        # retail demand elasticity 11: a 10% steady-state markup
        epsilon_D=11.0, epsilon_F=11.0,

        # banker discount factor 1/(1 + r) exactly
        beta_inter_D=1.0 / 1.003,  beta_inter_F=1.0 / 1.003,
        # leverage 8 (Greek financial accounts, 7.98); wedge 5 bp/yr (Euribor-Eurepo,
        # 2002Q2-2009Q4), so mu_bar = 0.0010
        leverage_target_D=8.0,            leverage_target_F=8.0,
        credit_spread_target_D=0.000125,  credit_spread_target_F=0.000125,

        # sovereign stocks are sized by _size_sovereign_stocks to these bank-asset shares
        # (EBA 2011, gross direct long positions, end-2010)
        own_sov_share_D=0.144,  # own sovereign / assets, six Greek banks
        own_sov_share_F=0.065,  # own sovereign / assets, twelve German banks
        foreign_share_D=0.127,  # German banks' share of Greek plus German banks' Greek holdings
        # Bunds / assets of the six Greek banks (EUR 0.41bn / 377bn); not in the paper, which
        # leaves D banks' holdings of F bonds open
        for_sov_share_D=0.0011,

        # coupon decay from average maturity (Eurostat); D's gives a 4.2-year duration
        delta_b_D=0.056,        delta_b_F=0.039,

        # risk factor: persistence 0.9 (tranquil-years estimate, bias-corrected); innovation
        # s.d. 0.41 (ten-year spreads rescaled by the CDS/bond volatility ratio)
        rho_s=0.9,
        sigma_s=0.41,

        # government purchases / GDP (Eurostat); annual Bohn coefficient 0.04 (Mendoza-Ostry),
        # a quarterly tax response gamma_tau = 0.01 per unit of debt above target
        G_D=0.204,              G_F=0.191,
        bohn_D=0.04,            bohn_F=0.04,

        # F/D nominal GDP 11.3; Greek import share of consumption 32% (2012); trade
        # elasticity 1.5 (Fitzgerald and Haller)
        size_D=1.0,             size_F=11.3,
        omega_home_D=0.68,      epsilon_trade=1.5,
    )


def _size_sovereign_stocks(cal):
    # SOVEREIGN STOCKS AND CROSS-BORDER BOOKS THAT DELIVER THE TARGET BANK-ASSET SHARES.
    # In the deterministic steady state every bank asset earns rk = r + wedge, so with Y = N = 1
    # capital is K = mc*alpha/(rk + delta), the working-capital loan is zeta*w and a perpetuity
    # costs delta_b/(rk + delta_b). Each bank's non-bond assets are then known, and each
    # share target is linear in the stocks. Exact at p = 1, which national clearing moves.
    def non_bond_assets(c):
        # CAPITAL PLUS THE WORKING-CAPITAL LOAN, PER CAPITA OF COUNTRY c.
        mc = (cal[f"epsilon_{c}"] - 1.0) / cal[f"epsilon_{c}"]
        rk = cal[f"r_dep_{c}_target"] + cal[f"credit_spread_target_{c}"]
        K = mc * cal[f"alpha_{c}"] / (rk + cal[f"delta_{c}"])
        w = mc * (1.0 - cal[f"alpha_{c}"]) / (1.0 + cal[f"zeta_wc_{c}"] * rk)
        return K + cal[f"zeta_wc_{c}"] * w

    def price(c, bond):
        # THE PERPETUITY PRICE OF bond TO A COUNTRY-c BANK.
        rk = cal[f"r_dep_{c}_target"] + cal[f"credit_spread_target_{c}"]
        return cal[f"delta_b_{bond}"] / (rk + cal[f"delta_b_{bond}"])

    sz = cal["size_F"] / cal["size_D"]
    # D bank: own D bonds and F bonds, per D capita
    assets_D = non_bond_assets("D") / (1.0 - cal["own_sov_share_D"] - cal["for_sov_share_D"])
    B_D = cal["own_sov_share_D"] * assets_D / (price("D", "D") * (1.0 - cal["foreign_share_D"]))
    b_F_D = cal["for_sov_share_D"] * assets_D / price("D", "F")
    # F bank: own F bonds and its slice of the D stock, per F capita
    held_D = price("F", "D") * cal["foreign_share_D"] * B_D / sz
    assets_F = (non_bond_assets("F") + held_D) / (1.0 - cal["own_sov_share_F"])
    b_F_F = cal["own_sov_share_F"] * assets_F / price("F", "F")
    cal.update(B_gov_D_ss=B_D, B_gov_F_ss=b_F_F + b_F_D / sz,
               b_D_F_ss=cal["foreign_share_D"] * B_D, b_F_D_ss=b_F_D)


def _bocola():
    # THE BOCOLA-BASED CALIBRATION OF THE 6 OCT 2026 REFERENCE RUN.
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
        # the risk factor s: p^d = 0.1% at rest, persistence 0.95 (Bocola's estimate)
        pd_rest=0.001,
        rho_s=0.95,
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

        # home bias; omega_home_F is derived in get_calibration so trade balances at p = 1
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
    return cal
