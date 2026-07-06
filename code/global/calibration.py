"""
Notes: 
* Shared bond denomination: all bonds are D-good claims (D is numeraire).
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
        a_min_D=0.0, a_max_D=150.0, n_a_D=250, a_curve_D=2.0,
        a_min_F=0.0, a_max_F=150.0, n_a_F=250, a_curve_F=2.0,

        # ── FIRMS (Cobb-Douglas + full price flexibility) ────────────────────
        epsilon_D=6.0, epsilon_F=6.0,   # demand elasticity → mc = (ε-1)/ε
        Z_ss_D=1.0,    Z_ss_F=1.0,

        # ── CAPITAL (Jermann 1998 adjustment cost) ────────────────────────────
        alpha_D=0.35,  alpha_F=0.35,    # capital share
        delta_D=0.025, delta_F=0.025,   # initial guess for SS solve; overwritten after SS solve
        ksi_D=0.50,    ksi_F=0.50,      # adjustment-cost curvature

        # ── FINANCIAL INTERMEDIARY ─────────────────────────────────────────────
        f_D=0.028,              f_F=0.028,
        r_dep_D_target=0.000,   r_dep_F_target=0.000,
        beta_inter_D=0.96,      beta_inter_F=0.96,
        lambda_K_D=0.30,        lambda_K_F=0.30,
        lambda_bD_D=0.04,       lambda_bD_F=0.04,
        lambda_bF_D=0.04,       lambda_bF_F=0.04,
        omega_ent_D=0.002,      omega_ent_F=0.002,

        # ── PORTFOLIO ADJUSTMENT COSTS (cross-border bonds) ──────────────────
        psi_bF_D=0.01,          psi_bD_F=0.01,
        b_F_D_ss=0.005,         b_D_F_ss=0.005,
        excess_return_F_D_ss=0.0,               # overwritten after SS solve
        excess_return_D_F_ss=0.0,               # overwritten after SS solve

        # ── GOVERNMENT BONDS ─────────────────
        delta_b_D=0.10,         delta_b_F=0.10,
        B_gov_D_ss=2.40,        B_gov_F_ss=2.40,

        # ── DEFAULT RISK ───────────────────────────────────────────────────────
        # Thresholds are debt-to-Y_ss ratios.  F is always safe.
        # b_ck_low/high are used by Cole-Kehoe crisis-zone logic (see solve_transition_ck).
        # b_ck_high also serves as the Bocola-Dovis fundamental default boundary (B̄).
        b_ck_low_D=0.55,        b_ck_low_F=99.0,
        b_ck_high_D=1.20,       b_ck_high_F=99.0,
        recovery_rate_D=0.40,   recovery_rate_F=0.0,   # Greek PSI-style haircut

        # Bocola-Dovis (2019): sunspot tightens GK IC for sovereign bonds
        # lbD_D_eff = lbD_D + psi_bd_D * xi_{t+1}  (D-bank and F-bank for D-bonds)
        # lbF_F_eff = lbF_F + psi_bd_F * xi_{t+1}  (F-bonds; set 0 — no F default)
        psi_bd_D=3.0,           psi_bd_F=0.0,

        # Outer CK fixed-point solver
        ck_max_iter=25,
        ck_tol=1e-5,
        ck_damping=0.5,

        # Outer BD fixed-point solver (separate budget; uses Anderson acceleration)
        bd_max_iter=50,
        bd_tol=1e-4,        # 0.004% of b_ss — adequate for quantitative IRF
        bd_anderson_m=3,    # Anderson window size

        # ── FISCAL ────────────────────────────────────────────────────────────
        phi_lamb_D=0.15,        phi_lamb_F=0.15,
        G_D=0.0,                G_F=0.0,

        # ── TRADE / CES BASKET ───────────────────────────────────────────────
        omega_home=0.85,       epsilon_trade=0.5,

        # ── SOLVER SETTINGS ───────────────────────────────────────────────────
        T=100,
        tol_hh=1e-9,
        tol_dist=1e-9,
        tol_mkt=1e-9,

        # ── INITIAL GUESSES FOR STEADY-STATE SOLVER ───────────────────────────
        rk_D_guess=0.010,       rk_F_guess=0.010,
        beta_guess_D=0.98,      beta_guess_F=0.98,
    )
    return cal
