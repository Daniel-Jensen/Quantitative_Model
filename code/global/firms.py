"Firm block: Cobb-Douglas production with full price flexibility."


def markup_ss(cal, country="D"):
    #CALCULATES THE STEADY STATE MARKUP (FLEX PRICES)
    return (cal[f"epsilon_{country}"] - 1) / cal[f"epsilon_{country}"]


def steady_state_firm(cal, Kap_ss, country="D"):
    # EASY PEASY COBB DOUGLAS PRODUCTION FUNCTION
    
    alpha  = cal[f"alpha_{country}"]
    delta  = cal[f"delta_{country}"]
    Z_ss   = cal[f"Z_ss_{country}"]
    frisch = cal[f"frisch_{country}"]
    G      = cal[f"G_{country}"]

    mc_ss  = markup_ss(cal, country)
    N_ss   = 1.0
    Y_ss   = Z_ss * Kap_ss ** alpha * N_ss ** (1 - alpha)
    # Working-capital wedge (Neumeyer-Perri): firms pre-finance zeta_wc of
    # the wage bill at r_wc = rdep + IC wedge.  At SS the single-λ wedge
    # equals the credit-spread target by construction (calibrate_bank_targets),
    # so r_wc_ss is a calibration constant — no bank fixed point needed here.
    zeta_wc  = cal.get(f"zeta_wc_{country}", 0.0)
    r_wc_ss  = (cal[f"r_dep_{country}_target"]
                + cal[f"credit_spread_target_{country}"])
    w_ss   = mc_ss * (1 - alpha) * Y_ss / N_ss / (1.0 + zeta_wc * r_wc_ss)
    mpk_ss = mc_ss * alpha * Y_ss / Kap_ss
    I_ss   = delta * Kap_ss
    C_ss   = Y_ss - I_ss - G

    # GHH static FOC at SS.
    chi = w_ss / (N_ss ** (1 / frisch))

    return dict(chi=chi, mc_ss=mc_ss, w_ss=w_ss, mpk_ss=mpk_ss,
                N_ss=N_ss, Y_ss=Y_ss, C_ss=C_ss, I_ss=I_ss)


def solve_firm_path(N_path, Kap_path, Z_path, cal, country="D"):
    # WE CAN SOLVE THE ENTIRE FIRM BLOCK IN ONE GO, GIVEN THE PATHS OF N, K, AND Z.

    """Firm quantities along a transition path (purely contemporaneous).

    Returns dict: Y, w, mpk, mc — all shape (T,).
    """
    alpha = cal[f"alpha_{country}"]
    mc    = markup_ss(cal, country)

    Y   = Z_path * Kap_path ** alpha * N_path ** (1 - alpha)
    w   = mc * (1 - alpha) * Y / N_path
    mpk = mc * alpha * Y / Kap_path

    return dict(Y=Y, w=w, mpk=mpk, mc=mc)
