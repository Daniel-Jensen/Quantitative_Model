# FIRM BLOCK: COBB-DOUGLAS PRODUCTION WITH FULLY FLEXIBLE PRICES.


def markup_ss(cal, country="D"):
    # FLEXIBLE-PRICE REAL MARGINAL COST mc = (eps-1)/eps.
    return (cal[f"epsilon_{country}"] - 1) / cal[f"epsilon_{country}"]


def steady_state_firm(cal, Kap_ss, country="D"):
    # STEADY-STATE FIRM BLOCK (Y, w, mpk, I, C, AND THE GHH chi THAT PINS N_ss=1).
    alpha  = cal[f"alpha_{country}"]
    delta  = cal[f"delta_{country}"]
    Z_ss   = cal[f"Z_ss_{country}"]
    frisch = cal[f"frisch_{country}"]
    G      = cal[f"G_{country}"]

    mc_ss  = markup_ss(cal, country)
    N_ss   = 1.0
    Y_ss   = Z_ss * Kap_ss ** alpha * N_ss ** (1 - alpha)
    zeta_wc = cal[f"zeta_wc_{country}"]
    r_wc_ss = cal[f"r_dep_{country}_target"] + cal[f"credit_spread_target_{country}"]
    w_ss   = mc_ss * (1 - alpha) * Y_ss / N_ss / (1.0 + zeta_wc * r_wc_ss)
    mpk_ss = mc_ss * alpha * Y_ss / Kap_ss
    I_ss   = delta * Kap_ss
    C_ss   = Y_ss - I_ss - G
    chi    = w_ss / (N_ss ** (1 / frisch))   # GHH static FOC

    return dict(chi=chi, w_ss=w_ss, mpk_ss=mpk_ss,
                Y_ss=Y_ss, C_ss=C_ss, I_ss=I_ss)


def solve_firm_path(N_path, Kap_prod_path, Z_path, cal, country="D"):
    # FIRM QUANTITIES ALONG THE PATH: Y, THE FRICTIONLESS WAGE, AND mpk.
    alpha = cal[f"alpha_{country}"]
    mc    = markup_ss(cal, country)

    # Kap_prod_path[t] is the stock PRODUCING at t (predetermined, Bocola eq. 6),
    # so mpk_t is the marginal product of the vintage banks bought at t-1
    Y   = Z_path * Kap_prod_path ** alpha * N_path ** (1 - alpha)
    w   = mc * (1 - alpha) * Y / N_path
    mpk = mc * alpha * Y / Kap_prod_path

    return dict(Y=Y, w=w, mpk=mpk, mc=mc)
