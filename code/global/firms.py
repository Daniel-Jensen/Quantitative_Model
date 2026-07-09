"Firm block: Cobb-Douglas production with full price flexibility."


def markup_ss(cal, country="D"):
    "Real marginal cost (inverse markup). Constant under full price flexibility."
    return (cal[f"epsilon_{country}"] - 1) / cal[f"epsilon_{country}"]


def steady_state_firm(cal, Kap_ss, country="D"):
    """Steady-state firm block. N_ss = 1 by normalisation; back-solves chi
    from the GHH labour-supply FOC chi·N^(1/frisch) = w/P_CES (P_CES=1 at SS).

    Returns: chi, mc_ss, w_ss, mpk_ss, N_ss, Y_ss, C_ss, I_ss.
    C_ss is residual Y − I − G (goods-market identity at SS).
    """
    alpha  = cal[f"alpha_{country}"]
    delta  = cal[f"delta_{country}"]
    Z_ss   = cal[f"Z_ss_{country}"]
    frisch = cal[f"frisch_{country}"]
    G      = cal[f"G_{country}"]

    mc_ss  = markup_ss(cal, country)
    N_ss   = 1.0
    Y_ss   = Z_ss * Kap_ss ** alpha * N_ss ** (1 - alpha)
    w_ss   = mc_ss * (1 - alpha) * Y_ss / N_ss
    mpk_ss = mc_ss * alpha * Y_ss / Kap_ss
    I_ss   = delta * Kap_ss
    C_ss   = Y_ss - I_ss - G

    # GHH static FOC at SS (P_CES normalised to 1 at symmetric SS):
    #   chi * N_ss^(1/frisch) = w_ss  →  chi = w_ss / N_ss^(1/frisch)
    chi = w_ss / (N_ss ** (1 / frisch))

    return dict(chi=chi, mc_ss=mc_ss, w_ss=w_ss, mpk_ss=mpk_ss,
                N_ss=N_ss, Y_ss=Y_ss, C_ss=C_ss, I_ss=I_ss)


def solve_firm_path(N_path, Kap_path, Z_path, cal, country="D"):
    """Firm quantities along a transition path (purely contemporaneous).

    Returns dict: Y, w, mpk, mc — all shape (T,).
    """
    alpha = cal[f"alpha_{country}"]
    mc    = markup_ss(cal, country)

    Y   = Z_path * Kap_path ** alpha * N_path ** (1 - alpha)
    w   = mc * (1 - alpha) * Y / N_path
    mpk = mc * alpha * Y / Kap_path

    return dict(Y=Y, w=w, mpk=mpk, mc=mc)
