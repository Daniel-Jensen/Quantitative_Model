# SINGLE-POINT PERIOD MAP FOR THE RECURSIVE SOLUTION (THE IMAGE OF _inner_economy).
# Evaluates the full economy at ONE grid point (d, S) in QUANTITY FORM (realized
# payoffs = quantities x current prices), so [K_D,K_F,P_D,P_F,B_D,s] is exactly
# sufficient -- the same algebra that makes Bocola's P sufficient for N. There is
# NO separate household-wealth state: a bank's gross deposit obligation P_X IS the
# household's gross deposit claim, so household carried wealth is P_X/P_CES_X.
# The SEVEN pointwise unknowns are the SAME as transition.py's stacked system --
# [N_D, N_F, Kp_D, Kp_F, rdep_D, rdep_F, p] -- and the seven residuals are its
# per-period image: 2 bank capital-Euler (occasionally binding), 2 labour, the
# household deposit Euler (pins rdep_D), deposit-UIP (pins rdep_F), goods-D.
# Deposits clear BY QUANTITY (A = dep/P_CES); the rate is what makes households
# willing to hold exactly what banks fund. The banker valuations alpha, Q_b and
# the multiplier mu are READ OFF the recursions given a FROZEN continuation (the
# previous time-iteration iterate).
#
# OCCASIONALLY-BINDING IC -- Bocola Prop. 1 CLOSED FORM: mu = max{1 - E[Om]*R*n /
# (lambda*divertable assets), 0}, explicit and bounded in [0,1); the capital Euler
# E[Om(R_K-R)] = lambda_K*mu is the residual pinning K'. Same KKT as transition.py's
# Fischer-Burmeister, in the recursive-stable form (user-approved 2026-07-30;
# transition.py / bank.py FB left untouched).
#
# Two-branch SDF convention (kept exactly as bank.py/risk_branch.py): the d'=0
# continuation discounts with beta_inter (bank.py's proxy); the d'=1 branch with
# the income-SDF Lam_d = beta_inter*(Y_d/Y_nd)^-sigma. Only the EXPECTATION is
# genuine now. pi=0 (no_default) nests the deterministic period map exactly.
#
# Tier-3 cuts (documented, reversible): cross-border bond splits at SS shares;
# B_F at SS; TPI off; the wc-wedge rate component uses current rdep.
import numpy as np

from firms import solve_firm_path, markup_ss
from capital import solve_capital_path
from trade import ces_price, import_demand, trade_balance
from state_grid import default_prob
from expectations import gh_nodes

IK_D, IK_F, IP_D, IP_F, IB_D, IS = range(6)

# the SEVEN pointwise Newton unknowns per regime (transition.py's stacked order)
SOLVE7 = ("N_D", "N_F", "Kp_D", "Kp_F", "rdep_D", "rdep_F", "p")
# valuations + household aggregates read off the recursions/closure (stored rules)
DERIVED4 = ("alpha_D", "alpha_F", "Q_bD", "Q_bF")


def _firm_capital(N, K, Kp, Z, cal, c):
    # ELEMENTWISE FIRM + CAPITAL BLOCK AT ONE POINT (untouched blocks).
    f = solve_firm_path(np.array([N]), np.array([K]), np.array([Z]), cal, c)
    cp = solve_capital_path(np.array([Kp]), K, 1.0, f["mpk"], cal, c,
                            Kap_lag_path=np.array([K]))
    return (f["Y"][0], f["w"][0], f["mpk"][0], cp["Q"][0], cp["rk"][0],
            cp["I"][0], cp["cap_profit"][0])


def _cont_capital(N_next, Kp, Kpp, Z, cal, c):
    # NEXT-PERIOD mpk' AND Q_K' FROM CONTINUATION RULES (for rk_next = t->t+1).
    # Clip the continuation next-capital to the Jermann-feasible band so a
    # transient-iterate Kp rule that dips too low does not raise inside the
    # untouched capital block (the guard is slack near the fixed point).
    delta, ksi = cal[f"delta_{c}"], cal[f"ksi_{c}"]
    floor_ratio = (1.0 - delta) - delta * ksi / (1.0 - ksi)
    Kpp = np.maximum(Kpp, 1.0001 * floor_ratio * Kp)
    f = solve_firm_path(N_next, np.full_like(N_next, Kp),
                        np.full_like(N_next, Z), cal, c)
    cp = solve_capital_path(Kpp, Kp, 1.0, f["mpk"], cal, c,
                            Kap_lag_path=np.full_like(Kpp, Kp))
    return f["mpk"], cp["Q"], f["Y"]


def point_residuals(S, d, x, cont, cal, ss, sproc, n_gh=7, no_default=False):
    # SEVEN UNIT-CONSISTENT RESIDUALS AT ONE POINT. x follows SOLVE7; cont is the
    # FROZEN continuation RuleSet (previous iterate). alpha/Q_b/mu AND the
    # household aggregates C/A are computed here and returned in `out` to be
    # stored as the next-iterate rules.
    N_D, N_F, Kp_D, Kp_F, rdep_D, rdep_F, p = x
    K_D, K_F, P_D, P_F, B_D, s = S
    bkD, bkF = ss["ss_bank_D"], ss["ss_bank_F"]

    # --- firms + capital (current) ----------------------------------------
    Y_D, w_D0, mpk_D, Q_D, _, I_D, capprof_D = _firm_capital(
        N_D, K_D, Kp_D, cal["Z_ss_D"], cal, "D")
    Y_F, w_F0, mpk_F, Q_F, _, I_F, capprof_F = _firm_capital(
        N_F, K_F, Kp_F, cal["Z_ss_F"], cal, "F")
    P_CES_D = ces_price(np.array([p]), cal, "D")[0]
    P_CES_F = ces_price(np.array([p]), cal, "F")[0]

    # --- government (per-period body, quantity form) ----------------------
    # Bocola fiscal rule tau = tau* + gamma_tau*B on the FIXED SS anchor: a default
    # lowers the (post-haircut) debt B*surv below b_gov_ss, so taxes FALL -- the
    # fiscal-relief leg of default. (No re-anchoring to the post-haircut stock: with
    # the correctly-sized debt the windfall is small and the bank-loss leg wins.)
    surv = 1.0 - float(d) * (1.0 - cal["recovery_rate_D"])
    anchor = ss["gs_D"]["b_gov_ss"]
    Tax_D = ss["gs_D"]["Tax_ss"] + cal["phi_lamb_D"] * (B_D * surv - anchor)
    Tax_F = ss["gs_F"]["Tax_ss"]
    coupon_D = cal["delta_b_D"] * B_D * surv

    # --- end-of-period stocks (payoffs at current prices; alpha/mu-free) ---
    shareF = cal["b_D_F_ss"] / cal["B_gov_D_ss"]
    b_D_D_lag = (1.0 - shareF) * B_D
    b_D_F_lag = shareF * B_D
    b_F_D_lag = cal["b_F_D_ss"]
    b_F_F_lag = cal["B_gov_F_ss"] - cal["b_F_D_ss"]
    db_D, db_F = cal["delta_b_D"], cal["delta_b_F"]

    eps, wq = gh_nodes(n_gh)
    pd = 0.0 if no_default else float(default_prob(s))
    s_next = ((1.0 - sproc["rho_s"]) * sproc["s_star"] + sproc["rho_s"] * s
              + sproc["sigma_s"] * eps)
    m = eps.size

    # current-period bank valuations = FROZEN previous-iterate rules at this state
    # (breaks the S'<->Q_b knot; coincide at the fixed point)
    Sm = np.atleast_2d(S)
    Q_bD_cur = float(cont.eval("Q_bD", d, Sm)[0])
    Q_bF_cur = float(cont.eval("Q_bF", d, Sm)[0])
    alpha_D_cur = float(cont.eval("alpha_D", d, Sm)[0])
    alpha_F_cur = float(cont.eval("alpha_F", d, Sm)[0])
    new_D = (cal["G_D"] + coupon_D - Tax_D) / Q_bD_cur
    Bp_D = (1.0 - cal["delta_b_D"]) * B_D * surv + new_D

    payD = surv * (db_D + (1.0 - db_D) * Q_bD_cur)
    payF = db_F + (1.0 - db_F) * Q_bF_cur
    X_D = ((mpk_D + (1.0 - cal["delta_D"]) * Q_D) * K_D
           + payD * b_D_D_lag + p * payF * b_F_D_lag)
    X_F = ((mpk_F + (1.0 - cal["delta_F"]) * Q_F) * K_F
           + payF * b_F_F_lag + payD * b_D_F_lag / p)
    ng_D, ng_F = X_D - P_D, X_F - P_F
    b_D_D_new = (1.0 - shareF) * Bp_D
    b_D_F_new = shareF * Bp_D
    # end-of-period portfolio valued at current PRICES (Q_b), not payoffs
    assets_D = Q_D * Kp_D + Q_bD_cur * b_D_D_new + p * Q_bF_cur * b_F_D_lag
    assets_F = Q_F * Kp_F + Q_bF_cur * b_F_F_lag + Q_bD_cur * b_D_F_new / p
    n_D = (1.0 - cal["f_D"]) * ng_D + cal["omega_ent_D"] * assets_D
    n_F = (1.0 - cal["f_F"]) * ng_F + cal["omega_ent_F"] * assets_F
    dep_D, dep_F = assets_D - n_D, assets_F - n_F
    Pp_D = (1.0 + rdep_D) * dep_D
    Pp_F = (1.0 + rdep_F) * dep_F

    # --- two-branch continuation objects (FROZEN cont) --------------------
    def _cont(d_next):
        # CONTINUATION RULE VALUES + NEXT-PERIOD CAPITAL RETURN, REGIME d_next.
        Sn = np.empty((m, 6))
        Sn[:, IK_D] = Kp_D; Sn[:, IK_F] = Kp_F
        Sn[:, IP_D] = Pp_D; Sn[:, IP_F] = Pp_F
        Sn[:, IB_D] = Bp_D
        Sn[:, IS] = s_next
        r = cont.eval_all(d_next, cont.grid.clip(Sn))
        # guard continuation outputs before they enter fractional powers / sqrt
        # (a deep default-regime iterate can push N, p, C negative -> NaN)
        for k in ("N_D", "N_F"):
            r[k] = np.maximum(r[k], 0.05)
        for k in ("C_D", "C_F"):
            r[k] = np.maximum(r[k], 1e-3)
        r["p"] = np.maximum(r["p"], 1e-2)
        mpkD_n, QKD_n, YD_n = _cont_capital(r["N_D"], Kp_D, r["Kp_D"],
                                            cal["Z_ss_D"], cal, "D")
        mpkF_n, QKF_n, YF_n = _cont_capital(r["N_F"], Kp_F, r["Kp_F"],
                                            cal["Z_ss_F"], cal, "F")
        rkD_n = (mpkD_n + (1.0 - cal["delta_D"]) * QKD_n) / Q_D - 1.0
        rkF_n = (mpkF_n + (1.0 - cal["delta_F"]) * QKF_n) / Q_F - 1.0
        return r, rkD_n, rkF_n, YD_n, YF_n

    r0, rkD_n0, rkF_n0, YD_n0, YF_n0 = _cont(0)
    if pd > 0.0:
        r1, rkD_n1, rkF_n1, YD_n1, YF_n1 = _cont(1)
    else:
        r1, rkD_n1, rkF_n1, YD_n1, YF_n1 = r0, rkD_n0, rkF_n0, YD_n0, YF_n0

    # branch SDFs (bank.py convention): d'=0 -> beta_inter; d'=1 -> income SDF
    biD, biF = cal["beta_inter_D"], cal["beta_inter_F"]
    Om0_D = biD * (cal["f_D"] + (1 - cal["f_D"]) * r0["alpha_D"])
    Om0_F = biF * (cal["f_F"] + (1 - cal["f_F"]) * r0["alpha_F"])
    Lam1_D = biD * (YD_n1 / np.maximum(YD_n0, 1e-9)) ** (-cal["sigma_D"])
    Lam1_F = biF * (YF_n1 / np.maximum(YF_n0, 1e-9)) ** (-cal["sigma_F"])
    Om1_D = Lam1_D * (cal["f_D"] + (1 - cal["f_D"]) * r1["alpha_D"])
    Om1_F = Lam1_F * (cal["f_F"] + (1 - cal["f_F"]) * r1["alpha_F"])

    w_nd, w_def = wq * (1.0 - pd), wq * pd

    def _E(v0, v1):
        # TWO-BRANCH EXPECTATION OVER eps_s NODES.
        return float(np.dot(w_nd, v0) + np.dot(w_def, v1))

    # --- banker FOC block (Bocola Prop. 1 CLOSED-FORM multiplier) -----------
    lKD, lKF = cal["lambda_K_D"], cal["lambda_K_F"]
    lbDD, lbFF = cal["lambda_bD_D"], cal["lambda_bF_F"]
    E_Om_D = _E(Om0_D, Om1_D)
    E_Om_F = _E(Om0_F, Om1_F)
    # divertable assets (same leverage term as the FB slack), frozen-price valued
    lev_D = max(lKD * Q_D * Kp_D + lbDD * Q_bD_cur * b_D_D_new
                + cal["lambda_bF_D"] * p * Q_bF_cur * b_F_D_lag, 1e-6)
    lev_F = max(lKF * Q_F * Kp_F + lbFF * Q_bF_cur * b_F_F_lag
                + cal["lambda_bD_F"] * Q_bD_cur * b_D_F_new / p, 1e-6)
    # closed-form mu in [0, mu_cap]: the lower max is Bocola's KKT switch; the
    # upper cap (guard mu<1) keeps alpha=E_Om R/(1-mu) finite when a default-regime
    # net worth n goes negative (deep post-haircut recession).
    _MU_CAP = 0.95
    mu_D = min(max(1.0 - E_Om_D * (1.0 + rdep_D) * n_D / lev_D, 0.0), _MU_CAP)
    mu_F = min(max(1.0 - E_Om_F * (1.0 + rdep_F) * n_F / lev_F, 0.0), _MU_CAP)
    # capital-Euler surplus E[Om(R_K - R)] - lambda_K*mu (the residual pinning K').
    # Normalise by the O(1) discount kernel E_Om so the residual is in return units
    # (dividing by lambda_K*mu_ss amplifies ~130x and wrecks the solver step).
    E_Om_rk_D = _E(Om0_D * (rkD_n0 - rdep_D), Om1_D * (rkD_n1 - rdep_D))
    E_Om_rk_F = _E(Om0_F * (rkF_n0 - rdep_F), Om1_F * (rkF_n1 - rdep_F))
    cap_eul_D = (E_Om_rk_D - lKD * mu_D) / E_Om_D
    cap_eul_F = (E_Om_rk_F - lKF * mu_F) / E_Om_F
    alpha_D_new = np.clip(E_Om_D * (1.0 + rdep_D) / (1.0 - mu_D), 0.3, 8.0)
    alpha_F_new = np.clip(E_Om_F * (1.0 + rdep_F) / (1.0 - mu_F), 0.3, 8.0)

    # bond Euler (D risky: haircut surv' on the d'=1 payoff; F safe), same mu
    payD_nd0 = db_D + (1 - db_D) * r0["Q_bD"]
    payD_nd1 = cal["recovery_rate_D"] * (db_D + (1 - db_D) * r1["Q_bD"])
    payF_nd0 = db_F + (1 - db_F) * r0["Q_bF"]
    payF_nd1 = db_F + (1 - db_F) * r1["Q_bF"]
    E_Om_payD = _E(Om0_D * payD_nd0, Om1_D * payD_nd1)
    E_Om_payF = _E(Om0_F * payF_nd0, Om1_F * payF_nd1)
    Q_bD_new = np.clip(E_Om_payD / (E_Om_D * (1.0 + rdep_D) + lbDD * mu_D),
                       0.2, 1.2)
    Q_bF_new = np.clip(E_Om_payF / (E_Om_F * (1.0 + rdep_F) + lbFF * mu_F),
                       0.2, 1.2)

    # --- working-capital wedge + dividends (untouched formulas) ------------
    r_wc_D = rdep_D + lKD * mu_D / E_Om_D
    r_wc_F = rdep_F + lKF * mu_F / E_Om_F
    zD, zF = cal["zeta_wc_D"], cal["zeta_wc_F"]
    w_D = w_D0 / (1.0 + zD * r_wc_D)
    w_F = w_F0 / (1.0 + zF * r_wc_F)
    wc_inc_D = zD * r_wc_D * w_D * N_D
    wc_inc_F = zF * r_wc_F * w_F * N_F
    mcD, mcF = markup_ss(cal, "D"), markup_ss(cal, "F")
    div_bank_D = cal["f_D"] * ng_D - cal["omega_ent_D"] * assets_D
    div_bank_F = cal["f_F"] * ng_F - cal["omega_ent_F"] * assets_F
    Div_D = (1 - mcD) * Y_D + capprof_D + div_bank_D + wc_inc_D
    Div_F = (1 - mcF) * Y_F + capprof_F + div_bank_F + wc_inc_F

    # --- household aggregate (BUDGET + DEPOSIT EULER, rep-agent GHH) ---------
    # Deposits clear BY QUANTITY: household composite deposits A = bank funding
    # dep / P_CES. Carried wealth = P/P_CES (the gross deposit claim = the bank's
    # gross obligation state). Consumption is the residual of the budget; the
    # deposit RATE is pinned by the household Euler (residual 5). beta_eff =
    # 1/(1+rdep_ss) makes the rep-agent Euler reproduce the HA aggregate at the SS;
    # hh_T anchors the SS budget so C = C_ss, A = A_ss exactly.
    frisch_D, frisch_F = cal["frisch_D"], cal["frisch_F"]
    sigD, sigF = cal["sigma_D"], cal["sigma_F"]
    A_D = dep_D / P_CES_D
    A_F = dep_F / P_CES_F
    inc_D = (w_D / P_CES_D) * N_D + (Div_D - Tax_D) / P_CES_D + ss.get("hh_T_D", 0.0)
    inc_F = (w_F / P_CES_F) * N_F + (Div_F - Tax_F) / P_CES_F + ss.get("hh_T_F", 0.0)
    C_D = np.clip(P_D / P_CES_D + inc_D - A_D, 0.3 * ss["C_D_ss"], 2.0 * ss["C_D_ss"])
    C_F = np.clip(P_F / P_CES_F + inc_F - A_F, 0.3 * ss["C_F_ss"], 2.0 * ss["C_F_ss"])

    # deposit Euler (GHH composite x = C - vN): (x)^-sigma = beta_eff*E[(1+r')x'^-sigma]
    vN_D = cal["chi_D"] * N_D ** (1 + 1 / frisch_D) / (1 + 1 / frisch_D)
    vN_F = cal["chi_F"] * N_F ** (1 + 1 / frisch_F) / (1 + 1 / frisch_F)
    vNp_D0 = cal["chi_D"] * r0["N_D"] ** (1 + 1 / frisch_D) / (1 + 1 / frisch_D)
    vNp_F0 = cal["chi_F"] * r0["N_F"] ** (1 + 1 / frisch_F) / (1 + 1 / frisch_F)
    vNp_D1 = cal["chi_D"] * r1["N_D"] ** (1 + 1 / frisch_D) / (1 + 1 / frisch_D)
    vNp_F1 = cal["chi_F"] * r1["N_F"] ** (1 + 1 / frisch_F) / (1 + 1 / frisch_F)
    xp_D0 = np.maximum(r0["C_D"] - vNp_D0, 1e-9)
    xp_F0 = np.maximum(r0["C_F"] - vNp_F0, 1e-9)
    xp_D1 = np.maximum(r1["C_D"] - vNp_D1, 1e-9)
    xp_F1 = np.maximum(r1["C_F"] - vNp_F1, 1e-9)
    rp_D0 = (1.0 + rdep_D) * P_CES_D / ces_price(r0["p"], cal, "D") - 1.0
    rp_F0 = (1.0 + rdep_F) * P_CES_F / ces_price(r0["p"], cal, "F") - 1.0
    rp_D1 = (1.0 + rdep_D) * P_CES_D / ces_price(r1["p"], cal, "D") - 1.0
    rp_F1 = (1.0 + rdep_F) * P_CES_F / ces_price(r1["p"], cal, "F") - 1.0
    beff_D = 1.0 / (1.0 + cal["r_dep_D_target"])
    beff_F = 1.0 / (1.0 + cal["r_dep_F_target"])
    E_mu_D = _E((1.0 + rp_D0) * xp_D0 ** (-sigD), (1.0 + rp_D1) * xp_D1 ** (-sigD))
    E_mu_F = _E((1.0 + rp_F0) * xp_F0 ** (-sigF), (1.0 + rp_F1) * xp_F1 ** (-sigF))
    xC_D = max(C_D - vN_D, 1e-9)
    xC_F = max(C_F - vN_F, 1e-9)
    euler_D = xC_D ** (-sigD) / (beff_D * E_mu_D) - 1.0
    euler_F = xC_F ** (-sigF) / (beff_F * E_mu_F) - 1.0

    IM_D = import_demand(np.array([p]), np.array([C_D]), np.array([P_CES_D]), cal, "D")[0]
    IM_F = import_demand(np.array([p]), np.array([C_F]), np.array([P_CES_F]), cal, "F")[0]
    NX_D, NX_F = (lambda a, b: (a[0], b[0]))(
        *trade_balance(np.array([p]), np.array([IM_D]), np.array([IM_F])))

    Ep_next = float(np.dot(wq, r0["p"]))

    # --- the SEVEN residuals (transition.py's per-period image) -------------
    slack_D = alpha_D_cur * n_D - lev_D
    slack_F = alpha_F_cur * n_F - lev_F
    res = np.array([
        cap_eul_D,                                               # 1 cap Euler D -> Kp_D
        cap_eul_F,                                               # 2 cap Euler F -> Kp_F
        (cal["chi_D"] * N_D ** (1 / frisch_D) - w_D / P_CES_D)    # 3 lab_D -> N_D
        / (w_D / P_CES_D),
        (cal["chi_F"] * N_F ** (1 / frisch_F) - w_F / P_CES_F)    # 4 lab_F -> N_F
        / (w_F / P_CES_F),
        euler_D,                                                 # 5 deposit Euler -> rdep_D
        (1.0 + rdep_D) - (1.0 + rdep_F) * Ep_next / p,           # 6 UIP -> rdep_F
        (Y_D - P_CES_D * C_D - I_D - NX_D - cal["G_D"])          # 7 goods_D -> p
        / ss["ss_firm_D"]["Y_ss"],
    ])
    out = dict(mu_D=mu_D, mu_F=mu_F, n_D=n_D, n_F=n_F, Y_D=Y_D, I_D=I_D,
               alpha_D=alpha_D_new, alpha_F=alpha_F_new,
               Q_bD=Q_bD_new, Q_bF=Q_bF_new, C_D=C_D, C_F=C_F, A_D=A_D, A_F=A_F,
               inc_D=inc_D, inc_F=inc_F, w_D=w_D, dep_D=dep_D, dep_F=dep_F,
               Pp_D=Pp_D, Pp_F=Pp_F, Bp_D=Bp_D, slack_D=slack_D, slack_F=slack_F)
    return res, out
