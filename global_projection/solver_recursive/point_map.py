# THE PERIOD MAP: THE WHOLE ECONOMY AT ONE GRID POINT, GIVEN THE CONTINUATION RULES.
from types import SimpleNamespace

import numpy as np

from global_projection.blocks.firms import solve_firm_path, markup_ss
from global_projection.blocks.capital import solve_capital_path
from global_projection.blocks.trade import ces_price, import_demand, trade_balance, size_ratio
from global_projection.solver_recursive.state_grid import default_prob

from global_projection.solver_recursive.state_grid import (IK_D, IK_F, IP_D, IP_F, IBDD, IBDF,
                                         IBFD, IV, IS, IZ, IM, IO, NSTATE)



def gh_nodes(n=7):
    # GAUSS-HERMITE NODES AND WEIGHTS FOR ONE STANDARD-NORMAL SHOCK (PROBABILISTS' RULE).
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return x, w / w.sum()


# smoothing scale for the guards below
_GUARD_EPS = 1e-5
# cap on mu, so alpha stays finite in deep default states
_MU_CAP = 0.95


def _smax(x, floor, eps=1e-3):
    # SMOOTH MAX(x, floor), SO THE JACOBIAN STAYS VALID.
    return floor + 0.5 * ((x - floor) + np.sqrt((x - floor) ** 2 + eps ** 2))


def _smin(x, cap, eps=1e-3):
    # SMOOTH MIN(x, cap).
    return cap - 0.5 * ((cap - x) + np.sqrt((cap - x) ** 2 + eps ** 2))


def _sclip(x, lo, hi, eps=None):
    # SMOOTH CLIP BETWEEN lo AND hi.
    return _smin(_smax(x, lo, eps or _GUARD_EPS), hi, eps or _GUARD_EPS)


def _fb(a, b, eps):
    # SMOOTHED FISCHER-BURMEISTER: ZERO WHEN a, b >= 0 AND a*b = eps^2/2.
    return a + b - np.sqrt(a * a + b * b + eps * eps)


def _regime_weights(wq, pd, reg):
    # PROBABILITY OF EACH (s' NODE, REGIME) CELL.
    return [wq * (pd if d_n else (1.0 - pd)) for d_n in reg]


def _expect(wgt, vals):
    # EXPECTATION OVER THE (s' NODE, REGIME) CELLS.
    return float(sum(np.dot(wgt[j], vals[j]) for j in range(len(wgt))))


def _ghh(C, N, chi, fr, floor, eps):
    # GHH CONSUMPTION COMPOSITE x = C - v(N), SMOOTHLY FLOORED.
    return _smax(C - chi * N ** (1 + 1 / fr) / (1 + 1 / fr), floor, eps)


def _jermann_floor(K, cal, c):
    # LOWEST FEASIBLE NEXT CAPITAL UNDER JERMANN ADJUSTMENT COSTS.
    delta, ksi = cal[f"delta_{c}"], cal[f"ksi_{c}"]
    floor_ratio = (1.0 - delta) - delta * ksi / (1.0 - ksi)
    return 1.0001 * floor_ratio * K


def _firm_capital(N, K, Kp, Z, cal, c):
    # FIRMS AND CAPITAL AT ONE POINT; Kp IS FLOORED SO THE MAP CAN ALWAYS BE EVALUATED.
    Kp = max(Kp, _jermann_floor(K, cal, c))
    f = solve_firm_path(np.array([N]), np.array([K]), np.array([Z]), cal, c)
    cp = solve_capital_path(np.array([Kp]), K, 1.0, f["mpk"], cal, c,
                            Kap_lag_path=np.array([K]))
    return (f["Y"][0], f["w"][0], f["mpk"][0], cp["Q"][0], cp["rk"][0],
            cp["I"][0], cp["cap_profit"][0])


def _cont_capital(N_next, Kp, Kpp, Z, cal, c):
    # NEXT PERIOD'S mpk AND Q_K FROM THE CONTINUATION RULES.
    Kpp = np.maximum(Kpp, _jermann_floor(Kp, cal, c))
    f = solve_firm_path(N_next, np.full_like(N_next, Kp),
                        np.full_like(N_next, Z), cal, c)
    cp = solve_capital_path(Kpp, Kp, 1.0, f["mpk"], cal, c,
                            Kap_lag_path=np.full_like(Kpp, Kp))
    return f["mpk"], cp["Q"], f["Y"]


def _read_state_and_regime(S, d, x, cont, cal, ss):
    # UNPACK THE STATE AND THE UNKNOWNS, AND SET THE REGIME.
    v = SimpleNamespace()
    (v.N_D, v.N_F, v.Kp_D, v.Kp_F, v.rdep_D, v.rdep_F, v.p,
     v.Q_bD, v.b_DF_new, v.Q_bF, v.b_FD_new, v.A_D, v.A_F, v.x_cb) = x
    (v.K_D, v.K_F, v.P_D, v.P_F, v.b_D_D_lag, v.b_D_F_lag, v.b_F_D_lag,
     v.V_dep, v.s, v.Z_D, v.M_lag, v.O_lag) = S
    v.S = S
    # gross debt: the banks' holdings plus the Eurosystem's
    v.B_D = v.b_D_D_lag + v.b_D_F_lag + v.M_lag
    v.reg, v.nreg = cont.reg, cont.n_regimes
    v.d = d
    v.d_reg = v.reg[d]
    # sz = size_F/size_D, used wherever the two countries are added up
    v.sz = size_ratio(cal)
    v.surv = 1.0 - float(v.d_reg) * (1.0 - cal["recovery_rate_D"])
    # the household's claim: V = W_D - P_D
    v.W_D = v.P_D + v.V_dep
    v.bkD, v.bkF = ss["ss_bank_D"], ss["ss_bank_F"]
    # floors for the GHH composite
    v.X_FLOOR_D, v.X_FLOOR_F = 0.05 * ss["C_D_ss"], 0.05 * ss["C_F_ss"]
    v.X_EPS = 1e-3 * ss["C_D_ss"]
    v.Sm = np.atleast_2d(S)
    _tpi_purchase_and_price(v, cal, ss)
    return v


def _tpi_floor(Q_bF, cal):
    # THE TPI FLOOR: THE D BOND PRICE AT WHICH y_D - y_F EQUALS THE CAP.
    db_D, db_F = cal["delta_b_D"], cal["delta_b_F"]
    return db_D / (db_D + db_F * (1.0 - Q_bF) / Q_bF + cal["tpi_cap_bp"] / 4e4)


def _tpi_purchase_and_price(v, cal, ss):
    # PURCHASE AND PRICE FROM ONE SMOOTH VARIABLE x (GARCIA-ZANGWILL).
    v.tpi = bool(cal["tpi_on"])
    v.Q_floor = _tpi_floor(v.Q_bF, cal)
    v.Q_bD_slot = v.Q_bD
    B = cal["B_gov_D_ss"]
    if v.tpi and not v.d_reg and cal["tpi_gz"]:
        # x > 0 buys at the floor; x < 0 is the price's gap above it
        eps = cal["tpi_eps"]
        v.m_cb = B * _smax(v.x_cb, 0.0, eps)
        v.Q_bD = v.Q_floor + ss["Q_bD_ss"] * _smax(-v.x_cb, 0.0, eps)
    else:
        # otherwise x is the purchase itself (TPI off, in default, or the solver's FB form)
        v.m_cb = B * v.x_cb


def _production(v, cal):
    # FIRMS AND CAPITAL AT THE CURRENT STATE.
    (v.Y_D, v.w_D0, v.mpk_D, v.Q_D, _, v.I_D, v.capprof_D) = _firm_capital(
        v.N_D, v.K_D, v.Kp_D, v.Z_D, cal, "D")
    (v.Y_F, v.w_F0, v.mpk_F, v.Q_F, _, v.I_F, v.capprof_F) = _firm_capital(
        v.N_F, v.K_F, v.Kp_F, cal["Z_ss_F"], cal, "F")
    v.P_CES_D = ces_price(np.array([v.p]), cal, "D")[0]
    v.P_CES_F = ces_price(np.array([v.p]), cal, "F")[0]


def _government(v, cal, ss):
    # THE BOHN TAX ON THE SURVIVING DEBT; THE DEFAULT REGIME RE-ANCHORS TO THE POST-HAIRCUT STOCK.
    if v.d_reg:
        anchor = cal["recovery_rate_D"] * ss["gs_D"]["b_gov_ss"]
        Tax_base = cal["G_D"] + cal["delta_b_D"] * anchor * (1.0 - ss["gs_D"]["Q_B_ss"])
    else:
        anchor, Tax_base = ss["gs_D"]["b_gov_ss"], ss["gs_D"]["Tax_ss"]
    v.Tax_D = Tax_base + ss["gs_D"]["gamma_tau"] * (v.B_D * v.surv - anchor)
    # Eurosystem P&L, remitted by capital key: D's share cuts issuance, F's cuts F taxes
    db_D = cal["delta_b_D"]
    v.Xi_D = v.surv * (db_D + (1.0 - db_D) * v.Q_bD)
    v.Pi_cb = v.Xi_D * v.M_lag - v.O_lag
    kD = cal["tpi_key_D"]
    v.Tax_F = ss["gs_F"]["Tax_ss"] - (1.0 - kD) * v.Pi_cb / (v.sz * v.p)
    v.coupon_D = cal["delta_b_D"] * v.B_D * v.surv
    v.new_D = (cal["G_D"] + v.coupon_D - v.Tax_D - kD * v.Pi_cb) / v.Q_bD
    v.Bp_D = (1.0 - cal["delta_b_D"]) * v.B_D * v.surv + v.new_D


def _balance_sheets(v, cont, cal):
    # BANK PAYOFFS, NET WORTH, DEPOSITS AND HOUSEHOLD CLAIMS.
    db_D, db_F = cal["delta_b_D"], cal["delta_b_F"]
    # the F bank holds the rest of the fixed F stock
    v.b_F_F_lag = cal["B_gov_F_ss"] - v.b_F_D_lag / v.sz
    # today's alpha and r_wc come from the frozen rules (exact at the nodes)
    v.alpha_D_cur = float(cont.eval("alpha_D", v.d, v.Sm)[0])
    v.alpha_F_cur = float(cont.eval("alpha_F", v.d, v.Sm)[0])
    v.r_wc_D_cur = float(cont.eval("r_wc_D", v.d, v.Sm)[0])
    v.r_wc_F_cur = float(cont.eval("r_wc_F", v.d, v.Sm)[0])
    # working-capital loans, a bank asset
    zD_, zF_ = cal["zeta_wc_D"], cal["zeta_wc_F"]
    v.L_wc_D = zD_ * (v.w_D0 / (1.0 + zD_ * v.r_wc_D_cur)) * v.N_D
    v.L_wc_F = zF_ * (v.w_F0 / (1.0 + zF_ * v.r_wc_F_cur)) * v.N_F

    payD_now = v.Xi_D
    payF_now = db_F + (1.0 - db_F) * v.Q_bF
    # the D bank is also repaid last period's TPI claim
    X_D = ((v.mpk_D + (1.0 - cal["delta_D"]) * v.Q_D) * v.K_D
           + payD_now * v.b_D_D_lag + v.p * payF_now * v.b_F_D_lag + v.O_lag)
    X_F = ((v.mpk_F + (1.0 - cal["delta_F"]) * v.Q_F) * v.K_F
           + payF_now * v.b_F_F_lag + payD_now * (v.b_D_F_lag / v.sz) / v.p)
    v.ng_D, v.ng_F = X_D - v.P_D, X_F - v.P_F
    # the floors stop a transient iterate giving a bank a negative book
    v.b_D_F_new = _smax(v.b_DF_new, 1e-4, 1e-5)
    # the TPI book, held to maturity, and the safe claim that paid for it
    v.M_cb_new = (1.0 - db_D) * v.surv * v.M_lag + v.m_cb
    v.Z_cb = v.Q_bD * v.M_cb_new
    b_DD = v.Bp_D - v.b_D_F_new - v.M_cb_new
    # with the TPI on, clearing is exact: b_DD = 0 is a genuine corner
    v.b_D_D_new = b_DD if v.tpi else _smax(b_DD, 1e-4, 1e-5)
    v.b_F_D_new = _smax(v.b_FD_new, 1e-4, 1e-5)
    v.b_F_F_new = _smax(cal["B_gov_F_ss"] - v.b_F_D_new / v.sz, 1e-4, 1e-5)
    # end-of-period assets at current prices
    v.assets_D = (v.Q_D * v.Kp_D + v.Q_bD * v.b_D_D_new + v.p * v.Q_bF * v.b_F_D_new
                  + v.L_wc_D + v.Z_cb)
    v.assets_F = (v.Q_F * v.Kp_F + v.Q_bF * v.b_F_F_new
                  + v.Q_bD * (v.b_D_F_new / v.sz) / v.p + v.L_wc_F)
    v.n_D = (1.0 - cal["f_D"]) * v.ng_D + cal["omega_ent_D"] * v.assets_D
    v.n_F = (1.0 - cal["f_F"]) * v.ng_F + cal["omega_ent_F"] * v.assets_F
    # Bocola's net-worth floor, for the deep default corners only
    nwf = cal.get("nw_floor_frac", 0.0)
    if nwf > 0.0:
        v.n_D = _smax(v.n_D, nwf * v.bkD["n_ss"])
        v.n_F = _smax(v.n_F, nwf * v.bkF["n_ss"])
    v.dep_D, v.dep_F = v.assets_D - v.n_D, v.assets_F - v.n_F
    # next period's deposit bill, net of the loan repaid at the lending rate (Bocola)
    v.Pp_D = (1.0 + v.rdep_D) * v.dep_D - (1.0 + v.r_wc_D_cur) * v.L_wc_D
    v.Pp_F = (1.0 + v.rdep_F) * v.dep_F - (1.0 + v.r_wc_F_cur) * v.L_wc_F
    # one union-wide deposit market; the gap is the cross-border position
    v.dep_union = v.dep_D + v.sz * v.p * v.dep_F
    v.save_union = v.A_D * v.P_CES_D + v.sz * v.p * v.A_F * v.P_CES_F
    v.nfa_dep_D = v.A_D * v.P_CES_D - v.dep_D
    # the household's claim mirrors the bank's obligation
    v.Wp_D = (1.0 + v.rdep_D) * v.A_D * v.P_CES_D - (1.0 + v.r_wc_D_cur) * v.L_wc_D
    # the cross-border position, carried at the D deposit rate
    v.Vp_dep = (1.0 + v.rdep_D) * v.nfa_dep_D
    # what the Eurosystem owes the D banks next period
    v.Op_cb = (1.0 + v.rdep_D) * v.Z_cb


def _next_states(v, sproc, eps):
    # THE NEXT-PERIOD STATE AT EACH s' QUADRATURE NODE.
    Z_next = (1.0 - sproc["rho_z"]) * sproc["z_star"] + sproc["rho_z"] * v.Z_D
    s_next = ((1.0 - sproc["rho_s"]) * sproc["s_star"] + sproc["rho_s"] * v.s
              + sproc["sigma_s"] * eps)
    Sn = np.empty((eps.size, NSTATE))
    Sn[:, IK_D] = v.Kp_D; Sn[:, IK_F] = v.Kp_F
    Sn[:, IP_D] = v.Pp_D; Sn[:, IP_F] = v.Pp_F
    Sn[:, IBDD] = v.b_D_D_new; Sn[:, IBDF] = v.b_D_F_new
    Sn[:, IBFD] = v.b_F_D_new
    Sn[:, IV] = v.Vp_dep
    Sn[:, IS] = s_next
    Sn[:, IZ] = Z_next
    Sn[:, IM] = v.M_cb_new; Sn[:, IO] = v.Op_cb
    return Sn, Z_next


def _regime_continuation(v, j, Sn, Z_next, cont, cal):
    # CONTINUATION RULES AND NEXT-PERIOD CAPITAL RETURNS IN REGIME j.
    r = cont.eval_all(j, cont.grid.clip(Sn))
    # keep continuation values positive before powers and roots
    for k in ("N_D", "N_F"):
        r[k] = np.maximum(r[k], 0.05)
    for k in ("C_D", "C_F"):
        r[k] = np.maximum(r[k], 1e-3)
    r["p"] = np.maximum(r["p"], 1e-2)
    mpkD_n, QKD_n, _ = _cont_capital(r["N_D"], v.Kp_D, r["Kp_D"], Z_next, cal, "D")
    mpkF_n, QKF_n, _ = _cont_capital(r["N_F"], v.Kp_F, r["Kp_F"], cal["Z_ss_F"], cal, "F")
    rkD_n = (mpkD_n + (1.0 - cal["delta_D"]) * QKD_n) / v.Q_D - 1.0
    rkF_n = (mpkF_n + (1.0 - cal["delta_F"]) * QKF_n) / v.Q_F - 1.0
    return r, rkD_n, rkF_n


def _continuation(v, cont, cal, sproc, n_gh, no_default):
    # QUADRATURE WEIGHTS AND THE CONTINUATION IN EACH REGIME.
    eps, wq = gh_nodes(n_gh)
    pd = 0.0 if no_default else float(default_prob(v.s))
    Sn, Z_next = _next_states(v, sproc, eps)
    v.wgt = _regime_weights(wq, pd, v.reg)
    v.R, v.RKD, v.RKF = [], [], []
    for j in range(v.nreg):
        # a regime with zero weight reuses regime 0, so pi = 0 nests the risk-free model exactly
        if j > 0 and not np.any(v.wgt[j] > 0.0):
            v.R.append(v.R[0]); v.RKD.append(v.RKD[0]); v.RKF.append(v.RKF[0])
            continue
        rj, rkDj, rkFj = _regime_continuation(v, j, Sn, Z_next, cont, cal)
        v.R.append(rj); v.RKD.append(rkDj); v.RKF.append(rkFj)


def _discount_kernels(v, cont, cal):
    # THE HOUSEHOLD SDF AND THE BANKER'S KERNEL OMEGA IN EACH REGIME.
    nreg = v.nreg
    frD, frF = cal["frisch_D"], cal["frisch_F"]
    sgD, sgF = cal["sigma_D"], cal["sigma_F"]
    x_cur_D = _ghh(float(cont.eval("C_D", v.d, v.Sm)[0]), v.N_D, cal["chi_D"], frD,
                   v.X_FLOOR_D, v.X_EPS)
    x_cur_F = _ghh(float(cont.eval("C_F", v.d, v.Sm)[0]), v.N_F, cal["chi_F"], frF,
                   v.X_FLOOR_D, v.X_EPS)
    v.XN_D = [_ghh(v.R[j]["C_D"], v.R[j]["N_D"], cal["chi_D"], frD, v.X_FLOOR_D, v.X_EPS)
              for j in range(nreg)]
    XN_F = [_ghh(v.R[j]["C_F"], v.R[j]["N_F"], cal["chi_F"], frF, v.X_FLOOR_D, v.X_EPS)
            for j in range(nreg)]
    Lam_D = [cal["beta_inter_D"] * (x_cur_D / v.XN_D[j]) ** sgD for j in range(nreg)]
    Lam_F = [cal["beta_inter_F"] * (x_cur_F / XN_F[j]) ** sgF for j in range(nreg)]
    v.Om_D = [Lam_D[j] * (cal["f_D"] + (1 - cal["f_D"]) * v.R[j]["alpha_D"])
              for j in range(nreg)]
    v.Om_F = [Lam_F[j] * (cal["f_F"] + (1 - cal["f_F"]) * v.R[j]["alpha_F"])
              for j in range(nreg)]
    v.E_Om_D = _expect(v.wgt, v.Om_D)
    v.E_Om_F = _expect(v.wgt, v.Om_F)


def _incentive_constraint(v, cal):
    # BOCOLA'S CLOSED-FORM MULTIPLIER, THE CAPITAL EULER AND THE FRANCHISE VALUE.
    lKD, lKF = cal["lambda_K_D"], cal["lambda_K_F"]
    lbDD, lbFF = cal["lambda_bD_D"], cal["lambda_bF_F"]
    # divertable assets; the TPI claim counts at the bond's lambda, so a swap leaves them unchanged
    lev_D = max(lKD * v.Q_D * v.Kp_D + lbDD * v.Q_bD * v.b_D_D_new
                + cal["lambda_bF_D"] * v.p * v.Q_bF * v.b_F_D_new
                + lKD * v.L_wc_D + lbDD * v.Z_cb, 1e-6)
    lev_F = max(lKF * v.Q_F * v.Kp_F + lbFF * v.Q_bF * v.b_F_F_new
                + cal["lambda_bD_F"] * v.Q_bD * (v.b_D_F_new / v.sz) / v.p
                + lKF * v.L_wc_F, 1e-6)
    v.lev_D, v.lev_F = lev_D, lev_F
    # mu = max(., 0) is the KKT switch; the cap is only a guard
    v.mu_D = float(_smin(max(1.0 - v.E_Om_D * (1.0 + v.rdep_D) * v.n_D / v.lev_D, 0.0),
                         _MU_CAP, _GUARD_EPS))
    v.mu_F = float(_smin(max(1.0 - v.E_Om_F * (1.0 + v.rdep_F) * v.n_F / v.lev_F, 0.0),
                         _MU_CAP, _GUARD_EPS))
    # the capital Euler, in return units
    E_Om_rk_D = _expect(v.wgt, [v.Om_D[j] * (v.RKD[j] - v.rdep_D) for j in range(v.nreg)])
    E_Om_rk_F = _expect(v.wgt, [v.Om_F[j] * (v.RKF[j] - v.rdep_F) for j in range(v.nreg)])
    v.cap_eul_D = (E_Om_rk_D - lKD * v.mu_D) / v.E_Om_D
    v.cap_eul_F = (E_Om_rk_F - lKF * v.mu_F) / v.E_Om_F
    # alpha = E[Om]R/(1-mu), capped so a divergent iterate stays finite
    alpha_cap = float(cal.get("alpha_cap", 40.0))
    v.alpha_D_new = float(_sclip(v.E_Om_D * (1.0 + v.rdep_D) / (1.0 - v.mu_D), 0.05, alpha_cap))
    v.alpha_F_new = float(_sclip(v.E_Om_F * (1.0 + v.rdep_F) / (1.0 - v.mu_F), 0.05, alpha_cap))
    # slack of the same constraint, so mu*slack = 0
    v.slack_D = v.alpha_D_cur * v.n_D - v.lev_D
    v.slack_F = v.alpha_F_cur * v.n_F - v.lev_F


def _bond_demands(v, cal):
    # BOTH BANKS' FOCS FOR BOTH SOVEREIGNS.
    nreg, R, wgt = v.nreg, v.R, v.wgt
    db_D, db_F = cal["delta_b_D"], cal["delta_b_F"]
    hc = [cal["recovery_rate_D"] if d_n else 1.0 for d_n in v.reg]
    payD_gross = [db_D + (1.0 - db_D) * R[j]["Q_bD"] for j in range(nreg)]
    payD = [hc[j] * payD_gross[j] for j in range(nreg)]
    payF = [db_F + (1.0 - db_F) * R[j]["Q_bF"] for j in range(nreg)]
    v.E_Om_payD = _expect(wgt, [v.Om_D[j] * payD[j] for j in range(nreg)])
    # the bond-price decomposition legs (diagnostics only)
    v.E_payD = _expect(wgt, payD)
    v.E_payD_nodef = _expect(wgt, [payD_gross[j] if not v.reg[j] else payD_gross[0]
                                   for j in range(nreg)])
    v.E_payF = _expect(wgt, payF)
    v.E_Om_payF = _expect(wgt, [v.Om_F[j] * payF[j] for j in range(nreg)])
    # the D bank's F-bond FOC, with its cross-border adjustment cost
    v.E_Om_payF_D = _expect(wgt, [v.Om_D[j] * payF[j] for j in range(nreg)])
    v.dmd_F_home = v.E_Om_F * (1.0 + v.rdep_F) + cal["lambda_bF_F"] * v.mu_F
    v.dmd_F_for = v.E_Om_D * (1.0 + v.rdep_D) + cal["lambda_bF_D"] * v.mu_D
    v.adj_D = 1.0 + cal["psi_bF_D"] * (v.b_F_D_new - cal["b_F_D_ss"]) / cal["B_gov_F_ss"]
    # the D bank's D-bond demand price
    v.dmd_D = v.E_Om_D * (1.0 + v.rdep_D) + cal["lambda_bD_D"] * v.mu_D
    # the F bank's D-bond FOC, with its cross-border adjustment cost
    v.E_Om_payD_F = _expect(wgt, [v.Om_F[j] * payD[j] for j in range(nreg)])
    v.dmd_F = v.E_Om_F * (1.0 + v.rdep_F) + cal["lambda_bD_F"] * v.mu_F
    v.adj_F = 1.0 + cal["psi_bD_F"] * (v.b_D_F_new - cal["b_D_F_ss"]) / cal["B_gov_D_ss"]


def _tpi_rule(v, cal, ss):
    # THE D BANK'S D-BOND FOC (A COMPLEMENTARITY PAIR WITH THE TPI ON) AND THE TPI RESIDUAL.
    v.foc_D = (v.E_Om_payD - v.dmd_D * v.Q_bD) / v.E_Om_D
    if not v.tpi:
        v.res_bondD, v.res_tpi = v.foc_D, v.x_cb
        return
    eps = cal["tpi_eps"]
    # with the TPI on, the bank may sell everything: b_DD >= 0 against foc_D <= 0
    v.res_bondD = _fb(-v.foc_D, v.b_D_D_new / cal["B_gov_D_ss"], eps)
    if v.d_reg:
        v.res_tpi = v.x_cb
    elif cal["tpi_gz"]:
        v.res_tpi = (v.Q_bD_slot - v.Q_bD) / ss["Q_bD_ss"]
    else:
        v.res_tpi = _fb(v.x_cb, (v.Q_bD - v.Q_floor) / ss["Q_bD_ss"], eps)


def _wages_and_dividends(v, cal):
    # THE WORKING-CAPITAL RATE r_wc = rdep + lambda_K*mu/E[Om], THE NET WAGE AND DIVIDENDS.
    lKD, lKF = cal["lambda_K_D"], cal["lambda_K_F"]
    v.r_wc_D = v.rdep_D + lKD * v.mu_D / v.E_Om_D
    v.r_wc_F = v.rdep_F + lKF * v.mu_F / v.E_Om_F
    zD, zF = cal["zeta_wc_D"], cal["zeta_wc_F"]
    v.w_D = v.w_D0 / (1.0 + zD * v.r_wc_D)
    v.w_F = v.w_F0 / (1.0 + zF * v.r_wc_F)
    wc_inc_D = zD * v.r_wc_D * v.w_D * v.N_D
    wc_inc_F = zF * v.r_wc_F * v.w_F * v.N_F
    mcD, mcF = markup_ss(cal, "D"), markup_ss(cal, "F")
    div_bank_D = cal["f_D"] * v.ng_D - cal["omega_ent_D"] * v.assets_D
    div_bank_F = cal["f_F"] * v.ng_F - cal["omega_ent_F"] * v.assets_F
    # wc_rebate = 1 would hand the financing income to households; the default 0 keeps it in the bank
    reb = float(cal.get("wc_rebate", 0.0))
    v.Div_D = (1 - mcD) * v.Y_D + v.capprof_D + div_bank_D + reb * wc_inc_D
    v.Div_F = (1 - mcF) * v.Y_F + v.capprof_F + div_bank_F + reb * wc_inc_F


def _households(v, cal, ss):
    # REPRESENTATIVE GHH HOUSEHOLDS: BUDGET, CONSUMPTION, DEPOSIT EULERS, TRADE.
    nreg, R = v.nreg, v.R
    frisch_D, frisch_F = cal["frisch_D"], cal["frisch_F"]
    sigD, sigF = cal["sigma_D"], cal["sigma_F"]
    # the firm's working-capital receipts, zeta*r_wc*w*N + L_wc, go to its owners each period; the
    # repayment (1+r_wc)*L_wc is charged to the household's claim next period (Wp). A constant
    # stood in for this flow before, which leaked its movement into the goods market off the SS.
    # wc_flow_D/F (default 1) weight the flow only for the solver's homotopy (_stage_wc).
    wcf_D = cal.get("wc_flow_D", 1.0) * (cal["zeta_wc_D"] * v.r_wc_D * v.w_D * v.N_D + v.L_wc_D)
    wcf_F = cal.get("wc_flow_F", 1.0) * (cal["zeta_wc_F"] * v.r_wc_F * v.w_F * v.N_F + v.L_wc_F)
    v.inc_D = ((v.w_D / v.P_CES_D) * v.N_D + (v.Div_D - v.Tax_D + wcf_D) / v.P_CES_D
               + ss.get("hh_T_D", 0.0))
    v.inc_F = ((v.w_F / v.P_CES_F) * v.N_F + (v.Div_F - v.Tax_F + wcf_F) / v.P_CES_F
               + ss.get("hh_T_F", 0.0))
    # F's claim: its bank's obligation less the cross-border position
    v.W_F = v.P_F - v.V_dep / (v.sz * v.p)
    # consumption, smoothly bounded
    v.C_D = float(_sclip(v.W_D / v.P_CES_D + v.inc_D - v.A_D,
                         0.15 * ss["C_D_ss"], 3.0 * ss["C_D_ss"], v.X_EPS))
    v.C_F = float(_sclip(v.W_F / v.P_CES_F + v.inc_F - v.A_F,
                         0.15 * ss["C_F_ss"], 3.0 * ss["C_F_ss"], v.X_EPS))
    # the deposit Eulers on the GHH composite
    xp_F = [_ghh(R[j]["C_F"], R[j]["N_F"], cal["chi_F"], frisch_F, v.X_FLOOR_F, v.X_EPS)
            for j in range(nreg)]
    rp_D = [(1.0 + v.rdep_D) * v.P_CES_D / ces_price(R[j]["p"], cal, "D") - 1.0
            for j in range(nreg)]
    rp_F = [(1.0 + v.rdep_F) * v.P_CES_F / ces_price(R[j]["p"], cal, "F") - 1.0
            for j in range(nreg)]
    beff_D = 1.0 / (1.0 + cal["r_dep_D_target"])
    beff_F = 1.0 / (1.0 + cal["r_dep_F_target"])
    E_mu_D = _expect(v.wgt, [(1.0 + rp_D[j]) * v.XN_D[j] ** (-sigD) for j in range(nreg)])
    E_mu_F = _expect(v.wgt, [(1.0 + rp_F[j]) * xp_F[j] ** (-sigF) for j in range(nreg)])
    xC_D = float(_ghh(v.C_D, v.N_D, cal["chi_D"], frisch_D, v.X_FLOOR_D, v.X_EPS))
    xC_F = float(_ghh(v.C_F, v.N_F, cal["chi_F"], frisch_F, v.X_FLOOR_F, v.X_EPS))
    v.euler_D = xC_D ** (-sigD) / (beff_D * E_mu_D) - 1.0
    v.euler_F = xC_F ** (-sigF) / (beff_F * E_mu_F) - 1.0
    IM_D = import_demand(np.array([v.p]), np.array([v.C_D]), np.array([v.P_CES_D]), cal, "D")[0]
    IM_F = import_demand(np.array([v.p]), np.array([v.C_F]), np.array([v.P_CES_F]), cal, "F")[0]
    NX_D, NX_F = trade_balance(np.array([v.p]), np.array([IM_D]), np.array([IM_F]), cal)
    v.NX_D = NX_D[0]
    # the F goods market, a Walras-redundant check (never a residual)
    v.goods_F = ((v.Y_F - v.P_CES_F * v.C_F - v.I_F - NX_F[0] - cal["G_F"])
                 / ss["ss_firm_F"]["Y_ss"])
    # E[p'] for deposit-UIP
    v.Ep_next = _expect(v.wgt, [R[j]["p"] for j in range(nreg)])


def _residual_vector(v, cal, ss):
    # THE 14 EQUILIBRIUM CONDITIONS, EACH IN O(1) UNITS.
    frisch_D, frisch_F = cal["frisch_D"], cal["frisch_F"]
    # deposit-UIP; union_nominal_rate (rdep_D = rdep_F) is only a falsification test
    if cal.get("union_nominal_rate", False):
        uip = (v.rdep_D - v.rdep_F) / (1.0 + cal["r_dep_D_target"])
    else:
        uip = ((1.0 + v.rdep_D) - (1.0 + v.rdep_F) * v.Ep_next / v.p
               + cal.get("kappa_nfa", 0.0) * v.nfa_dep_D / ss["ss_firm_D"]["Y_ss"])
    return np.array([
        v.cap_eul_D,  # 1 cap Euler D -> Kp_D
        v.cap_eul_F,  # 2 cap Euler F -> Kp_F
        (cal["chi_D"] * v.N_D ** (1 / frisch_D) - v.w_D / v.P_CES_D)  # 3 lab_D -> N_D
        / (v.w_D / v.P_CES_D),
        (cal["chi_F"] * v.N_F ** (1 / frisch_F) - v.w_F / v.P_CES_F)  # 4 lab_F -> N_F
        / (v.w_F / v.P_CES_F),
        v.euler_D,  # 5 D deposit Euler
        uip,  # 6 deposit-UIP
        (v.Y_D - v.P_CES_D * v.C_D - v.I_D - v.NX_D - cal["G_D"])  # 7 goods_D -> p
        / ss["ss_firm_D"]["Y_ss"],
        v.res_bondD,  # 8 D-bank D-bond -> Q_bD
        (v.E_Om_payD_F - v.dmd_F * v.Q_bD * v.adj_F) / v.E_Om_F,  # 9 F-bank D-bond -> b_DF
        v.euler_F,  # 10 F deposit Euler
        (v.save_union - v.dep_union) / ((1.0 + v.sz) * v.bkD["Dep_supply_ss"]),  # 11 union clearing
        (v.E_Om_payF - v.dmd_F_home * v.Q_bF) / v.E_Om_F,  # 12 F-bank F-bond -> Q_bF
        (v.E_Om_payF_D - v.dmd_F_for * v.Q_bF * v.adj_D) / v.E_Om_D,  # 13 D-bank F-bond -> b_FD
        v.res_tpi,  # 14 TPI rule -> m
    ])


def _outputs(v, cal):
    # EVERY OBJECT STORED AS A RULE OR READ BY THE REPORTING.
    lKD, lbDD, lbFF = cal["lambda_K_D"], cal["lambda_bD_D"], cal["lambda_bF_F"]
    return dict(E_payD=v.E_payD, E_payD_nodef=v.E_payD_nodef, E_Om_payD=v.E_Om_payD,
                lam_bD_mu_D=lbDD * v.mu_D,
                # F-side legs, for the D-F spread decomposition
                E_Om_F=v.E_Om_F, E_payF=v.E_payF, E_Om_payF=v.E_Om_payF,
                lam_bF_mu_F=lbFF * v.mu_F, mu_F_=v.mu_F,
                mu_D=v.mu_D, mu_F=v.mu_F, n_D=v.n_D, n_F=v.n_F, Y_D=v.Y_D, I_D=v.I_D,
                I_F=v.I_F, alpha_D=v.alpha_D_new, alpha_F=v.alpha_F_new,
                Q_bF=v.Q_bF, C_D=v.C_D, C_F=v.C_F,
                b_F_D_new=v.b_F_D_new, b_F_F_new=v.b_F_F_new,
                b_D_D_new=v.b_D_D_new, b_D_F_new=v.b_D_F_new, nfa_dep_D=v.nfa_dep_D,
                W_D=v.W_D, W_F=v.W_F, Wp_D=v.Wp_D, Vp_dep=v.Vp_dep, A_D=v.A_D, A_F=v.A_F,
                Q_bD=v.Q_bD, Bp_tot=v.Bp_D, dep_union=v.dep_union,
                save_union=v.save_union, B_D=v.B_D,
                inc_D=v.inc_D, inc_F=v.inc_F, w_D=v.w_D, dep_D=v.dep_D, dep_F=v.dep_F,
                Pp_D=v.Pp_D, Pp_F=v.Pp_F, Bp_D=v.Bp_D, slack_D=v.slack_D, slack_F=v.slack_F,
                # the TPI book and its flows
                m_cb=v.m_cb, x_cb=v.x_cb, M_cb_new=v.M_cb_new, Z_cb=v.Z_cb, Op_cb=v.Op_cb, Pi_cb=v.Pi_cb,
                Tax_F=v.Tax_F, goods_F=v.goods_F, Q_floor=v.Q_floor, foc_D=v.foc_D,
                # accounting legs for the output decomposition
                N_D=v.N_D, Kap_prod_D=v.K_D, Z_D=v.Z_D, Kp_D=v.Kp_D, P_CES_D=v.P_CES_D,
                E_Om_D=v.E_Om_D, r_wc_D=v.r_wc_D, wedge_sp_D=lKD * v.mu_D / v.E_Om_D,
                rdep_D=v.rdep_D, rdep_F=v.rdep_F, Div_D=v.Div_D, Tax_D=v.Tax_D, p=v.p,
                Y_F=v.Y_F, r_wc_F=v.r_wc_F, L_wc_D=v.L_wc_D, L_wc_F=v.L_wc_F,
                w_F=v.w_F, P_CES_F=v.P_CES_F,
                # diagnostics
                euler_F_resid=v.euler_F, euler_D_resid=v.euler_D,
                C_D_terms=(v.P_D / v.P_CES_D, v.inc_D, v.A_D))


def point_residuals(S, d, x, cont, cal, ss, sproc, n_gh=7, no_default=False):
    # THE PERIOD MAP AT ONE POINT: 14 RESIDUALS AND THE OBJECTS STORED AS RULES.
    v = _read_state_and_regime(S, d, x, cont, cal, ss)
    _production(v, cal)
    _government(v, cal, ss)
    _balance_sheets(v, cont, cal)
    _continuation(v, cont, cal, sproc, n_gh, no_default)
    _discount_kernels(v, cont, cal)
    _incentive_constraint(v, cal)
    _bond_demands(v, cal)
    _tpi_rule(v, cal, ss)
    _wages_and_dividends(v, cal)
    _households(v, cal, ss)
    return _residual_vector(v, cal, ss), _outputs(v, cal)
