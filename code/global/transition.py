# TWO-COUNTRY NONLINEAR TRANSITION-PATH SOLVER (7T STACKED MARKET CLEARING).
# Unknowns [N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p]. Residuals: 2 bank-IC
# complementarities (occasionally-binding mu, Fischer-Burmeister), 2 labour,
# union deposit clearing + deposit-UIP (real-rate parity), goods-market D (pins
# p). Goods-market F and the current account are Walras-redundant: monitored,
# never imposed. Capital is predetermined (Bocola eq. 6), so the stock producing
# at t was bought at t-1 and impact output moves through hours alone. Government
# debt is endogenous inside every residual: bonds priced from marginal
# conditions, debt forward-integrated under the Bohn tax, banks clearing against
# the true end-of-period stock. Only D is default-risky.
import numpy as np
from scipy.optimize import root

from solvers import fd_jacobian, newton_solve
from firms import solve_firm_path, markup_ss
from capital import solve_capital_path
from bank import bank_backward, bank_forward
from household import solve_backward_transition
from distribution import forward_paths
from trade import ces_price, import_demand, trade_balance
from government import govt_transition


def _inner_economy(N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                   Z_D_path, Z_F_path, ss, cal,
                   def_price_D=None, def_real_D=None,
                   init=None, risk_D=None,
                   tpi_D=None, s_tpi_D=None):
    # FULL INNER ECONOMY GIVEN THE 7T GUESSES: FIRMS, CAPITAL, BANKS, GOVT, HOUSEHOLDS, TRADE.
    T = len(p_path)
    if init is None:
        init = {}

    # predetermined capital (Bocola eq. 6): the stock producing at t is the one
    # carried INTO t, anchored at the initial stock
    K_init_D = init.get("Kap_lag_D", ss["Kap_D_ss"])
    K_init_F = init.get("Kap_lag_F", ss["Kap_F_ss"])
    Kap_prod_D = np.concatenate(([K_init_D], Kap_D[:-1]))
    Kap_prod_F = np.concatenate(([K_init_F], Kap_F[:-1]))

    firm_D = solve_firm_path(N_D, Kap_prod_D, Z_D_path, cal, country="D")
    firm_F = solve_firm_path(N_F, Kap_prod_F, Z_F_path, cal, country="F")

    cap_D = solve_capital_path(Kap_D, K_init_D, init.get("Q_lag_D", 1.0),
                               firm_D["mpk"], cal, "D", Kap_lag_path=Kap_prod_D)
    cap_F = solve_capital_path(Kap_F, K_init_F, init.get("Q_lag_F", 1.0),
                               firm_F["mpk"], cal, "F", Kap_lag_path=Kap_prod_F)

    P_CES_D = ces_price(p_path, cal, "D")
    P_CES_F = ces_price(p_path, cal, "F")

    bwd = bank_backward(
        cap_D["rk"], cap_F["rk"], rdep_D, rdep_F, p_path,
        cal, ss["ss_bank_D"], ss["ss_bank_F"],
        def_price_D=def_price_D, risk_D=risk_D,
        tpi_D=tpi_D, s_tpi_D=s_tpi_D,
    )
    cb_buy_D_path = bwd["cb_buy_D"]   # realized CB purchases (zeros if TPI off)

    # Working-capital wedge: firms pre-finance zeta x wage bill at
    # r_wc = rdep(-1) + lambda*mu/Omega, which lowers the received wage. The
    # financing income goes to dividends below (intra-period, never on the bank
    # balance sheet). zeta = 0 nests it off.
    rdep_prev_D = (init["bank_D"]["rdep_prev"] if "bank_D" in init
                   else cal["r_dep_D_target"])
    rdep_prev_F = (init["bank_F"]["rdep_prev"] if "bank_F" in init
                   else cal["r_dep_F_target"])
    zeta_D = cal["zeta_wc_D"]
    zeta_F = cal["zeta_wc_F"]
    r_wc_D = (np.concatenate([[rdep_prev_D], rdep_D[:-1]])
              + cal["lambda_K_D"] * bwd["mu_D"] / bwd["Omega_D"])
    r_wc_F = (np.concatenate([[rdep_prev_F], rdep_F[:-1]])
              + cal["lambda_K_F"] * bwd["mu_F"] / bwd["Omega_F"])
    if zeta_D != 0.0:
        firm_D["w"] = firm_D["w"] / (1.0 + zeta_D * r_wc_D)
    if zeta_F != 0.0:
        firm_F["w"] = firm_F["w"] / (1.0 + zeta_F * r_wc_F)
    wc_income_D = zeta_D * r_wc_D * firm_D["w"] * N_D
    wc_income_F = zeta_F * r_wc_F * firm_F["w"] * N_F

    gov_D = govt_transition(cal, ss["gs_D"], bwd["Q_bD"], def_real_D, "D",
                            b_gov0=init.get("b_gov0_D"),
                            b_anchor=init.get("b_anchor_D"),
                            recap_path=init.get("recap_D_path"))
    gov_F = govt_transition(cal, ss["gs_F"], bwd["Q_bF"], None, "F",
                            b_gov0=init.get("b_gov0_F"))

    # CB net cash flow: coupon + surviving continuation value of last period's
    # holding, minus this period's purchase cost. Rebated lump-sum below. The
    # historical TPI-1 audit bug was exactly the absence of this identity —
    # purchases subtracted from bank holdings with nobody recording the cost.
    def_real_D_arr = def_real_D if def_real_D is not None else np.zeros(T)
    surv_cb_D = 1.0 - def_real_D_arr * (1.0 - cal["recovery_rate_D"])
    cb_buy_D_lag = np.concatenate(([init.get("cb_buy_D_lag0", 0.0)], cb_buy_D_path[:-1]))
    rem_cb_D = (cal["delta_b_D"] * surv_cb_D * cb_buy_D_lag
                + bwd["Q_bD"] * surv_cb_D * (1 - cal["delta_b_D"]) * cb_buy_D_lag
                - bwd["Q_bD"] * cb_buy_D_path)

    # bond clearing against the true end-of-period stock: the domestic bank is
    # the residual holder, after the foreign bank's leg and any CB purchase
    b_D_D_path = gov_D["b_gov_eop"] - bwd["b_D_F"] - cb_buy_D_path
    b_F_F_path = gov_F["b_gov_eop"] - bwd["b_F_D"]
    if np.any(b_D_D_path <= 0) or np.any(b_F_F_path <= 0):
        raise RuntimeError("Domestic bond holdings non-positive: cross-border "
                           "FOC holdings (+ CB purchases) exceed outstanding "
                           "government stock.")

    fwd = bank_forward(
        Kap_D, Kap_F, cap_D["Q"], cap_F["Q"],
        cap_D["rk"], cap_F["rk"], rdep_D, rdep_F, p_path,
        b_D_D_path, b_F_F_path, bwd, cal, ss["ss_bank_D"], ss["ss_bank_F"],
        def_real_D=def_real_D,
        init_D=init.get("bank_D"), init_F=init.get("bank_F"),
        Q_bD_lag0=init.get("Q_bD_lag"), Q_bF_lag0=init.get("Q_bF_lag"),
        p_lag0=init.get("p_lag"),
        recap_D=init.get("recap_D_path"),
    )
    bk = {**bwd, **fwd}

    mc_D = markup_ss(cal, "D")
    mc_F = markup_ss(cal, "F")
    Div_D = (1 - mc_D) * firm_D["Y"] + cap_D["cap_profit"] + bk["div_D"] + wc_income_D
    Div_F = (1 - mc_F) * firm_F["Y"] + cap_F["cap_profit"] + bk["div_F"] + wc_income_F

    chi_D   = cal["chi_D"];   frisch_D = cal["frisch_D"]
    chi_F   = cal["chi_F"];   frisch_F = cal["frisch_F"]

    vN_D = chi_D * N_D ** (1 + 1 / frisch_D) / (1 + 1 / frisch_D)   # GHH labour disutility
    vN_F = chi_F * N_F ** (1 + 1 / frisch_F) / (1 + 1 / frisch_F)

    w_real_D = firm_D["w"] / P_CES_D
    w_real_F = firm_F["w"] / P_CES_F

    # CB rebate split by SS-GDP share. rem_cb_D is D-goods denominated. D's own
    # share stays inside D's budget loop; F's share is a genuine CROSS-BORDER
    # real transfer, so it converts via p AND must appear explicitly in both
    # goods residuals (omitting either side is the TPI-1/W-2 class of bug).
    Y_ss_D = ss["ss_firm_D"]["Y_ss"];  Y_ss_F = ss["ss_firm_F"]["Y_ss"]
    share_D = Y_ss_D / (Y_ss_D + Y_ss_F)
    rebate_D = share_D * rem_cb_D
    rebate_F = (1.0 - share_D) * rem_cb_D / p_path

    # individual income (T, n_e) in composite-good units
    y_D_path = ((w_real_D * N_D)[:, None] * ss["e_D"][None, :]
                + ((Div_D - gov_D["Tax"] + rebate_D) / P_CES_D)[:, None])
    y_F_path = ((w_real_F * N_F)[:, None] * ss["e_F"][None, :]
                + ((Div_F - gov_F["Tax"] + rebate_F) / P_CES_F)[:, None])

    # Fisher real deposit returns (the rate paid at t was locked at t-1); the
    # period -1 entry anchors from init or the SS
    P_CES_D_ext = np.concatenate([[init.get("P_lag_D", 1.0)], P_CES_D, [1.0]])
    P_CES_F_ext = np.concatenate([[init.get("P_lag_F", 1.0)], P_CES_F, [1.0]])
    rdep_D_full = np.concatenate([[rdep_prev_D], rdep_D])
    rdep_F_full = np.concatenate([[rdep_prev_F], rdep_F])
    r_D_path = (1.0 + rdep_D_full) * P_CES_D_ext[:-1] / P_CES_D_ext[1:] - 1.0
    r_F_path = (1.0 + rdep_F_full) * P_CES_F_ext[:-1] / P_CES_F_ext[1:] - 1.0

    use_fast = bool(cal["use_numba"])
    c_D_path, a_pol_D_path = solve_backward_transition(
        ss["a_grid_D"], ss["Pi_D"], r_D_path, y_D_path, ss["c_D_ss"],
        ss["beta_D_ss"], cal["sigma_D"], cal["a_min_D"], vN_path=vN_D,
        use_fast=use_fast,
    )
    c_F_path, a_pol_F_path = solve_backward_transition(
        ss["a_grid_F"], ss["Pi_F"], r_F_path, y_F_path, ss["c_F_ss"],
        ss["beta_F_ss"], cal["sigma_F"], cal["a_min_F"], vN_path=vN_F,
        use_fast=use_fast,
    )

    # D_start[t] (the dist entering t) is stored so a branch can launch from any date
    A_D_path, C_D_path, D_start_D = forward_paths(
        init.get("D_D", ss["D_D_ss"]), a_pol_D_path, c_D_path,
        ss["a_grid_D"], ss["Pi_D"], use_fast)
    A_F_path, C_F_path, D_start_F = forward_paths(
        init.get("D_F", ss["D_F_ss"]), a_pol_F_path, c_F_path,
        ss["a_grid_F"], ss["Pi_F"], use_fast)

    IM_D = import_demand(p_path, C_D_path, P_CES_D, cal, "D")
    IM_F = import_demand(p_path, C_F_path, P_CES_F, cal, "F")
    NX_D_path, NX_F_path = trade_balance(p_path, IM_D, IM_F)

    return dict(
        firm_D=firm_D, firm_F=firm_F,
        cap_D=cap_D, cap_F=cap_F,
        bk=bk,
        gov_D=gov_D, gov_F=gov_F,
        Div_D=Div_D, Div_F=Div_F,
        P_CES_D=P_CES_D, P_CES_F=P_CES_F,
        vN_D=vN_D, vN_F=vN_F,
        A_D=A_D_path, C_D=C_D_path,
        A_F=A_F_path, C_F=C_F_path,
        NX_D=NX_D_path, NX_F=NX_F_path,
        D_start_D=D_start_D, D_start_F=D_start_F,
        r_wc_D=r_wc_D, r_wc_F=r_wc_F,
        rem_cb_D=rem_cb_D,
    )


def make_residual(spec, verbose=False):
    # BUILD THE 7T MARKET-CLEARING RESIDUAL F(y) FROM A PICKLABLE SPEC.
    # The spec (not a closure over live objects) is what lets the Jacobian
    # workers rebuild an identical residual after a spawn.
    ss  = spec["ss"];   cal = spec["cal"]
    Z_D_path = spec["Z_D_path"];  Z_F_path = spec["Z_F_path"]
    def_price_D = spec.get("def_price_D")
    def_real_D  = spec.get("def_real_D")
    init = spec.get("init");  risk_D = spec.get("risk_D")
    tpi_D = spec.get("tpi_D");  s_tpi_D = spec.get("s_tpi_D")

    T = cal["T"]
    ncalls  = [0]
    chi_D   = cal["chi_D"];   frisch_D = cal["frisch_D"]
    chi_F   = cal["chi_F"];   frisch_F = cal["frisch_F"]
    n_ss_D  = ss["ss_bank_D"]["n_ss"]
    n_ss_F  = ss["ss_bank_F"]["n_ss"]
    mu_ss_D = ss["ss_bank_D"]["mu_ss"]
    mu_ss_F = ss["ss_bank_F"]["mu_ss"]
    Kap_scale = ss["Kap_D_ss"] + ss["Kap_F_ss"]
    Y_ss_D  = ss["ss_firm_D"]["Y_ss"]
    Y_ss_F  = ss["ss_firm_F"]["Y_ss"]
    share_F_gdp = Y_ss_F / (Y_ss_D + Y_ss_F)   # CB-rebate cross-border share
    G_D     = cal["G_D"]

    def residual(y):
        # STACKED RESIDUAL OF THE 7 MARKET-CLEARING CONDITIONS.
        ncalls[0] += 1
        N_D, N_F     = y[:T], y[T:2*T]
        Kap_D, Kap_F = y[2*T:3*T], y[3*T:4*T]
        rdep_D, rdep_F = y[4*T:5*T], y[5*T:6*T]
        p_path       = y[6*T:7*T]

        # domain guard: below these, fractional powers go NaN. The penalty 10.0
        # must MATCH the failure paths below (unequal walls bias hybr's gradient).
        if (np.any(p_path <= 0.05) or np.any(N_D <= 0.01) or np.any(N_F <= 0.01)
                or np.any(Kap_D <= 0.1) or np.any(Kap_F <= 0.1)):
            return np.full(7 * T, 10.0)

        try:
            out = _inner_economy(
                N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                Z_D_path, Z_F_path, ss, cal,
                def_price_D=def_price_D, def_real_D=def_real_D,
                init=init, risk_D=risk_D,
                tpi_D=tpi_D, s_tpi_D=s_tpi_D,
            )
        except (ValueError, RuntimeError, FloatingPointError) as e:
            if verbose:
                print(f"  call {ncalls[0]:3d}: FAILED ({e}); penalising")
            return np.full(7 * T, 10.0)

        bk     = out["bk"]
        firm_D = out["firm_D"]
        firm_F = out["firm_F"]

        # Occasionally-binding IC (Bocola): 0 <= mu _|_ slack >= 0, imposed via the
        # Fischer-Burmeister function phi(a,b) = a + b - sqrt(a^2+b^2), where
        # slack = alpha(n - n_IC). The two arguments are scaled independently
        # (FB's zero set is scaling-invariant); at the SS (mu = mu_ss > 0,
        # slack = 0) phi is identically 0 and smooth.
        mu_D_sc = bk["mu_D"] / mu_ss_D
        mu_F_sc = bk["mu_F"] / mu_ss_F
        slack_D = bk["alpha_D"] * (bk["n_D"] - bk["n_IC_D"]) / n_ss_D
        slack_F = bk["alpha_F"] * (bk["n_F"] - bk["n_IC_F"]) / n_ss_F
        cap_D_resid = mu_D_sc + slack_D - np.sqrt(mu_D_sc ** 2 + slack_D ** 2)
        cap_F_resid = mu_F_sc + slack_F - np.sqrt(mu_F_sc ** 2 + slack_F ** 2)

        # labour: GHH static FOC chi*N^(1/frisch) = w/P_CES
        lab_D_demand = firm_D["w"] / out["P_CES_D"]
        lab_F_demand = firm_F["w"] / out["P_CES_F"]
        lab_D_resid = (chi_D * N_D ** (1 / frisch_D) - lab_D_demand) / (lab_D_demand + 1e-12)
        lab_F_resid = (chi_F * N_F ** (1 / frisch_F) - lab_F_demand) / (lab_F_demand + 1e-12)

        # Union deposit market. Deposits stay own-good claims at national rates;
        # a frictionless union interbank replaces the two national clearings with
        # (i) one union-wide clearing in D-good units, whose absorption margin is
        # the cross-border deposit position, and (ii) real-rate parity
        # (1+rdep_D,t) = (1+rdep_F,t)*p_{t+1}/p_t — the flexible-price image of a
        # single nominal union rate plus national inflation differentials. Parity
        # makes the interbank pass-through zero-profit, hence no Walras leak.
        dep_union_resid = ((out["P_CES_D"] * out["A_D"] - bk["Dep_supply_D"])
                           + p_path * (out["P_CES_F"] * out["A_F"]
                                       - bk["Dep_supply_F"])) / Kap_scale
        p_next = np.append(p_path[1:], p_path[-1])   # terminal: p flat
        uip_resid = (1.0 + rdep_D) - (1.0 + rdep_F) * p_next / p_path

        # goods market D (pins p); the CB rebate's F-share leaves D's resources
        cb_transfer_D = share_F_gdp * out["rem_cb_D"]
        goods_D_resid = (firm_D["Y"] - out["P_CES_D"] * out["C_D"] - out["cap_D"]["I"]
                         - out["NX_D"] - G_D - cb_transfer_D) / Y_ss_D

        resid = np.concatenate([
            cap_D_resid, cap_F_resid,
            lab_D_resid, lab_F_resid,
            dep_union_resid, uip_resid,
            goods_D_resid,
        ])
        if not np.all(np.isfinite(resid)):
            return np.full(7 * T, 10.0)

        if verbose:
            walras_F = np.max(np.abs(firm_F["Y"] - out["P_CES_F"] * out["C_F"]
                                     - out["cap_F"]["I"] - out["NX_F"] - cal["G_F"]))
            print(f"  call {ncalls[0]:3d}: max|resid|={np.max(np.abs(resid)):.3e}"
                  f"  goods_F={walras_F:.3e}")

        return resid

    return residual


def solve_transition(ss, cal, Z_D_path, Z_F_path,
                     def_price_D=None, def_real_D=None,
                     verbose=True, maxiter=300, y0=None,
                     init=None, risk_D=None, jac_cache=None, accept_tol=None,
                     tpi_D=None, s_tpi_D=None):
    # SOLVE THE 7T SYSTEM: DAMPED NEWTON (jac_cache-REUSED), HYBR FALLBACK, NEWTON POLISH.
    T = cal["T"]
    assert len(Z_D_path) == T and len(Z_F_path) == T

    if y0 is None:
        y0 = np.concatenate([
            np.full(T, 1.0),                    # N_D
            np.full(T, 1.0),                    # N_F
            np.full(T, ss["Kap_D_ss"]),         # Kap_D
            np.full(T, ss["Kap_F_ss"]),         # Kap_F
            np.full(T, cal["r_dep_D_target"]),  # rdep_D
            np.full(T, cal["r_dep_F_target"]),  # rdep_F
            np.full(T, ss["p_ss"]),             # p
        ])

    spec = dict(ss=ss, cal=cal, Z_D_path=Z_D_path, Z_F_path=Z_F_path,
                def_price_D=def_price_D, def_real_D=def_real_D,
                init=init, risk_D=risk_D,
                tpi_D=tpi_D, s_tpi_D=s_tpi_D)
    residual = make_residual(spec, verbose=verbose)
    if accept_tol is None:      # branch probes pass 1e-9; else the cal default
        accept_tol = cal["tol_transition"]
    jc = jac_cache if jac_cache is not None else {}
    n_jobs = int(cal["n_jobs"])

    def build_jac(y, F):
        # FRESH PARALLEL FD JACOBIAN AT THE CURRENT ITERATE.
        return fd_jacobian(residual, y, F, spec=spec, n_jobs=n_jobs, verbose=verbose)

    y_sol, F_sol, ok = newton_solve(residual, y0, jac_cache=jc,
                                    accept_tol=accept_tol,
                                    build_jac=build_jac, verbose=verbose)
    resid_norm = np.max(np.abs(F_sol))

    if not ok:
        # hybr trust region for guesses outside the Newton basin, then a Newton
        # polish with a fresh Jacobian (hybr alone plateaus ~5e-11 on xtol)
        if verbose:
            print(f"  newton stalled at max|resid|={resid_norm:.3e}; falling back to hybr")
        sol = root(residual, y0, method="hybr",
                   options={"maxfev": max(maxiter * (7 * T + 1), 50000), "factor": 100.0})
        F_h = residual(sol.x)
        if np.max(np.abs(F_h)) < resid_norm:
            y_sol, F_sol = sol.x, F_h
            resid_norm = np.max(np.abs(F_h))
        if resid_norm > accept_tol:
            jc.clear()
            y_sol, F_sol, ok = newton_solve(residual, y_sol, F0=F_sol, jac_cache=jc,
                                            accept_tol=accept_tol,
                                            build_jac=build_jac, verbose=verbose)
            resid_norm = np.max(np.abs(F_sol))

    if resid_norm > accept_tol:
        raise RuntimeError(f"Transition path did not converge: max|resid|={resid_norm:.3e}")

    N_D, N_F     = y_sol[:T], y_sol[T:2*T]
    Kap_D, Kap_F = y_sol[2*T:3*T], y_sol[3*T:4*T]
    rdep_D, rdep_F = y_sol[4*T:5*T], y_sol[5*T:6*T]
    p_path       = y_sol[6*T:7*T]

    out = _inner_economy(N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                         Z_D_path, Z_F_path, ss, cal,
                         def_price_D=def_price_D, def_real_D=def_real_D,
                         init=init, risk_D=risk_D,
                         tpi_D=tpi_D, s_tpi_D=s_tpi_D)

    # complementarity monitor: at an exact FB solution mu >= 0 and slack >= 0
    # everywhere, so a violation means Newton accepted a spurious corner
    mu_min_D = float(np.min(out["bk"]["mu_D"]))
    mu_min_F = float(np.min(out["bk"]["mu_F"]))
    slack_min_D = float(np.min(out["bk"]["alpha_D"]
                               * (out["bk"]["n_D"] - out["bk"]["n_IC_D"])))
    slack_min_F = float(np.min(out["bk"]["alpha_F"]
                               * (out["bk"]["n_F"] - out["bk"]["n_IC_F"])))
    if min(mu_min_D, mu_min_F, slack_min_D, slack_min_F) < -1e-7:
        print(f"  [transition] WARNING: complementarity violated on the solved "
              f"path (min mu D/F={mu_min_D:+.2e}/{mu_min_F:+.2e}, "
              f"min slack D/F={slack_min_D:+.2e}/{slack_min_F:+.2e}).")

    zeros = np.zeros(T)
    return dict(
        mu_min_D=mu_min_D, mu_min_F=mu_min_F,
        slack_min_D=slack_min_D, slack_min_F=slack_min_F,
        # cross-border deposit position: D households' net claim on the union
        # interbank in D-goods (the union absorption margin, identically 0 pre-union)
        nfa_dep_D=out["P_CES_D"] * out["A_D"] - out["bk"]["Dep_supply_D"],
        N_D=N_D, N_F=N_F, Kap_D=Kap_D, Kap_F=Kap_F,
        rdep_D=rdep_D, rdep_F=rdep_F, p=p_path,
        Z_D=Z_D_path, Z_F=Z_F_path,
        **{k + "_D": v for k, v in out["firm_D"].items()},
        **{k + "_F": v for k, v in out["firm_F"].items()},
        I_D=out["cap_D"]["I"],    Q_D=out["cap_D"]["Q"],    rk_D=out["cap_D"]["rk"],
        I_F=out["cap_F"]["I"],    Q_F=out["cap_F"]["Q"],    rk_F=out["cap_F"]["rk"],
        A_D=out["A_D"], C_D=out["C_D"], A_F=out["A_F"], C_F=out["C_F"],
        vN_D=out["vN_D"], vN_F=out["vN_F"],
        r_wc_D=out["r_wc_D"], r_wc_F=out["r_wc_F"],
        D_start_D=out["D_start_D"], D_start_F=out["D_start_F"],
        NX_D=out["NX_D"], NX_F=out["NX_F"],
        Div_D=out["Div_D"], Div_F=out["Div_F"],
        Tax_D=out["gov_D"]["Tax"], Tax_F=out["gov_F"]["Tax"],
        P_CES_D=out["P_CES_D"], P_CES_F=out["P_CES_F"],
        # endogenous debt (end-of-period = bank-held stock)
        b_gov_D=out["gov_D"]["b_gov_eop"], b_gov_F=out["gov_F"]["b_gov_eop"],
        coupon_D=out["gov_D"]["coupon"], coupon_F=out["gov_F"]["coupon"],
        net_issuance_D=out["gov_D"]["net_issuance"],
        net_issuance_F=out["gov_F"]["net_issuance"],
        def_price_D=(def_price_D if def_price_D is not None else zeros),
        def_real_D=(def_real_D if def_real_D is not None else zeros),
        rem_cb_D=out["rem_cb_D"],
        s_tpi_D=(s_tpi_D if s_tpi_D is not None else zeros),
        **out["bk"],   # every bank path: alpha, mu, Omega, n, Q_b, rb, holdings,
                       # spreads, cb_buy_D, Q_bD_free, Q_floor_D
        y_vec=y_sol,   # solved unknowns (warm start for homotopy)
    )


def market_residuals(out, cal, ss=None):
    # MAX ABSOLUTE RESIDUALS OF THE IMPOSED AND WALRAS-REDUNDANT MARKETS.
    # cap_* is the complementarity product mu*slack (0 at an exact solution).
    # Pass ss whenever out carries a nonzero rem_cb_D (TPI active) so the CB
    # rebate's cross-border transfer enters both goods residuals exactly as it
    # does in the imposed goods_D_resid; omitting it shows a leak that isn't there.
    cb_transfer_D = 0.0
    cb_transfer_F = 0.0
    if ss is not None and "rem_cb_D" in out:
        Y_ss_D = ss["ss_firm_D"]["Y_ss"];  Y_ss_F = ss["ss_firm_F"]["Y_ss"]
        share_F = Y_ss_F / (Y_ss_D + Y_ss_F)
        rem_cb_D = np.asarray(out["rem_cb_D"])
        cb_transfer_D = share_F * rem_cb_D
        cb_transfer_F = share_F * rem_cb_D / np.asarray(out["p"])

    goods_D = np.max(np.abs(out["Y_D"] - out["P_CES_D"] * out["C_D"] - out["I_D"]
                            - out["NX_D"] - cal["G_D"] - cb_transfer_D))
    goods_F = np.max(np.abs(out["Y_F"] - out["P_CES_F"] * out["C_F"] - out["I_F"]
                            - out["NX_F"] - cal["G_F"] + cb_transfer_F))
    dep_union = np.max(np.abs((out["P_CES_D"] * out["A_D"] - out["Dep_supply_D"])
                              + out["p"] * (out["P_CES_F"] * out["A_F"]
                                            - out["Dep_supply_F"])))
    p_next = np.append(np.asarray(out["p"])[1:], np.asarray(out["p"])[-1])
    uip = np.max(np.abs((1.0 + np.asarray(out["rdep_D"]))
                        - (1.0 + np.asarray(out["rdep_F"])) * p_next / out["p"]))
    slack_D = np.asarray(out["alpha_D"]) * (np.asarray(out["n_D"]) - np.asarray(out["n_IC_D"]))
    slack_F = np.asarray(out["alpha_F"]) * (np.asarray(out["n_F"]) - np.asarray(out["n_IC_F"]))

    return dict(
        goods_D=goods_D, goods_F=goods_F, dep_union=dep_union, uip=uip,
        cap_D=np.max(np.abs(np.asarray(out["mu_D"]) * slack_D)),
        cap_F=np.max(np.abs(np.asarray(out["mu_F"]) * slack_F)),
        slack_min_D=float(np.min(slack_D)), slack_min_F=float(np.min(slack_F)),
        mu_min_D=float(np.min(out["mu_D"])), mu_min_F=float(np.min(out["mu_F"])),
        nfa_dep_D=float(np.max(np.abs(out["nfa_dep_D"]))) if "nfa_dep_D" in out else 0.0,
    )
