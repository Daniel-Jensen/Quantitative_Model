# **Two-country nonlinear transition-path solver (7T stacked market clearing).**
# Unknowns [N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p]; residuals: 2 bank-IC
# complementarities (occasionally-binding μ, Fischer-Burmeister), 2 labour,
# 2 deposit, goods-market D (pins p). Goods-market F and the current account
# are Walras-redundant (monitored, not imposed). Government debt is endogenous
# inside every residual: bonds priced from marginal conditions, debt forward-
# integrated with the Bohn tax, banks clear against the true end-of-period
# stock. Only D is default-risky (Bocola): def_price_D = exogenous PRICED
# default probability path π_t (into Q); def_real_D = REALIZED default.
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
                   def_price_D=None,
                   def_real_D=None,
                   init=None, risk_D=None):
    # **Full inner economy given the 7T guesses: firms → capital → banks → govt → households → trade.**
    T = len(p_path)
    if init is None:
        init = {}

    firm_D = solve_firm_path(N_D, Kap_D, Z_D_path, cal, country="D")
    firm_F = solve_firm_path(N_F, Kap_F, Z_F_path, cal, country="F")

    # quality0 < 1 only in default branches (GK capital-quality loss)
    cap_D = solve_capital_path(Kap_D, init.get("Kap_lag_D", ss["Kap_D_ss"]),
                               init.get("Q_lag_D", 1.0), firm_D["mpk"], cal, country="D",
                               quality0=init.get("quality0_D", 1.0))
    cap_F = solve_capital_path(Kap_F, init.get("Kap_lag_F", ss["Kap_F_ss"]),
                               init.get("Q_lag_F", 1.0), firm_F["mpk"], cal, country="F",
                               quality0=init.get("quality0_F", 1.0))

    P_CES_D = ces_price(p_path, cal, "D")
    P_CES_F = ces_price(p_path, cal, "F")

    # bank backward: prices, multipliers, cross-border FOC holdings
    bwd = bank_backward(
        cap_D["rk"], cap_F["rk"], rdep_D, rdep_F, p_path,
        cal, ss["ss_bank_D"], ss["ss_bank_F"],
        def_price_D=def_price_D, risk_D=risk_D,
    )

    # Working-capital wedge: firms pre-finance ζ×wage-bill at r_wc = rdep(-1)+λμ/Ω̃.
    # Lowers the received wage (spread→output channel); financing income goes to
    # dividends below (intra-period, never on the bank balance sheet). ζ=0 nests off.
    rdep_prev_D = (init["bank_D"]["rdep_prev"] if "bank_D" in init
                   else cal["r_dep_D_target"])
    rdep_prev_F = (init["bank_F"]["rdep_prev"] if "bank_F" in init
                   else cal["r_dep_F_target"])
    zeta_D = cal.get("zeta_wc_D", 0.0)
    zeta_F = cal.get("zeta_wc_F", 0.0)
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

    # government: forward-integrate debt (recap_path = default-branch bailout outlays)
    gov_D = govt_transition(cal, ss["gs_D"], bwd["Q_bD"], def_real_D, "D",
                            b_gov0=init.get("b_gov0_D"),
                            b_anchor=init.get("b_anchor_D"),
                            recap_path=init.get("recap_D_path"))
    gov_F = govt_transition(cal, ss["gs_F"], bwd["Q_bF"], None, "F",
                            b_gov0=init.get("b_gov0_F"))

    # bond clearing against the true end-of-period stock (domestic bank = residual holder)
    b_D_D_path = gov_D["b_gov_eop"] - bwd["b_D_F"]
    b_F_F_path = gov_F["b_gov_eop"] - bwd["b_F_D"]
    if np.any(b_D_D_path <= 0) or np.any(b_F_F_path <= 0):
        raise RuntimeError("Domestic bond holdings non-positive: cross-border "
                           "FOC holdings exceed outstanding government stock.")

    # bank forward: net worth, dividends, deposit supply from REALIZED returns (+ recap)
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
    sigma_D = cal["sigma_D"]; sigma_F  = cal["sigma_F"]

    vN_D = chi_D * N_D ** (1 + 1 / frisch_D) / (1 + 1 / frisch_D)   # GHH labour disutility
    vN_F = chi_F * N_F ** (1 + 1 / frisch_F) / (1 + 1 / frisch_F)

    w_real_D = firm_D["w"] / P_CES_D
    w_real_F = firm_F["w"] / P_CES_F
    e_D = ss["e_D"];  e_F = ss["e_F"]

    # individual income (T, n_e) in composite-good units
    y_D_path = (w_real_D * N_D)[:, None] * e_D[None, :] + ((Div_D - gov_D["Tax"]) / P_CES_D)[:, None]
    y_F_path = (w_real_F * N_F)[:, None] * e_F[None, :] + ((Div_F - gov_F["Tax"]) / P_CES_F)[:, None]

    # Fisher real deposit returns (rate predetermined at t-1); period -1 anchors from init or SS
    P_CES_D_ext = np.concatenate([[init.get("P_lag_D", 1.0)], P_CES_D, [1.0]])
    P_CES_F_ext = np.concatenate([[init.get("P_lag_F", 1.0)], P_CES_F, [1.0]])
    rdep_D_full = np.concatenate([[rdep_prev_D], rdep_D])
    rdep_F_full = np.concatenate([[rdep_prev_F], rdep_F])
    r_D_path = (1.0 + rdep_D_full) * P_CES_D_ext[:-1] / P_CES_D_ext[1:] - 1.0
    r_F_path = (1.0 + rdep_F_full) * P_CES_F_ext[:-1] / P_CES_F_ext[1:] - 1.0

    use_fast = bool(cal.get("use_numba", True))
    c_D_path, a_pol_D_path = solve_backward_transition(
        ss["a_grid_D"], ss["Pi_D"], r_D_path, y_D_path, ss["c_D_ss"],
        ss["beta_D_ss"], sigma_D, cal["a_min_D"], vN_path=vN_D,
        use_fast=use_fast,
    )
    c_F_path, a_pol_F_path = solve_backward_transition(
        ss["a_grid_F"], ss["Pi_F"], r_F_path, y_F_path, ss["c_F_ss"],
        ss["beta_F_ss"], sigma_F, cal["a_min_F"], vN_path=vN_F,
        use_fast=use_fast,
    )

    # distribution forward; C_t on start-of-period dist, A_t on end-of-period; D_start[t]
    # (the dist entering period t) stored so a default branch can launch from any date
    D_D = init.get("D_D", ss["D_D_ss"])
    D_F = init.get("D_F", ss["D_F_ss"])
    A_D_path, C_D_path, D_start_D = forward_paths(
        D_D, a_pol_D_path, c_D_path, ss["a_grid_D"], ss["Pi_D"], use_fast)
    A_F_path, C_F_path, D_start_F = forward_paths(
        D_F, a_pol_F_path, c_F_path, ss["a_grid_F"], ss["Pi_F"], use_fast)

    IM_D = import_demand(p_path, C_D_path, P_CES_D, cal, "D")
    IM_F = import_demand(p_path, C_F_path, P_CES_F, cal, "F")
    NX_D_path, NX_F_path = trade_balance(p_path, IM_D, IM_F)

    return dict(
        firm_D=firm_D, firm_F=firm_F,
        cap_D=cap_D, cap_F=cap_F,
        bk=bk,
        gov_D=gov_D, gov_F=gov_F,
        Tax_D=gov_D["Tax"], Tax_F=gov_F["Tax"],
        Div_D=Div_D, Div_F=Div_F,
        P_CES_D=P_CES_D, P_CES_F=P_CES_F,
        y_D=y_D_path, y_F=y_F_path,
        vN_D=vN_D, vN_F=vN_F,
        c_D=c_D_path, a_pol_D=a_pol_D_path,
        c_F=c_F_path, a_pol_F=a_pol_F_path,
        A_D=A_D_path, C_D=C_D_path,
        A_F=A_F_path, C_F=C_F_path,
        NX_D=NX_D_path, NX_F=NX_F_path,
        IM_D=IM_D, IM_F=IM_F,
        D_start_D=D_start_D, D_start_F=D_start_F,
        r_wc_D=r_wc_D, r_wc_F=r_wc_F,
    )


def make_residual(spec, verbose=False):
    # **Build the 7T market-clearing residual F(y) from a picklable spec (Jacobian workers reuse it).**
    ss  = spec["ss"];   cal = spec["cal"]
    Z_D_path = spec["Z_D_path"];  Z_F_path = spec["Z_F_path"]
    def_price_D = spec.get("def_price_D")
    def_real_D  = spec.get("def_real_D")
    init = spec.get("init");  risk_D = spec.get("risk_D")

    T = cal["T"]
    ncalls  = [0]
    chi_D   = cal["chi_D"];   frisch_D = cal["frisch_D"]
    chi_F   = cal["chi_F"];   frisch_F = cal["frisch_F"]
    n_ss_D  = ss["ss_bank_D"]["n_ss"]
    n_ss_F  = ss["ss_bank_F"]["n_ss"]
    mu_ss_D = ss["ss_bank_D"]["mu_ss"]
    mu_ss_F = ss["ss_bank_F"]["mu_ss"]
    Kap_D_ss = ss["Kap_D_ss"]
    Kap_F_ss = ss["Kap_F_ss"]
    Y_ss_D  = ss["ss_firm_D"]["Y_ss"]
    G_D     = cal["G_D"]

    def residual(y):
        # **Stacked residual of the 7 market-clearing conditions.**
        ncalls[0] += 1
        N_D, N_F     = y[:T], y[T:2*T]
        Kap_D, Kap_F = y[2*T:3*T], y[3*T:4*T]
        rdep_D, rdep_F = y[4*T:5*T], y[5*T:6*T]
        p_path       = y[6*T:7*T]

        # domain guard: below these fractional powers go NaN. Penalty 10.0 must
        # MATCH the failure paths below (unequal walls bias hybr's gradient).
        if (np.any(p_path <= 0.05) or np.any(N_D <= 0.01) or np.any(N_F <= 0.01)
                or np.any(Kap_D <= 0.1) or np.any(Kap_F <= 0.1)):
            return np.full(7 * T, 10.0)

        try:
            out = _inner_economy(
                N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                Z_D_path, Z_F_path, ss, cal,
                def_price_D=def_price_D,
                def_real_D=def_real_D,
                init=init, risk_D=risk_D,
            )
        except (ValueError, RuntimeError, FloatingPointError) as e:
            if verbose:
                print(f"  call {ncalls[0]:3d}: FAILED ({e}); penalising")
            return np.full(7 * T, 10.0)

        bk     = out["bk"]
        firm_D = out["firm_D"]
        firm_F = out["firm_F"]

        # occasionally-binding IC (Bocola): 0 ≤ μ ⊥ slack ≥ 0, imposed via the
        # Fischer-Burmeister function φ(a,b)=a+b−√(a²+b²). slack = α(n−n_IC)
        # = αn − λ·assets. Arguments are scaled independently (FB's zero set is
        # scaling-invariant); at the SS (μ=μ_ss>0, slack=0) φ≡0 and is smooth.
        mu_D_sc = bk["mu_D"] / mu_ss_D
        mu_F_sc = bk["mu_F"] / mu_ss_F
        slack_D = bk["alpha_D"] * (bk["n_D"] - bk["n_IC_D"]) / n_ss_D
        slack_F = bk["alpha_F"] * (bk["n_F"] - bk["n_IC_F"]) / n_ss_F
        cap_D_resid = mu_D_sc + slack_D - np.sqrt(mu_D_sc ** 2 + slack_D ** 2)
        cap_F_resid = mu_F_sc + slack_F - np.sqrt(mu_F_sc ** 2 + slack_F ** 2)

        # labour: GHH static FOC chi·N^(1/frisch) = w/P_CES
        lab_D_supply = chi_D * N_D ** (1 / frisch_D)
        lab_D_demand = firm_D["w"] / out["P_CES_D"]
        lab_D_resid  = (lab_D_supply - lab_D_demand) / (lab_D_demand + 1e-12)
        lab_F_supply = chi_F * N_F ** (1 / frisch_F)
        lab_F_demand = firm_F["w"] / out["P_CES_F"]
        lab_F_resid  = (lab_F_supply - lab_F_demand) / (lab_F_demand + 1e-12)

        # deposits: bank supply is nominal, household A is real → scale by P_CES.
        # pin_rdep (diagnostic, default off): hold the deposit rate fixed at
        # (r_D, r_F) and DROP the two deposit-market clearing conditions,
        # replacing them with rate pins. Nothing then forces household supply =
        # bank demand, so the un-cleared imbalance leaks into the monitored
        # goods_F / current-account identities (Walras-violating by design).
        pin = cal.get("pin_rdep", None)
        if pin is None:
            dep_D_resid = (out["P_CES_D"] * out["A_D"] - bk["Dep_supply_D"]) / Kap_D_ss
            dep_F_resid = (out["P_CES_F"] * out["A_F"] - bk["Dep_supply_F"]) / Kap_F_ss
        else:
            dep_D_resid = rdep_D - pin[0]
            dep_F_resid = rdep_F - pin[1]

        # goods market D (pins p)
        goods_D_resid = (firm_D["Y"] - out["P_CES_D"] * out["C_D"] - out["cap_D"]["I"]
                         - out["NX_D"] - G_D) / Y_ss_D

        resid = np.concatenate([
            cap_D_resid, cap_F_resid,
            lab_D_resid, lab_F_resid,
            dep_D_resid, dep_F_resid,
            goods_D_resid,
        ])
        if not np.all(np.isfinite(resid)):
            return np.full(7 * T, 10.0)

        if verbose:
            walras_D = np.max(np.abs(firm_D["Y"] - out["P_CES_D"] * out["C_D"] - out["cap_D"]["I"] - out["NX_D"] - G_D))
            walras_F = np.max(np.abs(firm_F["Y"] - out["P_CES_F"] * out["C_F"] - out["cap_F"]["I"] - out["NX_F"] - cal["G_F"]))
            print(f"  call {ncalls[0]:3d}: max|resid|={np.max(np.abs(resid)):.3e}"
                  f"  goods_F={walras_F:.3e}  goods_D={walras_D:.3e}")

        return resid

    return residual


def solve_transition(ss, cal, Z_D_path, Z_F_path,
                     def_price_D=None,
                     def_real_D=None,
                     verbose=True, maxiter=300, y0=None,
                     init=None, risk_D=None, jac_cache=None,
                     hybr_factor=100.0, accept_tol=None):
    # **Solve the 7T system: damped Newton (jac_cache-reused) → hybr fallback → Newton polish.**
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
                def_price_D=def_price_D,
                def_real_D=def_real_D,
                init=init, risk_D=risk_D)
    residual = make_residual(spec, verbose=verbose)
    # branch probes pass accept_tol=1e-9; else cal default 1e-10
    if accept_tol is None:
        accept_tol = cal.get("tol_transition", 1e-10)
    jc = jac_cache if jac_cache is not None else {}
    n_jobs = int(cal.get("n_jobs", 0))

    def build_jac(y, F):
        return fd_jacobian(residual, y, F, spec=spec, n_jobs=n_jobs,
                           verbose=verbose)

    y_sol, F_sol, ok = newton_solve(residual, y0, jac_cache=jc,
                                    accept_tol=accept_tol,
                                    build_jac=build_jac, verbose=verbose)
    resid_norm = np.max(np.abs(F_sol))

    if not ok:
        # hybr trust-region fallback for guesses outside the Newton basin, then a
        # Newton polish with a fresh Jacobian (hybr alone plateaus ~5e-11 on xtol)
        if verbose:
            print(f"  newton stalled at max|resid|={resid_norm:.3e}; "
                  "falling back to hybr")
        sol = root(residual, y0, method="hybr",
                   options={"maxfev": max(maxiter * (7 * T + 1), 50000),
                            "factor": hybr_factor})
        F_h = residual(sol.x)
        if np.max(np.abs(F_h)) < resid_norm:
            y_sol, F_sol = sol.x, F_h
            resid_norm = np.max(np.abs(F_h))
        if resid_norm > accept_tol:
            jc.clear()
            y_sol, F_sol, ok = newton_solve(residual, y_sol, F0=F_sol,
                                            jac_cache=jc,
                                            accept_tol=accept_tol,
                                            build_jac=build_jac,
                                            verbose=verbose)
            resid_norm = np.max(np.abs(F_sol))

    if resid_norm > accept_tol:
        raise RuntimeError(f"Transition path did not converge: max|resid|={resid_norm:.3e}")

    N_D, N_F     = y_sol[:T], y_sol[T:2*T]
    Kap_D, Kap_F = y_sol[2*T:3*T], y_sol[3*T:4*T]
    rdep_D, rdep_F = y_sol[4*T:5*T], y_sol[5*T:6*T]
    p_path       = y_sol[6*T:7*T]

    out = _inner_economy(N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                         Z_D_path, Z_F_path, ss, cal,
                         def_price_D=def_price_D,
                         def_real_D=def_real_D,
                         init=init, risk_D=risk_D)

    # complementarity monitor: at an exact FB solution μ ≥ 0 and slack ≥ 0
    # everywhere. Violations mean the Newton accepted a spurious corner.
    mu_min_D = float(np.min(out["bk"]["mu_D"]))
    mu_min_F = float(np.min(out["bk"]["mu_F"]))
    slack_min_D = float(np.min(out["bk"]["alpha_D"]
                               * (out["bk"]["n_D"] - out["bk"]["n_IC_D"])))
    slack_min_F = float(np.min(out["bk"]["alpha_F"]
                               * (out["bk"]["n_F"] - out["bk"]["n_IC_F"])))
    comp_tol = 1e-7
    if min(mu_min_D, mu_min_F) < -comp_tol or min(slack_min_D, slack_min_F) < -comp_tol:
        print(f"  [transition] WARNING: complementarity violated on the solved "
              f"path (min mu D/F={mu_min_D:+.2e}/{mu_min_F:+.2e}, "
              f"min slack D/F={slack_min_D:+.2e}/{slack_min_F:+.2e}).")

    zeros = np.zeros(T)
    return dict(
        slack_min_D=slack_min_D, slack_min_F=slack_min_F,
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
        b_gov_bop_D=out["gov_D"]["b_gov"], b_gov_bop_F=out["gov_F"]["b_gov"],
        coupon_D=out["gov_D"]["coupon"], coupon_F=out["gov_F"]["coupon"],
        net_issuance_D=out["gov_D"]["net_issuance"], net_issuance_F=out["gov_F"]["net_issuance"],
        def_price_D=(def_price_D if def_price_D is not None else zeros),
        def_real_D=(def_real_D if def_real_D is not None else zeros),
        mu_min_D=mu_min_D, mu_min_F=mu_min_F,
        **out["bk"],   # all bank paths (alpha, mu, Omega, n, Q_b, rb, holdings, spreads)
        y_vec=y_sol,   # solved unknowns (warm start for homotopy)
    )


