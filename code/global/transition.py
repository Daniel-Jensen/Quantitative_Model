"""Two-country nonlinear transition-path solver.

Outer unknowns (7T total):
  [N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p]

Seven market-clearing residuals per period:
  1. Capital market D: (n_IC_D − n_ACCUM_D) / n_ss_D          → pins Kap_D
  2. Capital market F: (n_IC_F − n_ACCUM_F) / n_ss_F          → pins Kap_F
  3. Labour market D:  (chi_D·N_D^(1/frisch) − w_D/P_CES_D) / (w_D/P_CES_D) → pins N_D
  4. Labour market F:  (chi_F·N_F^(1/frisch) − w_F/P_CES_F) / (w_F/P_CES_F) → pins N_F
  5. Deposit market D: (P_CES_D·A_D − Dep_supply_D) / Kap_D_ss → pins rdep_D
  6. Deposit market F: (P_CES_F·A_F − Dep_supply_F) / Kap_F_ss → pins rdep_F
  7. Goods market D:   (Y_D − P_CES_D·C_D − I_D − NX_D − G_D) / Y_D_ss → pins p

Walras-redundant (diagnostics only, not imposed):
  • Goods market F
  • Current account (external account D)

Government debt is ENDOGENOUS inside every residual evaluation (no outer
debt loop): bank_backward prices bonds from marginal conditions alone, the
debt stock is then forward-integrated from the budget identity with the Bohn
tax (government.govt_transition), and banks clear the bond market against the
true end-of-period stock:
  b_D_D[t] = b_gov_D_eop[t] − b_D_F[t]   (D-bank is the marginal/residual holder)
  b_F_F[t] = b_gov_F_eop[t] − b_F_D[t]
This keeps Walras exact when beliefs move the debt stock — new issuance at
depressed prices must be absorbed by constrained banks (the Cole-Kehoe /
Bocola amplification).

PRICED vs REALIZED default (Bocola 2016 pass-through experiment):
  def_price_* : probability of default priced into Q and expected-return FOCs.
  def_real_*  : realized haircut path hitting bond returns and government flows.
  Risk-only baseline: def_price = sunspot in the CK crisis zone, def_real = 0.

Household income (GHH, all in composite-good units):
  y_t(e) = (w_t/P_CES_t)·N_t·e + (Div_t − Tax_t)/P_CES_t
  where N_t is aggregate employment and e is the idiosyncratic productivity.
  The GHH labour-disutility vN_t = chi·N_t^(1+1/frisch)/(1+1/frisch) is
  subtracted from effective income in the EGM composite (passed via vN_path).
"""
import numpy as np
from scipy.optimize import root

from solvers import fd_jacobian, newton_solve
from firms import solve_firm_path, markup_ss
from capital import solve_capital_path
from bank import bank_backward, bank_forward
from household import solve_backward_transition
from distribution import forward_paths
from trade import ces_price, import_demand, trade_balance
from government import govt_transition, ck_default_prob


def _inner_economy(N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                   Z_D_path, Z_F_path, ss, cal,
                   def_price_D=None, def_price_F=None,
                   def_real_D=None, def_real_F=None,
                   init=None, risk_D=None):
    
    #GIVEN GUESSES CALCULATE EVERYTHING THATS INNER ECONOMY

    T = len(p_path)
    if init is None:
        init = {}

    # ── Firms
    firm_D = solve_firm_path(N_D, Kap_D, Z_D_path, cal, country="D")
    firm_F = solve_firm_path(N_F, Kap_F, Z_F_path, cal, country="F")

    # ── Capital (quality0 < 1 only in default branches: GK capital-quality loss)
    cap_D = solve_capital_path(Kap_D, init.get("Kap_lag_D", ss["Kap_D_ss"]),
                               init.get("Q_lag_D", 1.0), firm_D["mpk"], cal, country="D",
                               quality0=init.get("quality0_D", 1.0))
    cap_F = solve_capital_path(Kap_F, init.get("Kap_lag_F", ss["Kap_F_ss"]),
                               init.get("Q_lag_F", 1.0), firm_F["mpk"], cal, country="F",
                               quality0=init.get("quality0_F", 1.0))

    # ── Trade / CES (ces_price broadcasts over the whole path)
    P_CES_D = ces_price(p_path, cal, "D")
    P_CES_F = ces_price(p_path, cal, "F")

    # ── Bank backward pass: prices and cross-border FOC holdings
    bwd = bank_backward(
        cap_D["rk"], cap_F["rk"], rdep_D, rdep_F, p_path,
        cal, ss["ss_bank_D"], ss["ss_bank_F"],
        def_price_D=def_price_D, def_price_F=def_price_F, risk_D=risk_D,
    )

    # ── Working-capital wedge (Neumeyer-Perri) ───────────────────────────────
    # Firms pre-finance zeta_wc of the wage bill intra-period at
    # r_wc = rdep(−1) + λμ/Ω̃ (single-λ IC wedge, predetermined funding leg).
    # The wedge lowers the wage the worker receives — this is the
    # spread→output channel; the financing income is distributed to
    # households as dividends below (intra-period, never on the bank balance
    # sheet — Bocola-lite approximation, see calibration.py).  zeta = 0
    # reproduces the wedge-free model exactly.
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

    # ── Government: forward-integrate debt with Bohn tax (realized flows) ────
    # recap paths (default-branch bank recapitalization) are government
    # outlays here and equity injections in bank_forward below.
    gov_D = govt_transition(cal, ss["gs_D"], bwd["Q_bD"], def_real_D, "D",
                            b_gov0=init.get("b_gov0_D"),
                            b_anchor=init.get("b_anchor_D"),
                            recap_path=init.get("recap_D_path"))
    gov_F = govt_transition(cal, ss["gs_F"], bwd["Q_bF"], def_real_F, "F",
                            b_gov0=init.get("b_gov0_F"),
                            b_anchor=init.get("b_anchor_F"),
                            recap_path=init.get("recap_F_path"))

    # ── Bond market clearing against the TRUE end-of-period stock ────────────
    b_D_D_path = gov_D["b_gov_eop"] - bwd["b_D_F"]   # D-bank residual holder
    b_F_F_path = gov_F["b_gov_eop"] - bwd["b_F_D"]   # F-bank residual holder
    if np.any(b_D_D_path <= 0) or np.any(b_F_F_path <= 0):
        raise RuntimeError("Domestic bond holdings non-positive: cross-border "
                           "FOC holdings exceed outstanding government stock.")

    # ── Bank forward pass: net worth, dividends, deposit supply ──────────────
    fwd = bank_forward(
        Kap_D, Kap_F, cap_D["Q"], cap_F["Q"],
        cap_D["rk"], cap_F["rk"], rdep_D, rdep_F, p_path,
        b_D_D_path, b_F_F_path, bwd, cal, ss["ss_bank_D"], ss["ss_bank_F"],
        def_real_D=def_real_D, def_real_F=def_real_F,
        init_D=init.get("bank_D"), init_F=init.get("bank_F"),
        Q_bD_lag0=init.get("Q_bD_lag"), Q_bF_lag0=init.get("Q_bF_lag"),
        p_lag0=init.get("p_lag"),
        recap_D=init.get("recap_D_path"), recap_F=init.get("recap_F_path"),
    )
    bk = {**bwd, **fwd}

    # Why we have both bank forward and backward: backward is used to get the expected returns and the FOCs, 
    # while forward is used to get the actual realized flows and net worth of the banks. 
    # The backward pass uses the expected default probabilities to price the bonds and determine the optimal holdings, while the forward pass uses the realized default probabilities to update the banks' balance sheets and dividends.


    # ── Calculate bunch of things
    mc_D = markup_ss(cal, "D")
    mc_F = markup_ss(cal, "F")
    # wc_income: working-capital financing income, passed through to
    # households within the period (see wedge block above)
    Div_D = (1 - mc_D) * firm_D["Y"] + cap_D["cap_profit"] + bk["div_D"] + wc_income_D
    Div_F = (1 - mc_F) * firm_F["Y"] + cap_F["cap_profit"] + bk["div_F"] + wc_income_F

    chi_D   = cal["chi_D"];   frisch_D = cal["frisch_D"]
    chi_F   = cal["chi_F"];   frisch_F = cal["frisch_F"]
    sigma_D = cal["sigma_D"]; sigma_F  = cal["sigma_F"]

    vN_D = chi_D * N_D ** (1 + 1 / frisch_D) / (1 + 1 / frisch_D)
    vN_F = chi_F * N_F ** (1 + 1 / frisch_F) / (1 + 1 / frisch_F)

    w_real_D = firm_D["w"] / P_CES_D     # real wage in D-consumption units
    w_real_F = firm_F["w"] / P_CES_F
    e_D = ss["e_D"];  e_F = ss["e_F"]

    # y_path shape (T, n_e): individual income in composite-good units.
    y_D_path = (w_real_D * N_D)[:, None] * e_D[None, :] + ((Div_D - gov_D["Tax"]) / P_CES_D)[:, None]
    y_F_path = (w_real_F * N_F)[:, None] * e_F[None, :] + ((Div_F - gov_F["Tax"]) / P_CES_F)[:, None]

    # ── Household EGM backward ────────────────────────────────────────────────
    # Fisher-equation real returns: r_real[t] = (1+rdep[t-1]) * P_CES[t-1]/P_CES[t] - 1
    # Period -1 anchors (rate locked pre-path, price level pre-path) come from
    # init when the path starts mid-crisis; SS otherwise.  (rdep_prev_D/F are
    # defined in the working-capital block above.)
    P_CES_D_ext = np.concatenate([[init.get("P_lag_D", 1.0)], P_CES_D, [1.0]])
    P_CES_F_ext = np.concatenate([[init.get("P_lag_F", 1.0)], P_CES_F, [1.0]])
    rdep_D_full = np.concatenate([[rdep_prev_D], rdep_D])  # length T+1
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

    # ── Distribution forward ──────────────────────────────────────────────────
    # D_start[t] is the start-of-period distribution entering period t;
    # D_start[0] is the initial condition.  Needed to launch default branches
    # mid-path.  C_t aggregates the start-of-period distribution, A_t the
    # end-of-period one (forward_paths preserves this timing convention).
    D_D = init.get("D_D", ss["D_D_ss"])
    D_F = init.get("D_F", ss["D_F_ss"])
    A_D_path, C_D_path, D_start_D = forward_paths(
        D_D, a_pol_D_path, c_D_path, ss["a_grid_D"], ss["Pi_D"], use_fast)
    A_F_path, C_F_path, D_start_F = forward_paths(
        D_F, a_pol_F_path, c_F_path, ss["a_grid_F"], ss["Pi_F"], use_fast)

    # ── Trade (import_demand broadcasts over the whole path) ────────────────
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
    """Market-clearing residual y ↦ F(y) of the 7T stacked system, built from
    a PICKLABLE spec so that multiprocessing Jacobian workers can reconstruct
    the exact same function (solvers.fd_jacobian).  The body is the historical
    in-solver closure verbatim: same penalty walls, same normalizations."""
    ss  = spec["ss"];   cal = spec["cal"]
    Z_D_path = spec["Z_D_path"];  Z_F_path = spec["Z_F_path"]
    def_price_D = spec.get("def_price_D");  def_price_F = spec.get("def_price_F")
    def_real_D  = spec.get("def_real_D");   def_real_F  = spec.get("def_real_F")
    init = spec.get("init");  risk_D = spec.get("risk_D")

    T = cal["T"]
    ncalls  = [0]
    chi_D   = cal["chi_D"];   frisch_D = cal["frisch_D"]
    chi_F   = cal["chi_F"];   frisch_F = cal["frisch_F"]
    n_ss_D  = ss["ss_bank_D"]["n_ss"]
    n_ss_F  = ss["ss_bank_F"]["n_ss"]
    Kap_D_ss = ss["Kap_D_ss"]
    Kap_F_ss = ss["Kap_F_ss"]
    Y_ss_D  = ss["ss_firm_D"]["Y_ss"]
    G_D     = cal["G_D"]

    def residual(y):
        #GET ME THE RESIDUALS FOR THE 7T MARKET CLEARING CONDITIONS
        ncalls[0] += 1
        N_D, N_F     = y[:T], y[T:2*T]
        Kap_D, Kap_F = y[2*T:3*T], y[3*T:4*T]
        rdep_D, rdep_F = y[4*T:5*T], y[5*T:6*T]
        p_path       = y[6*T:7*T]

        # Domain guard: below these thresholds fractional powers go NaN.
        # Penalty height 10.0 must MATCH the other failure paths below —
        # unequal walls bias hybr's finite-difference gradient toward the
        # lower wall (it would rather step into NaN territory than out of
        # bounds).
        if (np.any(p_path <= 0.05) or np.any(N_D <= 0.01) or np.any(N_F <= 0.01)
                or np.any(Kap_D <= 0.1) or np.any(Kap_F <= 0.1)):
            return np.full(7 * T, 10.0)

        try:
            # calculate a full inner economy given the outer guess for the 7T unknowns 
            out = _inner_economy(
                N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                Z_D_path, Z_F_path, ss, cal,
                def_price_D=def_price_D, def_price_F=def_price_F,
                def_real_D=def_real_D, def_real_F=def_real_F,
                init=init, risk_D=risk_D,
            )
        except (ValueError, RuntimeError, FloatingPointError) as e:
            if verbose:
                # penalise infeasible guesses (same wall height as the guard)
                print(f"  call {ncalls[0]:3d}: FAILED ({e}); penalising")
            return np.full(7 * T, 10.0)

        bk     = out["bk"]
        firm_D = out["firm_D"]
        firm_F = out["firm_F"]

        # 1. Capital markets
        cap_D_resid = (bk["n_IC_D"] - bk["n_D"]) / n_ss_D
        cap_F_resid = (bk["n_IC_F"] - bk["n_F"]) / n_ss_F

        # 2. Labour markets (GHH static FOC: chi·N^(1/frisch) = w/P_CES)
        lab_D_supply = chi_D * N_D ** (1 / frisch_D)
        lab_D_demand = firm_D["w"] / out["P_CES_D"]
        lab_D_resid  = (lab_D_supply - lab_D_demand) / (lab_D_demand + 1e-12)

        lab_F_supply = chi_F * N_F ** (1 / frisch_F)
        lab_F_demand = firm_F["w"] / out["P_CES_F"]
        lab_F_resid  = (lab_F_supply - lab_F_demand) / (lab_F_demand + 1e-12)

        # 3. Deposit markets: bank supplies Dep_supply in nominal good units;
        # household A is in real composite units — multiply by P_CES to compare.
        dep_D_resid = (out["P_CES_D"] * out["A_D"] - bk["Dep_supply_D"]) / Kap_D_ss
        dep_F_resid = (out["P_CES_F"] * out["A_F"] - bk["Dep_supply_F"]) / Kap_F_ss

        # 4. Goods market D (pins p): Y_D = P_CES_D·C_D + I_D + NX_D + G_D
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
                     def_price_D=None, def_price_F=None,
                     def_real_D=None, def_real_F=None,
                     verbose=True, maxiter=300, y0=None,
                     init=None, risk_D=None, jac_cache=None,
                     hybr_factor=100.0, accept_tol=None):

    # SOLVE THE 2 COUNTRY TRANSITION PATH
    #
    # Strategy (solvers.py): damped Newton on an explicit finite-difference
    # Jacobian, which `jac_cache` (a caller-owned dict) carries ACROSS solves —
    # warm re-solves inside the CK/risk fixed points then cost a handful of
    # residual evaluations instead of hybr's unconditional 7T+1-call Jacobian
    # rebuild.  scipy hybr remains as the fallback for guesses far outside the
    # Newton basin, followed by a Newton polish (hybr alone stops on step size
    # and plateaus near max|resid| ≈ 5e-11).

    T = cal["T"]
    assert len(Z_D_path) == T and len(Z_F_path) == T

    # Empty arrays for initial guess if not provided--
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
                def_price_D=def_price_D, def_price_F=def_price_F,
                def_real_D=def_real_D, def_real_F=def_real_F,
                init=init, risk_D=risk_D)
    residual = make_residual(spec, verbose=verbose)
    # accept_tol override: default-branch probes accept 1e-9 (the tightest
    # test threshold) — quality-shocked branch systems can polish-stall a
    # few x above the 1e-10 default without any economic content.
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
        # Trust-region fallback for guesses outside the Newton basin, then a
        # Newton polish with a fresh Jacobian (replaces the old hybr polish
        # restarts and the cold-start krylov fallback, both of which stall:
        # hybr terminates on xtol, krylov diverges from cold starts here).
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
                         def_price_D=def_price_D, def_price_F=def_price_F,
                         def_real_D=def_real_D, def_real_F=def_real_F,
                         init=init, risk_D=risk_D)

    # μ monitor: the IC multiplier is a Lagrange multiplier — a negative value
    # on a SOLVED path means the always-binding IC is violated in equilibrium
    # (the imposed n_IC = n_ACCUM equality then manufactures a counterfactual
    # bank-recapitalization boom; see docs/sunspot_transition_study.md).
    # Checked here, NOT inside the residual, so Newton exploration is unaffected.
    mu_min_D = float(np.min(out["bk"]["mu_D"]))
    mu_min_F = float(np.min(out["bk"]["mu_F"]))
    if mu_min_D < 0 or mu_min_F < 0:
        print(f"  [transition] WARNING: IC multiplier negative on the solved "
              f"path (min mu_D={mu_min_D:+.4f}, min mu_F={mu_min_F:+.4f}) — "
              "always-binding IC violated; results in this region are suspect.")

    T_arr = cal["T"]
    zeros = np.zeros(T_arr)
    return dict(
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
        # Endogenous government debt (end-of-period = bank-held stock)
        b_gov_D=out["gov_D"]["b_gov_eop"], b_gov_F=out["gov_F"]["b_gov_eop"],
        b_gov_bop_D=out["gov_D"]["b_gov"], b_gov_bop_F=out["gov_F"]["b_gov"],
        coupon_D=out["gov_D"]["coupon"], coupon_F=out["gov_F"]["coupon"],
        net_issuance_D=out["gov_D"]["net_issuance"], net_issuance_F=out["gov_F"]["net_issuance"],
        # Default paths actually used (for plotting/tests)
        def_price_D=(def_price_D if def_price_D is not None else zeros),
        def_price_F=(def_price_F if def_price_F is not None else zeros),
        def_real_D=(def_real_D if def_real_D is not None else zeros),
        def_real_F=(def_real_F if def_real_F is not None else zeros),
        mu_min_D=mu_min_D, mu_min_F=mu_min_F,
        **out["bk"],   # all bank paths (alpha, mu, Omega, n, Q_b, rb, holdings, spreads)
        y_vec=y_sol,   # solved unknown vector (warm start for homotopy)
    )


def solve_transition_ck(ss, cal, Z_D_path, Z_F_path,
                        sunspot_D_path=None, sunspot_F_path=None,
                        def_real_D=None, def_real_F=None,
                        verbose=True, y0=None, jac_cache=None):
    """Cole-Kehoe sunspot solver (Bocola 2016 pass-through experiment).

    The sunspot xi_t ∈ [0,1] is the probability that lenders coordinate on
    the no-rollover equilibrium at t, conditional on the crisis zone being
    active (Bocola's exogenous default-risk process restricted to the CK
    zone).  It is PRICED into bonds:  def_price[t] = ck_default_prob(b[t]).
    In the risk-only baseline default is never REALIZED (def_real = 0):
    net worth falls purely from the mark-to-market revaluation of bonds.

    Because the debt stock is forward-integrated inside every residual
    evaluation (see _inner_economy), the only remaining fixed point is the
    crisis-zone indicator itself: def_price depends on b_gov/Y_ss crossing
    b_ck_low/b_ck_high.  When debt stays inside the crisis zone (the typical
    configuration: SS is in the zone), the map converges in one iteration.

    Returns the solve_transition dict plus 'sunspot_D', 'sunspot_F'.
    """
    T      = cal["T"]
    b_ss_D = cal["B_gov_D_ss"]
    b_ss_F = cal["B_gov_F_ss"]
    Y_ss_D = ss["ss_firm_D"]["Y_ss"]
    Y_ss_F = ss["ss_firm_F"]["Y_ss"]

    if sunspot_D_path is None:
        sunspot_D_path = np.zeros(T)
    if sunspot_F_path is None:
        sunspot_F_path = np.zeros(T)

    # Initial zone check at SS debt
    def_price_D = np.array([ck_default_prob(b_ss_D, Y_ss_D, cal, s, "D")
                            for s in sunspot_D_path])
    def_price_F = np.array([ck_default_prob(b_ss_F, Y_ss_F, cal, s, "F")
                            for s in sunspot_F_path])

    max_iter = cal.get("ck_max_iter", 25)
    tol      = cal.get("ck_tol", 1e-8)
    damp     = cal.get("ck_damping", 1.0)

    out = None
    for ck_it in range(max_iter):
        out = solve_transition(
            ss, cal, Z_D_path, Z_F_path,
            def_price_D=def_price_D, def_price_F=def_price_F,
            def_real_D=def_real_D, def_real_F=def_real_F,
            verbose=False, y0=y0, jac_cache=jac_cache,
        )
        y0 = out["y_vec"]   # warm start for the next zone iteration

        # Re-evaluate the crisis zone on the solved beginning-of-period debt
        def_price_D_new = np.array([
            ck_default_prob(out["b_gov_bop_D"][t], Y_ss_D, cal, sunspot_D_path[t], "D")
            for t in range(T)])
        def_price_F_new = np.array([
            ck_default_prob(out["b_gov_bop_F"][t], Y_ss_F, cal, sunspot_F_path[t], "F")
            for t in range(T)])

        err = max(np.max(np.abs(def_price_D_new - def_price_D)),
                  np.max(np.abs(def_price_F_new - def_price_F)))
        if verbose:
            print(f"  CK iter {ck_it + 1:2d}: zone err={err:.2e}"
                  f"  b_gov_D peak={np.max(out['b_gov_D']):.4f}"
                  f"  Q_bD[0]={out['Q_bD'][0]:.4f}")
        if err < tol:
            if verbose:
                print(f"  CK converged in {ck_it + 1} iteration(s).")
            break
        def_price_D = (1 - damp) * def_price_D + damp * def_price_D_new
        def_price_F = (1 - damp) * def_price_F + damp * def_price_F_new
    else:
        if verbose:
            print(f"  CK warning: zone indicator did not settle after {max_iter} iters "
                  f"(err={err:.2e}).")

    out["sunspot_D"] = sunspot_D_path
    out["sunspot_F"] = sunspot_F_path
    return out
