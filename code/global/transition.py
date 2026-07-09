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

from firms import solve_firm_path, markup_ss
from capital import solve_capital_path
from bank import bank_backward, bank_forward
from household import solve_backward_transition
from distribution import forward_iterate, aggregate_assets, aggregate_consumption
from trade import ces_price, import_demand, trade_balance
from government import govt_transition, ck_default_prob


def _inner_economy(N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                   Z_D_path, Z_F_path, ss, cal,
                   def_price_D=None, def_price_F=None,
                   def_real_D=None, def_real_F=None,
                   init=None, risk_D=None):
    """Given the outer guesses, solve all inner blocks for both countries.

    init   : optional dict of period-(-1)/period-0 initial conditions for
             paths that start mid-crisis (default branches, policy runs):
             D_D, D_F (household distributions), bank_D, bank_F (see
             bank_forward), b_gov0_D/F, Kap_lag_D/F, Q_lag_D/F,
             Q_bD_lag, Q_bF_lag, p_lag, P_lag_D/F.  Missing keys → SS.
    risk_D : optional Bocola risk-channel inputs for bank_backward
             (two-branch expectations over the D-default event).

    Returns a comprehensive dict of all time paths and residual inputs.
    """
    T = len(p_path)
    if init is None:
        init = {}

    # ── Firms ─────────────────────────────────────────────────────────────────
    firm_D = solve_firm_path(N_D, Kap_D, Z_D_path, cal, country="D")
    firm_F = solve_firm_path(N_F, Kap_F, Z_F_path, cal, country="F")

    # ── Capital ───────────────────────────────────────────────────────────────
    cap_D = solve_capital_path(Kap_D, init.get("Kap_lag_D", ss["Kap_D_ss"]),
                               init.get("Q_lag_D", 1.0), firm_D["mpk"], cal, country="D")
    cap_F = solve_capital_path(Kap_F, init.get("Kap_lag_F", ss["Kap_F_ss"]),
                               init.get("Q_lag_F", 1.0), firm_F["mpk"], cal, country="F")

    # ── Trade / CES ───────────────────────────────────────────────────────────
    P_CES_D = np.array([ces_price(p, cal, "D") for p in p_path])
    P_CES_F = np.array([ces_price(p, cal, "F") for p in p_path])

    # ── Bank backward pass: prices and cross-border FOC holdings ─────────────
    bwd = bank_backward(
        cap_D["rk"], cap_F["rk"], rdep_D, rdep_F, p_path,
        cal, ss["ss_bank_D"], ss["ss_bank_F"],
        def_price_D=def_price_D, def_price_F=def_price_F, risk_D=risk_D,
    )

    # ── Government: forward-integrate debt with Bohn tax (realized flows) ────
    gov_D = govt_transition(cal, ss["gs_D"], bwd["Q_bD"], def_real_D, "D",
                            b_gov0=init.get("b_gov0_D"),
                            b_anchor=init.get("b_anchor_D"))
    gov_F = govt_transition(cal, ss["gs_F"], bwd["Q_bF"], def_real_F, "F",
                            b_gov0=init.get("b_gov0_F"),
                            b_anchor=init.get("b_anchor_F"))

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
    )
    bk = {**bwd, **fwd}

    # ── Dividends ─────────────────────────────────────────────────────────────
    mc_D = markup_ss(cal, "D")
    mc_F = markup_ss(cal, "F")
    Div_D = (1 - mc_D) * firm_D["Y"] + cap_D["cap_profit"] + bk["div_D"]
    Div_F = (1 - mc_F) * firm_F["Y"] + cap_F["cap_profit"] + bk["div_F"]

    # ── GHH income and labour disutility ───────────────────────────────────────
    chi_D   = cal["chi_D"];   frisch_D = cal["frisch_D"]
    chi_F   = cal["chi_F"];   frisch_F = cal["frisch_F"]
    sigma_D = cal["sigma_D"]; sigma_F  = cal["sigma_F"]

    # Labour disutility v(N) = chi·N^(1+1/frisch)/(1+1/frisch)
    vN_D = chi_D * N_D ** (1 + 1 / frisch_D) / (1 + 1 / frisch_D)
    vN_F = chi_F * N_F ** (1 + 1 / frisch_F) / (1 + 1 / frisch_F)

    # Effective real income (real wage × idiosyncratic e + non-labour income)
    w_real_D = firm_D["w"] / P_CES_D     # real wage in D-consumption units
    w_real_F = firm_F["w"] / P_CES_F
    e_D = ss["e_D"];  e_F = ss["e_F"]

    # y_path shape (T, n_e): individual income in composite-good units.
    y_D_path = (w_real_D * N_D)[:, None] * e_D[None, :] + ((Div_D - gov_D["Tax"]) / P_CES_D)[:, None]
    y_F_path = (w_real_F * N_F)[:, None] * e_F[None, :] + ((Div_F - gov_F["Tax"]) / P_CES_F)[:, None]

    # ── Household EGM backward ────────────────────────────────────────────────
    # Fisher-equation real returns: r_real[t] = (1+rdep[t-1]) * P_CES[t-1]/P_CES[t] - 1
    # Period -1 anchors (rate locked pre-path, price level pre-path) come from
    # init when the path starts mid-crisis; SS otherwise.
    rdep_prev_D = (init["bank_D"]["rdep_prev"] if "bank_D" in init
                   else cal["r_dep_D_target"])
    rdep_prev_F = (init["bank_F"]["rdep_prev"] if "bank_F" in init
                   else cal["r_dep_F_target"])
    P_CES_D_ext = np.concatenate([[init.get("P_lag_D", 1.0)], P_CES_D, [1.0]])
    P_CES_F_ext = np.concatenate([[init.get("P_lag_F", 1.0)], P_CES_F, [1.0]])
    rdep_D_full = np.concatenate([[rdep_prev_D], rdep_D])  # length T+1
    rdep_F_full = np.concatenate([[rdep_prev_F], rdep_F])
    r_D_path = (1.0 + rdep_D_full) * P_CES_D_ext[:-1] / P_CES_D_ext[1:] - 1.0
    r_F_path = (1.0 + rdep_F_full) * P_CES_F_ext[:-1] / P_CES_F_ext[1:] - 1.0

    c_D_path, a_pol_D_path = solve_backward_transition(
        ss["a_grid_D"], ss["Pi_D"], r_D_path, y_D_path, ss["c_D_ss"],
        ss["beta_D_ss"], sigma_D, cal["a_min_D"], vN_path=vN_D,
    )
    c_F_path, a_pol_F_path = solve_backward_transition(
        ss["a_grid_F"], ss["Pi_F"], r_F_path, y_F_path, ss["c_F_ss"],
        ss["beta_F_ss"], sigma_F, cal["a_min_F"], vN_path=vN_F,
    )

    # ── Distribution forward ──────────────────────────────────────────────────
    A_D_path = np.empty(T)
    C_D_path = np.empty(T)
    A_F_path = np.empty(T)
    C_F_path = np.empty(T)

    D_D = init.get("D_D", ss["D_D_ss"])
    D_F = init.get("D_F", ss["D_F_ss"])

    # Start-of-period distributions (D_start[t] enters period t); D_start[0]
    # is the initial condition.  Needed to launch default branches mid-path.
    D_start_D = np.empty((T + 1,) + D_D.shape)
    D_start_F = np.empty((T + 1,) + D_F.shape)

    for t in range(T):
        D_start_D[t] = D_D
        D_start_F[t] = D_F
        C_D_path[t] = aggregate_consumption(D_D, c_D_path[t])
        C_F_path[t] = aggregate_consumption(D_F, c_F_path[t])
        D_D = forward_iterate(D_D, a_pol_D_path[t], ss["a_grid_D"], ss["Pi_D"])
        D_F = forward_iterate(D_F, a_pol_F_path[t], ss["a_grid_F"], ss["Pi_F"])
        A_D_path[t] = aggregate_assets(D_D, ss["a_grid_D"])
        A_F_path[t] = aggregate_assets(D_F, ss["a_grid_F"])
    D_start_D[T] = D_D
    D_start_F[T] = D_F

    # ── Trade ─────────────────────────────────────────────────────────────────
    IM_D = np.array([import_demand(p_path[t], C_D_path[t], P_CES_D[t], cal, "D")
                     for t in range(T)])
    IM_F = np.array([import_demand(p_path[t], C_F_path[t], P_CES_F[t], cal, "F")
                     for t in range(T)])
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
    )


def solve_transition(ss, cal, Z_D_path, Z_F_path,
                     def_price_D=None, def_price_F=None,
                     def_real_D=None, def_real_F=None,
                     verbose=True, maxiter=300, y0=None,
                     init=None, risk_D=None, try_krylov=True,
                     hybr_factor=100.0):
    """Solve the two-country transition path.

    Arguments
    ---------
    ss            : steady-state dict from solve_steady_state()
    cal           : calibration dict (mutated in-place with chi, excess_return)
    Z_D_path      : TFP path for D (length T)
    Z_F_path      : TFP path for F (length T)
    def_price_D/F : PRICED default probability paths or None (zeros)
    def_real_D/F  : REALIZED default (haircut) paths or None (zeros)
    y0            : optional warm-start for the 7T unknown vector (homotopy)
    init          : optional initial-conditions dict for mid-crisis starts
                    (see _inner_economy docstring); None → SS start
    risk_D        : optional Bocola risk-channel inputs (see bank_backward /
                    risk_branch.py); None → risk-neutral pricing
    hybr_factor   : hybr initial trust-region factor (scipy default 100).
                    Large trial steps can cross the alpha-explosion penalty
                    wall in bank_backward and poison the finite-difference
                    Jacobian, stalling hybr at the warm start; a small
                    factor (e.g. 0.1) keeps early steps feasible — used by
                    risk_branch's ladder rescue for flat (near-SS) starts.

    Returns
    -------
    dict with all 7T outer unknowns and all inner time paths, including the
    endogenous government debt paths b_gov_D/b_gov_F (end-of-period).
    """
    T = cal["T"]
    assert len(Z_D_path) == T and len(Z_F_path) == T

    Kap_D_ss = ss["Kap_D_ss"]
    Kap_F_ss = ss["Kap_F_ss"]
    N_ss     = 1.0

    if y0 is None:
        y0 = np.concatenate([
            np.full(T, N_ss),             # N_D
            np.full(T, N_ss),             # N_F
            np.full(T, Kap_D_ss),         # Kap_D
            np.full(T, Kap_F_ss),         # Kap_F
            np.full(T, cal["r_dep_D_target"]),  # rdep_D
            np.full(T, cal["r_dep_F_target"]),  # rdep_F
            np.full(T, ss["p_ss"]),       # p
        ])

    ncalls = [0]
    chi_D   = cal["chi_D"];   frisch_D = cal["frisch_D"]
    chi_F   = cal["chi_F"];   frisch_F = cal["frisch_F"]
    n_ss_D  = ss["ss_bank_D"]["n_ss"]
    n_ss_F  = ss["ss_bank_F"]["n_ss"]
    Y_ss_D  = ss["ss_firm_D"]["Y_ss"]
    G_D     = cal["G_D"]

    def residual(y):
        ncalls[0] += 1
        N_D, N_F     = y[:T], y[T:2*T]
        Kap_D, Kap_F = y[2*T:3*T], y[3*T:4*T]
        rdep_D, rdep_F = y[4*T:5*T], y[5*T:6*T]
        p_path       = y[6*T:7*T]

        # Domain guard: fractional powers of non-positive p/N/K produce NaNs
        # that do NOT raise — penalize explicitly before they poison hybr.
        if (np.any(p_path <= 0.05) or np.any(N_D <= 0.01) or np.any(N_F <= 0.01)
                or np.any(Kap_D <= 0.1) or np.any(Kap_F <= 0.1)):
            return np.full(7 * T, 10.0)

        try:
            out = _inner_economy(
                N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                Z_D_path, Z_F_path, ss, cal,
                def_price_D=def_price_D, def_price_F=def_price_F,
                def_real_D=def_real_D, def_real_F=def_real_F,
                init=init, risk_D=risk_D,
            )
        except (ValueError, RuntimeError, FloatingPointError) as e:
            if verbose:
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

    accept_tol = max(cal["tol_mkt"] * 100, 1e-6)
    sol = root(residual, y0, method="hybr",
               options={"maxfev": max(maxiter * (7 * T + 1), 5000),
                        "factor": hybr_factor})
    resid_norm = np.max(np.abs(residual(sol.x)))

    if resid_norm > accept_tol and try_krylov:
        sol2 = root(residual, y0, method="krylov",
                    options={"fatol": cal["tol_mkt"], "maxiter": maxiter})
        resid_norm2 = np.max(np.abs(residual(sol2.x)))
        if resid_norm2 < resid_norm:
            sol, resid_norm = sol2, resid_norm2

    if resid_norm > accept_tol:
        raise RuntimeError(f"Transition path did not converge: max|resid|={resid_norm:.3e}")

    N_D, N_F     = sol.x[:T], sol.x[T:2*T]
    Kap_D, Kap_F = sol.x[2*T:3*T], sol.x[3*T:4*T]
    rdep_D, rdep_F = sol.x[4*T:5*T], sol.x[5*T:6*T]
    p_path       = sol.x[6*T:7*T]

    out = _inner_economy(N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                         Z_D_path, Z_F_path, ss, cal,
                         def_price_D=def_price_D, def_price_F=def_price_F,
                         def_real_D=def_real_D, def_real_F=def_real_F,
                         init=init, risk_D=risk_D)

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
        **out["bk"],   # all bank paths (alpha, mu, Omega, n, Q_b, rb, holdings, spreads)
        y_vec=sol.x,   # solved unknown vector (warm start for homotopy)
    )


def solve_transition_ck(ss, cal, Z_D_path, Z_F_path,
                        sunspot_D_path=None, sunspot_F_path=None,
                        def_real_D=None, def_real_F=None,
                        verbose=True, y0=None):
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
            verbose=False, y0=y0,
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
