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

Household income (GHH, all in composite-good units):
  y_t(e) = (w_t/P_CES_t)·N_t·e + (Div_t − Tax_t)/P_CES_t
  where N_t is aggregate employment and e is the idiosyncratic productiviy.
  Dividing by P_CES_t converts nominal D-good flows to composite units consistent
  with the real wage w/P_CES and the EGM deposit stock (denominated in real units).
  The GHH labour-disutility vN_t = chi·N_t^(1+1/frisch)/(1+1/frisch) is
  subtracted from effective income in the EGM composite (passed via vN_path).
"""
import numpy as np
from scipy.optimize import root

from firms import solve_firm_path, markup_ss
from capital import solve_capital_path
from bank import solve_bank_paths
from household import solve_backward_transition
from distribution import forward_iterate, aggregate_assets, aggregate_consumption
from trade import ces_price, import_demand, trade_balance
from government import govt_transition, ck_default_prob, integrate_b_gov, bd_fundamental_default_prob


def _inner_economy(N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                   Z_D_path, Z_F_path, def_D_path, def_F_path, ss, cal,
                   b_gov_D_path=None, b_gov_F_path=None,
                   sunspot_D_path=None, sunspot_F_path=None):
    """Given the outer guesses, solve all inner blocks for both countries.

    Returns a comprehensive dict of all time paths and residual inputs.
    """
    # ── Firms ─────────────────────────────────────────────────────────────────
    firm_D = solve_firm_path(N_D, Kap_D, Z_D_path, cal, country="D")
    firm_F = solve_firm_path(N_F, Kap_F, Z_F_path, cal, country="F")

    # ── Capital ───────────────────────────────────────────────────────────────
    cap_D = solve_capital_path(Kap_D, ss["Kap_D_ss"], 1.0, firm_D["mpk"], cal, country="D")
    cap_F = solve_capital_path(Kap_F, ss["Kap_F_ss"], 1.0, firm_F["mpk"], cal, country="F")

    # ── Trade / CES ───────────────────────────────────────────────────────────
    P_CES_D = np.array([ces_price(p, cal, "D") for p in p_path])
    P_CES_F = np.array([ces_price(p, cal, "F") for p in p_path])

    T = len(p_path)
    if b_gov_D_path is None:
        b_gov_D_path = np.full(T, cal["B_gov_D_ss"])
    if b_gov_F_path is None:
        b_gov_F_path = np.full(T, cal["B_gov_F_ss"])

    # ── Bank paths (both banks simultaneously) ─────────────────────────────────
    # Bond market clearing uses the SS supply (fixed).  The endogenous b_gov
    # from the BD/CK outer loop feeds ONLY into govt_transition() (Bohn tax
    # rule and coupon accounting), NOT into the bank's IC clearing.  This
    # decouples the bank Newton problem from outer-loop iteration in b_gov,
    # ensuring the BD fixed-point map is contracting (Jacobian ≈ 0.85 < 1).
    bk = solve_bank_paths(
        Kap_D, Kap_F, cap_D["Q"], cap_F["Q"],
        cap_D["rk"], cap_F["rk"], rdep_D, rdep_F,
        p_path, cal["B_gov_D_ss"], cal["B_gov_F_ss"],
        cal, ss["ss_bank_D"], ss["ss_bank_F"],
        def_D_path=def_D_path, def_F_path=def_F_path,
        sunspot_D_path=sunspot_D_path, sunspot_F_path=sunspot_F_path,
    )

    # ── Government tax paths (Bohn rule on b_gov_path; CK survival in coupon) ──
    gov_D = govt_transition(cal, ss["gs_D"], bk["Q_bD"], def_D_path, b_gov_D_path, "D")
    gov_F = govt_transition(cal, ss["gs_F"], bk["Q_bF"], def_F_path, b_gov_F_path, "F",
                            p_path=p_path)

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
    # Non-labour income (Div_D - Tax_D) is in nominal D-goods; divide by P_CES to match
    # the real wage (w/P_CES) and the EGM deposit stock (real composite units).
    y_D_path = (w_real_D * N_D)[:, None] * e_D[None, :] + ((Div_D - gov_D["Tax"]) / P_CES_D)[:, None]
    y_F_path = (w_real_F * N_F)[:, None] * e_F[None, :] + ((Div_F - gov_F["Tax"]) / P_CES_F)[:, None]

    # ── Household EGM backward ────────────────────────────────────────────────
    r_D_path = np.concatenate(([cal["r_dep_D_target"]], rdep_D))
    r_F_path = np.concatenate(([cal["r_dep_F_target"]], rdep_F))

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

    D_D = ss["D_D_ss"]
    D_F = ss["D_F_ss"]

    for t in range(T):
        C_D_path[t] = aggregate_consumption(D_D, c_D_path[t])
        C_F_path[t] = aggregate_consumption(D_F, c_F_path[t])
        D_D = forward_iterate(D_D, a_pol_D_path[t], ss["a_grid_D"], ss["Pi_D"])
        D_F = forward_iterate(D_F, a_pol_F_path[t], ss["a_grid_F"], ss["Pi_F"])
        A_D_path[t] = aggregate_assets(D_D, ss["a_grid_D"])
        A_F_path[t] = aggregate_assets(D_F, ss["a_grid_F"])

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
    )


def solve_transition(ss, cal, Z_D_path, Z_F_path,
                     def_D_path=None, def_F_path=None,
                     b_gov_D_path=None, b_gov_F_path=None,
                     sunspot_D_path=None, sunspot_F_path=None,
                     verbose=True, maxiter=300):
    """Solve the two-country transition path.

    Arguments
    ---------
    ss            : steady-state dict from solve_steady_state()
    cal           : calibration dict (mutated in-place with chi, excess_return)
    Z_D_path      : TFP path for D (length T)
    Z_F_path      : TFP path for F (length T)
    def_D_path    : default rate path for D or None (zeros)
    def_F_path    : default rate path for F or None (zeros)
    b_gov_D_path  : beginning-of-period debt stock path for D, or None (B_gov_ss)
    b_gov_F_path  : beginning-of-period debt stock path for F, or None (B_gov_ss)

    Returns
    -------
    dict with all 7T outer unknowns and all inner time paths.
    """
    T = cal["T"]
    assert len(Z_D_path) == T and len(Z_F_path) == T

    if def_D_path is None:
        def_D_path = np.zeros(T)
    if def_F_path is None:
        def_F_path = np.zeros(T)

    Kap_D_ss = ss["Kap_D_ss"]
    Kap_F_ss = ss["Kap_F_ss"]
    N_ss     = 1.0

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

        try:
            out = _inner_economy(
                N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
                Z_D_path, Z_F_path, def_D_path, def_F_path, ss, cal,
                b_gov_D_path=b_gov_D_path, b_gov_F_path=b_gov_F_path,
                sunspot_D_path=sunspot_D_path, sunspot_F_path=sunspot_F_path,
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

        # 3. Deposit markets: bank supplies Dep_supply in nominal good units (D-goods/F-goods);
        # household A is in real composite units. Multiply by P_CES to convert to nominal before
        # comparing, so the residual is dimensionally homogeneous and Walras-consistent.
        dep_D_resid = (out["P_CES_D"] * out["A_D"] - bk["Dep_supply_D"]) / Kap_D_ss
        dep_F_resid = (out["P_CES_F"] * out["A_F"] - bk["Dep_supply_F"]) / Kap_F_ss

        # 4. Goods market D (pins p): Y_D = P_CES_D·C_D + I_D + NX_D + G_D
        # C_D is the composite-good index; expenditure in D-goods = P_CES_D·C_D.
        goods_D_resid = (firm_D["Y"] - out["P_CES_D"] * out["C_D"] - out["cap_D"]["I"]
                         - out["NX_D"] - G_D) / Y_ss_D

        resid = np.concatenate([
            cap_D_resid, cap_F_resid,
            lab_D_resid, lab_F_resid,
            dep_D_resid, dep_F_resid,
            goods_D_resid,
        ])

        if verbose:
            walras_D = np.max(np.abs(firm_D["Y"] - out["P_CES_D"] * out["C_D"] - out["cap_D"]["I"] - out["NX_D"] - G_D))
            walras_F = np.max(np.abs(firm_F["Y"] - out["P_CES_F"] * out["C_F"] - out["cap_F"]["I"] - out["NX_F"] - cal["G_F"]))
            print(f"  call {ncalls[0]:3d}: max|resid|={np.max(np.abs(resid)):.3e}"
                  f"  goods_F={walras_F:.3e}  goods_D={walras_D:.3e}")

        return resid

    accept_tol = max(cal["tol_mkt"] * 100, 1e-6)
    sol = root(residual, y0, method="hybr",
               options={"maxfev": max(maxiter * (7 * T + 1), 5000)})
    resid_norm = np.max(np.abs(residual(sol.x)))

    if resid_norm > accept_tol:
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
                         Z_D_path, Z_F_path, def_D_path, def_F_path, ss, cal,
                         b_gov_D_path=b_gov_D_path, b_gov_F_path=b_gov_F_path,
                         sunspot_D_path=sunspot_D_path, sunspot_F_path=sunspot_F_path)

    return dict(
        N_D=N_D, N_F=N_F, Kap_D=Kap_D, Kap_F=Kap_F,
        rdep_D=rdep_D, rdep_F=rdep_F, p=p_path,
        Z_D=Z_D_path, Z_F=Z_F_path,
        **{k + "_D": v for k, v in out["firm_D"].items()},
        **{k + "_F": v for k, v in out["firm_F"].items()},
        I_D=out["cap_D"]["I"],    Q_D=out["cap_D"]["Q"],    rk_D=out["cap_D"]["rk"],
        I_F=out["cap_F"]["I"],    Q_F=out["cap_F"]["Q"],    rk_F=out["cap_F"]["rk"],
        A_D=out["A_D"], C_D=out["C_D"], A_F=out["A_F"], C_F=out["C_F"],
        NX_D=out["NX_D"], NX_F=out["NX_F"],
        Div_D=out["Div_D"], Div_F=out["Div_F"],
        Tax_D=out["gov_D"]["Tax"], Tax_F=out["gov_F"]["Tax"],
        P_CES_D=out["P_CES_D"], P_CES_F=out["P_CES_F"],
        **out["bk"],   # all bank paths (alpha_D, mu_D, n_IC_D, n_D, ..., Q_bD, Q_bF, ...)
    )


def solve_transition_ck(ss, cal, Z_D_path, Z_F_path,
                        sunspot_D_path=None, sunspot_F_path=None, verbose=True):
    """Cole-Kehoe transition solver: outer fixed-point on (def_rate, b_gov).

    The CK circularity — bond price Q_bD needs def_rate needs b_gov needs Q_bD —
    is resolved by iterating the (def_rate, b_gov) pair around the existing 7T
    Newton solver until self-consistency is reached.

    Parameters
    ----------
    sunspot_D_path : (T,) AR(1) sunspot path ∈ [0,1] for country D.
                     0 = good equilibrium, 1 = full default probability.
    sunspot_F_path : (T,) same for F (usually zeros; F has no default).

    Returns
    -------
    Standard solve_transition dict plus 'b_gov_D', 'b_gov_F', 'def_D', 'def_F'.
    """
    T        = cal["T"]
    b_ss_D   = cal["B_gov_D_ss"]
    b_ss_F   = cal["B_gov_F_ss"]
    Y_ss_D   = ss["ss_firm_D"]["Y_ss"]
    Y_ss_F   = ss["ss_firm_F"]["Y_ss"]

    if sunspot_D_path is None:
        sunspot_D_path = np.zeros(T)
    if sunspot_F_path is None:
        sunspot_F_path = np.zeros(T)

    # Initial def_rate: evaluated at SS debt level
    def_D = np.array([ck_default_prob(b_ss_D, Y_ss_D, cal, s, "D") for s in sunspot_D_path])
    def_F = np.array([ck_default_prob(b_ss_F, Y_ss_F, cal, s, "F") for s in sunspot_F_path])
    b_gov_D = np.full(T, b_ss_D)
    b_gov_F = np.full(T, b_ss_F)

    max_iter = cal.get("ck_max_iter", 25)
    tol      = cal.get("ck_tol", 1e-5)
    damp     = cal.get("ck_damping", 0.5)

    for ck_it in range(max_iter):
        out = solve_transition(
            ss, cal, Z_D_path, Z_F_path,
            def_D_path=def_D, def_F_path=def_F,
            b_gov_D_path=b_gov_D, b_gov_F_path=b_gov_F,
            verbose=False,
        )

        # Forward-integrate new debt stock from budget identity.
        # integrate_b_gov returns END-of-period stocks [b1, b2, ..., bT].
        b_gov_D_eop = integrate_b_gov(
            b_ss_D, out["Q_bD"], def_D, out["Tax_D"], cal, "D")
        b_gov_F_eop = integrate_b_gov(
            b_ss_F, out["Q_bF"], def_F, out["Tax_F"], cal, "F")

        # Convert to BEGINNING-of-period stocks (lagged by one) for:
        # (a) Bohn rule in govt_transition  (b) CK threshold check.
        # b_gov_bop[t] = stock at start of period t = end of period t-1.
        b_gov_D_new = np.concatenate([[b_ss_D], b_gov_D_eop[:-1]])
        b_gov_F_new = np.concatenate([[b_ss_F], b_gov_F_eop[:-1]])

        # New default probabilities from Cole-Kehoe crisis-zone function
        def_D_new = np.array([
            ck_default_prob(b_gov_D_new[t], Y_ss_D, cal, sunspot_D_path[t], "D")
            for t in range(T)])
        def_F_new = np.array([
            ck_default_prob(b_gov_F_new[t], Y_ss_F, cal, sunspot_F_path[t], "F")
            for t in range(T)])

        err_def   = max(np.max(np.abs(def_D_new - def_D)),
                        np.max(np.abs(def_F_new - def_F)))
        err_bgov  = max(np.max(np.abs(b_gov_D_new - b_gov_D)) / b_ss_D,
                        np.max(np.abs(b_gov_F_new - b_gov_F)) / b_ss_F)
        err = max(err_def, err_bgov)

        if verbose:
            print(f"  CK iter {ck_it + 1:2d}: err={err:.2e}"
                  f"  (def={err_def:.2e}  b_gov={err_bgov:.2e})")

        # Damped update
        def_D   = (1.0 - damp) * def_D   + damp * def_D_new
        def_F   = (1.0 - damp) * def_F   + damp * def_F_new
        b_gov_D = (1.0 - damp) * b_gov_D + damp * b_gov_D_new
        b_gov_F = (1.0 - damp) * b_gov_F + damp * b_gov_F_new

        if err < tol:
            if verbose:
                print(f"  CK converged in {ck_it + 1} iterations.")
            break
    else:
        if verbose:
            print(f"  CK warning: did not converge after {max_iter} iterations (err={err:.2e}).")

    # Report end-of-period debt stocks for diagnostics/plotting
    out["b_gov_D"] = b_gov_D_eop
    out["b_gov_F"] = b_gov_F_eop
    out["def_D"]   = def_D
    out["def_F"]   = def_F
    return out


def _anderson_step(G_hist, F_hist, m=3, reg=1e-10):
    """Anderson(m) mixing step for fixed-point iteration g(x) = x.

    G_hist : list of g(x_k) values (most recent last).
    F_hist : list of f(x_k) = g(x_k) - x_k residuals (most recent last).

    Returns the Anderson-accelerated next iterate x_{k+1}.
    Solves: min_{c, 1'c=1} ||F @ c||_2  →  x_next = G @ c.
    """
    mk = min(len(G_hist), m)
    G = np.column_stack(G_hist[-mk:])   # (n, mk)
    F = np.column_stack(F_hist[-mk:])   # (n, mk)

    if mk == 1:
        return G[:, 0].copy()

    FtF = F.T @ F + reg * np.eye(mk)
    ones = np.ones(mk)
    try:
        v = np.linalg.solve(FtF, ones)
        c = v / (ones @ v)
    except np.linalg.LinAlgError:
        c = np.ones(mk) / mk
    return G @ c


def solve_transition_bd(ss, cal, Z_D_path, Z_F_path,
                        sunspot_D_path=None, sunspot_F_path=None, verbose=True):
    """Bocola-Dovis (2019) transition solver.

    Mechanism:
      sunspot xi_{t+1} → GK IC tightens (lbD_eff = lbD + psi_bd * xi)
        → Q_bD falls (bond price compressed)
        → government issues MORE face-value bonds to raise the same revenue
        → b_gov rises (beliefs → fundamentals)
        → if b_gov/Y_ss >= b_ck_high: fundamental default fires (surv < 1)
        → bank NW falls → investment falls → amplification

    Unlike CK, the sunspot does NOT directly cause default — it only lowers
    the bond price through the IC denominator.  Fundamental default only
    triggers if debt crosses b_ck_high (the Bocola-Dovis B̄ boundary).

    The outer loop uses Anderson(m) acceleration to converge the
    (b_gov, def_fund) fixed-point while the sunspot path is held fixed.

    Parameters
    ----------
    sunspot_D_path : (T,) path xi_t ∈ [0,1] — exogenous rollover-risk shock.
    sunspot_F_path : (T,) same for F (default zeros; F has no default).

    Returns
    -------
    Standard solve_transition dict plus:
      'b_gov_D', 'b_gov_F' : end-of-period debt paths
      'def_D',  'def_F'    : fundamental default rate paths (0/1)
      'sunspot_D', 'sunspot_F' : the exogenous sunspot paths
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

    # Initial guess: no fundamental default; debt at SS level
    def_fund_D = np.zeros(T)
    def_fund_F = np.zeros(T)
    b_gov_D    = np.full(T, b_ss_D)
    b_gov_F    = np.full(T, b_ss_F)

    max_iter = cal.get("bd_max_iter", 50)
    tol      = cal.get("bd_tol", 1e-4)
    m_aa     = cal.get("bd_anderson_m", 3)

    G_hist_D, F_hist_D = [], []   # Anderson history for b_gov_D
    G_hist_F, F_hist_F = [], []   # Anderson history for b_gov_F

    for bd_it in range(max_iter):
        out = solve_transition(
            ss, cal, Z_D_path, Z_F_path,
            def_D_path=def_fund_D, def_F_path=def_fund_F,
            b_gov_D_path=b_gov_D, b_gov_F_path=b_gov_F,
            sunspot_D_path=sunspot_D_path, sunspot_F_path=sunspot_F_path,
            verbose=False,
        )

        # Forward-integrate new debt stocks from budget identity
        b_gov_D_eop = integrate_b_gov(
            b_ss_D, out["Q_bD"], def_fund_D, out["Tax_D"], cal, "D")
        b_gov_F_eop = integrate_b_gov(
            b_ss_F, out["Q_bF"], def_fund_F, out["Tax_F"], cal, "F")

        # Beginning-of-period stocks (lagged one period) for Bohn rule and threshold check
        b_gov_D_new = np.concatenate([[b_ss_D], b_gov_D_eop[:-1]])
        b_gov_F_new = np.concatenate([[b_ss_F], b_gov_F_eop[:-1]])

        # Fundamental default: ONLY when b/Y >= b_ck_high (BD threshold B̄)
        def_fund_D_new = np.array([
            bd_fundamental_default_prob(b_gov_D_new[t], Y_ss_D, cal, "D")
            for t in range(T)])
        def_fund_F_new = np.array([
            bd_fundamental_default_prob(b_gov_F_new[t], Y_ss_F, cal, "F")
            for t in range(T)])

        err_def  = max(np.max(np.abs(def_fund_D_new - def_fund_D)),
                       np.max(np.abs(def_fund_F_new - def_fund_F)))
        err_bgov = max(np.max(np.abs(b_gov_D_new - b_gov_D)) / b_ss_D,
                       np.max(np.abs(b_gov_F_new - b_gov_F)) / b_ss_F)
        err = max(err_def, err_bgov)

        if verbose:
            bd_spread_peak = cal.get("psi_bd_D", 0.0) * np.max(sunspot_D_path)
            print(f"  BD iter {bd_it + 1:2d}: err={err:.2e}"
                  f"  b_gov_D_peak={np.max(b_gov_D_new):.4f}"
                  f"  bd_spread_peak={bd_spread_peak:.4f}"
                  f"  (def={err_def:.2e}  bgov={err_bgov:.2e})")

        if err < tol:
            if verbose:
                print(f"  BD converged in {bd_it + 1} iterations.")
            break

        # Anderson acceleration for b_gov paths (binary def_fund accepted directly)
        G_hist_D.append(b_gov_D_new.copy())
        F_hist_D.append((b_gov_D_new - b_gov_D).copy())
        G_hist_F.append(b_gov_F_new.copy())
        F_hist_F.append((b_gov_F_new - b_gov_F).copy())

        b_gov_D = _anderson_step(G_hist_D, F_hist_D, m=m_aa)
        b_gov_F = _anderson_step(G_hist_F, F_hist_F, m=m_aa)
        # Safety clip: debt cannot go below 50% of SS (numerical guard)
        b_gov_D = np.maximum(b_gov_D, 0.5 * b_ss_D)
        b_gov_F = np.maximum(b_gov_F, 0.5 * b_ss_F)

        def_fund_D = def_fund_D_new.copy()
        def_fund_F = def_fund_F_new.copy()
    else:
        if verbose:
            print(f"  BD warning: did not converge after {max_iter} iters (err={err:.2e}).")

    out["b_gov_D"]   = b_gov_D_eop
    out["b_gov_F"]   = b_gov_F_eop
    out["def_D"]     = def_fund_D
    out["def_F"]     = def_fund_F
    out["sunspot_D"] = sunspot_D_path
    out["sunspot_F"] = sunspot_F_path
    return out
