# BOCOLA (2016) RISK CHANNEL VIA REPRESENTATIVE POST-EVENT BRANCHES.
# Bankers weight a no-event continuation (the base path) against a default
# branch of probability pi_t and, when TPI is priced, a "backstop reneged"
# branch. Pairing HIGH marginal valuations Omega^d with LOW branch payoffs puts
# a covariance premium on bonds and capital (precautionary deleveraging).
# pi_t is EXOGENOUS (Bocola eqs. 11-12): an input path, never a function of
# debt. One representative branch (event at tau*=1) is reused at every date, and
# pi = 0 nests the risk-neutral model exactly. The feared default event is a
# PURE HAIRCUT to recovery_rate_D on the whole claim, so the default-state
# recession has to arise endogenously through bank balance sheets.
# Approximations: Lambda^nd = beta_inter on the base path; Lambda^d from a
# household-SDF proxy on aggregate income; the household deposit Euler is
# pi-blind (faithful to Bocola, where deposits are riskless too).
import time

import numpy as np

from transition import solve_transition

# Recap shares used ONLY as a numerical continuation ladder when the direct
# branch solve stalls: larger injections cushion the net-worth wipeout, which
# shrinks the jump from the shifted-base warm start. Intermediate solves are
# auxiliary and discarded — the branch that is returned always has zero recap.
_RECAP_LADDER = (0.5, 0.35, 0.2, 0.1, 0.05)


def extract_init_state(out, ss, cal, tau):
    # INITIAL-CONDITIONS DICT FOR A TRANSITION STARTING AT BASE PERIOD tau.
    s = tau - 1   # everything the branch inherits is a tau-1 object
    n_D_s = out["n_D"][s];  n_F_s = out["n_F"][s]
    p_s = out["p"][s]
    return dict(
        D_D=out["D_start_D"][tau], D_F=out["D_start_F"][tau],
        # portfolio shares on ACTUAL net worth (n_IC coincides only when the IC binds)
        bank_D=dict(
            n_prev=n_D_s,
            kappa_prev=out["Q_D"][s] * out["Kap_D"][s] / n_D_s,
            phi_bdom_prev=out["Q_bD"][s] * out["b_D_D"][s] / n_D_s,
            phi_bfor_prev=p_s * out["Q_bF"][s] * out["b_F_D"][s] / n_D_s,
            rdep_prev=out["rdep_D"][s],
        ),
        bank_F=dict(
            n_prev=n_F_s,
            kappa_prev=out["Q_F"][s] * out["Kap_F"][s] / n_F_s,
            phi_bdom_prev=out["Q_bF"][s] * out["b_F_F"][s] / n_F_s,
            phi_bfor_prev=out["Q_bD"][s] * out["b_D_F"][s] / (p_s * n_F_s),
            rdep_prev=out["rdep_F"][s],
        ),
        b_gov0_D=out["b_gov_D"][s], b_gov0_F=out["b_gov_F"][s],
        Kap_lag_D=out["Kap_D"][s], Kap_lag_F=out["Kap_F"][s],
        Q_lag_D=out["Q_D"][s], Q_lag_F=out["Q_F"][s],
        Q_bD_lag=out["Q_bD"][s], Q_bF_lag=out["Q_bF"][s],
        p_lag=p_s, P_lag_D=out["P_CES_D"][s], P_lag_F=out["P_CES_F"][s],
        # the CB's carried-over holding entering tau; without it a branch launched
        # from a state with an active CB position drops its legacy coupon income
        cb_buy_D_lag0=out["cb_buy_D"][s],
    )


def _shift_path(x, tau, T):
    # SHIFT A PATH FORWARD BY tau, PADDED WITH ITS FINAL VALUE.
    idx = np.minimum(np.arange(T) + tau, T - 1)
    return np.asarray(x)[idx]


def _shifted_y0(out, tau, T):
    # BRANCH NEWTON WARM START: THE BASE SOLUTION SHIFTED BY tau.
    blocks = [out["N_D"], out["N_F"], out["Kap_D"], out["Kap_F"],
              out["rdep_D"], out["rdep_F"], out["p"]]
    return np.concatenate([_shift_path(b, tau, T) for b in blocks])


def solve_default_branch(out, ss, cal, tau=1, verbose=False, y0=None,
                         jac_cache=None):
    # SOLVE THE ONE FEARED EVENT: A PURE HAIRCUT ON THE OUTSTANDING STOCK.
    # def_real_D[0] = 1 at recovery recovery_rate_D, with the Bohn rule
    # re-anchored to the surviving stock (taxing the pre-haircut stock would be
    # a ~31%-of-GDP artifact). Infeasible even after the continuation -> RAISES.
    T = cal["T"]
    init = extract_init_state(out, ss, cal, tau)

    rec_D = cal["recovery_rate_D"]
    Z_D_b = _shift_path(out["Z_D"], tau, T)
    Z_F_b = _shift_path(out["Z_F"], tau, T)

    def_real_D = np.zeros(T)
    def_real_D[0] = 1.0
    init["b_anchor_D"] = init["b_gov0_D"] * rec_D

    bond_exposure = out["Q_bD"][tau - 1] * out["b_D_D"][tau - 1]

    def _attempt(share, y_s):
        # ONE BRANCH SOLVE AT A GIVEN RECAP SHARE; None IF THE SOLVER STALLS.
        if share > 0.0:
            recap = np.zeros(T)
            recap[0] = share * (1.0 - rec_D) * bond_exposure
            init["recap_D_path"] = recap
        else:
            init.pop("recap_D_path", None)
        try:
            return solve_transition(
                ss, cal, Z_D_b, Z_F_b,
                def_real_D=def_real_D,
                verbose=False, y0=y_s, init=init,
                jac_cache=jac_cache, accept_tol=1e-9,
            )
        except (RuntimeError, ValueError):
            return None

    y_start = y0 if y0 is not None else _shifted_y0(out, tau, T)
    branch = _attempt(0.0, y_start)

    if branch is None:
        if verbose:
            print("  [risk_branch] direct solve stalled; recap-share continuation")
        y_c = y_start
        for share in _RECAP_LADDER:
            aux = _attempt(share, y_c)
            if aux is not None:
                y_c = aux["y_vec"]
        branch = _attempt(0.0, y_c)

    if branch is None:
        raise RuntimeError(
            f"default branch (pure haircut, recovery={rec_D}) infeasible even "
            "after the recap continuation.")

    branch["recap_D_path"] = np.zeros(T)
    if verbose:
        print(f"  [risk_branch] branch solved: n_D(0)/n_ss = "
              f"{branch['n_D'][0] / ss['ss_bank_D']['n_ss']:.3f}, "
              f"Q_bD(0) = {branch['Q_bD'][0]:.4f}, rk_D(0) = {branch['rk_D'][0]:+.4f}")
    return branch


def solve_tpi_branch(out, ss, cal, tau=1, verbose=False, y0=None, jac_cache=None):
    # SOLVE THE ONE FEARED TPI EVENT: THE CB BACKSTOP RENEGES AT tau.
    # Not a haircut — def_real_D stays 0 and the branch is a plain transition
    # (no further priced TPI or mechanical purchases inside the branch itself,
    # exactly as the default branch carries no further priced default risk).
    T = cal["T"]
    try:
        branch = solve_transition(
            ss, cal, _shift_path(out["Z_D"], tau, T), _shift_path(out["Z_F"], tau, T),
            verbose=False, y0=(y0 if y0 is not None else _shifted_y0(out, tau, T)),
            init=extract_init_state(out, ss, cal, tau),
            jac_cache=jac_cache, accept_tol=1e-9,
        )
    except (RuntimeError, ValueError):
        raise RuntimeError("TPI-reneged branch (backstop withdrawn) infeasible.")
    if verbose:
        print(f"  [risk_branch] TPI branch solved: Q_bD(0) = {branch['Q_bD'][0]:.4f}, "
              f"rk_D(0) = {branch['rk_D'][0]:+.4f}")
    return branch


def _income_sdf(branch, base, cal, T, label):
    # HOUSEHOLD-SDF PROXY Lambda^branch ON AGGREGATE INCOME, WITH A SIGN GATE.
    # Marginal value is high where branch output is low relative to the base
    # continuation at t+1. If the branch is NOT a recession the proxy would
    # produce a negative premium, so it falls back to beta with a warning.
    Y_b_D = branch["Y_D"][0]
    Y_b_F = branch["Y_F"][0]
    ratio_D = Y_b_D / _shift_path(base["Y_D"], 1, T)
    if np.any(ratio_D >= 1.0):
        print(f"  [risk_branch] WARNING: {label} income-SDF sign gate tripped — "
              f"branch Y_D(0) = {Y_b_D:.4f} is not below the base path. "
              "Lambda falls back to beta.")
        return np.full(T, cal["beta_inter_D"]), np.full(T, cal["beta_inter_F"])
    Lam_D = cal["beta_inter_D"] * ratio_D ** (-cal["sigma_D"])
    Lam_F = cal["beta_inter_F"] * (Y_b_F / _shift_path(base["Y_F"], 1, T)) ** (-cal["sigma_F"])
    return Lam_D, Lam_F


def make_risk_inputs(branch, base, ss, cal):
    # DEFAULT-BRANCH OBJECTS -> risk_D DICT FOR bank_backward.
    T = cal["T"]
    Lam_d_D, Lam_d_F = _income_sdf(branch, base, cal, T, "default")
    # Omega^d = Lambda^d * [f + (1-f)*alpha^d(0)] (Bocola kernel weights)
    Omega_d_D = Lam_d_D * (cal["f_D"] + (1 - cal["f_D"]) * branch["alpha_D"][0])
    Omega_d_F = Lam_d_F * (cal["f_F"] + (1 - cal["f_F"]) * branch["alpha_F"][0])

    return dict(
        Omega_d_D=Omega_d_D, Omega_d_F=Omega_d_F,
        rk_d_D=branch["rk_D"][0], rk_d_F=branch["rk_F"][0],
        Q_bD_d=branch["Q_bD"][0], Q_bF_d=branch["Q_bF"][0],
        p_d=branch["p"][0],
        surv_d=cal["recovery_rate_D"],   # haircut on the whole claim
    )


def make_tpi_inputs(branch, base, ss, cal):
    # TPI-RENEGED-BRANCH OBJECTS -> tpi_D DICT FOR bank_backward.
    T = cal["T"]
    Lam_tpi_D, Lam_tpi_F = _income_sdf(branch, base, cal, T, "TPI")
    Omega_tpi_D = Lam_tpi_D * (cal["f_D"] + (1 - cal["f_D"]) * branch["alpha_D"][0])
    Omega_tpi_F = Lam_tpi_F * (cal["f_F"] + (1 - cal["f_F"]) * branch["alpha_F"][0])

    return dict(
        Omega_tpi_D=Omega_tpi_D, Omega_tpi_F=Omega_tpi_F,
        rk_tpi_d_D=branch["rk_D"][0], rk_tpi_d_F=branch["rk_F"][0],
        Q_bD_tpi_d=branch["Q_bD"][0], Q_bF_tpi_d=branch["Q_bF"][0],
        p_tpi_d=branch["p"][0],
        surv_tpi_d=1.0,   # reneging repriced the claim; it never haircuts it
    )


def _branch_conv(old_in, new_in, keys, ss, cal):
    # CONVERGENCE GAUGE FOR ONE BRANCH: PRICE, RETURN, AND KERNEL GAPS.
    q_key, rk_key, om_key = keys
    return max(
        abs(np.asarray(old_in[q_key]) - new_in[q_key]) / ss["Q_bD_ss"],
        abs(np.asarray(old_in[rk_key]) - new_in[rk_key]),
        float(np.mean(np.abs(np.asarray(old_in[om_key]) - new_in[om_key])))
        / cal["beta_inter_D"],
    )


def solve_transition_risk(ss, cal, Z_D_path, Z_F_path, pi_D_path=None,
                          pi_tpi_D_path=None, s_tpi_D_path=None,
                          verbose=True, max_rounds=12, damp=0.5, tol=1e-3,
                          y0=None):
    # OUTER FIXED POINT: BASE PATH <-> DEFAULT BRANCH <-> TPI-RENEGED BRANCH.
    # Three independent "off" states, each nesting exactly:
    #   pi_D_path = 0      no priced default risk
    #   pi_tpi_D_path = 1  backstop never doubted (no TPI-reneged branch)
    #   s_tpi_D_path = 0   no mechanical CB purchases
    # The mechanical channel needs NO branch solve: bank_backward's price floor
    # reads only s_tpi_D, so "CB always buys, backstop never doubted" costs one
    # solve_transition call, the same as the risk-neutral path.
    T = cal["T"]
    pi_D_path = np.zeros(T) if pi_D_path is None else np.asarray(pi_D_path, dtype=float)
    pi_tpi_D_path = (np.ones(T) if pi_tpi_D_path is None
                     else np.asarray(pi_tpi_D_path, dtype=float))
    s_tpi_D_path = (np.zeros(T) if s_tpi_D_path is None
                    else np.asarray(s_tpi_D_path, dtype=float))

    def_live        = bool(np.any(pi_D_path))
    tpi_priced_live = bool(np.any(pi_tpi_D_path < 1.0))
    tpi_mech_live   = bool(np.any(s_tpi_D_path > 0.0))

    def _off_out(out):
        # TAG A SOLVE THAT NEEDED NO BRANCH AS A CONVERGED, BRANCHLESS RESULT.
        out["branch"] = None;      out["risk_D_inputs"] = None
        out["tpi_branch"] = None;  out["tpi_D_inputs"] = None
        out["pi_D"] = pi_D_path;   out["risk_converged"] = True
        return out

    if not def_live and not tpi_priced_live:
        return _off_out(solve_transition(
            ss, cal, Z_D_path, Z_F_path,
            def_price_D=(pi_D_path if tpi_mech_live else None),
            s_tpi_D=(s_tpi_D_path if tpi_mech_live else None),
            verbose=False, y0=y0))

    # one Jacobian cache per LIVE system kind, reused across re-solves of that kind
    jc_base   = {}
    jc_branch = {} if def_live else None
    jc_tpi    = {} if tpi_priced_live else None

    out = solve_transition(ss, cal, Z_D_path, Z_F_path,
                           def_price_D=pi_D_path, s_tpi_D=s_tpi_D_path,
                           verbose=False, y0=y0, jac_cache=jc_base)

    risk_in = None;  branch = None;      branch_y0 = None
    tpi_in  = None;  tpi_branch = None;  tpi_branch_y0 = None
    fixed_point_ok = True
    conv = np.inf
    for rd in range(1, max_rounds + 1):
        conv_def = conv_tpi = 0.0
        t_branch = t_tpi = 0.0

        if def_live:
            t0 = time.perf_counter()
            try:
                branch_new = solve_default_branch(out, ss, cal, tau=1, verbose=verbose,
                                                  y0=branch_y0, jac_cache=jc_branch)
            except RuntimeError as e:
                if branch is None:
                    raise
                print(f"  [risk_branch] round {rd}: default branch re-solve failed "
                      f"({e}); keeping the previous round's risk inputs.")
                fixed_point_ok = False
                branch_new = None
            if branch_new is not None:
                branch = branch_new
                t_branch = time.perf_counter() - t0
                branch_y0 = branch["y_vec"]
                new_in = make_risk_inputs(branch, out, ss, cal)
                if risk_in is None:
                    risk_in = new_in
                else:
                    risk_in = {k: (1 - damp) * np.asarray(risk_in[k]) + damp * np.asarray(new_in[k])
                               for k in new_in}
                conv_def = _branch_conv(risk_in, new_in,
                                        ("Q_bD_d", "rk_d_D", "Omega_d_D"), ss, cal)

        if tpi_priced_live:
            t0 = time.perf_counter()
            try:
                tpi_branch_new = solve_tpi_branch(out, ss, cal, tau=1, verbose=verbose,
                                                  y0=tpi_branch_y0, jac_cache=jc_tpi)
            except RuntimeError as e:
                if tpi_branch is None:
                    raise
                print(f"  [risk_branch] round {rd}: TPI branch re-solve failed "
                      f"({e}); keeping the previous round's TPI inputs.")
                fixed_point_ok = False
                tpi_branch_new = None
            if tpi_branch_new is not None:
                tpi_branch = tpi_branch_new
                t_tpi = time.perf_counter() - t0
                tpi_branch_y0 = tpi_branch["y_vec"]
                new_tpi_in = make_tpi_inputs(tpi_branch, out, ss, cal)
                if tpi_in is None:
                    tpi_in = new_tpi_in
                else:
                    tpi_in = {k: (1 - damp) * np.asarray(tpi_in[k]) + damp * np.asarray(new_tpi_in[k])
                              for k in new_tpi_in}
                conv_tpi = _branch_conv(tpi_in, new_tpi_in,
                                        ("Q_bD_tpi_d", "rk_tpi_d_D", "Omega_tpi_D"), ss, cal)

        conv = max(conv_def, conv_tpi)

        # risk_D must be non-None whenever tpi_priced_live (bank_backward needs
        # risk_mode for tpi_mode) even when def_live is False; pi = 0 then makes
        # the default term contribute exactly zero weight whatever it holds
        if risk_in is not None:
            risk_D = dict(pi=pi_D_path, **{k: np.asarray(v) if np.ndim(v) else float(v)
                                           for k, v in risk_in.items()})
        else:
            risk_D = dict(pi=pi_D_path,
                          Omega_d_D=np.full(T, cal["beta_inter_D"]),
                          Omega_d_F=np.full(T, cal["beta_inter_F"]),
                          rk_d_D=0.0, rk_d_F=0.0, Q_bD_d=1.0, Q_bF_d=1.0, p_d=1.0,
                          surv_d=cal["recovery_rate_D"])

        tpi_D = (dict(pi=pi_tpi_D_path, **{k: np.asarray(v) if np.ndim(v) else float(v)
                                           for k, v in tpi_in.items()})
                 if tpi_in is not None else None)

        t0 = time.perf_counter()
        out = solve_transition(
            ss, cal, Z_D_path, Z_F_path,
            def_price_D=pi_D_path,
            verbose=False, y0=out["y_vec"], risk_D=risk_D, tpi_D=tpi_D,
            s_tpi_D=s_tpi_D_path, jac_cache=jc_base,
        )

        if verbose:
            print(f"  risk round {rd}: conv={conv:.2e}  Q_bD[0]={out['Q_bD'][0]:.4f}"
                  f"  [def {t_branch:.0f}s, tpi {t_tpi:.0f}s, "
                  f"base {time.perf_counter() - t0:.0f}s]")
        if conv < tol and rd >= 2:
            if verbose:
                print(f"  risk channel converged in {rd} rounds.")
            break
    else:
        if verbose:
            print(f"  risk channel: max_rounds={max_rounds} reached (conv={conv:.2e}).")

    out["branch"] = branch
    out["risk_D_inputs"] = risk_in
    out["tpi_branch"] = tpi_branch
    out["tpi_D_inputs"] = tpi_in
    out["pi_D"] = pi_D_path
    out["risk_converged"] = bool(fixed_point_ok and conv < tol)
    return out


def bond_decomposition(out, ss, cal):
    # SPLIT THE D-BOND EXCESS RETURN INTO DEFAULT COMPENSATION + RISK PREMIUM + LIQUIDITY.
    # Exact per-period identity: payoff^nd/Q - 1 - rdep = defcomp + risk + lambda*mu/Omega.
    # The risk premium is 0 in risk-neutral mode and positive when Omega^d > Omega^nd.
    db  = cal["delta_b_D"]
    rec = cal["recovery_rate_D"]
    Q   = np.asarray(out["Q_bD"])
    Q_next  = np.append(Q[1:], ss["Q_bD_ss"])
    pi_next = np.append(np.asarray(out["def_price_D"])[1:], 0.0)

    payoff_nd = db + (1 - db) * Q_next
    if out.get("risk_D_inputs") is not None:
        Q_d = float(np.asarray(out["risk_D_inputs"]["Q_bD_d"]))
        surv_d = float(np.asarray(out["risk_D_inputs"].get("surv_d", rec)))
        payoff_d = surv_d * (db + (1 - db) * Q_d)
    else:
        payoff_d = rec * (db + (1 - db) * Q_next)   # surv-form equivalent

    Epay = (1 - pi_next) * payoff_nd + pi_next * payoff_d
    rdep = np.asarray(out["rdep_D"])
    ic   = np.asarray(out["ic_spread_bD_D"])
    bk   = ss["ss_bank_D"]
    ic_ss = bk["lambda_bD"] * bk["mu_ss"] / bk["Omega_ss"]

    y_sov = db / Q - db
    y_ss  = db / ss["Q_bD_ss"] - db

    return dict(
        total_yield=4e4 * (y_sov - y_ss),
        defcomp=4e4 * ((payoff_nd - Epay) / Q),
        risk=4e4 * (Epay / Q - 1 - rdep - ic),
        liquidity=4e4 * (ic - ic_ss),
        promised_excess=4e4 * (payoff_nd / Q - 1 - rdep - ic_ss),
    )
