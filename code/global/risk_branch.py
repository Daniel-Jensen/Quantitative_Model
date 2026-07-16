# **Bocola (2016) risk channel via a representative post-default branch.**
# Bankers weight two branches: no-default (prob 1-π, base path) and default
# (prob π, a haircut + recession). Pairing HIGH marginal valuations Ω^d with
# LOW default-state payoffs puts a covariance premium on bonds and capital
# (precautionary deleveraging). One representative branch (default at τ*=1) is
# reused at every date; π≡0 nests the risk-neutral model exactly. Approximations:
# Λ^nd≡beta_inter on the base path; household deposit Euler is π-blind.
import numpy as np

from government import ck_default_prob
from transition import solve_transition, solve_transition_ck


def extract_init_state(out, ss, cal, tau):
    # **Initial-conditions dict for a transition starting at base period tau (uses tau-1 objects).**
    s = tau - 1
    n_IC_D = out["n_IC_D"][s];  n_IC_F = out["n_IC_F"][s]
    p_s = out["p"][s]
    return dict(
        D_D=out["D_start_D"][tau], D_F=out["D_start_F"][tau],
        bank_D=dict(
            n_prev=out["n_D"][s],
            kappa_prev=out["Q_D"][s] * out["Kap_D"][s] / n_IC_D,
            phi_bdom_prev=out["Q_bD"][s] * out["b_D_D"][s] / n_IC_D,
            phi_bfor_prev=p_s * out["Q_bF"][s] * out["b_F_D"][s] / n_IC_D,
            rdep_prev=out["rdep_D"][s],
        ),
        bank_F=dict(
            n_prev=out["n_F"][s],
            kappa_prev=out["Q_F"][s] * out["Kap_F"][s] / n_IC_F,
            phi_bdom_prev=out["Q_bF"][s] * out["b_F_F"][s] / n_IC_F,
            phi_bfor_prev=out["Q_bD"][s] * out["b_D_F"][s] / (p_s * n_IC_F),
            rdep_prev=out["rdep_F"][s],
        ),
        b_gov0_D=out["b_gov_D"][s], b_gov0_F=out["b_gov_F"][s],
        Kap_lag_D=out["Kap_D"][s], Kap_lag_F=out["Kap_F"][s],
        Q_lag_D=out["Q_D"][s], Q_lag_F=out["Q_F"][s],
        Q_bD_lag=out["Q_bD"][s], Q_bF_lag=out["Q_bF"][s],
        p_lag=p_s, P_lag_D=out["P_CES_D"][s], P_lag_F=out["P_CES_F"][s],
    )


def _shift_path(x, tau, T):
    # **Shift a path forward by tau, padded with its final value.**
    idx = np.minimum(np.arange(T) + tau, T - 1)
    return np.asarray(x)[idx]


def _shifted_y0(out, tau, T, xi_K=0.0, rho_rebuild=0.975):
    # **Branch Newton warm start: base solution shifted by tau (Kap scaled by the quality profile).**
    # xi_K>0: scaling the Kap_D block toward the post-destruction path avoids an
    # implied one-quarter rebuild that would spike the Jermann price and stall
    # every probe on the penalty wall. rho_rebuild≈0.975 ≈ Jermann K half-life.
    blocks = [out["N_D"], out["N_F"], out["Kap_D"], out["Kap_F"],
              out["rdep_D"], out["rdep_F"], out["p"]]
    y = np.concatenate([_shift_path(b, tau, T) for b in blocks])
    if xi_K > 0.0:
        y[2*T:3*T] *= 1.0 - xi_K * rho_rebuild ** np.arange(T)
    return y


def solve_default_branch(out, ss, cal, tau=1, verbose=False, y0=None,
                         jac_cache=None):
    # **Solve ONE fixed feared event: haircut + capital-quality loss + government recap.**
    # scale = branch_haircut_scale (1.0 = full PSI). ξ_K destroys branch capital
    # (stops it being the safe haven). recap (on when recap_share_D>0) makes the
    # full event feasible. Infeasible → RAISE unless branch_use_ladder re-enables
    # the scale search. branch["rescue_mode"]/["haircut_scale"] record what ran.
    T = cal["T"]
    init = extract_init_state(out, ss, cal, tau)

    rec_D = cal["recovery_rate_D"]
    xi_K = cal.get("def_capital_quality_D", 0.0)
    recap_share = cal.get("recap_share_D", 0.0)

    # output cost of default (Arellano 2008) makes the branch a recession
    cost = cal.get("def_output_cost_D", 0.0)
    rho_c = cal.get("def_output_rho_D", 0.9)
    Z_D_b = _shift_path(out["Z_D"], tau, T) * (1.0 - cost * rho_c ** np.arange(T))
    Z_F_b = _shift_path(out["Z_F"], tau, T)

    init["quality0_D"] = 1.0 - xi_K   # GK capital-quality loss at h=0 (branch only)
    s = tau - 1
    bond_exposure_D = out["Q_bD"][s] * out["b_D_D"][s]   # for ex-ante recap sizing

    def _attempt(scale, recap_on, y_start, hybr_factor=100.0):
        # **One branch solve at a given haircut scale; None if infeasible.**
        def_real_D = np.zeros(T)
        def_real_D[0] = scale
        # re-anchor the Bohn rule to the post-haircut stock (else default is expansionary)
        init["b_anchor_D"] = init["b_gov0_D"] * (1.0 - scale * (1.0 - rec_D))
        if recap_on and recap_share > 0.0:
            recap = np.zeros(T)
            recap[0] = recap_share * scale * (1.0 - rec_D) * bond_exposure_D
            init["recap_D_path"] = recap
        else:
            init.pop("recap_D_path", None)
        try:
            return solve_transition(
                ss, cal, Z_D_b, Z_F_b,
                def_real_D=def_real_D,
                verbose=False, y0=y_start, init=init,
                jac_cache=jac_cache, hybr_factor=hybr_factor,
                accept_tol=1e-9,
            )
        except (RuntimeError, ValueError) as e:
            if verbose:
                print(f"    [branch] scale={scale:.3f} recap={recap_on} "
                      f"infeasible ({e})")
            return None

    def _climb(scales, y_start, hybr_factor=100.0):
        # **Opt-in ladder: try every scale, keep the LARGEST feasible event.**
        # Feasibility is non-monotone (small events boom-infeasible, large ones
        # crunch-infeasible), so never stop at the first failure.
        b, s_used, y = None, 0.0, y_start
        for scale in scales:
            cand = _attempt(scale, True, y, hybr_factor)
            if cand is None:
                continue
            b, s_used = cand, scale
            y = b["y_vec"]
        return b, s_used

    # one deterministic solve of the fixed event (recap on when calibrated)
    scale = float(cal.get("branch_haircut_scale", 1.0))
    recap_on = recap_share > 0.0
    y_start = y0 if y0 is not None else _shifted_y0(out, tau, T, xi_K)

    branch = _attempt(scale, recap_on, y_start)
    scale_used = scale
    mode = (("full" if scale >= 1.0 else f"scale{scale:.2f}")
            + ("+recap" if recap_on else ""))

    # opt-in fallback: the feasibility ladder (off by default → single deterministic solve)
    if branch is None and cal.get("branch_use_ladder", False):
        full_ladder = (0.075, 0.15, 0.3, 0.5, 0.75, 1.0)
        branch, scale_used = _climb(full_ladder, y_start)
        if branch is None:
            if verbose:
                print("    [ladder] all probes stalled; retrying with "
                      "small trust region (hybr_factor=0.1)")
            branch, scale_used = _climb(full_ladder,
                                        _shifted_y0(out, tau, T, xi_K),
                                        hybr_factor=0.1)
        if branch is not None:
            mode = (f"ladder({scale_used:.3f})"
                    + ("+recap" if recap_on else ""))

    if branch is None:
        raise RuntimeError(
            f"default branch infeasible at haircut scale={scale:.2f} "
            f"(recap={'on' if recap_on else 'off'}, share={recap_share:.2f}). "
            "Raise recap_share_D, lower branch_haircut_scale, or set "
            "branch_use_ladder=True to search for the largest feasible event.")

    branch["haircut_scale"] = scale_used
    branch["rescue_mode"] = mode
    branch["recap_D_path"] = init.get("recap_D_path", np.zeros(T)).copy() \
        if init.get("recap_D_path") is not None else np.zeros(T)
    if verbose and scale_used < 1.0:
        print(f"  [risk_branch] priced event = partial restructuring "
              f"(haircut {scale_used * (1 - rec_D):.0%}; full event "
              "infeasible even with recap)")
    b_post = init["b_gov0_D"] * (1.0 - scale_used * (1.0 - rec_D))
    if b_post / ss["ss_firm_D"]["Y_ss"] >= cal["b_ck_low_D"]:
        print(f"  [risk_branch] WARNING: post-default debt b/Y = "
              f"{b_post / ss['ss_firm_D']['Y_ss']:.2f} ≥ b_ck_low — branch "
              "not absorbing; risk inputs understate continuation risk.")
    if verbose:
        print(f"  [risk_branch] branch solved [{mode}]: n_D(0)/n_ss = "
              f"{branch['n_D'][0] / ss['ss_bank_D']['n_ss']:.3f}, "
              f"Q_bD(0) = {branch['Q_bD'][0]:.4f}, rk_D(0) = {branch['rk_D'][0]:+.4f}")
    return branch


def make_risk_inputs(branch, base, ss, cal):
    # **Branch objects → risk_D dict for bank_backward (Ω^d, branch prices/returns, surv_d).**
    # Ω^d = Λ^d·[(1-f)+f·α^d(0)]. sdf_mode: "income" (Euler loading on branch
    # output vs base t+1, sign-gated) | "empirical" (β·kappa_d) | "model" (GHH
    # composites, wrong-signed until the union-deposit fix).
    T = cal["T"]
    mode = cal.get("sdf_mode", "empirical")
    if mode == "income":
        # marginal value high where branch output is low, vs base continuation at t+1
        Y_d_D = branch["Y_D"][0]
        Y_d_F = branch["Y_F"][0]
        Y_nd_D = _shift_path(base["Y_D"], 1, T)
        Y_nd_F = _shift_path(base["Y_F"], 1, T)
        ratio_D = Y_d_D / Y_nd_D
        if np.any(ratio_D >= 1.0):
            print("  [risk_branch] WARNING: income-SDF sign gate tripped — "
                  f"branch Y_D(0) = {Y_d_D:.4f} is not below the base path "
                  "(default state not a recession).  Falling back to "
                  "sdf_mode='empirical' for this round.")
            kappa = cal.get("kappa_d", 1.0)
            Lam_d_D = np.full(T, cal["beta_inter_D"] * kappa)
            Lam_d_F = np.full(T, cal["beta_inter_F"] * kappa)
        else:
            Lam_d_D = cal["beta_inter_D"] * ratio_D ** (-cal["sigma_D"])
            Lam_d_F = cal["beta_inter_F"] * (Y_d_F / Y_nd_F) ** (-cal["sigma_F"])
    elif mode == "empirical":
        kappa = cal.get("kappa_d", 1.0)
        Lam_d_D = np.full(T, cal["beta_inter_D"] * kappa)
        Lam_d_F = np.full(T, cal["beta_inter_F"] * kappa)
    else:
        # model-consistent SDF from GHH composites x = C - v(N) (branch vs base t+1)
        x_d_D = branch["C_D"][0] - branch["vN_D"][0]
        x_d_F = branch["C_F"][0] - branch["vN_F"][0]
        x_nd_D = _shift_path(base["C_D"] - base["vN_D"], 1, T)
        x_nd_F = _shift_path(base["C_F"] - base["vN_F"], 1, T)
        Lam_d_D = cal["beta_inter_D"] * (x_d_D / x_nd_D) ** (-cal["sigma_D"])
        Lam_d_F = cal["beta_inter_F"] * (x_d_F / x_nd_F) ** (-cal["sigma_F"])

    alpha_d_D0 = branch["alpha_D"][0]
    alpha_d_F0 = branch["alpha_F"][0]
    Omega_d_D = Lam_d_D * ((1 - cal["f_D"]) + cal["f_D"] * alpha_d_D0)
    Omega_d_F = Lam_d_F * ((1 - cal["f_F"]) + cal["f_F"] * alpha_d_F0)

    return dict(
        Omega_d_D=Omega_d_D, Omega_d_F=Omega_d_F,
        rk_d_D=branch["rk_D"][0], rk_d_F=branch["rk_F"][0],
        Q_bD_d=branch["Q_bD"][0], Q_bF_d=branch["Q_bF"][0],
        p_d=branch["p"][0],
        surv_d=1.0 - branch["haircut_scale"] * (1.0 - cal["recovery_rate_D"]),
    )


def solve_transition_ck_risk(ss, cal, Z_D_path, Z_F_path,
                             sunspot_D_path=None, sunspot_F_path=None,
                             verbose=True, max_rounds=6, damp=0.5, tol=1e-3,
                             y0=None):
    # **Outer fixed point base path ↔ representative branch (Cole-Kehoe risk-only).**
    # Round 0: risk-off base. Round k: solve branch → damped risk inputs →
    # re-solve base with two-branch pricing → refresh zone indicator. Converges
    # on the branch objects; a final zone-consistency loop is a no-op in-zone.
    T = cal["T"]
    Y_ss_D = ss["ss_firm_D"]["Y_ss"]
    Y_ss_F = ss["ss_firm_F"]["Y_ss"]
    if sunspot_D_path is None:
        sunspot_D_path = np.zeros(T)
    if sunspot_F_path is None:
        sunspot_F_path = np.zeros(T)

    def _zone(b_bop, Y_ss, sun, ctry):
        # **Crisis-zone indicator re-evaluated on a solved debt path.**
        return np.array([ck_default_prob(b_bop[t], Y_ss, cal, sun[t], ctry)
                         for t in range(T)])

    chi_tilt = cal.get("chi_tilt", 1.0)

    # one Jacobian cache per system kind (base / branch); reused across re-solves
    jc_base = {}
    jc_branch = {}

    out = solve_transition_ck(ss, cal, Z_D_path, Z_F_path,
                              sunspot_D_path=sunspot_D_path,
                              sunspot_F_path=sunspot_F_path,
                              verbose=False, y0=y0, jac_cache=jc_base)
    def_price_D = np.asarray(out["def_price_D"])
    def_price_F = np.asarray(out["def_price_F"])
    def_price_D_used = def_price_D
    def_price_F_used = def_price_F

    # No priced risk anywhere → the branch is irrelevant (pi ≡ 0 nests the
    # risk-neutral model exactly), so return the risk-off path and skip the
    # branch solve, which from the healthy SS state needs the (off-by-default)
    # ladder as a homotopy.
    if not np.any(def_price_D) and not np.any(def_price_F):
        out["branch"] = None
        out["risk_D_inputs"] = None
        out["sunspot_D"] = sunspot_D_path
        out["sunspot_F"] = sunspot_F_path
        out["risk_converged"] = True
        return out

    import time as _time
    risk_in = None
    branch = None
    branch_y0 = None
    fixed_point_ok = True
    conv = np.inf
    for rd in range(1, max_rounds + 1):
        _t0 = _time.perf_counter()
        try:
            branch_new = solve_default_branch(out, ss, cal, tau=1,
                                              verbose=verbose,
                                              y0=branch_y0,
                                              jac_cache=jc_branch)
        except RuntimeError as e:
            if branch is None:
                raise
            print(f"  [risk_branch] round {rd}: branch re-solve failed ({e}); "
                  "keeping previous round's risk inputs.")
            fixed_point_ok = False
            break
        branch = branch_new
        _t_branch = _time.perf_counter() - _t0
        branch_y0 = branch["y_vec"]
        new_in = make_risk_inputs(branch, out, ss, cal)

        if risk_in is None:
            risk_in = new_in
        else:
            risk_in = {k: (1 - damp) * np.asarray(risk_in[k]) + damp * np.asarray(new_in[k])
                       for k in new_in}
        conv = max(
            abs(np.asarray(risk_in["Q_bD_d"]) - new_in["Q_bD_d"]) / ss["Q_bD_ss"],
            abs(np.asarray(risk_in["rk_d_D"]) - new_in["rk_d_D"]),
            float(np.mean(np.abs(np.asarray(risk_in["Omega_d_D"]) - new_in["Omega_d_D"])))
            / cal["beta_inter_D"],
        )

        pi = np.minimum(chi_tilt * def_price_D, 1.0)
        risk_D = dict(pi=pi, **{k: np.asarray(v) if np.ndim(v) else float(v)
                                for k, v in risk_in.items()})

        _t0 = _time.perf_counter()
        def_price_D_used = def_price_D
        def_price_F_used = def_price_F
        out = solve_transition(
            ss, cal, Z_D_path, Z_F_path,
            def_price_D=def_price_D, def_price_F=def_price_F,
            verbose=False, y0=out["y_vec"], risk_D=risk_D,
            jac_cache=jc_base,
        )
        _t_base = _time.perf_counter() - _t0

        def_price_D = _zone(out["b_gov_bop_D"], Y_ss_D, sunspot_D_path, "D")
        def_price_F = _zone(out["b_gov_bop_F"], Y_ss_F, sunspot_F_path, "F")

        if verbose:
            print(f"  risk round {rd}: conv={conv:.2e}  Q_bD[0]={out['Q_bD'][0]:.4f}"
                  f"  Q_bD_d={float(np.asarray(risk_in['Q_bD_d'])):.4f}"
                  f"  rk_d_D={float(np.asarray(risk_in['rk_d_D'])):+.4f}"
                  f"  [branch {_t_branch:.0f}s, base {_t_base:.0f}s]")
        if conv < tol and rd >= 2:
            if verbose:
                print(f"  risk channel converged in {rd} rounds.")
            break
    else:
        if verbose:
            print(f"  risk channel: max_rounds={max_rounds} reached (conv={conv:.2e}).")

    # final zone-consistency: re-solve if the refreshed indicator moved (no-op in-zone)
    zone_tol = cal.get("ck_tol", 1e-8)
    for _ in range(5):
        zone_err = max(np.max(np.abs(def_price_D - def_price_D_used)),
                       np.max(np.abs(def_price_F - def_price_F_used)))
        if zone_err <= zone_tol:
            break
        if verbose:
            print(f"  [risk_branch] zone indicator moved in the final round "
                  f"(err={zone_err:.2e}); re-solving for consistency.")
        pi = np.minimum(chi_tilt * def_price_D, 1.0)
        risk_D = dict(pi=pi, **{k: np.asarray(v) if np.ndim(v) else float(v)
                                for k, v in risk_in.items()})
        def_price_D_used = def_price_D
        def_price_F_used = def_price_F
        out = solve_transition(
            ss, cal, Z_D_path, Z_F_path,
            def_price_D=def_price_D, def_price_F=def_price_F,
            verbose=False, y0=out["y_vec"], risk_D=risk_D,
            jac_cache=jc_base,
        )
        def_price_D = _zone(out["b_gov_bop_D"], Y_ss_D, sunspot_D_path, "D")
        def_price_F = _zone(out["b_gov_bop_F"], Y_ss_F, sunspot_F_path, "F")
    else:
        print("  [risk_branch] WARNING: zone indicator did not settle in the "
              "final consistency step; returned path and def_price disagree.")
        fixed_point_ok = False

    out["branch"] = branch
    out["risk_D_inputs"] = risk_in
    out["def_price_D"] = def_price_D
    out["def_price_F"] = def_price_F
    out["sunspot_D"] = sunspot_D_path
    out["sunspot_F"] = sunspot_F_path
    out["risk_converged"] = bool(fixed_point_ok and conv < tol)
    return out


def bond_decomposition(out, ss, cal):
    # **Split the D-bond excess return into default compensation + risk premium + liquidity.**
    # Exact per-period identity: payoff^nd/Q-1-rdep = defcomp + risk + λμ/Ω̃.
    # Risk premium is 0 in risk-neutral mode, positive when Ω^d>Ω^nd depresses Q.
    T   = cal["T"]
    db  = cal["delta_b_D"]
    rec = cal["recovery_rate_D"]
    Q   = np.asarray(out["Q_bD"])
    Q_next   = np.append(Q[1:], ss["Q_bD_ss"])
    pi_next  = np.append(np.asarray(out["def_price_D"])[1:], 0.0)

    payoff_nd = db + (1 - db) * Q_next
    if "risk_D_inputs" in out:
        pi_next = np.minimum(cal.get("chi_tilt", 1.0) * pi_next, 1.0)
        Q_d = float(np.asarray(out["risk_D_inputs"]["Q_bD_d"]))
        surv_d = float(np.asarray(out["risk_D_inputs"].get("surv_d", rec)))
        payoff_d = surv_d * (db + (1 - db) * Q_d)
    else:
        payoff_d = rec * (db + (1 - db) * Q_next)   # surv-form equivalent

    Epay  = (1 - pi_next) * payoff_nd + pi_next * payoff_d
    rdep  = np.asarray(out["rdep_D"])
    ic    = np.asarray(out["ic_spread_bD_D"])
    bk    = ss["ss_bank_D"]
    ic_ss = bk["lambda_bD"] * bk["mu_ss"] / bk["Omega_ss"]

    defcomp = (payoff_nd - Epay) / Q
    risk    = Epay / Q - 1 - rdep - ic
    y_sov   = db / Q - db
    y_ss    = db / ss["Q_bD_ss"] - db

    return dict(
        total_yield=4e4 * (y_sov - y_ss),
        defcomp=4e4 * defcomp,
        risk=4e4 * risk,
        liquidity=4e4 * (ic - ic_ss),
        promised_excess=4e4 * (payoff_nd / Q - 1 - rdep - ic_ss),
    )
