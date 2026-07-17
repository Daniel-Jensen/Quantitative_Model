# **Bocola (2016) risk channel via a representative post-default branch.**
# Bankers weight two branches: no-default (prob 1-π, base path) and default
# (prob π, a haircut). Pairing HIGH marginal valuations Ω^d with LOW default-
# state payoffs puts a covariance premium on bonds and capital (precautionary
# deleveraging). π_t is EXOGENOUS (Bocola eqs. 11-12, s-shock analog): an
# input path, never a function of debt. One representative branch (default at
# τ*=1) is reused at every date; π≡0 nests the risk-neutral model exactly.
# The feared event is a PURE HAIRCUT (recovery = recovery_rate_D on the whole
# claim, Bocola's D); the default-state recession must arise endogenously via
# bank balance sheets. Diagnostic flags, all OFF (=0) at baseline:
# def_output_cost_D (Arellano output cost), def_capital_quality_D (GK ξ_K),
# recap_share_D (equity injection financed by issuance).
# Approximations: Λ^nd ≡ beta_inter on the base path; Λ^d from the household-
# SDF proxy on aggregate income; the household deposit Euler is π-blind
# (faithful to Bocola: deposits are riskless there too).
import numpy as np

from transition import solve_transition


def extract_init_state(out, ss, cal, tau):
    # **Initial-conditions dict for a transition starting at base period tau (uses tau-1 objects).**
    # Portfolio shares on ACTUAL net worth (n_IC only coincides when the IC binds).
    s = tau - 1
    n_D_s = out["n_D"][s];  n_F_s = out["n_F"][s]
    p_s = out["p"][s]
    return dict(
        D_D=out["D_start_D"][tau], D_F=out["D_start_F"][tau],
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
    # **Solve the ONE feared event: a pure haircut on the outstanding stock.**
    # def_real_D[0]=1 with recovery recovery_rate_D (Bocola's D=1-rec). The
    # Bohn rule re-anchors to the surviving stock (taxing the pre-haircut stock
    # would be a ~31%-of-GDP artifact). Diagnostic flags off at baseline:
    # output cost, ξ_K capital-quality loss, recap. If the direct solve stalls
    # (a ~54% net-worth wipeout is far outside the shifted-base Newton basin),
    # a NUMERICAL continuation chains warm starts down a recap-share ladder —
    # intermediate solves are auxiliary and discarded; the returned branch is
    # always the calibrated event. Infeasible after continuation → RAISES.
    T = cal["T"]
    init = extract_init_state(out, ss, cal, tau)

    rec_D       = cal["recovery_rate_D"]
    xi_K        = cal.get("def_capital_quality_D", 0.0)
    recap_share = cal.get("recap_share_D", 0.0)
    cost        = cal.get("def_output_cost_D", 0.0)
    rho_c       = cal.get("def_output_rho_D", 0.9)

    Z_D_b = _shift_path(out["Z_D"], tau, T)
    if cost > 0.0:
        Z_D_b = Z_D_b * (1.0 - cost * rho_c ** np.arange(T))
    Z_F_b = _shift_path(out["Z_F"], tau, T)

    if xi_K > 0.0:
        init["quality0_D"] = 1.0 - xi_K   # GK capital-quality loss at h=0

    def_real_D = np.zeros(T)
    def_real_D[0] = 1.0
    init["b_anchor_D"] = init["b_gov0_D"] * rec_D   # Bohn on the surviving stock

    s = tau - 1
    bond_exposure = out["Q_bD"][s] * out["b_D_D"][s]

    def _attempt(share, y_s):
        # **One branch solve at a given recap share; None if the solver stalls.**
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

    y_start = y0 if y0 is not None else _shifted_y0(out, tau, T, xi_K)
    branch = _attempt(recap_share, y_start)

    if branch is None:
        # continuation: larger recaps cushion the wipeout, shrinking the jump
        # from the shifted-base warm start; chain y down to the calibrated share
        if verbose:
            print("  [risk_branch] direct solve stalled; recap-share continuation")
        y_c = y_start
        for share in (0.5, 0.35, 0.2, 0.1, 0.05):
            if share <= recap_share:
                continue
            aux = _attempt(share, y_c)
            if aux is None:
                continue
            y_c = aux["y_vec"]
        branch = _attempt(recap_share, y_c)

    if branch is None:
        raise RuntimeError(
            f"default branch (pure haircut, recovery={rec_D}) infeasible even "
            "after recap continuation. Diagnostic flags available: "
            "def_output_cost_D, def_capital_quality_D, recap_share_D.")

    branch["recap_D_path"] = (init["recap_D_path"].copy()
                              if recap_share > 0.0 else np.zeros(T))
    if verbose:
        print(f"  [risk_branch] branch solved: n_D(0)/n_ss = "
              f"{branch['n_D'][0] / ss['ss_bank_D']['n_ss']:.3f}, "
              f"Q_bD(0) = {branch['Q_bD'][0]:.4f}, rk_D(0) = {branch['rk_D'][0]:+.4f}")
    return branch


def make_risk_inputs(branch, base, ss, cal):
    # **Branch objects → risk_D dict for bank_backward (Ω^d, branch prices/returns, surv_d).**
    # Ω^d = Λ^d·[f + (1-f)·α^d(0)] (Bocola kernel weights). Λ^d = household-SDF
    # proxy on aggregate income: marginal value high where branch output is low
    # relative to the base continuation at t+1. If the branch is NOT a recession
    # the proxy is sign-gated to Λ^d = β (no SDF-side premium) with a warning.
    T = cal["T"]
    Y_d_D  = branch["Y_D"][0]
    Y_d_F  = branch["Y_F"][0]
    Y_nd_D = _shift_path(base["Y_D"], 1, T)
    Y_nd_F = _shift_path(base["Y_F"], 1, T)
    ratio_D = Y_d_D / Y_nd_D
    if np.any(ratio_D >= 1.0):
        print("  [risk_branch] WARNING: income-SDF sign gate tripped — "
              f"branch Y_D(0) = {Y_d_D:.4f} is not below the base path "
              "(default state not a recession). Λ^d falls back to β.")
        Lam_d_D = np.full(T, cal["beta_inter_D"])
        Lam_d_F = np.full(T, cal["beta_inter_F"])
    else:
        Lam_d_D = cal["beta_inter_D"] * ratio_D ** (-cal["sigma_D"])
        Lam_d_F = cal["beta_inter_F"] * (Y_d_F / Y_nd_F) ** (-cal["sigma_F"])

    alpha_d_D0 = branch["alpha_D"][0]
    alpha_d_F0 = branch["alpha_F"][0]
    Omega_d_D = Lam_d_D * (cal["f_D"] + (1 - cal["f_D"]) * alpha_d_D0)
    Omega_d_F = Lam_d_F * (cal["f_F"] + (1 - cal["f_F"]) * alpha_d_F0)

    return dict(
        Omega_d_D=Omega_d_D, Omega_d_F=Omega_d_F,
        rk_d_D=branch["rk_D"][0], rk_d_F=branch["rk_F"][0],
        Q_bD_d=branch["Q_bD"][0], Q_bF_d=branch["Q_bF"][0],
        p_d=branch["p"][0],
        surv_d=cal["recovery_rate_D"],   # haircut on the whole claim
    )


def solve_transition_risk(ss, cal, Z_D_path, Z_F_path, pi_D_path=None,
                          verbose=True, max_rounds=12, damp=0.5, tol=1e-3,
                          y0=None):
    # **Outer fixed point base path ↔ representative branch at exogenous π_t.**
    # Round 0: risk-neutral base (π priced linearly). Round k: solve branch →
    # damped risk inputs → re-solve base with two-branch pricing. π_D ≡ 0
    # returns the risk-neutral path immediately (exact nesting).
    T = cal["T"]
    if pi_D_path is None:
        pi_D_path = np.zeros(T)
    pi_D_path = np.asarray(pi_D_path, dtype=float)

    # one Jacobian cache per system kind (base / branch); reused across re-solves
    jc_base = {}
    jc_branch = {}

    out = solve_transition(ss, cal, Z_D_path, Z_F_path,
                           def_price_D=pi_D_path,
                           verbose=False, y0=y0, jac_cache=jc_base)

    if not np.any(pi_D_path):
        out["branch"] = None
        out["risk_D_inputs"] = None
        out["pi_D"] = pi_D_path
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

        risk_D = dict(pi=pi_D_path, **{k: np.asarray(v) if np.ndim(v) else float(v)
                                       for k, v in risk_in.items()})

        _t0 = _time.perf_counter()
        out = solve_transition(
            ss, cal, Z_D_path, Z_F_path,
            def_price_D=pi_D_path,
            verbose=False, y0=out["y_vec"], risk_D=risk_D,
            jac_cache=jc_base,
        )
        _t_base = _time.perf_counter() - _t0

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

    out["branch"] = branch
    out["risk_D_inputs"] = risk_in
    out["pi_D"] = pi_D_path
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
    if out.get("risk_D_inputs") is not None:
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
