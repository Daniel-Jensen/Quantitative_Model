"""Bocola (2016) risk channel via a representative default branch.

The risk channel: bankers discount with the household SDF
Λ(S',S) = β·u_c'/u_c and carry state-contingent continuation values α(S').
News that default is more likely pairs HIGH marginal valuations (Λ^d, α^d —
the default state is a recession with depressed bank net worth) with LOW
asset returns → a covariance premium on BOTH sovereign bonds and capital
claims, i.e. precautionary deleveraging even when funding is cheap (Bocola's
LTRO result; the channel explains up to 45% of the crisis impact on firms'
borrowing costs in his estimation).

Perfect-foresight implementation (R2-lite): at each base date t the banker's
expectations average two branches —
  no-default (prob 1−π_{t+1}) : continue on the base path;
  default    (prob π_{t+1})   : haircut 1−recovery realized at t+1, economy
                                jumps to a post-default transition.
One REPRESENTATIVE branch is solved (default at τ*=1, from the impact
state); its h=0 objects (α^d, rk^d, Q^d, p^d, consumption composites) are
reused at every base date.  The branch is absorbing: post-default debt
(1−D)·b/Y falls below b_ck_low, so the branch carries no further default
risk.  Approximations (documented): Λ^nd ≡ beta_inter on the base path (the
channel lives in π·(Ω^d − Ω^nd)); branch state dependence across base dates
is second order; household expectations stay PERFECT-FORESIGHT on the base
path — the deposit Euler does not weight the default branch by π (deposits
are safe in rate terms even in the branch, so this omits only the
precautionary-savings response to default-state income risk; risk pricing
lives entirely in the bank block).  Extension: per-date branches.

Outer fixed point: base path ↔ branch (the branch starts from the base
impact state; the base backward pass uses branch objects), damped iteration
with warm starts.
"""
import numpy as np

from government import ck_default_prob
from transition import solve_transition, solve_transition_ck


# ─────────────────────────────────────────────────────────────────────────────
# State extraction and branch solve
# ─────────────────────────────────────────────────────────────────────────────

def extract_init_state(out, ss, cal, tau):
    """Initial-conditions dict for a transition starting at base period tau
    (uses only period tau−1 objects of the solved base path `out`)."""
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
    """x shifted forward by tau, padded with its final value."""
    idx = np.minimum(np.arange(T) + tau, T - 1)
    return np.asarray(x)[idx]


def _shifted_y0(out, tau, T):
    """Warm start for the branch Newton: base solution shifted by tau."""
    blocks = [out["N_D"], out["N_F"], out["Kap_D"], out["Kap_F"],
              out["rdep_D"], out["rdep_F"], out["p"]]
    return np.concatenate([_shift_path(b, tau, T) for b in blocks])


def solve_default_branch(out, ss, cal, tau=1, verbose=False, y0=None,
                         target_scale=None):
    """Post-default PF transition: the FULL haircut (1 − recovery) realized
    at branch period 0, starting from the base-path state entering period tau.

    The haircut feasibility ladder was removed (2026-07-13): with
    recovery_rate_D = 0.90 the event costs banks ~10% of their sovereign
    exposure and never threatens equity, so the largest-feasible-event
    search is moot.  haircut_scale ≡ 1.0.  (`target_scale` kept for caller
    compatibility; unused.)
    """
    T = cal["T"]
    init = extract_init_state(out, ss, cal, tau)

    rec_D = cal["recovery_rate_D"]

    # Canonical output cost of default (Arellano 2008): the default state is
    # a recession — this is what makes Λ^d, α^d high while returns are low.
    cost = cal.get("def_output_cost_D", 0.0)
    rho_c = cal.get("def_output_rho_D", 0.9)
    Z_D_b = _shift_path(out["Z_D"], tau, T) * (1.0 - cost * rho_c ** np.arange(T))
    Z_F_b = _shift_path(out["Z_F"], tau, T)

    def_real_D = np.zeros(T)
    def_real_D[0] = 1.0
    # No fiscal windfall: re-anchor the Bohn rule to the post-haircut stock,
    # otherwise φ·(b − b_ss) turns the haircut into large tax cuts and the
    # default state becomes expansionary (wrong-signed risk premium).
    init["b_anchor_D"] = init["b_gov0_D"] * rec_D

    # Warm start: previous round's branch solution if available, else the
    # shifted base path.  Retry with a small trust region if hybr stalls on
    # the alpha-explosion penalty wall (poisoned finite-difference Jacobian).
    y_start = y0 if y0 is not None else _shifted_y0(out, tau, T)
    branch = None
    for hybr_factor in (100.0, 0.1):
        try:
            branch = solve_transition(
                ss, cal, Z_D_b, Z_F_b,
                def_real_D=def_real_D,
                verbose=False, y0=y_start, init=init,
                try_krylov=False,
                hybr_factor=hybr_factor,
            )
            break
        except (RuntimeError, ValueError) as e:
            if verbose:
                print(f"    [branch] full-haircut solve failed "
                      f"(hybr_factor={hybr_factor}): {e}")
    if branch is None:
        raise RuntimeError("default branch infeasible at the full haircut — "
                           "bank equity wiped out (recap not modelled)")
    scale_used = 1.0
    branch["haircut_scale"] = scale_used
    # Absorbing-branch check at the ACTUAL event size
    b_post = init["b_gov0_D"] * (1.0 - scale_used * (1.0 - rec_D))
    if b_post / ss["ss_firm_D"]["Y_ss"] >= cal["b_ck_low_D"]:
        print(f"  [risk_branch] WARNING: post-default debt b/Y = "
              f"{b_post / ss['ss_firm_D']['Y_ss']:.2f} ≥ b_ck_low — branch "
              "not absorbing; risk inputs understate continuation risk.")
    if verbose:
        print(f"  [risk_branch] branch solved: n_D(0)/n_ss = "
              f"{branch['n_D'][0] / ss['ss_bank_D']['n_ss']:.3f}, "
              f"Q_bD(0) = {branch['Q_bD'][0]:.4f}, rk_D(0) = {branch['rk_D'][0]:+.4f}")
    return branch


# ─────────────────────────────────────────────────────────────────────────────
# Risk inputs for bank_backward
# ─────────────────────────────────────────────────────────────────────────────

def make_risk_inputs(branch, base, ss, cal):
    """Branch objects → risk_D dict for bank_backward (without 'pi').

    Ω^d_X[t] = Λ^d_X[t] · [(1−f_X) + f_X·α^d_X(0)], with the SDF ratio taken
    across branches at t+1.  Three sdf_mode variants (calibration.py):
      "income"    : Λ^d = beta_inter · (Y^d(0)/Y^nd_{t+1})^(−σ) — Euler-
                    consistent loading on aggregate OUTPUT.  Prices the
                    default state as a bad time whenever it is a recession
                    (Y^d < Y^nd), which holds even while the deposit-market
                    comovement problem makes branch CONSUMPTION rise — the
                    robust rep-agent proxy for the HA economy.  Sign-gated:
                    if the branch is not a recession the gate falls back to
                    the empirical mode loudly (never price default as good).
      "empirical" : Λ^d = beta_inter · kappa_d (free loading, Bocola's
                    asset-price-disciplined route).
      "model"     : Λ^d from GHH consumption composites x = C − v(N).
                    WRONG-SIGNED until the union-deposit-market fix (branch
                    C rises); kept for post-fix use.
    chi_tilt ≥ 1 (calibration, default 1) tilts the default probability
    pessimistically (EZ-lite dial; off by default).
    """
    T = cal["T"]
    mode = cal.get("sdf_mode", "empirical")
    if mode == "income":
        # Euler-consistent GDP loading: marginal value is high in states
        # where aggregate output is low.  Branch h=0 output vs the base
        # no-default continuation at t+1.
        Y_d_D = branch["Y_D"][0]
        Y_d_F = branch["Y_F"][0]
        Y_nd_D = _shift_path(base["Y_D"], 1, T)     # Y^nd_{t+1}
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
            # F banker loads on F output in the same D-default event; the
            # contagion recession there is small, so this leg stays near 1.
            Lam_d_F = cal["beta_inter_F"] * (Y_d_F / Y_nd_F) ** (-cal["sigma_F"])
    elif mode == "empirical":
        # Empirically disciplined default-state SDF loading (Bocola's own
        # route: E[Λ̂|default] proxied from asset prices, not model C).
        kappa = cal.get("kappa_d", 1.0)
        Lam_d_D = np.full(T, cal["beta_inter_D"] * kappa)
        Lam_d_F = np.full(T, cal["beta_inter_F"] * kappa)
    else:
        # Model-consistent SDF: branch h=0 composites vs the base no-default
        # continuation at t+1.  WARNING: wrong-signed while the comovement
        # problem makes branch consumption rise (see calibration.py note).
        x_d_D = branch["C_D"][0] - branch["vN_D"][0]
        x_d_F = branch["C_F"][0] - branch["vN_F"][0]
        x_nd_D = _shift_path(base["C_D"] - base["vN_D"], 1, T)   # x^nd_{t+1}
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
        # Survival factor of the PRICED event (haircut may be partial if the
        # full haircut is infeasible without bank recapitalization)
        surv_d=1.0 - branch["haircut_scale"] * (1.0 - cal["recovery_rate_D"]),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Outer loop: base path ↔ representative branch
# ─────────────────────────────────────────────────────────────────────────────

def solve_transition_ck_risk(ss, cal, Z_D_path, Z_F_path,
                             sunspot_D_path=None, sunspot_F_path=None,
                             verbose=True, max_rounds=6, damp=0.5, tol=1e-3,
                             y0=None):
    """Cole-Kehoe risk-only experiment WITH the Bocola risk channel.

    Rounds:
      0. risk-off base path (existing solve_transition_ck);
      k. representative default branch from the base impact state →
         risk inputs (damped) → re-solve the base path with two-branch
         pricing (warm-started) → update the crisis-zone def_price
         (both countries).
    Converges on the branch objects (Q^d, rk^d, mean Ω^d).  A final
    consistency step re-solves the base path if the refreshed zone
    indicator differs from the one used in the last solve (debt crossing
    a zone threshold in the last round) — a no-op when debt stays
    interior to the crisis zone.

    Returns the base-path dict plus 'branch' (the default branch),
    'risk_D_inputs', 'sunspot_D' and 'sunspot_F'.
    """
    T = cal["T"]
    Y_ss_D = ss["ss_firm_D"]["Y_ss"]
    Y_ss_F = ss["ss_firm_F"]["Y_ss"]
    if sunspot_D_path is None:
        sunspot_D_path = np.zeros(T)
    if sunspot_F_path is None:
        sunspot_F_path = np.zeros(T)

    def _zone(b_bop, Y_ss, sun, ctry):
        """Crisis-zone indicator re-evaluated on a solved debt path."""
        return np.array([ck_default_prob(b_bop[t], Y_ss, cal, sun[t], ctry)
                         for t in range(T)])

    chi_tilt = cal.get("chi_tilt", 1.0)

    # Round 0: risk-off base (also settles the crisis-zone indicator).
    # Pass y0 from a pre-solved risk-off run to skip the cold Newton start
    # (large sunspots need a homotopy, which main.py already performs).
    out = solve_transition_ck(ss, cal, Z_D_path, Z_F_path,
                              sunspot_D_path=sunspot_D_path,
                              sunspot_F_path=sunspot_F_path,
                              verbose=False, y0=y0)
    def_price_D = np.asarray(out["def_price_D"])
    def_price_F = np.asarray(out["def_price_F"])
    # Indicators actually used in the most recent base solve (round 0's CK
    # loop converges on its own indicator, so out is consistent here).
    def_price_D_used = def_price_D
    def_price_F_used = def_price_F

    import time as _time
    risk_in = None
    branch = None
    branch_y0 = None
    scale_star = None   # haircut_scale ≡ 1.0 (ladder removed); kept as warm-start hint
    fixed_point_ok = True
    conv = np.inf
    for rd in range(1, max_rounds + 1):
        _t0 = _time.perf_counter()
        try:
            branch_new = solve_default_branch(out, ss, cal, tau=1,
                                              verbose=verbose,
                                              y0=branch_y0,
                                              target_scale=scale_star)
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
        scale_star = branch["haircut_scale"]
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
        )
        _t_base = _time.perf_counter() - _t0

        # Refresh crisis-zone indicators on the new debt paths
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

    # Final zone consistency: the last base solve used the previous round's
    # indicator; if the refreshed indicator differs (debt crossed a zone
    # threshold in the last round), re-solve until path and indicator agree.
    # Interior to the crisis zone (the typical configuration) this never runs.
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


# ─────────────────────────────────────────────────────────────────────────────
# Diagnostics: exact per-period D-bond spread decomposition
# ─────────────────────────────────────────────────────────────────────────────

def bond_decomposition(out, ss, cal):
    """Decompose the D-bond promised excess return (annualized bps, deviation
    from SS) into default compensation + risk premium + liquidity premium.

    Exact identity per period:
      payoff^nd/Q − 1 − rdep  =  (payoff^nd − E[payoff])/Q      [default comp.]
                               + (E[payoff]/Q − 1 − rdep − λμ/Ω̃) [risk premium]
                               + λμ/Ω̃                            [liquidity]
    where E is the physical expectation over the default event.  In
    risk-neutral mode the risk premium is 0 by the pricing equation; with the
    Bocola channel it is positive because Ω^d > Ω^nd depresses Q.

    Returns dict of (T,) arrays in annualized bps (deviations from SS where
    applicable): total_yield, defcomp, risk, liquidity, plus raw arrays.
    """
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
