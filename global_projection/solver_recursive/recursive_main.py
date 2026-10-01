# RECURSIVE GLOBAL SOLUTION DRIVER (DELIVERABLE D): TIME ITERATION.
# Orchestrates the global recursive equilibrium of the frozen two-country model:
#   1. build the Smolyak state box + cold-start rules from the SS;
#   2. calibrate the rep-agent household anchors so the SS is an exact rest point;
#   3. TIME ITERATION: each sweep freezes the previous iterate as the continuation
#      and solves the SEVEN market-clearing unknowns pointwise (the per-period
#      image of transition.py). The banker valuations alpha/Q_b AND the household
#      aggregates C/A are read off inside the point map (the rep-agent deposit
#      Euler closure -- Tier-1 KS ladder), stored, damped and refit;
#   4. simulate from the ergodic mean and read the s-shock IRFs (deliverable E).
# Shock/experiment definition lives here; console output goes through prints.
# The economic blocks are untouched. Pointwise solves are serial scipy-root, so
# no multiprocessing spawn hazard; the __main__ guard is kept regardless.
import numpy as np
from scipy.optimize import root

from solver_recursive.decision_rules import RuleSet, SOLVE7, DERIVED
from solver_recursive.point_map import point_residuals
from solver_recursive.state_grid import build_state_box, s_process_params, NSTATE


def ss_state(ss, cal, sproc):
    # THE SS POINT IN THE 10-STATE VECTOR
    # [K_D, K_F, P_D, P_F, b_DD, b_DF, b_FD, V_dep, s, Z_D]. b_DD/b_DF are the
    # two banks' carried holdings of the D sovereign and b_FD is the D bank's carried
    # holding of the F sovereign (both splits endogenous); V_dep is the carried
    # cross-border deposit position W_D - P_D, the margin national clearing suppressed.
    # It is ZERO at the symmetric SS, where each household's claim equals its own bank's
    # obligation. The CB backstop carries no state, so the SS point is the same vector
    # whether or not the facility is on -- which is what makes phi = 0 nest exactly.
    b_DF = cal["b_D_F_ss"]
    return np.array([ss["Kap_D_ss"], ss["Kap_F_ss"],
                     ss["ss_bank_D"]["P_state_ss"], ss["ss_bank_F"]["P_state_ss"],
                     cal["B_gov_D_ss"] - b_DF, b_DF, cal["b_F_D_ss"], 0.0,
                     sproc["s_star"], cal["Z_ss_D"]])


def ss_x(ss, cal):
    # THE SS VALUES OF THE THIRTEEN UNKNOWNS, in decision_rules.SOLVE order.
    return np.array([1.0, 1.0, ss["Kap_D_ss"], ss["Kap_F_ss"],
                     cal["r_dep_D_target"], cal["r_dep_F_target"], ss["p_ss"],
                     ss["Q_bD_ss"], cal["b_D_F_ss"],
                     ss["Q_bF_ss"], cal["b_F_D_ss"],
                     ss["A_D_ss"], ss["A_F_ss"]])


def calibrate_household_anchors(cal, ss, sproc, tol=1e-13, max_it=12):
    # SET ss["hh_T_D/F"] SO THE BUDGET GIVES C = C_ss AT THE SS (with A = dep/P_CES
    # = A_ss, deposits clearing by quantity).
    # C is linear in hh_T, so ONE zero-anchor evaluation would pin it -- except when
    # that trial evaluation lands on a GUARD. At hh_T = 0 the household is short the
    # whole working-capital repayment (1+r_wc)*L_wc, which puts C below the 0.3*C_ss
    # floor in _sclip; anchoring off the clipped value then leaves a permanent gap of
    # exactly that size and the SS stops being a rest point. Iterating to a fixed point
    # is self-correcting whatever guards are active, and costs a handful of evaluations.
    ss["hh_T_D"] = 0.0
    ss["hh_T_F"] = 0.0
    grid = build_state_box(ss, cal, mu=1)
    rules = RuleSet.from_ss(grid, ss, cal)
    S0, x0 = ss_state(ss, cal, sproc), ss_x(ss, cal)
    for _ in range(max_it):
        _, out = point_residuals(S0, 0, x0, rules, cal, ss, sproc, no_default=True)
        gap_D = ss["C_D_ss"] - out["C_D"]
        gap_F = ss["C_F_ss"] - out["C_F"]
        ss["hh_T_D"] += gap_D
        ss["hh_T_F"] += gap_F
        if max(abs(gap_D), abs(gap_F)) < tol:
            break
    return ss["hh_T_D"], ss["hh_T_F"]


def solve_point(S, d, cont, cal, ss, sproc, x0, no_default=False, n_gh=7,
                x_ss=None, no_cb=False):
    # SOLVE THE MARKET-CLEARING UNKNOWNS AT ONE POINT (FROZEN CONTINUATION).
    # hybr from the warm start (the common case: one cheap solve near the fixed
    # point); a single fallback from the SS guess only if that misses.
    def f(x):
        try:
            return point_residuals(S, d, x, cont, cal, ss, sproc, n_gh=n_gh,
                                   no_default=no_default, no_cb=no_cb)[0]
        except (ValueError, RuntimeError, FloatingPointError):
            return np.full(len(SOLVE7), 10.0)

    sol = root(f, x0, method="hybr", tol=1e-12)
    best = (sol.x, np.max(np.abs(sol.fun)))
    if best[1] > 1e-9 and x_ss is not None:
        sol2 = root(f, x_ss, method="hybr", tol=1e-12)
        if np.max(np.abs(sol2.fun)) < best[1]:
            best = (sol2.x, np.max(np.abs(sol2.fun)))
    # evaluate at the best root, falling back to x0; if BOTH raise (a genuinely
    # infeasible point) return a failure sentinel (fn=1e3) rather than crash -- the
    # sweep then retains/masks it. Needed for off-box warm starts (e.g. EDS points).
    for xt in (best[0], x0):
        try:
            _, out = point_residuals(S, d, xt, cont, cal, ss, sproc, n_gh=n_gh,
                                     no_default=no_default, no_cb=no_cb)
            return xt, out, best[1]
        except (ValueError, RuntimeError, FloatingPointError):
            continue
    return x0, None, 1e3


def _sweep(rules, cont, cal, ss, sproc, regimes, no_default, n_gh,
           keep_tol=1e-3, no_cb=False):
    # ONE TIME-ITERATION SWEEP: SOLVE EVERY POINT, RETURN NEW RULE VALUE ARRAYS.
    # A point whose solve does not clear (fn > keep_tol) RETAINS the previous
    # iterate's values -- a failed corner must never poison the continuation.
    n = rules.grid.n
    new = {k: {d: np.empty(n) for d in regimes} for k in STORE()}
    wt = {d: np.ones(n) for d in regimes}     # per-point fit weight (0 = failed corner)
    x_ss = ss_x(ss, cal)     # was a duplicated literal; it silently kept 7 entries
    worst, n_fail = 0.0, 0
    for d in regimes:
        for i in range(n):
            S = rules.grid.points[i]
            x0 = np.array([rules.vals[k][d][i] for k in SOLVE7])
            x, out, fn = solve_point(S, d, cont, cal, ss, sproc, x0,
                                     no_default=no_default, n_gh=n_gh, x_ss=x_ss,
                                     no_cb=no_cb)
            worst = max(worst, fn if np.isfinite(fn) else 1e3)
            if (not np.isfinite(fn)) or fn > keep_tol:   # retain old values, mask fit
                n_fail += 1
                wt[d][i] = 0.0
                for k in STORE():
                    new[k][d][i] = rules.vals[k][d][i]
            else:
                for j, k in enumerate(SOLVE7):
                    new[k][d][i] = x[j]
                for k in DERIVED:
                    new[k][d][i] = out[k]
    return new, worst, n_fail, wt


def STORE():
    # THE 15 STORED RULE NAMES (SOLVE7 + DERIVED), one place.
    return SOLVE7 + DERIVED


def p_block_rotation(ss, cal, sproc, eps=1e-3, mu=1, mu_vec=None, probe_kw=None):
    # EIGENBASIS OF THE (P_D, P_F) TRANSITION JACOBIAN AT THE SS -- BOCOLA'S V.
    # Returns (rot, centre, J, |eig J|) for build_state_box: rot maps NATURAL states
    # to the coordinates in which the deposit-obligation block is diagonal.
    # The probe grid MUST carry the same mu/mu_vec and bands as the grid the rotation
    # is for: J depends on which states are live. On the isotropic mu=1 box capital is
    # a real state, the continuation absorbs the capital response and rho(|J|) = 0.92;
    # on the production anisotropic box K_D/K_F have a single node, the rules are flat
    # in capital, that channel is switched off and rho(|J|) = 1.96. Probing the wrong
    # grid therefore designs a box for dynamics the solver will not have.
    # WHY: the measured Jacobian d(P_D',P_F')/d(P_D,P_F) at the SS is
    #   [[ 1.076, -0.310], [ 2.016, -1.252]]
    # -- eigenvalues 0.766 and -0.942, so the dynamics are stable, but the matrix is
    # strongly non-normal (singular values 2.61 / 0.28). Its entrywise absolute value
    # has spectral radius 1.96, and |J| b <= b then has NO positive solution b, i.e.
    # NO axis-aligned box centred on the SS is one-step invariant, at any bandwidth.
    # On the eigenbasis the map is diagonal with |lambda| < 1, so every box IS.
    grid = build_state_box(ss, cal, mu=mu, mu_vec=mu_vec, **(probe_kw or {}))
    rules = RuleSet.from_ss(grid, ss, cal)
    S0, x0 = ss_state(ss, cal, sproc), ss_x(ss, cal)
    _, o0, _ = solve_point(S0, 0, rules, cal, ss, sproc, x0, no_default=True,
                           n_gh=5, x_ss=x0)
    base = np.array([o0["Pp_D"], o0["Pp_F"]])
    J = np.zeros((2, 2))
    for j, idx in enumerate((2, 3)):
        S = S0.copy(); S[idx] *= (1.0 + eps)
        _, o, _ = solve_point(S, 0, rules, cal, ss, sproc, x0, no_default=True,
                              n_gh=5, x_ss=x0)
        J[:, j] = (np.array([o["Pp_D"], o["Pp_F"]]) - base) / (S0[idx] * eps)
    w, V = np.linalg.eig(J)
    if np.iscomplexobj(V):                 # complex pair -> use the real Schur basis
        V = np.linalg.qr(np.column_stack([V.real[:, 0], V.imag[:, 0]]))[0]
    V = np.real(V)
    rot = np.eye(NSTATE)
    rot[np.ix_((2, 3), (2, 3))] = np.linalg.inv(V)
    return rot, S0.copy(), J, np.abs(w)


def time_iteration(rules, cal, ss, sproc, regimes=(0, 1),
                   no_default=False, damp=0.5, tol=1e-7, max_it=60,
                   n_gh=7, verbose=False, no_cb=False):
    # TIME-ITERATE THE RULES TO A FIXED POINT (IN PLACE).
    # Returns (converged, iters, worst_point_residual, n_fail). n_fail is part of the
    # contract because the exit test needs BOTH a settled rule AND every point
    # clearing: a sweep can look converged on `change` while a quarter of the grid is
    # frozen on its previous values, and the caller must be able to see that.
    # cal["fit_mask"]/cal["fit_ridge"] (optional) switch the coefficient fit from the
    # exact square collocation to a ROBUST weighted ridge-LS fit: failed corners get a
    # small weight fit_fail_weight (default 0.1) rather than 0 -- they still anchor the
    # fit (keeping it full-rank and sane) but their reason-2 poisoning is cut ~10x, and
    # the ridge damps the high-degree wiggle. HARD masking (weight 0) collapses the fit
    # when a sweep has many failures. Absent -> exact fit, behaviour unchanged.
    fit_mask = bool(cal.get("fit_mask", False))
    fit_ridge = float(cal.get("fit_ridge", 0.0))
    fit_fw = float(cal.get("fit_fail_weight", 0.1))
    rules.n_gh = int(n_gh)                        # stamp: readers must match the solve
    worst, n_fail = np.inf, len(regimes) * rules.grid.n
    for it in range(max_it):
        cont = rules.copy()                       # frozen continuation
        new, worst, n_fail, wt = _sweep(rules, cont, cal, ss, sproc, regimes,
                                        no_default, n_gh, no_cb=no_cb)
        change = 0.0
        for k in STORE():
            for d in regimes:
                old = rules.vals[k][d]
                upd = damp * new[k][d] + (1.0 - damp) * old
                change = max(change, np.max(np.abs(upd - old))
                             / (np.max(np.abs(old)) + 1e-8))
                weights = (np.where(wt[d] > 0.5, 1.0, fit_fw) if fit_mask else None)
                rules.set_values(k, d, upd, weights=weights, ridge=fit_ridge)
        if verbose:
            print(f"    [time-it {it + 1:2d}] max|F_point|={worst:.2e}  "
                  f"rel rule change={change:.2e}  fails={n_fail}/"
                  f"{len(regimes) * rules.grid.n}")
        if change < tol and worst < 1e-6:
            return True, it + 1, worst, n_fail
    return False, max_it, worst, n_fail
