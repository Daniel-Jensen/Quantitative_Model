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

from decision_rules import RuleSet, SOLVE7, DERIVED
from point_map import point_residuals
from state_grid import build_state_box, s_process_params


def ss_state(ss, cal, sproc):
    # THE STEADY-STATE POINT IN THE 6-STATE VECTOR [K_D,K_F,P_D,P_F,B_D,s].
    return np.array([ss["Kap_D_ss"], ss["Kap_F_ss"],
                     (1 + cal["r_dep_D_target"]) * ss["ss_bank_D"]["Dep_supply_ss"],
                     (1 + cal["r_dep_F_target"]) * ss["ss_bank_F"]["Dep_supply_ss"],
                     cal["B_gov_D_ss"], sproc["s_star"]])


def ss_x(ss, cal):
    # THE SS VALUES OF THE SEVEN UNKNOWNS.
    return np.array([1.0, 1.0, ss["Kap_D_ss"], ss["Kap_F_ss"],
                     cal["r_dep_D_target"], cal["r_dep_F_target"], ss["p_ss"]])


def calibrate_household_anchors(cal, ss, sproc):
    # SET ss["hh_T_D/F"] SO THE BUDGET GIVES C = C_ss AT THE SS (with A = dep/P_CES
    # = A_ss, deposits clearing by quantity). C depends linearly on hh_T, so one
    # zero-anchor evaluation pins it.
    ss["hh_T_D"] = 0.0
    ss["hh_T_F"] = 0.0
    grid = build_state_box(ss, cal, mu=1)
    rules = RuleSet.from_ss(grid, ss, cal)
    _, out = point_residuals(ss_state(ss, cal, sproc), 0, ss_x(ss, cal),
                             rules, cal, ss, sproc, no_default=True)
    ss["hh_T_D"] = ss["C_D_ss"] - out["C_D"]
    ss["hh_T_F"] = ss["C_F_ss"] - out["C_F"]
    return ss["hh_T_D"], ss["hh_T_F"]


def solve_point(S, d, cont, cal, ss, sproc, x0, no_default=False, n_gh=7,
                x_ss=None):
    # SOLVE THE 7 MARKET-CLEARING UNKNOWNS AT ONE POINT (FROZEN CONTINUATION).
    # hybr from the warm start (the common case: one cheap solve near the fixed
    # point); a single fallback from the SS guess only if that misses.
    def f(x):
        try:
            return point_residuals(S, d, x, cont, cal, ss, sproc,
                                   n_gh=n_gh, no_default=no_default)[0]
        except (ValueError, RuntimeError, FloatingPointError):
            return np.full(7, 10.0)

    sol = root(f, x0, method="hybr", tol=1e-12)
    best = (sol.x, np.max(np.abs(sol.fun)))
    if best[1] > 1e-9 and x_ss is not None:
        sol2 = root(f, x_ss, method="hybr", tol=1e-12)
        if np.max(np.abs(sol2.fun)) < best[1]:
            best = (sol2.x, np.max(np.abs(sol2.fun)))
    x = best[0]
    try:
        _, out = point_residuals(S, d, x, cont, cal, ss, sproc,
                                 n_gh=n_gh, no_default=no_default)
    except (ValueError, RuntimeError, FloatingPointError):
        x = x0
        _, out = point_residuals(S, d, x, cont, cal, ss, sproc,
                                 n_gh=n_gh, no_default=no_default)
    return x, out, best[1]


def _sweep(rules, cont, cal, ss, sproc, regimes, no_default, n_gh,
           keep_tol=1e-3):
    # ONE TIME-ITERATION SWEEP: SOLVE EVERY POINT, RETURN NEW RULE VALUE ARRAYS.
    # A point whose solve does not clear (fn > keep_tol) RETAINS the previous
    # iterate's values -- a failed corner must never poison the continuation.
    n = rules.grid.n
    new = {k: {d: np.empty(n) for d in regimes} for k in STORE()}
    x_ss = np.array([1.0, 1.0, ss["Kap_D_ss"], ss["Kap_F_ss"],
                     cal["r_dep_D_target"], cal["r_dep_F_target"], ss["p_ss"]])
    worst, n_fail = 0.0, 0
    for d in regimes:
        for i in range(n):
            S = rules.grid.points[i]
            x0 = np.array([rules.vals[k][d][i] for k in SOLVE7])
            x, out, fn = solve_point(S, d, cont, cal, ss, sproc, x0,
                                     no_default=no_default, n_gh=n_gh, x_ss=x_ss)
            worst = max(worst, fn if np.isfinite(fn) else 1e3)
            if (not np.isfinite(fn)) or fn > keep_tol:   # retain old values
                n_fail += 1
                for k in STORE():
                    new[k][d][i] = rules.vals[k][d][i]
            else:
                for j, k in enumerate(SOLVE7):
                    new[k][d][i] = x[j]
                for k in DERIVED:
                    new[k][d][i] = out[k]
    return new, worst, n_fail


def STORE():
    # THE 15 STORED RULE NAMES (SOLVE7 + DERIVED), one place.
    return SOLVE7 + DERIVED


def time_iteration(rules, cal, ss, sproc, regimes=(0, 1),
                   no_default=False, damp=0.5, tol=1e-7, max_it=60,
                   n_gh=7, verbose=False):
    # TIME-ITERATE THE RULES TO A FIXED POINT (IN PLACE).
    worst = np.inf
    for it in range(max_it):
        cont = rules.copy()                       # frozen continuation
        new, worst, n_fail = _sweep(rules, cont, cal, ss, sproc, regimes,
                                    no_default, n_gh)
        change = 0.0
        for k in STORE():
            for d in regimes:
                old = rules.vals[k][d]
                upd = damp * new[k][d] + (1.0 - damp) * old
                change = max(change, np.max(np.abs(upd - old))
                             / (np.max(np.abs(old)) + 1e-8))
                rules.set_values(k, d, upd)
        if verbose:
            print(f"    [time-it {it + 1:2d}] max|F_point|={worst:.2e}  "
                  f"rel rule change={change:.2e}  fails={n_fail}/"
                  f"{len(regimes) * rules.grid.n}")
        if change < tol and worst < 1e-6:
            return True, it + 1, worst
    return False, max_it, worst
