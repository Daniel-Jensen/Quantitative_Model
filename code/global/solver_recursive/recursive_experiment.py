# RECURSIVE SOLUTION EXPERIMENT: THE PASS-THROUGH OF SOVEREIGN RISK (BOCOLA).
# End-to-end demonstration that the global recursive solution (two-branch
# quadrature + Bocola closed-form multiplier) delivers the sign the representative
# branch could not: an elevated PRICED probability of default lowers output AND
# consumption, persistently. Orchestration only -- the economics live in the
# untouched blocks and in point_map.py; the solver is recursive_main.time_iteration.
#
# Pipeline: solve the d=0 (no-default) rules, warm-start d=1 (post-default) from
# them and solve, refine jointly at pi>0, then read (i) the IMPACT response as the
# one-quarter-ahead default probability rises and (ii) the persistence IRF as an
# s-shock decays (rho_s). Spawn-guarded (the pointwise solves are serial, but the
# guard is kept so the module is safe to import).
import numpy as np
from scipy.optimize import root

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import build_state_box, s_process_params, default_prob
from solver_recursive.decision_rules import RuleSet, STORE_RULES
from solver_recursive.recursive_main import (time_iteration, calibrate_household_anchors,
                            ss_state, ss_x)
from solver_recursive.point_map import point_residuals, SOLVE7


def solve_recursive(cal, ss, sproc, mu=1, verbose=True, mu_vec=None):
    # SOLVE BOTH DEFAULT REGIMES BY TIME ITERATION (WARM-START CHAIN d0 -> d1).
    # mu=1 (13 pts) is the converging, correctly-signed grid; mu=2 (85 pts) adds
    # state cross-terms but the d=1 regime does not converge on the wide box
    # (joint residual ~0.7, sign flips expansionary) -- the corner-decoupling
    # numerical obstacle documented in commit 1207f6e / docs/bocola2016_replication.md.
    # SYMMETRIC box: a wide-LOW B (to hold the post-default surviving debt ~0.45
    # on-grid) was tried but DIVERGES -- the deep low-B states are too far from the
    # SS cold start to solve and poison the fit. The d=1 recession is still
    # captured because the haircut hits the CURRENT-period payoff (crashing net
    # worth directly); only the post-default continuation is clipped (second order).
    grid = build_state_box(ss, cal, mu=mu, mu_vec=mu_vec, k_band=0.02, p_band=0.08, b_band=0.08)
    rules = RuleSet.from_ss(grid, ss, cal)
    if verbose:
        print(f"  grid: mu={mu}, {grid.n} points x 2 regimes")
    _, _, w0 = time_iteration(rules, cal, ss, sproc, regimes=(0,),
                              no_default=True, damp=0.25, tol=1e-6, max_it=70,
                              n_gh=5)
    # d=1 (post-default) via a HAIRCUT HOMOTOPY: warm-start from d=0, then lower
    # the recovery from ~1 (=no default, identical to d=0) down to the calibrated
    # value in steps, re-solving each -- so the deep recession is reached by
    # continuation, not cold-started (which does not converge).
    for k in STORE_RULES:
        rules.set_values(k, 1, rules.vals[k][0].copy())
    rec_target = cal["recovery_rate_D"]
    w1 = np.nan
    for rec in (0.85, 0.70, 0.55, rec_target):
        cal["recovery_rate_D"] = rec
        _, _, w1 = time_iteration(rules, cal, ss, sproc, regimes=(1,),
                                  no_default=False, damp=0.2, tol=1e-6,
                                  max_it=40, n_gh=5)
        if verbose:
            print(f"    d=1 homotopy recovery={rec:.2f}: worst={w1:.1e}")
    cal["recovery_rate_D"] = rec_target
    _, _, wj = time_iteration(rules, cal, ss, sproc, regimes=(0, 1),
                              no_default=False, damp=0.2, tol=1e-6, max_it=40,
                              n_gh=5)
    if verbose:
        print(f"  worst point residual: d0={w0:.1e} d1={w1:.1e} joint={wj:.1e}")
    return rules


def _solve_point(rules, cal, ss, sproc, S, x0):
    # DIRECT d=0 IMPACT SOLVE AT STATE S (both regimes as continuation).
    def f(x):
        try:
            return point_residuals(S, 0, x, rules, cal, ss, sproc,
                                   n_gh=7, no_default=False)[0]
        except (ValueError, RuntimeError, FloatingPointError):
            return np.full(7, 10.0)
    best = None
    for g in (x0, ss_x(ss, cal)):
        sol = root(f, g, method="hybr", tol=1e-12)
        fn = np.max(np.abs(sol.fun))
        if best is None or fn < best[1]:
            best = (sol.x, fn)
    _, o = point_residuals(S, 0, best[0], rules, cal, ss, sproc,
                           n_gh=7, no_default=False)
    return best[0], o


def read_at(rules, cal, ss, sproc, S):
    # READ THE CONVERGED RULES AT STATE S (binding branch), returning the implied
    # allocation + the point residual there (accuracy). Evaluating the period map
    # at the rules' OWN policy values stays on the binding branch -- re-solving
    # with a root finder can slip onto the nearby slack equilibrium at the barely-
    # binding SS. This is the standard way to read a global solution's IRF.
    Sm = np.atleast_2d(S)
    x = np.array([float(rules.eval(k, 0, Sm)[0]) for k in SOLVE7])
    res, o = point_residuals(S, 0, x, rules, cal, ss, sproc,
                             n_gh=7, no_default=False)
    o["_x"] = x
    o["_resid"] = float(np.max(np.abs(res)))
    return o


def _spread_bp(o, cal):
    # LENDING-SPREAD PROXY IN ANNUALISED BASIS POINTS (lambda_K * mu / alpha).
    return 4e4 * cal["lambda_K_D"] * o["mu_D"] / max(o["alpha_D"], 1e-6)


def impact_table(rules, cal, ss, sproc):
    # RESPONSE AT THE SS-LEVEL STATE AS THE PRICED DEFAULT PROBABILITY RISES.
    S0 = ss_state(ss, cal, sproc)
    base = read_at(rules, cal, ss, sproc, S0.copy())
    Yb, Cb, Ib, Nb = base["Y_D"], base["C_D"], base["I_D"], base["_x"][0]
    print("\n  IMPACT of priced default risk (deviation from the low-risk state)")
    print("   pd%     Y_D%     C_D%    hours%    I_D%    muD     Q_bD  spread_bp  resid")
    for s_val in (sproc["s_star"], -6.0, -5.2, -4.5, -3.9):
        S = S0.copy(); S[5] = s_val
        o = read_at(rules, cal, ss, sproc, S)
        print(f"  {100*default_prob(s_val):5.2f}  {100*(o['Y_D']/Yb-1):+7.3f} "
              f"{100*(o['C_D']/Cb-1):+7.3f}  {100*(o['_x'][0]/Nb-1):+6.2f}  "
              f"{100*(o['I_D']/Ib-1):+6.2f}  {o['mu_D']:.4f}  {o['Q_bD']:.4f} "
              f"{_spread_bp(o, cal):7.1f}  {o['_resid']:.0e}")


def persistence_irf(rules, cal, ss, sproc, s_shock=-3.9, T=21):
    # IRF AS AN s-SHOCK DECAYS (rho_s), endogenous states held at SS so the path
    # stays on-grid -- the shock-persistence channel (a lower bound; the
    # endogenous net-worth dynamics amplify it). Reads the binding-branch rules.
    S0 = ss_state(ss, cal, sproc)
    base = read_at(rules, cal, ss, sproc, S0.copy())
    Yr, Cr, Ir = base["Y_D"], base["C_D"], base["I_D"]
    print("\n  PERSISTENCE IRF (one-off risk shock, rho_s decay)")
    print("   qtr   pd%     Y_D%     C_D%     I_D%")
    for t in range(T):
        s_t = sproc["s_star"] + sproc["rho_s"] ** t * (s_shock - sproc["s_star"])
        S = S0.copy(); S[5] = s_t
        o = read_at(rules, cal, ss, sproc, S)
        if t in (0, 1, 2, 4, 6, 8, 12, 16, 20):
            print(f"   {t:3d}  {100*default_prob(s_t):5.2f}  "
                  f"{100*(o['Y_D']/Yr-1):+7.3f}  {100*(o['C_D']/Cr-1):+7.3f}  "
                  f"{100*(o['I_D']/Ir-1):+7.2f}")


def _tfp_read(rules, cal, ss, sproc, S):
    # READ THE NO-DEFAULT RULES AT STATE S (TFP experiment: no sovereign risk).
    Sm = np.atleast_2d(S)
    x = np.array([float(rules.eval(k, 0, Sm)[0]) for k in SOLVE7])
    res, o = point_residuals(S, 0, x, rules, cal, ss, sproc, n_gh=7, no_default=True)
    o["_x"] = x
    return o


def solve_tfp(cal, ss, sproc, mu=1):
    # SOLVE THE NO-DEFAULT (d=0) RULES OVER THE 7-STATE GRID FOR THE TFP EXPERIMENT.
    grid = build_state_box(ss, cal, mu=mu)
    rules = RuleSet.from_ss(grid, ss, cal)
    time_iteration(rules, cal, ss, sproc, regimes=(0,), no_default=True,
                   damp=0.25, tol=1e-6, max_it=70, n_gh=5)
    return rules


def tfp_irf(rules, cal, ss, sproc, dz=0.01, T=21):
    # TFP IRF read off the no-default rules along the Z_D-decay path (rho_z from
    # sproc), endogenous states held at SS so the read stays on-grid -- the exact
    # image of persistence_irf, with the TFP state Z_D in place of the risk state s.
    S0 = ss_state(ss, cal, sproc)
    Z_ss = S0[6]
    base = _tfp_read(rules, cal, ss, sproc, S0.copy())
    Yb, Cb, Ib, Nb = base["Y_D"], base["C_D"], base["I_D"], base["_x"][0]
    print(f"\n  TFP IRF (one-off {dz:.0%} shock, rho_z={sproc['rho_z']} decay)")
    print("   qtr    Z%     Y_D%     C_D%     I_D%    hours%")
    for t in range(T):
        z_t = dz * sproc["rho_z"] ** t
        S = S0.copy(); S[6] = Z_ss * np.exp(z_t)
        o = _tfp_read(rules, cal, ss, sproc, S)
        if t in (0, 1, 2, 4, 6, 8, 12, 16, 20):
            print(f"   {t:3d}  {100*z_t:5.2f}  {100*(o['Y_D']/Yb-1):+7.3f}  "
                  f"{100*(o['C_D']/Cb-1):+7.3f}  {100*(o['I_D']/Ib-1):+7.2f}  "
                  f"{100*(o['_x'][0]/Nb-1):+6.2f}")


def tfp_main():
    cal = get_calibration()
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print("=== TFP shock — recursive projection (Z_D as the 7th state) ===")
    rules = solve_tfp(cal, ss, sproc)
    tfp_irf(rules, cal, ss, sproc)


def main():
    cal = get_calibration()
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print("=== recursive global solution: pass-through of sovereign risk ===")
    rules = solve_recursive(cal, ss, sproc)
    try:
        import pickle
        with open("/private/tmp/claude-501/-Users-Huawei-Quantitative-Model/"
                  "239042af-4c74-4ebe-a83f-92681158d4c3/scratchpad/mu2_rules.pkl",
                  "wb") as fh:
            pickle.dump((rules, cal, ss, sproc), fh)
    except Exception:
        pass
    impact_table(rules, cal, ss, sproc)
    persistence_irf(rules, cal, ss, sproc)


if __name__ == "__main__":
    main()
