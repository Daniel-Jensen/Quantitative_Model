# M2, WARM-STARTED (2026-10-28): THE REST-POINT CALIBRATION AT APPENDIX B. calibrate_srp.py COLD-STARTS EVERY
# EVALUATION, WHICH FAILS HERE; THIS DRIVER SEEDS EACH EVALUATION FROM THE LAST ACCEPTED SOLUTION (SS-RELATIVE),
# SOLVES THE HARD SYSTEM (GENUINE IC MAX, HARD BOOKS, HARD FLOOR) WITH THE CALM GUARD, AND FALLS BACK THROUGH THE
# WORKING SMOOTHINGS. PHASE A RAISES THE D SPREAD INSTRUMENT UNTIL THE D BANK BINDS AT REST (THE MOMENT IS FLAT AT
# 0 WHILE SLACK, SO BROYDEN CANNOT START THERE); PHASE B IS calibrate_stochastic's BROYDEN ON ALL SIX.
#   python3 -u m2_warm.py <hard checkpoint>   (e.g. m1_calm_hard_mu_eps_book_eps_floor_eps.pkl)
import os, pickle, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np

HARD = dict(mu_eps=0.0, book_eps=0.0, floor_eps=0.0)
CUR = {}
# the grid refinement the evaluations solve on (None: the coarse grid; (IS, 9, "cross"): Bocola's fine grid)
REFINE = None


def seed_rel(prev, cal_p, ss_p, cal_k, ss_k):
    # THE NEW ECONOMY'S RULES ON THE SAME GRID DESIGN: ITS SS PLUS THE LAST SOLUTION'S DEVIATION, NODE BY NODE.
    from global_projection.solver_recursive.decision_rules import RuleSet, STORE_RULES, to_fit, from_fit
    from global_projection.solver_recursive import recursive_experiment as RE
    ss_old = RuleSet.from_ss(prev.grid, ss_p, cal_p, n_regimes=2)
    new = RuleSet.from_ss(RE.normal_box(ss_k, cal_k, mu=1, refine=REFINE), ss_k, cal_k, n_regimes=2)
    for k in STORE_RULES:
        for d in range(new.n_regimes):
            dev = to_fit(k, prev.vals[k][d]) - to_fit(k, ss_old.vals[k][d])
            new.set_values(k, d, from_fit(k, to_fit(k, new.vals[k][d]) + dev))
    new.n_gh = prev.n_gh
    return new


def solve_at(RE, x, verbose=False):
    # THE HARD SYSTEM AT INSTRUMENTS x, SEEDED FROM CUR; CUR MOVES ONLY ON SUCCESS. RETURNS (ok, max|F|).
    import global_projection.solver_recursive.calibrate_stochastic as CS
    from diag_homotopy import economy
    from m1_coarse import stage, calm_guard
    relative_seed = seed_rel
    v = dict(CUR["vals"], **HARD, **{p[0]: float(xi) for p, xi in zip(CS.PAIRS, x)})
    routes = (("hard", v), ("working", dict(v, **RE.HARD_WORK)))
    rules, cal, ss = CUR["rules"], CUR["cal"], CUR["ss"]
    for name, vv in routes:
        cal_k, ss_k, sp_k = economy(CUR["base"], vv)
        trial = relative_seed(rules, cal, ss, cal_k, ss_k)
        ok, w = stage(RE, trial, cal_k, ss_k, sp_k, f"m2 {name}", maxit=30)
        g = calm_guard(RE, trial, cal_k, ss_k, sp_k, f"m2 {name}", maxit=30)
        if g is not None:
            ok, w = g
        if ok and name == "working":
            # from the working smoothings, the genuine conditions polished in
            cal_h, ss_h, sp_h = economy(CUR["base"], v)
            hard = relative_seed(trial, cal_k, ss_k, cal_h, ss_h)
            ok, w = stage(RE, hard, cal_h, ss_h, sp_h, "m2 hard polish", maxit=30)
            trial, cal_k, ss_k, sp_k = hard, cal_h, ss_h, sp_h
        if ok:
            CUR.update(rules=trial, cal=cal_k, ss=ss_k, sp=sp_k, vals=v)
            return True, w
    return False, w


def evaluate(RE, x):
    # ONE WARM EVALUATION: THE REST-POINT MOMENTS AND THE VERDICT, IN calibrate_stochastic's FORMAT.
    import global_projection.solver_recursive.calibrate_stochastic as CS
    from global_projection.solver_recursive.state_grid import STATE_NAMES
    t0 = time.perf_counter()
    ok, w = solve_at(RE, x)
    if not ok:
        return np.full(len(CS.PAIRS), np.nan), dict(solve_ok=False, solve_worst=float(w))
    S, o = CS.rest_point_read(CUR["rules"], CUR["cal"], CUR["ss"], CUR["sp"])
    m = CS.rest_point_moments(o)
    u = CUR["rules"].grid.to_unit(np.atleast_2d(S))[0]
    out = [k for k, ui in zip(STATE_NAMES, u) if abs(ui) > 1.0]
    print(f"  [m2] solved {w:.1e} in {time.perf_counter() - t0:.0f}s; rest point mu_D {o['mu_D']:.5f} mu_F "
          f"{o['mu_F']:.5f}; outside the box: {out or 'none'}", flush=True)
    with open(HERE / "runs" / "m2w_current.pkl", "wb") as fh:
        pickle.dump({k: CUR[k] for k in ("rules", "cal", "ss", "vals")}, fh)
    return (np.array([m[p[2]] for p in CS.PAIRS]),
            dict(solve_ok=True, solve_worst=float(w), read_resid=float(o["_resid"]),
                 p_ss=float(CUR["ss"]["p_ss"]), p_rest=float(o["p"])))


def bind_phase(RE, x, tries=8):
    # RAISE THE D SPREAD INSTRUMENT (x2 PER RUNG, HALVED ON A FAILED SOLVE) UNTIL THE D BANK BINDS AT REST.
    import global_projection.solver_recursive.calibrate_stochastic as CS
    i = [p[0] for p in CS.PAIRS].index("credit_spread_target_D")
    target, hi = CUR["cal"]["srp_spread_D"], CS.PAIRS[i][5]
    good, fac = x.copy(), 2.0
    for _ in range(tries):
        x = good.copy()
        x[i] = min(good[i] * fac, hi)
        print(f"\n=== bind phase: credit_spread_target_D {4e4 * x[i]:.1f} bp/yr", flush=True)
        m, info = evaluate(RE, x)
        if not info["solve_ok"]:
            fac = 1.0 + 0.5 * (fac - 1.0)
            print(f"  [bind] solve failed; factor {fac:.3f}", flush=True)
            continue
        good = x
        print(f"  [bind] rest-point D spread {4e4 * m[i]:.2f} bp (target {4e4 * target:.1f})", flush=True)
        if m[i] >= 0.5 * target:
            return good
    return good


def main():
    from calibration.global_projection import get_calibration
    import global_projection.solver_recursive.calibrate_stochastic as CS
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import s_process_params
    from diag_path import FLOOR
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    base = get_calibration()
    base.update(nw_floor_frac=FLOOR, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3, jacobian="complex")
    CUR.update(rules=sv["rules"], cal=sv["cal"], ss=sv["ss"], sp=s_process_params(sv["cal"]),
               vals=dict(sv["vals"], **HARD), base=base)
    x = np.array([CUR["cal"][p[0]] for p in CS.PAIRS])
    state = HERE / "runs" / "m2w_state.json"
    if not state.exists():
        x = bind_phase(RE, x)
    # phase B: calibrate_stochastic's Broyden, its evaluation replaced by the warm one
    CS.get_calibration = lambda: dict(CUR["cal"], **{p[0]: float(xi) for p, xi in zip(CS.PAIRS, x)})
    CS.evaluate = lambda xx, settings, verbose=True: evaluate(RE, xx)
    xs, m = CS.calibrate(dict(), state, max_evals=16)
    print(f"\n[m2] {'CONVERGED' if m is not None else 'not converged'}: "
          + "  ".join(f"{p[0]}={v:.8g}" for p, v in zip(CS.PAIRS, xs)), flush=True)


if __name__ == "__main__":
    main()
