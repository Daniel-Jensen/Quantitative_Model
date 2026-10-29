# M2, OCCASIONALLY BINDING (2026-10-28 DECISION): THE D BANK'S IC IS SLACK AT THE REST POINT AND f_D IS CALIBRATED SO
# IT STARTS TO BIND AT THE CRISIS THRESHOLD p^d* (OB_PD, PER QUARTER, THE OTHER STATES AT REST). THE s BOX IS EXTENDED
# UPWARD SO THE THRESHOLD IS INSIDE IT; THE OTHER FIVE REST-POINT TARGETS (LEVERAGE D/F, THE F SPREAD, EXPOSURES D/F)
# ARE MATCHED BY calibrate_stochastic's BROYDEN WITH f_D HELD; THE D SPREAD INSTRUMENT STAYS AT ITS SS VALUE.
#   python3 -u m2_ob.py <hard checkpoint> [--fine] [--f F_D] [--broyden]
import json, os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np

OB_PD = 0.04
OB_TOL = 0.002
# the alternative f_D target (2026-10-29 decision): the lending-spread rise at p^d 2% from the rest point, appendix B's
# GR-DE overdraft-spread row (108 bp, an upper bound), exact reads
SPREAD_PD, SPREAD_BP, SPREAD_TOL = 0.02, 108.0, 10.0
# the largest move in f_D per warm fine solve
F_STEP = 0.02
# the extended box's top, p^d per quarter, walked in rungs (the 2010-11 average is 9.1%)
TOP_RUNGS = (0.06, 0.12)
LOGIT = lambda p: float(np.log(p / (1.0 - p)))


def slack(RE, rules, cal, ss, sp, S):
    # THE D SLACK RATIO E[Om](1+r)n/(divertable assets) AT STATE S, EXACT READ (1 = BINDING).
    from measure import omega_legs
    b = omega_legs(rules, cal, ss, sp, S, RE.read_exact)
    return b["E_Om"] * (1.0 + b["o"]["rdep_D"]) * b["n_D"] / b["lev_D"]


def threshold(RE, M):
    # THE p^d AT WHICH THE D IC STARTS TO BIND, THE OTHER STATES AT THE REST POINT; ABOVE THE BOX, THE LINEAR
    # EXTRAPOLATION OF THE SLACK IN p^d FROM THE TOP TWO POINTS (FLAGGED).
    from global_projection.solver_recursive.state_grid import IS
    r, cal, ss, sp = M.CUR["rules"], M.CUR["cal"], M.CUR["ss"], M.CUR["sp"]
    S0 = RE.stochastic_rest_point(r, cal, ss, sp, verbose=False)
    top = 1.0 / (1.0 + np.exp(-r.grid.hi[IS]))
    pd_of = lambda s: 1.0 / (1.0 + np.exp(-s))

    def g(pd):
        S = S0.copy(); S[IS] = LOGIT(pd)
        return slack(RE, r, cal, ss, sp, S) - 1.0
    lo, hi = pd_of(S0[IS]), 0.98 * top
    g_lo, g_hi = g(lo), g(hi)
    if g_lo <= 0.0:
        return lo, "binding at rest"
    if g_hi > 0.0:
        g_m = g(0.8 * hi)
        slope = (g_hi - g_m) / (hi - 0.8 * hi)
        return (hi - g_hi / slope if slope < 0 else np.inf), f"above the box top {100 * top:.1f}% (extrapolated)"
    for _ in range(30):
        mid = np.sqrt(lo * hi)
        if g(mid) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi / lo - 1.0 < 1e-3:
            break
    return np.sqrt(lo * hi), "inside the box"


def spread_rise(RE, M, pd=SPREAD_PD):
    # THE LENDING-SPREAD RISE r_wc - r_dep (bp/yr) AT p^d = pd FROM THE REST POINT, THE OTHER STATES THERE (EXACT).
    from global_projection.solver_recursive.state_grid import IS
    r, cal, ss, sp = M.CUR["rules"], M.CUR["cal"], M.CUR["ss"], M.CUR["sp"]
    S0 = RE.stochastic_rest_point(r, cal, ss, sp, verbose=False)
    S = S0.copy(); S[IS] = LOGIT(pd)
    o0, o1 = RE.read_exact(r, cal, ss, sp, S0), RE.read_exact(r, cal, ss, sp, S)
    return 4e4 * ((o1["r_wc_D"] - o1["rdep_D"]) - (o0["r_wc_D"] - o0["rdep_D"]))


def extend_box(RE, M, top):
    # THE s BOX's TOP RAISED TO p^d = top: A NEW GRID SEEDED BY EVALUATING THE CURRENT RULES THERE, SOLVED HARD.
    from global_projection.solver_recursive.decision_rules import RuleSet
    from m1_coarse import stage, calm_guard
    RE.BOX_KW = dict(RE.BOX_KW, s_hi=LOGIT(top))
    cal, ss, sp = M.CUR["cal"], M.CUR["ss"], M.CUR["sp"]
    new = RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1, refine=M.REFINE), 2), M.CUR["rules"])
    ok, w = stage(RE, new, cal, ss, sp, f"box top {100 * top:.0f}%", maxit=40)
    g = calm_guard(RE, new, cal, ss, sp, f"box top {100 * top:.0f}%", maxit=40)
    if g is not None:
        ok, w = g
    if ok:
        M.CUR["rules"] = new
    return ok, w


def main():
    from calibration.global_projection import get_calibration
    import global_projection.solver_recursive.calibrate_stochastic as CS
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import s_process_params
    from diag_path import FLOOR
    import m2_warm as M
    from global_projection.solver_recursive.state_grid import IS
    import m1_coarse
    args = sys.argv[2:]
    fine = "--fine" in args
    by_spread = "--target" in args and args[args.index("--target") + 1] == "spread2"
    if fine:
        M.REFINE = (IS, 9, "cross")
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    base = get_calibration()
    base.update(nw_floor_frac=FLOOR, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3, jacobian="complex")
    M.CUR.update(rules=sv["rules"], cal=sv["cal"], ss=sv["ss"], sp=s_process_params(sv["cal"]),
                 vals=dict(sv["vals"], **M.HARD), base=base)
    if "--f" in args:
        M.CUR["vals"]["f_D"] = float(args[args.index("--f") + 1])
    tag0 = "fine" if fine else "coarse"
    log = HERE / "runs" / f"m2_ob_{tag0}_state.json"
    hist = []

    def save(tag):
        # THE CURRENT SOLUTION AND THE HISTORY.
        with open(HERE / "runs" / f"m2_ob_{tag0}_{tag}.pkl", "wb") as fh:
            pickle.dump({k: M.CUR[k] for k in ("rules", "cal", "ss", "vals")}, fh)
        log.write_text(json.dumps(hist, indent=1))

    # 1. the box extended upward, rung by rung, at the starting f_D
    spec = args[args.index("--rungs") + 1] if "--rungs" in args else None
    rungs = TOP_RUNGS if spec is None else (() if spec == "none" else tuple(float(v) for v in spec.split(",")))
    if spec == "none":
        # resuming on an extended box: its top, for the seeds of later evaluations
        from global_projection.solver_recursive.state_grid import IS as _IS
        RE.BOX_KW = dict(RE.BOX_KW, s_hi=float(M.CUR["rules"].grid.hi[_IS]))
    for top in rungs:
        ok, w = extend_box(RE, M, top)
        print(f"[box] top p^d {100 * top:.0f}%: ok {ok}  max|F| {w:.2e}", flush=True)
        if ok:
            save(f"box{100 * top:.0f}")
        if not ok:
            print("[box] the extension failed; stopping", flush=True)
            return
    save("box")
    # 2. f_D by the secant method on the threshold, the other instruments held
    x = np.array([M.CUR["cal"][p[0]] for p in CS.PAIRS])
    pts = []
    f = float(M.CUR["vals"].get("f_D", M.CUR["cal"]["f_D"]))
    for k in range(8):
        M.CUR["vals"]["f_D"] = f
        ok, w = M.solve_at(RE, x)
        if not ok:
            print(f"[f] solve failed at f_D {f:.4f}", flush=True)
            if not pts:
                return
            f = 0.5 * (f + pts[-1][0])
            continue
        p, where = threshold(RE, M)
        sp2 = spread_rise(RE, M)
        # the moment the secant matches: the 2% spread rise (--target spread2) or the binding threshold
        mom, tgt, tol = (sp2, SPREAD_BP, SPREAD_TOL) if by_spread else (p, OB_PD, OB_TOL)
        pts.append((f, mom))
        hist.append(dict(stage="f", f_D=f, threshold=p, where=where, spread_rise_2pct=sp2, solve=w))
        print(f"[f] f_D {f:.4f}: IC binds from p^d {100 * p:.2f}% ({where}); spread rise at p^d 2% {sp2:.1f} bp "
              f"(target {'%.0f bp' % SPREAD_BP if by_spread else '%.1f%% threshold' % (100 * OB_PD)})", flush=True)
        save("f")
        if abs(mom - tgt) <= tol:
            break
        if "--search" not in args:
            break
        if len(pts) == 1:
            # measured on the fine 12%-top box (2026-10-29): a higher f lowers the 2% spread rise (f 0.06: 547 bp,
            # 0.08: 452 bp) and moves the threshold up
            f = f + F_STEP if mom > tgt else f - F_STEP
        else:
            (f1, p1), (f2, p2) = pts[-2], pts[-1]
            f_new = f2 + (tgt - p2) * (f2 - f1) / (p2 - p1) if np.isfinite(p1 + p2) and p2 != p1 else f2 + F_STEP
            # a trust region: one warm fine solve does not survive a larger move in f (0.06 -> 0.152 stalled)
            f = f2 + float(np.clip(f_new - f2, -F_STEP, F_STEP))
        f = float(np.clip(f, 0.02, 0.5))
    # the impact rows at crisis-level shocks, fitted and exact
    m1_coarse.PDS = (0.02, 0.04, 0.06, 0.091)
    m1_coarse.reads(RE, M.CUR["rules"], M.CUR["cal"], M.CUR["ss"], M.CUR["sp"], f"OB {tag0}, f_D {M.CUR['vals']['f_D']:.3f}")
    if "--broyden" not in args:
        return
    # 3. the other five rest-point targets by Broyden, f_D held (the D spread pair left out)
    CS.PAIRS = tuple(pr for pr in CS.PAIRS if pr[0] != "credit_spread_target_D")
    CS.TOL = np.array([pr[3] for pr in CS.PAIRS])
    x = np.array([M.CUR["cal"][pr[0]] for pr in CS.PAIRS])
    CS.get_calibration = lambda: dict(M.CUR["cal"], **{pr[0]: float(xi) for pr, xi in zip(CS.PAIRS, x)})
    CS.evaluate = lambda xx, settings, verbose=True: M.evaluate(RE, xx)
    xs, m = CS.calibrate(dict(), HERE / "runs" / f"m2_ob_{tag0}_broyden.json", max_evals=12)
    p, where = threshold(RE, M)
    hist.append(dict(stage="broyden", converged=m is not None, threshold=p, where=where,
                     x={pr[0]: float(v) for pr, v in zip(CS.PAIRS, xs)}, f_D=float(M.CUR["vals"]["f_D"])))
    save("final")
    print(f"\n[ob] {'CONVERGED' if m is not None else 'not converged'}; f_D {M.CUR['vals']['f_D']:.4f}; IC binds from "
          f"p^d {100 * p:.2f}% ({where}); " + "  ".join(f"{pr[0]}={v:.6g}" for pr, v in zip(CS.PAIRS, xs)), flush=True)


if __name__ == "__main__":
    main()
