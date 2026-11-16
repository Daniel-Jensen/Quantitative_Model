# THE TPI CAP WALK WITH GAUSS-SEIDEL NODE SWEEPS (2026-11-02). THE COARSE WALK STALLS ONE RUNG AFTER THE CAP FIRST BINDS
# (1150 -> 1137.5 bp AT FLOOR 0.5), YET THE BINDING NODE (THE p^d-12% TOP) SOLVED ALONE HAS A CLEAN ROOT AT EVERY CAP DOWN
# TO 200 bp (91% BOUGHT, THE D BANKS AT THEIR CORNER): THE GLOBAL NEWTON, NOT THE ECONOMICS, FAILS. HERE EACH CAP RUNG
# FIRST SWEEPS THE NODES ONE BY ONE (EACH SOLVED ALONE, KEPT IF ITS RESIDUAL FALLS), THEN POLISHES WITH THE GLOBAL NEWTON.
# OPTIONS: --sweeps N (0: NEWTON ONLY), --eps E (THE TPI SMOOTHING DURING THE WALK, ANNEALED TO cal["tpi_eps"] AT THE END,
# AS THE REFINED SOLVE DID ON 2026-10-05).
# --target T [--step S --min-step M]: AN ADAPTIVE WALK TO T INSTEAD OF A FIXED CAP LIST; --maxit N PER RUNG. A RUNG PICKLE
# (tpi_gs_<tag>_<cap>.pkl) RESUMES A WALK. --gz: THE WALK IN THE GARCIA-ZANGWILL FORM; --tag NAMES THE OUTPUTS.
# --no-anneal: THE WALK's SMOOTHING IS KEPT AT THE END.
# --presolve: BEFORE EACH RUNG's NEWTON, THE NODES THE STEP MOVES (TPI RESIDUAL > 1e-6) ARE SOLVED ALONE (switching()).
# --gaps [--margin M]: EACH RUNG's CAP KEPT AT LEAST M bp FROM EVERY NODE's SWITCHING CAP (next_cap()).
# --rung-tol T: AN INTERMEDIATE RUNG (ABOVE --target) IS ACCEPTED AT max|F| <= T; THE TARGET KEEPS THE FULL TOLERANCE.
# --active: THE NODES A RUNG SWITCHES ON START ON THE BUYING SIDE OF THEIR KINK (activate()).
# --reuse K: EACH RUNG FIRST RUNS NEWTON WITH THE JACOBIAN REFRESHED EVERY K ITERATIONS, THEN THE FULL NEWTON.
#   python3 -u tpi_gs_walk.py <pickle> [caps, e.g. 1137.5,1000,800] [--sweeps N] [--eps E] [--target T] [--maxit N]
import json, os, pickle, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np

CAPS = (1137.5, 1100.0, 1000.0, 800.0, 600.0, 400.0, 200.0)
SWEEPS = 6
NODE_TOL = 1e-10


def sweep(RE, r, c, ss, sp):
    # ONE GAUSS-SEIDEL PASS OVER THE NODES (BOTH REGIMES): EACH SOLVED ALONE, KEPT IF ITS RESIDUAL FALLS. RETURNS THE
    # LARGEST NODE RESIDUAL AFTER THE PASS AND THE NUMBER OF NODES MOVED.
    regimes, worst, moved = (0, 1), 0.0, 0
    for i in range(r.grid.n):
        z_now = RE._node_values(r, i, regimes)
        f_now = float(np.max(np.abs(RE._node_residual(r, c, ss, sp, i, regimes, False, z_now))))
        if f_now <= NODE_TOL:
            worst = max(worst, f_now)
            continue
        f_new = RE._node_solve(r, c, ss, sp, i, regimes, False, z_now)
        if np.isfinite(f_new) and f_new < f_now:
            moved += 1
        else:
            RE._node_residual(r, c, ss, sp, i, regimes, False, z_now)
            f_new = f_now
        worst = max(worst, f_new)
    return worst, moved


def switching(RE, r, c, ss, sp, tol=1e-6):
    # THE NODES THE CAP STEP MOVES OFF THEIR ROOT (THE TPI RESIDUAL ABOVE tol: THE ONES ABOUT TO BUY OR BUYING MORE), EACH
    # SOLVED ALONE FROM ITS CURRENT VALUES AND KEPT IF ITS RESIDUAL FALLS: THE KINK OF A NODE THAT SWITCHES TO BUYING IS
    # RESOLVED LOCALLY, SO THE GLOBAL NEWTON ONLY RECONCILES THE SMOOTH COUPLING. RETURNS (WORST NODE RESIDUAL, NODES MOVED).
    from global_projection.solver_recursive import collocation as C
    F = C.make_residual(r, c, ss, sp, (0, 1), False, RE.N_GH, None, ("normal",))
    f = np.abs(F(C.pack(r, (0, 1), ("normal",)))).reshape(2, r.grid.n, C.N_RES)
    itpi = C.RES_NAMES.index("tpi")
    worst, moved = 0.0, 0
    for i in np.argsort(-f[0, :, itpi]):
        if f[0, i, itpi] <= tol:
            break
        z_now = RE._node_values(r, i, (0, 1))
        f_now = float(np.max(np.abs(RE._node_residual(r, c, ss, sp, i, (0, 1), False, z_now))))
        f_new = RE._node_solve(r, c, ss, sp, i, (0, 1), False, z_now)
        if np.isfinite(f_new) and f_new < f_now:
            moved += 1
        else:
            RE._node_residual(r, c, ss, sp, i, (0, 1), False, z_now)
            f_new = f_now
        worst = max(worst, f_new)
    return worst, moved


def thresholds(RE, r, c, ss, sp, B):
    # THE CAPS AT WHICH THE NODES NOT YET BUYING START TO BUY: THEIR SOVEREIGN SPREADS (bp/yr) AT THE CURRENT RULES. A NODE
    # AT ITS SWITCH BUYS A SMOOTHING-SIZED SLIVER (~tpi_eps of B), SO "NOT YET BUYING" IS BELOW 1e-4 OF B.
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    g, out = r.grid, []
    for i in range(g.n):
        x = np.array([r.vals[k][0][i] for k in SOLVE])
        _, o = point_residuals(g.points[i], 0, x, r, c, ss, sp, RE.N_GH)
        if o["m_cb"] <= 1e-4 * B:
            out.append(4e4 * (c["delta_b_D"] * (1 / o["Q_bD"] - 1) - c["delta_b_F"] * (1 / o["Q_bF"] - 1)))
    return np.array(out)


def next_cap(last, step, thr, target, margin):
    # THE NEXT RUNG: THE CAP NEAREST last - step THAT SITS AT LEAST margin FROM EVERY NODE's SWITCHING CAP (A NODE EXACTLY
    # AT ITS SWITCH LEAVES THE COMPLEMENTARITY's TWO ARGUMENTS AT ZERO, WHERE NEWTON PLATEAUS AT ~1e-5); SEARCHED OUTWARD IN
    # 0.25 bp, BELOW last AND NOT BELOW target.
    cand = max(last - step, target)
    near = lambda x: bool(np.any(np.abs(thr - x) < margin))
    if cand <= target or not near(cand):
        return cand
    for d in np.arange(0.25, 2.0 * step + 4.0 * margin, 0.25):
        for x in (cand - d, cand + d):
            if target <= x < last - 0.5 and not near(x):
                return float(x)
    return cand


def activate(RE, r, c, ss, sp, B, last, cap, slope=0.005):
    # AN ACTIVE-SET WARM START: EVERY NODE NOT YET BUYING WHOSE SWITCHING CAP LIES ABOVE THE NEW CAP IS MOVED TO THE BUYING
    # SIDE OF ITS KINK (PURCHASE slope*(threshold - cap) OF B, ~0.5% PER bp AT THE 1002.5 bp SOLUTION; PRICE AT THE NEW
    # FLOOR) AND NOTHING ELSE, SO NEWTON's LINEARISATION STARTS ON THE RIGHT SIDE OF EVERY SWITCH. RETURNS THE NODES MOVED.
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    g, moved = r.grid, []
    xs, qs = r.vals["x_cb"][0].copy(), r.vals["Q_bD"][0].copy()
    for i in range(g.n):
        x = np.array([r.vals[k][0][i] for k in SOLVE])
        _, o = point_residuals(g.points[i], 0, x, r, c, ss, sp, RE.N_GH)
        if o["m_cb"] > 1e-4 * B:
            continue
        t = 4e4 * (c["delta_b_D"] * (1 / o["Q_bD"] - 1) - c["delta_b_F"] * (1 / o["Q_bF"] - 1))
        if cap < t < last + 1.0:
            xs[i], qs[i] = slope * (t - cap), o["Q_floor"]
            moved.append(i)
    if moved:
        r.set_values("x_cb", 0, xs)
        r.set_values("Q_bD", 0, qs)
    return moved


def main():
    from tpi_floor_sweep import economy, book
    from global_projection.solver_recursive import recursive_experiment as RE
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    args = sys.argv[2:]
    opt = lambda k, d: type(d)(args[args.index(k) + 1]) if k in args else d
    caps = tuple(float(v) for v in args[0].split(",")) if args and not args[0].startswith("--") else CAPS
    sweeps, eps, maxit = opt("--sweeps", SWEEPS), opt("--eps", 0.0), opt("--maxit", 20)
    presolve, active = "--presolve" in args, "--active" in args
    reuse = opt("--reuse", 1)
    rung_tol, target_cap = opt("--rung-tol", 0.0), opt("--target", 0.0)
    # a floor-sweep pickle (its TPI rules) or a rung pickle written below (its rules, cap and smoothing: a resume);
    # RuleSet.copy() drops tpi_cap_solved, so the last solved cap is tracked here
    if "on" in sv:
        r, last = sv["on"].copy(), float(sv["on"].tpi_cap_solved)
    else:
        r, last = sv["rules"].copy(), float(sv["cap"])
        eps = eps or float(sv["tpi_eps"])
    floor = sv["floor"]
    cal, ss, sp = economy(floor)
    # a start written for another box (bocola_size.py: pd_box_top 0.062) carries its calibration changes
    extra = sv.get("cal_extra", {})
    cal.update(extra)
    c = dict(cal, tpi_on=True, tpi_cap_bp=last)
    hard = RE._soften(c)
    eps_target = c["tpi_eps"]
    if eps > eps_target:
        c["tpi_eps"] = eps
    # --gz: the walk roots the purchase rule in the Garcia-Zangwill form (one smooth variable per node) instead of the
    # Fischer-Burmeister pair, whose two near-zero arguments stall Newton at a node's switch point (2026-11-03)
    gz_walk = "--gz" in args
    RE._tpi_form(r, c, ss, gz=gz_walk)
    tag = (args[args.index("--tag") + 1] if "--tag" in args else sv.get("tag")
           or Path(sys.argv[1]).stem + (f"_eps{eps:g}" if eps else "") + (f"_gs{sweeps}" if sweeps != SWEEPS else ""))
    hist = []

    def rung(cap):
        # ONE CAP RUNG: THE SWEEPS, THEN THE GLOBAL NEWTON; THE RULES MOVE ONLY ON SUCCESS, EACH SOLVED RUNG IS PICKLED.
        nonlocal r, last
        t0 = time.perf_counter()
        good = r.copy()
        c["tpi_cap_bp"] = cap
        w_s, n = np.nan, 0
        if active:
            ids = activate(RE, r, c, ss, sp, cal["B_gov_D_ss"], last, cap)
            print(f"  [gs {cap:.1f} bp] active set: {len(ids)} nodes moved to buying {ids}", flush=True)
        if presolve:
            w_s, n = switching(RE, r, c, ss, sp)
            print(f"  [gs {cap:.1f} bp] presolve: {n} nodes solved alone, worst node residual {w_s:.1e}", flush=True)
        for k in range(sweeps):
            w_s, n = sweep(RE, r, c, ss, sp)
            print(f"  [gs {cap:.1f} bp] sweep {k + 1}: worst node residual {w_s:.1e}, {n} nodes moved", flush=True)
            if w_s <= NODE_TOL or n == 0:
                break
        ok_r, w_r = False, np.inf
        if reuse > 1:
            # a chord pass first (the Jacobian refreshed every --reuse iterations: seconds per iteration, not minutes);
            # where it stalls above the rung tolerance, the full Newton continues from its iterate
            from global_projection.solver_recursive.collocation import solve_collocation
            ok_r, _, w_r = solve_collocation(r, c, ss, sp, regimes=(0, 1), no_default=False, maxit=maxit, verbose=True,
                                             label=f" TPI cap={cap:.1f}bp (reuse {reuse})", jac_every=reuse)
        if ok_r or (cap > target_cap and w_r <= rung_tol):
            ok, w = ok_r, w_r
        else:
            ok, _, w = RE._stage(r, c, ss, sp, (0, 1), False, f"TPI cap={cap:.1f}bp", True, warm=0, maxit=maxit)
        # an intermediate rung only seeds the next one, so it is accepted at --rung-tol; the target needs the full tolerance
        if not ok and cap > target_cap and np.isfinite(w) and w <= rung_tol:
            print(f"  [gs] {cap:.1f} bp accepted as a path rung at max|F| {w:.1e} (rung tolerance {rung_tol:.0e})", flush=True)
            ok = True
        rec = dict(cap=cap, tpi_eps=c["tpi_eps"], sweeps_worst=float(w_s), newton_ok=bool(ok), newton_worst=float(w),
                   seconds=time.perf_counter() - t0)
        if ok:
            last = cap
            # the footprint in the read (GZ) form, on copies: _tpi_form switches cal["tpi_gz"] in place
            gz, cg = r.copy(), dict(c)
            gz.tpi_cap_solved = cap
            RE._tpi_form(gz, cg, ss, gz=True)
            rec["book"] = book(RE, gz, cg, ss, sp)
            with open(HERE / "runs" / f"tpi_gs_{tag}_{cap:.0f}.pkl", "wb") as fh:
                pickle.dump(dict(rules=r, floor=floor, cap=cap, tpi_eps=c["tpi_eps"], tag=tag, cal_extra=extra,
                                 smoothing={k: c[k] for k in hard}), fh)
        else:
            r = good
        print(f"[gs] cap {cap:.1f} bp (tpi_eps {c['tpi_eps']:g}): Newton ok {ok} ({w:.1e}); footprint {rec.get('book')} "
              f"({rec['seconds']:.0f}s)", flush=True)
        hist.append(rec)
        (HERE / "runs" / f"tpi_gs_{tag}.json").write_text(json.dumps(hist, indent=1, default=lambda z: z.item()))
        return ok

    if "--target" in args:
        # the adaptive walk: the step halves after a failed rung and grows by --grow after a solved one (to 100 bp)
        target, step, min_step = opt("--target", 200.0), opt("--step", 50.0), opt("--min-step", 12.5)
        # a failed rung costs several solved ones on the fine grid (2026-11-02: 2 h against 20 min), so --grow 1 holds the step
        grow = opt("--grow", 1.5)
        gaps, margin = "--gaps" in args, opt("--margin", 2.0)
        while last > target and step >= min_step:
            cap = max(last - step, target)
            if gaps:
                thr = thresholds(RE, r, c, ss, sp, cal["B_gov_D_ss"])
                cap = next_cap(last, step, thr, target, margin)
                print(f"  [gs] next rung {cap:.2f} bp (step {step:.2f}); thresholds ahead "
                      f"{', '.join(f'{t:.1f}' for t in sorted(thr[thr < last], reverse=True)[:6])}", flush=True)
            if rung(cap):
                step = min(grow * step, 100.0)
            else:
                step *= 0.5
        print(f"[gs] walk ended at {last:.1f} bp (target {target:.0f} bp, step {step:.1f} bp)", flush=True)
    else:
        for cap in caps:
            if not rung(cap):
                print(f"[gs] stopped: cap {cap:.1f} bp not solved; last solved {last:.1f} bp", flush=True)
                break
    # the TPI smoothing annealed to its calibrated value at the last solved cap (--no-anneal: kept, for a later resume)
    c["tpi_cap_bp"] = last
    for e in [] if "--no-anneal" in args else [e for e in (1e-3,) if eps_target < e < c["tpi_eps"]] + ([eps_target] if c["tpi_eps"] > eps_target else []):
        c["tpi_eps"] = e
        ok, _, w = RE._stage(r, c, ss, sp, (0, 1), False, f"TPI cap={last:.1f}bp tpi_eps={e:g}", True, warm=0, maxit=maxit)
        print(f"[gs] anneal at {last:.1f} bp: tpi_eps {e:g} ok {ok} ({w:.1e})", flush=True)
        hist.append(dict(cap=last, anneal_eps=e, newton_ok=bool(ok), newton_worst=float(w)))
        (HERE / "runs" / f"tpi_gs_{tag}.json").write_text(json.dumps(hist, indent=1, default=lambda z: z.item()))
        if not ok:
            break
        with open(HERE / "runs" / f"tpi_gs_{tag}_{last:.0f}_eps{e:g}.pkl", "wb") as fh:
            pickle.dump(dict(rules=r, floor=floor, cap=last, tpi_eps=e, tag=tag, cal_extra=extra, smoothing={k: c[k] for k in hard}), fh)


if __name__ == "__main__":
    main()
