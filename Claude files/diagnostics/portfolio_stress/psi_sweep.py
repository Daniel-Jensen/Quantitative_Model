# STRESS TEST: HOW MUCH ONE CALIBRATION DIAL DRIVES THE SOVEREIGN-RISK IRF OF THE GLOBAL MODEL.
#
#   python3 "Claude files/diagnostics/portfolio_stress/psi_sweep.py"                       # psi, quick grid
#   python3 "Claude files/diagnostics/portfolio_stress/psi_sweep.py" --values 0.5 2 8 --s-refine 5
#   python3 "Claude files/diagnostics/portfolio_stress/psi_sweep.py" --param kappa_nfa --values 0.0025 0.01 0.04
#   python3 "Claude files/diagnostics/portfolio_stress/psi_sweep.py" --set deposit_market=national   # domestic-funding-fix only
#   python3 "Claude files/diagnostics/portfolio_stress/psi_sweep.py" --compare                       # rank every saved sweep
#
# Each value is a full global re-solve (warm-started from its neighbour when the SS is unchanged),
# read as Bocola's differenced IRF from the model's own rest point. Three measurements per value:
# the IRF itself, the impact at the REFERENCE rules (the within-period share of the dial's effect),
# and the accuracy guards (solver verdict, box clipping, fitted-vs-exact gap).
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "results" / "STRESS"
PSI = "psi_bF_D,psi_bD_F"
LEVELS = ("Y_D", "C_D", "I_D", "n_D", "n_F", "Q_bD", "Y_F", "b_D_F_new", "b_F_D_new")
RATES = ("rdep_D", "r_wc_D", "wedge_sp_D")
METRICS = ("Y_impact", "Y_trough", "Y_cum", "I_cum", "C_impact", "n_D_impact", "n_F_impact",
           "Q_bD_impact", "spread_impact_bp", "rdep_impact_bp", "sov_impact_bp",
           "share_F_q4_pp", "nfa_q8_pctY", "retrench_elast")


def _parse(v):
    # A COMMAND-LINE VALUE AS BOOL, FLOAT OR STRING.
    if v in ("True", "False"):
        return v == "True"
    try:
        return float(v)
    except ValueError:
        return v


def setup(keys, value, overrides):
    # CALIBRATION, SS AND HOUSEHOLD ANCHORS WITH THE DIAL SET.
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    cal = get_calibration()
    cal["nw_floor_frac"] = 0.15  # run.py's NW_FLOOR
    cal.update(overrides)
    for k in keys:
        cal[k] = value
    ss = solve_steady_state(cal)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    return cal, ss, sproc


def _yield(Q, db):
    # HM PERPETUITY FLOW YIELD, PER QUARTER.
    return db * (1.0 - Q) / Q


def _row(o, cal):
    # THE OBJECTS RECORDED EACH QUARTER.
    r = {k: float(o[k]) for k in LEVELS + RATES}
    r["share_F"] = float(o["b_D_F_new"] / o["Bp_D"])
    # D's net external position: the deposit leg plus the net cross-border sovereign book
    r["nfa"] = float(o["nfa_dep_D"] + o["p"] * o["Q_bF"] * o["b_F_D_new"]
                     - o["Q_bD"] * o["b_D_F_new"])
    r["sov"] = _yield(r["Q_bD"], cal["delta_b_D"]) - _yield(float(o["Q_bF"]), cal["delta_b_F"])
    return r


def irf(rules, cal, ss, sproc, pd, T):
    # BOCOLA'S DIFFERENCED IRF FROM THE REST POINT, WITH THE LARGEST BOX CLIP ON THE WAY.
    from global_projection.solver_recursive.recursive_experiment import (
        read_at, advance, stochastic_rest_point, s_from_pd)
    from global_projection.solver_recursive.state_grid import IS
    S0 = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    paths, clip = {}, 0.0
    for label, s0 in (("ref", None), ("shock", s_from_pd(pd))):
        S = S0.copy()
        if s0 is not None:
            S[IS] = s0
        rows = []
        for _ in range(T):
            o = read_at(rules, cal, ss, sproc, S)
            rows.append(_row(o, cal))
            Sn = advance(o, S, sproc)
            S = rules.grid.clip(Sn)[0]
            clip = max(clip, float(np.max(np.abs(Sn - S) / np.maximum(np.abs(Sn), 1e-12))))
        paths[label] = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    ref, sh = paths["ref"], paths["shock"]
    Yss = ss["ss_firm_D"]["Y_ss"]
    out = {k: 100 * (sh[k] / ref[k] - 1) for k in LEVELS}
    out.update({k: 4e4 * (sh[k] - ref[k]) for k in RATES + ("sov",)})
    out["share_F"] = 100 * (sh["share_F"] - ref["share_F"])
    out["nfa"] = 100 * (sh["nfa"] - ref["nfa"]) / Yss
    return out, S0, clip


def metrics(r):
    # SCALAR SUMMARIES OF ONE IRF.
    Y = r["Y_D"]
    m = dict(Y_impact=Y[0], Y_trough=Y.min(), Y_trough_q=int(Y.argmin()), Y_cum=Y.sum(),
             I_cum=r["I_D"].sum(), C_impact=r["C_D"][0], n_D_impact=r["n_D"][0],
             n_F_impact=r["n_F"][0], Q_bD_impact=r["Q_bD"][0],
             spread_impact_bp=r["wedge_sp_D"][0], rdep_impact_bp=r["rdep_D"][0],
             sov_impact_bp=r["sov"][0], share_F_q4_pp=r["share_F"][4], nfa_q8_pctY=r["nfa"][8])
    # model analogue of the retrenchment moment: pp of D debt shed by F banks per 100 bp of spread
    m["retrench_elast"] = m["share_F_q4_pp"] / max(abs(m["sov_impact_bp"]) / 100, 1e-9)
    return {k: float(v) for k, v in m.items()}


def impact_gap(rules, cal, ss, sproc, S0, pd):
    # FITTED VS EXACTLY CLEARED IMPACT OUTPUT (THE IDENTIFICATION GAP AT THE KKT KINK).
    from global_projection.solver_recursive.recursive_experiment import read_at, read_exact, s_from_pd
    from global_projection.solver_recursive.state_grid import IS
    S = S0.copy()
    S[IS] = s_from_pd(pd)
    fit = read_at(rules, cal, ss, sproc, S)["Y_D"] / read_at(rules, cal, ss, sproc, S0)["Y_D"]
    ex = read_exact(rules, cal, ss, sproc, S)["Y_D"] / read_exact(rules, cal, ss, sproc, S0)["Y_D"]
    return 100 * (fit - 1), 100 * (ex - 1)


def _try(fn):
    # A SOLVE THAT RAISES (NaN ITERATE, SINGULAR STEP) COUNTS AS A FAILED SOLVE.
    try:
        rules = fn()
        return rules if rules.solve_ok else None
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as e:
        print(f"    solve raised {type(e).__name__}: {e}", flush=True)
        return None


def solve(keys, v, v_seed, overrides, s_refine, warm, substeps=(1, 2, 4)):
    # ONE GLOBAL SOLVE: WARM FROM THE NEIGHBOUR, THEN IN 2 AND 4 LOG-STEPS, THEN COLD.
    from global_projection.solver_recursive.recursive_experiment import solve_recursive
    if warm is not None:
        for k in substeps:
            rules = warm
            for u in np.geomspace(v_seed, v, k + 1)[1:]:
                cal, ss, sproc = setup(keys, float(u), overrides)
                rules = _try(lambda: solve_recursive(cal, ss, sproc, s_refine=s_refine,
                                                     base=rules.copy()))
                if rules is None:
                    break
            if rules is not None:
                return rules, f"warm/{k}"
    cal, ss, sproc = setup(keys, v, overrides)
    rules = _try(lambda: solve_recursive(cal, ss, sproc, s_refine=s_refine, base_out=[]))
    return rules, "cold"


def _same_ss(a, b):
    # TRUE WHEN TWO STEADY STATES SHARE THE GRID BOX (SO A WARM START IS VALID).
    return all(np.allclose(a[k], b[k]) for k in ("Kap_D_ss", "C_D_ss", "C_F_ss", "Q_bD_ss"))


def walk_order(values, ref):
    # THE REFERENCE FIRST, THEN OUTWARD IN EACH DIRECTION, EACH STEP SEEDED BY ITS NEIGHBOUR.
    vals = sorted(set(values) | {ref})
    up = [v for v in vals if v > ref]
    down = [v for v in vals if v < ref][::-1]
    return [(ref, None)] + [(v, (up[i - 1] if i else ref)) for i, v in enumerate(up)] \
        + [(v, (down[i - 1] if i else ref)) for i, v in enumerate(down)]


def elasticities(res, ref):
    # d ln|m| / d ln(dial) AT THE REFERENCE, FROM ITS TWO NEIGHBOURS; AND THE SWEEP RANGE.
    vals = sorted(res)
    i = vals.index(ref)
    lo, hi = vals[max(i - 1, 0)], vals[min(i + 1, len(vals) - 1)]
    out = {}
    for k in METRICS:
        m = {v: res[v]["metrics"][k] for v in vals}
        e = np.nan
        if hi != lo and lo > 0 and m[lo] * m[hi] > 0:
            e = (np.log(abs(m[hi])) - np.log(abs(m[lo]))) / (np.log(hi) - np.log(lo))
        rng = (max(m.values()) - min(m.values())) / max(abs(m[ref]), 1e-12)
        out[k] = dict(elasticity=float(e), range_rel=float(rng))
    return out


def run(args):
    # THE SWEEP.
    keys = args.param.split(",")
    overrides = dict(kv.split("=", 1) for kv in args.set)
    overrides = {k: _parse(v) for k, v in overrides.items()}
    from calibration.global_projection import get_calibration
    ref = float(args.ref if args.ref is not None else get_calibration()[keys[0]])
    tag = f"{keys[0]}_{overrides.get('deposit_market', 'default')}_s{args.s_refine}"
    OUT.mkdir(parents=True, exist_ok=True)
    res, rules_by, ss_by, failed = {}, {}, {}, []
    for v, seed in walk_order(args.values, ref):
        t0 = time.perf_counter()
        cal, ss, sproc = setup(keys, v, overrides)
        warm = (rules_by.get(seed) if seed in rules_by and _same_ss(ss, ss_by[seed])
                else None)
        rules, how = solve(keys, v, seed, overrides, args.s_refine, warm)
        if rules is None:
            # a value the ladder cannot solve is itself a result: record it, keep sweeping
            failed.append(v)
            print(f"[{keys[0]} = {v:g}] NO SOLUTION (warm walk and cold ladder both failed)",
                  flush=True)
            continue
        r, S0, clip = irf(rules, cal, ss, sproc, args.pd, args.T)
        fit, ex = impact_gap(rules, cal, ss, sproc, S0, args.pd)
        rules_by[v], ss_by[v] = rules, ss
        res[v] = dict(value=v, start=how, solve_ok=bool(rules.solve_ok),
                      solve_worst=float(rules.solve_worst), max_box_clip=clip,
                      impact_fitted=fit, impact_exact=ex, metrics=metrics(r),
                      irf={k: a.tolist() for k, a in r.items()},
                      seconds=time.perf_counter() - t0)
        m = res[v]["metrics"]
        print(f"[{keys[0]} = {v:g}] {how} ok={rules.solve_ok} clip={clip:.1e}  "
              f"Y0 {m['Y_impact']:+.4f}%  trough {m['Y_trough']:+.4f}% (q{m['Y_trough_q']})  "
              f"cumY {m['Y_cum']:+.3f}  shareF(q4) {m['share_F_q4_pp']:+.3f}pp  "
              f"gap fit/exact {fit:+.4f}/{ex:+.4f}  ({res[v]['seconds']:.0f}s)", flush=True)
        _save(OUT / f"{tag}.json", keys, ref, overrides, args, res, failed)

    if ref not in res:
        print(f"  the reference {ref:g} did not solve; no decomposition. Failed: {failed}")
        return
    # the dial's effect at FIXED reference rules: what the period map alone does
    from global_projection.solver_recursive.recursive_experiment import (
        read_exact, stochastic_rest_point, s_from_pd)
    from global_projection.solver_recursive.state_grid import IS
    cal_r, ss_r, sp_r = setup(keys, ref, overrides)
    S0 = stochastic_rest_point(rules_by[ref], cal_r, ss_r, sp_r, verbose=False)
    S1 = S0.copy()
    S1[IS] = s_from_pd(args.pd)
    for v in res:
        cal, ss, sproc = setup(keys, v, overrides)
        y = (read_exact(rules_by[ref], cal, ss, sproc, S1)["Y_D"]
             / read_exact(rules_by[ref], cal, ss, sproc, S0)["Y_D"])
        res[v]["impact_fixed_rules_exact"] = 100 * (y - 1)
    for v in res:
        d_full = res[v]["impact_exact"] - res[ref]["impact_exact"]
        d_fix = res[v]["impact_fixed_rules_exact"] - res[ref]["impact_fixed_rules_exact"]
        res[v]["within_period_share"] = float(d_fix / d_full) if abs(d_full) > 1e-9 else np.nan
    el = elasticities(res, ref)
    _save(OUT / f"{tag}.json", keys, ref, overrides, args, res, failed, el)
    report(keys[0], ref, res, el, failed)
    plot(OUT / f"{tag}.png", keys[0], ref, res)


def _save(path, keys, ref, overrides, args, res, failed, el=None):
    # WRITE THE SWEEP SO FAR.
    path.write_text(json.dumps(dict(param=keys, ref=ref, overrides=overrides, pd=args.pd,
                                    T=args.T, s_refine=args.s_refine,
                                    results={f"{v:g}": r for v, r in sorted(res.items())},
                                    failed=failed, elasticities=el), indent=1, default=float))


def report(name, ref, res, el, failed):
    # CONSOLE TABLE: EVERY VALUE, THEN THE ELASTICITIES AT THE REFERENCE.
    print(f"\n  {name}: sweep against the reference {ref:g}")
    hdr = ("value", "ok", "Y0%", "trough%", "cumY", "cumI", "n_D0%", "spr0bp", "rdep0bp",
           "shF4pp", "nfa8%Y", "fixed-rule share")
    print("  " + "".join(f"{h:>10}" for h in hdr))
    for v in sorted(res):
        r, m = res[v], res[v]["metrics"]
        print("  " + "".join(f"{x:>10}" for x in (
            f"{v:g}", str(r["solve_ok"])[0], f"{m['Y_impact']:+.4f}", f"{m['Y_trough']:+.4f}",
            f"{m['Y_cum']:+.3f}", f"{m['I_cum']:+.3f}", f"{m['n_D_impact']:+.3f}",
            f"{m['spread_impact_bp']:+.1f}", f"{m['rdep_impact_bp']:+.1f}",
            f"{m['share_F_q4_pp']:+.3f}", f"{m['nfa_q8_pctY']:+.3f}",
            f"{r.get('within_period_share', np.nan):.2f}")))
    print(f"\n  elasticity d ln|m| / d ln {name} at {ref:g}, and the sweep range as a share of |m(ref)|")
    for k, e in el.items():
        print(f"    {k:18s} {e['elasticity']:+8.3f}   range {e['range_rel']:7.2%}")
    if failed:
        print(f"  no solution at {name} = {', '.join(f'{v:g}' for v in failed)}")


def plot(path, name, ref, res):
    # SIX IRF PANELS, ONE LINE PER VALUE.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = (("Y_D", "output D, %"), ("I_D", "investment D, %"), ("n_D", "D bank net worth, %"),
              ("wedge_sp_D", "credit spread, bp/yr"), ("share_F", "F-bank share of D debt, pp"),
              ("nfa", "D net external position, % of Y"))
    fig, ax = plt.subplots(2, 3, figsize=(12, 6.5), sharex=True)
    cmap = plt.get_cmap("viridis")
    vals = sorted(res)
    for j, v in enumerate(vals):
        for a, (k, _) in zip(ax.flat, panels):
            a.plot(res[v]["irf"][k], color=cmap(j / max(len(vals) - 1, 1)),
                   lw=2.2 if v == ref else 1.2, label=f"{name} = {v:g}")
    for a, (_, t) in zip(ax.flat, panels):
        a.set_title(t, fontsize=10)
        a.axhline(0, color="0.6", lw=0.6)
    ax.flat[0].legend(fontsize=7)
    for a in ax[1]:
        a.set_xlabel("quarter")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"  figure -> {path}")


def compare():
    # RANK EVERY SAVED SWEEP BY THE ELASTICITY OF CUMULATIVE OUTPUT (THE TORNADO, AS A TABLE).
    rows = []
    for f in sorted(OUT.glob("*.json")):
        d = json.loads(f.read_text())
        if d.get("elasticities"):
            e = d["elasticities"]
            rows.append((f.stem, e["Y_cum"]["elasticity"], e["Y_cum"]["range_rel"],
                         e["Y_impact"]["elasticity"], e["Y_trough"]["elasticity"]))
    rows.sort(key=lambda r: -abs(np.nan_to_num(r[1])))
    print(f"  {'sweep':40s} {'el. cumY':>9} {'range cumY':>11} {'el. Y0':>8} {'el. trough':>11}")
    for r in rows:
        print(f"  {r[0]:40s} {r[1]:+9.3f} {r[2]:11.1%} {r[3]:+8.3f} {r[4]:+11.3f}")


def main():
    # COMMAND LINE.
    ap = argparse.ArgumentParser()
    ap.add_argument("--param", default=PSI, help="comma-separated keys set to the same value")
    ap.add_argument("--values", type=float, nargs="+", default=[0.25, 0.5, 1.0, 4.0, 8.0, 32.0])
    ap.add_argument("--ref", type=float, default=None, help="reference value (default: calibrated)")
    ap.add_argument("--set", nargs="*", default=[], help="extra calibration overrides key=value")
    ap.add_argument("--s-refine", type=int, default=0, help="0 = coarse (quick) grid, 5 = full")
    ap.add_argument("--pd", type=float, default=0.02, help="headline shock, p^d per quarter")
    ap.add_argument("--T", type=int, default=25)
    ap.add_argument("--compare", action="store_true")
    args = ap.parse_args()
    compare() if args.compare else run(args)


if __name__ == "__main__":
    # the collocation Jacobian pool re-imports this file in every worker
    main()
