# THE TPI AT BOCOLA's (2016) SCALE (2026-11-03). HIS REPLICATION (replication/bocola2016) SOLVES ON s IN [-11.41, -2.71], SO
# p^d FROM 0.001% TO 6.2% A QUARTER, WITH 9 NODES PER AXIS (ISOTROPIC SMOLYAK LEVEL 3, ~389 POINTS IN 6 STATES), AND HIS IRF
# SHOCK RAISES p^d TO 2.5%. A LEVEL-3 GRID IN OUR 12 STATES (~2650 POINTS, ~106k UNKNOWNS) IS NOT FEASIBLE; THE FINE CROSS
# GRID KEEPS HIS 9 NODES ON THE RISK AXIS AND HIS SYSTEM SIZE (~3000 UNKNOWNS), SO HERE ITS RISK BOX IS CUT TO HIS 6.2% TOP.
#   stage 1: the no-TPI baseline on that box, seeded from the solved 12%-top rules        -> runs/bocola_base.pkl
#   stage 2: its economics: rest point, the p^d where the D IC binds, the sovereign spread at every s level of the grid
#            (where a 200 bp cap binds), and the impact of the 2.5% shock                 -> runs/bocola_base.json
#   stage 3: the TPI start for tpi_gs_walk.py (the cap just above the highest node spread) -> runs/tpi_start_bocola.pkl
#   python3 -u bocola_size.py [--from-stage N]
import json, os, pickle, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np

PD_TOP = 0.062
SHOCK_PD = 0.025
FLOOR = 0.5


def economy():
    # run.py's CALIBRATION WITH BOCOLA's RISK BOX TOP.
    from tpi_floor_sweep import economy as eco
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    cal, ss, sp = eco(FLOOR)
    cal["pd_box_top"] = PD_TOP
    calibrate_household_anchors(cal, ss, sp)
    return cal, ss, sp


def sov(o, cal):
    # THE SOVEREIGN SPREAD OVER F, bp/yr.
    return 4e4 * (cal["delta_b_D"] * (1 / o["Q_bD"] - 1) - cal["delta_b_F"] * (1 / o["Q_bF"] - 1))


def stage1(cal, ss, sp):
    # THE NO-TPI BASELINE ON THE 6.2% BOX, SEEDED FROM THE 12% RULES (THE NEW BOX LIES INSIDE THE OLD: INTERPOLATION).
    import run
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive.state_grid import IS
    fine = pickle.load(open(run.WARM_START, "rb"))
    grid = RE.normal_box(ss, cal, mu=1, refine=(IS, 9, "cross"))
    print(f"[bocola] box s [{grid.lo[IS]:.3f}, {grid.hi[IS]:.3f}] (p^d {100 / (1 + np.exp(-grid.lo[IS])):.4f}% to "
          f"{100 / (1 + np.exp(-grid.hi[IS])):.2f}%), {grid.n} points; the 12% box had s top {fine.grid.hi[IS]:.3f}", flush=True)
    seed = RE._seed_from(RuleSet(grid, 2), fine)
    t0 = time.perf_counter()
    base = RE.solve_recursive(dict(cal), ss, sp, base=seed, s_refine=0)
    print(f"[bocola] baseline: ok {base.solve_ok}, worst {base.solve_worst:.1e} ({time.perf_counter() - t0:.0f}s)", flush=True)
    with open(HERE / "runs" / "bocola_base.pkl", "wb") as fh:
        pickle.dump(dict(rules=base, floor=FLOOR, pd_box_top=PD_TOP), fh)
    return base


def stage2(base, cal, ss, sp):
    # THE NO-TPI ECONOMICS ON BOCOLA's BOX (EXACT READS).
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    from global_projection.solver_recursive.point_map import _tpi_floor
    S0 = RE.stochastic_rest_point(base, cal, ss, sp, verbose=False)
    o0 = RE.read_exact(base, cal, ss, sp, S0)
    sprd = lambda o: 4e4 * (o["r_wc_D"] - o["rdep_D"])
    read = lambda pd: RE.read_exact(base, cal, ss, sp, np.where(np.arange(S0.size) == IS, RE.s_from_pd(pd), S0))
    # the p^d where the D IC starts to bind (bisection on mu_D, the other states at rest)
    lo, hi = 1 / (1 + np.exp(-S0[IS])), 0.98 * PD_TOP
    thr = None if read(hi)["mu_D"] <= 0 else hi
    while thr is not None and hi / lo - 1 > 1e-3:
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if read(mid)["mu_D"] <= 0 else (lo, mid)
        thr = hi
    out = dict(rest=dict(mu_D=float(o0["mu_D"]), sov=sov(o0, cal), E_Om=float(o0["E_Om_D"]), Q_bD=float(o0["Q_bD"])),
               ic_binds_from=thr)
    print(f"[bocola] rest point: mu_D {o0['mu_D']:.4f}, sovereign spread {sov(o0, cal):.0f} bp, E[Om] {o0['E_Om_D']:.3f}; "
          f"the D IC binds from p^d {'never in the box' if thr is None else f'{100 * thr:.2f}%'}", flush=True)
    for pd in (0.02, SHOCK_PD, 0.04, 0.05):
        o = read(pd)
        row = dict(mu_D=float(o["mu_D"]), sov=sov(o, cal), spread=sprd(o) - sprd(o0),
                   **{f"d{k}": 100 * (o[k] / o0[k] - 1) for k in ("Q_bD", "n_D", "n_F", "Y_D", "C_D", "I_D")})
        out[f"impact_{100 * pd:.1f}"] = row
        print(f"[bocola] impact at p^d {100 * pd:.1f}%: sovereign {row['sov']:.0f} bp, lending spread {row['spread']:+.0f} bp, "
              f"Q_bD {row['dQ_bD']:+.2f}%, n_D {row['dn_D']:+.2f}%, n_F {row['dn_F']:+.3f}%, Y_D {row['dY_D']:+.3f}%, "
              f"C_D {row['dC_D']:+.3f}%, I_D {row['dI_D']:+.2f}%, mu_D {row['mu_D']:.3f}", flush=True)
    # the sovereign spread at every s level of the grid, every node (where a 200 bp cap binds, and how many nodes at once)
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    g, levels = base.grid, {}
    for i in range(g.n):
        x = np.array([base.vals[k][0][i] for k in SOLVE])
        _, o = point_residuals(g.points[i], 0, x, base, dict(cal, tpi_on=False), ss, sp, RE.N_GH)
        levels.setdefault(round(float(g.points[i][IS]), 8), []).append(sov(o, cal))
    out["levels"] = []
    for s, spreads in sorted(levels.items()):
        a = np.array(spreads)
        out["levels"].append(dict(pd=float(1 / (1 + np.exp(-s))), nodes=len(a), sov_min=float(a.min()), sov_max=float(a.max())))
        print(f"[bocola] s level p^d {100 / (1 + np.exp(-s)):7.3f}%: {len(a):2d} nodes, sovereign spread {a.min():6.0f}"
              f"-{a.max():6.0f} bp", flush=True)
    (HERE / "runs" / "bocola_base.json").write_text(json.dumps(out, indent=1, default=lambda z: z.item()))
    return out


def stage3(base, out):
    # THE TPI WALK's START: THE BASELINE IS THE TPI SOLUTION AT ANY CAP ABOVE EVERY NODE's SPREAD.
    top = max(l["sov_max"] for l in out["levels"])
    base.tpi_cap_solved = float(np.ceil(top / 25.0) * 25.0 + 25.0)
    with open(HERE / "runs" / "tpi_start_bocola.pkl", "wb") as fh:
        pickle.dump(dict(on=base, floor=FLOOR, tag="bocola", cal_extra=dict(pd_box_top=PD_TOP)), fh)
    print(f"[bocola] TPI start written: cap {base.tpi_cap_solved:.0f} bp (highest node spread {top:.0f} bp)", flush=True)


def main():
    args = sys.argv[1:]
    first = int(args[args.index("--from-stage") + 1]) if "--from-stage" in args else 1
    cal, ss, sp = economy()
    base = stage1(cal, ss, sp) if first <= 1 else pickle.load(open(HERE / "runs" / "bocola_base.pkl", "rb"))["rules"]
    out = stage2(base, cal, ss, sp) if first <= 2 else json.loads((HERE / "runs" / "bocola_base.json").read_text())
    stage3(base, out)


if __name__ == "__main__":
    main()
