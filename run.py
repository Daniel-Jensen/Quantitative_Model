"""Unified entry point: the linear (SSJ) model, the global model, or both, then a comparison.

    python3 run.py                    run MODEL with QUICK (both set below) -- what a
                                      "Run" click in the editor does
    python3 run.py ssj|global|both    override MODEL
    python3 run.py global --quick     coarse grid: a preview, NOT converged (= QUICK = True)
    python3 run.py both --plots-only  redraw every figure from the saved data, no solving
    python3 run.py compare            linear-vs-global tables and overlays only

THE PIPELINE, for every selected model:
    1. compute  -> results/<MODEL>/data/     every computed object, as plain data
    2. plot     -> results/<MODEL>/figures/  drawn from the SAVED data, after all computing
    3. compare  -> results/COMPARISON/       when both models ran
Console output of each model is also written to results/<MODEL>/run.log.

The two solvers need different interpreters (CLAUDE.md "Environment"): the GLOBAL model
runs in this process (plain python3), the SSJ model in a subprocess under SSJ_PYTHON.
"""
import argparse
import contextlib
import os
import pickle
import shutil
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")  # figures are files; never open a window

import results_io
from results_io import RESULTS, load, save

ROOT = Path(__file__).resolve().parent

# ================================ CONFIGURATION ================================
MODEL = "global"  # "ssj", "global" or "both"; the command line overrides it
QUICK = False  # True: the coarse grid only, a quick preview (not converged); --quick forces it
# the SSJ model's interpreter (the base env's liblapack is broken)
SSJ_PYTHON = os.environ.get("SSJ_PYTHON", "/opt/anaconda3/envs/ssj/bin/python")

# global model settings; its parameters live in calibration/global_projection.py
NW_FLOOR = 0.5  # the funded recapitalisation floor in default, a share of n_ss (Bocola's 0.15 until 2026-10-26)
MU = 1  # Smolyak level of the TFP grid
RISK_MU_VEC = None  # per-state Smolyak levels for the risk grid; None = mu 1
S_REFINE = 9  # nodes in s: 9 = Bocola's (75 points with recursive_experiment.S_REFINE_MODE "cross"), 5, 0 = coarse only
ROTATE_P = False  # eigenbasis box for P: off, it measures worse
ACCURACY_T = 1200  # simulated periods for the Euler-error report
DECOMP_T = 25  # quarters in the decompositions
TFP_SHOCK = 0.01  # one-off TFP shock to Z_D, decaying at rho_z = 0.9
# one-off risk shock: p^d jumps from its rest value (cal["pd_rest"]) to this per quarter, then decays. 0.06
# (2026-11-06, author's choice): a crisis-sized shock, so the D banks' constraint binds on impact (it binds from ~5%
# on the 12% box; at the old 2% the D IC is slack along the whole path and there is no recession to repair)
RISK_SHOCK_PD = 0.06
# the top of the risk box, p^d per quarter. 0.12 = the appendix-B box, WARM_START's own: the only box on which the
# D IC binds on impact (at the 6% shock: GDP -0.14%, lending spread +197 bp). 0.062 = Bocola (2016)'s range (his s in
# [-11.41, -2.71]), where it never binds; 0.02 = the illustrative 200 bp cap of 2026-11-05/06
# (Claude files/reports/tpi_root_cause_2026-11-05.md §6-8)
PD_BOX_TOP = 0.12
# the box half-width of the Eurosystem's D-bond book, a share of the D banks' SS book: 0.32 = WARM_START's own band,
# +-27.9% of the D debt, which holds the 6% shock's relief purchase (11.8%). 1.0 (the D banks' whole book) was the
# illustrative 200 bp cap's band on the 2% box; wider than 1.0 the bank's next book turns negative at the M+ edge
# (diagnostics/tpi_root_cause/m_band_edge.py). The band of the Eurosystem's claim keeps its old width
TPI_M_BAND = 0.32
# solve the TPI and compare it with no TPI (on again 2026-11-04, with the cap walk below)
RUN_TPI = True
# how the TPI acts (2026-11-06, author's choice): "relief" = at the shock, unanticipated, the Eurosystem buys just enough
# D bonds to switch the D banks' constraint off and holds them to maturity, on the no-TPI rules (no TPI solve;
# tpi_experiment.intervention). "cap" = the spread cap TPI_CAP solved as an anticipated rule (the walk below); a
# 200 bp cap has no solution at crisis shocks, where the spread left after relief is expected default loss
TPI_MODE = "relief"
JACOBIAN = "complex"  # the collocation Jacobian: "complex" (the complex step, exact) or "fd" (forward differences)
# the solved appendix-B rules (fine grid, 75 points, 12% box) the risk solve starts from; None = the cold-start
# ladder, which does not reach appendix B from the SS. A warm solve re-solves these rules at the current calibration
# on their own grid (on another PD_BOX_TOP, evaluated at the new box's points first), so S_REFINE is not used
WARM_START = ROOT / "data" / "global_warm_start.pkl"
TPI_CAP = 200.0  # TPI spread cap, bp/yr; "rest" = defend the no-TPI rest-point spread (~56 bp)
# the TPI cap walk (recursive_experiment.TPI_WALK): first step and its growth, smallest and largest step (bp), the
# distance each rung keeps from every node's switching cap (bp), the tolerance of a rung above the target, Newton
# iterations per rung, and how often the Jacobian is refreshed in each rung's first pass
TPI_WALK = dict(step=6.0, grow=1.2, min_step=0.5, max_step=100.0, margin=3.0, path_tol=1e-4, maxit=10, reuse=3,
                # a rung whose fast pass ends within rescue x path_tol gets the full Newton
                rescue=100.0,
                # the cap in a few large smoothed rungs (2026-11-05), each started from its nodes above the cap solved
                # alone; then the smoothing of the complementarities walked down from eps_start (one jump to 200 bp
                # stalled at 1.2e-2 on the 1.0 box); False = the cap walk above
                direct=True, eps_start=1e-2, ladder=(270.0, 245.0, 225.0, 210.0), halvings=3)
# each accepted rung is saved here, so a stopped TPI walk resumes where it stopped (None: no checkpoints); delete the
# folder to restart the walk from the top
TPI_CHECKPOINTS = ROOT / "results" / "GLOBAL_tpi_walk"
# the solved TFP and risk rules are stored here, keyed on the calibration and the solve settings, so a run with
# nothing changed reuses them instead of re-solving (~1 h); None = always solve; delete the folder to force a re-solve
SOLVE_CACHE = ROOT / "results" / "GLOBAL_cache"
# what a global run prints first: the steps, their typical times, where everything goes, how to resume
GLOBAL_GUIDE = """
  WHAT THIS RUN DOES (full guide: Claude files/docs/RUN_GLOBAL.md)
    Step 1/7  steady state ................................................ ~1 min
    Step 2/7  TFP-shock solve ............................. ~10 min (seconds if cached)
    Step 3/7  sovereign-risk solve, no TPI .................. ~50 min (seconds if cached)
    Step 4/7  rest point, IRFs, calibration table, decompositions, accuracy ... ~20 min
    Step 5/7  save the no-TPI results ................ results/GLOBAL/data (safe from here on)
    Step 6/7  TPI: TPI_MODE "relief" = the relief purchase at the shock vs no TPI ......... ~10 min
              ("cap": the cap solved as a rule, 1-3 h, saved in results/GLOBAL_tpi_walk and reused)
    Step 7/7  figures ............................................... results/GLOBAL/figures
  The console is logged to results/GLOBAL/run.log. Keep the Mac awake (lid open, on power): a sleeping
  Mac pauses the run. Settings are at the top of run.py (PD_BOX_TOP, RISK_SHOCK_PD, RUN_TPI, TPI_MODE, TPI_CAP, TPI_WALK).
"""


# ============================== THE GLOBAL MODEL ==============================

def compute_global(quick=False):
    """Every GLOBAL computation, step by step; each result is saved to results/GLOBAL/data."""
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES, IS, default_prob
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.recursive_experiment import (
        solve_tfp, tfp_irf, solve_recursive, report_rest_point, impact_table,
        persistence_irf, dynamic_irf, s_from_pd)
    from global_projection.solver_recursive.output_decomposition import (
        simulate, s_decay_path, decompose_output, active_channels, decompose_bond_price,
        exact_bond_reads, BOND_CHANNELS)
    from global_projection.solver_recursive.accuracy import accuracy_report
    from global_projection.solver_recursive.calibrate_stochastic import (
        rest_point_read, rest_point_moments, calibration_values)
    from global_projection.reporting.prints import (banner, print_ss_table,
                                                    print_output_decomposition,
                                                    print_bond_decomposition,
                                                    print_calibration_table)
    from global_projection.reporting.export import export_irfs
    from global_projection.solver_recursive import tpi_experiment
    from global_projection.reporting.prints import print_tpi_report

    data, _, log = results_io.model_dirs(RESULTS / "GLOBAL", fresh=True)
    settings = dict(quick=quick, S_REFINE=0 if quick else S_REFINE,
                    NW_FLOOR=NW_FLOOR, MU=MU,
                    RISK_MU_VEC=RISK_MU_VEC, ROTATE_P=ROTATE_P, ACCURACY_T=ACCURACY_T,
                    DECOMP_T=DECOMP_T, TFP_SHOCK=TFP_SHOCK, RISK_SHOCK_PD=RISK_SHOCK_PD,
                    RUN_TPI=RUN_TPI, TPI_MODE=TPI_MODE, TPI_CAP=TPI_CAP, JACOBIAN=JACOBIAN, PD_BOX_TOP=PD_BOX_TOP,
                    TPI_M_BAND=TPI_M_BAND,
                    TPI_WALK=TPI_WALK, TPI_CHECKPOINTS=str(TPI_CHECKPOINTS) if TPI_CHECKPOINTS else None,
                    WARM_START=str(WARM_START) if WARM_START and WARM_START.exists() else None)
    t0 = time.perf_counter()
    _keep_awake()
    with _logged(log, "w"):
        print(GLOBAL_GUIDE)
        # 1. calibration and the deterministic steady state
        banner("Step 1/7 — steady state (two-country HANK-GK monetary union)")
        cal = get_calibration()
        cal["nw_floor_frac"] = NW_FLOOR
        cal["jacobian"] = JACOBIAN
        cal["pd_box_top"] = PD_BOX_TOP
        RE.TPI_WALK.update(TPI_WALK)
        RE.TPI_CHECKPOINTS = TPI_CHECKPOINTS
        # o_band keeps the claim's half-width at 0.25 x 0.32 of the D banks' book whatever the book's band
        RE.BOX_KW.update(m_band=TPI_M_BAND, o_band=0.08 / TPI_M_BAND)
        ss = solve_steady_state(cal)
        sproc = s_process_params(cal)
        calibrate_household_anchors(cal, ss, sproc)
        print(f"  solved in {time.perf_counter() - t0:.0f}s")
        print_ss_table(ss, cal)

        # the solved rules both global solves start from (WARM_START), if any
        warm = None
        if settings["WARM_START"]:
            with open(settings["WARM_START"], "rb") as fh:
                warm = pickle.load(fh)
            print(f"  warm start: {settings['WARM_START']} ({warm.grid.n} points)")
            grid = RE.normal_box(ss, cal, mu=1, refine=(IS, S_REFINE, "cross"))
            if abs(grid.hi[IS] - warm.grid.hi[IS]) > 1e-9:
                # solved on another risk box: the rules evaluated at this box's points seed the solve
                print(f"  the warm rules' box tops at p^d {100 * default_prob(warm.grid.hi[IS]):.2f}%; seeded onto "
                      f"PD_BOX_TOP {100 * PD_BOX_TOP:.2f}% ({grid.n} points)")
                warm = RE._seed_from(RuleSet(grid, 2), warm)

        # 2. TFP shock
        key = _cache_key(cal, settings)
        _step_done(1, t0)
        banner("Step 2/7 — TFP shock, global collocation (Z_D as the TFP state)")
        t1 = time.perf_counter()
        rules_tfp = _cached("tfp", key, lambda: solve_tfp(cal, ss, sproc, mu=MU, seed=warm))
        if not rules_tfp.solve_ok:
            print(f"  NOTE: the TFP solve did not converge (max|F| {rules_tfp.solve_worst:.1e}); its IRF is not usable")
        tfp = tfp_irf(rules_tfp, cal, ss, sproc, dz=TFP_SHOCK, exact=True)

        # 3. sovereign-risk pass-through: solve, find the rest point, read the IRFs
        _step_done(2, t1, t0)
        banner("Step 3/7 — sovereign-risk pass-through, global collocation (12 states, Newton)")
        t1 = time.perf_counter()
        base = []  # the coarse no-TPI baseline, reused by the TPI
        rules = _cached("risk", key, lambda: solve_recursive(cal, ss, sproc, mu_vec=RISK_MU_VEC, rotate=ROTATE_P,
                                                             s_refine=0 if warm else settings["S_REFINE"],
                                                             base=warm, base_out=base))
        if not base:
            base.append(rules)
        _step_done(3, t1, t0)
        banner("Step 4/7 — rest point, IRFs, calibration table, decompositions, accuracy")
        t1 = time.perf_counter()
        # the TPI starts from the coarse baseline, or, warm, from the solved rules on their own grid
        tpi_base, tpi_refine = (rules, 0) if warm else (base[0], settings["S_REFINE"])
        # every readout clears the period map exactly at each visited state (2026-11-06, author's decision): at the 6%
        # crisis shock the fitted rules' residual (4-6e-3) doubled the impact recession against the exact read
        S_rest = report_rest_point(rules, cal, ss, sproc)  # every IRF starts here
        impact_table(rules, cal, ss, sproc)
        persistence = persistence_irf(rules, cal, ss, sproc, pd_shock=RISK_SHOCK_PD, exact=True)
        risk = dynamic_irf(rules, cal, ss, sproc, pd_shock=RISK_SHOCK_PD, T=25, exact=True)

        # 3b. the calibration table: every target and validation moment at the rest point
        banner("Calibration table — data against the stochastic rest point")
        _, o_rest = rest_point_read(rules, cal, ss, sproc, exact=True)
        srp_moments = rest_point_moments(o_rest)
        # its transmission rows are the p^d 2%/qtr impact (appendix B's 108 bp row), whatever the headline shock
        irf_2 = (risk if abs(RISK_SHOCK_PD - 0.02) < 1e-12 else
                 dynamic_irf(rules, cal, ss, sproc, pd_shock=0.02, T=25, rest_verbose=False, exact=True))
        calib = calibration_values(o_rest, cal, sproc, irf_2)
        print_calibration_table(calib, srp_moments, cal)

        # 4. decompositions of the same shock
        banner("Decompositions — which channels produce the response")
        s_path = s_decay_path(sproc, s_from_pd(RISK_SHOCK_PD), DECOMP_T)
        sim = simulate(rules, cal, ss, sproc, s_path, S_init=S_rest, exact=True)
        ref = simulate(rules, cal, ss, sproc, [sproc["s_star"]] * DECOMP_T, S_init=S_rest, exact=True)
        output_dec = decompose_output(sim, ref, cal)
        print_output_decomposition(output_dec, active_channels(output_dec))
        # Bocola's split of the bond's excess return, on period maps cleared exactly along both paths
        bond_dec = decompose_bond_price(exact_bond_reads(rules, cal, ss, sproc, sim),
                                        exact_bond_reads(rules, cal, ss, sproc, ref), cal)
        print_bond_decomposition(bond_dec, BOND_CHANNELS)

        # 5. solution accuracy on the ergodic set
        banner("Solution accuracy — Euler errors on the ergodic set")
        accuracy = accuracy_report(rules, cal, ss, sproc, S_rest, T=ACCURACY_T,
                                   label="sovereign-risk rules")

        # 6. save every no-TPI result before the long TPI solve, with the TPI's cap already set
        cal["tpi_cap_bp"] = (tpi_experiment.rest_point_row(rules, cal, ss, sproc)["sov_bp"]
                             if TPI_CAP == "rest" else float(TPI_CAP))
        _step_done(4, t1, t0)
        banner(f"Step 5/7 — saving the no-TPI results -> {data}")
        save(data / "settings", settings)
        save(data / "calibration", cal)
        save(data / "steady_state", ss)
        save(data / "tfp_irf", tfp)
        save(data / "rest_point", dict(zip(STATE_NAMES, S_rest)))
        save(data / "persistence_irf", persistence)
        save(data / "risk_irf", risk)
        save(data / "calibration_table", dict(values=calib, moments=srp_moments))
        save(data / "output_decomposition", output_dec)
        save(data / "bond_decomposition", bond_dec)
        save(data / "accuracy", accuracy)
        # the solved rules, and the coarse baseline the TPI starts from
        for name, r in (("rules_tfp", rules_tfp), ("rules_risk", rules),
                        ("rules_base_coarse", tpi_base)):
            with open(data / f"{name}.pkl", "wb") as fh:
                pickle.dump(r, fh)
        export_irfs(risk, tfp, sproc, RISK_SHOCK_PD, TFP_SHOCK, data / "comparison_irfs.json",
                    grid_note=" [--quick: coarse grid]" if quick else "")

        # 7. the TPI: the same economy with the spread cap switched on
        print("  step 5/7 done: every no-TPI result is saved; an interrupted TPI walk below loses none of it")
        if RUN_TPI and TPI_MODE == "relief":
            t1 = time.perf_counter()
            banner(f"Step 6/7 — TPI relief intervention at the {100 * RISK_SHOCK_PD:.0f}% shock "
                   "(unanticipated, held to maturity, no-TPI rules)")
            tpi = tpi_experiment.intervention(cal, ss, sproc, rules, pd_shock=RISK_SHOCK_PD)
            print_tpi_report(tpi, BOND_CHANNELS)
            save(data / "tpi", tpi)
            _step_done(6, t1, t0)
        elif RUN_TPI:
            t1 = time.perf_counter()
            banner(f"Step 6/7 — TPI backstop, cap walk to {cal['tpi_cap_bp']:.1f} bp/yr over the F bond"
                   + (" (the no-TPI rest-point spread)" if TPI_CAP == "rest" else ""))
            tpi, rules_tpi = tpi_experiment.run(cal, ss, sproc, rules, tpi_base,
                                                s_refine=tpi_refine,
                                                mu_vec=RISK_MU_VEC, pd_shock=RISK_SHOCK_PD)
            print_tpi_report(tpi, BOND_CHANNELS)
            save(data / "tpi", tpi)
            with open(data / "rules_tpi.pkl", "wb") as fh:
                pickle.dump(rules_tpi, fh)
            _step_done(6, t1, t0)
        print(f"  {len(list(data.iterdir()))} files written")
        print(f"\nTOTAL  {time.perf_counter() - t0:.0f}s")


def plot_global():
    """Every GLOBAL figure, drawn from results/GLOBAL/data into results/GLOBAL/figures."""
    from global_projection.reporting import plots
    from global_projection.reporting.export import export_calibration_tex
    from global_projection.reporting.prints import banner
    from global_projection.solver_recursive.output_decomposition import (active_channels,
                                                                         BOND_CHANNELS)

    data, figures, log = results_io.model_dirs(RESULTS / "GLOBAL")
    plots.OUTDIR = str(figures)
    with _logged(log, "a"):
        banner("Step 7/7 — figures, drawn from the saved results")
        s = load(data / "settings")
        shock = f"p-d shock to {100 * s['RISK_SHOCK_PD']:.2f}%/qtr"
        output_dec = load(data / "output_decomposition")
        paths = [
            plots.plot_tfp_irf(load(data / "tfp_irf"),
                               note=f"{s['TFP_SHOCK']:.0%} shock to Z_D"),
            plots.plot_risk_irf(load(data / "risk_irf"), note=f"{shock}, states evolving"),
            plots.plot_output_decomposition(
                [(shock, output_dec)], active_channels(output_dec),
                note="d log Y_D split by the production function and the GHH labour FOC; "
                     "the two wedge legs are a symmetric (Shapley) split, so the channels "
                     "sum to the total identically"),
            plots.plot_bond_decomposition(load(data / "bond_decomposition"),
                                          list(BOND_CHANNELS),
                                          note="Bocola (2016): liquidity premium $\\lambda\\mu/E[\\Omega]$ "
                                               "and risk premium $-Cov(\\Omega, R_b)/E[\\Omega]$"),
        ]
        if s.get("RUN_TPI"):
            tpi = load(data / "tpi")
            paths += [plots.plot_tpi_irf(tpi, note=shock),
                      plots.plot_tpi_mechanism(tpi, list(BOND_CHANNELS))]
        # the paper's calibration table, filled from the saved rest-point moments
        paths.append(export_calibration_tex(load(data / "calibration_table")["values"],
                                            str(figures / "calibration_table.tex")))
        for p in paths:
            print(f"  figure -> {p}")


# ================================ THE SSJ MODEL ================================

def compute_ssj(quick=False):
    """Every SSJ computation, saved to results/SSJ/data (subprocess under SSJ_PYTHON).

    quick has no SSJ counterpart: the whole linear pipeline takes ~15 min.
    """
    _, _, log = results_io.model_dirs(RESULTS / "SSJ")
    _run_ssj_stage("compute", log, "w")


def plot_ssj():
    """Every SSJ figure, drawn from results/SSJ/data into results/SSJ/figures."""
    _, _, log = results_io.model_dirs(RESULTS / "SSJ")
    _run_ssj_stage("plot", log, "a")


def _run_ssj_stage(stage, log, mode):
    if shutil.which(SSJ_PYTHON) is None:
        sys.exit(f"SSJ interpreter not found: {SSJ_PYTHON}\n"
                 f"  set SSJ_PYTHON in run.py (or the environment), e.g. SSJ_PYTHON=python3")
    cmd = [SSJ_PYTHON, "-m", "linear_ssj.main", "--stage", stage,
           "--results", str(RESULTS / "SSJ")]
    print(f"\n>>> {' '.join(cmd)}\n", flush=True)
    env = dict(os.environ, PYTHONUNBUFFERED="1")
    with open(log, mode) as fh:
        proc = subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, env=env)
        for line in proc.stdout:  # to the terminal and the run log
            sys.stdout.write(line)
            fh.write(line)
        if proc.wait():
            sys.exit(proc.returncode)


# =================================== RUNNER ===================================

class _Tee:
    # a stdout that writes to the terminal and a log file at once
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for s in self.streams:
            s.write(text)

    def flush(self):
        for s in self.streams:
            s.flush()


def _cache_key(cal, settings):
    """A short hash of everything the TFP and risk solves depend on: the calibration and the solve settings."""
    import hashlib
    import json
    scalars = {k: (round(v, 12) if isinstance(v, float) else v) for k, v in cal.items()
               if isinstance(v, (bool, int, float, str))}
    solve = {k: settings.get(k) for k in ("quick", "S_REFINE", "NW_FLOOR", "MU", "RISK_MU_VEC", "ROTATE_P",
                                          "JACOBIAN", "PD_BOX_TOP", "TPI_M_BAND", "WARM_START")}
    return hashlib.md5(json.dumps(dict(cal=scalars, solve=solve), sort_keys=True, default=str).encode()).hexdigest()[:12]


def _cached(name, key, solve):
    """The rules `name` from SOLVE_CACHE when solved before with this key; otherwise solve() and store them."""
    path = None if SOLVE_CACHE is None else Path(SOLVE_CACHE) / f"{name}_{key}.pkl"
    if path is not None and path.exists():
        with open(path, "rb") as fh:
            rules = pickle.load(fh)
        print(f"  reusing the {name} rules solved before with these settings ({path}); delete "
              f"{SOLVE_CACHE} to re-solve")
        return rules
    rules = solve()
    if path is not None and getattr(rules, "solve_ok", False):
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as fh:
            pickle.dump(rules, fh)
    return rules


def _step_done(k, t_step, t0=None):
    """One line after each step of the global run: its time and the time so far."""
    now = time.perf_counter()
    total = "" if t0 is None else f", {(now - t0) / 60:.0f} min since the start"
    print(f"\n  step {k}/7 done in {(now - t_step) / 60:.1f} min{total}")


def _keep_awake():
    """On a Mac, stop idle sleep while this process runs (closing the lid still sleeps it)."""
    if sys.platform == "darwin" and shutil.which("caffeinate"):
        subprocess.Popen(["caffeinate", "-ims", "-w", str(os.getpid())])


@contextlib.contextmanager
def _logged(log, mode):
    with open(log, mode) as fh, contextlib.redirect_stdout(_Tee(sys.stdout, fh)):
        yield


COMPUTE = {"SSJ": compute_ssj, "GLOBAL": compute_global}
PLOT = {"SSJ": plot_ssj, "GLOBAL": plot_global}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("model", nargs="?", default=MODEL,
                   choices=("ssj", "global", "both", "compare"))
    p.add_argument("--quick", action="store_true",
                   help="global model: coarse grid (fast, not converged)")
    p.add_argument("--plots-only", action="store_true",
                   help="skip computing; redraw the figures from the saved data")
    a = p.parse_args()
    models = {"ssj": ["SSJ"], "global": ["GLOBAL"], "both": ["SSJ", "GLOBAL"],
              "compare": []}[a.model]

    if not a.plots_only:  # 1. compute every selected model
        for m in models:
            COMPUTE[m](quick=a.quick or QUICK)
    for m in models:  # 2. then draw every figure from the saved data
        PLOT[m]()
    if a.model in ("both", "compare"):  # 3. compare the two solutions
        import compare
        compare.main()


if __name__ == "__main__":
    main()
