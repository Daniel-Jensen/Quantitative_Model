"""Unified entry point: the linear (SSJ) model, the global model, or both, then a comparison.

    python3 run.py                    run MODEL (set below)
    python3 run.py ssj|global|both    override MODEL
    python3 run.py global --quick     coarse grid: a preview, NOT converged
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

os.environ.setdefault("MPLBACKEND", "Agg")   # figures are files; never open a window

import results_io
from results_io import RESULTS, load, save

ROOT = Path(__file__).resolve().parent

# ================================ CONFIGURATION ================================
MODEL = "both"            # "ssj", "global" or "both" -- the command-line argument overrides it
# the SSJ model's interpreter: the ssj conda env (the base env has a broken liblapack
# symlink, CLAUDE.md); the SSJ_PYTHON environment variable overrides it
SSJ_PYTHON = os.environ.get("SSJ_PYTHON", "/opt/anaconda3/envs/ssj/bin/python")

# GLOBAL model settings -- how it is solved and which experiments run. Its economic
# parameters live in calibration/global_projection.py. Measured cost: one dense
# finite-difference Jacobian per Newton step, so a solve scales as (points x regimes)^2;
# ~100 min end to end at S_REFINE = 5.
NW_FLOOR = 0.15           # Bocola's net-worth floor (fraction of n_ss): keeps the deep default corners feasible
MU = 1                    # Smolyak level of the TFP grid (no risk dimension to resolve there)
RISK_MU_VEC = None        # per-state Smolyak levels for the coarse risk grid; None = isotropic mu = 1
S_REFINE = 5              # dense Chebyshev nodes in s for the risk solve: 5 = 95 points (converged for
                          # every reported object), 9 = Bocola's resolution (~4 h), 0 = coarse only
ROTATE_P = False          # P-block eigenbasis box: right in theory, measures worse (see solve_recursive)
ACCURACY_T = 1200         # simulated periods for the Euler-error report
DECOMP_T = 25             # quarters in the decompositions
TFP_SHOCK = 0.01          # one-off TFP shock: +1% to Z_D, decaying at rho_z = 0.9
RISK_SHOCK_PD = 0.0198    # one-off risk shock: p^d jumps 0.10% -> 1.98%/qtr, decays at rho_s
RUN_TPI = True            # the TPI backstop against no TPI, at calibration tpi_cap_bp: one more
                          # solve (the cap homotopy on the coarse grid, then the refined grid)


# ============================== THE GLOBAL MODEL ==============================
# Solved GLOBALLY as recursive decision rules on a Smolyak grid by one collocation
# Newton (global_projection/solver_recursive). Read the two spreads correctly: the BANK
# CREDIT spread lambda_K*mu/alpha is the wedge the leverage constraint puts on capital and,
# through r_wc = rdep + lambda*mu/E[Om], the only channel into output under GHH; the
# SOVEREIGN spread y_D - y_F comes out of the bond Euler and is what the figures plot.

def compute_global(quick=False):
    """Every GLOBAL computation, step by step; each result is saved to results/GLOBAL/data."""
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.recursive_experiment import (
        solve_tfp, tfp_irf, solve_recursive, report_rest_point, impact_table,
        persistence_irf, dynamic_irf, s_from_pd)
    from global_projection.solver_recursive.output_decomposition import (
        simulate, s_decay_path, decompose_output, active_channels, decompose_bond_price,
        BOND_CHANNELS)
    from global_projection.solver_recursive.accuracy import accuracy_report
    from global_projection.reporting.prints import (banner, print_ss_table,
                                                    print_output_decomposition,
                                                    print_bond_decomposition)
    from global_projection.reporting.export import export_irfs
    from global_projection.solver_recursive import tpi_experiment
    from global_projection.reporting.prints import print_tpi_report

    data, _, log = results_io.model_dirs(RESULTS / "GLOBAL", fresh=True)
    settings = dict(quick=quick, S_REFINE=0 if quick else S_REFINE,
                    NW_FLOOR=NW_FLOOR, MU=MU,
                    RISK_MU_VEC=RISK_MU_VEC, ROTATE_P=ROTATE_P, ACCURACY_T=ACCURACY_T,
                    DECOMP_T=DECOMP_T, TFP_SHOCK=TFP_SHOCK, RISK_SHOCK_PD=RISK_SHOCK_PD,
                    RUN_TPI=RUN_TPI)
    t0 = time.perf_counter()
    with _logged(log, "w"):
        # 1. calibration and the deterministic steady state
        banner("Two-country HANK-GK monetary union: steady state")
        cal = get_calibration()
        cal["nw_floor_frac"] = NW_FLOOR
        ss = solve_steady_state(cal)
        sproc = s_process_params(cal)
        calibrate_household_anchors(cal, ss, sproc)
        print(f"  solved in {time.perf_counter() - t0:.0f}s")
        print_ss_table(ss, cal)

        # 2. TFP shock (solved at pi == 0, where the deterministic SS IS the rest point)
        banner("TFP shock — global collocation (Z_D as the TFP state)")
        rules_tfp = solve_tfp(cal, ss, sproc, mu=MU)
        tfp = tfp_irf(rules_tfp, cal, ss, sproc, dz=TFP_SHOCK)

        # 3. sovereign-risk pass-through: solve, find the model's own rest point, read IRFs
        banner("Sovereign-risk pass-through — global collocation (12-state, Newton)")
        base = []                                          # the coarse no-TPI baseline,
        rules = solve_recursive(cal, ss, sproc, mu_vec=RISK_MU_VEC, rotate=ROTATE_P,
                                s_refine=settings["S_REFINE"], base_out=base)  # reused by the TPI
        S_rest = report_rest_point(rules, cal, ss, sproc)   # every IRF is read against this
        impact_table(rules, cal, ss, sproc)
        persistence = persistence_irf(rules, cal, ss, sproc, pd_shock=RISK_SHOCK_PD)
        risk = dynamic_irf(rules, cal, ss, sproc, pd_shock=RISK_SHOCK_PD, T=25)

        # 4. decompositions of the same shock, off the same rules, around the rest point
        banner("Decompositions — which channels produce the response")
        s_path = s_decay_path(sproc, s_from_pd(RISK_SHOCK_PD), DECOMP_T)
        sim = simulate(rules, cal, ss, sproc, s_path, S_init=S_rest)
        ref = simulate(rules, cal, ss, sproc, [sproc["s_star"]] * DECOMP_T, S_init=S_rest)
        output_dec = decompose_output(sim, ref, cal)
        print_output_decomposition(output_dec, active_channels(output_dec))
        bond_dec = decompose_bond_price(sim, ref, cal)
        print_bond_decomposition(bond_dec, BOND_CHANNELS)

        # 5. solution accuracy on the ergodic set, simulated from the rest point
        banner("Solution accuracy — Euler errors on the ergodic set")
        accuracy = accuracy_report(rules, cal, ss, sproc, S_rest, T=ACCURACY_T,
                                   label="sovereign-risk rules")

        # 6. the TPI backstop: the same economy with the Eurosystem's spread cap switched on
        tpi = rules_tpi = None
        if RUN_TPI:
            banner(f"TPI backstop — spread cap {cal['tpi_cap_bp']:.0f} bp/yr over the F bond")
            tpi, rules_tpi = tpi_experiment.run(cal, ss, sproc, rules, base[0],
                                                s_refine=settings["S_REFINE"],
                                                mu_vec=RISK_MU_VEC, pd_shock=RISK_SHOCK_PD)
            print_tpi_report(tpi, BOND_CHANNELS)

        # 7. save every computed object
        banner(f"Saving -> {data}")
        save(data / "settings", settings)
        save(data / "calibration", cal)
        save(data / "steady_state", ss)
        save(data / "tfp_irf", tfp)
        save(data / "rest_point", dict(zip(STATE_NAMES, S_rest)))
        save(data / "persistence_irf", persistence)
        save(data / "risk_irf", risk)
        save(data / "output_decomposition", output_dec)
        save(data / "bond_decomposition", bond_dec)
        save(data / "accuracy", accuracy)
        if tpi is not None:
            save(data / "tpi", tpi)
        for name, r in (("rules_tfp", rules_tfp), ("rules_risk", rules), ("rules_tpi", rules_tpi)):
            if r is not None:
                with open(data / f"{name}.pkl", "wb") as fh:   # the solved decision rules
                    pickle.dump(r, fh)
        export_irfs(risk, tfp, sproc, RISK_SHOCK_PD, TFP_SHOCK, data / "comparison_irfs.json",
                    grid_note=" [--quick: coarse grid]" if quick else "")
        print(f"  {len(list(data.iterdir()))} files written")
        print(f"\nTOTAL  {time.perf_counter() - t0:.0f}s")


def plot_global():
    """Every GLOBAL figure, drawn from results/GLOBAL/data into results/GLOBAL/figures."""
    from global_projection.reporting import plots
    from global_projection.reporting.prints import banner
    from global_projection.solver_recursive.output_decomposition import (active_channels,
                                                                         BOND_CHANNELS)

    data, figures, log = results_io.model_dirs(RESULTS / "GLOBAL")
    plots.OUTDIR = str(figures)
    with _logged(log, "a"):
        banner("Figures — drawn from the saved results")
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
                                          note="the D bank first-order condition, split leg by leg"),
        ]
        if s.get("RUN_TPI"):
            paths.append(plots.plot_tpi_irf(load(data / "tpi"), note=shock))
        for p in paths:
            print(f"  figure -> {p}")


# ================================ THE SSJ MODEL ================================
# Linearised in sequence space (sequence_jacobian): calibration -> steady state -> IC-delta
# and depreciation calibration -> Jacobian + IRFs -> TPI. The steps live in
# linear_ssj/main.py, run under SSJ_PYTHON; its stages mirror the global ones.

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
        for line in proc.stdout:          # to the terminal AND the run log
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

    if not a.plots_only:                      # 1. compute every selected model
        for m in models:
            COMPUTE[m](quick=a.quick)
    for m in models:                          # 2. then draw every figure from saved data
        PLOT[m]()
    if a.model in ("both", "compare"):        # 3. compare the two solutions
        import compare
        compare.main()


if __name__ == "__main__":
    main()
