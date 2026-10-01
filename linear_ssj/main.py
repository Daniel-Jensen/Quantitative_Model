"""
Two-Country MU HANK — the sequence-space (SSJ) pipeline, in two stages.

    /opt/anaconda3/envs/ssj/bin/python -m linear_ssj.main [--stage compute|plot|all] [--results DIR]
or, with the interpreter picked for you:  python3 run.py ssj

compute  calibration -> steady state -> IC-delta -> depreciation -> Jacobian + IRFs -> TPI,
         then every result is saved as plain data to DIR/data
plot     every figure, redrawn from DIR/data into DIR/figures
DIR defaults to results/SSJ.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root

import results_io
from calibration.ssj                            import get_calibration
from linear_ssj.solve.steady_state             import solve_steady_state
from linear_ssj.solve.ic_delta_calibration     import calibrate_ic_delta
from linear_ssj.solve.depreciation_calibration import calibrate_depreciation
from linear_ssj.model.full_model               import build_and_solve
from linear_ssj.solve.tpi                      import run_tpi
from linear_ssj.reporting.export               import export_irfs
from linear_ssj.reporting.irf_plots            import generate_irf_plots
from linear_ssj.reporting.tpi_plots            import generate_tpi_plots

DEFAULT_RESULTS = results_io.RESULTS / 'SSJ'
# not saved: the TPI Jacobian (large, recomputable, no figure reads it) and the steady
# state, which is saved on its own as data/steady_state
TPI_NOT_SAVED = ('G_tpi', 'ss_final')


def _step(n, title):
    print("=" * 60)
    print(f"Step {n}: {title}")
    print("=" * 60)


def compute(results):
    data, _, _ = results_io.model_dirs(results, fresh=True)
    print(f"Results directory: {results}\n")

    _step(1, "Calibration")
    calibration_start = get_calibration()
    print(f"  {len(calibration_start)} parameters loaded.\n")

    _step(2, "Steady State — initial solve + portfolio targeting")
    ss_results = solve_steady_state(calibration_start)
    print()

    _step(3, "IC Delta Calibration (back-solve divertable fraction)")
    ss_results = calibrate_ic_delta(ss_results)
    print()

    _step(4, "Depreciation Calibration + Final Steady-State Re-solve")
    ss_results = calibrate_depreciation(ss_results)
    print()

    _step(5, "Full Dynamic Model + Baseline IRFs")
    model_results = build_and_solve(ss_results)
    print()

    _step(6, "TPI Experiment (Jacobian + closed-loop IRFs)")
    tpi_results = run_tpi(model_results)
    print()

    _step(7, f"Saving results -> {data}")
    results_io.save(data / 'steady_state', model_results['ss_final'])
    results_io.save(data / 'irfs', {k: model_results[k] for k in
                                    ('irfs_Z_D', 'irfs_def_D', 'dShock_def_D', 'T')})
    results_io.save(data / 'tpi', {k: v for k, v in tpi_results.items()
                                   if k not in TPI_NOT_SAVED})
    export_irfs(model_results, data / 'comparison_irfs.json')
    print(f"  {len(list(data.iterdir()))} files written\n")


def plot(results):
    data, figures, _ = results_io.model_dirs(results)

    _step(8, f"Baseline IRF Plots (from {data})")
    generate_irf_plots(results_io.load(data / 'irfs'), figures)
    print()

    _step(9, "TPI Plots")
    tpi_results = results_io.load(data / 'tpi')
    tpi_results['ss_final'] = results_io.load(data / 'steady_state')
    generate_tpi_plots(tpi_results, figures)
    print()

    print("=" * 60)
    print(f"Done — all figures saved to: {figures}")
    print("=" * 60)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--stage', choices=('compute', 'plot', 'all'), default='all')
    p.add_argument('--results', type=Path, default=DEFAULT_RESULTS)
    a = p.parse_args()
    if a.stage in ('compute', 'all'):
        compute(a.results)
    if a.stage in ('plot', 'all'):
        plot(a.results)


if __name__ == '__main__':
    main()
