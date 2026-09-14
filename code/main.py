"""
Two-Country MU HANK — main orchestrator.

Run from the repo root with the ssj conda environment:
    /opt/anaconda3/envs/ssj/bin/python code/main.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from calibration              import get_calibration
from steady_state             import solve_steady_state
from ic_delta_calibration     import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model               import build_and_solve
from irf_plots                import generate_irf_plots
from tpi                      import run_tpi
from tpi_plots                import generate_tpi_plots

OUTPUT_DIR = Path(__file__).parent.parent / 'outputs'


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}\n")

    print("=" * 60)
    print("Step 1: Calibration")
    print("=" * 60)
    calibration_start = get_calibration()
    print(f"  {len(calibration_start)} parameters loaded.\n")

    print("=" * 60)
    print("Step 2: Steady State — initial solve + portfolio targeting")
    print("=" * 60)
    ss_results = solve_steady_state(calibration_start)
    print()

    print("=" * 60)
    print("Step 3: IC Delta Calibration (back-solve divertable fraction)")
    print("=" * 60)
    ss_results = calibrate_ic_delta(ss_results)
    print()

    print("=" * 60)
    print("Step 4: Depreciation Calibration + Final Steady-State Re-solve")
    print("=" * 60)
    ss_results = calibrate_depreciation(ss_results)
    print()

    print("=" * 60)
    print("Step 5: Full Dynamic Model + Baseline IRFs")
    print("=" * 60)
    model_results = build_and_solve(ss_results)
    print()

    print("=" * 60)
    print("Step 6: Baseline IRF Plots")
    print("=" * 60)
    generate_irf_plots(model_results, OUTPUT_DIR)
    print()

    print("=" * 60)
    print("Step 7: TPI Experiment (Jacobian + closed-loop IRFs)")
    print("=" * 60)
    tpi_results = run_tpi(model_results)
    print()

    print("=" * 60)
    print("Step 8: TPI Plots")
    print("=" * 60)
    generate_tpi_plots(tpi_results, OUTPUT_DIR)
    print()

    print("=" * 60)
    print(f"Done — all figures saved to: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    main()
