# ENTRY POINT: CHEBYSHEV-SMOLYAK PROJECTION SOLVER FOR THE WHOLE MODEL.
# The two-country model is solved GLOBALLY as recursive decision rules on a Smolyak
# grid (solver_recursive/) -- there is NO perfect-foresight / representative-branch
# machinery. main() solves the steady state once, then runs three projection
# experiments, each a full time-iteration solve: a TFP shock (Z_D as the 7th state),
# the Bocola sovereign-risk pass-through, and the OMT/TPI activation comparison
# (0/50/100%). Computationally heavy by design (~20-30 min end to end).
import time

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import s_process_params
from solver_recursive.recursive_main import calibrate_household_anchors
from solver_recursive.recursive_experiment import (
    solve_tfp, tfp_irf, solve_recursive, impact_table, persistence_irf)
from solver_recursive import tpi_recursive_experiment
from reporting.prints import banner, print_ss_table


def main():
    t0 = time.perf_counter()
    cal = get_calibration()

    banner("Two-country HANK-GK monetary union: steady state")
    ss = solve_steady_state(cal)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_ss_table(ss, cal)

    banner("TFP shock — recursive projection (Z_D as the 7th state)")
    rules_tfp = solve_tfp(cal, ss, sproc)
    tfp_irf(rules_tfp, cal, ss, sproc)

    banner("Sovereign-risk pass-through — recursive projection")
    rules_risk = solve_recursive(cal, ss, sproc)
    impact_table(rules_risk, cal, ss, sproc)
    persistence_irf(rules_risk, cal, ss, sproc)

    banner("OMT/TPI activation — recursive projection (0/50/100%)")
    tpi_recursive_experiment.run(cal, ss, sproc)

    print(f"\nTOTAL  {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
