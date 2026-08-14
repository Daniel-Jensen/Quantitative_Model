# ENTRY POINT: CHEBYSHEV-SMOLYAK PROJECTION SOLVER FOR THE WHOLE MODEL.
# The two-country model is solved GLOBALLY as recursive decision rules on a Smolyak
# grid (solver_recursive/) -- there is NO perfect-foresight / representative-branch
# machinery. main() solves the steady state once, then runs the projection
# experiments at MU=2, each a full time-iteration solve: a TFP shock (Z_D as the 7th
# state) and the Bocola sovereign-risk pass-through, plus the OMT/TPI activation
# overlay (RUN_TPI, at MU_TPI over the selectable TPI_ACTIVATIONS probabilities).
# Computationally heavy by design (mu=2, ~2-3h end to end).
import time

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import s_process_params
from solver_recursive.recursive_main import calibrate_household_anchors
from solver_recursive.recursive_experiment import (
    solve_tfp, tfp_irf, solve_recursive, impact_table, persistence_irf)
from solver_recursive import tpi_recursive_experiment
from reporting.prints import banner, print_ss_table

# MU=1 IS THE PRODUCTION GRID, at the Bocola-faithful calibration (leverage 4, 200bp,
# mu_ss~0.02; see config/calibration.py). Bocola's own posterior mu^bg~0.001 (~8bp) sits
# ON the constraint kink, which our pointwise solver cannot cross; GK11's 200bp is the
# robustly-binding proxy. mu=2 does NOT converge at this realistic spread -- the barely-
# binding constraint lets the mu=2 grid slip to the slack/expansionary branch. Forcing
# mu=2 by "binding harder" only worked by pushing the SS spread to an implausible 720bp
# (and even then the d=1 fit stayed non-convergent) -- a regression, now removed. The
# real mu=2 fix is Bocola's structural device: fit the franchise value alpha as a POLICY
# via its forward Euler (residual_nodefault.m), NOT our pointwise least-root fixed point,
# so the franchise fold never arises. That is a solver change, not a calibration knob;
# see the Bocola-2016 comparison report. MU_TPI matches MU.
MU = 1
RUN_TPI = False                         # append the OMT/TPI activation overlay
TPI_ACTIVATIONS = (0.00, 0.50, 1.00)    # priced activation probabilities to compare
MU_TPI = 1


def main():
    t0 = time.perf_counter()
    cal = get_calibration()   # standard Bocola-faithful SS: leverage 4, 200bp, mu_ss~0.02
                              # (no harder-binding override; the 720bp detour is removed)

    banner("Two-country HANK-GK monetary union: steady state")
    ss = solve_steady_state(cal)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_ss_table(ss, cal)

    banner("TFP shock — recursive projection (Z_D as the 7th state)")
    rules_tfp = solve_tfp(cal, ss, sproc, mu=MU)
    tfp_irf(rules_tfp, cal, ss, sproc)

    banner("Sovereign-risk pass-through — recursive projection")
    rules_risk = solve_recursive(cal, ss, sproc, mu=MU)
    impact_table(rules_risk, cal, ss, sproc)
    persistence_irf(rules_risk, cal, ss, sproc)

    if RUN_TPI:
        pct = "/".join(f"{round(a * 100):.0f}%" for a in TPI_ACTIVATIONS)
        banner(f"OMT/TPI activation — recursive projection ({pct})")
        tpi_recursive_experiment.run(cal, ss, sproc, mu=MU_TPI, activations=TPI_ACTIVATIONS)

    print(f"\nTOTAL  {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
