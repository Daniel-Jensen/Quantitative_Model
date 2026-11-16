# PROBE (2026-10-26): THE d0 RESIDUAL AT THE SS SEED, THEN AFTER EACH WARM TIME-ITERATION SWEEP, NODE BY NODE.
import os, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
import numpy as np
sys.path.insert(0, str(HERE))
import diag_d0 as DD


def main():
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors, time_iteration
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive import recursive_experiment as RE
    DD.STATE = STATE_NAMES
    cal = get_calibration(); cal["nw_floor_frac"] = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal)
    RE._soften(cal); cal["mu_eps"] = RE.MU_EPS_START
    cal["wc_flow_D"] = cal["wc_flow_F"] = 0.0; cal["fx_val"] = 0.0
    calibrate_household_anchors(cal, ss, sp)
    rules = RuleSet.from_ss(RE.normal_box(ss, cal, mu=1), ss, cal, n_regimes=2); rules.n_gh = RE.N_GH
    DD.worst_residuals(rules, cal, ss, sp, "SS seed")
    for k in range(1, 13):
        time_iteration(rules, cal, ss, sp, regimes=(0,), no_default=True, damp=0.5, tol=1e-4,
                       max_it=1, n_gh=RE.N_GH, verbose=False)
        if k in (1, 2, 4, 8, 12):
            DD.worst_residuals(rules, cal, ss, sp, f"after {k} sweeps", k=4)


if __name__ == "__main__":
    main()
