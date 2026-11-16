# T1 PROBE (2026-10-26): THE d0 LADDER AT pi = 0 WITH FEWER WARM TIME-ITERATION SWEEPS BEFORE THE FIRST NEWTON.
#   python3 -u diag_warm.py <warm sweeps> <floor>
import os, sys, pickle
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np
import diag_d0 as DD


def main():
    warm, floor = int(sys.argv[1]), float(sys.argv[2])
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive import recursive_experiment as RE
    DD.STATE = STATE_NAMES
    d = list(RE._stage.__defaults__); d[1] = warm; RE._stage.__defaults__ = tuple(d)  # backend, warm, maxit, sets
    cal = get_calibration(); cal["nw_floor_frac"] = floor
    if len(sys.argv) > 3:  # solve the SS consistent with this smoothing (bank.ic_slack_factor)
        cal["mu_eps"] = float(sys.argv[3])
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    rules = RuleSet.from_ss(RE.normal_box(ss, cal, mu=1), ss, cal, n_regimes=2); rules.n_gh = RE.N_GH
    RE._soften(cal); cal["mu_eps"] = RE.MU_EPS_START
    print(f"SS lambda {cal['lambda_K_D']:.5f}  mu_ss {ss['ss_bank_D']['mu_ss']:.5f}", flush=True)
    ok, it, w = RE._stage_wc(rules, cal, ss, sp, (0,), True, f"d0 warm={warm}", True)
    print(f"\n_stage_wc warm={warm} floor={floor}: ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
    DD.worst_residuals(rules, cal, ss, sp, f"warm {warm}")
    zs = DD.kink_table(rules, cal, ss, sp, f"warm {warm}")
    with open(HERE / "runs" / f"diag_warm{warm}_floor{floor}.pkl", "wb") as fh:
        pickle.dump(dict(rules=rules, cal=cal, ss=ss), fh)
    if not ok:
        return
    base = rules.copy()
    for eps in (7e-3, 3e-3, 1e-3, 0.0):
        r = base.copy(); cal["mu_eps"] = eps; calibrate_household_anchors(cal, ss, sp)
        ok, it, w = RE._stage(r, cal, ss, sp, (0,), True, f"d0 eps={eps:.0e}", True, warm=0, maxit=20)
        print(f"\n== from 1e-2 straight to mu_eps {eps:.0e}: ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
        DD.worst_residuals(r, cal, ss, sp, f"eps {eps:.0e}", k=5)
        DD.kink_table(r, cal, ss, sp, f"eps {eps:.0e}")


if __name__ == "__main__":
    main()
