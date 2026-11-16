# A COARSE PRODUCTION SOLVE WITH CALIBRATION OR RISK-PROCESS OVERRIDES, SAVED FOR diag_slack.py (2026-10-22).
# usage: python3 -u probe_solve.py TAG key=value[,...] [pd_rest=0.0002] [--no-haircut-ladder] [--refine]
#   key=value pairs override the calibration; pd_rest sets the quarterly default probability at rest (s_star)
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def main():
    # SOLVE AND SAVE.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    import global_projection.solver_recursive.recursive_experiment as RE
    import run
    tag, cal_over, sp_over = sys.argv[1], {}, {}
    for a in sys.argv[2:]:
        if a.startswith("--") or "=" not in a:
            continue
        for kv in a.split(","):
            k, v = kv.split("=")
            if k == "pd_rest":
                sp_over["s_star"] = float(np.log(float(v) / (1.0 - float(v))))
            else:
                cal_over[k] = float(v)
    if "--no-haircut-ladder" in sys.argv:
        RE.HAIRCUT_LADDER = ()
    t0 = time.perf_counter()
    cal = get_calibration()
    cal["nw_floor_frac"] = run.NW_FLOOR
    cal.update(cal_over)
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    sp.update(sp_over)
    calibrate_household_anchors(cal, ss, sp)
    print(f"[config] cal {cal_over} sp {sp_over} floor {cal['nw_floor_frac']}", flush=True)
    ref = 5 if "--refine" in sys.argv else 0
    rules = RE.solve_recursive(cal, ss, sp, s_refine=ref, verbose=True)
    print(f"[VERDICT] {tag} ok={rules.solve_ok} worst={rules.solve_worst:.2e} ({time.perf_counter() - t0:.0f}s)", flush=True)
    with open(os.path.join(HERE, "runs", f"{tag}.pkl"), "wb") as fh:
        pickle.dump(dict(rules=rules.copy(), nw_floor_frac=cal["nw_floor_frac"], mu_eps=cal["mu_eps"],
                         cal_over=cal_over, sp_over=sp_over, ok=bool(rules.solve_ok),
                         worst=float(rules.solve_worst)), fh)


if __name__ == "__main__":
    main()
