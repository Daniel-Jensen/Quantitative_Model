# THE PRODUCTION RISK SOLVE (solve_recursive, WITH THE INSTRUMENT WALK) AT OVERRIDDEN INSTRUMENTS, SAVED; --tfp ALSO
# SOLVES THE TFP RULES THERE. Writes runs/<TAG>.pkl next to this file (keep work out of /tmp: a reboot wipes it).
# usage: python3 -u walk_solve.py TAG key=value[,key=value...] [--tfp]   e.g. credit_spread_target_D=0.0035
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def main():
    # SOLVE, SAVE, AND OPTIONALLY SOLVE THE TFP RULES.
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    import global_projection.solver_recursive.recursive_experiment as RE
    import run
    tag = sys.argv[1]
    over = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in sys.argv[2].split(",")} if len(sys.argv) > 2 and "=" in sys.argv[2] else {}
    os.makedirs(os.path.join(HERE, "runs"), exist_ok=True)
    t0 = time.perf_counter()
    cal = get_calibration()
    cal["nw_floor_frac"] = run.NW_FLOOR
    cal.update(over)
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    print(f"[config] {over}, S_REFINE {RE.S_REFINE} ({RE.S_REFINE_MODE})", flush=True)
    base = []
    rules = RE.solve_recursive(cal, ss, sp, s_refine=RE.S_REFINE, verbose=True, base_out=base)
    print(f"[VERDICT risk] solve_ok={rules.solve_ok} worst={rules.solve_worst:.2e} points {rules.grid.n} "
          f"({time.perf_counter() - t0:.0f}s)", flush=True)
    with open(os.path.join(HERE, "runs", f"{tag}.pkl"), "wb") as fh:
        pickle.dump(dict(rules=rules.copy(), base=base[0].copy() if base else None, over=over,
                         mu_eps=cal["mu_eps"]), fh)
    if "--tfp" in sys.argv:
        RE.solve_tfp(cal, ss, sp)
        print(f"[VERDICT tfp] ({time.perf_counter() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
