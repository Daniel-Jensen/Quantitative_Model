# CHECK (2026-10-30): DOES run.py's TFP SOLVE (solve_tfp: THE NO-DEFAULT COARSE RULES FROM THE SS) CONVERGE AT APPENDIX B?
import os, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)


def main():
    import run
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.recursive_experiment import solve_tfp, tfp_irf
    cal = get_calibration(); cal.update(nw_floor_frac=run.NW_FLOOR, jacobian=run.JACOBIAN)
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    t = time.perf_counter()
    import pickle
    warm = pickle.load(open(run.WARM_START, "rb")) if run.WARM_START and run.WARM_START.exists() else None
    rules = solve_tfp(cal, ss, sp, mu=run.MU, seed=warm)
    print(f"[tfp] solve_ok {getattr(rules, 'solve_ok', '?')}  ({time.perf_counter() - t:.0f}s)", flush=True)
    irf = tfp_irf(rules, cal, ss, sp, dz=run.TFP_SHOCK)
    print("[tfp] impact Y %+.3f%%  C %+.3f%%  I %+.3f%%" % (irf["Y"][0], irf["C"][0], irf["I"][0]), flush=True)


if __name__ == "__main__":
    main()
