# IS THE TPI OPERATIONAL AT APPENDIX B? (2026-10-30) ON THE FINE GRID THE CAP WALK STARTS BINDING AT THE p^d-12% BOX TOP
# EVEN AT 700 bp (Q_bD 0.64 AGAINST A FLOOR OF 0.73) AND CREEPS AT ~10 MIN PER ITERATION. HERE: A COARSE BASELINE SEEDED
# FROM THE SOLVED FINE RULES, THEN THE COARSE TPI CAP WALK FROM ABOVE THE BOX-TOP SPREAD DOWN TO THE TARGET.
#   python3 -u tpi_check.py [cap start bp] [cap step bp]
import os, pickle, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)


def main():
    import run
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive import recursive_experiment as RE
    cal = get_calibration(); cal.update(nw_floor_frac=run.NW_FLOOR, jacobian=run.JACOBIAN)
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    fine = pickle.load(open(run.WARM_START, "rb"))
    RE.TPI_CAP_START = float(sys.argv[1]) if len(sys.argv) > 1 else 1500.0
    RE.TPI_CAP_STEP = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
    # the coarse no-TPI baseline, seeded from the fine rules and solved with the exact conditions
    t = time.perf_counter()
    coarse = RE.solve_recursive(dict(cal), ss, sp, base=RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1), 2), fine),
                                s_refine=0)
    print(f"[base] coarse no-TPI: ok {coarse.solve_ok}  worst {coarse.solve_worst:.1e}  ({time.perf_counter() - t:.0f}s)",
          flush=True)
    # the TPI on the coarse grid only (no refinement): does the cap walk reach the target?
    t = time.perf_counter()
    c = dict(cal, tpi_on=True, tpi_cap_bp=float(run.TPI_CAP))
    on = RE.solve_recursive(c, ss, sp, base=coarse, s_refine=0)
    print(f"[tpi] coarse TPI: ok {on.solve_ok}  worst {on.solve_worst:.1e}  cap solved {on.tpi_cap_solved}  "
          f"(target {run.TPI_CAP}; {time.perf_counter() - t:.0f}s)", flush=True)
    with open(HERE / "runs" / "tpi_check.pkl", "wb") as fh:
        pickle.dump(dict(coarse=coarse, on=on), fh)


if __name__ == "__main__":
    main()
