# RE-RUN run.py's TPI COMPARISON (STEP 6's compare, THEN STEP 7's FIGURES) ON THE SAVED RULES, WITHOUT RE-SOLVING: AFTER A
# CHANGE TO A READER (2026-11-06: the rest-point row now reads exactly). Writes results/GLOBAL/data/tpi and the figures.
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

if __name__ == "__main__":
    import run
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive import tpi_experiment
    from global_projection.solver_recursive.output_decomposition import BOND_CHANNELS
    from global_projection.reporting.prints import print_tpi_report
    from results_io import save, load
    data = ROOT / "results" / "GLOBAL" / "data"
    cal = get_calibration()
    cal.update(nw_floor_frac=run.NW_FLOOR, jacobian=run.JACOBIAN, pd_box_top=run.PD_BOX_TOP)
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    cal["tpi_cap_bp"] = float(run.TPI_CAP)
    rules_off = pickle.load(open(data / "rules_risk.pkl", "rb"))
    rules_on = pickle.load(open(data / "rules_tpi.pkl", "rb"))
    old = load(data / "tpi")
    tpi = tpi_experiment.compare(cal, ss, sproc, rules_on, rules_off, pd_shock=run.RISK_SHOCK_PD)
    # the solve's own verdicts are not the comparison's to change
    tpi.update(cap_solved=old["cap_solved"], solve_ok=old["solve_ok"])
    print_tpi_report(tpi, BOND_CHANNELS)
    save(data / "tpi", tpi)
    run.plot_global()
