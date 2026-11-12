# ONE COLLOCATION RESIDUAL ON THE 2% BOX (THE WARM RULES SEEDED AS run.py DOES), TIMED: IS THE MODEL SLOWER THERE?
import pickle
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from calibration.global_projection import get_calibration
from global_projection.steady_state import solve_steady_state
from global_projection.solver_recursive.state_grid import s_process_params, IS
from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
from global_projection.solver_recursive.decision_rules import RuleSet
from global_projection.solver_recursive import recursive_experiment as RE
from global_projection.solver_recursive.collocation import make_residual, pack

if __name__ == "__main__":
    cal = get_calibration()
    cal.update(nw_floor_frac=0.5, jacobian="complex", pd_box_top=0.02)
    RE.BOX_KW.update(m_band=0.65, o_band=0.08 / 0.65)
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    RE._soften(cal)
    warm = pickle.load(open(ROOT / "data" / "global_warm_start.pkl", "rb"))
    rules = RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1, refine=(IS, 9, "cross")), 2), warm)
    F = make_residual(rules, cal, ss, sproc, regimes=(0, 1), n_gh=rules.n_gh or 5)
    x = pack(rules, (0, 1))
    for _ in range(2):
        t = time.perf_counter(); f = F(x); dt = time.perf_counter() - t
        print(f"one F on the 2% box: {dt:.2f} s, max|F| at the seed {abs(f).max():.2e}")
