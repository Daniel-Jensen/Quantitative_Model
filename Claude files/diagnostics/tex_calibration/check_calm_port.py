# CHECK (2026-10-28): THE PORTED recursive_experiment._calm_branch ON THE APPENDIX-B d0 CHECKPOINT (CRISIS AT NODES
# 9 AND 22) MUST REPRODUCE branch_select.py's CONVERGED d0 (2.5e-12, NO CRISIS NODE LEFT).
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))


def main():
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import s_process_params
    sv = pickle.load(open(HERE / "runs" / "d0_appendixB_bohn006.pkl", "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    cal.update(nw_floor_frac=0.5, jacobian="complex")
    sp = s_process_params(cal)
    print("crisis nodes before:", RE._crisis_nodes(rules, cal, ss, sp, True), flush=True)
    res = RE._calm_branch(rules, cal, ss, sp, (0,), True, "d0", True)
    print("result:", res, "| crisis nodes after:", RE._crisis_nodes(rules, cal, ss, sp, True), flush=True)
    assert res is not None and res[0] and res[1] <= 1e-9, res
    print("check_calm_port: PASSED")


if __name__ == "__main__":
    main()
