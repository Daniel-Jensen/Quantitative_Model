# THE FASTER FISCAL RULE BY CONTINUATION (2026-11-06, author's choice): THE APPENDIX-B NO-TPI RULES ON THE 12% BOX
# (data/global_warm_start.pkl, solved at bohn_annual 0.06) RE-SOLVED AS bohn_annual WALKS TO 0.25. EACH RUNG IS WARM,
# HARD FIRST, ELSE THE WORKING SMOOTHINGS AND A HARD POLISH; A FAILED RUNG IS HALVED. THE STEADY STATE IS REBUILT AT EACH
# RUNG (the tax response gamma_tau = bohn_annual / 4 lives in its government block).
#   python3 -u bohn_homotopy.py   (run detached; the collocation workers re-import this file)
import copy
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

START = 0.06
RUNGS = (0.10, 0.15, 0.20, 0.25)
MAX_HALVINGS = 3
OUT = HERE / "runs" / "rules_12pct_bohn0.250.pkl"


def economy_at(base, bohn):
    # THE CALIBRATION, STEADY STATE AND RISK PROCESS AT ONE BOHN COEFFICIENT.
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    cal = dict(base, bohn_annual=bohn)
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    return cal, ss, sp


def solve_rung(RE, rules, cal, ss, sp, label):
    # ONE RUNG FROM rules: HARD, ELSE WORKING SMOOTHINGS THEN THE HARD POLISH; RETURNS (solved rules or None, max|F|).
    trial = copy.deepcopy(rules)
    ok, _, w = RE._stage(trial, cal, ss, sp, (0, 1), False, f"{label} hard", True, warm=0, maxit=12)
    if not ok:
        trial = copy.deepcopy(rules)
        ok, _, w = RE._stage(trial, dict(cal, **RE.HARD_WORK), ss, sp, (0, 1), False, f"{label} working", True,
                             warm=0, maxit=12)
        if ok:
            ok, _, w = RE._stage(trial, cal, ss, sp, (0, 1), False, f"{label} polish", True, warm=0, maxit=12)
    return (trial if ok else None), w


if __name__ == "__main__":
    from crisis_shock import setup12
    from key_test import max_resid
    from global_projection.solver_recursive import recursive_experiment as RE
    base, _, _, rules = setup12()
    cal0, ss0, sp0 = economy_at(base, START)
    print(f"[bohn] warm rules at {START}: max|F| {max_resid(rules, cal0, ss0, sp0):.2e}", flush=True)
    last, todo, halvings = START, list(RUNGS), 0
    while todo:
        b = todo[0]
        cal, ss, sp = economy_at(base, b)
        print(f"[bohn] rung {b:.4f}: start max|F| {max_resid(rules, cal, ss, sp):.2e}", flush=True)
        got, w = solve_rung(RE, rules, cal, ss, sp, f"bohn {b:.4f}")
        print(f"[bohn] rung {b:.4f}: {'ok' if got is not None else 'FAILED'}, max|F| {w:.2e}", flush=True)
        if got is None:
            if halvings >= MAX_HALVINGS:
                print(f"[bohn] stopped at {last:.4f}", flush=True)
                sys.exit(1)
            todo.insert(0, 0.5 * (last + b))
            halvings += 1
            continue
        rules, last = got, b
        todo.pop(0)
        with open(HERE / "runs" / f"rules_12pct_bohn{b:.3f}.pkl", "wb") as fh:
            pickle.dump(rules, fh)
    print(f"[bohn] solved at {last:.4f}: {OUT}", flush=True)
