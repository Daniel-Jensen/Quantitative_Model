# THE REMITTANCE TEST BY CONTINUATION (2026-11-06): THE NO-TPI ECONOMY ON THE APPENDIX-B 12% BOX WALKED FROM
# tpi_key_D 0.071 TO 1 IN RUNGS (THE DIRECT JUMP STALLED AT max|F| 0.14: F's DEFAULT-BRANCH EULER AT THE BOOK EDGES),
# EACH RUNG WARM, HARD FIRST, THEN THE WORKING SMOOTHINGS AND A HARD POLISH; A FAILED RUNG IS HALVED. THEN key_test's
# EXPERIMENTS ON BOTH RULE SETS.   python3 -u key_homotopy.py   (run detached; workers re-import this file)
import copy
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

RUNGS = (0.25, 0.45, 0.7, 1.0)
MAX_HALVINGS = 2


def solve_rung(RE, rules, cal, ss, sproc, key):
    # ONE RUNG FROM rules (MODIFIED IN PLACE ON SUCCESS): HARD, ELSE WORKING SMOOTHINGS THEN THE HARD POLISH.
    c = dict(cal, tpi_key_D=key)
    trial = copy.deepcopy(rules)
    ok, its, w = RE._stage(trial, c, ss, sproc, (0, 1), False, f"key {key:.3f} hard", True, warm=0, maxit=12)
    if not ok:
        trial = copy.deepcopy(rules)
        cw = dict(c, **RE.HARD_WORK)
        ok, its, w = RE._stage(trial, cw, ss, sproc, (0, 1), False, f"key {key:.3f} working", True, warm=0, maxit=12)
        if ok:
            ok, its, w = RE._stage(trial, c, ss, sproc, (0, 1), False, f"key {key:.3f} polish", True, warm=0, maxit=12)
    print(f"[key] rung {key:.3f}: ok {ok}, max|F| {w:.2e}", flush=True)
    return (trial if ok else None), w


if __name__ == "__main__":
    from crisis_shock import setup12
    from global_projection.solver_recursive import recursive_experiment as RE
    import key_test
    cal, ss, sproc, rules = setup12()
    last, todo, halvings = 0.071, list(RUNGS), 0
    while todo:
        key = todo[0]
        got, w = solve_rung(RE, rules, cal, ss, sproc, key)
        if got is None:
            if halvings >= MAX_HALVINGS:
                print(f"[key] stopped: rung {key:.3f} failed after {halvings} halvings; last solved key {last:.3f}",
                      flush=True)
                sys.exit(1)
            todo.insert(0, 0.5 * (last + key))
            halvings += 1
            continue
        rules, last = got, key
        todo.pop(0)
        with open(HERE / "runs" / f"rules_12pct_key{key:.3f}.pkl", "wb") as fh:
            pickle.dump(rules, fh)
    with open(key_test.OUT, "wb") as fh:
        pickle.dump(rules, fh)
    print("[key] key_D 1 solved; experiments follow", flush=True)
    sys.argv = [sys.argv[0]]
    exec(compile(open(HERE / "key_test.py").read(), str(HERE / "key_test.py"), "exec"), {"__name__": "__main__"})
