# M1, CONTINUED (2026-10-27): THE JOINT STAGE FROM ITS CHECKPOINT WITH MORE NEWTON STEPS (VERBOSE), THEN THE REST OF M1.
#   python3 -u m1_resume.py <m1 joint checkpoint> [maxit]
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import m1_coarse as M
from diag_path import FLOOR


def main():
    from calibration.global_projection import get_calibration
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive import recursive_experiment as RE
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    maxit = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    rules, cal, ss, vals = sv["rules"], sv["cal"], sv["ss"], dict(sv["vals"])
    sp = s_process_params(cal)
    base = get_calibration(); base.update(nw_floor_frac=FLOOR, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3)
    ok, it, w = RE._stage(rules, cal, ss, sp, (0, 1), False, "joint eps 1e-2 (resumed)", True, warm=0, maxit=maxit)
    print(f"[joint eps 1e-2 resumed] ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
    M.save("joint_eps1e-2_resumed", rules=rules, cal=cal, ss=ss, vals=vals)
    for k in range(3):
        if ok:
            break
        # the rest point sits at the box edges: re-centre the bank stocks there (build_state_box's at) and re-solve
        from global_projection.solver_recursive.decision_rules import RuleSet
        from global_projection.solver_recursive.state_grid import _box_centre
        S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)
        at = {key: float(S0[RE._STOCK_INDEX[key]]) for key in RE.BOX_STOCKS}
        c0 = _box_centre(ss, cal)
        print(f"[recentre {k + 1}] " + ", ".join(f"{key} {100 * (v / c0[key] - 1):+.2f}%" for key, v in at.items()), flush=True)
        new = RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1, at=at), 2), rules)
        new.box = dict(mu=1, mu_vec=None, rot=None, centre=None, box_kw=dict(RE.BOX_KW), aligned=True, at=at)
        rules = new
        ok, it, w = RE._stage(rules, cal, ss, sp, (0, 1), False, f"joint re-centred {k + 1}", True, warm=0, maxit=maxit)
        print(f"[joint re-centred {k + 1}] ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
        M.save(f"joint_recentred{k + 1}", rules=rules, cal=cal, ss=ss, vals=vals)
    if ok:
        M.finish(RE, base, rules, cal, ss, sp, vals)
    else:
        M.reads(RE, rules, cal, ss, sp, "joint eps 1e-2 (unconverged)")


if __name__ == "__main__":
    main()
