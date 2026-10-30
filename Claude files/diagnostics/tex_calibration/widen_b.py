# THE LAST BOX STEP (2026-10-29): THE BOND-BOOK BAND WIDENED (b_band 0.12 -> RUNGS) ON THE CALIBRATED FINE APPENDIX-B
# SOLUTION, SO THE CRISIS IRF PATH (b_DD PAST THE EDGE FROM QUARTER 6, b_FD BELOW IT ON IMPACT) STAYS INSIDE THE GRID;
# EACH RUNG A WARM HARD SOLVE WITH THE CALM GUARD, THEN THE REST-POINT TARGETS RE-READ.
#   python3 -u widen_b.py <checkpoint> v1,v2,... [box key, default b_band; earlier widened keys key=value,...]
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np


def main():
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import s_process_params, IS
    from global_projection.solver_recursive.decision_rules import RuleSet
    import global_projection.solver_recursive.calibrate_stochastic as CS
    from m1_coarse import stage, calm_guard
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    cal["jacobian"] = "complex"
    sp = s_process_params(cal)
    s_hi = float(rules.grid.hi[IS])
    key = sys.argv[3] if len(sys.argv) > 3 else "b_band"
    # bands an earlier run already widened, kept (e.g. b_band=0.25)
    keep = dict((kv.split("=")[0], float(kv.split("=")[1])) for kv in sys.argv[4].split(",")) if len(sys.argv) > 4 else {}
    for b in map(float, sys.argv[2].split(",")):
        RE.BOX_KW = dict(RE.BOX_KW, s_hi=s_hi, **keep, **{key: b})
        new = RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1, refine=(IS, 9, "cross")), 2), rules)
        ok, w = stage(RE, new, cal, ss, sp, f"{key} {b:.2f}", maxit=40)
        g = calm_guard(RE, new, cal, ss, sp, f"{key} {b:.2f}", maxit=40)
        if g is not None:
            ok, w = g
        print(f"[b] {key} {b:.2f}: ok {ok}  max|F| {w:.2e}", flush=True)
        if not ok:
            return
        rules = new
        S, o = CS.rest_point_read(rules, cal, ss, sp)
        m = CS.rest_point_moments(o)
        u = rules.grid.to_unit(np.atleast_2d(S))[0]
        print("[b] rest: " + "  ".join(f"{p[2]} {m[p[2]]:.5g}" for p in CS.PAIRS) + f"  (u K_D {u[0]:+.2f}, "
              f"b_DD {u[4]:+.2f})", flush=True)
        with open(HERE / "runs" / f"final_{key}{100 * b:.0f}.pkl", "wb") as fh:
            pickle.dump(dict(sv, rules=rules, cal=cal), fh)


if __name__ == "__main__":
    main()
