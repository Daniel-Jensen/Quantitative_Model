# M1 VARIANT (2026-10-28): THE JOINT STAGE ON AN SS-CENTRED BOX WITH A WIDER CAPITAL BAND (THE REST POINT SITS AT
# K_D +2.07% AGAINST THE +-2% BAND), SEEDED FROM THE SS-BOX JOINT CHECKPOINT; THEN THE REST OF M1.
#   M1_TAG=_k3 python3 -u m1_kband.py <joint checkpoint> <k_band>
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import m1_coarse as M
from diag_path import FLOOR


def main():
    from calibration.global_projection import get_calibration
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive import recursive_experiment as RE
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    kb = float(sys.argv[2])
    bb = float(sys.argv[3]) if len(sys.argv) > 3 else RE.BOX_KW["b_band"]
    RE.BOX_KW = dict(RE.BOX_KW, k_band=kb, b_band=bb)
    print(f"box: {RE.BOX_KW}", flush=True)
    rules, cal, ss, vals = sv["rules"], sv["cal"], sv["ss"], dict(sv["vals"])
    sp = s_process_params(cal)
    base = get_calibration(); base.update(nw_floor_frac=FLOOR, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3)
    new = RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1), 2), rules)
    new.box = dict(mu=1, mu_vec=None, rot=None, centre=None, box_kw=dict(RE.BOX_KW), aligned=True, at=None)
    ok, w = M.stage(RE, new, cal, ss, sp, f"joint k_band {kb} b_band {bb}", maxit=60)
    M.save("joint_eps1e-2", rules=new, cal=cal, ss=ss, vals=vals)
    if ok:
        M.finish(RE, base, new, cal, ss, sp, vals)
    else:
        M.reads(RE, new, cal, ss, sp, f"joint k_band {kb} b_band {bb} (unconverged)")


if __name__ == "__main__":
    main()
