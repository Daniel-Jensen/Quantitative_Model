# T1.6 PROBE (2026-10-26): THE SAME d0 ECONOMY (THE PATH'S HOME-BIAS ENDPOINT, FLOOR 6e-8 ON THE V5 BOX) ON A BOX WITH
# THE BANK-STOCK BANDS SCALED BY k, SEEDED BY EVALUATING THE SOLVED RULES AT THE NEW NODES. DOES THE NEAR-FOLD GO?
#   python3 -u probe_box.py <checkpoint> <k>
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np


def main():
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive import recursive_experiment as RE
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    k = float(sys.argv[2])
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    sp = s_process_params(cal)
    RE.BOX_KW = {key: (v * k if key in ("k_band", "p_band_D", "p_band_F", "b_band", "w_band") else v)
                 for key, v in RE.BOX_KW.items()}
    print(f"box bands x{k}: {RE.BOX_KW}", flush=True)
    new = RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1), 2), rules)
    ok, it, w = RE._stage(new, cal, ss, sp, (0,), True, f"d0 box x{k}", True, warm=0, maxit=30)
    print(f"[box x{k}] ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
    with open(HERE / "runs" / f"probe_box_x{k}.pkl", "wb") as fh:
        pickle.dump(dict(rules=new, cal=cal, ss=ss, vals=sv["vals"], box_kw=RE.BOX_KW), fh)


if __name__ == "__main__":
    main()
