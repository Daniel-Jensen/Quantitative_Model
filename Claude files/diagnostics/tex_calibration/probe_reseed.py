# PROBE (2026-10-26): THE RELATIVE RE-SEED ON THE SAVED PATH SOLUTION: THE SAME ECONOMY (MUST GIVE ~1e-10),
# THEN TINY HOME-BIAS MOVES; THE SEED'S RESIDUAL AND WHERE IT SITS, BEFORE ANY NEWTON STEP.
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np
import diag_d0 as DD
from diag_homotopy import OLD, economy
from diag_path import relative_seed


def main():
    from calibration.global_projection import get_calibration
    from global_projection.solver_recursive.state_grid import STATE_NAMES
    from global_projection.solver_recursive import recursive_experiment as RE
    DD.STATE = STATE_NAMES
    sv = pickle.load(open(HERE / "runs" / "diag_path_foreign_share.pkl", "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    base = get_calibration(); base.update(nw_floor_frac=0.15, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3)
    new_hb = float(base["omega_home_D"])
    for frac in (0.0, 1e-3, 1e-2):
        vals = dict(sv["vals"], omega_home_D=OLD["omega_home_D"] + frac * (new_hb - OLD["omega_home_D"]))
        cal_k, ss_k, sp_k = economy(base, vals)
        seed = relative_seed(rules, cal, ss, cal_k, ss_k)
        print(f"\nhome-bias fraction {frac}: p_ss {ss_k['p_ss']:.6f} (was {ss['p_ss']:.6f})  "
              f"hh_T_D {ss_k.get('hh_T_D', 0):+.2e}", flush=True)
        DD.worst_residuals(seed, cal_k, ss_k, sp_k, f"seed at fraction {frac}", k=6)


if __name__ == "__main__":
    main()
