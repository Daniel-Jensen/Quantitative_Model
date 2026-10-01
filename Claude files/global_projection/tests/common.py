# SHARED TEST FIXTURES: SOLVE THE STEADY STATE ONCE PER PROCESS.
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402
from config.calibration import get_calibration  # noqa: E402
from config.steady_state import solve_steady_state  # noqa: E402

_CACHE = {}


def get_ss():
    # CALIBRATION + STEADY STATE, SOLVED ONCE AND CACHED FOR THE PROCESS.
    if "ss" not in _CACHE:
        cal = get_calibration()
        _CACHE["cal"] = cal
        _CACHE["ss"] = solve_steady_state(cal, verbose=False)
    return _CACHE["cal"], _CACHE["ss"]


def ss_input_paths(cal, ss):
    # CONSTANT STEADY-STATE INPUT PATHS FOR THE BANK BLOCK (LENGTH T).
    T = cal["T"]
    return dict(
        Kap_D=np.full(T, ss["Kap_D_ss"]), Kap_F=np.full(T, ss["Kap_F_ss"]),
        Q_D=np.ones(T), Q_F=np.ones(T),
        rk_D=np.full(T, ss["rk_D_ss"]), rk_F=np.full(T, ss["rk_F_ss"]),
        rdep_D=np.full(T, cal["r_dep_D_target"]),
        rdep_F=np.full(T, cal["r_dep_F_target"]),
        p_path=np.full(T, ss["p_ss"]),
    )
