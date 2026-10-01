# SHARED TEST FIXTURES: SOLVE THE STEADY STATE ONCE PER PROCESS.
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from calibration.global_projection import get_calibration  # noqa: E402
from global_projection.steady_state import solve_steady_state  # noqa: E402

_CACHE = {}


def get_ss():
    # CALIBRATION + STEADY STATE, SOLVED ONCE AND CACHED FOR THE PROCESS.
    if "ss" not in _CACHE:
        cal = get_calibration()
        _CACHE["cal"] = cal
        _CACHE["ss"] = solve_steady_state(cal, verbose=False)
    return _CACHE["cal"], _CACHE["ss"]


