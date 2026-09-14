"""Thin shim so the audit probe can reuse diagnostics/regimes/lottery_math's
pole scan without importing the whole regimes package (which pulls the model)."""
import numpy as np

POLE_SAFETY_FRACTION = 0.75


def pole_scan(A_cb, T, lo=0.25, hi=60.0, n=240, cond_max=1.0e4):
    """First gamma at which (I - gamma A_cb) becomes ill-conditioned.

    Same construction as diagnostics/regimes/lottery_math.closed_loop_pole and
    code/tpi.py's inline guard. Returns None when no pole is found on [lo, hi].
    """
    I = np.eye(T)
    for g in np.linspace(lo, hi, n):
        if np.linalg.cond(I - g * A_cb) > cond_max:
            return {"gamma_pole": float(g),
                    "gamma_safe_max": float(POLE_SAFETY_FRACTION * g)}
    return None
