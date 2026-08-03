"""The dY decomposition must reproduce the market_clearing_D identity exactly."""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, os.path.join(ROOT, "diagnostics", "regimes"))
sys.path.insert(0, HERE)


def _synthetic():
    """An IRF/SS pair that satisfies the identity by construction.

    Y_D = P_CES_D*C_D + I_D + G_D + Phi_D + T_D + NX_D, linearised:
    dY = P_ss*dC + C_ss*dP + dI + dG + dPhi + dT + dNX
    """
    rng = np.random.default_rng(0)
    T = 20
    ss = {"P_CES_D_ss": 1.3, "C_D_ss": 0.7}
    irf = {k: rng.normal(size=T) * 1e-3
           for k in ("C_D", "P_CES_D", "I_D", "G_D", "Phi_D", "T_D", "NX_D")}
    irf["Y_D"] = (ss["P_CES_D_ss"] * irf["C_D"] + ss["C_D_ss"] * irf["P_CES_D"]
                  + irf["I_D"] + irf["G_D"] + irf["Phi_D"] + irf["T_D"] + irf["NX_D"])
    return irf, ss


def test_decomposition_closes_on_synthetic_input():
    from e2_dy_decomposition import decompose_dY
    irf, ss = _synthetic()
    components, residual = decompose_dY(irf, ss)
    assert np.max(np.abs(residual)) < 1e-15


def test_components_sum_to_dY():
    from e2_dy_decomposition import decompose_dY
    irf, ss = _synthetic()
    components, _ = decompose_dY(irf, ss)
    total = sum(components.values())
    assert np.allclose(total, irf["Y_D"], atol=1e-15)


def test_missing_term_is_detected_not_absorbed():
    """Dropping a real term must show up in the residual, never be silently absorbed."""
    from e2_dy_decomposition import decompose_dY
    irf, ss = _synthetic()
    irf["NX_D"] = irf["NX_D"] + 1e-4   # perturb one component only
    _, residual = decompose_dY(irf, ss)
    assert np.max(np.abs(residual)) > 1e-6


def test_absent_optional_term_is_treated_as_zero():
    """T_D may legitimately be absent from the cache (T0=T1=0)."""
    from e2_dy_decomposition import decompose_dY
    irf, ss = _synthetic()
    irf["Y_D"] = irf["Y_D"] - irf["T_D"]
    del irf["T_D"]
    components, residual = decompose_dY(irf, ss)
    assert np.max(np.abs(residual)) < 1e-15
    assert np.all(components["macropru_tax"] == 0.0)
