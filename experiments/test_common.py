"""Tests for experiments/common.py and the regime_model cache-key contract."""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, os.path.join(ROOT, "diagnostics", "regimes"))
sys.path.insert(0, HERE)


def test_cache_path_reflects_calibration_override():
    """cache_path must read the LIVE calibration, not an import-time constant.

    Regression guard: regime_model computed CAL_FINGERPRINT at import. Under an
    override applied afterwards, E3's cache would be written to the baseline
    filename and silently clobber it.
    """
    import regime_model
    from common import calibration_override

    base = regime_model.cache_path(8.5)
    with calibration_override(writeoff_enabled_D=1.0):
        overridden = regime_model.cache_path(8.5)
    after = regime_model.cache_path(8.5)

    assert base != overridden, "override did not change the cache filename"
    assert base == after, "override leaked out of the context manager"


def test_cache_path_carries_schema_version():
    import regime_model
    assert f"_v{regime_model.CACHE_SCHEMA}_" in regime_model.cache_path(8.5)


def test_required_outputs_cover_the_dy_identity():
    """E2's identity needs these in the cache or it cannot close."""
    import regime_model
    for name in ("Y_D", "C_D", "P_CES_D", "I_D", "NX_D", "Phi_D", "def_rate_D"):
        assert name in regime_model.REQUIRED, f"{name} missing from REQUIRED"


def test_calibration_override_restores_on_exception():
    from calibration import get_calibration
    from common import calibration_override

    original = get_calibration()["writeoff_enabled_D"]
    with pytest.raises(RuntimeError):
        with calibration_override(writeoff_enabled_D=1.0):
            raise RuntimeError("boom")
    assert get_calibration()["writeoff_enabled_D"] == original


def test_bp_ann_converts_quarterly_rate_to_annual_basis_points():
    from common import BP_ANN, bp_ann
    assert BP_ANN == 4.0e4
    assert bp_ann(np.array([0.00376]))[0] == pytest.approx(150.4, abs=0.1)


def test_pct_of_ss_divides_by_the_steady_state_level():
    """Level deviations are only percentages where the SS level is 1."""
    from common import pct_of_ss
    out = pct_of_ss(np.array([-0.072270]), 2.138)
    assert out[0] == pytest.approx(-3.380, abs=0.01)
