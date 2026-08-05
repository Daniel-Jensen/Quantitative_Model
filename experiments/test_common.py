"""Tests for experiments/common.py and the regime_model cache-key contract.

The first three tests below (cache_path override/leak, schema version, REQUIRED
coverage) exercise diagnostics/regimes/regime_model.py rather than common.py — they
live here, not co-located with regime_model.py, because they are regression guards
for the common.calibration_override contract that this file is otherwise
responsible for.
"""
import json
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
    """The monkeypatch must be undone even when the block raises.

    Assert on the MODULE ATTRIBUTE, not on a name imported beforehand: a local
    `from calibration import get_calibration` binds the original function object
    and would pass even if the patch leaked. That is what the first version of
    this test did, and it guarded nothing.
    """
    import calibration
    from common import calibration_override

    before = calibration.get_calibration
    with pytest.raises(RuntimeError):
        with calibration_override(writeoff_enabled_D=1.0):
            raise RuntimeError("boom")
    assert calibration.get_calibration is before, "override leaked after an exception"
    assert calibration.get_calibration()["writeoff_enabled_D"] == 0.0


def test_calibration_override_rejects_unknown_key():
    """A typo'd override key (e.g. psi_lambda_b_D for psi_lambda_B_D) must be loud.

    Silently adding a junk key and leaving the intended parameter at its default
    is exactly the "silently ran a different model" failure mode this project has
    already been burned by (the retired audit_artifacts/ harness).
    """
    import calibration
    from common import calibration_override

    with pytest.raises(KeyError):
        with calibration_override(psi_lambda_b_D=999.0):
            calibration.get_calibration()


def test_bp_ann_converts_quarterly_rate_to_annual_basis_points():
    from common import BP_ANN, bp_ann
    assert BP_ANN == 4.0e4
    assert bp_ann(np.array([0.00376]))[0] == pytest.approx(150.4, abs=0.1)


def test_pct_of_ss_divides_by_the_steady_state_level():
    """Level deviations are only percentages where the SS level is 1."""
    from common import pct_of_ss
    out = pct_of_ss(np.array([-0.072270]), 2.138)
    assert out[0] == pytest.approx(-3.380, abs=0.01)


def test_write_results_round_trips_a_numpy_array_as_a_list(tmp_path, monkeypatch):
    """IRF paths are numpy arrays; json.dump(default=float) chokes on those (it
    only handles length-1 arrays), so write_results needs its own array handling."""
    import common
    monkeypatch.setattr(common, "RESULTS_DIR", str(tmp_path))

    path = common.write_results("array_case", {"irf": np.array([1.0, 2.0, 3.0])})

    with open(path) as fh:
        loaded = json.load(fh)
    assert loaded["irf"] == [1.0, 2.0, 3.0]


def test_write_results_rejects_nan(tmp_path, monkeypatch):
    """A NaN payload is a modelling failure (or a 0/0), not a value to persist
    silently as the non-JSON-standard token `NaN`. Must fail loudly here rather
    than travel into a results table."""
    import common
    monkeypatch.setattr(common, "RESULTS_DIR", str(tmp_path))

    with pytest.raises(ValueError):
        common.write_results("nan_case", {"bad": float("nan")})
