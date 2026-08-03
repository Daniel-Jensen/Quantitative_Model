"""Shared plumbing for the experiments package.

Everything here reads the PRODUCTION calibration and the PRODUCTION equation
files via diagnostics/regimes/regime_model.py. Nothing in this package carries
its own copy of a parameter value — that is what made the retired
audit_artifacts/ harness test a different model than code/main.py for weeks.
"""
import contextlib
import datetime
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (os.path.join(ROOT, "code"), os.path.join(ROOT, "diagnostics", "regimes")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

RESULTS_DIR = os.path.join(HERE, "results")
FIGURES_DIR = os.path.join(HERE, "figures")

# Quarterly rate deviation -> annualised basis points.
BP_ANN = 4.0e4


@contextlib.contextmanager
def calibration_override(**overrides):
    """Run a block with get_calibration() returning a modified dict.

    Every consumer (steady_state, regime_model, ...) calls get_calibration()
    rather than taking parameters, so patching the module attribute is the
    minimal-change way to run a variant. regime_model imports the symbol INSIDE
    its functions, so it picks the patch up at call time.
    """
    import calibration
    original = calibration.get_calibration

    def patched():
        cal = original()
        cal.update(overrides)
        return cal

    calibration.get_calibration = patched
    try:
        yield
    finally:
        calibration.get_calibration = original


def bp_ann(path):
    """Quarterly rate deviation -> annualised basis points."""
    return np.asarray(path) * BP_ANN


def pct_of_ss(path, ss_level):
    """Level deviation -> percent of steady state.

    SSJ IRFs are LEVEL deviations, so x100 is a percentage only where the SS level
    is ~1. Y_D_ss ~ 1 passes; n_inter_D_ss=2.138 and K_D_ss=10.8 do NOT. Two Stage A
    figure panels and a main.py print were mislabelled "%" on that basis (2.1x on
    net worth, 10x on capital) before the 2026-07-31 units fix.
    """
    return np.asarray(path) * 100.0 / float(ss_level)


def provenance():
    """Stamp every result with the model it came from. Read live, never hardcoded."""
    from calibration import BANK_SCOPE, EBA_CALIBRATION, get_calibration
    from regime_model import CACHE_SCHEMA, _calibration_fingerprint

    cal = get_calibration()
    try:
        sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                      cwd=ROOT, text=True).strip()
    except subprocess.CalledProcessError:
        sha = "unknown"
    return {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "git_sha": sha,
        "cal_fingerprint": _calibration_fingerprint(),
        "cache_schema": CACHE_SCHEMA,
        "BANK_SCOPE": BANK_SCOPE,
        "EBA_CALIBRATION": bool(EBA_CALIBRATION),
        "psi_lambda_B_D": float(cal["psi_lambda_B_D"]),
        "mv_rule_D": float(cal["mv_rule_D"]),
        "recovery_rate_D": float(cal["recovery_rate_D"]),
        "writeoff_enabled_D": float(cal["writeoff_enabled_D"]),
        "zeta_writeoff_D": float(cal["zeta_writeoff_D"]),
        "phi_lamb_D": float(cal["phi_lamb_D"]),
        "delta_b_D": float(cal["delta_b_D"]),
    }


def write_results(name, payload):
    """Persist a machine-readable result so downstream tooling never parses markdown."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"{name}.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=float)
    return path
