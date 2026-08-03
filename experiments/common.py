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

    Patches the MODULE ATTRIBUTE calibration.get_calibration. This reaches
    regime_model.build_caches, _calibration_fingerprint and _live_psilam because
    each does a function-local `from calibration import get_calibration`, resolved
    at call time. It then propagates onward through the dict those functions pass
    to solve_steady_state(cal), which takes the calibration as an argument rather
    than fetching it.

    FOOTGUN: a module-level `from calibration import get_calibration` binds the
    original function object at import and will NOT see the override. Any new
    experiment must either call through one of the functions above or import the
    module and call calibration.get_calibration() at use time.
    """
    import calibration
    original = calibration.get_calibration

    def patched():
        cal = original()
        unknown = set(overrides) - set(cal)
        if unknown:
            raise KeyError(
                f"calibration_override got unknown key(s): {sorted(unknown)}. "
                f"A typo here would silently leave the parameter at its default and "
                f"produce a wrong-but-plausible number — the exact failure mode the "
                f"retired audit_artifacts/ harness had. Check spelling against "
                f"code/calibration.py.")
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
    # This runs after a solve that can take up to 20 minutes, so a missing git
    # binary (FileNotFoundError/OSError) or any other git failure must degrade to
    # "unknown" rather than crash and lose the result. stderr is silenced so a
    # not-a-repo checkout doesn't print "fatal: not a git repository" to the console.
    try:
        sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                      cwd=ROOT, text=True,
                                      stderr=subprocess.DEVNULL, timeout=10).strip()
    except (subprocess.SubprocessError, OSError):
        sha = "unknown"
    # A provenance stamp taken with uncommitted edits is otherwise indistinguishable
    # from a clean run at the same SHA — a real gap for numbers headed into the
    # paper. None (not False) on failure: we must not claim "clean" when we could
    # not check.
    try:
        dirty_out = subprocess.check_output(["git", "status", "--porcelain"],
                                            cwd=ROOT, text=True,
                                            stderr=subprocess.DEVNULL, timeout=10)
        git_dirty = bool(dirty_out.strip())
    except (subprocess.SubprocessError, OSError):
        git_dirty = None
    return {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "git_sha": sha,
        "git_dirty": git_dirty,
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


def _json_default(obj):
    """numpy -> JSON. Arrays become lists; numpy scalars become Python numbers."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(f"{type(obj).__name__} is not JSON-serialisable: {obj!r}")


def write_results(name, payload):
    """Persist a machine-readable result so downstream tooling never parses markdown.

    allow_nan=False on purpose: a NaN in a results file is either a real modelling
    failure or a division by an empty denominator, and both should surface here
    rather than travel silently into a table. Encode a deliberate 'not applicable'
    as None (JSON null) instead — e.g. loading at the passive regime, where there
    are no purchases and the expected-loss denominator is zero.
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"{name}.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=_json_default, allow_nan=False)
    return path
