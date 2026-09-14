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


def load_cache(psilam=None):
    """Load the response-matrix cache for the live calibration."""
    from regime_model import load_cache as _load, _live_psilam
    return _load(_live_psilam() if psilam is None else psilam)


def cache_outputs(cache):
    """Names of every output cached with a cb_buy_D column."""
    return sorted({k.split("__")[0] for k in cache if k.endswith("__cb_buy_D")})


def irf_from_cache(cache, cb_path, eps):
    """Closed-loop IRF: response = M_def @ eps + M_cb @ cb_path.

    Mirrors run_regimes.irf_all so the two agree by construction.
    """
    out = {}
    for name in cache_outputs(cache):
        out[name] = (cache[f"{name}__shock_def_D"] @ eps
                     + cache[f"{name}__cb_buy_D"] @ cb_path)
    out["cb_buy_D"] = np.asarray(cb_path)
    return out


def named_regime_gammas(cache):
    """passive / medium / aggressive, with gamma SOLVED for peak-spread compression.

    Solved rather than chosen so the regimes keep their meaning across
    recalibrations, instead of a round number drifting into a different policy
    stance. Targets are spec section 7: 25% (medium), 50% (aggressive).

    ⚠ **"aggressive" no longer means 50% since the 2026-08-18 GK structural refactor.**
    The closed loop has a POLE at gamma ~ 27.3, and the maximum compression reachable
    below it is ~46.6%. The 50% target is only met on the far side of the singularity,
    which is a different branch, not a stronger version of the same policy. Rather than
    chase it there (wrong) or drop the regime (would break eight paper figures), the
    aggressive regime falls back to `lottery_math.POLE_SAFETY_FRACTION * pole` = 0.75 x
    pole, and this function prints the compression it actually achieves (~40.3%). 0.75
    rather than "as close to the pole as possible": the loading schedule is monotone in
    gamma only up to ~0.85 x pole, and at 0.98 x pole the discounted consumption gains
    reach +12% of steady-state consumption — that is the singularity talking, not the
    policy. Every downstream table already reports each regime's gamma and
    peak spread, so the artefacts stay self-describing; but any prose calling the
    aggressive regime "50% compression" is now WRONG and must say ~46.6% or, better,
    "maximum feasible". See `lottery_math.gamma_for_compression` and
    `lottery_math.closed_loop_pole`.
    """
    from lottery_math import (gamma_for_compression, closed_loop_pole, closed_loop,
                              peak, CompressionInfeasible, POLE_SAFETY_FRACTION)
    A_def = cache["spread_rb__shock_def_D"]
    A_cb = cache["spread_rb__cb_buy_D"]
    eps = np.asarray(cache["dShock_def_D"])
    assert float(A_cb[0, 0]) < 0.0, (
        f"A_cb[0,0]={float(A_cb[0,0]):+.4e} >= 0: CB purchases WIDEN the spread, so "
        "compression targeting is infeasible. This is the ms-regime SA-1 pathology, "
        "which must be absent on main — investigate before reporting anything.")
    p0 = peak(closed_loop(A_def, A_cb, eps, 0.0)[0])
    out = {"passive": 0.0}
    for name, target in (("medium", 0.25), ("aggressive", 0.50)):
        try:
            out[name] = float(gamma_for_compression(A_def, A_cb, eps, target=target))
        except CompressionInfeasible as exc:
            pole = closed_loop_pole(A_cb)
            assert pole is not None, exc          # infeasible with no pole = a real bug
            out[name] = POLE_SAFETY_FRACTION * float(pole)
            got = 1.0 - peak(closed_loop(A_def, A_cb, eps, out[name])[0]) / p0
            print(f"  [named_regime_gammas] {name}: {100*target:.0f}% compression is "
                  f"UNREACHABLE below the closed-loop pole at gamma = {pole:.2f} "
                  f"(max reachable there ~46.6%). Falling back to "
                  f"POLE_SAFETY_FRACTION x pole = {out[name]:.3f}, which achieves "
                  f"{100*got:.2f}% and keeps a margin from the singularity. Do not "
                  f"describe this regime as {100*target:.0f}% compression.")
    return out


def regime_irfs(cache):
    """{regime_name: (gamma, irf_dict)} for the three named regimes."""
    from lottery_math import closed_loop
    A_def = cache["spread_rb__shock_def_D"]
    A_cb = cache["spread_rb__cb_buy_D"]
    eps = np.asarray(cache["dShock_def_D"])
    out = {}
    for name, gamma in named_regime_gammas(cache).items():
        _spread, cb = closed_loop(A_def, A_cb, eps, gamma)
        out[name] = (gamma, irf_from_cache(cache, cb, eps))
    return out
