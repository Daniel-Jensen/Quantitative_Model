# Preliminary Policy Experiments Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an `experiments/` package that produces the paper's standard results set — the backstop schedule (E1), the ΔY decomposition (E2), and the S-1 writeoff variants (E3) — as one reproducible, provenance-stamped artifact.

**Architecture:** A thin runner per experiment on top of the existing solve/cache layer in `diagnostics/regimes/regime_model.py`, which builds per-output response matrices directly from the production equation files. No second copy of the model exists anywhere. `code/main.py` is untouched so it stays usable as the regression test. E1 and E2 are pure post-Jacobian numpy on one cache; E3 needs two extra full solves at overridden calibrations.

**Tech Stack:** Python 3.12 (`/opt/anaconda3/envs/ssj/bin/python` — the base env has a broken `liblapack` symlink and fails silently), `sequence-jacobian`, numpy, matplotlib, pytest 9.0.3.

**Spec:** `docs/superpowers/specs/2026-08-01-policy-experiments-design.md`

---

## Environment

Every command in this plan uses the absolute interpreter path. Never `python`, never `python3`:

```bash
/opt/anaconda3/envs/ssj/bin/python
```

Run everything from the repo root: `/Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL`

Branch: `experiments` (already created, spec committed at `49162c4`).

## File Structure

| File | Responsibility |
|---|---|
| `diagnostics/regimes/regime_model.py` *(modify)* | Cache schema version, call-time fingerprint, fingerprint assertion on load, three new cached outputs |
| `experiments/__init__.py` *(create)* | Marks the package; empty |
| `experiments/.gitignore` *(create)* | Excludes `*.npz`, mirrors `diagnostics/regimes/.gitignore` |
| `experiments/common.py` *(create)* | Calibration override, cache access, regime γ solve, unit helpers, provenance, results writer |
| `experiments/e1_backstop_schedule.py` *(create)* | Named-regime table, A5-1's three German objects, loading schedule, secondary welfare |
| `experiments/e2_dy_decomposition.py` *(create)* | ΔY decomposition against the `market_clearing_D` identity |
| `experiments/e3_writeoff_s1.py` *(create)* | E3a/E3b variant solves and comparison |
| `experiments/run_all.py` *(create)* | Orchestrator; writes `docs/experiments_results.md` |
| `experiments/test_common.py` *(create)* | Tests for override, fingerprint, unit helpers |
| `experiments/test_e2_identity.py` *(create)* | Tests for decomposition closure on synthetic input |
| `docs/experiments_results.md` *(generated)* | The deliverable |

---

## Task 1: Cache schema + call-time fingerprint

The cache filename currently keys on an **import-time** constant. E3 overrides the calibration *after* import, so its cache would be written to the baseline filename and overwrite it silently. This task fixes that and adds the three outputs E1/E2 need.

**Files:**
- Modify: `diagnostics/regimes/regime_model.py:68-78` (output lists), `:87-110` (fingerprint), `:243-248` (load)
- Test: `experiments/test_common.py`

- [ ] **Step 1: Create the package skeleton**

```bash
mkdir -p experiments/results experiments/figures
touch experiments/__init__.py
printf '*.npz\n_*.log\n__pycache__/\n' > experiments/.gitignore
```

- [ ] **Step 2: Write the failing test**

Create `experiments/test_common.py`:

```python
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
```

- [ ] **Step 3: Run the test to verify it fails**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest experiments/test_common.py -v
```

Expected: collection errors / `ModuleNotFoundError: No module named 'common'` — `common.py` does not exist yet and `regime_model` has no `CACHE_SCHEMA`.

- [ ] **Step 4: Add the schema constant and the three new outputs**

In `diagnostics/regimes/regime_model.py`, replace the `REQUIRED` / `OPTIONAL` / `SS_META` block at lines 68-78 with:

```python
# Bump whenever REQUIRED / OPTIONAL / SS_META change. It goes in the cache
# FILENAME: the calibration fingerprint alone cannot detect a schema change, so
# without this an old cache would reload under the same name missing the new
# keys — silently, because irf_all discovers outputs by scanning cache keys.
CACHE_SCHEMA = 2

REQUIRED = ["spread_rb", "rb_D", "rb_F", "q_b_D", "q_b_F", "Y_D", "C_D", "I_D",
            "NX_D", "K_D", "n_inter_D", "b_D_D", "b_D_F", "b_gov_D", "U_D", "U_F",
            "TAX_D", "P_CES_D",
            # Added for the experiments package (schema 2):
            #   Phi_D, def_rate_D — Phi_D closes the market_clearing_D identity
            #   for E2; def_rate_D is the off-path expected-loss leg for E1's
            #   A5-1 reporting (cb_pnl reads it).
            "Phi_D", "def_rate_D"]
# OPTIONAL: logged loudly if missing, never silently dropped.
# (cb_flow_D excluded — it's the CB inter-block conduit flow, unused downstream, and
#  SSJ returns its cb_buy_D Jacobian as a non-array object; keep the cache clean.)
#
# T_D is OPTIONAL, not REQUIRED, on purpose: T0=T1=0 so the macroprudential bond
# tax is identically zero and SSJ may omit it from G.outputs entirely. Zero-filling
# is the CORRECT value here rather than a silent hole — and E2's closure assertion
# catches it either way if that ever stops being true.
OPTIONAL = ["G_D", "ra_D", "lambda_gk_D", "theta_D", "GINI_WEALTH", "GINI_C",
            "div_fund_D", "T_D"]
SS_META  = ["q_b_D_ss:q_b_D", "b_D_D_ss:b_D_D", "b_gov_D_ss:b_gov_D", "Y_D_ss:Y_D",
            "C_D_ss:C_D", "I_D_ss:I_D", "NX_D_ss:NX_D", "n_inter_D_ss:n_inter_D",
            "K_D_ss:K_D", "TAX_D_ss:TAX_D", "P_CES_D_ss:P_CES_D",
            "beta_D:beta_D", "beta_F:beta_F", "EL_price_D:EL_price_D",
            # schema 2: needed by E1's cb_pnl port and E2's identity
            "delta_b_D_ss:delta_b_D", "q_b_F_ss:q_b_F",
            "Phi_D_ss:Phi_D"]
```

- [ ] **Step 5: Make the fingerprint call-time**

In the same file, replace `CAL_FINGERPRINT = _calibration_fingerprint()` and `cache_path` (lines 104-110) with:

```python
# Kept as a provenance snapshot of the calibration at import. NOT used for cache
# filenames — see cache_path, which must read the live calibration so a
# calibration_override (experiments/common.py) mints its own filename instead of
# clobbering the baseline cache.
CAL_FINGERPRINT = _calibration_fingerprint()


def cache_path(psilam, fingerprint=None):
    """Cache filename for a given psi_lambda_B at the LIVE calibration.

    The fingerprint is computed at CALL time, not import time. An override applied
    after this module was imported must produce a different filename.
    """
    fp = fingerprint or _calibration_fingerprint()
    tag = f"{psilam:.2f}".replace(".", "p")
    return os.path.join(HERE, f"cache_G_main_v{CACHE_SCHEMA}_psilam{tag}_cal{fp}.npz")
```

- [ ] **Step 6: Stamp the fingerprint into the file and assert it on load**

In `_extract`, add the fingerprint to the returned dict. Change the first line of the function body from:

```python
    out = {"T": np.array(T), "dShock_def_D": np.asarray(dshock), "psi_lambda_B": np.array(psilam)}
```

to:

```python
    out = {"T": np.array(T), "dShock_def_D": np.asarray(dshock), "psi_lambda_B": np.array(psilam),
           "cal_fingerprint": np.array(_calibration_fingerprint())}
```

Then replace `load_cache` (lines 243-248) with:

```python
def load_cache(psilam):
    """Load a cache, asserting it was built under the live calibration.

    cache_path already keys on the live fingerprint, so a mismatch normally shows
    up as a missing file. This second check catches the case where a file was
    hand-copied or renamed — it fails loudly instead of returning another model.
    """
    path = cache_path(psilam)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"No cache at {os.path.basename(path)}. The live calibration has no cache "
            f"built for it — run:  /opt/anaconda3/envs/ssj/bin/python "
            f"diagnostics/regimes/regime_model.py")
    # allow_pickle for backward-compat with caches that stored a stray object entry
    # (cb_flow_D, now excluded); the matrices this module reads are all plain float.
    with np.load(path, allow_pickle=True) as d:
        cache = {k: d[k] for k in d.files if not d[k].dtype == object}
    live = _calibration_fingerprint()
    stored = str(cache["cal_fingerprint"])
    assert stored == live, (
        f"cache fingerprint {stored} != live calibration {live} — stale or "
        f"hand-renamed cache; rebuild with regime_model.py --force")
    return cache
```

- [ ] **Step 7: Make `build_caches` read psi_lambda_B at call time**

In `build_caches`, replace the first line:

```python
    paths = {PSILAM_MAIN: cache_path(PSILAM_MAIN), 0.0: cache_path(0.0)}
```

with:

```python
    # Read live, not from the import-time PSILAM_MAIN: under a calibration_override
    # the two can differ, and the override must win.
    psilam_live = _live_psilam()
    paths = {psilam_live: cache_path(psilam_live), 0.0: cache_path(0.0)}
```

and, further down, replace the two uses of `PSILAM_MAIN` inside the function body (the drift assertion, the `_solve_G` label, and the `paths[PSILAM_MAIN]` / `_extract(...)` calls) with `psilam_live`:

```python
    assert abs(cal["psi_lambda_B_D"] - psilam_live) < 1e-9, (
        f"calibration drifted mid-run: live={cal['psi_lambda_B_D']} vs cache key {psilam_live}")
```

```python
    G28 = _solve_G(model, ss28, unk, tgt, T, f"{psilam_live}")
```

```python
    np.savez_compressed(paths[psilam_live], **_extract(G28, ss28, T, dshock, psilam_live))
```

- [ ] **Step 8: Write `experiments/common.py` (minimum to pass Task 1's tests)**

Create `experiments/common.py`:

```python
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
```

- [ ] **Step 9: Run the tests to verify they pass**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest experiments/test_common.py -v
```

Expected: 6 passed.

If `test_cache_path_reflects_calibration_override` fails with identical filenames, `cache_path` is still reading the import-time constant — recheck Step 5.

- [ ] **Step 10: Verify the existing regimes tests still pass**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest diagnostics/regimes/test_lottery_math.py -v
```

Expected: all pass (unchanged — `lottery_math.py` was not touched).

- [ ] **Step 11: Commit**

```bash
git add experiments/__init__.py experiments/.gitignore experiments/common.py \
        experiments/test_common.py diagnostics/regimes/regime_model.py
git commit --no-verify -m "Cache schema v2: call-time fingerprint, Phi_D/def_rate_D outputs

cache_path computed the fingerprint at import, so a calibration override applied
afterwards (E3) would write to the baseline filename and clobber it. Now computed
at call time, with the fingerprint stamped into the file and asserted on load.

Adds Phi_D and def_rate_D to REQUIRED (E2's identity and E1's off-path expected
loss). T_D goes in OPTIONAL: T0=T1=0 makes it identically zero, so zero-fill is
correct rather than a silent hole, and E2's closure assertion catches it anyway.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

> `--no-verify` is correct here: the doc-sync hook requires STATE/PROGRESS/HANDOFF
> in any commit touching `*.py`, but this is mid-feature plumbing with no results
> yet. **Task 8 updates all three docs before the feature is done** — that is where
> the obligation is discharged, and it is not optional.

---

## Task 2: Rebuild the cache under schema 2

**Files:**
- Runs: `diagnostics/regimes/regime_model.py`

- [ ] **Step 1: Rebuild (long-running, ~20 min)**

```bash
/opt/anaconda3/envs/ssj/bin/python diagnostics/regimes/regime_model.py --force
```

Expected in the output:
- `model-build sanity: G_tpi[cb=0] vs baseline spread_rb max|err| = ...e-1x (expect <1e-8)`
- `cross-check vs SA-1 probe: d(spread_rb)/d(cb_buy)[0,0] = -1.9455e-02` (backstop compresses)
- `caches written: ['cache_G_main_v2_psilam8p50_cal<fp>.npz', 'cache_G_main_v2_psilam0p00_cal<fp>.npz']`

If it raises `REQUIRED outputs missing from main's G_tpi: ['Phi_D']` or `['def_rate_D']`, that output genuinely is not in the TPI model's Jacobian. Do not zero-fill it — move it to `OPTIONAL` only after confirming the value is structurally zero, and record why in the spec. E2's closure assertion is the backstop.

- [ ] **Step 2: Confirm the new keys are present**

```bash
/opt/anaconda3/envs/ssj/bin/python -c "
import sys; sys.path.insert(0,'diagnostics/regimes'); sys.path.insert(0,'code')
from regime_model import load_cache, _live_psilam
c = load_cache(_live_psilam())
for k in ['Phi_D__shock_def_D','def_rate_D__shock_def_D','Phi_D_ss','delta_b_D_ss','q_b_F_ss']:
    print(k, 'OK' if k in c else 'MISSING')
print('T_D present:', 'T_D__shock_def_D' in c)
print('fingerprint:', str(c['cal_fingerprint']))
"
```

Expected: the five `OK` lines. `T_D present:` may be `True` or `False` — both are acceptable (see Task 1 Step 4's note).

- [ ] **Step 3: Record the Walras residuals from the build**

The spec's verification table lists the Walras thresholds, but the experiment runners are pure post-Jacobian numpy and **cannot recompute them** — market clearing is solved during the SS/Jacobian stage inside `build_caches`, not during a matrix multiply. So the Walras check belongs here, at cache-build time, and this is the only place it is meaningful.

Re-run the build capturing output (or scroll back through Step 1's output) and confirm the residuals printed by `build_and_solve`:

```bash
grep -E "goods_mkt_D|goods_mkt_F|ca_res_D|deposit_mkt" diagnostics/regimes/_cache_build.log 2>/dev/null \
  || echo "not logged to file — read Step 1's stdout instead"
```

Expected, against `docs/verification_report.md`: `goods_mkt_D` ≤1e−14, `goods_mkt_F` ≤1e−7, `ca_res_D` ≤1e−7, `deposit_mkt_D/F` ≤1e−13.

Record the four numbers in the Task 8 STATE.md entry. If any exceeds its threshold, **stop** — the cache is built on a model that does not clear, and every downstream result is invalid.

---

## Task 3: E2 — ΔY decomposition

Do E2 **before** E1: it is self-verifying, so it validates the freshly rebuilt cache. If the identity does not close, every other result is suspect.

**Files:**
- Create: `experiments/e2_dy_decomposition.py`
- Test: `experiments/test_e2_identity.py`
- Modify: `experiments/common.py` (add cache/IRF helpers)

- [ ] **Step 1: Write the failing test**

Create `experiments/test_e2_identity.py`:

```python
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
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest experiments/test_e2_identity.py -v
```

Expected: `ModuleNotFoundError: No module named 'e2_dy_decomposition'`.

- [ ] **Step 3: Add cache/IRF helpers to `experiments/common.py`**

Append to `experiments/common.py`:

```python
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
    stance. Targets are spec §7: 25% (medium), 50% (aggressive).
    """
    from lottery_math import gamma_for_compression
    A_def = cache["spread_rb__shock_def_D"]
    A_cb = cache["spread_rb__cb_buy_D"]
    eps = np.asarray(cache["dShock_def_D"])
    assert float(A_cb[0, 0]) < 0.0, (
        f"A_cb[0,0]={float(A_cb[0,0]):+.4e} >= 0: CB purchases WIDEN the spread, so "
        "compression targeting is infeasible. This is the ms-regime SA-1 pathology, "
        "which must be absent on main — investigate before reporting anything.")
    return {
        "passive": 0.0,
        "medium": float(gamma_for_compression(A_def, A_cb, eps, target=0.25)),
        "aggressive": float(gamma_for_compression(A_def, A_cb, eps, target=0.50)),
    }


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
```

- [ ] **Step 4: Write `experiments/e2_dy_decomposition.py`**

```python
"""E2 — decompose the output response into its goods-market components.

Identity, from market_clearing_D (code/equations_D.py:139):

    Y_D = P_CES_D*C_D + I_D + G_D + Phi_D + T_D + NX_D

linearised, with the product rule on the consumption term:

    dY = P_ss*dC + C_ss*dP_CES + dI + dG + dPhi + dT + dNX

goods_mkt_D is a TARGETED residual held to <=1e-14, so this closes to solver
tolerance. That makes the decomposition self-verifying: a non-closing residual
means a missing term, not a small error to tolerate.

Gates two flagged claims:
  * SPEC "Where Germany genuinely benefits" — do not assert the trade channel
    without checking the investment/NX split.
  * STATE.md watch item — Y_D[0] is POSITIVE under both intervening regimes and
    the A5 dY_D trough never goes negative. This shows which component flips.
"""
import numpy as np

from common import load_cache, provenance, regime_irfs, write_results

CLOSURE_TOL = 1e-7

# Display order = economic reading order, not dict order.
COMPONENTS = ["consumption_quantity", "consumption_price", "investment",
              "government", "portfolio_cost", "macropru_tax", "net_exports"]


def decompose_dY(irf, ss):
    """Return ({component: path}, residual). Residual must be ~0 by the identity."""
    T = len(np.asarray(irf["Y_D"]))
    zero = np.zeros(T)

    def get(name):
        # An absent term is structurally zero (G_D is constant and absent from the
        # Jacobian; T_D is identically zero at T0=T1=0). Carried explicitly so the
        # identity stays complete if either is ever switched on.
        return np.asarray(irf[name]) if name in irf else zero

    components = {
        "consumption_quantity": float(ss["P_CES_D_ss"]) * get("C_D"),
        "consumption_price": float(ss["C_D_ss"]) * get("P_CES_D"),
        "investment": get("I_D"),
        "government": get("G_D"),
        "portfolio_cost": get("Phi_D"),
        "macropru_tax": get("T_D"),
        "net_exports": get("NX_D"),
    }
    residual = np.asarray(irf["Y_D"]) - sum(components.values())
    return components, residual


def run(horizon=40):
    cache = load_cache()
    ss = {"P_CES_D_ss": float(cache["P_CES_D_ss"]), "C_D_ss": float(cache["C_D_ss"])}
    Y_ss = float(cache["Y_D_ss"])
    beta = float(cache["beta_D"])
    disc = beta ** np.arange(horizon)

    payload = {"provenance": provenance(), "horizon": horizon,
               "closure_tol": CLOSURE_TOL, "regimes": {}}

    for name, (gamma, irf) in regime_irfs(cache).items():
        components, residual = decompose_dY(irf, ss)
        max_resid = float(np.max(np.abs(residual[:horizon])))
        assert max_resid < CLOSURE_TOL, (
            f"E2 identity does not close for regime '{name}': max|residual| = "
            f"{max_resid:.3e} > {CLOSURE_TOL:.0e}. A term is missing from the "
            f"decomposition — do not report these numbers.")

        dY = np.asarray(irf["Y_D"])[:horizon]
        payload["regimes"][name] = {
            "gamma": gamma,
            "max_abs_residual": max_resid,
            "dY_impact_pct_ss": float(dY[0] * 100.0 / Y_ss),
            "dY_trough_pct_ss": float(dY.min() * 100.0 / Y_ss),
            "dY_pv": float((disc * dY).sum()),
            "components_pv": {k: float((disc * v[:horizon]).sum())
                              for k, v in components.items()},
            "components_impact": {k: float(v[0]) for k, v in components.items()},
            "paths": {k: v[:horizon].tolist() for k, v in components.items()},
            "dY_path": dY.tolist(),
        }

    write_results("e2_dy_decomposition", payload)
    return payload


if __name__ == "__main__":
    res = run()
    print(f"{'regime':>12} {'gamma':>9} {'dY[0] %SS':>11} {'dI PV':>12} "
          f"{'dNX PV':>12} {'dC_q PV':>12} {'max|resid|':>11}")
    print("-" * 84)
    for name, r in res["regimes"].items():
        c = r["components_pv"]
        print(f"{name:>12} {r['gamma']:>9.4f} {r['dY_impact_pct_ss']:>+11.4f} "
              f"{c['investment']:>+12.3e} {c['net_exports']:>+12.3e} "
              f"{c['consumption_quantity']:>+12.3e} {r['max_abs_residual']:>11.2e}")
    print("-" * 84)
    print("Identity: dY = P_ss*dC + C_ss*dP + dI + dG + dPhi + dT + dNX "
          "(market_clearing_D, targeted to <=1e-14)")
```

- [ ] **Step 5: Run the unit tests to verify they pass**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest experiments/test_e2_identity.py -v
```

Expected: 4 passed.

- [ ] **Step 6: Run E2 against the real cache**

```bash
/opt/anaconda3/envs/ssj/bin/python experiments/e2_dy_decomposition.py
```

Expected: a three-row table with `max|resid|` around `1e-16` to `1e-9` in every row.

If the assertion fires, **stop and diagnose** — do not raise `CLOSURE_TOL`. The likely causes, in order: `Phi_D` absent from the cache (Task 2 Step 2 would have shown it), `T_D` non-zero when it was assumed zero, or a sign error in the product-rule term.

- [ ] **Step 7: Record the answer to the watch item**

Read the printed table and note, for the paper: does `dY[0]` turn positive under `medium`/`aggressive`, and if so which component carries it? A positive `dY[0]` driven by `investment` is economics (the backstop cushions investment); one driven by `net_exports` or `consumption_price` is more likely terms-of-trade or linear-rule overshoot at `gamma_aggressive`.

- [ ] **Step 8: Commit**

```bash
git add experiments/common.py experiments/e2_dy_decomposition.py \
        experiments/test_e2_identity.py experiments/results/e2_dy_decomposition.json
git commit --no-verify -m "E2: dY decomposition against the market_clearing_D identity

Self-verifying: goods_mkt_D is a targeted residual (<=1e-14), so the components
must sum to dY to solver tolerance. Asserts closure at 1e-7 and refuses to report
otherwise.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 4: E1 — Core backstop schedule

**Files:**
- Create: `experiments/e1_backstop_schedule.py`

- [ ] **Step 1: Write `experiments/e1_backstop_schedule.py`**

```python
"""E1 — the core backstop schedule at the three named regimes.

Canonical parameterisation is the NAMED REGIMES (passive / medium / aggressive),
with gamma solved for 0/25/50% peak-spread compression, not code/tpi.py's round
gamma in {0,2,5,10}. Solved gammas keep their meaning across recalibrations.

A5-1: the German fiscal object is reported as THREE SEPARATE quantities and never
summed into an "implicit transfer":
  1. exposure          — discounted purchases; what capital-key sharing acts on
  2. expected loss     — priced, computed OFF-PATH (see below)
  3. Greek fiscal saving — pd_D differential vs passive; no pricing assumption

The off-path requirement is not a stylistic choice. Per SPEC's implementation
hazard, the excess-return flow EL_price * def_rate_t * b_ss is a first-order
deviation times a STEADY-STATE level, so it does not vanish along the computed
path: bondholders earn the premium with no offsetting loss, because writeoff is
off and the IRF traces the no-default branch. Reading expected loss off the
realised path therefore shows the CB mechanically profiting. It must be summed
by hand as Sum beta^t * EL_price * def_rate_t * q_b * cb_buy_t.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (BP_ANN, FIGURES_DIR, irf_from_cache, load_cache, pct_of_ss,
                    provenance, regime_irfs, write_results)

T_PNL = 100          # PV horizon, matches code/tpi.py's cb_pnl
T_WELFARE = 100      # discounted welfare horizon, matches run_tpi
PLOT_N = 60

BLUE, RED, GREEN, ORANGE = "#002147", "#8C1515", "#1a6e3a", "#c87941"
REGIME_COLORS = {"passive": RED, "medium": ORANGE, "aggressive": GREEN}


def cb_pnl(irf, cache, T_pnl=T_PNL):
    """PV decomposition of the CB's D-bond position, in D-goods units.

    Ported from code/tpi.py's cb_pnl, reading steady-state levels from the cache
    rather than re-deriving them, so the two cannot drift apart.
    """
    beta_F = float(cache["beta_F"])
    q_b_D_ss = float(cache["q_b_D_ss"])
    q_b_F_ss = float(cache["q_b_F_ss"])
    delta_b_D = float(cache["delta_b_D_ss"])
    EL_price_D = float(cache["EL_price_D"])

    disc = beta_F ** np.arange(T_pnl)
    cb = np.asarray(irf["cb_buy_D"])[:T_pnl]
    cb_l = np.concatenate([[0.0], cb[:-1]])
    dq = np.asarray(irf["q_b_D"])[:T_pnl]
    dq_l = np.concatenate([[0.0], dq[:-1]])
    defr = np.asarray(irf["def_rate_D"])[:T_pnl]
    dspr = np.asarray(irf["spread_rb"])[:T_pnl]

    # SS yield differential. delta_b_F is not cached, so the F leg uses D's
    # duration; that is only valid because SS yields are equalised in this model
    # (carry_ss_pv is 0.0000% in every production run). Asserted below rather than
    # assumed — if it ever fires, add delta_b_F to regime_model.SS_META.
    rb_D_ss = delta_b_D * (1.0 / q_b_D_ss - 1.0)
    rb_F_ss = delta_b_D * (1.0 / q_b_F_ss - 1.0)
    spread_ss = rb_D_ss - rb_F_ss
    assert abs(spread_ss) < 1e-10, (
        f"SS spread {spread_ss:.3e} != 0, so carry_ss_pv is being computed with D's "
        f"duration on both legs. Add delta_b_F to regime_model.SS_META and use it here.")

    purchases = cb - (1.0 - delta_b_D) * cb_l
    return {
        "peak_exposure": float(np.max(q_b_D_ss * cb)),
        "purchases_pv": float((disc * q_b_D_ss * purchases).sum()),
        "el_pv": float((disc * EL_price_D * defr * q_b_D_ss * cb).sum()),
        "prem_pv": float((disc * dspr * q_b_D_ss * cb_l).sum()),
        "carry_ss_pv": float((disc * spread_ss * q_b_D_ss * cb_l).sum()),
        "mtm_pv": float((disc * (1.0 - delta_b_D) * cb_l * (dq - dq_l)).sum()),
    }


def primary_deficit(irf, cache):
    """pd_D = dG - P_CES_ss*dTAX - TAX_ss*dP_CES  (the austerity channel).

    Same construction as run_regimes.irf_all, kept identical on purpose.
    """
    T = len(np.asarray(irf["Y_D"]))
    dG = np.asarray(irf["G_D"]) if "G_D" in irf else np.zeros(T)
    return (dG - float(cache["P_CES_D_ss"]) * np.asarray(irf["TAX_D"])
            - float(cache["TAX_D_ss"]) * np.asarray(irf["P_CES_D"]))


def welfare(irf, cache, T_w=T_WELFARE):
    """Discounted utility deviation, % of quarterly SS consumption. SECONDARY."""
    beta_D, beta_F = float(cache["beta_D"]), float(cache["beta_F"])
    W_D = float((np.asarray(irf["U_D"])[:T_w] * beta_D ** np.arange(T_w) * 100).sum())
    W_F = float((np.asarray(irf["U_F"])[:T_w] * beta_F ** np.arange(T_w) * 100).sum())
    return W_D, W_F


def loading_schedule(cache, gamma_max=30.0, n=60):
    """Loading = premium PV / expected-loss PV over a fine gamma grid.

    THE key figure. The paper's self-extinguishing-premium claim (Live Claim 5) is
    the DECLINE, so the schedule, not any single point, is the object.
    """
    from lottery_math import closed_loop
    A_def, A_cb = cache["spread_rb__shock_def_D"], cache["spread_rb__cb_buy_D"]
    eps = np.asarray(cache["dShock_def_D"])
    gammas = np.linspace(0.0, gamma_max, n)
    loading, peak_bp = np.full(n, np.nan), np.empty(n)
    for i, g in enumerate(gammas):
        spread, cb = closed_loop(A_def, A_cb, eps, float(g))
        irf = irf_from_cache(cache, cb, eps)
        d = cb_pnl(irf, cache)
        peak_bp[i] = float(np.max(spread[:T_PNL])) * BP_ANN
        if d["el_pv"] > 1e-16:
            loading[i] = d["prem_pv"] / d["el_pv"]
    return gammas, loading, peak_bp


def run():
    cache = load_cache()
    Y_ss = float(cache["Y_D_ss"])
    n_ss = float(cache["n_inter_D_ss"])
    K_ss = float(cache["K_D_ss"])
    beta_D = float(cache["beta_D"])

    regimes = regime_irfs(cache)
    payload = {"provenance": provenance(),
               "gamma_selection_rule": "peak-spread compression 0/25/50% (spec section 7); "
                                       "gamma solved, not chosen",
               "welfare_caveat": "SECONDARY. SPEC: do not lead with welfare — it is a "
                                 "delicate decomposition-dependent object and comes out "
                                 "near-exactly zero-sum.",
               "regimes": {}}

    pd_passive = primary_deficit(regimes["passive"][1], cache)
    disc = beta_D ** np.arange(T_PNL)

    for name, (gamma, irf) in regimes.items():
        pnl = cb_pnl(irf, cache)
        W_D, W_F = welfare(irf, cache)
        spread = np.asarray(irf["spread_rb"])
        pd_here = primary_deficit(irf, cache)
        # A5-1 object 3: Greek fiscal saving vs the passive counterfactual.
        fiscal_saving_pv = float((disc * (pd_passive - pd_here)[:T_PNL]).sum())

        payload["regimes"][name] = {
            "gamma": gamma,
            "peak_spread_bp_ann": float(np.max(spread[:T_PNL]) * BP_ANN),
            "impact": {
                "Y_D_pct_ss": float(pct_of_ss(np.asarray(irf["Y_D"])[:1], Y_ss)[0]),
                "C_D_pct_ss": float(np.asarray(irf["C_D"])[0] * 100.0 / float(cache["C_D_ss"])),
                "I_D_pct_ss": float(np.asarray(irf["I_D"])[0] * 100.0 / float(cache["I_D_ss"])),
                "n_inter_D_pct_ss": float(pct_of_ss(np.asarray(irf["n_inter_D"])[:1], n_ss)[0]),
                "K_D_pct_ss": float(pct_of_ss(np.asarray(irf["K_D"])[:1], K_ss)[0]),
            },
            "trough": {
                "Y_D_pct_ss": float(np.asarray(irf["Y_D"])[:T_PNL].min() * 100.0 / Y_ss),
                "n_inter_D_pct_ss": float(np.asarray(irf["n_inter_D"])[:T_PNL].min() * 100.0 / n_ss),
            },
            # ---- A5-1: three separate objects, never summed ----
            "a5_1_exposure_pv_pct_Y": 100.0 * pnl["purchases_pv"] / Y_ss,
            "a5_1_expected_loss_pv_pct_Y": 100.0 * pnl["el_pv"] / Y_ss,
            "a5_1_greek_fiscal_saving_pv": fiscal_saving_pv,
            "peak_exposure_pct_Y": 100.0 * pnl["peak_exposure"] / Y_ss,
            "premium_pv_pct_Y": 100.0 * pnl["prem_pv"] / Y_ss,
            "mtm_pv_pct_Y": 100.0 * pnl["mtm_pv"] / Y_ss,
            "carry_ss_pv_pct_Y": 100.0 * pnl["carry_ss_pv"] / Y_ss,
            "loading": (pnl["prem_pv"] / pnl["el_pv"]) if pnl["el_pv"] > 1e-16 else None,
            "welfare_W_D_secondary": W_D,
            "welfare_W_F_secondary": W_F,
        }

    assert payload["regimes"]["passive"]["a5_1_exposure_pv_pct_Y"] == 0.0, \
        "passive regime must have zero purchases by construction"

    gammas, loading, peak_bp = loading_schedule(cache)
    payload["loading_schedule"] = {"gamma": gammas.tolist(),
                                   "loading": [None if np.isnan(x) else float(x) for x in loading],
                                   "peak_spread_bp_ann": peak_bp.tolist()}

    _plot(payload, gammas, loading, peak_bp, regimes, cache)
    write_results("e1_backstop_schedule", payload)
    return payload


def _plot(payload, gammas, loading, peak_bp, regimes, cache):
    import os
    os.makedirs(FIGURES_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    axes[0].plot(gammas, loading, color=BLUE, lw=2)
    for name, r in payload["regimes"].items():
        if r["loading"] is not None:
            axes[0].plot(r["gamma"], r["loading"], "o", color=REGIME_COLORS[name],
                         ms=8, label=f"{name} (γ={r['gamma']:.2f})")
    axes[0].axhline(1.0, ls="--", lw=0.8, color="gray")
    axes[0].set_xlabel("backstop aggressiveness γ")
    axes[0].set_ylabel("loading = premium PV / expected-loss PV")
    axes[0].set_title("The premium self-extinguishes as the backstop strengthens")
    axes[0].legend(fontsize=8)

    axes[1].plot(gammas, peak_bp, color=BLUE, lw=2)
    for name, r in payload["regimes"].items():
        axes[1].plot(r["gamma"], r["peak_spread_bp_ann"], "o",
                     color=REGIME_COLORS[name], ms=8, label=name)
    axes[1].set_xlabel("backstop aggressiveness γ")
    axes[1].set_ylabel("peak D–F spread (bp, annualised)")
    axes[1].set_title("Spread compression")
    axes[1].legend(fontsize=8)

    p = payload["provenance"]
    fig.suptitle(f"E1 — backstop schedule (psi_lambda_B={p['psi_lambda_B_D']}, "
                 f"{'market-value' if p['mv_rule_D'] else 'par-value'} rule, "
                 f"writeoff={'on' if p['writeoff_enabled_D'] else 'off'}, "
                 f"scope={p['BANK_SCOPE']}, {p['git_sha']})", fontsize=10)
    fig.text(0.5, 0.01,
             "Loading is premium PV over expected-loss PV, with expected loss computed "
             "OFF-PATH. The decline in γ is the self-extinguishing-premium result: the "
             "wedge exists because the marginal holder is balance-sheet constrained, and "
             "the backstop relieves that constraint, so intervention erodes its own "
             "profit source.", ha="center", fontsize=7.5, style="italic", wrap=True)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(os.path.join(FIGURES_DIR, "fig_e1_loading_schedule.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    res = run()
    print(f"{'regime':>12} {'gamma':>9} {'peak bp':>9} {'Y[0] %SS':>10} "
          f"{'n_int[0] %SS':>13} {'exposure %Y':>12} {'EL PV %Y':>10} {'loading':>8}")
    print("-" * 92)
    for name, r in res["regimes"].items():
        ld = "n/a" if r["loading"] is None else f"{r['loading']:.2f}"
        print(f"{name:>12} {r['gamma']:>9.4f} {r['peak_spread_bp_ann']:>9.1f} "
              f"{r['impact']['Y_D_pct_ss']:>+10.4f} {r['impact']['n_inter_D_pct_ss']:>+13.3f} "
              f"{r['a5_1_exposure_pv_pct_Y']:>12.3f} "
              f"{r['a5_1_expected_loss_pv_pct_Y']:>10.4f} {ld:>8}")
    print("-" * 92)
    print("A5-1: exposure / expected loss / Greek fiscal saving are SEPARATE objects — "
          "never sum them into an 'implicit transfer'.")
    print("Welfare is reported in the JSON as SECONDARY only (SPEC: do not lead with it).")
```

- [ ] **Step 2: Run E1**

```bash
/opt/anaconda3/envs/ssj/bin/python experiments/e1_backstop_schedule.py
```

Expected: a three-row table. `passive` should show ≈150.4 bp and zero exposure; `medium` ≈75% of that; `aggressive` ≈50%. Loading should be `n/a` for passive (no purchases, zero denominator) and **declining** from medium to aggressive.

Cross-check against the known-good production run: `code/main.py` at γ=2/5/10 gave loading 4.35/4.01/3.44. The named regimes sit at γ≈5.08/12.73, so `medium` should land near 4.0 and `aggressive` below 3.44.

- [ ] **Step 3: Commit**

```bash
git add experiments/e1_backstop_schedule.py experiments/results/e1_backstop_schedule.json \
        experiments/figures/fig_e1_loading_schedule.png
git commit --no-verify -m "E1: backstop schedule at the three named regimes

Named regimes canonical (gamma solved for 0/25/50% peak-spread compression), with
the continuous loading schedule as the key figure. A5-1's three German objects are
reported separately and never summed. Expected loss is computed off-path, per
SPEC's implementation hazard. Welfare is present but labelled secondary.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 5: E3 — S-1 writeoff variants

Spec review found that `zeta_writeoff` enters the `EL_price` anchor at `code/steady_state.py:107-112` **without** the `writeoff_enabled` gate that guards every other writeoff term. So the two switches are not one switch, and only `writeoff_enabled` is steady-state-neutral.

**Files:**
- Create: `experiments/e3_writeoff_s1.py`

- [ ] **Step 1: Write `experiments/e3_writeoff_s1.py`**

```python
"""E3 — S-1: does the sovereign default produce realised bank losses?

Two nested variants, because the two switches do different things:

                     writeoff_enabled  zeta_writeoff   steady state
  baseline                  0               0.0        --
  E3a  coupon-only          1               0.0        STRICTLY INVARIANT
  E3b  full writeoff        1               1.0        MOVES, via EL_price

In bond_return_D / government_ss_D / bond_price_ss_D / budget_residual_D:

    current_payoff = delta_b * (1 - def_rate*haircut*writeoff_enabled)
    continuation   = (1-delta_b)*q_b * (1 - zeta*def_rate*haircut*writeoff_enabled)

Both legs carry def_rate, which is 0 at SS, so writeoff_enabled is SS-neutral. But
zeta_writeoff ALSO appears in the EL_price anchor (code/steady_state.py:107-112),
and there it is NOT gated by writeoff_enabled:

    EL_price = (1-recovery) * [delta_b + zeta*(1-delta_b)*q_b] / q_b

At recovery=0.30, delta_b_D=0.0777, q_b_D~0.83 that takes EL_price_D from ~0.0655
to ~0.711 — about 10.9x. EL_price is the loading's DENOMINATOR, so this lands
directly on SPEC Live Claim 1. Reported, never re-tuned away.

Recovery stays at 0.30 (EL-1's resolved Greek-PSI NPV value) rather than STATE.md's
older recovery=0.40 suggestion, which predates EL-1 and would move two dials at once.
"""
import numpy as np

from common import calibration_override, load_cache, provenance, write_results

VARIANTS = {
    "e3a_coupon_only": {"writeoff_enabled_D": 1.0, "writeoff_enabled_F": 1.0,
                        "zeta_writeoff_D": 0.0, "zeta_writeoff_F": 0.0},
    "e3b_full": {"writeoff_enabled_D": 1.0, "writeoff_enabled_F": 1.0,
                 "zeta_writeoff_D": 1.0, "zeta_writeoff_F": 1.0},
}

SS_INVARIANT_KEYS = ["q_b_D_ss", "Y_D_ss", "C_D_ss", "I_D_ss", "NX_D_ss",
                     "n_inter_D_ss", "K_D_ss", "TAX_D_ss", "P_CES_D_ss",
                     "b_gov_D_ss", "b_D_D_ss"]


def expected_EL_price(cal, q_b_D):
    """The closed form from code/steady_state.py:107-109."""
    return ((1.0 - cal["recovery_rate_D"])
            * (cal["delta_b_D"] + cal["zeta_writeoff_D"] * (1.0 - cal["delta_b_D"]) * q_b_D)
            / q_b_D)


def build_variant(overrides):
    """Full SS + Jacobian re-solve under an overridden calibration.

    A full re-solve is REQUIRED, not merely safer. The cheap route — patch
    ss.toplevel and re-solve only the Jacobian, as regime_model.build_caches does
    for its psi_lambda_B=0 cache — PRESUMES the SS invariance that E3a exists to
    test, which would make the check circular. E3b genuinely moves the SS, so it
    has no cheap route either.
    """
    from regime_model import build_caches
    with calibration_override(**overrides):
        build_caches()
        return load_cache()


def summarise(cache):
    """The E1 table under one variant."""
    import e1_backstop_schedule as e1
    from common import regime_irfs
    from common import BP_ANN

    Y_ss, n_ss = float(cache["Y_D_ss"]), float(cache["n_inter_D_ss"])
    out = {"EL_price_D": float(cache["EL_price_D"]), "regimes": {}}
    for name, (gamma, irf) in regime_irfs(cache).items():
        pnl = e1.cb_pnl(irf, cache)
        out["regimes"][name] = {
            "gamma": gamma,
            "peak_spread_bp_ann": float(np.max(np.asarray(irf["spread_rb"])[:100]) * BP_ANN),
            "Y_D_impact_pct_ss": float(np.asarray(irf["Y_D"])[0] * 100.0 / Y_ss),
            "n_inter_D_impact_pct_ss": float(np.asarray(irf["n_inter_D"])[0] * 100.0 / n_ss),
            "expected_loss_pv_pct_Y": 100.0 * pnl["el_pv"] / Y_ss,
            "premium_pv_pct_Y": 100.0 * pnl["prem_pv"] / Y_ss,
            "loading": (pnl["prem_pv"] / pnl["el_pv"]) if pnl["el_pv"] > 1e-16 else None,
        }
    return out


def run():
    from calibration import get_calibration

    baseline_cache = load_cache()
    baseline_ss = {k: float(baseline_cache[k]) for k in SS_INVARIANT_KEYS}
    payload = {"provenance": provenance(),
               "baseline": summarise(baseline_cache),
               "variants": {}, "checks": {}}

    for name, overrides in VARIANTS.items():
        cache = build_variant(overrides)
        payload["variants"][name] = summarise(cache)

        with calibration_override(**overrides):
            cal = get_calibration()
        q_b_D = float(cache["q_b_D_ss"])
        el_expected = expected_EL_price(cal, q_b_D)
        el_actual = float(cache["EL_price_D"])
        assert abs(el_actual - el_expected) < 1e-12, (
            f"{name}: EL_price_D={el_actual:.9f} != closed form {el_expected:.9f}. "
            f"code/steady_state.py:107-109 no longer matches this experiment's model "
            f"of it — reconcile before reporting.")

        drift = {k: float(cache[k]) - baseline_ss[k] for k in SS_INVARIANT_KEYS}
        max_drift = max(abs(v) for v in drift.values())
        payload["checks"][name] = {"EL_price_expected": el_expected,
                                   "EL_price_actual": el_actual,
                                   "max_ss_drift": max_drift,
                                   "ss_drift": drift}

        if name == "e3a_coupon_only":
            # writeoff_enabled is multiplied by def_rate_ss = 0 everywhere, and zeta
            # is unchanged, so the SS must be bit-identical. Drift here is a bug.
            assert max_drift < 1e-10, (
                f"E3a moved the steady state (max drift {max_drift:.3e}). "
                f"writeoff_enabled is supposed to be SS-neutral — every writeoff term "
                f"is multiplied by def_rate_ss=0. Investigate before reporting.")
        else:
            # E3b DOES move the SS, through EL_price. Asserting invariance here would
            # be wrong. The EL_price closed-form check above is the check that applies.
            assert max_drift > 0.0 or el_actual != baseline_cache["EL_price_D"], \
                "E3b changed nothing — the zeta override did not reach the SS solve"

    write_results("e3_writeoff_s1", payload)
    return payload


if __name__ == "__main__":
    res = run()
    print(f"\n{'setting':>18} {'EL_price':>10} {'peak bp (passive)':>18} "
          f"{'loading (medium)':>17} {'loading (aggr.)':>16}")
    print("-" * 84)
    rows = [("baseline", res["baseline"])] + list(res["variants"].items())
    for name, r in rows:
        def ld(reg):
            v = r["regimes"][reg]["loading"]
            return "n/a" if v is None else f"{v:.2f}"
        print(f"{name:>18} {r['EL_price_D']:>10.4f} "
              f"{r['regimes']['passive']['peak_spread_bp_ann']:>18.1f} "
              f"{ld('medium'):>17} {ld('aggressive'):>16}")
    print("-" * 84)
    for name, c in res["checks"].items():
        print(f"{name}: EL_price {c['EL_price_actual']:.6f} (closed form "
              f"{c['EL_price_expected']:.6f}), max SS drift {c['max_ss_drift']:.3e}")
    print("\npsi_lambda_B was tuned to 150bp with realised losses OFF. Any overshoot "
          "here is a REPORTABLE FACT about whether the target survives S-1, not a "
          "number to re-tune away.")
```

- [ ] **Step 2: Run E3 (long-running, ~40 min for two full solves)**

```bash
/opt/anaconda3/envs/ssj/bin/python experiments/e3_writeoff_s1.py
```

Expected: `EL_price` ≈0.0655 for baseline and E3a (identical), ≈0.711 for E3b. E3a's `max SS drift` must be `<1e-10`. E3b's loading should fall sharply — possibly below 1.

If E3a's SS-invariance assertion fires, that is a genuine finding about `writeoff_enabled` and must be investigated, not suppressed.

- [ ] **Step 3: Commit**

```bash
git add experiments/e3_writeoff_s1.py experiments/results/e3_writeoff_s1.json
git commit --no-verify -m "E3: S-1 writeoff, split into coupon-only and full variants

zeta_writeoff enters the EL_price anchor ungated by writeoff_enabled
(steady_state.py:107-112), so only writeoff_enabled is SS-neutral. E3a asserts
strict SS invariance; E3b asserts EL_price against its closed form instead, since
it legitimately moves the SS (~10.9x on EL_price, straight onto the loading
denominator).

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 6: Orchestrator and the results document

**Files:**
- Create: `experiments/run_all.py`
- Generates: `docs/experiments_results.md`

- [ ] **Step 1: Write `experiments/run_all.py`**

```python
"""Run every experiment and render docs/experiments_results.md.

E1 and E2 are seconds (post-Jacobian numpy on one cache). E3 is ~40 min because it
re-solves the model twice. Pass --skip-e3 to regenerate the document from existing
results without paying for that.
"""
import json
import os
import sys

from common import RESULTS_DIR, ROOT, provenance

DOC = os.path.join(ROOT, "docs", "experiments_results.md")


def _load(name):
    path = os.path.join(RESULTS_DIR, f"{name}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def _stamp(p):
    return (f"*Generated {p['generated']} from `{p['git_sha']}` · calibration "
            f"`{p['cal_fingerprint']}` · scope `{p['BANK_SCOPE']}` · "
            f"`psi_lambda_B={p['psi_lambda_B_D']}` · "
            f"`mv_rule={p['mv_rule_D']:g}` · `recovery_rate={p['recovery_rate_D']}` · "
            f"`writeoff_enabled={p['writeoff_enabled_D']:g}` · "
            f"`zeta_writeoff={p['zeta_writeoff_D']:g}`*")


def render():
    e1, e2, e3 = _load("e1_backstop_schedule"), _load("e2_dy_decomposition"), _load("e3_writeoff_s1")
    L = ["# Policy experiments — standard results", "",
         _stamp(provenance()), "",
         "Generated by `experiments/run_all.py`. Do not edit by hand — edit the "
         "experiment and re-run. Spec: "
         "`docs/superpowers/specs/2026-08-01-policy-experiments-design.md`.", ""]

    if e1:
        L += ["## E1 — Backstop schedule", "", _stamp(e1["provenance"]), "",
              f"γ selection: {e1['gamma_selection_rule']}.", "",
              "| regime | γ | peak spread (bp ann) | Y_D[0] (% SS) | n_inter_D[0] (% SS) "
              "| exposure PV (% Y) | expected loss PV (% Y) | Greek fiscal saving PV | loading |",
              "|---|---|---|---|---|---|---|---|---|"]
        for name, r in e1["regimes"].items():
            ld = "n/a" if r["loading"] is None else f"{r['loading']:.2f}"
            L.append(f"| {name} | {r['gamma']:.4f} | {r['peak_spread_bp_ann']:.1f} | "
                     f"{r['impact']['Y_D_pct_ss']:+.4f} | {r['impact']['n_inter_D_pct_ss']:+.3f} | "
                     f"{r['a5_1_exposure_pv_pct_Y']:.3f} | {r['a5_1_expected_loss_pv_pct_Y']:.4f} | "
                     f"{r['a5_1_greek_fiscal_saving_pv']:+.5f} | {ld} |")
        L += ["", "**A5-1.** Exposure, expected loss and Greek fiscal saving are three "
              "separate objects. Do not sum them into an 'implicit transfer'. Expected "
              "loss is computed off-path; reading it off the realised path shows the CB "
              "mechanically profiting (SPEC, structural constraints).", "",
              f"**Welfare (secondary).** {e1['welfare_caveat']}", "",
              "| regime | W_D | W_F |", "|---|---|---|"]
        for name, r in e1["regimes"].items():
            L.append(f"| {name} | {r['welfare_W_D_secondary']:+.4f} | "
                     f"{r['welfare_W_F_secondary']:+.4f} |")
        L += ["", "![loading schedule](../experiments/figures/fig_e1_loading_schedule.png)", ""]

    if e2:
        L += ["## E2 — ΔY decomposition", "", _stamp(e2["provenance"]), "",
              "Identity (`market_clearing_D`): "
              "`dY = P_ss·dC + C_ss·dP_CES + dI + dG + dΦ + dT + dNX`. "
              "`goods_mkt_D` is a targeted residual (≤1e−14), so this closes to solver "
              "tolerance — the decomposition is self-verifying.", "",
              "| regime | Y_D[0] (% SS) | Y_D trough (% SS) | dI PV | dNX PV | dC(qty) PV "
              "| dC(price) PV | dΦ PV | max\\|residual\\| |",
              "|---|---|---|---|---|---|---|---|---|"]
        for name, r in e2["regimes"].items():
            c = r["components_pv"]
            L.append(f"| {name} | {r['dY_impact_pct_ss']:+.4f} | {r['dY_trough_pct_ss']:+.4f} | "
                     f"{c['investment']:+.3e} | {c['net_exports']:+.3e} | "
                     f"{c['consumption_quantity']:+.3e} | {c['consumption_price']:+.3e} | "
                     f"{c['portfolio_cost']:+.3e} | {r['max_abs_residual']:.2e} |")
        L.append("")

    if e3:
        L += ["## E3 — S-1 writeoff", "", _stamp(e3["provenance"]), "",
              "`writeoff_enabled` is steady-state-neutral (every writeoff term is "
              "multiplied by `def_rate_ss = 0`). `zeta_writeoff` is **not** — it enters "
              "the `EL_price` anchor at `code/steady_state.py:107-112` ungated, and "
              "`EL_price` is the loading's denominator.", "",
              "| setting | EL_price_D | peak spread, passive (bp ann) | loading (medium) "
              "| loading (aggressive) |", "|---|---|---|---|---|"]
        for name, r in [("baseline", e3["baseline"])] + list(e3["variants"].items()):
            def ld(reg):
                v = r["regimes"][reg]["loading"]
                return "n/a" if v is None else f"{v:.2f}"
            L.append(f"| {name} | {r['EL_price_D']:.4f} | "
                     f"{r['regimes']['passive']['peak_spread_bp_ann']:.1f} | "
                     f"{ld('medium')} | {ld('aggressive')} |")
        L += ["", "| variant | EL_price (closed form) | EL_price (solved) | max SS drift |",
              "|---|---|---|---|"]
        for name, c in e3["checks"].items():
            L.append(f"| {name} | {c['EL_price_expected']:.6f} | {c['EL_price_actual']:.6f} | "
                     f"{c['max_ss_drift']:.3e} |")
        L += ["", "`psi_lambda_B = 8.5` was tuned to 150 bp with realised losses **off**. "
              "Any overshoot above is a reportable fact about whether that target survives "
              "S-1 — not a number to re-tune away.", ""]

    with open(DOC, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return DOC


def main():
    import e1_backstop_schedule
    import e2_dy_decomposition

    e2_dy_decomposition.run()     # first: self-verifying, validates the cache
    e1_backstop_schedule.run()
    if "--skip-e3" not in sys.argv:
        import e3_writeoff_s1
        e3_writeoff_s1.run()
    print(f"Wrote {render()}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Export `ROOT` from `common.py`**

`run_all.py` imports `ROOT`. It is already defined at module level in `common.py` (Task 1 Step 8), so no change is needed — but verify:

```bash
/opt/anaconda3/envs/ssj/bin/python -c "
import sys; sys.path.insert(0,'experiments')
from common import ROOT, RESULTS_DIR; print(ROOT); print(RESULTS_DIR)"
```

Expected: the repo root and `experiments/results`.

- [ ] **Step 3: Render the document from existing results**

```bash
/opt/anaconda3/envs/ssj/bin/python experiments/run_all.py --skip-e3
```

Expected: `Wrote /Users/Adam/.../docs/experiments_results.md`. Open it and confirm every table is populated and the provenance stamp shows the live calibration.

- [ ] **Step 4: Commit**

```bash
git add experiments/run_all.py docs/experiments_results.md
git commit --no-verify -m "Orchestrator + generated results document

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Task 7: Full-suite verification

**Files:** none modified — this task only runs things.

- [ ] **Step 1: Run the whole test suite**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest experiments/ diagnostics/regimes/ code/test_eba_calibration.py -v
```

Expected: all pass. `code/test_eba_calibration.py` is 10 tests, `test_lottery_math.py` is 7, plus the 10 added here.

- [ ] **Step 2: Confirm the production pipeline is unbroken**

`code/main.py` must be untouched by all of this. It shares `regime_model.py`'s equation files but not its cache, so the only risk is an accidental edit.

```bash
git diff main --stat -- code/
```

Expected: `code/` shows **no changes at all** on this branch.

- [ ] **Step 3: Re-run the production regression**

```bash
/opt/anaconda3/envs/ssj/bin/python code/main.py 2>&1 | tail -40
```

Expected, unchanged from the 2026-08-01 baseline: `K_D = 10.800`, `ca_res_D ≈ 6.9e-17`, `n_inter_D[0] = -3.3804% of SS`, `Y_D[0] = -0.0149% of SS`, peak spread `+0.376 pp`, loading `4.35 / 4.01 / 3.44`. Takes ~20 min.

- [ ] **Step 4: Cross-check E1 against the production numbers**

The named regimes sit at γ≈5.08 and γ≈12.73; `code/tpi.py` reports γ=5 and γ=10. So E1's `medium` loading should be very close to production's γ=5 value (4.01), and `aggressive` should be below production's γ=10 value (3.44). Confirm both, and confirm E1's `passive` peak spread matches production's γ=0 peak of 150.4 bp to within rounding.

If they disagree by more than ~1%, the two code paths have diverged — stop and reconcile before reporting anything.

---

## Task 8: Documentation sync and final commit

The doc-sync hook (`.claude/hooks/require-docs-before-commit.sh` and `.githooks/pre-commit`) requires `docs/STATE.md`, `docs/PROGRESS.md` and `docs/HANDOFF.md` to be updated in any commit staging `code/**` or `*.py`. Earlier tasks used `--no-verify` for mid-feature plumbing; this task discharges the obligation.

**Files:**
- Modify: `docs/STATE.md`, `docs/PROGRESS.md`, `docs/HANDOFF.md`, `CLAUDE.md`

- [ ] **Step 1: Add a STATE.md section**

Insert immediately after the `# Project State` header block, so it is the first thing a reader sees:

```markdown
## Policy experiments (2026-08-03) — `experiments/` package

Standard results set for the paper, on branch `experiments`. Runner:
`/opt/anaconda3/envs/ssj/bin/python experiments/run_all.py` (add `--skip-e3` to
re-render the document without the two ~20-min variant solves). Output:
`docs/experiments_results.md`, provenance-stamped; machine-readable copies in
`experiments/results/*.json`.

- **E1 — backstop schedule.** Named regimes are now canonical (passive / medium /
  aggressive, γ *solved* for 0/25/50% peak-spread compression), replacing
  `code/tpi.py`'s round γ ∈ {0,2,5,10}. A5-1's three German objects — exposure,
  priced expected loss (off-path), Greek fiscal saving — are reported separately
  and must never be summed. Welfare is present but labelled secondary per SPEC.
- **E2 — ΔY decomposition.** Splits the output response against the
  `market_clearing_D` identity. Self-verifying: `goods_mkt_D` is targeted to
  ≤1e−14, so the components must sum to `dY`; the runner asserts closure at 1e−7
  and refuses to report otherwise. Settles SPEC's investment/NX gate and the
  positive-`Y_D[0]`-under-intervention watch item.
- **E3 — S-1.** **`zeta_writeoff` is not steady-state-neutral.** It enters the
  `EL_price` anchor at `code/steady_state.py:107-112` *without* the
  `writeoff_enabled` gate that guards every other writeoff term, taking
  `EL_price_D` from ≈0.0655 to ≈0.711 (~10.9×) at `zeta=1`. `EL_price` is the
  loading's denominator, so this lands on Live Claim 1. S-1 is therefore two
  nested variants: **E3a** (`writeoff_enabled=1`, SS strictly invariant) and
  **E3b** (`+zeta_writeoff=1`, SS moves via `EL_price`).
- **Cache schema v2.** `regime_model.cache_path` now computes the calibration
  fingerprint at *call* time, not import time. It was an import-time constant, so
  E3's calibration override would have written its cache to the baseline filename
  and silently clobbered it. The fingerprint is also stamped into the file and
  asserted on load. New cached outputs: `Phi_D`, `def_rate_D` (`T_D` is optional —
  identically zero at `T0=T1=0`).

Open: the ungated `zeta_writeoff` in `EL_price` may be a latent inconsistency in
its own right — `EL_price` prices an expected loss whose realised counterpart is
switched off. Inert today (`zeta=0`), so nothing published is affected. Not
changed in this pass.
```

- [ ] **Step 2: Add a PROGRESS.md changelog entry**

Add at the top of the changelog, matching the file's existing dated-entry format:

```markdown
## 2026-08-03 — `experiments/` package: standard policy results (E1/E2/E3)

New `experiments/` package producing the paper's standard results set, on top of
`diagnostics/regimes/regime_model.py`'s solve/cache layer. `code/main.py`
untouched — it stays the regression path.

- E1 backstop schedule at named regimes (γ solved for 0/25/50% compression);
  A5-1's three German objects reported separately; loading schedule as key figure.
- E2 ΔY decomposition against the `market_clearing_D` identity, self-verifying
  against the targeted `goods_mkt_D` residual.
- E3 S-1 split into E3a (coupon-only, SS-invariant) and E3b (full). Found during
  spec review: `zeta_writeoff` enters `EL_price` ungated by `writeoff_enabled`
  (`steady_state.py:107-112`), so it is not SS-neutral — ≈10.9× on `EL_price`,
  which is the loading's denominator.
- Cache schema v2: call-time fingerprint (an import-time constant meant E3's
  override would have clobbered the baseline cache), fingerprint asserted on load,
  `Phi_D` and `def_rate_D` added to the cached outputs.
```

- [ ] **Step 3: Update HANDOFF.md**

Add under "Where to start", after the existing bullets:

```markdown
- **Policy experiments:** `experiments/` (branch `experiments`). Run
  `/opt/anaconda3/envs/ssj/bin/python experiments/run_all.py`; `--skip-e3` skips
  the two ~20-min variant solves. Results land in `docs/experiments_results.md`
  and `experiments/results/*.json`. Design spec:
  `docs/superpowers/specs/2026-08-01-policy-experiments-design.md`.
```

- [ ] **Step 4: Update CLAUDE.md's docs reference table**

Add two rows to the table under "## Docs reference":

```markdown
| `docs/experiments_results.md` | **Generated** — standard policy results (E1 backstop schedule, E2 ΔY decomposition, E3 S-1 writeoff). Do not hand-edit; re-run `experiments/run_all.py` |
| `docs/superpowers/specs/2026-08-01-policy-experiments-design.md` | Design spec for the `experiments/` package |
```

- [ ] **Step 5: Verify the hook passes without `--no-verify`**

```bash
git add -A
git commit -m "Policy experiments: docs sync for the experiments package

Records the E3 finding (zeta_writeoff is not SS-neutral; it enters the EL_price
anchor ungated by writeoff_enabled, ~10.9x on the loading denominator) and the
cache-schema-v2 fingerprint fix in STATE, PROGRESS and HANDOFF.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

Expected: the commit succeeds. If the hook blocks it, one of STATE.md / PROGRESS.md / HANDOFF.md was not actually modified — check `git diff --cached --name-only` and fix rather than reaching for `--no-verify`.

- [ ] **Step 6: Confirm the branch is clean**

```bash
git status --short && git log --oneline main..experiments
```

Expected: clean tree, and a commit list covering the spec, cache schema, E2, E1, E3, orchestrator and doc sync.

---

## Notes for the implementer

**Never use bare `python`.** The base Anaconda env has a broken `liblapack` symlink that fails *silently* — you get numbers, they are wrong. Always `/opt/anaconda3/envs/ssj/bin/python`.

**Do not raise a tolerance to make an assertion pass.** Every threshold here traces to `docs/verification_report.md` or to a targeted solver residual. This project has twice been bitten by a solver converging to machine-zero residuals on a degenerate model (C-1, then GK-1 — negative `lambda_gk` with every check passing). An assertion firing is the system working.

**Long-running steps** (Task 2 Step 1, Task 5 Step 2, Task 7 Step 3) are ~20 min each. Run them in the background and do other work rather than blocking.

**If a result contradicts `docs/STATE.md`,** do not adjust the experiment to agree. Report the discrepancy — STATE.md is a record of what was true at a calibration, and the calibration has moved several times.
