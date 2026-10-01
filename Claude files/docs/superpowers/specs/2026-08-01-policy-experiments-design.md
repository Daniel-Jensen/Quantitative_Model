# Preliminary policy experiments — standard results for the paper

**Date:** 2026-08-01
**Branch:** `experiments` (from `main` @ `607b804`)
**Status:** Approved (design), building
**Calibration at design time:** EBA live, `BANK_SCOPE="broad"`, `psi_lambda_B=8.5`,
fingerprint `de195df2`

## Objective

Produce the **standard results set** the paper's claims rest on, as one reproducible,
citable artifact — not as numbers scattered across run logs. Three experiments:

1. **E1 — core backstop schedule.** The γ-schedule results: spread compression, the
   loading (self-extinguishing premium), and incidence. SPEC Live Claims 1 and 5.
2. **E2 — ΔY decomposition (I vs NX).** Splits the output response into its
   goods-market components. Gates SPEC's trade-channel claim.
3. **E3 — S-1 (writeoff).** The one genuinely structural open author decision. Spec
   review split it into two nested variants, `writeoff_enabled` and `zeta_writeoff`,
   because only the first is steady-state-neutral.

Explicitly **out of scope this pass** (author decision, 2026-08-01): distributional
incidence by decile (DIST-1), and the wide robustness grid over `delta_b` / `theta` /
`phi_own` / `def_scale` / `phi_lamb`. Deferred, not rejected.

## Baseline verified before design

Full `code/main.py` run on `607b804`, exit 0, **20m21s wall clock**. Reproduces
`docs/STATE.md` exactly:

| Check | Value |
|---|---|
| `K_D` / `K_F` (target 10.8) | 10.800 / 10.832 |
| `ca_res_D` | 6.85e−17 |
| `goods_mkt_D` / `goods_mkt_F` | −4.25e−07 / −4.19e−07 |
| `n_inter_D[0]` | −3.3804% of SS |
| `Y_D[0]` | −0.0149% of SS |
| peak spread (γ=0) | +0.376 pp = 150.4 bp ann |
| loading at γ=2/5/10 | 4.35 / 4.01 / 3.44 |

The 20-minute figure supersedes the "~3 min per Jacobian" note in CLAUDE.md, which
timed the Jacobian solve alone, not the pipeline.

## Architecture

New `experiments/` package. Production `code/main.py` is **untouched** — CLAUDE.md
makes it the regression test, and a 20-minute variant solve inside it would make that
unusable.

```
experiments/
  common.py                cache access, named-regime γ solve, units helpers
  e1_backstop_schedule.py
  e2_dy_decomposition.py
  e3_writeoff_s1.py
  run_all.py               orchestrator → docs/experiments_results.md
  results/*.json           machine-readable, provenance-stamped
  figures/*.png
```

**Solve/cache layer is reused, not reimplemented.** `diagnostics/regimes/regime_model.py`
already builds and caches per-output response matrices from the production equation
files, keyed by a hash of the whole live calibration. `experiments/common.py` imports
it. No second copy of the model exists anywhere in this package — this is the specific
failure that made the retired `audit_artifacts/` harness test a different model than
`main.py` for weeks.

**Units discipline.** SSJ IRFs are *level* deviations. `×100` is a percentage only
where the SS level is ≈1. `Y_D_ss≈1` passes; `n_inter_D_ss=2.138` and `K_D_ss=10.8` do
not. Every reported percentage divides by its own SS level and is labelled `% of SS`.
Quarterly rate deviations annualise as `×4×1e4` for bp.

### Cache schema change (required)

E2 needs `Phi_D` (portfolio adjustment cost) and `T_D` (macroprudential bond tax) in
the cache; neither is currently cached. Both go in `REQUIRED`, and `Phi_D_ss`/`T_D_ss`
in `SS_META`.

**Hazard:** `cache_path()` keys the filename on the *calibration* fingerprint only, so
after this change an existing cache file would be reloaded under its old name while
missing the new keys — silently, since `irf_all` discovers outputs by scanning cache
keys. Fix: add a schema version constant to the fingerprint input, so any change to
the cached output list mints a new filename. This is the same class of bug the
fingerprint was introduced to prevent.

### Calibration overrides (E3) and the import-time fingerprint

E3 runs the model at a modified calibration, but `get_calibration()` takes no
arguments — every consumer (`steady_state`, `regime_model`, …) calls it directly. The
minimal-change mechanism is a context manager in `experiments/common.py` that patches
`calibration.get_calibration` to return an overridden dict, which all downstream
callers then pick up at solve time.

**Hazard, and it is the dangerous one.** `regime_model` computes `CAL_FINGERPRINT` and
`PSILAM_MAIN` as **module-level constants at import time**. An override applied after
that import would not change them, so E3's cache would be written to the *baseline*
filename — silently colliding with and overwriting the baseline cache. That is exactly
the silent-stale-cache failure the fingerprint exists to prevent, reintroduced one
level up.

Fix, both parts required:

1. `cache_path()` computes the fingerprint **at call time**, not from an import-time
   constant, so an override always mints the correct filename.
2. `load_cache()` asserts the fingerprint of the live calibration matches the file it
   is loading, so a mismatch fails loudly instead of returning the wrong model.

`experiments/common.py` additionally applies overrides *before* importing
`regime_model`, so both mechanisms agree even if (1) is ever regressed.

## E1 — Core backstop schedule

**Canonical parameterisation: named regimes** (author decision, 2026-08-01), replacing
`code/tpi.py`'s round γ ∈ {0,2,5,10}.

| Regime | Definition |
|---|---|
| passive | γ = 0 (fixed anchor) |
| medium | γ solved for 25% peak-spread compression |
| aggressive | γ solved for 50% peak-spread compression |

γ is *solved*, not chosen, so it re-derives itself after any recalibration instead of
silently changing meaning. `gamma_for_compression` in `diagnostics/regimes/lottery_math.py`
already implements this. At the design-time calibration these are γ = 0 / 5.0798 /
12.7260.

Round-γ results stay available underneath for continuity with existing figures and with
SPEC's currently-quoted numbers, but are not the reporting basis.

**Reported per regime:** peak spread (bp ann); impact and trough of `Y_D`, `C_D`,
`I_D`, `NX_D`, `n_inter_D` (% of SS); discounted CB purchases.

**A5-1: the German fiscal object is reported as three separate quantities, never
summed into an "implicit transfer".**

1. **Exposure** — discounted purchases `Σ β^t q_b cb_t`. What capital-key sharing acts on.
2. **Priced expected loss** — `Σ β^t · EL_price · def_rate_t · q_b · cb_buy_t`,
   computed **off-path**. Per SPEC's implementation hazard, this must never be read off
   the realised path, which mechanically shows the CB profiting: the excess-return flow
   is a first-order deviation times a steady-state level and so does not vanish.
3. **Greek fiscal saving** — `pd_D` differential vs the passive counterfactual. The
   cleanest headline; needs no pricing assumption.

Discount rate `beta_F` (creditor side) is retained, and the choice is recorded
explicitly in the results doc rather than left as an unexamined default — A5-1 asks
for exactly this.

**Key figure: the continuous loading schedule** over a fine γ grid, with the three
named regimes marked. Loading = premium PV / expected-loss PV. The paper's
self-extinguishing-premium claim (Live Claim 5) is the *decline*, so the schedule, not
any single point, is the object.

**Welfare is secondary, and labelled so.** ΔW_D / ΔW_F is computed and reported, with
SPEC's "do not lead with this" caveat attached inline in the results doc. At the
design-time calibration it is a near-exact zero-sum transfer (γ=10: +0.0497 / −0.0522),
which is a decomposition-sensitive object.

## E2 — ΔY decomposition (I vs NX)

`market_clearing_D` (`code/equations_D.py:139`) gives the exact identity

```
Y_D = P_CES_D·C_D + I_D + G_D + Phi_D + T_D + NX_D
```

Linearised, with the product rule on the consumption term:

```
dY = P_CES_ss·dC + C_ss·dP_CES + dI + dG + dPhi + dT + dNX
```

**The decomposition verifies itself.** `goods_mkt_D` is a *targeted* residual held to
≤1e−14 (CLAUDE.md acceptance thresholds), so the components must sum to `dY` to solver
tolerance. `run_all` asserts closure at ≤1e−7 and fails loudly otherwise — a
non-closing decomposition means a missing term, not a small error to tolerate.

Known term values at this calibration: `dG_D = 0` (`G_D` is constant, absent from the
Jacobian — already noted in `run_regimes.py`); `T_D ≡ 0` since `T0=T1=0`. Both are
carried explicitly anyway so the identity is complete and stays correct if either is
switched on later.

**Reported:** contribution paths over the plot horizon, and PV shares, per regime.

**Two flagged gates this closes:**

- SPEC "Where Germany genuinely benefits": *"Do not lead with ΔY. A small headline
  output number can be only small because two large channels (investment contraction,
  NX cushion) are netting out — and they land on different households. Check the
  current model's investment/NX decomposition before asserting this."*
- STATE.md watch item: `Y_D[0]` is **positive** under both intervening regimes and the
  A5 `dY_D` trough never goes negative. The decomposition shows *which component* flips
  sign, which distinguishes linear-rule overshoot at γ_aggressive=12.7 from economics.
  If the answer is overshoot, that is a reportable limit on the intervening-regime
  output paths, not a result to publish.

## E3 — S-1 (writeoff)

Recovery stays at **0.30** throughout, EL-1's resolved Greek-PSI NPV value, rather than
STATE.md's older `recovery=0.40, zeta=1.0` suggestion — that line predates EL-1's
resolution, and moving two dials at once would confound the comparison.

### The two switches are not one switch (found during spec review, 2026-08-01)

STATE.md's suggested resolution bundles `writeoff_enabled=1` with `zeta_writeoff=1.0`.
Reading the code, these do different things and only one of them is SS-neutral.

In `bond_return_D` / `government_ss_D` / `bond_price_ss_D` / `budget_residual_D`
(`code/equations_D.py`), the writeoff enters two legs:

```
current_payoff = delta_b * (1 - def_rate*haircut*writeoff_enabled)          # coupon leg
continuation   = (1-delta_b)*q_b * (1 - zeta*def_rate*haircut*writeoff_enabled)  # principal leg
```

Both legs carry `def_rate`, which is 0 at SS, so **`writeoff_enabled` is strictly
SS-neutral.** But `zeta_writeoff` *also* appears in the `EL_price` anchor at
`code/steady_state.py:107-112`, and there it is **not gated by `writeoff_enabled`**:

```
EL_price = (1-recovery) * [delta_b + zeta*(1-delta_b)*q_b] / q_b
```

So `zeta_writeoff` moves the steady state whether or not writeoff is enabled. At
design-time values (`recovery=0.30`, `delta_b_D=0.0777`, `q_b_D≈0.83`):

| `zeta_writeoff` | `EL_price_D` |
|---|---|
| 0.0 (live) | ≈ 0.0655 |
| 1.0 | ≈ 0.711 |

**≈10.9×.** `EL_price` is the denominator of the loading, so this lands directly on
Live Claim 1. It is a first-order change to a headline object, not a robustness wiggle.

(The comment at `steady_state.py:106` — "Enters only the bond FOCs → SS-neutral" — is
about where `EL_price` *enters*, which is correct. It should not be read as saying
`zeta_writeoff` is SS-neutral. Flagged here; not changed in this pass.)

### Consequently, two nested variants

| | `writeoff_enabled` | `zeta_writeoff` | SS |
|---|---|---|---|
| baseline | 0 | 0.0 | — |
| **E3a** coupon-only writeoff | 1 | 0.0 | **strictly invariant** |
| **E3b** full writeoff | 1 | 1.0 | moves via `EL_price` |

This separates two questions S-1 had conflated: *does default produce realised losses
at all* (E3a), and *how much of the bond's value is written off* (E3b).

**E3a assertion: SS strictly identical to baseline.** Every writeoff term is
multiplied by `def_rate_ss = 0`, and `zeta` is unchanged, so any SS drift is a bug and
must halt the experiment.

**E3b assertion: `EL_price` matches the closed form above** to 1e−12, and no *other*
SS quantity moves except through it. The SS-invariance assertion does **not** apply
here and asserting it would be wrong.

**Both variants require a full SS + Jacobian re-solve.** The cheap route — patch
`ss.toplevel` and re-solve only the Jacobian, as `regime_model.build_caches` does for
its `psi_lambda_B=0` cache — is not available: it *presumes* the SS invariance that
E3a exists to test, making the check circular. E3b genuinely moves the SS, so it has
no cheap route either.

**Reported:** the E1 table under all three settings, side by side.

**Expected finding, stated in advance so it is not mistaken for a bug:**
`psi_lambda_B=8.5` was tuned to hit 150bp with realised losses off. Turning them on
gives the default loading another channel, so the peak spread should overshoot 150bp,
and E3b's ~10.9× `EL_price` should compress the loading sharply — plausibly below 1,
which would invert Live Claim 1's over-compensation result *under that setting*. Both
outcomes are reportable facts about whether the calibration target and the headline
claim survive S-1. They must be reported as such and **not silently re-tuned away**;
whether to re-tune `psi_lambda_B` is a separate author decision this result informs.

## Verification

Every runner asserts rather than warns; a failure raises.

| Check | Threshold | Source |
|---|---|---|
| `goods_mkt_D` | ≤1e−14 | CLAUDE.md |
| `goods_mkt_F` | ≤1e−7 | CLAUDE.md |
| `ca_res_D` | ≤1e−7 | CLAUDE.md |
| `deposit_mkt_D/F` | ≤1e−13 | CLAUDE.md |
| E2 identity closure | ≤1e−7 | this spec |
| passive regime purchases | exactly 0 | `run_regimes.py` precedent |
| `n_inter_D[0]`, `Y_D[0]` signs (γ=0) | both negative | CLAUDE.md step 4 |
| E3a SS invariance | SS identical to baseline | this spec |
| E3b `EL_price` | matches closed form, ≤1e−12 | this spec |
| cache fingerprint on load | matches live calibration | this spec |

## Outputs

`docs/experiments_results.md`, generated by `run_all.py`, carrying a provenance stamp
— calibration fingerprint, git SHA, date, and the live values of `psi_lambda_B`,
`mv_rule`, `recovery_rate`, `writeoff_enabled`, `BANK_SCOPE` — at the head of every
table. Read from the live calibration, never hardcoded: `run_regimes.py` shipped a
hardcoded `"market-value rule"` suptitle while running at `mv_rule=0`, which is the
mistake this stamp exists to prevent.

Machine-readable `experiments/results/*.json` alongside, so downstream figure or paper
tooling never re-parses markdown.

## Compute budget

| Step | Cost |
|---|---|
| cache rebuild (schema change, 2 psi_lambda_B points) | ~20 min, background |
| E1 | seconds (post-Jacobian numpy) |
| E2 | seconds (post-Jacobian numpy) |
| E3a (full SS + Jacobian, `writeoff_enabled=1`) | ~20 min, background |
| E3b (full SS + Jacobian, `+ zeta_writeoff=1`) | ~20 min, background |

~60 min total, nearly all unattended. E3 grew from one variant to two when spec review
separated the switches; the cheap patch-the-SS route is unavailable for either (see E3).

## Doc-sync obligation

`.claude/hooks/require-docs-before-commit.sh` and `.githooks/pre-commit` block any
commit staging `code/**` or `*.py` unless `docs/STATE.md`, `docs/PROGRESS.md` and
`docs/HANDOFF.md` are updated in the same commit. `experiments/*.py` is covered by the
`*.py` rule, so the three living docs must be updated in the same commit as the code.
