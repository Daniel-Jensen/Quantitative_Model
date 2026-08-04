# PROGRESS — Changelog

Reverse-chronological development log for the two-country MU-HANK sovereign-risk
model (Greece 2010–12 / ECB TPI). **Newest first.** Reverse-engineered from the
git history (135 commits, first commit 2026-05-14), `docs/STATE.md`,
`docs/audit.md`, `docs/eba_calibration.md`, and `docs/HANDOFF.md` — the detailed
derivations live in those docs; this file is the timeline.

Convention: `[hash]` is a commit on `main`; dates are commit dates. Doc-sync
policy: every commit touching model/code (`code/**`, `*.py`) gets an entry here.
Enforced by `.claude/hooks/require-docs-before-commit.sh` (Claude Code commits)
and `.githooks/pre-commit` (terminal commits; enable with
`git config core.hooksPath .githooks`).

---

## 2026-08-04 — S-1 resolved (`writeoff_enabled=0`); first-draft figures and tables  [this commit]

**S-1 resolved by author decision: the paper keeps `writeoff_enabled=0`**, the
pure risk-premium framing. E3 is retained as an appendix robustness result, and
its content becomes a stated caveat rather than an open question — the
over-compensation claim is *conditional on no realised principal writedown*.

`experiments/paper_outputs.py` emits five captioned figures (`experiments/paper/`)
and three tables (`docs/paper_draft_results.md`): a calibration/identification
ledger separating measured from targeted from free parameters, a moment-match
table, and the main results table. Every number is derived live from the solved SS
or the cached response matrices; only source citations are literal text. Captions
are baked into each PNG rather than living in the LaTeX, so they survive reuse.

**Corrected: the default-loading split is 3.1% fundamental expected loss / 96.9%
collateral friction** (`EL_price_D=0.056134`, `psi_spread_D=1.737724`), not the
10.9%/89% recorded in CLAUDE.md and STATE.md at the pre-EBA calibration. This is
the quantitative core of the constrained-seller argument and 96.9% is a materially
stronger version of it.

**New finding — the backstop damps the oscillation, it does not shift the spread
path down.** Building the transmission figure surfaced that regime differences are
concentrated at impact: by quarter four the net-worth and investment paths have
converged, and the spread ordering *reverses* (t=8: passive +0.5bp vs aggressive
+15.5bp; t=12: passive −10.1bp vs aggressive +5.9bp). Passive overshoots downward
later; intervention decays monotonically. An earlier caption claimed the backstop
"cushions each link roughly proportionally", which is true only on impact — it was
corrected rather than shipped.

Figure palette re-validated with the dataviz validator: the project's `#8C1515`
failed the lightness band and `#002147` failed the chroma floor (reads gray). The
paper set uses `#1B6CA8 / #A62B22 / #c87941 / #1a6e3a`, all passing.

## 2026-08-03 — `experiments/` package complete; production regression bit-identical

Consolidation of the E1/E2/E3 work. `CLAUDE.md` gains an `experiments/`
architecture section, three doc-reference rows, and a rewritten **S-1** row
carrying E3's numbers — S-1 is no longer an open flag but a quantified decision
that determines whether the paper's central claim holds.

**`CLAUDE.md`'s `EL_price_D = 0.0717` was stale and is corrected to 0.056134.**
That figure predates the EBA `delta_b = 0.0777` / `q_b = 0.969`. It is the TPI
loading's denominator, so the error propagated into every loading figure quoted
from memory rather than re-derived. E3 made it load-bearing.

**Production regression re-run after all experiment work** (`code/main.py`,
exit 0, 18m46s): **bit-identical to the pre-work baseline** —
`goods_mkt_D = -4.2493506589857954e-07`, `ca_res_D = 6.852157730108388e-17`,
`K_D = 10.800`/`K_F = 10.832`, `n_inter_D[0] = -3.3804%`, `Y_D[0] = -0.0149%`,
`b_gov_D[499] = 1.4e-05`, loading 4.35/4.01/3.44. `git diff main -- code/` is
empty. The package reads the production equation files but writes nothing into
them, so `code/main.py` remains the regression path. Full suite: 31 passed.

**Independent cross-validation of E1.** Two separately written code paths agree:
`code/main.py` gives loading 4.01 at γ=5, E1 gives 4.00 at γ=5.0798; production
3.44 at γ=10, E1 3.17 at γ=12.726 (correctly lower). Peak spreads hit the
compression targets exactly.

**Two items now blocking paper text, both author decisions:**
1. **S-1** — `zeta_writeoff=1` inverts Live Claim 1 (loading 0.37/0.28, below 1).
2. **A5-1's third object is misnamed** — the reported
   `Σ β^t (pd_passive − pd_intervention)` is negative *because* the backstop
   relaxes austerity, so negative means Greece is better off, the opposite of what
   "fiscal saving" implies. Flip the sign or rename it before it is quoted.

## 2026-08-03 — E3: S-1 writeoff. Full writeoff INVERTS Live Claim 1

`experiments/e3_writeoff_s1.py`. S-1 resolved into two nested variants, because
`writeoff_enabled` and `zeta_writeoff` do different things and only the first is
steady-state-neutral.

**E3a (coupon-only, `zeta=0`) is negligible**: peak spread 150.3 → 149.1 bp,
loading 4.00 → 3.93. Measured SS drift exactly **0.000e+00**, confirming
`writeoff_enabled` is strictly SS-neutral. The coupon is only ~7.8% of the bond,
so haircutting it alone barely registers.

**E3b (full, `zeta=1`) inverts SPEC Live Claim 1.** `EL_price_D` goes 0.056134 →
0.701743 (**12.5×**, matching the closed form to 1e−12) and the loading collapses
from 4.00/3.17 to **0.37/0.28 — below 1**. The CB becomes *under*-compensated,
receiving ~30% of the actuarially fair expected loss, where the paper's central
claim is over-compensation. **The "monetary-financing objection fails on the
model's own terms" argument does not survive `zeta_writeoff = 1`.**

The mechanism is attributable to the denominator alone: premium income barely
moves (`prem_PV` +9%) while priced expected loss goes ×11.8. It is a repricing of
the expected loss, not a change in what the CB earns.

`psi_lambda_B = 8.5` also stops hitting its 150 bp anchor under E3b (168.9 bp,
+12.4%). Reported, not re-tuned away.

**Correction to the design spec.** It predicted E3b "moves the steady state via
`EL_price`". More precisely: `EL_price` changes value 12.5× but **no SS allocation
moves** (drift 0.000e+00 across eleven quantities) — it multiplies `def_rate`,
which is 0 at SS, so it is allocation-neutral while still changing the linearised
bond FOC and hence every dynamic result.

**Unanticipated: under E3b the named-regime construction breaks.** Peak spread
stops being monotone in γ, so compression targeting has no unique solution. Two
violations on a 40-point grid over γ∈[0,15]: a trivial one at γ≈0.385 and a large
spike at γ≈3.46 (144.4 → 166.6 bp), after which the curve resumes falling to
82.7 bp at γ=15. The isolated spike sits where `I − γ·A_cb` is plausibly
near-singular, so it reads as a linear-algebra pathology rather than economics —
but it means compression-targeted regimes are undefined under full writeoff. E3
therefore evaluates every variant at the **baseline's** γ held fixed, so the model
changes without the policy also changing.

Two bugs were caught by assertions written before the code was first run. The
`gamma_for_compression` monotonicity guard surfaced the finding above rather than
silently bisecting to a meaningless γ. And `expected_EL_price` was initially handed
a stale calibration because `run()` did `from calibration import get_calibration`
at the top, binding the original function *before* the override context opened —
the exact footgun documented in `calibration_override`'s docstring after the
earlier code review. The closed-form check caught it; the fix is to import the
module and resolve at use time.

## 2026-08-03 — orchestrator: `run_all.py` → `docs/experiments_results.md`

`experiments/run_all.py` runs every experiment and renders the generated results
document. E2 runs first because it is self-verifying, so it validates the cache
before anything else reports numbers off it. `--skip-e3` re-renders without paying
for E3's two model re-solves; `--render-only` rebuilds the document from results
already on disk.

Every table carries a provenance stamp read **live** from the calibration —
including a **working-tree-dirty flag**, so a document generated with uncommitted
edits can no longer be mistaken for a clean run at that SHA. The stamp exists
because `run_regimes.py` once shipped a hardcoded "market-value rule" caption while
actually running the par rule.

The document is generated and carries a do-not-hand-edit warning; it also states
plainly when an experiment's results are missing rather than rendering a partial
table silently.

## 2026-08-03 — E1: backstop schedule; cache schema 3 (`delta_b_F`)

`experiments/e1_backstop_schedule.py`. Named regimes canonical, γ **solved** for
0/25/50% peak-spread compression (0 / 5.0798 / 12.7260). Reports the regime table,
A5-1's three German objects separately, the loading schedule, and welfare labelled
secondary. Figure: `experiments/figures/fig_e1_loading_schedule.png`.

**Every cross-check against the independent `code/main.py` pipeline passes.**
Loading 4.00 at `medium` (γ=5.08) vs production's 4.01 at γ=5; 3.17 at
`aggressive` (γ=12.73), correctly below production's 3.44 at γ=10. Peak spreads
hit the compression targets exactly (150.3 → 112.7 → 75.2 bp). `n_inter_D[0]`
(−3.380 / −2.167 / −1.099 % SS) and `Y_D[0]` (−0.0149 / +0.0111 / +0.0338)
reproduce `docs/STATE.md`'s regime table to every printed digit.

**Live Claim 5 confirmed on a fine grid.** The loading schedule is monotone
decreasing at **all 59 finite grid points**, 4.51 (γ=0.51) → 2.07 (γ=30) —
stronger evidence for the self-extinguishing premium than the three points
previously on record.

**Cache schema 3: `delta_b_F_ss` added, and the reason is a bug this caught.**
`cb_pnl` computes each carry leg's SS yield as `delta_b·(1/q_b_ss − 1)`. The first
draft used `delta_b_D` on *both* legs, but `delta_b_F = 0.056779 ≠ delta_b_D =
0.077701` — the two countries' bank books have different measured maturity
ladders. That put the SS spread at −9.2e−04 rather than its true ~1e−17 and would
have silently contaminated `carry_ss_pv`. An assertion written into `cb_pnl`
before the code was ever run caught it. `carry_ss_pv` now comes out 1.2e−16 /
3.2e−16, i.e. numerically zero as it should be. Cache rebuilt (schema 3); E2
re-ran with identical γ and `dY` values.

**Open, needs an author decision: A5-1's third object is misnamed.** The code
reports `Σ β^t (pd_passive − pd_intervention)`, which comes out **negative**
(−0.0015 medium, −0.0047 aggressive) because the backstop lets Greece run a
*larger* primary deficit — it relaxes required austerity. A negative number
therefore means Greece is better off, the opposite of what "Greek fiscal saving"
implies. Either flip the sign or rename it ("austerity relief, PV"). Magnitudes
are unaffected; the label must not ship as-is.

## 2026-08-03 — E2: ΔY decomposition against the `market_clearing_D` identity

`experiments/e2_dy_decomposition.py` + `experiments/test_e2_identity.py`, and five
cache/IRF helpers appended to `experiments/common.py` (`load_cache`,
`cache_outputs`, `irf_from_cache`, `named_regime_gammas`, `regime_irfs`). The
named regimes are canonical: γ is **solved** for 0/25/50% peak-spread compression
(0 / 5.0798 / 12.7260), not chosen as round numbers, so the regimes keep their
meaning across recalibrations.

**Self-verifying by construction.** The decomposition is the linearised
`market_clearing_D` identity `dY = P_ss·dC + C_ss·dP + dI + dG + dΦ + dT + dNX`,
and `goods_mkt_D` is a *targeted* solver residual, so the components must sum to
`dY` to solver tolerance. Achieved **5.8e−17 / 1.1e−16 / 1.5e−16** against a 1e−7
assertion that halts rather than warns.

**Finding — the headline output number is the residue of two much larger
offsetting channels.** Passive → aggressive, `dY[0]` moves +4.87e−04 while
investment moves +2.16e−03 and net exports −1.85e−03, each ~4× the headline and
opposite in sign. This confirms `docs/SPEC.md`'s standing caution ("a small
headline output number can be *only* small because two large channels are netting
out — and they land on different households") **as a measured property of this
calibration**, and settles the gate it placed on the trade-channel claim.

**The `Y_D[0] > 0` watch item is answered.** `dY[0]` = −0.0149 / +0.0111 /
+0.0338 % of SS. The proximate driver of the sign flip is **consumption quantity,
not investment**: at `medium`, investment is still negative on impact and output
is positive only because consumption outweighs it. Consumption is already positive
on impact at `passive` with no backstop, so it is not manufactured by the policy
rule — which argues against pure linear-rule overshoot. But the magnitudes
(0.01–0.03% of SS) are small differences of much larger terms and should not be
leaned on. Full table in `docs/STATE.md`.

`Phi_D` and `G_D` are **verified** zero rather than merely uncached: `Phi_D` has
no Jacobian column (the portfolio adjustment cost is quadratic about its anchor,
so its level deviation is second-order), and `G_D` is absent from `G_tpi.outputs`
because government spending is constant. The identity closes *because* both are
genuinely zero. Cache rebuilt under schema 2 (7m27s); `G_tpi[cb=0]` vs baseline
`max|err| = 0.00e+00`.

Also noted: `EL_price_D` is **0.056134** at the live calibration, not the `0.0717`
still quoted in older doc sections and CLAUDE.md — that figure predates the EBA
`delta_b`/`q_b`. It is the TPI loading's denominator, so it must be re-derived
wherever quoted, not copied. Not fixed in this commit.

## 2026-08-03 — `experiments/` package: cache schema v2 (call-time fingerprint)

First commit of the new `experiments/` package (branch `experiments`), which will
produce the paper's standard policy results on top of the regimes cache layer.
`code/main.py` untouched.

**The fix that had to come first.** `regime_model.cache_path` built its filename
from `CAL_FINGERPRINT`, computed at module **import**. Experiment E3 needs to solve
the model under a calibration override applied at run time; with an import-time
fingerprint its cache would have been written to the baseline filename and silently
overwritten it. Now computed at call time, stamped into the `.npz` as
`cal_fingerprint`, and asserted against the live calibration on load, with a clear
`FileNotFoundError` naming the rebuild command when no cache matches.

`CACHE_SCHEMA = 2` now appears in the filename: the calibration fingerprint alone
cannot detect a change to the cached *output list*, so without it an old cache would
reload under an unchanged name missing the new keys — invisibly, since `irf_all`
discovers outputs by scanning cache keys. Added `Phi_D` and `def_rate_D` to
`REQUIRED` (`Phi_D` closes the `market_clearing_D` identity for E2; `def_rate_D` is
the off-path expected-loss leg for E1). `T_D` went to `OPTIONAL` deliberately —
`T0=T1=0` makes it identically zero, so zero-filling is correct rather than a silent
hole, and E2's closure assertion catches it if that ever changes.

Also `build_caches` now reads `psi_lambda_B` live rather than from the import-time
`PSILAM_MAIN` constant, so an override wins there too.

New: `experiments/common.py` (calibration-override context manager, unit helpers,
provenance stamp, results writer) and `experiments/test_common.py` (regression
guards including one that the override changes the cache filename and does not
leak out of the context manager).

**Hardened after code review (same day).** `calibration_override` now raises
`KeyError` on an override key not present in the calibration dict — previously a
typo (e.g. `psi_lambda_b_D` for `psi_lambda_B_D`) would silently add a junk key
while leaving the real parameter at its default, producing a wrong-but-plausible
number without any error: exactly the failure mode this whole package exists to
close off, one level up from the cache-fingerprint fix above. `write_results` now
serialises numpy arrays/scalars properly (`json.dump(default=float)` raised on any
multi-element array — the normal shape of an IRF payload) and passes
`allow_nan=False`, so a `NaN` in a result is a loud `ValueError` at write time
rather than a token that travels silently into a table and that strict JSON
parsers reject anyway. Also fixed: the exception-restore test asserted on a
name bound at import (`from calibration import get_calibration`), which cannot
observe the module attribute the context manager patches and so passed
regardless of whether the restore worked — rewritten to assert on
`calibration.get_calibration` directly, and verified to fail when the `finally:`
restore is removed. `cache_path`'s `fingerprint` parameter is now actually used by
`load_cache` instead of being computed and discarded twice. 17/17 tests pass
(`experiments/` + `diagnostics/regimes/test_lottery_math.py`).

## 2026-07-31 — regimes re-run at the broad scope; `PSILAM_BREAKDOWN` re-derived; units bug fixed

Follow-up to the broad-scope commit below: the policy-regime diagnostics had never
been run at the live calibration, and their hard guard blocked it.

**`PSILAM_BREAKDOWN` re-derived: 2.5 → 15.0.** New
`diagnostics/psilam_breakdown_sweep.py` sweeps 16 points above the live 8.5,
solving the SS once (`psi_spread` is exactly linear in `psi_lambda_B`) and
re-solving the Jacobian per point with both dials moved together. Peak spread
150.3 → 223.8 → 273.6 → 625.3 → 1034.5 → 8903.8 bp at 8.5 / 14 / 20 / 25 / 26 / 27,
then **sign-flips at 28** — a pole between 27 and 28. The A7 >1000bp flag first
fires at 26; the *first* pathology is earlier, `n_inter_D[0]` shrinking over
`[14,18]` while the spread still rises. The guard is set from that, not from the
pole, so it keeps real margin. The old 2.5 was CT1-scope and would have blocked
the live calibration outright.

**Regimes re-run, all three figures regenerated.** Caches rebuilt (the calibration
fingerprint forced it). `A_cb=-1.889e-2` — backstop still compresses, SA-1 absent.
`gamma_aggressive=12.726` / `gamma_medium=5.080`; peak spread 75.2 / 112.7 /
**150.3** bp, so `run_regimes.py`'s 120–180bp band now passes. A6 amplifier
invariance holds in both Stage A (3.78/4.67/5.44 at `psi_lambda_B=0`) and the
lottery (2.51/3.09/3.59). 18/18 tests pass.

**New: `diagnostics/regimes/certainty_equivalence.py`** — answers whether the
regime lottery is degenerate under first-order certainty equivalence. It is not,
but the non-degeneracy is ~9% of what the Stage B figure shows: comparing the
lottery to one *known* CB at `gamma_bar` (same silence-until-`k` timing) gives
`LOT − CE = +1.477bp` at ergodic beliefs against a 16.6bp total belief-shift
effect. The wedge is exactly `A_cb[0,:] @ (cb^e − cb_CE)` (verified to 6e−15 bp);
mechanism is that uncertainty **back-loads** the expected backstop while the
date-0 spread weights near-term purchases most. See `docs/STATE.md`.

**Units bug fixed (pre-existing, affected published numbers).** SSJ IRFs are
*level* deviations, so `×100` is a percent only where the SS level is ≈1.
`Y_D_ss≈1` passes; `n_inter_D_ss=2.138` and `K_D_ss=10.8` do not. `main.py`
printed `n_inter_D[0]×100` as `%` — the widely-quoted **−7.227%** is the level
deviation and the true impact is **−3.380% of SS** — and the Stage A figure
titled two panels `(%)` on the same basis (2.1× on net worth, **10×** on
capital). Both fixed; level deviations retained alongside for continuity.
Consequence: `PT-1`'s pass-through moment is **−2.25%/100bp**, not ≈−4.5%; still
inside the Acharya–Drechsler–Schnabl-implied range but at its low end.

**Watch item recorded:** `Y_D[0]` is positive under both intervening regimes and
the A5 `dY_D` trough never goes negative — output never falls at all under the
backstop. At `gamma_aggressive=12.7` that is more likely linear-rule overshoot
than economics; flagged in STATE.md, not reported as a result.

Full `code/main.py` re-verified (exit 0): residuals, `b_gov_D[499]=+1.4e−05`,
`ρ_b=0.845`, peak spread 0.376pp, TPI loading 4.35/4.01/3.44 all unchanged.

---

## 2026-07-31 — broad capital-funding-sector scope: the EBA calibration goes LIVE

The last blocker was the **scope of `n_inter`**, not any parameter. EBA CT1 is the
capital of the *sovereign-exposed stress-test sample*; the model's `n_inter` is
the net worth of the agent intermediating the **whole capital stock**. Using CT1
pins `omega_K` tiny and makes the accelerator gain ~`1/n_inter`.

**New `BANK_SCOPE` in `code/calibration.py`** (`"broad"` live, `"ct1"` kept for
comparison). Under `"broad"` the intermediary is the entire capital-funding
sector and its net worth follows from the measured leverage and the balance
sheet, `N = (Q*K + sovereign)/theta`, so **`omega_K = 1` by construction** and the
passive-fund device disappears entirely.

| | CT1 scope | **broad (live)** | pre-EBA placeholder |
|---|---|---|---|
| `n_inter_D/F` | 0.408 / 0.175 | **2.138 / 1.627** | 3.0 / 3.0 |
| `phi_own_D/F` | 2.390 / 2.758 | **0.456 / 0.296** | 0.25 / 0.25 |
| `omega_K_D/F` | 0.117 / 0.067 | **1.0 / 1.0** | 1.0 |
| `theta_D/F` | 5.51 / 6.94 (measured) | **5.51 / 6.94** | 4.0 |

Kept measured: `theta`, the sovereign book, `delta_b` (ladder), `K/Y`. Given up:
`n_inter` as observed CT1, and `phi_own = 2.39` as a *model* parameter — 2.39 is
concentration *within the stress-tested slice*, not within the whole
capital-funding sector, and only the latter is what the model's `phi_own` means.
**Load-bearing assumption:** applying the EBA sample's `theta` to the whole
sector. This is now the only such assumption left in the bank block.

Consequences: `Delta` returns to **0.2/0.4** (the 0.85 bound was a CT1 artifact;
at `phi_own=0.456`, `f*theta = 0.661 > 0.367`), the fiscal rule to
`phi_lamb=0.15` / `mv_rule=0`, and `psi_lambda_B` retunes to **8.5**.

**Verified end-to-end** (`code/main.py`, exit 0):

| Check | Result |
|---|---|
| over-identifying `K` | `K_D = 10.800`, `K_F = 10.832` (target 10.8) |
| IC residual | −8.9e−16 (D) / 0.0 (F) |
| Walras | `ca_res_D = 6.9e−17`; all block residuals < 1e−8 |
| stability | `b_gov_D[499] = +1.4e−05`, `ρ_b = 0.845` |
| `n_inter_D[0]` | **−7.227%** ✓ |
| `Y_D[0]` | **−0.0149%** ✓ — **Y-1 RESOLVED** |
| `rk_D`, `rk_F` | both exactly 0.010000 — **RK-1 RESOLVED** |
| peak spread (γ=0) | 0.376pp = **150.4bp**, on target |
| TPI loading | **4.35 / 4.01 / 3.44** at γ=2/5/10, declining |

The paper's self-extinguishing-premium claim survives on a properly identified
calibration.

`psi_lambda_B = 8.5` is higher than the historical 0.31 / 1.18 / 3.0 because the
broad scope's `phi_own = 0.456` is well below CT1's 2.39, so more of the default
loading comes from the friction. The mapping is smooth and monotone with **no
breakdown region** — 5.4 / 23.1 / 58.6 / 111.0 / 142.5 / 157.8 bp at
0 / 1 / 3 / 6 / 8 / 9 — so the old "breakdown above ~1.5-2.0" warning (specific to
CT1-thin net worth) does not apply.

**Method bug found and fixed.** `psi_spread_D` is derived from `psi_lambda_B`
inside `_apply_ss_anchors`, so a sweep **must re-solve the SS per point**.
Patching the flag onto an already-solved SS leaves `psi_spread` stale and
*inverts the apparent sign* of the spread response — an earlier sweep did exactly
that and reported the spread falling in `psi_lambda_B`. Those numbers are void;
all figures above come from full per-point re-solves.

---

## 2026-07-31 — fix `omega_K`: capital fund holds a fixed QUANTITY, not a fixed share  [301ffd2]

`omega_K` as a **fixed share** was the defect. The passive fund held
`(1-omega_K)*K` at all times, so it mechanically *mirrored* bank deleveraging —
a 1% fall in the bank's book dragged the other ~88% of the capital stock with
it. That is the `1/omega_K` lever, and it is an assumption nobody would defend
stated plainly: non-bank capital holders do not shrink in lockstep with bank
equity.

**Fix: `fund_rule_D/F`** (`code/calibration.py`, committed at **1**).

| `fund_rule` | fund holds | bank holds | `dK/dN` |
|---|---|---|---|
| 0 (legacy) | `(1-omega_K)*K` | `omega_K*K` | `theta/omega_K` |
| **1** | **constant `K_fund`** | `K - K_fund` | **`theta`** |

With `K_fund = (1-omega_K)*K_ss` the rules are **identical in steady state** —
the change is purely dynamic. Verified: `lambda_gk_D=+0.9271`, `Omega_D=+4.62`,
`K_D=10.800` under both. Gain 47.1 → 5.5. At `omega_K=1` (pre-EBA), `K_fund=0`
and the rules coincide exactly, so the placeholder calibration is untouched by
construction — confirmed by a full `main.py` run reproducing
`n_inter_D[0]=-3.0009%`, `Y_D[0]=-0.0261%`, `ρ_b=0.8451`, `b_gov_D[499]=-1.3e-5`,
IC residual `+8.9e-16`.

Touched four sites per country (`smart_steady`, `intermediation_IC` kappa,
`capital_fund`, `k_balance_sheet`) plus the `steady_state.py` over-identifying
print, which is now exact for both countries (`K_D = K_F = 10.800`).

**Helps substantially, does not stabilise alone:** `b_gov[499]` +3.98e2 → +2.25e2
at `psi_lambda_B=0`, and −2.08e4 → −1.23e3 at `psi_lambda_B=1` (17×).

**Three compounding amplifiers** (at `fund_rule=1`, `psi_lambda_B=0`):

| bank block | `phi_own` | `b_gov[499]` |
|---|---|---|
| pre-EBA (`n=3.0, theta=4.0, omega_K=1`) | 0.25 | **−4.4e−08 stable** |
| EBA (`n=0.41, theta=5.51`) | 0.25 | +1.85 |
| EBA | 2.39 (measured) | +5.32e+02 |
| pre-EBA | 2.39 | solver failed to converge |

1. the fixed-share fund (`1/omega_K` ≈ 8.5×) — **now fixed**;
2. thin measured net worth (`n_inter` 3.0 → 0.41): ~1e−8 → 1.85;
3. measured concentration (`phi_own` 0.25 → 2.39): 1.85 → 532.

(2) and (3) are both *measured*, so neither can be tuned away. `EBA_CALIBRATION`
stays **False**.

**Also ruled out this session** (all with evidence in `docs/eba_calibration.md`):
the sovereign-risk schedule — `def_scale` 0.25 → **0.00** makes it *worse*
(−2.1e4 → −1.6e5), so flattening it (e.g. a bounded `tanh`) cannot help; note
also that at first order only the local slope `a*b` matters, so a bounded form is
observationally equivalent in the linearised IRFs. Separately, measuring
`B_supply` (2.4 → 1.116) steepened that schedule's local slope ~45% (0.0799 →
0.1158) as a side effect, since lower debt sits further up the concave curve —
real, but not the driver.

---

## 2026-07-31 — fix the collateral mapping; EBA steady state now correct, dynamics still explosive  [c6230a2]

Follow-up to the rebuild below. The previous entry concluded the feasible `Delta`
set was empty. **That was wrong, and the reason is instructive:** the emptiness
came entirely from an undocumented convention, not from the economics.

`ic_delta_calibration._ic_delta` back-solved `Delta` from one equation in two
unknowns and closed the system with a hardcoded
`ratio = Delta_cross/Delta_own = 2.0`. With `Delta_cross <= 1` that caps
`Delta_own <= 0.5`, against the `> ~0.73` GK well-posedness requires. It was a
back-solve closure masquerading as a consistency check — which is exactly why the
inherited `0.2/0.4` pair "passed" it.

**Fixed.** The convention is gone. `Delta_own`/`Delta_cross` are free structural
parameters and the module now checks the IC **residual** directly
(`ic_residual`, tol 1e-8, plus a positive-divertable-leverage check). Verified on
the pre-EBA calibration: residual `+8.9e-16` (D) / `+1.8e-15` (F).

**`Delta = 0.85/0.90` resolves the steady state:**

| `Delta_own`/`Delta_cross` | `lambda_gk_D` | `lambda_gk_F` | `Omega_D` | `Omega_F` | `K_D` | `K_F` |
|---|---|---|---|---|---|---|
| **0.85 / 0.90** | **+0.927** | **+0.960** | +4.62 | +5.98 | 10.80 | 10.65 |
| 0.90 / 0.95 | +0.488 | +0.456 | +2.49 | +2.91 | 10.80 | 10.65 |

`lambda_gk_D = +0.927` is essentially the pre-EBA `+0.923`, so the amplification
block keeps its strength with measured concentration. Not a fudge: leverage of
5.5× on a 43%-sovereign book is inconsistent with sovereigns being good
collateral — if they were, the bank would lever further and `theta=5.5` would not
bind. Also the right story for 2010–12 Greece (collapsing GGB eligibility, rising
ECB haircuts). Note `Delta=0.80/0.90` lands just past the `lambda_gk` **pole**
(`lambda_gk_F=-12.45`): the closed-form bound ignores endogenous `rn`, so the real
frontier is a fixed point — sweep, don't solve.

**STILL BLOCKING — dynamic instability (GK-2).** With the SS correct,
`b_gov_D[499] ~ 1e2–1e3` vs a ~1e-5 target. Measured amplification is
`theta*phi_own = 13.17` against `4.0*0.25 = 1.0` for the placeholder. Diagnosed:

| Test | Result |
|---|---|
| `psi_lambda_B = 0` | still explosive, `b_gov[499] = -2038` |
| `phi_lamb` 0.6 → 25 | **flat** (peak ~1.1e7bp at 0.6, 1.5 *and* 25) — not the fiscal mode |
| `mv_rule = 1` | does not fix it |
| `chi1` 0 → 0.5 | peak spread **1.1e7bp → 6.0bp**, `b_gov[499]` −2038 → +70 |
| `chi1` ∈ [0.2, 5.0] | `b_gov[499]` stays 70–560 — amplitude damped, root not removed |

`chi1` (Auclert intermediary adjustment cost, currently 0) is the strongest lever
and restores a sane spread response, but no tested value removes the root. Three
open routes in `docs/eba_calibration.md` → *Dynamic instability*.

**Also fixed:** `steady_state.py` did `from calibration import EBA_CALIBRATION`,
binding the flag at import time so a sweep flipping it would silently keep the old
portfolio targets — the same stale-binding trap as the regimes cache key and
`PSILAM_MAIN`. Now imports the module and resolves at call time. The GK guard also
reports both countries' `lambda_gk`/`Omega` on failure, not just the failing one.

`EBA_CALIBRATION` stays **False**. Full `main.py` re-verified: `n_inter_D[0]=-3.0009%`,
`Y_D[0]=-0.0261%`, `ρ_b=0.8451`, `b_gov_D[499]=-1.3e-5`, all residuals < 1e-8.

---

## 2026-07-31 — EBA calibration rebuilt and identified; measured moments found infeasible in the GK block  [050e54a]

Rebuilt the EBA 2011 moment set to be **identified rather than back-solved**, then
found that the measured moments are structurally incompatible with the model's
collateral constraint. The rebuild is sound; the incompatibility is the finding.

**Live calibration unchanged.** `EBA_CALIBRATION = False` (`code/calibration.py`)
selects the pre-EBA values bit-exactly — verified parameter-by-parameter and by a
full `main.py` run reproducing `n_inter_D[0]=-3.0009%`, `Y_D[0]=-0.0261%`,
`ρ_b=0.8451`, `b_gov_D[499]=-1.3e-5`, TPI peaks 0.468/0.330/0.236/0.163 pp.

**Newly measured** (`code/eba_calibration.py` → `data/eba_moments.json`):

| Parameter | Was | Now | From |
|---|---|---|---|
| `delta_b_D/F` | 0.10, no EBA counterpart | 0.0777 / 0.0568 | maturity ladder (`MATURITY_CODE` 125–155, previously unread) repriced at the end-2010 market yield |
| `theta_D/F` | 4.0 assumed | 5.51 / 6.94 | (corp ex-CRE + CRE + sovereign) EAD / CT1 |
| `omega_K_D/F` | back-solved plug | 0.117 / 0.067 | corp+CRE EAD ÷ K at `K/Y_ann=2.7` |

- **`delta_b` retires the F-1 duration blocker.** The old "port 0.036/0.038 (7y)"
  target measured the sovereign's *whole outstanding stock*; `delta_b` governs the
  *bank-held book at the yields banks faced*, whose modified duration is 3.12y
  (GGB) / 4.22y (Bund). Close to the old 0.10 — `mv_rule=1` + `phi_lamb=0.60` not
  needed.
- **`omega_K` kept, not dropped.** Banks fund ~12%/7% of the capital stock, so
  `omega_K=1` is counterfactual. Now measured, with `K` an output: the balance
  sheet delivers `K_D = 10.800` vs the 10.8 target.
- **Amplification moment is Acharya–Steffen, not CT1 depletion.** Ladder ×
  observed 2011 yield moves (GGB 12.01%→21.14%) gives mechanical
  **−5.73%/100bp** of CT1, −39.8% realised over 2011. The pre-EBA calibration
  generates −0.61%/100bp, **~10× too weak** — almost entirely `phi_own`
  (0.25 vs 2.39), not duration. This is what `psi_lambda_B=3.0` was standing in
  for, and it reframes S-1's "89% collateral friction" as a calibration artifact.
  The 2011 adverse scenario is *deliberately rejected* (it excluded banking-book
  sovereign default) and guarded by a test.

**BLOCKING FINDING — GK feasibility.** The block is well-posed only if
`f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross`. At the measured
moments this is violated by −1.26 (D) / −1.42 (F):

```
lambda_gk_D = -0.0869   Omega_D = -0.3013
lambda_gk_F = -0.0723   Omega_F = -0.3217
```

Negative IC multiplier and negative banker franchise value — yet the solver
converged, `ca_res_D = -7.7e-17`, the IC-δ check passed exactly, stability passed,
and the TPI loading still declined. Every IRF meaningless. C-1's failure mode
reborn.

The bound is `Delta_own > ~0.73`, but `_ic_delta` hardcodes
`ratio = Delta_cross/Delta_own = 2.0`, which with `Delta_cross <= 1` caps
`Delta_own <= 0.5`. **Feasible set empty.** `f` would need > 0.349 (literature
0.03–0.12), `theta` > 16.03 (measured 5.51). Tested: `Delta = 0.80/0.90` →
`lambda_gk_F = -12.45`, `Omega_F = -75.91`, just past the pole. Escaping is a
modelling decision; clearing the pole needs `Delta_own ≈ 0.85–0.95`, i.e.
sovereign bonds nearly worthless as collateral — removing the channel the doom
loop runs on.

**New guard: `steady_state.assert_gk_well_posed`**, called from
`_apply_ss_anchors`, so it fires on every solved SS in both `steady_state.py` and
`depreciation_calibration.py`. Exports `gk_feasibility_margin` / `min_Delta_own`.
The most valuable artifact here — it makes this failure class impossible to commit
silently.

**Also fixed** (all pre-existing, surfaced while working):
- `diagnostics/regimes/` cache filenames keyed only on `psi_lambda_B`, so the
  `psilam=0` cache would be **silently reused across different models**. Now
  carries a SHA of the whole calibration. Stale pre-EBA caches deleted.
- `PSILAM_BREAKDOWN` 4.0 → 2.5, with its net-worth dependence documented.
- `run_regimes.py`'s figure suptitle hardcoded `"market-value rule"` while the
  committed calibration is the par rule; now reads live `mv_rule`. (The earlier
  fix covered the JSON provenance string only.)
- `data/README.md`: full worksheet-4/5 code map and the maturity-bucket decode.

**Tests:** `code/test_eba_calibration.py` 10/10 pass, including ladder-exhausts-total,
duration-below-maturity, `delta_b` round-trip, the balance-sheet identity, the
mechanical-MTM magnitude, and a guard that the adverse scenario is never used.

**NOT done:** the `psi_lambda_B` retune to 150bp. A sweep was run but on the
degenerate model, so it is **void and discarded**; redo it after the GK
feasibility question is settled.

---

## 2026-07-31 - fix A6 lottery invariance test (measurement-window error)

The A6 amplifier-invariance check in `uncertain_regime.py` ranked the **full-sample**
peak spread across the three lottery branches. That peak is the common pre-`k` spread —
no branch has acted before revelation (the same pre-`k` identity asserted in sec 10.3) —
so it is identical across branches *by construction*, and `pk0[0] < pk0[1] < pk0[2]` was
comparing floating-point noise. Measured both ways:

| `psi_lambda_B` | full-sample peaks | strict `<` | gap |
|---|---|---|---|
| 3.0 (calibrated) | all `184.992662` bp | **False** | `-5.1e-13` bp |
| 0 | all `9.302930` bp | **True** | `+2.3e-14` bp |

Same expression, opposite verdicts, decided by the last ULP — the reported "YES" was
luck, not a result.

**Fix:** rank on the **post-revelation window `t>=k`**, where branches have actually
diverged; check **both** amplifier settings (that is what makes it an invariance test
rather than a single reading); and require separation above a `1e-3` bp margin so
numerical noise can never produce a verdict.

| `psi_lambda_B` | aggressive | medium | passive | ordered |
|---|---|---|---|---|
| as calibrated | 76.98 | 117.87 | 160.19 | YES |
| 0 | 4.91 | 5.89 | 6.43 | YES |

**A6 invariance genuinely holds in the lottery** (separation 0.548bp at `psi_lambda_B=0`,
~550x the margin) — the earlier "false pass" was a measurement-window error, not a
failure of the economics. The deterministic Stage A A6 (9.33/8.54/7.14bp) is unaffected
and independently confirms it.

No other Stage B output changed: `E_pi[W_D]=-14.3343`, impact spread 151.46 -> 187.19bp,
the `k` sweep and the welfare table are all identical. Full suite re-verified, exit 0:
`test_lottery_math.py`, Stage A, Stage B-lite.

Not touched: the Stage B table's `lottery` column is 185.0bp for all three branches, and
that is *correct and deliberate* — before revelation nobody can tell the branches apart,
which is exactly what the uncertainty premium (lottery minus known-delayed) prices.

---

## 2026-07-31 - doc-sync enforced by a git hook; regimes feature unblocked and verified  [this commit]

**Doc-sync is now actually enforced.** Correction to yesterday's entry: the policy *was*
already implemented as a Claude Code PreToolUse hook
(`.claude/hooks/require-docs-before-commit.sh`, tracked); the "not installed" claim came
from checking `.git/hooks/`. Two real gaps found:

- A PreToolUse hook cannot see terminal commits at all. Added **`.githooks/pre-commit`**,
  the git-native twin and now the primary enforcement: enable once per clone with
  `git config core.hooksPath .githooks`. It inspects the git *index*, so unlike the
  command-string-matching PreToolUse hook it cannot false-positive. Verified: exits 1
  with the block message when code is staged without the docs.
- The PreToolUse hook's `if` filter is prefix-matched, so it never fired on compound
  commands (`cd X; git commit ...`). Removing the filter was tried and **reverted** - it
  makes the script scan every Bash command string and it then denied an unrelated
  documentation edit. The filter stays; the git hook is the real gate, and the
  PreToolUse hook is early feedback for the simple-command case.
- Stale `PROCESS.md` reference in the script header corrected to `PROGRESS.md`.

**`run_regimes.py` provenance string no longer hardcoded.** It read `"mv_rule=1"`
regardless of calibration, silently mislabelling every number in
`regimes_calibration.json`. Now reads `psi_lambda_B`, `mv_rule`, `recovery_rate` live.

**Policy-regime feature now runs end-to-end at the pre-EBA calibration** - it did not
before. `gamma_for_compression`'s scan range narrowed 60 -> 25: bisection's validity
condition is monotonicity on the *bracketing interval*, and `peak(gamma)` falls
monotonically 187.2 -> 34.2bp on [0,25] then ticks up 1.1bp at gamma=30 (saturation at
81% compression, not economics), which aborted the whole run.

Verified, all exit 0: `test_lottery_math.py`, Stage A `run_regimes.py`, Stage B-lite
`uncertain_regime.py`.

- Stage A: `A_cb=-2.406e-2` (compresses; SA-1 absent), passive peak 187.2bp,
  `gamma_aggressive=5.0813` (50% compression), `gamma_medium=1.5730` (25%). Peak spread
  187.2 / 140.4 / 93.6 bp passive / medium / aggressive.
- Stage A A6 (deterministic, `psi_lambda_B=0`): 9.33 / 8.54 / 7.14 bp - amplifier
  invariance genuinely holds, clean separation.
- Stage B-lite: sec 10.3 assertions PASS; uncertainty premium +69.0 / +34.4 / -2.2 bp
  (aggressive / medium / passive); impact spread rises with the passive belief weight
  (151.5 -> 187.2bp), sign computed not targeted.

**Open - A6 invariance at the LOTTERY stage is a false pass.** All three branch peaks are
numerically identical (`9.302980` bp, gap `0.000000`), yet
`pk0[0] < pk0[1] < pk0[2]` returns True on last-ULP floating-point noise. Structural, not
a coding slip: with `k=2` the peak spread falls in the common pre-`k` window where no
branch has acted, so the branches coincide there by construction. Stage A's deterministic
A6 is the informative test; the lottery A6 line should be reworded or dropped rather than
reported as a result.

**Also open, none EBA-related:** `psi_lambda_B=3.0` gives 187.2bp, outside
`run_regimes.py`'s own 120-180bp sanity band (it logs "investigate"); `delta_b=0.10`
still needs `mv_rule=1` AND `phi_lamb=0.60` jointly; `beliefs.json` dates from
2026-07-23, before the calibration revert.

---

## 2026-07-30 — revert calibration to pre-EBA; drop audit_artifacts; F-1 hard-break measured  [this commit]

**Calibration reverted to the pre-EBA values** in force at `abcbb6e` (the last commit
before the EBA work began at `eade414`). Verified parameter-by-parameter against that
commit: the only additions are `omega_K_D/F=1.0`, which is the structural no-op
(capital fund empty, `div_fund=0`, so the pre-EBA balance sheet is recovered exactly).

| Parameter | EBA | now |
|---|---|---|
| `psi_lambda_B_D/F` | 1.1793 | **3.0** |
| `n_inter_D/F` | 0.408 / 0.175 | **3.0 / 3.0** |
| `omega_K_D/F` | 0.0602 / 0.0190 | **1.0 / 1.0** |
| `B_supply`/`b_gov`/`b_gov_ss` | 1.19 / 0.591 | **2.4** |
| `phi_lamb_D/F` | 0.60 | **0.15** |
| `mv_rule_D/F` | 1.0 | **0.0** (par) |
| `phi_bF_D_ss` / `phi_bD_F_ss` | 0.018 / 0.069 | **0.25 / 0.25** |
| portfolio targets (`steady_state.py`) | 2.39/0.018/0.069/2.76 | **0.25/0.15/0.15/0.25** |

**All structural fixes retained** — C-1 (multi-asset `lambda_gk`), W-1/W-2/W-3, T-2, A-2,
TPI-1, the `omega_K` generalisation and the capital-key conduit are untouched. Only
parameter values moved. `recovery_rate_D/F=0.30` (EL-1) was **kept**, not reverted: 0.00
asserts a counterfactual 100% loss-given-default.

**Verified end-to-end** (`code/main.py`, exit 0): `n_inter_D[0]=-3.0009%`,
`Y_D[0]=-0.0261%` (both negative ✓); `b_gov_D[499]=-1.3e-5` (default) / `+7.8e-5` (TFP);
`rho_b=0.8451`; IC-δ exact at 0.2000/0.4000; `max|ca_res_D|=6.3e-8`,
`max|goods_mkt_F|=1.1e-9`. TPI peak spread monotone in γ: 0.468→0.330→0.236→0.163 pp.

**Finding F-1 sharpened — the `[0.15,0.18]` zone is a hard break, not a mild one.**
Measured directly by switching `mv_rule` on the otherwise-pre-EBA calibration:

| `phi_lamb` | `mv_rule` | `n_inter_D[0]` | `Y_D[0]` | `b_gov[499]` | spread |
|---|---|---|---|---|---|
| 0.15 | 1 | **-1554.0%** | **+0.170%** ✗ | 1.6e-2 | 124.0bp |
| 0.60 | 1 | -5.89% | -0.024% ✓ | 0.0 | 219.8bp |
| 0.15 | 0 | -3.00% | -0.026% ✓ | -1.3e-5 | 187.2bp |

`mv_rule=1` and `phi_lamb=0.15` are **not a usable pair**. Porting empirical duration
(`delta_b=0.036/0.038`) is therefore a two-parameter move, not one.

**Default-loading decomposition measured.** `EL_price_D=0.0717` vs `psi_spread_D=0.8385`
(= `lambda_gk*psi_lambda_B/(beta_inter*Omega)` = 1.8031*3.0/(0.9975*6.4670)): the
fundamental expected loss is only **10.9%** of the total default loading; the GK
collateral friction is **89%**. Quantifies S-1 — the model is an amplification story,
not a credit-loss story. Consequence: `recovery_rate` 0.00→0.30 moves the loading by
3.3% partial / 6.0% in GE (peak spread 199.2→187.2bp ann).

**Units correction:** `spread_rb` is a *quarterly* rate deviation; the 150bp target is
annualised. Annualise with ×4×1e4 (`BP_ANN` in `run_regimes.py`). At `psi_lambda_B=3.0`
the model gives **187.2bp ann**, i.e. ~25% *over* target — not under.

**`audit_artifacts/` deleted (30 files).** `run_audit.py` carried its own hardcoded copy
of the calibration instead of importing `get_calibration()`, so it silently tested a
different model than `code/main.py` — its `ca_res_D` "regression failure" (1.479e-7 vs
1e-7) is reproduced exactly by the committed baseline artifact and predates any change
here. The `ms-regime` branch had already deleted it for the same reason (`2fa1b55`).
`code/main.py` is now the only regression path. Doc references updated; historical
provenance citations in finding write-ups left verbatim under a deprecation note.

**`diagnostics/regimes/` — two fixes and one open blocker.** `PSILAM_MAIN` was hardcoded
to `1.1793` while serving as **both** the guard and the **cache filename key**, so it
pointed at a stale cache built under a different model; it now reads the live
calibration. The `psi_lambda_B<1.5` guard was EBA-specific (thin net worth) and is now
`PSILAM_BREAKDOWN=4.0` for pre-EBA net worth. Open: `gamma_for_compression` scans
`linspace(0,60,61)` and requires *global* monotonicity; at the current calibration
`peak(γ)` is monotone on `[0,25]` and ticks up 1.1bp at γ=30 (81% compression,
saturation), so the run still aborts. Targets sit at γ≈1.6 (25%) and γ≈5.1 (50%) — far
inside the monotone region; narrowing `hi` to 25 would fix it (bisection only needs
monotonicity on its bracketing interval). **`run_regimes.py:75` hardcodes `"mv_rule=1"`
in its JSON provenance string — now false.** Note the regimes feature has *no*
structural `mv_rule` dependence; `A_cb=-2.406e-2` (compression) holds at `mv_rule=0`, so
SA-1 does not recur.

**`docs/PROCESS.md` deleted** — superseded by PROGRESS.md at `0c99013`; dangling
references in CLAUDE.md/HANDOFF.md/PROGRESS.md fixed. The doc-sync pre-commit hook is
documented but **not installed** in `.git/hooks`.

---

## 2026-07-24 — phi_lamb doc correction, Finding F-2 (ring), doc-sync hook  [803ecd2, +this commit]
- **Corrected a stale `phi_lamb`**: the committed value is **0.60** (~Bohn), not the
  `0.15` STATE.md's calibration table/IRF header claimed nor the "0.30" in
  `calibration.py`'s comment. STATE.md line 37 and the code itself already had 0.60
  (the run's `ρ_b=0.373` inverts to 0.60). `phi_lamb` governs the debt/fiscal mode
  only — well above the F-1 near-unit-root zone.
- **Finding F-2 (financial-accelerator ring).** The damped ~25q oscillation in the
  asset-price / bank-net-worth IRFs (`n_inter_D`, `q_b_D`, `q_b_F`, `C_D`; absent
  from `Y_D`/`w_D`) is the GK IC/leverage accelerator: **|λ|=0.954, period ≈25q
  (~6.25 yr), half-life ≈14.6q (~3.6 yr), R²=1.0.** A 100× PAC sweep
  (`psi_bF_D=psi_bD_F`, 0.05→5.0) moves |λ| by <0.006 and leaves the period fixed;
  committed `PAC=0.5` sits at the modulus *minimum*. The ring is intrinsic to the
  amplification block (pinned to the 150bp spread target) → **no free cosmetic
  fix**; a structural financial cycle to describe, not patch. New:
  `audit_artifacts/pac_sweep.py`, `pac_sweep_results.json`.
- Re-ran `code/main.py` end-to-end (clean: SS `goods_mkt_D≈-4.8e-7`,
  `ca_res_D≈2.3e-16`; `b_gov_D[499]≈2e-6`). Confirmed default-shock
  `n_inter_D[0]=-2.83%`, peak spread `+0.392pp`, TPI loading 3.59/3.03/2.47 and
  spread compression 0.392→0.244pp (−38%) over γ=0→10.
- **Added the pre-commit doc-sync hook** (`.claude/settings.json`,
  `.claude/hooks/require-docs-before-commit.sh`) and **this changelog**
  (`PROGRESS.md`); the hook's required-doc set is `STATE.md / PROGRESS.md /
  HANDOFF.md`.

## 2026-07-23 — Policy regimes (exogenous backstop aggressiveness)  [55f2031, 9e123a8]
- Rebuilt the policy-regime feature on `main`: three exogenous ECB backstop regimes
  over the TPI feedback coefficient γ (`TPI_t = γ·(spread−spread_ss)`), plus a
  Stage-B belief lottery over an ex-ante-unknown CB type. All post-Jacobian numpy on
  cached `G_tpi`; production `main.py` untouched (code in `diagnostics/regimes/`).
- **The capital-key backstop compresses spreads** (`d(spread_rb)/d(cb_buy)=−1.95e-2`);
  aggressive/medium γ hit 78/117 bp peak vs 156 bp passive.
- Provenance: the opposite "SA-1 spread-*widening*" result from the retired
  `ms-regime` branch was an artifact of a superseded model (single-country conduit,
  par rule, `psi_lambda_B=2.8` in the breakdown region) — **void on `main`.**

## 2026-07-22 — C-1 fixed at root; psi_lambda_B → 1.1793; EL-1/PT-1 resolved  [53e1783, e73659d, af02685, 0c70882]
- **C-1 fixed at its root** [53e1783]: `steady_auxilliary_D/F` now solve `lambda_gk`
  from the multi-asset IC directly; `Delta_bD_D/F=0.2/0.4` are genuine hardcoded
  inputs (was a degenerate back-solve, `Δ_cross=1.45>1`). The "explosive EBA doom
  loop" was superseded, not merely improved — `b_gov_D[499]`≈1e-5 near-stationary.
  Fixed 3 downstream drift bugs found while verifying (stale audit harness, a
  diagnostic sign/scale error, a real TPI conduit accounting leak). Re-tested
  Finding F-1 with a validated **order-selected Prony** eigenvalue estimator (the
  earlier energy-ratio proxy was itself an overfitting artifact).
- **psi_lambda_B recalibrated** [e73659d]: 0.31 → **1.1284** to hit the 2010 GR-DE
  ~150 bp spread target (0.31 was only a bug-dodge, undershooting ~3×). Found a
  non-monotonic linear-approximation breakdown region above ψ≈2 — the old literature
  values 2.8/3.0 now sit inside it and must not be restored.
- **EL-1 + PT-1 resolved** [0c70882]: `recovery_rate` 0.00 → **0.30** (NPV Greek-PSI
  framing; Zettelmeyer–Trebesch–Gulati). `psi_lambda_B` re-tuned 1.1284 → **1.1793**
  to restore the 150 bp target after `EL_price` shrank. Pass-through (PT-1) validated
  at ≈−4.5%/100 bp, inside the Acharya–Drechsler–Schnabl range. **TPI loading rose to
  3.59/3.03/2.47** at γ=2/5/10 (over-compensated, declining — self-extinguishing
  premium holds).
- Retired `docs/FRAMING_HANDOFF.md` into SPEC/STATE/HANDOFF/CLAUDE [af02685].

## 2026-07-21 — EBA-2011 calibration goes live  [eade414, e003284, f43f3a7, dd704eb, 1ce5312]
- Ported EBA 2011 stress-test bank-sovereign concentration (`phi_bD_D=2.39`,
  `phi_bF_F=2.76`; thin bank net worth `n_inter_D=0.408`) and added a passive
  capital-intermediation fund (`omega_K_D=0.0601`, `omega_K_F=0.0190`) to reconcile
  thin EBA net worth with the aggregate capital stock (Walras-neutral; no-op at
  `omega_K=1`).
- Correct-signed doom loop achieved at `psi_lambda_B=0.31`, `mv_rule=1`,
  `phi_lamb=0.6` [e003284]. Documented the endogenous-`omega_K` negative result and
  the C-1 root cause analytically (ruling out EBA leverage as `theta`).

## 2026-07-13–17 — Fundamental EL channel, ECB conduit, EBA spec  [c6f5707, 995a957, 1e82e22, 738311d, abcbb6e, b25cd1f, c981f19, a6961ab, 7d560ed, a129238]
- Added the fundamental expected-loss channel (`EL_price`) to the sovereign bond FOC
  (macro-pru-fix) [c6f5707], independent of `psi_lambda_B`; disciplined `psi_lambda_B`
  by a moment sweep and a `psi_lambda_B=0` output experiment [995a957, 1e82e22,
  738311d].
- **ECB balance sheet as a capital-key conduit** [abcbb6e]: CB net cash flow split
  `kappa_cb_F=0.929` to the F treasury / rest to D; two-leg (carry + credit) P&L;
  off-path expected loss hand-computed (never read off the linear path). Figure 8
  insurance-loading schedule + captions baked into all TPI figures [b25cd1f]
  (loading then 4.86→4.06→3.22 at the ψ=0.31 calibration).
- EBA-calibration design spec + a moments module decoding the 2011 stress-test CSV
  into calibration targets [a6961ab, c981f19]; 2010 base-year decision and the
  `omega_K` fork resolution [7d560ed].
- **Retired `model_v12.ipynb`; designated `code/main.py` as the production
  pipeline** [a129238].

## 2026-06-23–29 — Market-value fiscal rule (Finding F-1)  [d4fa259, f25dcf6]
- Added the switchable **market-value fiscal rule** (`mv_rule`): the Bohn rule can
  react to the mark-to-market debt gap `q_b·b_gov(-1) − mv_gov_ss` instead of par —
  it "sees" the current spread, restoring stationarity with empirical long-duration
  bonds under the risk-premium framing. Documented as Finding F-1 across nav docs.

## 2026-06-22 — Modular reorganization (PR #28)  [7037d72, ab52722, e6affd0]
- Reorganized the monolithic `model_v12.ipynb` into modular Python files
  (`code/equations_{D,F,global}.py`, `main.py`, `steady_state.py`,
  `ic_delta_calibration.py`, `depreciation_calibration.py`, `full_model.py`,
  `tpi.py`, plot modules); removed stale notebooks/folders. `main.py` becomes the
  single source of truth. Fixed stale branch/PR references after the audit merge.

## 2026-06-10–12 — Forensic audit + six structural fixes (PR #27)  [396cbd9, 4c810e1, 52f17d5, a0ddc18, 1e68440]
- Hostile-referee forensic audit: the model converged but was **internally
  inconsistent**. Six fixes, each verified numerically (`docs/audit.md`):
  - **T-2** (critical) — deposit rate re-dated to `Rgross=(1+rdep(-1))·P(-1)/P` (was a
    period-t unknown on the t−1 deposit stock, ~9× quarterly GDP). The doom loop had
    run *backwards*: bank net worth and GDP *rose* after a sovereign-default shock.
    Required `phi_lamb` 0.02→0.15.
  - **W-2/W-3** — p-conversion of F-bank's D-good-denominated bond book
    (`bank_return_F`, `divert_bond_foc_F`); `goods_mkt_F` had leaked up to **2% of
    GDP**, corrupting all cross-country spillovers by 39–124%.
  - **W-1** — capital-timing: keep `Y=F(K_t)`, pay `mpk·ΔK` to the capital producer so
    **CA=ΔNFA** holds at first order [4c810e1].
  - **TPI-1** — added the CB budget constraint (remittance); unbacked closed-loop flows
    had inflated welfare gains ~40% at γ=10.
  - **A-2** — aligned `m` vs Φ/T bookkeeping between SS and dynamics (needed before any
    `chi1≠0` experiment).
- Acceptance (**passed**): max|ca_res_D|,|goods_mkt_F| ≤1e-7 on all shocks incl. TPI;
  `n_inter_D` and `Y_D` fall on a default shock; system stationary. Post-fix
  notebook-era TPI: ΔW_D +1.88 / ΔW_F −1.90 at γ=10 (later superseded). Left open at
  audit: **S-1** (writeoff regime), **C-1** (`Δ_cross>1`).

## 2026-06-03–10 — TPI experiment, Walras-leak fixes, bond-pricing inversion  [f810776, 6aaa9ec, 25caad4, 754ea86, ecbae65, 09cf30e]
- Added the **TPI closed-loop spread rule** with welfare figures (PR #25) [754ea86];
  endogenous default rule + lump-sum tax fiscal rule [ecbae65].
- Fixed a Walras leak by adding `cap_profit` to the resource constraint [f810776];
  **inverted bond pricing** so `q_b` is the forward-looking PV and `rb` the implied
  YTM [6aaa9ec]. Precomputed the Jacobian `G` and got IRFs via matrix multiply
  [25caad4]. Added `spread_rb`/welfare aggregates; "sovereign risk without actual
  write-off" via `writeoff_enabled` [09cf30e].
- Discovered and reverted a "perverse GDP rise under default" fix [985b0b6→670c0ce] —
  the precursor symptom later diagnosed as T-2.

## 2026-05-21–28 — Multi-asset GK IC, endogenous default, macroprudential charge  [353c053, 5d7227a, a0c0828, 19df483, 26a223d, afd3bf8, b15811d, 6214e84]
- Built the **multi-asset Gertler-Karadi incentive constraint** (three-ν,
  asset-specific diversion rates) and a canonical-GK steady state (Ω, `f` convention,
  `theta` retune, `lambda_gk` calibration).
- Added **endogenous default** via a smooth-power debt-gap response (`def_scale`), an
  endogenous sovereign risk weight (lagged default + rating shock), and a Basel-style
  `mp_wedge` sovereign capital charge. Merged several topic branches (PRs #5–14).
- Forked `model_v12` as the no-central-bank baseline.

## 2026-05-14–20 — Foundations  [4aa2155, PR #1–3, 2327f4c, 19df483, 9385c05]
- Initial repo, README, and two-country model assembled from merged predecessor
  notebooks (v10/v11).
- Government default (PR #2); bond pricing; cross-border excess-return anchors; a
  spread-miscalibration fix [9385c05] (had targeted `rdep_F−rdep_D` not `rb−rdep`,
  linearizing around a non-stationary point → ~60/12 bp errors).
- nbstripout/nbdime notebook hygiene (PR #3); `zeta_writeoff` dial; endogenous-default
  channel via `def_scale`.

---

## Predecessor models (pre-git / `OLD models/`)
- `model_v11.ipynb` — predecessor with free bond trade between intermediaries.
- `model_v12` — added NK labour, portfolio adjustment costs, and the TPI extension;
  reorganized into `code/main.py` (PR #28) and then removed.

*(This section folded in from the former `docs/PROCESS.md` "Version history".)*

---

**Sources:** `git log`, `docs/STATE.md`, `docs/audit.md`, `docs/eba_calibration.md`,
`docs/HANDOFF.md`.
