# EBA-anchored calibration

Source: EBA **2011 EU-wide stress test** disclosure (`data/DATA_DISCLOSURE.CSV`),
base **31 Dec 2010 actual**; nominal GDP from Eurostat `nama_10_gdp` (2010);
sovereign yields from FRED via `Empirics/outputs/spreads_fred.csv`.
Country map: **D = Greece** (banks GR030–GR035, issuer GR), **F = Germany**
(banks DE0xx, issuer DE). Moments produced by `code/eba_calibration.py` →
`data/eba_moments.json`. Code map validated: `DE017 = Deutsche Bank`
(2010 RWA €346.6 bn, CT1 €30.4 bn, 8.76% — matches published).

---

# REBUILD (2026-07-31) — supersedes everything below

The 2026-07-22 EBA build was reverted on 2026-07-30 because it was **not
identified**: it matched a set of EBA moments by inventing a free parameter to
absorb the resulting over-identification. This rebuild fixes that. Everything
from "Parameter → moment map" down is the historical record of the old build.

## What was wrong before

| # | Failure | Fix |
|---|---|---|
| 1 | **`omega_K` was a free plug.** `theta=4.0` was *assumed*; `omega_K` was then back-solved as `N(θ−φ_own−φ_cross)/(Q·K_target)` to force `K=10.8`. A parameter with no data counterpart absorbing an assumption. | `theta` is now **measured** (CT1 / GK-eligible assets) and `omega_K` is **measured** (corporate+CRE EAD / K). `K` becomes an *output* and the 10.8 comparison becomes a real over-identifying check. |
| 2 | **`delta_b` had no EBA counterpart at all.** Only the `MATURITY_CODE=999` total row was read; the target `0.036/0.038` came from the sovereign's *whole outstanding stock* average maturity — the wrong object. | The **maturity ladder** (`MATURITY_CODE` 125–155) is now read and repriced. |
| 3 | **`psi_lambda_B` absorbed a mis-calibrated mechanical channel** and its moment-match was non-monotone (219bp@2.0, 97bp@2.6, 853bp@2.8). | The mechanical MTM channel is now measured, so `psi_lambda_B` only has to supply genuine amplification. It remains unidentified — see the ledger. |
| 4 | **Amplification risked being matched to adverse-scenario CT1 depletion.** | Explicitly rejected and guarded by a test. The 2011 adverse scenario excluded sovereign default in the banking book (this is why Dexia passed and failed months later); its depletion understates sovereign pass-through *by construction*. |

## The maturity ladder → `delta_b`

The EBA 2011 sovereign template reports gross direct long exposure in seven
residual-maturity buckets (upper bounds 3M/1Y/2Y/3Y/5Y/10Y/15Y). They sum
exactly to the reported total (asserted in `validate()`), so the last bucket is
effectively "10Y and over".

| Book | Face €m | Wtd-avg residual maturity | Modified duration | → `delta_b` |
|---|---|---|---|---|
| GR banks / GGB | 54,447 | 5.13y | **3.12y** @ 12.01% | **0.0777** |
| DE banks / Bund | 315,313 | 4.86y | **4.22y** @ 2.91% | **0.0568** |
| DE banks / GGB (contagion) | 7,934 | 7.41y | 4.05y | 0.0593 |
| GR banks / Bund | 411 | 7.27y | 6.45y | 0.0363 |

The gap between maturity and duration is the whole point: at Greece's 12%
end-2010 yield, a 5.13y-maturity ladder has only **3.12y** of duration. Duration,
not maturity, is what governs the mark-to-market loss `delta_b` exists to
generate. Durations are inverted through the Hatchondo–Martinez perpetuity
(`−dlog q_b/dr = 1/(r_q + delta_b)` quarters).

**This retires the `delta_b` / F-1 blocker.** The standing item was "port
`delta_b=0.036/0.038` (7y/6.5y) from `bank-cal`", which F-1 showed requires
`mv_rule=1` *and* `phi_lamb=0.60` jointly, with a hard break at `phi_lamb=0.15`.
The measured value is **0.0777/0.0568**, much closer to the old `0.10` than to
`0.036` — so the par rule (`mv_rule=0`) is not pushed into F-1's explosive
region and the two-parameter move is not needed. The old target was measuring
the wrong object.

## `theta` and `omega_K` from the same balance sheet

`theta` multiplies only the GK book, so the denominator is **GK-eligible
assets** = corporate (ex-CRE) + commercial real estate + sovereign EAD, own
country. Residential mortgages are excluded (the model has no housing stock).

| | Corp ex-CRE | CRE | Sovereign | GK assets | CT1 | **θ** | (CT1/TA memo) |
|---|---|---|---|---|---|---|---|
| D (GR) | 57,254 | 13,410 | 54,858 | 125,522 | 22,778 | **5.51** | 14.9 |
| F (DE) | 352,593 | 117,679 | 323,247 | 793,519 | 114,317 | **6.94** | 32.9 |

The total-assets version (14.9 / 32.9) stays rejected — it includes interbank,
reserves and retail, and was previously verified not to converge in steady state.

`omega_K` = (corporate + CRE) / K with `K = 2.7 × annual GDP`: **0.117 (D)**,
**0.067 (F)**. Banks fund roughly a tenth of the productive capital stock; the
rest is equity, retained earnings and non-bank credit. Setting `omega_K=1` would
assert banks fund *all* capital, which is the counterfactual that generated the
original over-identification — so `omega_K` is **kept, not dropped**, but it is
now measured rather than back-solved.

**Over-identifying check.** With θ, `n_inter`, `omega_K` and the φ ratios all
measured independently, the balance sheet `omega_K·Q·K + q_b·bonds = θ·N`
implies `K_D = K_F = 10.80` — against the conventional `K/Y_annual = 2.7`
target of 10.8. `steady_state.py` prints this every run.

## The amplification moment: exposure ladder × observed yield moves

Acharya–Steffen ("The Greatest Carry Trade Ever") construction — bank-level
mark-to-market loss from the observed sovereign repricing, **not** the stress
test's own scenario output.

GGB 10y: 12.01% (Dec-2010) → 21.14% (Dec-2011), **+913bp**. Bund: 2.91% → 1.93%.

| | Mechanical `dNW/CT1` per 100bp | Realised 2011 |
|---|---|---|
| D (GR banks, GGB book) | **−5.73%** | **−39.8% of CT1** |
| F (DE banks, Bund book) | −11.92% | +12.3% of CT1 (Bund rally) |

This is the headline diagnostic of the rebuild:

- The measured mechanical channel is **−5.73%/100bp**. The pre-EBA calibration
  (`phi_own=0.25`, `delta_b=0.10`) could generate only **−0.61%/100bp** — an
  order of magnitude too weak. Setting `phi_own=2.39` with the *old* `delta_b`
  already gives −5.85%, so essentially the entire gap was the concentration
  ratio, not the duration.
- That explains S-1's "expected loss is 10.9% of the default loading, collateral
  friction 89%" split as a **calibration artifact**: `psi_lambda_B=3.0` was
  standing in for a mechanical mark-to-market loss that had been mis-calibrated
  by ~10×. With the mechanical channel measured, `psi_lambda_B` only has to
  supply true amplification.
- Only **11.1%** of GR banks' GGB book was fair-valued (AFS+FVO+trading); the
  rest sat at amortised cost. The −39.8% figure is therefore an *economic*, not
  an accounting, loss — which is the right concept, because the model's
  `n_inter` is economic net worth and marks the whole book through `q_b`.

## GK feasibility — the measured concentration BOUNDS `Delta_own`

**The most consequential finding of the rebuild.** Feeding the measured moments
into the model at the inherited `Delta_bD_D=0.2 / Delta_bF_D=0.4` produces a
steady state with

```
lambda_gk_D = -0.0869    Omega_D = -0.3013
lambda_gk_F = -0.0723    Omega_F = -0.3217
```

— a **negative IC multiplier and negative banker franchise value**. The solver
still converges, every Walras residual is machine-zero (`ca_res_D = -7.7e-17`),
the IC-δ consistency check passes exactly, the stability check passes
(`b_gov_D[499]=2e-6`), and the TPI loading schedule still declines. **All of it
is meaningless** — the banker's continuation value is negative. This is the C-1
failure mode in a new guise: silent degeneracy that passes every check the
pipeline previously ran.

From `steady_auxilliary_D/F`, `lambda_gk > 0` requires a positive denominator:

```
f * theta  >  (1 - Delta_own) * phi_own + (1 - Delta_cross) * phi_cross
```

The banker's franchise value (left) must cover the non-divertable "good
collateral" part of the sovereign book (right). At the measured moments:

| | `f*theta` | required | margin |
|---|---|---|---|
| D | 0.6613 | 1.9231 | **−1.2618 VIOLATED** |
| F | 0.8330 | 2.2482 | **−1.4153 VIOLATED** |

**This is an identification result, not just a bug.** The measured concentration
`phi_own = 2.39 / 2.76` puts a hard **lower bound** on `Delta_own`:

| `Delta_cross` | min `Delta_own` (D) | min `Delta_own` (F) |
|---|---|---|
| 0.4 | 0.7279 | 0.7131 |
| 0.6 | 0.7264 | 0.7081 |
| 0.8 | 0.7249 | 0.7030 |
| 0.9 | 0.7241 | 0.7005 |

The bound is essentially independent of `Delta_cross` (cross-holdings are tiny).
So EBA data does not point-identify `Delta_own`, but it **rules out the entire
region below ~0.73 — including the 0.2 that was hardcoded.** A parameter the
previous ledger called wholly unidentified turns out to be bounded by data.

The other two levers are out of range:

- `f` would need **> 0.349** (GK literature: 0.03–0.12).
- `theta` would need **> 16.03** (measured 5.51; even the rejected CT1/total-assets
  14.9 falls short, and that variant was already verified not to converge).

### What had pinned `Delta` at 0.2/0.4: a hidden `ratio = 2.0` convention

`ic_delta_calibration._ic_delta` back-solved `Delta` from one equation in two
unknowns, closing the system with a hardcoded
**`ratio = Delta_cross / Delta_own = 2.0`**. That was an undocumented structural
convention masquerading as a consistency check — it is exactly why the inherited
`0.2/0.4` pair "passed". With `Delta_cross <= 1` it also caps `Delta_own <= 0.5`,
below the ~0.73 the measured moments require.

**Fixed 2026-07-31.** The convention is gone. `Delta_own` and `Delta_cross` are
free structural parameters, and the module now checks the thing that matters —
the IC residual

```
value  ==  lambda_gk * [theta - (1-Delta_own)*phi_own - (1-Delta_cross)*phi_cross]
```

a genuine one-equation residual with no free closure (`ic_residual`, tolerance
1e-8, plus a positive-divertable-leverage check).

### Resolution: `Delta = 0.85 / 0.90`

Sweeping `Delta_own` with `Delta_cross = Delta_own + 0.05`:

| `Delta_own`/`Delta_cross` | `lambda_gk_D` | `lambda_gk_F` | `Omega_D` | `Omega_F` | `K_D` | `K_F` |
|---|---|---|---|---|---|---|
| 0.85 / 0.90 | **+0.927** | **+0.960** | +4.62 | +5.98 | 10.80 | 10.65 |
| 0.90 / 0.95 | +0.488 | +0.456 | +2.49 | +2.91 | 10.80 | 10.65 |

**0.85/0.90 is committed for the EBA branch.** `lambda_gk_D = +0.927` is
essentially identical to the pre-EBA `+0.923`, so the amplification block keeps
its previous strength while the concentration becomes measured — the cleanest
possible basis for comparing the two calibrations.

Note the closed-form bound (~0.73) understates the requirement because it ignores
the endogenous banker return `rn`, which enters as `D_target/(beta_inter*(1+rn))`.
`Delta = 0.80/0.90` lands just past the `lambda_gk` **pole** (the denominator
crosses zero, so `lambda_gk` runs `-inf -> +inf`) with `lambda_gk_F = -12.45`.
The real frontier is a fixed point, not a formula — sweep, don't solve.

**Why this is a correction and not a fudge.** Measured leverage is only 5.5× on a
book that is ~43% sovereign. You cannot simultaneously claim sovereigns are
excellent collateral: if they were, the bank would lever further and `theta=5.5`
would not be the binding constraint. High concentration at low leverage *implies*
bonds are nearly as divertable as capital. That is also the right story for
2010–12 Greece — collapsing GGB collateral eligibility, rising ECB haircuts.

## Dynamic instability — the remaining blocker

With the steady state correct, the **dynamics are still explosive**:
`b_gov_D[499] ~ 1e2-1e3` against a ~1e-5 target. This is a *separate* problem
from the collateral mapping, and it is structural.

The measured moments give a financial-accelerator gain

```
theta * phi_own = 5.51 * 2.39 = 13.17     (placeholder: 4.0 * 0.25 = 1.0)
```

— a ~13× stronger leverage loop. Diagnosed 2026-07-31; it is **not** the fiscal
mode and **not** the collateral friction:

| Test | Result |
|---|---|
| `psi_lambda_B = 0` (friction fully off) | still explosive, `b_gov[499] = -2038` |
| `phi_lamb` 0.6 → 25 | **flat**: peak spread ~1.1e7bp at 0.6, 1.5 *and* 25 |
| `mv_rule = 1` | does not fix it |
| `chi1` 0 → 0.5 | peak spread **1.1e7bp → 6.0bp**, `b_gov[499]` −2038 → +70 |
| `chi1` ∈ [0.2, 5.0] | `b_gov[499]` stays 70–560 — amplitude damped, root not removed |

`chi1` (the Auclert intermediary capital adjustment cost, currently **0**) is by
far the strongest lever and makes the spread response sane again, but no value
tested removes the unstable root.

**Open routes:**

1. **Damp the accelerator structurally** — `chi1` is the natural home but needs
   pairing with something else (slower `theta` adjustment, or a leverage rule
   with inertia).
2. **Re-scope the bank block** so `phi_own` is not 2.39 — model only the
   sovereign-exposed sub-book, or let `n_inter` be broader than stress-test CT1.
   Keeps the mechanism, weakens the "EBA-measured" claim.
3. **Replace the GK IC** with a constraint tolerating high concentration
   (value-at-risk / risk-weighted rather than linear divertability).

**Committed: `EBA_CALIBRATION = False`** — the pre-EBA calibration, which solves.
The switch turns the whole measured moment set on in one line, and
`assert_gk_well_posed` guarantees the steady-state failure mode can never return
silently.

**The headline.** Measured Greek bank-sovereign concentration implies a
financial-accelerator gain ~13× the placeholder's, and this model is linearly
unstable there regardless of the fiscal rule or the collateral friction. The
2026-07-22 build hid the steady-state half of this by back-solving `omega_K`
around an assumed `theta`; the 2026-07-30 revert avoided both halves by
abandoning EBA anchoring. Neither is a reason the concentration is wrong — it is
measured — so the model, not the data, is what has to give.

**Guarded in code.** `steady_state.assert_gk_well_posed` now runs inside
`_apply_ss_anchors`, i.e. on every solved steady state in both
`steady_state.py` and `depreciation_calibration.py`. It raises with the
feasibility inequality and the suggested levers. `gk_feasibility_margin` and
`min_Delta_own` are exported for sweeps. **This guard is the single most
valuable artifact of the rebuild** — it makes the C-1 class of failure
impossible to commit silently.

## Identification ledger

Emitted into `data/eba_moments.json` under `identification`, and mirrored here.

**Identified by EBA 2011:**

| Parameter | Moment |
|---|---|
| `n_inter_D/F` | CT1 / own quarterly nominal GDP → 0.4075 / 0.1748 |
| `phi_bD_D_ss`, `phi_bF_F_ss` | own-sovereign gross long / CT1 → 2.390 / 2.758 |
| `phi_bD_F_ss`, `phi_bF_D_ss` | cross-border sovereign / CT1 → 0.069 / 0.018 |
| `theta_D/F` | (corp ex-CRE + CRE + sovereign) EAD / CT1 → 5.51 / 6.94 |
| `delta_b_D/F` | ladder modified duration at the end-2010 market yield → 0.0777 / 0.0568 |
| `B_supply_D/F` | bank-held sovereign / own quarterly GDP → 1.116 / 0.483 |

**Identified jointly with one standard macro target:**

| Parameter | Moment |
|---|---|
| `omega_K_D/F` | measured corp+CRE EAD ÷ K, with `K/Y_annual = 2.7` → 0.117 / 0.067 |

**Still NOT identified — the honest list:**

| Parameter | Status |
|---|---|
| `Delta_bD_D`, `Delta_bF_F` (own) | **Partially identified — bounded below, not point-identified.** GK feasibility at the measured `theta`/`phi_own`/`f` forces `Delta_own > ~0.73 (D) / ~0.71 (F)`; see "GK feasibility" above. Committed at 0.80. The level above the bound is an author decision. |
| `Delta_bF_D`, `Delta_bD_F` (cross) | Still unidentified. Committed at 0.90 to preserve `Delta_own < Delta_cross`. The feasibility bound barely constrains these (cross-holdings are ~1% of the book), so they are close to free. |
| `psi_lambda_B_D/F` | Tuned to the 150bp GR–DE spread target. The mechanical channel is now measured, so this parameter does far less work than before — but its level is still one moment, one parameter. **Proper identification needs bank equity returns regressed on the EBA exposure cross-section (Acharya–Steffen). Those returns are not in this repo.** The exposure cross-section *is* (per-bank columns in the disclosure), so this is a data-acquisition task, not a modelling one. |
| `def_scale_D/F` | 0.25, hand-set. Exceeds the 2011 GR crisis peak (0.12–0.23). |
| `f_D/F` | 0.12, GK literature. `bank-cal` has 0.03. Not an EBA object. |
| `phi_lamb_D/F` | 0.60, ~Bohn (1998). Chosen for stationarity under the measured doom loop, not moment-matched. |
| `recovery_rate_D/F` | 0.30, Greek-PSI NPV literature (EL-1). Not EBA. |
| `theta` own-country vs total book | Own-country (5.51/6.94) is used. The total-book alternative (7.21/13.30) is equally defensible for F, which lends heavily outside the euro area. Reported under `leverage_alternatives`. |

**Deliberately rejected:**

- **Adverse-scenario (105) CT1 depletion** — the 2011 exercise excluded
  sovereign default in the banking book, so its capital depletion is not a
  measure of sovereign-stress pass-through. Guarded by
  `test_adverse_scenario_not_used`.
- **`theta` from total assets** — 14.9 / 32.9; includes assets `theta` does not
  multiply; previously verified not to converge.

## Conventions that are choices, not measurements

Recorded so they can be argued with rather than discovered later:

1. **Discount rate for the ladder duration.** Durations are computed at the
   **31-Dec-2010 market yield**, so the model reproduces the *realised* MTM
   sensitivity. Evaluating at the model's own SS yield (1% annual) instead gives
   longer durations (4.88y GGB / 4.67y Bund → `delta_b` 0.0489/0.0511), because
   the model's SS has no sovereign risk premium. The ~40% gap between the two is
   a real modelling choice.
2. **Bucket midpoints.** Each bucket is represented by its midpoint; the last
   ("10Y–15Y") uses 12.5y.
3. **Coupons.** GGB 4.7%, Bund 3.5% (legacy stock outstanding at end-2010).
   `coupon_sensitivity_modified_duration_b_D_D` in the JSON reports the effect.
4. **Full book vs fair-valued book.** The MTM moment marks the *whole* book,
   matching the model's economic-net-worth concept, not 2011 accounting.

## Unclosed check

`omega_K` implies bank credit to non-financial corporations of ~12% (D) / ~7%
(F) of the capital stock. This has an independent counterpart in ECB BSI
statistics that has **not** been pulled. If it disagrees materially, `omega_K`
and hence `theta` need revisiting.

---

# Historical: the 2026-07-22 build (superseded by the rebuild above)

## Parameter → moment map

| Parameter | Old | EBA value | Moment (EUR m, 2010) |
|---|---|---|---|
| `n_inter_D` | 3.0 | **0.408** | CT1_GR 22,778 / quarterly GDP 55,898 |
| `n_inter_F` | 3.0 | **0.175** | CT1_DE 114,317 / quarterly GDP 653,815 |
| `phi_bD_D_ss` (GR banks' GR book / capital) | 0.25 | **2.39** | 54,447 / 22,778 |
| `phi_bF_F_ss` (DE banks' Bund / capital) | 0.25 | **2.76** | 315,313 / 114,317 |
| `phi_bD_F_ss` (DE banks' GR / capital, contagion) | 0.15 | **0.069** | 7,934 / 114,317 |
| `phi_bF_D_ss` (GR banks' Bund / capital) | 0.15 | **0.018** | 411 / 22,778 |
| `B_supply_D`, `b_gov_D` (bank-held GR debt) | 2.4 | **1.19** | 62,381 / qGDP (27.9% ann. GDP) |
| `B_supply_F`, `b_gov_F` (bank-held Bund) | 2.4 | **0.591** | 323,247 / qGDP (12.1% ann. GDP) |

`phi` = q·(sovereign book)/(bank net worth). The own-holding ratios (2.39×, 2.76×)
are the doom-loop vulnerability; cross-holdings (0.069, 0.018) are the thin direct
contagion channel — a **10–30× reduction** from the previous 0.25 placeholder.

**Design note (b1):** `B_supply` is the *bank-held* sovereign stock (~28%/12% of GDP),
NOT headline debt/GDP (~150% GR). The non-bank/official/ECB residual is intended to
enter via the SMP conduit (phase 2). PSI recovery (Stage 3) is deferred.

## Structural change: `omega_K` capital-intermediation fund

EBA over-identifies the bank block: matching the 2.39× sovereign/capital ratio needs
thin net worth (`n_inter_D=0.408`), but the model ties net worth to intermediating the
*entire* capital stock. Resolution (user-approved option A): banks hold `omega_K` of
capital, a **passive deposit-funded capital fund** holds `(1-omega_K)K` and rebates its
spread `(rk-rdep)` to households as `div_fund` (`capital_fund_D/F`, into `income_D/F`).
`omega_K` set to preserve `K≈10.8` (K/annualY≈2.7): `omega_K = N·(θ-φ_own-φ_cross)/(Q·K_target)`
→ **0.060 (D), 0.019 (F)**.

**Verification:** the change is Walras-clean to machine precision at SS **and** along
IRFs (`ca_res_D`, `goods_mkt_F` ~1e-8; solved targets ~1e-16). Verified as a no-op at
`omega_K=1` (reproduces the pre-EBA SS exactly).

## Stability: what works, what's open

**Update (2026-07-22): C-1 is now fixed at its root (see below), and the explosive
root reported below turned out to be a symptom of C-1, not of the EBA concentration
itself.** With the multi-asset `lambda_gk` fix in place, the same calibration
(`psi_lambda_B=0.31`, `mv_rule=1`, `phi_lamb=0.60`, unchanged) is **near-stationary**:
`b_gov_D[499]` is `~1e-5` on both the TFP and default shocks (was 3.4–79.5 across the
variants tested below), and `ρ_b (partial-eq.) = 0.373` (well under the 0.95 target).
The section below is kept as the historical record of the pre-fix diagnosis; see
"C-1 fix (2026-07-22)" further down for what changed and why.

**`psi_lambda_B` was subsequently recalibrated 0.31 → 1.1284 later the same day**
(see "Remaining open items" below) -- 0.31 was never itself a calibration target,
only a value that avoided the C-1 degeneracy. Re-verified stationary at the new
value too (`b_gov_D[499] ~1e-6/1e-7`, `ρ_b = 0.373`, unchanged) -- the stability
finding above is not an artifact of the specific (stale) `psi_lambda_B`.

With `psi_lambda_B` rescaled **3.0 → 0.31** (doom-loop amplification ~ `psi·phi_own`;
`phi_own` rose ~10×, so `psi` falls ~10× to keep comparable amplification and keep
`Delta_eff = Delta + psi·def_rate` below 1), plus `mv_rule=1` and `phi_lamb=0.60`:

- ✅ **Doom-loop signs correct**: a default shock gives spread ↑, `q_b` ↓, `rb` ↓,
  bank net worth ↓ (previously inverted by the C-1 `Delta→1` collateral flip).
- ✅ **Walras clean** on all IRFs.
- ⚠️ ~~Residual ~1.2%/period explosive root~~ **RESOLVED by the C-1 fix** (`K_D`,
  `b_gov_D` drift up to t=499). `phi_lamb` reduces the amplitude (`K_D[499]` 9.5→3.4
  as `phi_lamb` 0.3→0.6) but not the root -- *this diagnosis turned out to be
  incomplete; the root was the C-1 degeneracy, not the concentration ratio itself.*

### Why C-1 forces `Delta→1` (analytical; explains the collateral flip)

C-1 is **not** a calibration input — it is an artifact of an inconsistency between the SS
and dynamic blocks. `steady_auxilliary` builds `lambda_gk` from the *single-asset* GK
formula (`lambda_gk = f/(θ·(1/(β_inter(1+rn)) − (1−f)))`), which forces the identity
`value/lambda_gk = θ`. The *dynamic* IC is multi-asset with divertability weights `Delta`.
Substituting the identity into `ic_delta_calibration._ic_delta` collapses the back-solve to
a function of the portfolio ratio alone:

```
Delta_own = (phi_own + phi_cross) / (phi_own + ratio·phi_cross)        [ratio = 2.0, hardcoded]
```

Verified exactly against all three observed cases:

| case | phi_own | phi_cross | predicted `Delta_cross` | observed |
|---|---|---|---|---|
| baseline (pre-EBA) | 0.250 | 0.150 | 1.4545 | **1.45** (C-1 in CLAUDE.md) |
| EBA D-bank | 2.390 | 0.018 | 1.9852 | **1.9852** |
| EBA F-bank | 2.760 | 0.069 | 1.9524 | **1.9524** |

So `Delta_cross > 1` whenever `Delta_own > 0.5` — i.e. essentially always. The pre-EBA
calibration merely sat far enough from the boundary (`Delta_own=0.73`) to limp; EBA's
realistic concentration (`phi_own=2.39`) drives `Delta_own→0.99`, and `Delta_eff = Delta +
psi·def_rate` then crosses 1, flipping the sign of the collateral channel. **This is exactly
why `psi_lambda_B` had to be rescaled 3.0 → 0.31** to restore correct signs.

*Permanent fix -- DONE (2026-07-22):* `steady_auxilliary_D/F` (`code/equations_D.py`,
`code/equations_F.py`) now solve `lambda_gk` from the multi-asset IC directly. The
"fixed point" collapses to closed form because `value` is *linear* in `lambda_gk`
(via `Omega = f + (1-f)*lambda_gk*theta`):

```
D_target = theta - (1-Delta_own)*phi_own - (1-Delta_cross)*phi_cross
lambda_gk = f / (D_target/(beta_inter*(1+rn)) - (1-f)*theta)
```

This is the *original* single-asset formula with `theta` replaced by `D_target` in
the first term only -- `D_target == theta` (recovering the old formula exactly)
iff `Delta_own == Delta_cross == 1`, confirming the old formula silently assumed
full (`Delta=1`) bond divertability. `Delta_bD_D=0.2` / `Delta_bF_D=0.4` (D) and
`Delta_bF_F=0.2` / `Delta_bD_F=0.4` (F) are now genuine hardcoded calibration
inputs (already sitting as placeholders in `code/calibration.py`, previously
overwritten by the degenerate back-solve). `code/ic_delta_calibration.py`'s
back-solve is no longer authoritative -- it now re-derives the *same* equation as
a consistency check and asserts it recovers the hardcoded Delta (verified to
recover 0.20000000000000157 / 0.40000000000000313 etc., i.e. exact to
floating-point). `smart_steady_D/F` were extended to also return `phi_bD_D`/
`phi_bF_D` (D) and `phi_bF_F`/`phi_bD_F` (F) so `steady_auxilliary_*` can consume
them -- SSJ wires `@simple` block outputs to same-named inputs automatically, but
note its output-name parser is a single-line regex over the `return` statement
(`utilities/function.py::output_list`), so multi-line `return (...)` tuples
silently break block wiring with an opaque `ValueError: ... is output twice`. Keep
all `@simple` return statements on one line.

**Consequence for stability:** this also resolved the "explosive root" reported
below as intrinsic to the EBA calibration -- see the update note at the top of
this section. `run_audit.py`'s IC/Walras/sign checks (frozen pre-EBA calibration)
and `code/main.py`'s live EBA pipeline both confirmed the fix at the time.
(`audit_artifacts/` removed 2026-07-30 — that frozen-calibration divergence is
exactly why; `code/main.py`'s IC-δ consistency check is now the live check.)

*Ruled out:* using EBA **total-asset** leverage as `θ` (GR 16.56, DE 42.62). `θ` multiplies
only the GK book (capital+sovereign), whereas EBA total assets include low-yield
loans/reserves; `θ=16.56` on the model's fixed `rk−rdep=0.74%` spread implies ~52% annual
banker ROE and the SS **does not converge** (verified). Do not retry.

### Formerly "open issue": the EBA sovereign doom loop is intrinsically explosive here -- RESOLVED (2026-07-22)

*(Kept verbatim below as the historical record of the diagnosis at the time; the
conclusion it reaches turned out to be wrong -- see the correction after.)*

> Full stationarity is **not** achieved. `phi_lamb` (fiscal feedback) damps the
> amplitude (`K_D[499]` 9.5→3.4 as `phi_lamb` 0.3→0.6) but leaves a residual
> ~1.2%/period explosive root. Two capital-structure variants were tested;
> **neither** removes the root:
>
> - **Fixed `omega_K`** (bank holds `omega_K·K`): `K = (θN−bonds)/(omega_K·Q)`
>   levers `K` `1/omega_K≈16–67×` on thin net worth. *Least unstable of the two*
>   (default-shock `b_gov[499]=0.17` at `phi_lamb=0.6`). **This is the current
>   committed state.**
> - **Endogenous `omega_K` = fixed fund quantity** (bank holds `Kbank=K−Kfund_ss`):
>   tested (commit history) and made it **worse** (`b_gov[499]=79.5`). Reverted.
>
> **Conclusion:** the instability is not the capital plumbing — it is the **strong
> EBA sovereign doom loop itself** (`phi_bD_D=2.39`). Both variants are correctly
> signed and Walras-clean; both mildly explode.

**Correction:** that conclusion was wrong. The explosive root was **the C-1
degeneracy** (`lambda_gk` built from the single-asset formula, implicitly
`Delta=1`, while `Delta_eff` was being pushed near/over 1 by the `psi_lambda_B`
collateral-friction term) -- not the 2.39× concentration ratio itself. With the
multi-asset `lambda_gk` fix (previous section) and **no other calibration change**
(`psi_lambda_B=0.31`, `mv_rule=1`, `phi_lamb=0.60` all unchanged), `code/main.py`'s
live EBA pipeline gives:

```
b_gov_D[499] on TFP shock:      -0.000010   (was 3.4-79.5 across variants above)
b_gov_D[499] on default shock:   0.000002
rho_b (partial-eq.) = 0.373     (target < 0.95)
```

i.e. **near-stationary to numerical noise**, not merely "less explosive." The 2.39×
concentration ratio is not, by itself, destabilizing; it was the internally
inconsistent collateral value (bonds silently priced as if `Delta=1`, a single-asset
artifact) amplifying itself through the leverage/IC loop. The 4 "options for the
author" from the original diagnosis (relax the 2.39× target, microfound the fund,
etc.) are **moot** -- none were needed. This does *not* mean the SMP-conduit
stabiliser story is wrong, only that it is no longer required to explain why the
baseline doom loop is stable; the conduit experiment can now be evaluated as a
genuine "how does welfare/spread respond to CB purchases" question on a system that
is already well-behaved without it, rather than as an emergency stabiliser for a
broken one.

### TPI conduit accounting leak -- found and fixed alongside C-1 (2026-07-22)

Running the TPI (ECB capital-key conduit) experiment on the EBA calibration surfaced
a second, independent bug: `code/tpi.py` builds its own copy of the full dynamic
model (`ha_full_tpi`) rather than reusing `full_model.py`'s, and that copy predates
the `omega_K` capital-fund commits -- it was missing the `capital_fund_D`/
`capital_fund_F` blocks entirely (same "shared equations evolved, a downstream
duplicate didn't" pattern as the `run_audit.py` drift below). Effect: `div_fund_D`/
`div_fund_F` (the passive fund's rebate into household income) was silently frozen
at its steady-state level in the TPI model instead of responding to the shock, so
the TPI's goods-market/external-account accounting quietly diverged from the
baseline model's by exactly that omitted response.

Symptom (pre-fix, EBA calibration): `max|ca_res_D|` up to **6.2e-2** and
`max|goods_mkt_F|` up to **1.3e-3** at every `gamma` -- five to seven orders of
magnitude above the ~1e-8 floor the rest of the model holds to. The `G_tpi[cb=0]`
vs. baseline-Jacobian sanity check (should match to `<1e-8` since zeroing CB
purchases should exactly reproduce the no-TPI model) differed by `2.54e-3` --
**the same order of magnitude as the reported welfare gains themselves**. Any TPI
welfare number computed on the EBA calibration before this fix is not reliable --
the leak was large enough to be a first-order contributor to, not just noise
around, the result.

Fix: added `capital_fund_D`/`capital_fund_F` to `tpi.py`'s imports and to
`ha_full_tpi`'s block list, in the same position as `full_model.py`. Post-fix:
`G_tpi[cb=0]` vs baseline Jacobian `max|err| = 0.00e+00` (was `2.54e-3`);
`max|ca_res_D|` down to `~2-5e-8`, `max|goods_mkt_F|` down to `~6-7e-10` at every
gamma -- back in line with the rest of the model. The `loading` schedule (premium
PV / EL PV) is still declining in gamma post-fix (0.83 → 0.74 → 0.62 at
gamma=2/5/10), consistent with the self-extinguishing-premium narrative, but the
*levels* differ from anything reported before this fix (and before the EBA
calibration switch) -- treat all TPI welfare/loading numbers from before
2026-07-22 as superseded.

### Remaining open items (author decisions, not yet acted on)

- ~~**`psi_lambda_B` re-exploration**~~ **DONE (2026-07-22).** Re-ran the moment-
  matching exercise (`audit_artifacts/psilam_moment_sweep_postC1.py`, same
  external target as the pre-fix `2.8` calibration: 2010 GR-DE spread ~150bp on
  a 1pp default shock) on today's model. Result: `psi_lambda_B=1.1284`
  (verified: 151.3bp), replacing `0.31` (never a target, undershot 150bp by
  >3x) as the committed value. Neither `2.8` nor the original `3.0` transfers:
  the spread-vs-`psi_lambda_B` response is smooth/monotonic only up to
  `psi_lambda_B~1.5-2.0`, then turns wildly non-monotonic (219bp at 2.0, 97bp
  at 2.6, 853bp at 2.8, 353bp at 3.0, 11478bp at 5.0) -- both legacy values now
  sit inside a linear-approximation-breakdown region on this model, not a
  region of valid economic moments. (EBA's thinner bank net worth pulls this
  breakdown much earlier than the old calibration's -- there it only appeared
  around `psi_lambda_B` 4-5.) Side effect: this restores the TPI loading above
  1 (2.54/2.14/1.74 at gamma=2/5/10, still declining), resolving the "loading
  <1" concern raised in a paper-direction hostile review the same day -- at
  roughly a third the magnitude of the stale "~7x" figure previously in
  `docs/FRAMING_HANDOFF.md` (now retired; see `docs/SPEC.md`'s theoretical
  framing section, updated accordingly). **Re-tuned again same day to
  `psi_lambda_B=1.1793`** after resolving `recovery_rate` below (`EL_price`
  fell, pulling the spread response with it; re-verified 150.02bp).
- ~~**EL-1 (`recovery_rate` placeholder)**~~ **DONE (2026-07-22).** Set to
  `0.30`, an NPV-recovery estimate for the actual March 2012 Greek PSI
  (Zettelmeyer-Trebesch-Gulati "Autopsy": 59-65% NPV loss; contemporaneous bank
  estimates: 73-78%; 0.30 recovery / 70% haircut sits centrally). `EL_price`
  fell from 0.1025 to 0.0717 as a result, which required the `psi_lambda_B`
  re-tune noted above. TPI loading rose to 3.59/3.03/2.47 at gamma=2/5/10 (a
  smaller `EL_price` denominator makes the same premium look like a bigger
  multiple of fair compensation).
- ~~**PT-1 (pass-through validation)**~~ **DONE (2026-07-22).** Checked the
  model's ≈-4.5%/100bp bank-net-worth-to-spread pass-through against
  Acharya-Drechsler-Schnabl (2014 JF, Table 6): their bank-equity-return-on-
  sovereign-CDS coefficient (-0.096, post-bailout, full controls) implies
  -1.8% to -8.6%/100bp depending on the baseline CDS level assumed for the
  log-to-level conversion. The model's figure sits within this range under
  every reasonable baseline -- same order of magnitude, not a stray number.
- **`rk_F` depreciation-calibration miss:** `code/depreciation_calibration.py`
  targets `rk_F=0.01` by setting `delta_F` from the *pre-re-solve* `K_F`/`Y_F`
  (a one-shot, not iterated-to-convergence, calibration), then re-solves SS with
  the new `delta_F`. Post-EBA (thin `n_inter_F`, `omega_K_F=0.019`), this one-shot
  approximation misses: actual `rk_F` comes out at `0.0133`, not `0.0100`
  (`rk_D` hits its target exactly; only F misses). Likely pre-existing and
  independent of the C-1 fix, but not verified against a pre-EBA baseline.
- **Small `Y_D[0]` sign anomaly on the default shock:** `Y_D[0]` comes out
  *positive* (`+0.0124%` at baseline calibration, `+6.17e-4%`/`+1.24e-4%` in the
  EBA runs) even though `n_inter_D[0]` is correctly negative. Economically
  plausible as a portfolio-substitution effect (bonds `Delta=0.4` are now good
  collateral at baseline and become relatively worse on impact, so capital demand
  can rise enough to offset the direct hit to `Y_D` at t=0) -- this tension
  between "substitution" and "deleveraging" channels was already flagged in an
  earlier diagnostic (`diagnostics/` "substitution vs deleveraging" experiment,
  commit `1e82e22`). Small in magnitude (two orders below `n_inter_D`), but worth
  checking before reporting `Y_D` impact signs in the paper.

## Reproduce

```bash
/opt/anaconda3/envs/ssj/bin/python code/eba_calibration.py     # moments
/opt/anaconda3/envs/ssj/bin/python code/main.py                # full pipeline
```
