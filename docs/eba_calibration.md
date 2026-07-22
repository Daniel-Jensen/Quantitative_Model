# EBA-anchored calibration

Source: EBA **2011 EU-wide stress test** disclosure (`data/DATA_DISCLOSURE.CSV`),
base **31 Dec 2010 actual**; nominal GDP from Eurostat `nama_10_gdp` (2010).
Country map: **D = Greece** (banks GR030–GR035, issuer GR), **F = Germany**
(banks DE0xx, issuer DE). Moments produced by `code/eba_calibration.py` →
`data/eba_moments.json`. Code map validated: `DE017 = Deutsche Bank`
(2010 RWA €346.6 bn, CT1 €30.4 bn, 8.76% — matches published).

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
and `code/main.py`'s live EBA pipeline both confirm the fix; see
`audit_artifacts/run_audit.py`'s `check()` output for the regression harness.

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
