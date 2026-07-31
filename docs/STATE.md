# Project State

**Branch:** `main` | **Date:** 2026-07-22 | **Status:** EBA-anchored calibration + C-1 structural fix (see below); post-forensic-audit baseline still applies underneath

## CURRENT CALIBRATION (2026-07-30) — supersedes every table below

**The calibration was reverted to its pre-EBA values.** Everything below this section
describes the EBA-anchored calibration and is now **historical**. Structural fixes
(C-1, W-1/W-2/W-3, T-2, A-2, TPI-1, `omega_K`, capital-key conduit) are all retained —
only parameter values moved.

| Parameter | Value | Note |
|---|---|---|
| `psi_lambda_B_D/F` | 3.0 | pre-EBA default; gives 187.2bp ann vs 150bp target (~25% over) |
| `n_inter_D/F` | 3.0 (=0.75×4) | bank net worth, 0.75 of annual GDP |
| `omega_K_D/F` | 1.0 | capital fund empty (`div_fund=0`) — pre-EBA balance sheet exactly |
| `phi_lamb_D/F` | 0.15 | minimum stabilising value under the par rule |
| `mv_rule_D/F` | 0 (par) | see Finding F-1 below — **1 is unusable at `phi_lamb=0.15`** |
| `recovery_rate_D/F` | 0.30 | EL-1 **retained**; live only via `EL_price` while `writeoff_enabled=0` |
| `delta_b_D/F` | 0.10 | 2.5yr; empirical 7yr/6.5yr still not ported (needs `mv_rule=1`) |
| portfolio targets | 0.25/0.15/0.15/0.25 | `steady_state.py`; EBA's 2.39/0.018/0.069/2.76 retired |

**Verified end-to-end** (`code/main.py`, exit 0): `n_inter_D[0]=-3.0009%`,
`Y_D[0]=-0.0261%`; `b_gov_D[499]=-1.3e-5`; `rho_b=0.8451`; IC-δ exact;
`max|ca_res_D|=6.3e-8`, `max|goods_mkt_F|=1.1e-9`; TPI monotone in γ.

**Finding F-1 sharpened.** The `phi_lamb∈[0.15,0.18]` zone identified under `mv_rule=1`
is a **hard break**, not the "narrow, mild zone" described below. Measured directly:
`mv_rule=1` + `phi_lamb=0.15` → `n_inter_D[0]=-1554%`, `Y_D[0]=+0.170%` (perverse sign),
`b_gov_D[499]=1.6e-2`. `mv_rule=1` needs `phi_lamb=0.60` (→ `-5.89%` / `-0.024%` / `0.0`).

**Default-loading split.** `EL_price_D=0.0717` vs `psi_spread_D=0.8385` → fundamental
expected loss is **10.9%** of the default loading, GK collateral friction **89%**.

**Units.** `spread_rb` is a *quarterly* rate deviation; annualise ×4×1e4 for comparison
with the 150bp target.

**Policy-regime feature (2026-07-31): runs end-to-end at this calibration.** Stage A,
Stage B-lite and the unit tests all pass (exit 0) after three fixes: `PSILAM_MAIN` now
reads the live calibration (it was hardcoded while doubling as the cache key),
`PSILAM_BREAKDOWN=4.0` replaces the EBA-specific `<1.5` guard, and
`gamma_for_compression`'s scan narrowed 60->25. Headline: `A_cb=-2.406e-2` (backstop
compresses, SA-1 absent), `gamma_aggressive=5.0813` / `gamma_medium=1.5730`, peak spread
187.2/140.4/93.6bp passive/medium/aggressive; Stage A A6 amplifier invariance holds
cleanly (9.33/8.54/7.14bp at `psi_lambda_B=0`).

**A6 at the lottery stage — fixed 2026-07-31, and the invariance is real.** The old
check ranked the *full-sample* peak, which is the common pre-`k` spread (no branch has
acted yet), so it was identical across branches by construction and the strict
inequality compared floating-point noise: it read YES at `psi_lambda_B=0` on a
`+2e-14`bp gap and NO at the calibrated `psi_lambda_B` on a `-5e-13`bp gap — same
expression, opposite verdicts. Now ranked on the **post-revelation window `t>=k`** and
checked at **both** amplifier settings, with a `1e-3`bp separation margin so noise can
never produce a verdict:

| `psi_lambda_B` | aggressive | medium | passive | ordered |
|---|---|---|---|---|
| as calibrated (3.0) | 76.98 | 117.87 | 160.19 | YES |
| 0 (fundamental floor) | 4.91 | 5.89 | 6.43 | YES |

Separation at `psi_lambda_B=0` is 0.548bp, ~550x the margin. **A6 amplifier invariance
holds in the lottery as well as in deterministic Stage A** — the previous "false pass"
verdict was a measurement-window error, not a failure of the economics.

> **Note:** `audit_artifacts/` was removed 2026-07-30 — the harness carried its own
> hardcoded copy of the calibration instead of importing `get_calibration()`, so it
> silently tested a different model than `code/main.py`. Every `audit_artifacts/*`
> path cited below is **historical provenance** (what was run at the time), not a
> live command. Scripts recoverable from git history at `0c99013`.

## EBA calibration + C-1 fix (2026-07-22) -- current state supersedes the sections below

Everything from here to "What is complete (post-audit)" describes the state as of
the 2026-06-22 forensic audit. Since then (2026-07-13 to 2026-07-22), the
calibration was re-anchored to EBA 2011 stress-test data (see `docs/eba_calibration.md`
for the full derivation) and a structural bug (C-1) was found and fixed. Summary:

- **Calibration is now EBA-anchored, not the placeholder values below.**
  `n_inter_D/F`, `phi_bD_D_ss` etc., bank-held debt stocks, and bank concentration
  ratios all come from the EBA 2011 disclosure (`data/DATA_DISCLOSURE.CSV`), not
  the round-number placeholders in the table further down. See `docs/eba_calibration.md`
  for the parameter -> moment map.
- **A passive capital-intermediation fund (`omega_K`) was added**: banks hold only
  `omega_K` of the physical capital stock (`omega_K_D=0.0601`, `omega_K_F=0.0190`);
  the rest sits in a deposit-funded passive fund that rebates its spread to
  households (`capital_fund_D/F`, `div_fund_D/F`). Needed because EBA's thin bank
  net worth (`n_inter_D=0.408`) can't plausibly intermediate the entire capital
  stock. Verified Walras-neutral and a no-op at `omega_K=1`.
- **C-1 (`Delta_cross > 1`, degenerate IC) is fixed at its root, not calibrated
  around.** The single-asset `lambda_gk` formula in `steady_auxilliary_D/F`
  implicitly assumed full (`Delta=1`) bond divertability; substituted into the
  back-solve, this forced `Delta_cross > 1` whenever `Delta_own > 0.5` -- i.e.
  essentially always at realistic (EBA) concentration. `steady_auxilliary_D/F`
  now solve `lambda_gk` from the multi-asset IC directly (closed form, since
  `value` is linear in `lambda_gk`), taking `Delta_bD_D=0.2`/`Delta_bF_D=0.4` (D)
  and `Delta_bF_F=0.2`/`Delta_bD_F=0.4` (F) as genuine hardcoded calibration
  inputs instead of back-solved outputs. Verified to recover those exact values
  (to floating-point) via `code/ic_delta_calibration.py`'s consistency check and
  independently via `audit_artifacts/run_audit.py`. See `docs/eba_calibration.md`
  "Why C-1 forces Delta->1" for the full derivation.
- **The "explosive EBA doom loop" finding in `docs/eba_calibration.md` was
  superseded by the C-1 fix, not just improved.** With C-1 fixed and *no other
  calibration change* (`psi_lambda_B=0.31`, `mv_rule=1`, `phi_lamb=0.60` all
  unchanged from the pre-fix explosive runs), `b_gov_D[499]` is `~1e-5` on both
  shocks (was 3.4-79.5 across the variants previously tested) -- near-stationary,
  not merely damped. The explosive root was the C-1 degeneracy amplifying itself
  through the leverage/IC loop, not an intrinsic feature of the 2.39x EBA
  concentration ratio as previously concluded.
- **A second, independent bug was found and fixed alongside C-1**: `code/tpi.py`
  builds its own copy of the dynamic model and had not been updated for the
  `omega_K` capital-fund commits -- it was missing `capital_fund_D/F`, so
  `div_fund_D/F` was silently frozen at its SS level instead of responding to
  shocks in the TPI experiment. This produced a real (not machine-precision)
  accounting leak on every `gamma` (`max|ca_res_D|` up to 6.2e-2, comparable in
  size to the reported welfare gains). Fixed by adding the missing blocks; the
  `G_tpi[cb=0]` vs. baseline-Jacobian sanity check now matches to `0.00e+00`
  (was `2.54e-3`). **Any TPI welfare number from before 2026-07-22 is superseded.**
- **`audit_artifacts/run_audit.py` had drifted out of sync** with the shared
  equation files across three separate additions (`omega_K`, `mv_rule`,
  `EL_price` from the macro-pru-fix) and could not run at all without patching.
  Fixed: its calibration fixture now has the missing (neutral-default) params,
  its duplicated `_apply_ss_anchors` now imports the one true definition from
  `steady_state.py`, and it gained real regression assertions (Walras thresholds
  from this file's own acceptance table, the default-shock sign check from
  "Typical iteration" step 4, and a C-1 consistency check) that fail loudly with
  a collected summary instead of silently passing or crashing on unrelated
  `KeyError`s.
- **Two smaller items are flagged, not yet resolved** (see `docs/eba_calibration.md`
  "Remaining open items"): a `rk_F` depreciation-calibration miss (0.0133 vs.
  target 0.01, likely a pre-existing one-shot-calibration imprecision) and a small
  (two orders of magnitude below bank net worth) positive `Y_D[0]` on the default
  shock, plausibly a portfolio-substitution effect rather than a bug.
- **`psi_lambda_B` recalibrated 0.31 → 1.1284 (2026-07-22, later same day).**
  0.31 was never a calibration target -- it was chosen only to dodge the (now
  fixed) C-1 degeneracy and undershot the paper's external anchor (2010 GR-DE
  spread ~150bp on a 1pp default shock) by more than 3x (44bp). Re-ran the
  moment-matching exercise on today's model (`audit_artifacts/psilam_moment_sweep_postC1.py`)
  and found the spread-vs-`psi_lambda_B` response is smooth/monotonic only up
  to about 1.5-2.0, then turns wildly non-monotonic (219bp at 2.0, 97bp at 2.6,
  853bp at 2.8, 353bp at 3.0, 11478bp at 5.0) -- a linear-approximation
  breakdown, not real moments. **Both the pre-fix literature value (2.8, per
  `docs/SPEC.md`'s theoretical framing) and the original round-number default
  (3.0) now sit inside that broken region on this model and must not be
  restored as-is.** `psi_lambda_B=1.1284` verified directly to hit 151.3bp
  (smooth neighbourhood: 147bp at 1.10, 154bp at 1.15) and is now the
  committed value. Side effect: this also resolves a "loading <1" concern
  raised in a same-day paper-direction hostile review (the model's TPI
  compensation appeared to have inverted below the actuarially-fair level at
  the pre-recalibration `psi_lambda_B=0.31`) -- loading is back above 1
  (2.54/2.14/1.74 at γ=2/5/10) and still declining, restoring the paper's
  central over-compensation claim, at roughly a third the magnitude of the
  stale "~7x" figure in `docs/SPEC.md`'s theoretical framing, which has been
  updated accordingly.

## Policy regimes (exogenous backstop aggressiveness) — added 2026-07-23

Three exogenous ECB backstop regimes over the TPI feedback coefficient γ on the
endogenous D–F yield spread (`TPI_t = γ·(spread_t − spread_ss)`, dormant at SS),
plus a Stage-B lottery in which banks price an ex-ante-unknown CB type revealed at
date k. All post-Jacobian numpy on cached `G_tpi` response matrices — production
`main.py` untouched. Code/outputs in `diagnostics/regimes/` (`run_regimes.py` =
Stage A, `uncertain_regime.py` = Stage B, `lottery_math.py`+tests = the math core,
`regime_model.py` = the cache layer; see `regimes_log.md`). Beliefs disciplined by
the 3-state Markov-switching intervention-regime estimation in `Empirics/` (via
`beliefs_from_empirics.py`).

**The backstop compresses spreads on this model, as intended.** A unit CB purchase
gives `d(spread_rb)/d(cb_buy) = −1.95e-2` (compression): the ECB **capital-key
conduit** (`kappa_cb_F=0.929`) funds a D-bond purchase 92.9% through F's treasury,
so D banks shed bonds (`d(b_D_D)=−0.72`) and the periphery–core spread narrows.
Stage A (spec §7 spread-compression targeting): aggressive/medium γ = 17.05/5.72
hit 78/117 bp peak spread vs 156 bp passive, cushioning bank net worth; the A6
amplifier-invariance ranking **survives** at `psi_lambda_B=0`. Stage B: the impact
spread rises with the belief weight on a passive CB (fear of no backstop is priced
before any policy acts), with a positive regime-uncertainty premium (+52/+23 bp on
the aggressive/medium branches); the §10.3 structural identities hold exactly.

**Provenance caveat (important).** An earlier build of this feature on the retired
`ms-regime` branch reported the OPPOSITE — that CB purchases *widen* the impact
spread ("Finding SA-1"), forcing an output-protection workaround and an inverted
Stage-B sign. That result was an artifact of `ms-regime`'s **superseded** model
(single-country conduit, par-value fiscal rule, `psi_lambda_B=2.8` — now known to
sit in this model's linear-approximation-breakdown region). It does **not** hold on
main; the capital-key conduit + market-value rule resolve it. Any SA-1/SB-1
"spread-widening" or "output-protection" statement from `ms-regime` is void here.

## Historical state (as of 2026-06-22 forensic audit) -- see above for what changed since

## Current status

Six structural/accounting bugs were found and fixed in the 2026-06-11 forensic audit (W-1, W-2, W-3, T-2, A-2, TPI-1). See `docs/audit.md` for the full ranked finding list and `docs/verification_report.md` for verified fix status. All six are **merged into `main`** via PR #27 (`AB-audit`, merged 2026-06-11). PR #26 (`audit` branch) was closed as superseded. `main` was subsequently reorganised into modular Python files via PR #28 (merged 2026-06-22).

Core equations: `code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py`. `code/main.py` is the production entry point (added by PR #28; orchestrates `code/full_model.py`, `code/steady_state.py`, `code/tpi.py`, …). The legacy notebook `code/model_v12.ipynb` has been removed. TPI/IRF figures are produced by `code/tpi_plots.py` and `code/irf_plots.py` (no committed `plots/` dir).

`main` now contains all six fixes (PR #27) plus the modular-file reorganisation (PR #28). Use `main` — or a feature branch off it — for all new work. The `audit` / `AB-audit` branches are historical and should not be reused.

## What is complete (post-audit)

- Household deposit choice and GHH preferences for D (Greece) and F (Germany) households.
- Bank steady-state and intermediation blocks: capital, bond returns, fees, GK Bellman (P1) and IC constraint (P3/lambda_gk).
- Production, capital adjustment, and capital producer profit — W-1 fixed: production uses `Y=F(K_t)`, capital producer receives `mpk·(K−K(−1))` so all capital income is allocated.
- Deposit return predetermined correctly: `Rgross = (1+rdep(−1))·P(−1)/P` — T-2 fix. Funding legs in `bank_return_*` and FOCs in `intermediation_P1_*`/`divert_*` use ex-ante deposit rate.
- F-bank bond returns converted to F-goods via `p(−1)/p` in `bank_return_F` — W-2 fix (the dominant leak; was causing ~2% of F GDP goods_mkt_F residual on a 1% TFP shock).
- Cross-border bond FOC in F-bank uses `p/p(+1)` for expected return — W-3 fix (optimality condition; does not affect Walras but required for internally consistent portfolio choice).
- Smart steady-state blocks: `m = n·(1−(1−f)·(1+rn))` without spurious `+Phi+T` — A-2 fix (required for any `chi1≠0` calibration, e.g. bank-cal's chi1=0.5).
- Global goods market, external account, bond clearing, and portfolio adjustment cost blocks.
- Domestic and foreign bond pricing, yields, spreads, and Hatchondo-Martinez geometric-decay perpetuity default mechanics.
- TPI extension: CB budget closed via `budget_residual_D_tpi` remittance — TPI-1 fix. Before the fix, unbacked CB flows inflated welfare gains by ~40% at γ=10.
- **ECB balance sheet as a capital-key conduit** (ecb-balance-sheet branch): the TPI-1 wiring remitted the *entire* CB cash flow to the Greek treasury (self-financing; no German side). Now `cb_flow_D = (1+rb_actual_D)·q_b(-1)·cb_buy(-1) − q_b·cb_buy` splits `kappa_cb_F=0.929` to the F treasury (`budget_residual_F_tpi`, converted `/p`) and the rest to D; the F share of the CB book enters the external account like `b_D_F` (`external_account_D_tpi`), so the CB hole cancels identically and `ca_res_D`/`goods_mkt_F` stay at baseline leak levels at every γ. Carry/credit legs and off-path expected loss are post-processed in `cb_pnl` (`code/tpi.py`), discounted at `beta_F`; `audit_artifacts/tpi_test.py` imports the production blocks and asserts `CB_CONDUIT_TEST: PASS` at ≤1e−7.
- EBA bilateral sovereign exposures in calibration cell: b_D_D/asset=24.47%, b_F_F/asset=25.79%, b_F_D/asset=0.18%, b_D_F/asset=0.65%.

## Walras accounting (post-fix, verified)

| Residual | 1% TFP-D shock | 1pp default-D shock |
|----------|----------------|---------------------|
| goods_mkt_D (targeted) | ≤1e−16 | ≤1e−16 |
| goods_mkt_F (untargeted) | ≤8e−10 | ≤1e−9 |
| ca_res_D = CA−ΔNFA (untargeted) | ≤5.8e−8 | ≤3.5e−8 |
| deposit_mkt_D/F | ≤4e−15 | ≤4e−15 |

Pre-fix peaks for reference: goods_mkt_F 2.0e−2 (~2% of F GDP); ca_res_D 1.5e−4. All cross-country spillover and welfare results from the pre-fix model are first-order invalid and must be regenerated from `main`.

## IRF summary (historical, pre-EBA — phi_lamb=0.15 era; SUPERSEDED)

> **Superseded (2026-07-24).** The committed fiscal-feedback coefficient is
> `phi_lamb=0.60` (see the C-1-fix section at line ~37 and `code/calibration.py`),
> **not** the `0.15` this section's title assumed — that was a stale table value,
> now corrected. The impact numbers below also predate the current EBA /
> `psi_lambda_B` / `recovery_rate` calibration; the 2026-07-24 production run gives
> `n_inter_D[0]=−2.83%`, `Y_D[0]=+0.032%`, peak D–F spread `+0.392pp`,
> `b_gov_D[499]≈2e−6`.

**1pp default shock to D (ρ=0.8):**
- `n_inter_D[0] = −3.5%` (falls), `Y_D[0] = −2.5e−4` (falls) — both signs correct post-T-2-fix; were positive/perverse pre-fix.
- `n_inter_F[0] ≈ −0.33%` — contagion small, sign correct.
- Spread widens on impact; doom loop is live with correct sign.

**TPI (ECB capital-key conduit, ecb-balance-sheet branch, 2026-07-16):**
- ΔW_D = +1.38/+2.75/+4.13, ΔW_F = −1.36/−2.70/−4.03 (% quarterly SS consumption, 100q) at γ=2/5/10. Near-transfer, slightly positive-sum (+0.10 at γ=10). Spread now compresses with γ: peak 0.409pp (γ=0) → 0.161pp (γ=10).
- ECB P&L (PV at β_F over 100q, % of quarterly SS Y_D, γ=10): peak exposure 1.57%, purchases PV 4.20%, expected-loss leg 0.0070%, default-premium leg 0.0227%, MTM leg +0.0253%; SS-carry 0 (SS yields equalised).
- **Loading (premium PV / EL PV) declines monotonically in γ: 4.86 → 4.06 → 3.22** — the self-extinguishing premium (`docs/SPEC.md` "Theoretical framing," Live Claim 5) confirmed in-model at this (now superseded) calibration; `loading_arr` over `gammas_fine` in `tpi_results` is the key-figure schedule. **Superseded 2026-07-22**: current committed calibration gives 2.54/2.14/1.74 at γ=2/5/10 (see the EBA/C-1-fix section at the top of this file).
- Germany (κ=0.929): bears EL PV 0.0065% Y_D, receives premium PV 0.0210% at γ=10 (memo at full EA key 26.1%: EL 0.0018%).
- Conduit residuals at every γ: `ca_res_D` ≤ 7.4e−8, `goods_mkt_F` ≤ 6.2e−10. The pre-conduit "zero-sum transfer, no spread compression" numbers (ΔW_D +1.88/ΔW_F −1.90) were from the notebook-era model and are superseded.

## Calibration summary (current, main)

| Parameter | Value | Source / note |
|-----------|-------|---------------|
| `phi_lamb_D/F` | 0.60 | ~Bohn (1998) fiscal-feedback magnitude. Hand-set for stationarity under EBA doom-loop amplification (not moment-matched); enters only `tax_rule_D/F`, and under `mv_rule=1` reacts to the market-value debt gap (so it couples to `q_b`). Well above the F-1 near-unit-root zone `[0.15,0.18]`; stability confirmed (2026-07-24 run: `b_gov_D[499]≈2e−6`, `ρ_b=0.373`). Governs the debt/fiscal mode only — **not** the ~25q financial-accelerator ring (that mode is set by the PAC `psi_bF_D/psi_bD_F`). **Corrected 2026-07-24 from a stale `0.15` entry.** |
| `def_scale_D` | 0.25 | Strong amplification. Exceeds GR 2011 crisis peak (0.12–0.23 from spread-debt slope calibration). |
| `delta_b_D/F` | 0.10 | 2.5yr avg maturity. Empirically too short; bank-cal has 0.036/0.038 matching GR/DE 2011 ~7yr/6.5yr. |
| `theta_D/F` | 4.0 | GK leverage; conservative vs 2011 historical 10–25×. |
| `psi_lambda_B_D/F` | 3.0 | State-dependent divertability; primary amplification dial; no direct empirical counterpart. |
| `f_D/F` | 0.12 | Bank exit rate; bank-cal has 0.03 (standard GK range). |
| `kappa_cb_F` | 0.929 | F share of CB conduit cash flows: two-country renormalised euro-area capital key (BuBa 26.1 / BoG 2.0). TPI layer only; SS-neutral. Memo reporting also at full EA key 26.1%. |
| `Delta_cross` | 1.4545 | Back-solved (`_ic_delta`, ratio=2.0); degenerate >1. See C-1. |
| `recovery_rate_D` | 0.00 | No realized losses; inert while writeoff_enabled=0. |
| `writeoff_enabled_D/F` | 0.0 | Default produces no balance-sheet losses. See S-1. |
| `chi1_D/F` | 0.0 | Intermediation adjustment cost off. A-2 fix makes chi1≠0 safe. |
| `frisch` | 0.5 | Frisch elasticity. |
| nDep_D/F | 500/500 | Household deposit grid points. |
| income rho_z/sigma/nZ | 0.90/0.30/15 | D and F income process (Markov approximation). |

## Open issues

| ID | Description | Status |
|----|-------------|--------|
| C-1 | ~~`Delta_cross=1.45 > 1`: divertable fraction exceeds 1, making the multi-asset IC constraint degenerate.~~ | **RESOLVED (2026-07-22).** `lambda_gk` now solved from the multi-asset IC directly; `Delta_bD_D/F=0.2/0.4` are genuine hardcoded inputs, verified to bind exactly. See "EBA calibration + C-1 fix" above and `docs/eba_calibration.md`. |
| S-1 | `writeoff_enabled=0`: default shock produces zero realized bank losses. `recovery_rate` and `zeta_writeoff` are set but inert. Model is currently a pure risk-premium loop, not a balance-sheet doom loop. | Author decision. Resolution: set `writeoff_enabled=1` with `recovery=0.40, zeta=1.0` (GR 2012 ~50% haircut → ~0.4–0.5 recovery). |
| RK-1 | `rk_F` depreciation calibration misses its target (0.0133 vs. 0.0100) post-EBA; the one-shot `delta_F` calibration in `depreciation_calibration.py` isn't iterated to convergence against the re-solved SS. `rk_D` hits its target exactly; only F misses. | Open (2026-07-22). Not yet investigated against a pre-EBA baseline. |
| Y-1 | Small positive `Y_D[0]` on the default shock (`+0.01%` to `+0.03%` across calibrations) even though `n_inter_D[0]` is correctly negative. Plausibly a portfolio-substitution effect (bonds are now good collateral at `Delta=0.4` and become relatively worse on impact) partially offsetting the direct hit — a tension already flagged in an earlier "substitution vs deleveraging" diagnostic. | Open (2026-07-22). Small magnitude; worth checking before reporting `Y_D` impact signs. |
| EL-1 | ~~`recovery_rate_D/F=0.00` is a placeholder~~ **RESOLVED (2026-07-22).** Set to `0.30` (NPV-recovery framing for the actual March 2012 Greek PSI): Zettelmeyer, Trebesch & Gulati ("The Greek Debt Restructuring: An Autopsy", PIIE WP13-8) find actual investor NPV losses of 59-65%, considerably below the ~75% commonly quoted; contemporary bank estimates (Credit Suisse, Morgan Stanley) put NPV haircuts at 73-78%. 0.30 recovery (70% haircut) sits centrally in this range. NPV framing chosen over the face-value framing (53.5% face cut -> ~46.5% recovery) because `EL_price` multiplies a payoff already priced at market terms (`q_b`), matching the NPV concept. `EL_price` fell from 0.1025 to 0.0717 as a result; `psi_lambda_B` re-tuned from 1.1284 to 1.1793 to restore the 150bp target exactly (verified: 150.02bp). Full pipeline reverified clean at the new calibration. | **Resolved.** Side effect: TPI loading rose to 3.59/3.03/2.47 at γ=2/5/10 (was 2.54/2.14/1.74) — a smaller `EL_price` denominator makes the same premium income look like a larger multiple of fair compensation. Update `docs/SPEC.md`'s theoretical framing numbers accordingly. |
| PT-1 | ~~Bank-net-worth-to-spread pass-through unvalidated~~ **CHECKED (2026-07-22).** Model gives ≈−4.5%/100bp at the calibrated point. Acharya-Drechsler-Schnabl (2014 JF, Table 6 col. 9, post-bailout, full controls) report bank equity return regressed on Δlog(sovereign CDS): coefficient −0.096 (s.e. 0.017, 1% sig.). Converting their log-CDS elasticity to a per-100bp-of-level moment requires a baseline CDS level (log vs. level don't translate 1:1): using their post-bailout sample mean (108.5bp) gives ≈−6.3%/100bp; using their reported typical weekly move (11.3% log-change at that mean, linearly rescaled to 100bp) gives ≈−8.6%/100bp; using a higher, more Greece-specific crisis baseline (~500bp) gives ≈−1.8%/100bp. **The model's −4.5%/100bp sits within this literature-implied range (−1.8% to −8.6%/100bp) under every reasonable baseline choice** — same order of magnitude, not a stray number. Altavilla-Pagano-Simonelli (2017, Rev. Finance) corroborates the *mechanism* (bank sovereign exposure amplifies transmission of sovereign stress to bank outcomes) but its headline results are on lending, not equity/net-worth, so it wasn't used for a second point estimate. | **Resolved as "within plausible range," not exact-match** — appropriate given the model has no free parameter tuned to hit this moment; it emerges from the calibrated collateral/leverage block. Re-check if `psi_lambda_B` or the bank-capital calibration changes again. |
| DIST-1 | Distributional resolution is Gini-only this iteration; per-decile deferred. `GINI_C` is specifically the wrong statistic for the Greek crisis — measured inequality barely moved (already highly unequal pre-crisis; worst-affected households dropped out of surveys). | Open. Any incidence claim in the paper must state its distributional resolution explicitly. |
| A5-1 | The TPI's German fiscal cost/benefit should be reported as three distinct objects — **exposure** (discounted purchases, what capital-key sharing acts on), **priced expected loss** (hand-computed off-path per `docs/SPEC.md`'s "structural constraints"), and **Greek fiscal saving** (`pd_D` differential vs. a passive counterfactual, the cleanest headline, no pricing assumption needed) — not conflated into a single "implicit transfer." Discount rate should be German risk-free for a cross-border object (currently `beta_F`, the creditor side — check this is deliberate and justified in the paper, not just a default). | Open; check `code/tpi.py`'s `cb_pnl` reports all three cleanly before writing this section of the paper. |
| X-1 | Dead-code imports in notebook cell 7: blocks no longer in the model remain in the import list. | Minor cleanup; no numerical effect. |

## Next priorities

**Superseded by the 2026-07-22 EBA calibration + C-1 fix (see top of file):**
EBA bilateral exposures are ported and live (item 1's exposures/Delta=0.2/0.4
portion is done); `delta_b`/`f` are **not yet** ported to the bank-cal empirical
values (still `0.10`/`0.12`) — that part of item 1 remains open. Figures were
regenerated 2026-07-22 under the fixed, EBA-anchored calibration (item 4 done for
now, but re-do again after any further calibration change).

1. **Port remaining bank-cal calibration values**: `delta_b=0.036/0.038` (empirical
   7yr/6.5yr duration), `f=0.03` (bank exit rate). Given C-1 is now fixed and the
   baseline doom loop is stable without it, re-test Finding F-1's explosive-root
   concern (below) under the *fixed* model before assuming it still applies —
   it was diagnosed on the pre-C-1-fix, non-EBA calibration.
2. **Decide S-1**: set `writeoff_enabled=1` to give default realized losses, or
   keep pure risk-premium framing and state it explicitly in the paper.
3. ~~Re-explore `psi_lambda_B`~~ **DONE (2026-07-22).** Recalibrated to 1.1284
   (150bp spread target), replacing the never-validated 0.31. See the entry
   above this table for the full rationale and the discovered non-monotonic
   breakdown region above `psi_lambda_B~2`.
4. **Investigate RK-1 and Y-1** (open issues table above) before reporting `rk_F`
   or `Y_D` impact-sign results in the paper.
5. **Re-generate all figures** from `main` after any of the above — the current
   `outputs/` figures (2026-07-22) reflect the EBA calibration + C-1 fix + TPI
   conduit-leak fix, but not any further calibration change.

## Finding F-1: market-value fiscal rule (duration ↔ fiscal-stability tension)

**Historical diagnosis below predates the C-1 fix** (kept verbatim; the doom loop's
amplification was resting on a degenerate collateral value at the time). **Re-tested
2026-07-22** on the now-fixed model, using today's actual committed calibration
(EBA-anchored, `psi_lambda_B=0.31`, `def_scale=0.25`, `writeoff=0`) with only
`delta_b` (duration), `phi_lamb`, and `mv_rule` varied. This went through three
rounds before landing on a trustworthy answer -- worth keeping all three, since the
methodological correction is as important as the result:

**Round 1 (coarse grid, energy-ratio proxy, `audit_artifacts/philamb_sweep_postC1.py`).**
`mv_rule=0` (par value): explosive at every `phi_lamb` in `{0.05,...,0.60}` tested,
but mildly (modulus 1.002-1.005, vs. 1.005-1.015 pre-fix; `bgov[499]` ~1e-6, vs.
9.5-79.5 pre-fix). `mv_rule=1` (market value): STABLE at 0.05, EXPLOSIVE at 0.10 and
0.15, STABLE again at 0.25-0.60 -- a non-monotonic "explosive island" sitting exactly
where the pre-fix finding's stable plateau (`phi_lamb` in `[0.07,0.12]`) was reported.

**Round 2 (finer grid, naive Prony/eigenvalue estimator, `philamb_sweep_postC1_fine.py`).**
The two-window energy-ratio proxy used in round 1 is known to be foolable by
oscillatory near-unit-circle modes, so a proper linear-prediction (Prony) estimate
was built: fit an order-`p` AR recursion to the IRF tail (t=150-500), take the
eigenvalues of the resulting companion matrix. Validated first against synthetic
signals (recovered known poles to ~1e-15) and an adversarial two-mode case (correctly
resolved a true dominant pole an energy-ratio proxy would have underestimated by 7x)
-- then applied at `p=8` (with `p=4` as a cross-check) to the real IRFs. **This gave
an implausible modulus of ~3.2** at every grid point, flatly contradicting the
observed `bgov[499]` levels (~1e-6) a true modulus of 3.2 could never produce over
500 periods.

**Root cause (`audit_artifacts/philamb_order_selection_check.py`).** R^2 of the AR
fit already saturates to machine precision (`1-1e-9`) at order `p=2-3`; every order
beyond that fits pure numerical noise into spurious high-modulus "ghost" poles that
contribute negligibly to the actual signal but have huge companion-matrix
eigenvalues -- textbook AR/Prony overfitting once model order exceeds what the data
supports. `p=4` (round 2's "conservative" cross-check) was *already* past this onset
for at least one grid point -- it wasn't a safe fallback either.

**Round 3 (corrected: automatic order selection, `philamb_sweep_postC1_fine_v2.py`).**
Refit every cell at `p=1..4`, taking the modulus at the *smallest* order whose R^2
clears `0.99999` (falling back to `p=3`, flagged, if none does -- never triggered).
Result, `mv_rule=1`, fine grid over `phi_lamb` in `{0.05, 0.075, ..., 0.25}`:

```
phi_lamb:  0.05    0.075   0.10    0.125   0.15    0.175   0.20    0.225   0.25
trusted:   0.997   0.982   0.998   0.996   0.999   1.003   0.990   0.990   0.990
verdict:   STABLE  STABLE  STABLE  STABLE  ~1.000  mild-X  STABLE  STABLE  STABLE
```

**Verdict: the "explosive island" was mostly an artifact of the energy-ratio
proxy.** The properly order-selected estimate is comfortably stable (modulus
0.98-0.997) across nearly the *entire* `phi_lamb ∈ [0.05, 0.25]` range under
`mv_rule=1` -- including 0.10 and 0.125, which round 1's proxy called explosive.
There is a real (not noise-floor -- this system's Jacobian conditioning implies a
~1e-10 floor, far below the ~0.003 excess seen here), but narrow and mild, near-unit-
root zone right around `phi_lamb ≈ 0.15-0.175` (modulus 0.999-1.003). That is much
smaller and narrower than either the round-1 or round-2 estimates suggested, and it
sits at a different location than the original pre-C-1-fix "stable plateau"
(`[0.07,0.12]`) claimed, but "avoid `phi_lamb` in roughly `[0.15,0.18]`, most other
values in `[0.05,0.25]` are fine" is a defensible, much more precise statement than
either "there's a stable plateau at 0.10" (pre-fix finding) or "there's an explosive
island at 0.10-0.15" (round-1 proxy). `mv_rule=0` (par value) was not re-swept at
finer resolution or order-selected precision -- round 1's mild-explosive-everywhere
reading for it should be treated with the same "proxy could be off" caveat pending
that redo. See `audit_artifacts/philamb_sweep_postC1_fine_v2_results.json` for the
full per-cell R^2/order/modulus breakdown, and `philamb_order_selection_check.py` for
the overfitting diagnostic.

*Historical diagnosis (pre-C-1-fix, kept for record):*

> The Bohn rule responds to the **par/face-value** lagged debt gap (`tax_rule_*`). With empirical long-duration bonds (`delta_b=0.036/0.038`, 7yr/6.5yr) this is **explosive at every `phi_lamb ∈ [0.02, 0.50]`** under both the balance-sheet (write-off ON) and risk-premium channels — the debt dynamics collapse to a near-unit-root ~250-quarter cycle (dominant modulus ≈ 1.005–1.015) that fiscal feedback cannot damp. The original model is stable at `phi_lamb=0.15` only *because* its duration is short (`delta_b=0.10`). So empirical duration, a literature `phi_lamb`, and a stable live doom loop are jointly infeasible with the par-value rule. (Sweeps: `audit_artifacts/philamb_sweep*.py`.)
>
> **Resolution — market-value rule.** Reacting to the mark-to-market debt gap `q_b·b_gov(-1) − q_b_ss·b_gov_ss` (it sees the current spread) opens a stable plateau at `phi_lamb ∈ [0.07, 0.12]` (modulus down to 0.983 at 0.12; `phi_lamb≈0.10` robustly interior) with empirical duration and a live, correctly-signed doom loop — **but only in the risk-premium framing** (`psi_lambda_B=1.0`, `def_scale=0.10`, write-off OFF). With write-off ON it is a *false victory*: `|λ|<1` but the default shock is perverse (spread narrows, bank net worth and output rise). So adopting the market-value rule **forces the risk-premium framing — it couples to the S-1 decision** (keep write-off OFF).

**Status.** Implemented as a switchable option: `mv_rule_D/F` in calibration (`0`=par, default, behaviour unchanged; `1`=market value). `mv_gov_ss_D/F` is set from the solved SS in `build_and_solve`. Adopting `delta_b` (empirical duration) as the baseline calibration remains an author decision. Per the round-3 (order-selected) re-test above, `phi_lamb` in roughly `[0.05, 0.125]` or `[0.20, 0.25]` reads as genuinely stable under `mv_rule=1`; avoid the narrow `[0.15, 0.18]` zone. Do not reuse the pre-fix `[0.07, 0.12]` "stable plateau" language as if it were re-confirmed — it wasn't specifically re-tested at the old paper's exact framing (risk-premium, `psi_lambda_B=1.0`, `def_scale=0.10`), only at today's committed EBA amplification (`psi_lambda_B=0.31`, `def_scale=0.25`), which happens to also land in a stable region there but for different reasons. **Committed value (2026-07-24): `phi_lamb=0.60`** — above the entire fine-swept range `[0.05,0.25]`. Round-1 (coarse) read `mv_rule=1` as stable at `0.25–0.60`, and the 2026-07-24 production run confirms it (`b_gov_D[499]≈2e−6`; debt mode well-damped, `ρ_b=0.373`). The `[0.15,0.18]` caveat therefore does **not** bind the committed calibration — it matters only if `phi_lamb` is ever lowered toward the literature range. Note this is the *debt/fiscal* dominant mode; the ~25q financial-accelerator ring in the asset-price IRFs is a separate, faster complex pair governed by the PAC (`psi_bF_D/psi_bD_F`), not `phi_lamb` — see `audit_artifacts/pac_sweep.py`.

## Finding F-2: financial-accelerator ring (2026-07-24)

The baseline IRFs carry a damped oscillatory "ring", concentrated in the
asset-price / bank-net-worth variables (`n_inter_D`, `q_b_D`, `q_b_F`, `C_D`) and
**absent** from the real block (`Y_D`, `w_D` decay monotonically). Quantified at
the committed calibration with a Prony / companion-eigenvalue extractor
(`audit_artifacts/pac_sweep.py`):

- **Complex pair: |λ| = 0.954, period ≈ 25q (~6.25 yr), half-life ≈ 14.6q
  (~3.6 yr).** R²=1.00000; consistent across `n_inter_D`/`q_b_F`/`q_b_D`, Prony
  orders 6 and 8, and the empirical peak-spacing (~24q). Comfortably damped —
  **not** a stability concern (each swing falls to ~30% of the last; ~99% gone by
  q100).
- **It is the GK IC/leverage financial accelerator, not a portfolio-friction
  artifact.** A PAC sweep (`psi_bF_D=psi_bD_F` over 0.05→5.0, a 100× range) moves
  |λ| by <0.006 (0.954→0.960) and leaves the period fixed at ~25q; the dependence
  is mildly U-shaped with the committed `PAC=0.5` sitting at the modulus
  *minimum*. So the ring is intrinsic to the amplification block
  (`psi_lambda_B`/`def_scale`/`theta`/`f`), which is pinned to the 150bp spread
  target — **there is no free cosmetic fix.** Three modes, three owners:
  debt/fiscal → `phi_lamb` (fast-damped, `ρ_b=0.373`); cross-border position →
  PAC; the ~6yr financial cycle → the accelerator. Neither the fiscal nor the
  portfolio dial touches the third.
- **Takeaway:** treat the ring as a structurally-identified ~6yr financial cycle
  to *describe* in the paper's propagation discussion, not a bug to patch. Full
  per-cell R²/order/modulus table: `audit_artifacts/pac_sweep_results.json`.
