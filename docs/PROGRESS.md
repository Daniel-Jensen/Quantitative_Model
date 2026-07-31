# PROGRESS — Changelog

Reverse-chronological development log for the two-country MU-HANK sovereign-risk
model (Greece 2010–12 / ECB TPI). **Newest first.** Reverse-engineered from the
git history (135 commits, first commit 2026-05-14), `docs/STATE.md`,
`docs/audit.md`, `docs/eba_calibration.md`, and `docs/HANDOFF.md` — the detailed
derivations live in those docs; this file is the timeline.

Convention: `[hash]` is a commit on `main`; dates are commit dates. Doc-sync
policy: every commit touching model/code (`code/**`, `*.py`) gets an entry here.
(The pre-commit hook that was meant to enforce this is documented but not
installed in `.git/hooks` — treat it as convention.)

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
