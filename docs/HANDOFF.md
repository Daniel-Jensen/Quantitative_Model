# Handoff Notes

## Open problem: foreign banks do not retrench

On a 1pp default shock the F bank **increases** its Greek holdings (`b_D_F` rises;
peak around t≈13) instead of cutting them. That is contrary to the 2010–12 record
and to the closest published analogue — Bi, Foerster & Traum (FRBSF WP 2025-10),
whose Foreign intermediary reduces Italian holdings by ~2% of GDP while its own
economy still expands. It also undercuts the instrument's rationale: the ECB is
meant to be the risk-neutral buyer of last resort, with private core banks
retrenching and reinforcing the doom loop TPI addresses.

The country-size asymmetry (2026-08-07) shrank the impact response by an order of
magnitude — to a near-miss — but did not flip the sign, and from t≈4 the
mark-to-market term dominates outright. The FOC decomposition of `d b_D_F` off
`divert_portfolio_adj` isolates the four contributions (MTM, terms of trade,
`−d rdep_F`, risk premium) and is the right diagnostic for any candidate fix.

Ruled out by inspection: convex/threshold terms in `def_rate`. `def_rate_ss = 0`,
so anything quadratic has exactly zero first-order effect in the linearised
solution. Candidate routes, none yet chosen: a home-bias risk-pricing wedge in
`prem_DF`; a larger cross-border `psi_lambda_B` on `Delta_bD_eff` in
`intermediation_IC_F`; porting BFT's CES/Krenz bond aggregator (`sigma_b = -2`,
`gamma_b` from domestically-held shares); or a moral-suasion/forced-absorption
device on the D side, which is closest to BFT's actual mechanism.

**Note on targets.** The numbered "Live Claims", and the declining loading
schedule in particular, are **not** success criteria and must not be used to
evaluate model changes. Report them as outputs. The real gates are the calibration
moments (150bp peak GR–DE spread) and the correctness checks (Walras residuals,
GK well-posedness, impact signs).

## Where to start

- **FIRST: regenerate the downstream artefacts. They are stale.** The country-size
  asymmetry landed on 2026-08-07 (`size_F = 11.697`) with `psi_lambda_B` re-tuned
  **2.92 → 3.01**; see `docs/STATE.md` → *Country-size asymmetry*. `code/main.py`
  is verified against the new calibration, but E1–E4 and every paper artefact
  still reflect the old one — **and everything TPI-related predates the
  `rem_cb_F` conduit fix, so those numbers are wrong by an amount that grows with
  γ.** Run **in this order** — the ordering is load-bearing, the experiments never
  re-solve the model:
  ```
  /opt/anaconda3/envs/ssj/bin/python diagnostics/regimes/regime_model.py --force
  /opt/anaconda3/envs/ssj/bin/python experiments/run_all.py
  /opt/anaconda3/envs/ssj/bin/python experiments/e4_distribution.py
  /opt/anaconda3/envs/ssj/bin/python experiments/paper_outputs.py
  ```
  Affected: `docs/experiments_results.md`, `docs/paper_draft_results.md`, the
  eight tracked `experiments/paper/fig0*.png`.

- **PAPER EDIT REQUIRED — the constrained-seller number changed.** The default
  loading split is now **8.36% fundamental / 91.64% collateral friction**
  (`EL_price_D = 0.056134`, `psi_spread_D = 0.615358`), a ratio of 10.96:1. It
  was 3.4% / 96.6% (28.6:1). (`psi_spread` is exactly linear in `psi_lambda_B`
  at a fixed SS — verified 2026-08-07 when the 2.92 → 3.01 retune moved
  `psi_spread_F` 0.465088 → 0.479423, precisely the 3.01/2.92 ratio — so this
  rescaled from the 2.92 figures without a re-solve.) The claim survives in
  direction but "essentially
  all of the spread was a constrained-seller phenomenon" must become "roughly
  nine tenths of it". `experiments/paper_outputs.py`'s
  `fig04_spread_decomposition` caption is derived at run time and will pick this
  up automatically once regenerated — but the *prose in the paper* will not.

- **Sweeping `psi_lambda_B`: re-solve the pipeline, do not patch the SS.** The SS
  is genuinely `psi_lambda_B`-neutral (bit-identical `goods_mkt_D`, `K_D`,
  `beta_D` at every value), but patching `psi_lambda_B_D/F` + `psi_spread_D/F`
  onto a solved SS and re-solving only the Jacobian is still **wrong** — it
  predicted 150.33bp at `psi = 2.73` where the pipeline gives 139.60, because
  `intermediation_IC_D`'s `Delta_bD_eff` collateral channel ignores the patch.
  A full re-solve is ~2 minutes. The `rho_def` bisection was thrown away and
  redone over exactly this.

- **`rho_def` and `rho_Z` now live in `code/calibration.py`**, section *Shock
  processes*, not in `code/full_model.py`. `rho_Z` stays at 0.80 — the
  Markov-switching estimate disciplines the sovereign-risk shock only.

- **MOSTLY CLOSED: `Y_D` negative for only ONE quarter (issue I-1).**
  `rho_def = 0.9408` took cumulative 40-quarter `Y` from −0.049 to **−2.542** and
  the count of negative-`Y` quarters in the first 40 from 5 to **37**. Residual
  defect: a small positive blip at q2–q4 (+0.0115, +0.0264, +0.0111 — all under
  +0.03% of SS) before `Y` goes negative again at q5 and stays there through q39.
  If a deeper, monotone bust is wanted, the next hypothesis is the `n_inter_D`
  rebound (+1.09 by q3, peaking **+2.94 at q8**), **not** another capital
  friction. The original I-1 write-up below is retained because its two rejected
  hypotheses must not be re-tested.

  *Historical (pre-`rho_def` fix), retained for the rejected hypotheses:* see
  `docs/STATE.md` → *Open issue I-1*. On the default shock `Y_D` was −0.5064,
  −0.0026, **+0.0929**, … then a positive hump; Bi–Foerster–Traum stay negative
  ~20 quarters. **Two frictions were tested and both rejected — do not
  re-test either.**
  1. `chi1` (intermediary capital adjustment cost): raising it makes both the
     trough *and* the rebound bigger. Stays 0.
  2. `omega_I` (investment-flow cost `S(I/I(-1))`, added 2026-08-06): the block
     is **live in `capital_adj_D/F` but calibrated to 0**. Sweeping 0/2/5/10
     shrinks the contraction toward zero rather than lengthening it — at
     `omega_I >= 5`, `Y_D[0]` goes **positive**. Full table in STATE.md.

  Both fail identically: they reallocate the impact between `I` and `C`
  (`C_D[0]` goes −0.51 → +0.23 as `omega_I` rises) without deepening the
  aggregate contraction. **The persistence problem is not a missing investment
  friction.** The `n_inter_D` rebound to **+3.6% by q5** — which gets *larger*,
  not smaller, at every positive `omega_I` — is the more promising next
  hypothesis. Start there, not on another capital-adjustment parameter.

  If a positive `omega_I` is ever adopted, `psi_lambda_B` must be re-tuned:
  peak spread drifts 150.1 → 163–168 bp off the 150 bp target.

- **Regenerating the `omega_I = 0` equivalence reference.** `/tmp/nkpc_irfs_nominal.npz`
  is **stale** (predates the `psi_lambda_B` 8.5 → 7.85 re-tune; differs by 1.56).
  Regenerate with `/opt/anaconda3/envs/ssj/bin/python code/dump_irfs.py OUT.npz`
  *before* making a change, and compare after. The current change passes at
  1.08e-13 over all 45 arrays.

- **The model is sticky-price with nominal deposit contracts.** The `add-nkpc`
  workstream (`docs/superpowers/plans/2026-08-05-nominal-rigidities.md`) is
  **COMPLETE** — Tasks 1–16, all committed, all results regenerated. Read
  *Nominal rigidities (`add-nkpc`) — complete* further down this file before
  touching anything, and `docs/STATE.md`'s top section for the full tables.
  The two things not to rediscover the hard way are `solve_jacobian_padded()`
  (SSJ cannot solve this system without it) and the regime-cache rebuild
  ordering; both are written up in that section.
- Working branch: `add-nkpc` (to be merged to `main`). Production entry point: `code/main.py` (orchestrates
  `calibration.py`, `steady_state.py`, `ic_delta_calibration.py`,
  `depreciation_calibration.py`, `full_model.py`, `tpi.py`, `irf_plots.py`,
  `tpi_plots.py`). The legacy notebook `code/model_v12.ipynb` has been removed.
- Core equations: `code/equations_D.py`, `code/equations_F.py`,
  `code/equations_global.py`.
- Read `docs/STATE.md` first for current model status, calibration table, and
  open issues.
- Read `docs/SPEC.md` for the paper's theoretical framing and research goals
  (merged in from the now-retired `docs/FRAMING_HANDOFF.md`).
- Read `docs/eba_calibration.md` for the EBA-2011 derivation, the identification
  ledger, and the three structural fixes of 2026-07-31 (collateral mapping,
  `omega_K` fund rule, `n_inter` scope). **This is the live calibration.**
- **Policy experiments:** `experiments/` on branch `experiments` — the paper's
  standard results set (E1 backstop schedule, E2 ΔY decomposition, E3 S-1
  writeoff). Spec: `docs/superpowers/specs/2026-08-01-policy-experiments-design.md`;
  plan: `docs/superpowers/plans/2026-08-03-policy-experiments.md`. **COMPLETE
  2026-08-03 — E1, E2, E3 and the orchestrator all landed.** Results:
  `docs/experiments_results.md` (generated). Production regression re-run after
  all of it and **bit-identical** to the pre-work baseline; `git diff main --
  code/` is empty. Full suite 31 passed.

  **Two author decisions now block paper text — see the two bold items below.**
  The schema-3 cache is built (`cache_G_main_v3_*.npz`); rebuild with
  `/opt/anaconda3/envs/ssj/bin/python diagnostics/regimes/regime_model.py --force`
  after any calibration change. Run everything with `experiments/run_all.py` (`--skip-e3` to skip the two re-solves, `--render-only` to rebuild the doc). Results land in `docs/experiments_results.md`. Run E2 alone with
  `/opt/anaconda3/envs/ssj/bin/python experiments/e2_dy_decomposition.py`.
  `experiments/common.py` was hardened after code review the same day:
  `calibration_override` now rejects an unrecognised override key instead of
  silently running a mistyped calibration, and `write_results` refuses to write
  `NaN`.

  **E2's headline finding, which changes how ΔY should be reported:** report the
  decomposition, never the headline ΔY — the channels **cancel** and land on
  different households. (Under the flexible-price model the output response was
  additionally a small *residue* of channels ~4× larger; under sticky prices the
  largest channel is 0.25× the headline, so the magnitude ordering has reversed
  but the instruction has not.) See `docs/STATE.md` for the table.

  **E1's headline:** the loading schedule is monotone decreasing at all 59 finite
  grid points (**4.43 → 1.49** over γ ∈ [0.51, 30.00] on the sticky-price model),
  confirming Live Claim 5 on a fine grid. Every cross-check against
  `code/main.py` passes. Run with
  `/opt/anaconda3/envs/ssj/bin/python experiments/e1_backstop_schedule.py`.

  **First-draft material is ready.** `experiments/paper_outputs.py` →
  **eight** captioned figures in `experiments/paper/` + **four** tables in
  `docs/paper_draft_results.md` (calibration/identification ledger, moment match,
  main results, distributional incidence). Regenerate with
  `/opt/anaconda3/envs/ssj/bin/python experiments/paper_outputs.py`. The decile /
  quintile cache is built separately and rarely:
  `/opt/anaconda3/envs/ssj/bin/python experiments/e4_distribution.py` (~4 min).

  **DIST-1 addressed; no Ginis.** Incidence is reported by **income quintile**,
  and the binning choice is load-bearing. Income bins have mass invariant to the
  shock (verified `max|Δmass| ≈ 1e−19`), so the response is purely behavioural.
  **Wealth** bins do not: masses move 2–3% and the net per-capita number is a
  residue of two nearly-cancelling terms (bottom decile PV: −41.6 consumption vs
  −44.4 mass, netting +2.8). Never describe the wealth cut as household behaviour.

  **Incidence result** (regenerated on the sticky-price model 2026-08-06 — see
  Table 4 of `docs/paper_draft_results.md`, which is authoritative): the crisis is
  progressive. PV consumption **+0.4250%** for the lowest income quintile against
  **−0.9073%** for the highest, monotone in between; the backstop's protection
  runs the same way (**+2.01** vs **+1.34**). The flex-price figures previously
  quoted here (+0.95 / −0.59, gains +0.40 / +0.07) are superseded.

  > **RESOLVED (Task 17, 2026-08-06).** `experiments/paper_outputs.py` no longer
  > has a hardcoded `CAPTIONS` dict. Each figure now builds its own caption from
  > the arrays it plots and hands it to `save()`, which registers it — so
  > `fig08_deciles`'s caption reads the same `pv` object Table 4 does and cannot
  > drift from it. **Figure captions are quotable again**, but quote them from a
  > freshly regenerated `docs/paper_draft_results.md`, not from memory.

  **S-1 RESOLVED 2026-08-04: `writeoff_enabled=0` stays** — the pure risk-premium
  framing. E3 becomes an appendix robustness result and a *stated caveat*: the
  over-compensation claim is conditional on no realised principal writedown.

  **Default-loading split is 3.4% / 96.6%** (fundamental expected loss /
  collateral friction) at the live `psi_lambda_B = 7.85`. `EL_price_D = 0.056134`
  is invariant to `psi_lambda_B`; `psi_spread_D` is *linear* in it, so the
  8.5 → 7.85 sticky-price re-tune moved it 1.737724 → 1.604839 and the split
  3.1%/96.9% → **3.4%/96.6%** (Task 17). The 10.9% / 89% in older sections is
  pre-EBA. Either way 96.6% is a strong version of the constrained-seller claim —
  use it, but **re-derive it whenever `psi_lambda_B` moves.**

  **New: the backstop damps the oscillation, it does not lower the spread path.**
  Cushioning is concentrated at impact; by ~q4 the paths converge and the spread
  ordering reverses (t=8: passive +0.5bp vs aggressive +15.5bp). Do not claim
  uniform compression over the whole path.

  **E3's numbers, for the appendix — read before writing the TPI section.**
  Full writeoff (`writeoff_enabled=1`, `zeta_writeoff=1`) takes `EL_price_D` from
  0.056134 to 0.701743 (12.5×) and takes the loading from 3.82/2.90 to
  **2.46/0.26**. On the sticky-price model the inversion is **partial: medium
  holds above 1 at 2.46, and only aggressive falls below at 0.26** — the
  flex-price model had both below 1 (0.37/0.28). Coupon-only writeoff (`zeta=0`)
  is negligible by contrast (3.77/2.87). So S-1 still decides whether the paper's
  over-compensation result holds *under strong intervention*, but no longer
  overturns it across the schedule.

  **Blocking a paper claim — A5-1's third object is misnamed.** The code reports
  `Σ β^t (pd_passive − pd_intervention)`, which is **negative** because the
  backstop lets Greece run a larger primary deficit, i.e. relaxes austerity. So
  negative = Greece better off, the opposite of what "Greek fiscal saving"
  implies. **Flip the sign or rename it ("austerity relief, PV") before this
  number appears anywhere.** Author decision; magnitudes (0.0015 / 0.0047 PV) are
  unaffected.

  **Careful with `EL_price_D`:** it is **0.056134** at the live calibration, not
  the `0.0717` still quoted in older doc sections and CLAUDE.md. It is the TPI
  loading's denominator — re-derive it, don't copy it.

> **The EBA calibration is LIVE and verified** (established 2026-07-31; the
> *steady-state* content below is current, the *dynamics* were re-measured on the
> sticky-price model 2026-08-06 and are given in the `add-nkpc` section).
> `EBA_CALIBRATION = True`, `BANK_SCOPE = "broad"` in `code/calibration.py`.
>
> Measured: `theta` 5.51/6.94 (GK-eligible assets / CT1), `delta_b` 0.0777/0.0568
> (sovereign maturity ladder repriced at the end-2010 market yield), the sovereign
> book, `K/Y`. Implied: `n_inter` 2.138/1.627 = `(Q*K + sovereign)/theta`, and
> `phi_own` 0.456/0.296. `omega_K = 1` — the passive-fund device is gone.
> Free/tuned: **`psi_lambda_B = 7.85`** (150bp target; was 8.5 until the
> sticky-price re-tune of 2026-08-06), `Delta = 0.2/0.4`, `phi_lamb = 0.15`,
> `mv_rule = 0`.
>
> Verified end-to-end, steady state (unchanged by the sticky-price work):
> `K_D=10.800`/`K_F=10.832` (target 10.8), IC residual −8.9e−16,
> `ca_res_D=6.9e−17`, `rk_D=rk_F=0.010000` (**RK-1 resolved**).
> Dynamics, **current** (sticky prices + nominal deposits, `psi_lambda_B=7.85`):
> `b_gov_D[499]=4.6e−05`, `n_inter_D[0]=−4.296% of SS`, `Y_D[0]=−0.5064%`
> (**Y-1 resolved**), peak spread 150.0bp, TPI loading 3.82/2.90 declining.
> *(The flex-price values were `n_inter_D[0]=−3.380%`, `Y_D[0]=−0.0149%`, peak
> spread 150.4bp, loading 4.35/4.01/3.44.)*
>
> Getting here took three fixes, all documented in `docs/eba_calibration.md`:
> (1) the hidden `ratio=Delta_cross/Delta_own=2.0` closure in
> `ic_delta_calibration`; (2) `omega_K` as a fixed share (`fund_rule=1` makes the
> fund a fixed quantity, SS-identical); (3) the `n_inter` scope — CT1 of the
> stress-test sample is not the whole capital-funding sector.
>
> **Set `BANK_SCOPE="ct1"` to reproduce the CT1 variant**, which is explosive
> (`b_gov[499] ~ 1e2-1e3`) and needs `Delta_own > 0.73`.
>
> **Policy regimes RE-RUN at this calibration (2026-07-31)** — Stage A, Stage B-lite,
> the new certainty-equivalence decomposition and 18/18 tests all exit 0; all three
> figures regenerated. `A_cb=-1.889e-2` (backstop compresses, SA-1 absent),
> `gamma` 12.726/5.080, peak spread 75.2/112.7/**150.3**bp, so the 120–180bp sanity
> band now passes. A6 invariance holds in Stage A and in the lottery.
> `PSILAM_BREAKDOWN` re-derived **2.5 → 15.0** (pole located between 27 and 28;
> guard set from the first pathology at 14–18, not the pole). Doc-sync is enforced
> by `.githooks/pre-commit` — enable once per clone:
> `git config core.hooksPath .githooks`.
>
> **Units fix (2026-07-31), affects published numbers.** SSJ IRFs are *level*
> deviations; `×100` is a percent only where the SS level is ≈1. `n_inter_D_ss=2.138`
> and `K_D_ss=10.8` are not, so `main.py`'s `n_inter_D[0]` print and two Stage A
> figure panels were mislabelled `%` (2.1× and **10×**). Fixed. **Do not quote the
> old −7.227%** — it is the level deviation; the impact is −3.380% of SS. `PT-1`'s
> pass-through is consequently −2.25%/100bp, not ≈−4.5%.
>
> **Open items:** (1) `beliefs.json` dates from 2026-07-23 (estimated MS chain on the
> FRED peripheral–Bund composite; calibration-independent). (2) `Y_D[0]` is positive
> under both intervening regimes and the A5 `dY_D` trough never goes negative —
> output never falls under the backstop. **This is now an order of magnitude larger
> (+0.2008 / +0.8721) and `n_inter_D[0]` has gone positive too (+0.924), i.e. the
> aggressive backstop produces an impact boom** — see the watch item in the
> `add-nkpc` section. Still plausibly linear-rule overshoot; diagnose before
> reporting intervening-regime paths. (3) The `theta`-for-the-whole-sector
> assumption is the one load-bearing judgement left in the bank block; an ECB BSI
> cross-check on bank credit to NFCs would test it. (4) S-1 **RESOLVED 2026-08-04**
> (`writeoff_enabled=0` stays, pure risk-premium framing).

## Quick start

```bash
conda activate ssj   # or use /opt/anaconda3/envs/ssj/bin/python explicitly
cd /path/to/QUANTITATIVE_MODEL
/opt/anaconda3/envs/ssj/bin/python code/main.py
```

Install deps if needed:
```bash
pip install sequence-jacobian numpy scipy matplotlib nbstripout nbdime
nbstripout --install && nbdime config-git --enable
```

Regression test after any equation change — the full pipeline is the regression test:
```bash
/opt/anaconda3/envs/ssj/bin/python code/main.py
```
(`audit_artifacts/` was removed 2026-07-30; it tested a hardcoded calibration, not `get_calibration()`.)

## Latest session (2026-07-24)

- **Ran `code/main.py` end-to-end** (clean: SS `goods_mkt_D≈-4.8e-7`,
  `ca_res_D≈2.3e-16`; `b_gov_D[499]≈2e-6`). Fresh impacts on the default shock:
  `n_inter_D[0]=-2.83%`, `Y_D[0]=+0.032%` (Y-1), peak D–F spread `+0.392pp`. TPI:
  spread compresses 0.392→0.244pp (−38%) over γ=0→10; loading 3.59/3.03/2.47 at
  γ=2/5/10 (matches SPEC Live Claim 1). NB: TPI welfare deltas came out small and
  D-negative under the committed calibration — differs from STATE.md's old
  line-165 (2026-07-16, pre-EBA) numbers; it's the delicate decomposition-dependent
  object SPEC says not to lead with.
- **Corrected a stale `phi_lamb` value in the docs.** Committed value is `0.60`
  (~Bohn); STATE.md's calibration table and IRF-summary header wrongly said `0.15`
  (STATE.md line 37 and `calibration.py` already had 0.60; the table + the
  `calibration.py` comment lagged). Fixed in STATE.md and the `calibration.py`
  comment.
- **PAC sweep → Finding F-2 (STATE.md).** The ~25q ring in the asset-price IRFs is
  the IC/leverage financial accelerator (**|λ|=0.954, ~6yr period, 3.6yr
  half-life**), NOT a portfolio-friction artifact — a 100× PAC change barely moves
  it. No free cosmetic fix; it's a structural financial cycle to describe.
  (`audit_artifacts/pac_sweep.py`.)
- **Added a pre-commit doc-sync hook** (`.claude/settings.json` +
  `.claude/hooks/require-docs-before-commit.sh`): blocks committing model/code
  changes unless STATE.md, PROGRESS.md, HANDOFF.md are updated in the same commit.
- **Created `docs/PROGRESS.md`** — a comprehensive changelog reverse-engineered from
  the 135-commit git history + STATE/audit/EBA docs (the retired `PROCESS.md`'s old "Version
  history" moved here). PROGRESS is now the hook's required changelog; PROCESS stays
  the (rarely-changing) workflow doc.

## Current model state (2026-07-22)

Calibration is EBA-2011-anchored (real GR/DE bank-sovereign concentration data,
not placeholder round numbers), C-1 (the degenerate multi-asset collateral
constraint) is fixed at its root, and `psi_lambda_B` is calibrated to the
paper's external spread target rather than left at a value chosen only to
dodge a bug. Three downstream drift bugs surfaced and were fixed while
verifying the C-1 fix (stale audit harness, a diagnostic sign/scale error, a
real TPI conduit accounting leak). Full details: `docs/eba_calibration.md`,
`docs/STATE.md`.

**Calibration (main, `code/calibration.py`):**
- `phi_bD_D_ss=2.39`, `phi_bF_F_ss=2.76` (own-book concentration, EBA 2011),
  `phi_bF_D_ss=0.018`, `phi_bD_F_ss=0.069` (cross-holdings)
- `n_inter_D=0.408`, `n_inter_F=0.175` (EBA CT1/quarterly-GDP bank capital)
- `omega_K_D=0.0601`, `omega_K_F=0.0190` (passive capital-fund split, added to
  reconcile EBA-thin bank net worth with a plausible aggregate capital stock)
- `Delta_bD_D=0.2`/`Delta_bF_D=0.4` (D), `Delta_bF_F=0.2`/`Delta_bD_F=0.4` (F) —
  genuine hardcoded collateral parameters (C-1 fix), not a degenerate back-solve
- `psi_lambda_B_D/F=1.1793` (recalibrated 2026-07-22 to the 150bp spread target,
  re-tuned from 1.1284 when EL-1 resolved; **do not raise above ~1.5-2.0 without
  re-checking stability** — the model enters a linear-approximation-breakdown
  region there)
- `mv_rule_D/F=1` (market-value fiscal rule), `phi_lamb_D/F=0.60` (~Bohn; governs
  the debt/fiscal mode, not the F-2 financial-accelerator ring)
- `delta_b_D/F=0.10` (2.5yr duration; empirical target 0.036/0.038 not yet
  ported), `f_D/F=0.12` (bank exit rate; bank-cal target 0.03 not yet ported)
- `writeoff_enabled=0` (pure risk-premium framing; S-1, still an author
  decision), `recovery_rate_D/F=0.30` (EL-1 resolved 2026-07-22; NPV Greek-PSI
  recovery, not the old 0.00 placeholder)

**Verified (2026-07-22):**
- Walras clean: `goods_mkt_D/F` ~1e-7, `ca_res_D` ~1e-16, TPI conduit
  `max|ca_res_D|` ~1e-7 across all `gamma`
- Stability: `b_gov_D[499]` ~1e-6/1e-7, `ρ_b(partial-eq.)=0.373` (target <0.95)
- Signs: default shock gives `n_inter_D[0]≈-2.83%` (correct); `Y_D[0]` is a
  small, persistent positive anomaly (Y-1, open)
- TPI: loading (premium PV / expected-loss PV) = 3.59/3.03/2.47 at
  `gamma=2/5/10` — over-compensated, declining (self-extinguishing premium
  claim holds)

## Open issues (author decisions required)

See `docs/STATE.md`'s "Open issues" table for the full, current list (C-1
resolved; S-1, RK-1, Y-1, EL-1, PT-1, DIST-1, A5-1 open). The two most
consequential for the paper right now:

- **EL-1 / recovery rate**: `recovery_rate_D/F=0` should be the actual Greek
  PSI recovery (~25–35%), not zero. It directly sets `EL_price`, one of the two
  pillars of the TPI over-compensation claim.
- **PT-1 / pass-through validation**: the bank-net-worth-to-spread
  pass-through moved ~5× under the EBA recalibration and has never been
  checked against bank-equity/sovereign-spread event studies.

## Next session priorities

1. **Validate PT-1** against event studies (Altavilla–Pagano–Simonelli;
   Acharya–Drechsler–Schnabl) at the new magnitude.
2. **Resolve EL-1** (confirm or justify `recovery_rate_D/F`).
3. **Update `docs/SPEC.md`'s "Theoretical framing" numbers** as PT-1/EL-1 are
   resolved and as the model changes further — several are explicitly flagged
   there as needing re-verification.
4. **Decide S-1** (`writeoff_enabled`), coupled to whatever EL-1 resolves to.
5. **Port remaining bank-cal calibration values**: `delta_b=0.036/0.038`,
   `f=0.03`. Re-test Finding F-1 (`docs/STATE.md`) again after porting — the
   most recent re-test used today's `delta_b=0.10`/`f=0.12` only for the
   duration dimension.
6. **Investigate RK-1 and Y-1** before reporting `rk_F` or `Y_D` impact-sign
   results in the paper.
7. **Re-generate all figures** from `main` after any of the above.

## Important file locations

| File | Purpose |
|------|---------|
| `code/main.py` | Production pipeline (calibration → SS → Jacobian → IRFs → TPI) |
| `code/equations_D.py` | Country D blocks (C-1 fix: `steady_auxilliary_D`) |
| `code/equations_F.py` | Country F blocks (C-1 fix: `steady_auxilliary_F`) |
| `code/equations_global.py` | Global clearing + portfolio costs |
| `code/tpi.py` | TPI/ECB conduit experiment; `cb_pnl` for off-path P&L |
| `docs/STATE.md` | Current model status, calibration table, open issues |
| `docs/SPEC.md` | Research goals, modelling choices, **and paper theoretical framing** |
| `docs/eba_calibration.md` | EBA parameter→moment map; C-1 fix; `psi_lambda_B` recalibration |
| `docs/audit.md` | Master audit log (ranked findings, fix history) |
| `docs/verification_report.md` | Post-fix verification with numerical evidence |
| `docs/bank_cal_review.md` | bank-cal branch analysis; remaining calibration porting roadmap |
| `docs/walras_forensics.md` | Analytical Walras derivation; all leaks proven |
| ~~`audit_artifacts/*`~~ | Removed 2026-07-30 — regression harness, `psi_lambda_B` sweep, F-1/F-2 estimators. Results retained in `docs/STATE.md`; scripts in git history at `0c99013`. |
| `code/tpi_plots.py`, `code/irf_plots.py` | Figure-generation scripts (regenerate from `main`) |
| Overleaf | https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c |

## Nominal rigidities (`add-nkpc`) — complete (Tasks 1–16, 2026-08-05/06)

**The workstream is done.** Sixteen tasks, all committed, all downstream results
regenerated. Full tables in `docs/STATE.md`'s top section; the changelog entry is
in `docs/PROGRESS.md`.

### What the model now is

- **Sticky prices.** Rotemberg price Phillips curves
  `pi = beta*pi(+1) + kappa_p*(mu_p*mc − 1)` in both countries (`price_nkpc_D/F`),
  a markup wedge `w = mu_p*mc*(1−alpha)*Y/N` in `labor_demand_D/F`, and a markup
  rent `profit = (1 − mu_p*mc)*(1−alpha)*Y` (`firm_profit_D/F`) distributed to
  households **in proportion to productivity `e`** via `income_D/F`. **Wages stay
  flexible** — `labor_market_D/F` is untouched, and is *allowed* to be untouched
  precisely because the rent is routed proportional to `e` rather than lump-sum.
- **Nominal closure with no policy rate.** `terms_of_trade` turns the
  monetary-union identity `p/p(-1) = (1+pi_F)/(1+pi_D)` into a residual on the
  existing unknown `p`, pinning the inflation differential; `union_inflation`
  (`omega_pi_D*pi_D + (1−omega_pi_D)*pi_F = 0`) pins the level as the
  `phi_pi → ∞` limit of an ECB rule on union PPI, stated as an abstraction. At
  `omega_pi_D = 0.071` (renormalised capital key) **93% of any terms-of-trade
  move is Greek deflation, 7% German inflation** — the internal-devaluation
  pattern. No financial contract in the model carries a policy rate, so no Fisher
  relation is needed to close it.
- **Nominal deposits.** `i_dep_D/F` is the nominal rate and the solver unknown.
  `rdep_D/F` **keeps its name** as the derived ex-ante real rate (so
  `intermediation_P1`, `divert_bond_foc` and `divert_portfolio_adj` were never
  touched); `rdep_expost_D/F` is the realised rate carrying the inflation
  surprise, consumed by `bank_return_D/F` and `capital_fund_D/F`. That is the
  Fisher channel. `rdep_expost` carries its own `(-1)` — do not double-lag it.
- **Sovereign bonds stay real.** A deliberate asymmetry that maximises banks'
  Fisher exposure. **It must be stated as such in the paper.**
- **27×27 solver system** (was 23×23): `+mc_D, pi_D, mc_F, pi_F` unknowns,
  `+nkpc_p_res_D/F, tot_res, union_pi_res` targets. One block-list definition,
  `full_model.build_block_list()`.
- **The steady state is bit-identical to pre-change** — markups are
  subsidy-neutralised (`mc_ss = 1/mu_p`) and `pi_ss = 0`, so every new object is
  exactly zero at SS.

### The numbers

Impact on the 1pp default shock, % of own SS level. Both columns sit on the same
150bp peak-spread moment, so this is like-for-like:

| | flex, real deposits (`psi_lambda_B=8.5`) | sticky + nominal (7.85) |
|---|---|---|
| peak spread | 150.4 bp | 150.0 bp |
| `Y_D[0]` | −0.0149 | **−0.5064** |
| `C_D[0]` | **+0.2164** | **−0.5103** |
| `I_D[0]` | −0.7718 | −1.0114 |
| `n_inter_D[0]` | −3.3804 | −4.2962 |

Price stickiness alone does most of it (`Y_D[0]` → −0.4923, `C_D[0]` → −0.4904);
nominal deposits add a Fisher amplification ~11× larger on bank net worth than on
output (`n_inter_D[0]` −4.0140 → −4.6155), which is the correct signature.
`psi_lambda_B` was then re-tuned 8.5 → 7.85 to put peak spread back on 150bp.

E1 regime table (regenerated):

| regime | γ | peak spread bp | `Y_D[0]` | `C_D[0]` | `I_D[0]` | `n_inter_D[0]` | loading |
|---|---|---|---|---|---|---|---|
| passive | 0 | 150.1 | −0.5064 | −0.5103 | −1.0114 | −4.296 | n/a |
| medium | 3.2515 | 112.6 | +0.2008 | +0.5285 | −0.2934 | −1.649 | 3.82 |
| aggressive | 9.0163 | 75.1 | +0.8721 | +1.5143 | +0.3977 | **+0.924** | 2.90 |

**Live Claim 5 survives** (loading monotone decreasing, 4.43 → 1.49 over 59 grid
points on γ ∈ [0.51, 30.00], above 1 throughout). **Live Claim 1 survives**
(3.82 / 2.90). E2's identity closes at 3.5e−17 against its 1e−07 assertion.

### Open items

1. **The one-quarter-spike caveat — do not drop it from the write-up.** Output
   and consumption are both positive from quarter 1, and flexible-price
   consumption is in fact *more* persistently negative from quarter 2 on. Nominal
   deposits deepen the impact quarter but do not lengthen the recession
   (`C_D[1]` is essentially unmoved: +0.1141 → +0.1144). Bi-Foerster-Traum's
   output stays negative ~20 quarters. **The honest claim is that the model fixes
   the impact quarter, not that it resolves the investment-bust counterfactual.**
2. **WATCH ITEM: the aggressive backstop now produces an impact boom.**
   `n_inter_D[0] = +0.924` where it was −1.099, with `Y_D[0] = +0.8721` and
   `C_D[0] = +1.5143`. That is a much stronger intervention effect than before and
   a referee will press on it. It may still be linear-rule overshoot at
   `γ = 9.02`, but it now reaches bank net worth, not just output. Diagnose before
   reporting intervening-regime paths.
3. **Candidate follow-ons.** **Nominal sovereign bonds** (would give the sovereign
   an inflation-erosion channel and flip the sign of the bank's net Fisher
   exposure) and a **Sims-Wu loan-in-advance constraint** (Bi-Foerster-Traum's
   persistence device — the one-quarter spike is the symptom it would address).
4. **Prose-vs-table agreement is only partly guarded.** `paper_outputs.py`'s
   `CAPTIONS` dict is fixed (Task 17): it is empty at import and filled by
   `save()` from each figure's own arrays, so a caption cannot outlive the numbers
   it describes, and `main()` asserts fig01's caption and Table 3 agree on impact
   net worth. **`experiments/run_all.py` still has no such assertion**, and there
   is no pytest covering either — the checks are runtime asserts inside the
   generators. A rendered-prose test remains a genuine follow-on. The generic
   lesson stands: **never write a number, a direction, or the word "monotone" into
   a caption as a literal** — derive it, and let a sign flip rewrite the sentence.

### Two things not to rediscover the hard way

1. **SSJ 1.0.0 drops H_Z rows for targets reachable from no shock**, so stock
   `Block.solve_jacobian` cannot solve this system — it returns a 23-row H_Z
   against a 27×27 H_U and numpy raises `size 11500 is different from 13500`.
   Everything routes through `full_model.solve_jacobian_padded()`, which restores
   those rows as zeros (exact, since `dH/dZ` at fixed unknowns is identically zero
   when the shock never appears in the equation). All nine call sites were
   converted in Task 9b, and the invariant
   `grep -rn "\.solve_jacobian(" --include="*.py" code experiments diagnostics |
   grep -v solve_jacobian_padded` must stay **empty**. A 25×25 rewrite was
   considered and rejected — same defect, smaller numbers.
2. **Rebuild the regime cache BEFORE running `experiments/run_all.py`.** The
   experiments never re-solve the model; they read
   `diagnostics/regimes/regime_model.py`'s cached Jacobians. The cache is keyed on
   a hash of the whole live calibration
   (`regime_model._calibration_fingerprint`), so a stale cache can never be picked
   up silently *by name* — but running the experiments first will happily
   re-report the old model. Current tag: `psilam7p85_cal685f7838`. And **E4 is not
   wired into `run_all.py`**: `experiments/e4_distribution.py` is a separate entry
   point feeding `experiments/paper_outputs.py`. Regenerating E1–E3 does not
   regenerate E4 or the eight tracked paper figures.

### Test entry points

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v   # 17 tests, ~1s
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py code/test_eba_calibration.py experiments/ -v   # 40 passed
```

## Run environment

```
/opt/anaconda3/envs/ssj/bin/python   ← always use this
```

Base env has a broken `liblapack` symlink. Each Jacobian solve at current
calibration (T=500) takes ~3 min.
