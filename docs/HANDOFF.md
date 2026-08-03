# Handoff Notes

## Where to start

- Working branch: `main`. Production entry point: `code/main.py` (orchestrates
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
  plan: `docs/superpowers/plans/2026-08-03-policy-experiments.md`. **In progress
  — E1 and E2 have landed; E3 and the orchestrator are next.**
  The schema-3 cache is built (`cache_G_main_v3_*.npz`); rebuild with
  `/opt/anaconda3/envs/ssj/bin/python diagnostics/regimes/regime_model.py --force`
  after any calibration change. Run everything with `experiments/run_all.py` (`--skip-e3` to skip the two re-solves, `--render-only` to rebuild the doc). Results land in `docs/experiments_results.md`. Run E2 alone with
  `/opt/anaconda3/envs/ssj/bin/python experiments/e2_dy_decomposition.py`.
  `experiments/common.py` was hardened after code review the same day:
  `calibration_override` now rejects an unrecognised override key instead of
  silently running a mistyped calibration, and `write_results` refuses to write
  `NaN`.

  **E2's headline finding, which changes how ΔY should be reported:** the output
  response is the small residue of an investment channel and a net-export channel
  each ~4× larger and opposite in sign. Report the decomposition, never the
  headline ΔY. See `docs/STATE.md` for the table.

  **E1's headline:** the loading schedule is monotone decreasing at all 59 finite
  grid points (4.51 → 2.07 over γ=0.5→30), confirming Live Claim 5 on a fine grid.
  Every cross-check against `code/main.py` passes. Run with
  `/opt/anaconda3/envs/ssj/bin/python experiments/e1_backstop_schedule.py`.

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

> **Current state (2026-07-31). The EBA calibration is LIVE and verified.**
> `EBA_CALIBRATION = True`, `BANK_SCOPE = "broad"` in `code/calibration.py`.
>
> Measured: `theta` 5.51/6.94 (GK-eligible assets / CT1), `delta_b` 0.0777/0.0568
> (sovereign maturity ladder repriced at the end-2010 market yield), the sovereign
> book, `K/Y`. Implied: `n_inter` 2.138/1.627 = `(Q*K + sovereign)/theta`, and
> `phi_own` 0.456/0.296. `omega_K = 1` — the passive-fund device is gone.
> Free/tuned: `psi_lambda_B = 8.5` (150bp target), `Delta = 0.2/0.4`,
> `phi_lamb = 0.15`, `mv_rule = 0`.
>
> Verified end-to-end: `K_D=10.800`/`K_F=10.832` (target 10.8), IC residual
> −8.9e−16, `ca_res_D=6.9e−17`, `b_gov_D[499]=1.4e−05`,
> `n_inter_D[0]=−3.380% of SS` (level dev −7.227 — see the units fix below),
> `Y_D[0]=−0.0149%` (**Y-1 resolved**), `rk_D=rk_F=0.010000` (**RK-1 resolved**),
> peak spread 150.4bp, TPI loading 4.35/4.01/3.44 declining.
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
> output never falls under the backstop. At `gamma_aggressive=12.7` this is likely
> linear-rule overshoot; check before reporting intervening-regime output paths.
> (3) The `theta`-for-the-whole-sector assumption is the one load-bearing judgement
> left in the bank block; an ECB BSI cross-check on bank credit to NFCs would test
> it. (4) S-1 (`writeoff_enabled=0`) still an author decision.

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

## Run environment

```
/opt/anaconda3/envs/ssj/bin/python   ← always use this
```

Base env has a broken `liblapack` symlink. Each Jacobian solve at current
calibration (T=500) takes ~3 min.
