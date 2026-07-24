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
- Read `docs/eba_calibration.md` for the EBA-2011-anchored calibration
  derivation, the C-1 structural fix, and the `psi_lambda_B` recalibration.

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

Regression test after any equation change (~6 min):
```bash
/opt/anaconda3/envs/ssj/bin/python audit_artifacts/run_audit.py
```

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
  the 135-commit git history + STATE/audit/EBA docs (`PROCESS.md`'s old "Version
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
| `docs/PROCESS.md` | Workflow, debugging, regression test |
| `docs/audit.md` | Master audit log (ranked findings, fix history) |
| `docs/verification_report.md` | Post-fix verification with numerical evidence |
| `docs/bank_cal_review.md` | bank-cal branch analysis; remaining calibration porting roadmap |
| `docs/walras_forensics.md` | Analytical Walras derivation; all leaks proven |
| `audit_artifacts/run_audit.py` | Full Walras/sign/C-1 regression pipeline |
| `audit_artifacts/psilam_moment_sweep_postC1.py` | `psi_lambda_B` moment-matching sweep |
| `audit_artifacts/philamb_sweep_postC1_fine_v2.py` | Finding F-1 stability re-test (order-selected Prony estimator) |
| `audit_artifacts/pac_sweep.py` | Finding F-2: ring eigenvalue vs PAC (Prony complex-pair extractor) |
| `code/tpi_plots.py`, `code/irf_plots.py` | Figure-generation scripts (regenerate from `main`) |
| Overleaf | https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c |

## Run environment

```
/opt/anaconda3/envs/ssj/bin/python   ← always use this
```

Base env has a broken `liblapack` symlink. Each Jacobian solve at current
calibration (T=500) takes ~3 min.
