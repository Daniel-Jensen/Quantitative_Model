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
- Read `docs/eba_calibration.md` for the C-1 structural fix (still live) and the
  EBA-2011 calibration derivation — **historical**: the calibration was reverted
  to pre-EBA values on 2026-07-30.

> **Current state (2026-07-31).** Calibration is still **pre-EBA**
> (`psi_lambda_B=3.0`, `n_inter=3.0`, `omega_K=1.0`, `phi_lamb=0.15`, `mv_rule=0`,
> cross-holdings 0.25) plus EL-1's `recovery_rate=0.30` — now selected by
> `EBA_CALIBRATION = False` in `code/calibration.py`. Verified bit-exact against
> the 2026-07-30 values and by a full `main.py` run. Spread 187.2bp ann vs the
> 150bp target.
>
> **The EBA moment set was rebuilt 2026-07-31 and is identified** — maturity
> ladder → `delta_b`, GK-eligible assets → `theta`, measured EAD → `omega_K`,
> Acharya–Steffen MTM for amplification (the 2011 adverse scenario is rejected).
> Flipping `EBA_CALIBRATION = True` turns the whole set on in one line.
>
> **Steady state now works; dynamics do not.** Two problems, one solved:
>
> 1. **SOLVED — collateral mapping.** The GK block needs
>    `f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross`, i.e.
>    `Delta_own > ~0.73`. `Delta=0.2` gave `lambda_gk_D=-0.087`, `Omega_D=-0.301`
>    while the solver converged cleanly. The culprit was `_ic_delta`'s hidden
>    `ratio=Delta_cross/Delta_own=2.0` closure, now removed — `Delta` is free and
>    the IC residual is checked directly. At **0.85/0.90**: `lambda_gk_D=+0.927`
>    (vs pre-EBA +0.923), `Omega_D=+4.62`, `K_D=10.80`.
> 2. **OPEN — dynamic instability.** `b_gov_D[499] ~ 1e2–1e3`. Amplification is
>    `theta*phi_own = 13.17` vs 1.0 for the placeholder. Not fiscal (flat in
>    `phi_lamb` to 25), not the friction (present at `psi_lambda_B=0`). `chi1`
>    0→0.5 cuts peak spread 1.1e7bp→6.0bp but no value removes the root.
>
> Read `docs/STATE.md` → *EBA REBUILD* and `docs/eba_calibration.md` →
> *GK feasibility* / *Dynamic instability* before touching this.
> `steady_state.assert_gk_well_posed` makes the steady-state failure loud.
>
> **Policy-regime feature runs end-to-end** (Stage A + Stage B-lite + unit tests, all
> exit 0). Doc-sync is enforced by `.githooks/pre-commit` - enable once per clone:
> `git config core.hooksPath .githooks`.
>
> A6 lottery invariance was fixed 2026-07-31 (ranked post-revelation, checked at both
> amplifier settings, with a noise margin) — it now genuinely holds.
>
> **Open items:** (1) `psi_lambda_B=3.0` gives 187.2bp, outside `run_regimes.py`'s own
> 120-180bp band - retuning to the 150bp target is unfinished. (2) `delta_b=0.10` still
> short of the empirical 7yr/6.5yr - porting needs `mv_rule=1` **and** `phi_lamb=0.60`
> together. (3) `beliefs.json` predates the calibration revert (2026-07-23).
> (4) EBA calibration to be revisited on a new branch.

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
