# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Two-country heterogeneous-agent New Keynesian model with Gertler-Karadi financial intermediaries and sovereign debt, calibrated to the 2010–2012 Greek sovereign debt crisis. Application: the ECB's Transmission Protection Instrument (TPI). Primary output is a research paper (Overleaf: https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c).

## Environment

Always use `/opt/anaconda3/envs/ssj/bin/python`. The base Anaconda environment has a broken `liblapack` symlink that causes silent numerical failures.

```bash
conda activate ssj
/opt/anaconda3/envs/ssj/bin/python code/main.py
```

Install dependencies if needed:
```bash
pip install sequence-jacobian numpy scipy matplotlib nbstripout nbdime
nbstripout --install && nbdime config-git --enable
```

## Running and testing

**Structural regression test** — run after any equation change; prints max Walras residuals across all shocks (~6 min total):
```bash
/opt/anaconda3/envs/ssj/bin/python audit_artifacts/run_audit.py
```

**Acceptance thresholds** (from `docs/verification_report.md`):
- `goods_mkt_D` ≤ 1e−14
- `goods_mkt_F` ≤ 1e−7
- `ca_res_D` ≤ 1e−7
- `deposit_mkt_D/F` ≤ 1e−13

**Targeted audit scripts:**
```bash
/opt/anaconda3/envs/ssj/bin/python audit_artifacts/fix_test.py        # W-1/W-2 Walras repair
/opt/anaconda3/envs/ssj/bin/python audit_artifacts/tpi_test.py        # TPI CB accounting
/opt/anaconda3/envs/ssj/bin/python audit_artifacts/philamb_test.py    # phi_lamb stability sweep
/opt/anaconda3/envs/ssj/bin/python audit_artifacts/bankcal_stability_test.py  # low-amplification probe
```

Each Jacobian solve at current calibration (T=500) takes ~3 min.

## Architecture

The model is implemented in the `sequence_jacobian` (SSJ) library. Blocks are defined as `@simple` or `@het` decorated Python functions in three equation files, then assembled and solved by the modular pipeline (`code/main.py`).

### Equation files (edit these; the pipeline imports them)

- `code/equations_D.py` — Country D (Greece): household EGM het block (`hh_D`), deposit return, bank steady-state and intermediation, production, capital, government fiscal, bond pricing/default
- `code/equations_F.py` — Country F (Germany): symmetric analogues of all D blocks
- `code/equations_global.py` — global goods market, external account, bond clearing, portfolio adjustment costs, trade balance, bond yield formula

### Production pipeline (run this)

- `code/main.py` — orchestrator: calibration → steady state → IC-δ / depreciation calibration → Jacobian + baseline IRFs → TPI experiment → figures. Runs the whole model end-to-end.
- `code/calibration.py`, `code/steady_state.py`, `code/ic_delta_calibration.py`, `code/depreciation_calibration.py`, `code/full_model.py` — the calibration/solve stages `main.py` calls.
- `code/tpi.py`, `code/tpi_plots.py`, `code/irf_plots.py` — TPI experiment and figure generation.

The legacy `code/model_v12.ipynb` has been removed; the modular pipeline above (added in PR #28) is the source of truth. `docs/equation_reconstruction.md` cites notebook cells 2–21 for historical provenance only.

### Routines

- `routines/grids.py` — deposit and income grids; supports both standard Rouwenhorst Markov chains and GMAR discrete-time process (loaded from `Discretisation/Outputs/`)
- `routines/income.py`, `routines/calculate_gini.py` — income process and distributional statistics

### Audit artifacts

- `audit_artifacts/run_audit.py` — full regression pipeline (the canonical post-fix verification tool)
- `audit_artifacts/*.py` — targeted tests for individual bugs (W-1/W-2, TPI-1, phi_lamb sweep)
- `audit_artifacts/*.json` — result logs from each audit run

## Key modelling choices

These are deliberate design decisions — do not "fix" them without checking `docs/SPEC.md`:

- **`Y = F(K_t)` (current-period capital):** production uses same-period capital stock; capital producer receives `mpk·(K−K(-1))` to close capital income accounting (W-1 fix). The alternative `K(-1)` timing eliminates this term but is equally valid.
- **Predetermined deposit rate:** `Rgross = (1+rdep(-1))·P(-1)/P`. Deposit contracts are non-contingent — the rate is locked at t−1. Using `rdep` (a period-t unknown) instead was T-2, the critical doom-loop sign inversion.
- **Hatchondo-Martinez perpetuity:** bond coupon decays at rate `1−delta_b`; duration ≈ 1/delta_b quarters. This is what generates MTM capital losses on bank balance sheets.
- **Walras redundancy:** `ca_res_D` and `goods_mkt_F` are *dropped* from the solver target system (not a bug). Post-fix they hold to machine tolerance; monitoring them is the primary regression check.
- **p-conversion in F-bank returns:** F-bank's D-bond book is denominated in D-goods; returns must be converted via `p(-1)/p` to F-goods before entering the F-goods budget constraint (W-2 fix). Missing this causes `goods_mkt_F` to leak up to 2% of GDP.

## Branch convention

- `main` — **use this for all new work**. Contains all six structural fixes (W-1, W-2, W-3, T-2, A-2, TPI-1, merged via PR #27) plus the modular-file reorganisation (PR #28).
- `audit` / `AB-audit` — historical audit branches. `AB-audit` was merged into `main` (PR #27); `audit` (PR #26) was closed as superseded. Do not reuse.
- `bank-cal` — old calibration branch predating structural fixes. **Do not merge.** Port calibration values only (see `docs/bank_cal_review.md`).

## Current model state and open issues

See `docs/STATE.md` for the full calibration table. Key tensions:

| Issue | Description |
|-------|-------------|
| **C-1** | **RESOLVED (2026-07-22).** Was: `Delta_cross=1.45>1`, back-solved divertable fraction exceeds 1, multi-asset IC degenerate. Fixed at its root: `steady_auxilliary_D/F` now solve `lambda_gk` from the multi-asset IC directly; `Delta_bD_D/F=0.2/0.4` are genuine hardcoded inputs, verified to bind exactly. See `docs/eba_calibration.md`. |
| **S-1** | `writeoff_enabled=0`: default shock produces no realized bank losses. Model is currently a pure risk-premium loop. Enabling writeoff (`writeoff_enabled=1`, `recovery_rate=0.40`) gives the balance-sheet doom loop. Author decision pending — **now coupled to F-1: with the market-value fiscal rule, writeoff must stay OFF (risk-premium framing), else the default response is perverse.** Also coupled to `recovery_rate_D/F` currently being a placeholder (0.00, not the actual Greek PSI ~25-35%) — see `docs/STATE.md` issue EL-1. |
| **Calibration** | EBA 2011 bank-sovereign concentration (`phi_bD_D_ss=2.39` etc.) and `psi_lambda_B=1.1284` (data-disciplined to the 150bp spread target, 2026-07-22) are now live — see `docs/eba_calibration.md`. `delta_b_D/F=0.10` (2.5yr) is still empirically too short; target is `0.036/0.038` (7yr/6.5yr GR/DE), not yet ported. **`psi_lambda_B` must not be raised above ~1.5-2.0 without re-checking stability — the model enters a linear-approximation-breakdown region there on the current calibration (both the pre-EBA value 2.8 and the original default 3.0 now sit inside it).** |
| **F-1** | Re-tested 2026-07-22 on the C-1-fixed, EBA-anchored model with a validated (order-selected Prony/eigenvalue) stability estimator — the original energy-ratio-proxy estimate was itself found to be an overfitting artifact partway through. Result: the market-value rule (`mv_rule_D/F`, default 0=par, currently committed at 1) is stable across nearly all of `phi_lamb∈[0.05,0.25]` with empirical long-duration bonds, except a narrow, mild zone around `phi_lamb≈0.15-0.18` — do not reuse the pre-fix "`phi_lamb≈0.10`" plateau claim, it sits inside a range not re-confirmed at the old paper's exact risk-premium parameterization. See `docs/STATE.md` Finding F-1 for all three re-test rounds. |

## Typical iteration

1. Edit equation files (`equations_D.py`, `equations_F.py`, `equations_global.py`).
2. Re-run the pipeline: `/opt/anaconda3/envs/ssj/bin/python code/main.py` (calibration → steady state → Jacobian → IRFs → TPI).
3. Inspect residuals: `goods_mkt_D`, `goods_mkt_F`, `ca_res_D`, `deposit_mkt_D/F` — all ≤ 1e−7.
4. Verify default shock: `n_inter_D[0]` and `Y_D[0]` must both fall (positive = timing bug).
5. Run `audit_artifacts/run_audit.py` to confirm no regression.
6. Update `docs/STATE.md` after any calibration or structural change.
7. Commit the changed `.py` files.

## Docs reference

| File | Contains |
|------|----------|
| `docs/STATE.md` | Current calibration table, Walras residuals, open issues, next priorities |
| `docs/SPEC.md` | Research goals, functional requirements, modelling choices, calibration targets, **and the paper's theoretical framing/narrative** (merged in from the retired `docs/FRAMING_HANDOFF.md`) |
| `docs/eba_calibration.md` | EBA 2011 parameter→moment map; the C-1 structural fix; `psi_lambda_B` recalibration and its breakdown-region warning |
| `docs/PROCESS.md` | Workflow, debugging steps, EBA verification assertions |
| `docs/HANDOFF.md` | Quick-start, session priorities, important file locations |
| `docs/audit.md` | Master audit log: all findings ranked by severity, fix history, open hypotheses |
| `docs/walras_forensics.md` | Analytical derivation of all three Walras leaks and their proofs |
| `docs/bank_cal_review.md` | bank-cal branch analysis; calibration porting roadmap |
| `docs/verification_report.md` | Post-fix numerical verification with residual tables |
