# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Two-country heterogeneous-agent New Keynesian model with Gertler-Karadi financial intermediaries and sovereign debt, calibrated to the 2010–2012 Greek sovereign debt crisis. Application: the ECB's Transmission Protection Instrument (TPI). Primary output is a research paper (Overleaf: https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c).

## Environment

Always use `/opt/anaconda3/envs/ssj/bin/python`. The base Anaconda environment has a broken `liblapack` symlink that causes silent numerical failures.

```bash
conda activate ssj
jupyter notebook code/model_v12.ipynb
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

The model is implemented in the `sequence_jacobian` (SSJ) library. Blocks are defined as `@simple` or `@het` decorated Python functions in three equation files, then assembled and solved in the notebook.

### Equation files (edit these; notebook imports them)

- `code/equations_D.py` — Country D (Greece): household EGM het block (`hh_D`), deposit return, bank steady-state and intermediation, production, capital, government fiscal, bond pricing/default
- `code/equations_F.py` — Country F (Germany): symmetric analogues of all D blocks
- `code/equations_global.py` — global goods market, external account, bond clearing, portfolio adjustment costs, trade balance, bond yield formula

### Active notebook

- `code/model_v12.ipynb` — calibration cell, steady-state solve, Jacobian computation, IRFs (TFP + default shocks), TPI policy experiment, welfare calculation

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
| **C-1** | `Delta_cross=1.45>1`: back-solved divertable fraction exceeds 1; multi-asset IC is degenerate. Preferred resolution: hardcode `Delta_D=0.2, Delta_F=0.4` per bank-cal branch. |
| **S-1** | `writeoff_enabled=0`: default shock produces no realized bank losses. Model is currently a pure risk-premium loop. Enabling writeoff (`writeoff_enabled=1`, `recovery_rate=0.40`) gives the balance-sheet doom loop. Author decision pending — **now coupled to F-1: with the market-value fiscal rule, writeoff must stay OFF (risk-premium framing), else the default response is perverse.** |
| **Calibration** | `delta_b_D/F=0.10` (2.5yr) is empirically too short; target is `0.036/0.038` (7yr/6.5yr GR/DE). Porting from bank-cal is the next major task (see `docs/bank_cal_review.md`). **But empirical long duration is explosive under the par-value fiscal rule at every `phi_lamb` (F-1) — it requires `mv_rule=1`.** |
| **F-1** | Par-value Bohn rule is explosive with empirical long-duration bonds at all `phi_lamb`. The switchable **market-value rule** (`mv_rule_D/F`, default 0=par) restores stationarity at `phi_lamb≈0.10` in the risk-premium framing. See `docs/STATE.md` Finding F-1. |

## Typical iteration

1. Edit equation files (`equations_D.py`, `equations_F.py`, `equations_global.py`).
2. Restart notebook kernel and re-run calibration → steady-state → Jacobian cells.
3. Inspect residuals: `goods_mkt_D`, `goods_mkt_F`, `ca_res_D`, `deposit_mkt_D/F` — all ≤ 1e−7.
4. Verify default shock: `n_inter_D[0]` and `Y_D[0]` must both fall (positive = timing bug).
5. Run `audit_artifacts/run_audit.py` to confirm no regression.
6. Update `docs/STATE.md` after any calibration or structural change.
7. Commit cleaned notebook (nbstripout strips outputs automatically).

## Docs reference

| File | Contains |
|------|----------|
| `docs/STATE.md` | Current calibration table, Walras residuals, open issues, next priorities |
| `docs/SPEC.md` | Research goals, functional requirements, modelling choices, calibration targets |
| `docs/PROCESS.md` | Workflow, debugging steps, EBA verification assertions |
| `docs/HANDOFF.md` | Quick-start, session priorities, important file locations |
| `docs/audit.md` | Master audit log: all findings ranked by severity, fix history, open hypotheses |
| `docs/walras_forensics.md` | Analytical derivation of all three Walras leaks and their proofs |
| `docs/bank_cal_review.md` | bank-cal branch analysis; calibration porting roadmap |
| `docs/verification_report.md` | Post-fix numerical verification with residual tables |
