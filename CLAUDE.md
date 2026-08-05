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

**Structural regression test** — the full pipeline is the regression test. Run after any
equation change and inspect the printed residuals:
```bash
/opt/anaconda3/envs/ssj/bin/python code/main.py
```

> The former `audit_artifacts/` harness (`run_audit.py` + targeted scripts and JSON logs)
> was removed on 2026-07-30. It carried its own hardcoded copy of the calibration rather
> than importing `get_calibration()`, so it silently tested a *different* model than
> `code/main.py` and its results were misleading. Recover from git history if needed.

**Acceptance thresholds** (from `docs/verification_report.md`):
- `goods_mkt_D` ≤ 1e−14
- `goods_mkt_F` ≤ 1e−7
- `ca_res_D` ≤ 1e−7
- `deposit_mkt_D/F` ≤ 1e−13

**Targeted audit scripts:** removed with `audit_artifacts/` (2026-07-30). The findings they
produced are recorded in `docs/audit.md` and `docs/STATE.md`; the scripts themselves are in
git history (last present at `0c99013`).

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

### Policy experiments (`experiments/`, added 2026-08-03)

The paper's standard results set. **`code/` is deliberately untouched by this package** so
`code/main.py` stays usable as the regression test.

- `experiments/run_all.py` — runs everything, renders `docs/experiments_results.md`.
  `--skip-e3` avoids E3's two model re-solves (~11 min); `--render-only` rebuilds the
  document from results already on disk.
- `experiments/e1_backstop_schedule.py` — named regimes (γ **solved** for 0/25/50%
  peak-spread compression), A5-1's three German objects reported separately, loading
  schedule, welfare labelled secondary.
- `experiments/e2_dy_decomposition.py` — ΔY against the `market_clearing_D` identity;
  self-verifying, asserts closure at 1e−7.
- `experiments/e3_writeoff_s1.py` — the S-1 writeoff variants.
- `experiments/common.py` — cache access, `calibration_override`, unit helpers, provenance.

It runs on `diagnostics/regimes/regime_model.py`'s cached Jacobian response matrices, which
are built from the production equation files — **no copy of the model or the calibration
lives in this package**, which is the failure that made the retired `audit_artifacts/`
harness silently test a different model for weeks. Rebuild the cache after any calibration
change: `/opt/anaconda3/envs/ssj/bin/python diagnostics/regimes/regime_model.py --force`.

**Two gotchas worth knowing before extending it.** `calibration_override` patches the
*module attribute*, so a module-level `from calibration import get_calibration` binds the
original and silently misses the override — import the module and resolve at use time.
And percentages must divide by their own SS level (`common.pct_of_ss`): `n_inter_D_ss=2.138`
and `K_D_ss=10.8` are not ≈1, and a past bug mislabelled exactly those by 2.1× and 10×.

### Routines

- `routines/grids.py` — deposit and income grids; supports both standard Rouwenhorst Markov chains and GMAR discrete-time process (loaded from `Discretisation/Outputs/`)
- `routines/income.py`, `routines/calculate_gini.py` — income process and distributional statistics

### Audit artifacts

Removed 2026-07-30 (see *Running and testing*). `code/main.py` is now the only regression
path; findings live in `docs/audit.md` and `docs/STATE.md`.

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
| **S-1** | **RESOLVED (author decision, 2026-08-04): `writeoff_enabled=0` stays.** The paper commits to the **pure risk-premium framing** — `def_rate` is a genuine probability, agents price the expected loss, and the IRF traces the no-default branch. This is a standard risk-premium-shock device, *not* "default is impossible", and must be stated as such. E3 (`experiments/e3_writeoff_s1.py`) quantified the alternative before the decision: `writeoff_enabled=1` alone is SS-neutral and negligible (loading 4.00→3.93), but adding `zeta_writeoff=1` takes `EL_price_D` 0.0561→0.7017 (12.5×) and collapses the loading to **0.37/0.28 — below 1**, inverting Live Claim 1. It also breaks the named-regime construction (peak spread stops being monotone in γ). Retained as an appendix robustness result: the over-compensation claim is **conditional on no realised principal writedown**, and that conditionality must be stated in the paper, not buried. `recovery_rate_D/F=0.30` (EL-1, Greek PSI NPV framing) stays live through `EL_price`. **`EL_price_D` is 0.056134 at the live calibration, not the 0.0717 previously recorded here** — that predates the EBA `delta_b=0.0777`/`q_b=0.969`. Re-derive it wherever quoted; it is the loading's denominator. |
| **GK-1** | **RESOLVED (2026-07-31) — collateral mapping.** The GK block is well-posed only if `f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross`. At measured EBA moments `Delta_own=0.2` violated it by −1.26/−1.42, giving **negative** `lambda_gk`/`Omega` while the solver converged with machine-zero residuals (C-1's silent-degeneracy mode). Cause: `_ic_delta`'s hidden `ratio=Delta_cross/Delta_own=2.0` back-solve closure, which capped `Delta_own<=0.5` against a required `>~0.73`. Removed; `Delta` is now free and the IC **residual** is checked directly. `Delta=0.85/0.90` → `lambda_gk_D=+0.927` (pre-EBA: +0.923). Guarded by `steady_state.assert_gk_well_posed` on every solved SS. |
| **GK-2** | **RESOLVED (2026-07-31) — `n_inter` scope.** Three compounding amplifiers made the CT1-scope EBA calibration explosive. Fixed in order: the hidden `ratio=2.0` closure (GK-1); `omega_K` as a *fixed share* (new `fund_rule=1` → fund holds a fixed quantity, `dK/dN = theta` not `theta/omega_K`, steady state identical); and finally the **scope of `n_inter`** — CT1 is the stress-test sample, not the agent intermediating the whole capital stock. New **`BANK_SCOPE="broad"`**: `n_inter = (Q*K + sovereign)/theta`, `omega_K = 1`, fund device gone. Model is stable and on target: `b_gov_D[499]=1.4e-05`, spread 150.4bp, `Y_D[0]=-0.0149%` (Y-1 resolved), `rk_D=rk_F=0.010000` (RK-1 resolved), TPI loading 4.35/4.01/3.44 declining. |
| **EBA switch** | `EBA_CALIBRATION` in `code/calibration.py`, default **False** (pre-EBA values, bit-exact, solves). `True` turns on the rebuilt measured moment set; SS is then correct but dynamics are explosive (GK-2). The moment set itself (`code/eba_calibration.py` → `data/eba_moments.json`) is rebuilt, identified, and tested (10/10). |
| **Calibration** | **Reverted to pre-EBA values 2026-07-30** (`psi_lambda_B=3.0`, `n_inter=3.0`, `omega_K=1.0`, `phi_lamb=0.15`, `mv_rule=0`, cross-holdings `0.25`), keeping all structural fixes and the EL-1 `recovery_rate=0.30`. The EBA 2011 anchoring (`phi_bD_D_ss=2.39`, `psi_lambda_B=1.1793`) is **no longer live** — `docs/eba_calibration.md` is now historical. Spread response is **187.2bp annualised** per 1pp default shock vs the paper's 150bp target (~25% over). `delta_b=0.10` (2.5yr) is still empirically short; `0.036/0.038` (7yr/6.5yr) requires `mv_rule=1` (see F-1). **The `psi_lambda_B<1.5` breakdown warning was EBA-specific** (thin net worth `n_inter=0.408`); at pre-EBA `n_inter=3.0` the documented breakdown region is ~4-5 and 3.0 runs clean. |
| **F-1** | `mv_rule_D/F` **committed at 0 (par)**. The near-unit-root zone `phi_lamb≈0.15-0.18` that F-1 identified under `mv_rule=1` is **not mild — it is a hard break**, measured directly 2026-07-30: `mv_rule=1` at the pre-EBA `phi_lamb=0.15` gives `n_inter_D[0]=-1554%`, `Y_D[0]=+0.17%` (perverse sign), `b_gov_D[499]=1.6e-2`. It needs `phi_lamb=0.60` to stay healthy (`n_inter_D[0]=-5.89%`, `Y_D[0]=-0.024%`, `b_gov_D[499]=0.0`). **`mv_rule=1` and `phi_lamb=0.15` are not a usable pair** — porting empirical duration is a two-parameter move. See `docs/STATE.md` Finding F-1. |

## Typical iteration

1. Edit equation files (`equations_D.py`, `equations_F.py`, `equations_global.py`).
2. Re-run the pipeline: `/opt/anaconda3/envs/ssj/bin/python code/main.py` (calibration → steady state → Jacobian → IRFs → TPI).
3. Inspect residuals: `goods_mkt_D`, `goods_mkt_F`, `ca_res_D`, `deposit_mkt_D/F` — all ≤ 1e−7.
4. Verify default shock: `n_inter_D[0]` and `Y_D[0]` must both fall (positive = timing bug).
5. Confirm the IC-δ consistency check and Walras residuals printed by `main.py` are unchanged.
6. Update the living docs after any calibration or structural change — **STATE.md, PROGRESS.md (changelog entry), HANDOFF.md** (not just CLAUDE.md). This is **enforced** by two hooks that block the commit otherwise:
   - `.claude/hooks/require-docs-before-commit.sh` — PreToolUse gate, fires when Claude Code runs the commit.
   - `.githooks/pre-commit` — git-native twin, covers terminal commits. **Enable once per clone: `git config core.hooksPath .githooks`.**

   Both fire only when the commit stages `code/**` or any `*.py`; doc-only commits pass. Bypass a false positive with `git commit --no-verify`. Keep the required-doc set in the two files in sync. (`docs/PROCESS.md` was retired 2026-07-30, superseded by PROGRESS.md.)
7. Commit the changed `.py` files, with the doc updates in the same commit.

## Docs reference

| File | Contains |
|------|----------|
| `docs/STATE.md` | Current calibration table, Walras residuals, open issues, next priorities |
| `docs/PROGRESS.md` | Changelog — dated development timeline (git history + findings); one entry per code commit (convention; hook not installed) |
| `docs/SPEC.md` | Research goals, functional requirements, modelling choices, calibration targets, **and the paper's theoretical framing/narrative** (merged in from the retired `docs/FRAMING_HANDOFF.md`) |
| `docs/eba_calibration.md` | **REBUILT 2026-07-31.** Identified EBA parameter→moment map (maturity ladder→`delta_b`, GK-eligible assets→`theta`, measured EAD→`omega_K`, Acharya–Steffen MTM), the **identification ledger** (identified / bounded / still-free / deliberately-rejected), and the **GK feasibility** finding (GK-1). The 2026-07-22 build is retained below it as history. |
| `docs/HANDOFF.md` | Quick-start, session priorities, important file locations |
| `docs/audit.md` | Master audit log: all findings ranked by severity, fix history, open hypotheses |
| `docs/walras_forensics.md` | Analytical derivation of all three Walras leaks and their proofs |
| `docs/bank_cal_review.md` | bank-cal branch analysis; calibration porting roadmap |
| `docs/verification_report.md` | Post-fix numerical verification with residual tables |
| `docs/experiments_results.md` | **GENERATED — do not hand-edit.** Standard policy results: E1 backstop schedule, E2 ΔY decomposition, E3 S-1 writeoff. Regenerate with `experiments/run_all.py` (`--skip-e3` skips the two model re-solves, `--render-only` rebuilds from results on disk) |
| `docs/superpowers/specs/2026-08-01-policy-experiments-design.md` | Design spec for the `experiments/` package |
| `docs/superpowers/plans/2026-08-03-policy-experiments.md` | Implementation plan for the same |
