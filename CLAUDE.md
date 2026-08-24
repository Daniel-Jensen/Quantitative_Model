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

**Fast unit tests** — the sticky-price/nominal-deposit blocks, ~1 s, no model solve. Run
these first; they catch a wiring mistake in a second rather than in twelve minutes:
```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v          # 17 tests
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py code/test_eba_calibration.py experiments/ -v   # 40 passed
```

**Acceptance thresholds** (from `docs/verification_report.md`):
- `goods_mkt_D` ≤ 1e−6 — **corrected 2026-08-18.** The old `1e−14` was never met on this
  calibration: the pre-refactor `main` (91ac778) prints `−4.2493e−07` and the current
  branch prints `−4.2281e−07`. Measured directly against a clean worktree at 91ac778, so
  this is a documentation fix, not a regression. The genuine machine-zero residual is
  `ca_res_D`.
- `goods_mkt_F` ≤ 1e−6 (same story; `−4.18e−07`)
- `ca_res_D` ≤ 1e−13 (actually ~1e−16)
- `deposit_mkt_D/F` ≤ 1e−13
- `nkpc_p_res_D/F`, `tot_res`, `union_pi_res` — exactly `0.000000e+00` at SS
- **GK portfolio FOCs** — `report_gk_steady_state` raises if `nu_i/nu_K ≠ Delta_i_eff` on any
  of the four legs (own legs to 1e−9, cross legs to 5e−4), or if any `Delta_*_eff` leaves
  [0,1]. This runs on every solved SS; it is the check that the sovereign spread is coming
  from intermediary optimality and not from a wedge.

**Targeted audit scripts:** removed with `audit_artifacts/` (2026-07-30). The findings they
produced are recorded in `docs/audit.md` and `docs/STATE.md`; the scripts themselves are in
git history (last present at `0c99013`).

Each Jacobian solve at current calibration (T=500) takes ~3 min.

## Architecture

The model is implemented in the `sequence_jacobian` (SSJ) library. Blocks are defined as `@simple` or `@het` decorated Python functions in three equation files, then assembled and solved by the modular pipeline (`code/main.py`).

### Equation files (edit these; the pipeline imports them)

- `code/equations_D.py` — Country D (Greece): household EGM het block (`hh_D`), deposit return, bank steady-state and intermediation, production, capital, government fiscal, bond pricing/default
- `code/equations_F.py` — Country F (Germany): symmetric analogues of all D blocks
- `code/equations_global.py` — global goods market, external account, bond clearing, the cross-border GK portfolio FOC (`gk_cross_border_foc`), trade balance, bond yield formula, terms of trade, union inflation

### Production pipeline (run this)

- `code/main.py` — orchestrator: calibration → steady state → IC-δ / depreciation calibration → Jacobian + baseline IRFs → TPI experiment → figures. Runs the whole model end-to-end.
- `code/calibration.py`, `code/steady_state.py`, `code/ic_delta_calibration.py`, `code/depreciation_calibration.py`, `code/full_model.py` — the calibration/solve stages `main.py` calls.
- `code/tpi.py`, `code/tpi_plots.py`, `code/irf_plots.py` — TPI experiment and figure generation.

**`full_model.build_block_list()` is the single model definition.** `full_model.py`,
`code/tpi.py` and `diagnostics/regimes/regime_model.py` all call it; the TPI layer supplies
its four `_tpi` swaps (`budget_residual_D/F`, `external_account_D`, `domestic_bond_clearing`)
through `tpi.tpi_overrides()`. Add a block once, there. The three used to hardcode their own
copies of the `sj.create_model([...])` list.

The legacy `code/model_v12.ipynb` has been removed; the modular pipeline above (added in PR #28) is the source of truth. `docs/equation_reconstruction.md` cites notebook cells 2–21 for historical provenance only.

### ⚠ `solve_jacobian_padded()` — never call `Block.solve_jacobian` directly

SSJ 1.0.0's `CombinedBlock._jacobian` seeds from the shock list and ends with
`total_Js[original_outputs & total_Js.outputs, :]`, only visiting blocks whose inputs
intersect that list. **A target that is a pure function of the solver's own unknowns is
therefore silently dropped from H_Z.** All four sticky-price targets (`nkpc_p_res_D/F`,
`tot_res`, `union_pi_res`) are exactly that, so stock SSJ returns a 23-row H_Z against the
27×27 H_U and `np.linalg.solve` dies with `size 11500 is different from 13500`.

`full_model.solve_jacobian_padded()` restores the missing rows as zeros. This is **exact, not
an approximation** — `dH/dZ` at fixed unknowns is identically zero when the shock symbol never
appears in the equation — and it otherwise mirrors `Block.solve_jacobian` line-for-line,
printing the padded row names on every solve so the padding can never go silent.

**Every Jacobian call site in the repo routes through it.** The invariant:

```bash
grep -rn "\.solve_jacobian(" --include="*.py" code experiments diagnostics | grep -v solve_jacobian_padded
```

must stay **empty**. A 25×25 rewrite folding the four targets into existing equations was
considered and rejected: it hits the identical defect with smaller numbers.

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
- `experiments/e3_writeoff_s1.py` — S-1 and the payoff specification. Rebased 2026-08-18: baseline is now `zeta_writeoff = 1`, so the variants are `e3a_realised_writeoff` (`writeoff_enabled = 1` — the pure S-1 test) and `e3b_coupon_only_pricing` (`zeta = 0` — the §12 Arm-3 diagnostic). Both are exactly SS-neutral; `zeta` is allocation-neutral but dynamically decisive.
- `experiments/common.py` — cache access, `calibration_override`, unit helpers, provenance.

**`run_all.py` covers E1–E3 only. E4 is a separate entry point.**
`experiments/e4_distribution.py` (~4 min) builds `cache_e4_deciles.npz`, which feeds
`experiments/paper_outputs.py` → the eight tracked `experiments/paper/fig0*.png` and
`docs/paper_draft_results.md`. Regenerating E1–E3 does **not** regenerate any of that. Run
both, or you will ship paper artefacts built on the previous model — which is exactly what
nearly happened in the sticky-price workstream.

It runs on `diagnostics/regimes/regime_model.py`'s cached Jacobian response matrices, which
are built from the production equation files — **no copy of the model or the calibration
lives in this package**, which is the failure that made the retired `audit_artifacts/`
harness silently test a different model for weeks. Rebuild the cache after any calibration
change, and **rebuild it *before* running the experiments** (they never re-solve the model,
so the reverse order silently re-reports the old one):
`/opt/anaconda3/envs/ssj/bin/python diagnostics/regimes/regime_model.py --force`.

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
- **Predetermined deposit rate:** the rate is locked at t−1 and deposit contracts are non-contingent. Since 2026-08-06 the contracted rate is **nominal**: `i_dep_D/F` is the solver unknown, `rdep_D/F` keeps its name as the derived **ex-ante** real rate (t → t+1), and `rdep_expost_D/F` is the realised rate carrying the inflation surprise. Using a period-t rate instead was T-2, the critical doom-loop sign inversion. Note `rdep_expost` carries its own `(-1)` internally — writing `rdep_expost_D(-1)` double-lags it.
- **Hatchondo-Martinez perpetuity:** bond coupon decays at rate `1−delta_b`; duration ≈ 1/delta_b quarters. This is what generates MTM capital losses on bank balance sheets.

**Sovereign pricing is structural — no wedges (2026-08-18, `gk-structural-foc` stages 2–5).**
Read `docs/STATE.md` → *GK structural refactor* before touching the bank block.

- **One source of truth for the bond payoff.** `bond_return_D/F` emits three objects:
  `rb_exp_D` (expected payoff over the default distribution — the ONLY return the pricing
  equations read), `rb_actual_D` (the realised return on the branch the IRF traces, gated by
  `writeoff_enabled_D = 0`), and `EL_load_D` (a pure diagnostic, read by `code/tpi.py`'s CB
  P&L and by nothing else). Never price off `rb_actual`.
- **`zeta_writeoff_D/F = 1`.** A default writes down the perpetuity's continuation value as
  well as its coupon. Coupon-only pricing under-states the loss on a 12.9-quarter claim by
  ~12.6×. `zeta_writeoff` governs what is PRICED; `writeoff_enabled` governs what is
  REALISED. They are independent and both are needed.
- **`psi_spread_D/F` and `EL_price_D/F` are DELETED, not recalibrated.** So are
  `divert_bond_foc_D/F`, `divert_portfolio_adj`, `domestic_bond_foc_D/F`,
  `portfolio_adj_cost`, `bond_price_ss_D/F` and the `excess_return_*_ss` anchors.
  `code/test_nkpc_blocks.py::test_no_ad_hoc_sovereign_spread_wedge_anywhere` AST-scans
  `code/*.py` and fails if any of those names reappears in live code. There must be no
  equation of the form `spread += parameter * default_probability`.
- **The spread comes from the GK portfolio FOC.** `gk_bond_foc_D/F` impose
  `nu_own/nu_K = Delta_own_eff`, which with `intermediation_P1_D/F` is
  `rb_exp(+1) − rdep = Delta_eff · (rk(+1) − rdep)`. `q_b_D` and `q_b_F` are **SS unknowns**
  pinned by these two residuals (`rb_D_res`, `rb_F_res`), not by `bond_price_ss_*`.
- **Cross-border legs get the same FOC plus a portfolio adjustment cost.**
  `gk_cross_border_foc` (in `equations_global.py`) states
  `nu_cross/nu_K = Delta_cross_eff + psi·(b − b_ss)`, divided through by `SDF_banker·Omega_p1`
  so `psi_bF_D`/`psi_bD_F` keep their calibrated units. Four proportionality conditions
  cannot hold against two bond prices, so the own legs pin PRICES and the cross legs pin
  QUANTITIES. The `psi` costs load on the bond stock, carry no `def_rate`, and are zero at
  the calibrated position — they are not spread wedges.
- **All four `Delta` are 0.20.** Cross-border `Delta_bF_D`/`Delta_bD_F` moved 0.40 → 0.20
  because GK optimality at a riskless SS with `rk_D = rk_F` and `rdep_D = rdep_F = 0` forces
  it. Holding 0.40 leaves an 80bp/yr constant cross-border wedge.
- **`psi_lambda_B_D/F` is SLATED FOR DELETION (decided 2026-08-24), not a diagnostic arm.**
  It is 0 in the live calibration and will never again be set otherwise: the parameter and
  the risk-sensitive branch of `collateral_quality_D/F` are to be removed outright, leaving
  `Delta_*_eff = Delta_*` as constants. The Greek episode gives no independent observable for
  a sovereign-specific haircut *elasticity*, and keeping a switched-off dial invites the
  reader to ask what it does. This does not remove collateral: the IC still binds and
  `Delta_bD_D = 0.20` still makes sovereign paper LESS divertable than capital — i.e. BETTER
  collateral, which is why it earns a fifth of the capital premium. (The former wording here
  said "worse collateral"; that contradicted the code comment in `intermediation_IC_D` and
  the IC algebra `theta_tgt = value/lambda_gk + (1−Delta)·phi_b`, and is corrected.)
  **The paper's introduction now states that the model contains no device making a
  balance-sheet friction a function of the default probability.** Until the deletion lands,
  `code/` and the paper disagree; the deletion is the fix, not a caveat in the text. The old
  `psi_lambda_B = 3.01` counterfactual is void and must never be described as "the
  non-fundamental share of the spread".
- **Do not write "x% fundamental / y% non-fundamental".** This is a linearised equilibrium
  model; the channels operate jointly. Describe mechanisms: direct expected-default-loss
  pricing, intermediary balance-sheet amplification, and (if enabled) risk-sensitive
  collateral amplification.
- **Walras redundancy:** `ca_res_D` and `goods_mkt_F` are *dropped* from the solver target system (not a bug). Post-fix they hold to machine tolerance; monitoring them is the primary regression check.
- **p-conversion in F-bank returns:** F-bank's D-bond book is denominated in D-goods; returns must be converted via `p(-1)/p` to F-goods before entering the F-goods budget constraint (W-2 fix). Missing this causes `goods_mkt_F` to leak up to 2% of GDP.

**Nominal rigidities (added 2026-08-06 on `add-nkpc`; see `docs/SPEC.md` for the full
rationale and `docs/STATE.md` for numbers):**

- **Rotemberg price NKPC, subsidy-neutralised.** `pi = beta*pi(+1) + kappa_p*(mu_p*mc − 1)` in both countries; markup wedge `w = mu_p*mc*(1−alpha)*Y/N` in `labor_demand_D/F`. **Wages stay flexible** — `labor_market_D/F` is untouched. `mc_ss = 1/mu_p` neutralises the markup, so `mu_p*mc = 1`, `profit_ss = 0`, `pi_ss = 0` and **the steady state is bit-identical to the flex-price model**. `kappa_p → ∞` recovers flexible prices exactly — the standing equivalence gate.
- **Markup rent proportional to `e`, not lump-sum.** `firm_profit_D/F` route `(1 − mu_p*mc)(1−alpha)Y` through `income_D/F` in proportion to productivity, which makes labour-plus-profit income per unit of `e` exactly `(1−alpha)Y·e` — identical to the flex-price model, so the wedge bites only on hiring, never on household income. A lump-sum rebate was rejected as countercyclical (it would manufacture a progressive incidence result as an artifact of the rebate rule).
- **No policy rate.** The union-inflation normalisation `omega_pi_D*pi_D + (1−omega_pi_D)*pi_F = 0` (`phi_pi → ∞` limit of an ECB rule on union PPI) pins the inflation level; `p/p(-1) = (1+pi_F)/(1+pi_D)` pins the differential off the existing unknown `p`. No contract in the model carries a policy rate, so no Fisher relation is needed. `omega_pi_D = 0.071` is the renormalised capital key, **not** GDP weights — GDP weights would erase the 93/7 Greek-deflation split.
- **Nominal deposits against REAL sovereign bonds** — a deliberate asymmetry that maximises banks' Fisher exposure (nominal debtors, real creditors). **Must be stated as a modelling choice in the paper.** Nominal sovereign bonds are a candidate extension, not a correction.
- **27×27 solver system** (was 23×23): `+mc_D, pi_D, mc_F, pi_F` unknowns (and `rdep_D/F → i_dep_D/F`), `+nkpc_p_res_D/F, tot_res, union_pi_res` targets.

## Branch convention

- `add-nkpc` — **the current working branch.** Sticky prices + nominal deposit contracts (Tasks 1–16, 2026-08-05/06), to be merged to `main`.
- `main` — contains all six structural fixes (W-1, W-2, W-3, T-2, A-2, TPI-1, merged via PR #27) plus the modular-file reorganisation (PR #28).
- `audit` / `AB-audit` — historical audit branches. `AB-audit` was merged into `main` (PR #27); `audit` (PR #26) was closed as superseded. Do not reuse.
- `bank-cal` — old calibration branch predating structural fixes. **Do not merge.** Port calibration values only (see `docs/bank_cal_review.md`).

## Current model state and open issues

See `docs/STATE.md` for the full calibration table. Key tensions:

| Issue | Description |
|-------|-------------|
| **C-1** | **RESOLVED (2026-07-22).** Was: `Delta_cross=1.45>1`, back-solved divertable fraction exceeds 1, multi-asset IC degenerate. Fixed at its root: `steady_auxilliary_D/F` now solve `lambda_gk` from the multi-asset IC directly; `Delta_bD_D/F=0.2/0.4` are genuine hardcoded inputs, verified to bind exactly. See `docs/eba_calibration.md`. |
| **S-1** | **SUPERSEDED 2026-08-18 by the GK structural refactor.** The *framing* decision stands — `writeoff_enabled = 0`, the IRF traces the no-default branch, `def_rate` is a genuine probability and agents price the expected loss. What changed is that **`zeta_writeoff = 1` is now the baseline**, so the priced loss covers the perpetuity's principal/continuation value and not just the current coupon. E3's old finding that this drives the TPI loading below 1 is therefore no longer a robustness variant — it *is* the baseline: **loading 0.520 / 0.504 / 0.482 at γ = 2/5/10**, and near-invariant across all three arms. Any paper text claiming ECB over-compensation must be rewritten, not caveated. `EL_price_D` is deleted; the diagnostic loading is the endogenous `EL_load_D = 0.7014` (was 0.056134 under coupon-only pricing). Full-writeoff no longer breaks the named-regime construction — peak spread is monotone in γ in the baseline. |
| **GK-1** | **RESOLVED (2026-07-31) — collateral mapping.** The GK block is well-posed only if `f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross`. At measured EBA moments `Delta_own=0.2` violated it by −1.26/−1.42, giving **negative** `lambda_gk`/`Omega` while the solver converged with machine-zero residuals (C-1's silent-degeneracy mode). Cause: `_ic_delta`'s hidden `ratio=Delta_cross/Delta_own=2.0` back-solve closure, which capped `Delta_own<=0.5` against a required `>~0.73`. Removed; `Delta` is now free and the IC **residual** is checked directly. Guarded by `steady_state.assert_gk_well_posed` on every solved SS. **UPDATED 2026-08-18:** the cross-border `Delta` is no longer free — imposing the cross-border portfolio FOC forces `Delta_bF_D = Delta_bD_F = Delta_own = 0.20` at a riskless SS with `rk_D = rk_F` and `rdep_D = rdep_F = 0`. **Live: all four `Delta = 0.20`, `lambda_gk_D = 2.1087`, `Omega_D = 10.3462`.** Anyone computing the collateral channel from the old 0.85/0.90 or from `Delta_cross = 0.40` is wrong. |
| **GK-2** | **RESOLVED (2026-07-31) — `n_inter` scope.** Three compounding amplifiers made the CT1-scope EBA calibration explosive. Fixed in order: the hidden `ratio=2.0` closure (GK-1); `omega_K` as a *fixed share* (new `fund_rule=1` → fund holds a fixed quantity, `dK/dN = theta` not `theta/omega_K`, steady state identical); and finally the **scope of `n_inter`** — CT1 is the stress-test sample, not the agent intermediating the whole capital stock. New **`BANK_SCOPE="broad"`**: `n_inter = (Q*K + sovereign)/theta`, `omega_K = 1`, fund device gone. Model is stable and on target. *(The dynamic numbers originally recorded here — spread 150.4bp, `Y_D[0]=-0.0149%`, loading 4.35/4.01/3.44 — were flex-price; current values are in the Calibration row below.)* `rk_D=rk_F=0.010000` (RK-1 resolved), Y-1 resolved. |
| **EBA switch** | `EBA_CALIBRATION` in `code/calibration.py` is **`True` and LIVE** since 2026-07-31, with `BANK_SCOPE="broad"`. GK-2's explosive dynamics are fixed; the `False` branch keeps the pre-EBA values as a fallback. The moment set (`code/eba_calibration.py` → `data/eba_moments.json`) is rebuilt, identified, and tested (10/10). |
| **Calibration** | **LIVE: EBA + `BANK_SCOPE="broad"` + sticky prices + nominal deposits + the GK structural refactor (2026-08-18).** **`psi_lambda_B_D/F = 0`** — the 150bp-moment tuning history is void; it was conditional on a payoff that under-priced default 12.6×. `zeta_writeoff_D/F = 1`, `writeoff_enabled_D/F = 0`, all four `Delta = 0.20`, `recovery_rate = 0.30` (EL-1), `EL_load_D = 0.7014` (endogenous, replaces the deleted `EL_price_D`). `n_inter = 2.138/1.627`, `omega_K = 1.0`, `phi_lamb = 0.15`, `mv_rule = 0`, `delta_b = 0.0777/0.0568`. Solved SS: `q_b_D = 0.974906`, `q_b_F = 0.965974`, both yields 80.0bp annualised, **SS spread exactly 0**, `beta_D = 0.999534`, `rdep_D = 0`, `rk_D = rk_F = 0.010000`. Sticky-price block unchanged: `mu_p_D/F = 1.20`, `mc_D/F = 1/1.20`, `kappa_p_D/F = 0.0871`, `pi_D/F = 0.0`, `omega_pi_D = 0.071`, `i_dep_D/F = 0.000`. Dynamics on a 1pp default shock, **none of it tuned to a moment**: peak spread **205.87bp**, `Y_D[0] = −1.974%`, `C_D[0] = −2.511%`, `n_inter_D[0] = −11.407%`, German yield **−16.2bp** (flight to quality), `b_DD[0] = +2.116%` with `K_D[0] = −0.067%` (balance-sheet crowding out), `b_DF[0] = −2.101%` (German banks retrench). TPI loading **0.520 / 0.504 / 0.482**, declining. |
| **F-1** | `mv_rule_D/F` **committed at 0 (par)**. The near-unit-root zone `phi_lamb≈0.15-0.18` that F-1 identified under `mv_rule=1` is **not mild — it is a hard break**, measured directly 2026-07-30: `mv_rule=1` at the pre-EBA `phi_lamb=0.15` gives `n_inter_D[0]=-1554%`, `Y_D[0]=+0.17%` (perverse sign), `b_gov_D[499]=1.6e-2`. It needs `phi_lamb=0.60` to stay healthy (`n_inter_D[0]=-5.89%`, `Y_D[0]=-0.024%`, `b_gov_D[499]=0.0`). **`mv_rule=1` and `phi_lamb=0.15` are not a usable pair** — porting empirical duration is a two-parameter move. See `docs/STATE.md` Finding F-1. |

## Typical iteration

1. Edit equation files (`equations_D.py`, `equations_F.py`, `equations_global.py`). New blocks go in `full_model.build_block_list()` — one place, not three.
2. Run the fast tests first: `/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v` (17 tests, ~1 s). Then re-run the pipeline: `/opt/anaconda3/envs/ssj/bin/python code/main.py` (calibration → steady state → Jacobian → IRFs → TPI).
3. Inspect residuals: `goods_mkt_D`, `goods_mkt_F`, `ca_res_D`, `deposit_mkt_D/F` — all ≤ 1e−7 — plus the four sticky-price targets `nkpc_p_res_D`, `nkpc_p_res_F`, `tot_res`, `union_pi_res`, which must be exactly `0.000000e+00` at SS (they are zero by construction at `mu_p*mc = 1`, `pi = 0`; anything else means the subsidy neutralisation broke).
4. Verify default shock: `n_inter_D[0]` and `Y_D[0]` must both fall (positive = timing bug). Current values: `−11.4073%` and `−1.9742%` of SS. Also check the GK FOC table `main.py` now prints — all four legs must read `OK`.
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
| `docs/cb_mechanism.md` | **CANONICAL (2026-08-19).** What the TPI central bank is, how it transmits, the 2×2 sovereign-holdings matrix, and the reporting rules that follow. Supersedes ad-hoc CB descriptions elsewhere. Key results: the spread identity `rb_exp_D(+1) − rdep_D = 0.20·(rk_D(+1) − rdep_D)` means TPI's spread and investment effects are **one** effect; ~72% of the CB book is bought from **German** banks; concentration relief is a net-worth not a quantity effect; the "closed-loop pole" is a T=500 truncation artefact. Evidence in `diagnostics/cb_audit/` |
| `docs/experiments_results.md` | **GENERATED — do not hand-edit.** Standard policy results: E1 backstop schedule, E2 ΔY decomposition, E3 S-1 writeoff. Regenerate with `experiments/run_all.py` (`--skip-e3` skips the two model re-solves, `--render-only` rebuilds from results on disk) |
| `docs/paper_draft_results.md` | **GENERATED — do not hand-edit.** First-draft tables and figures. Emitted by `experiments/paper_outputs.py`, which needs `experiments/e4_distribution.py`'s cache first. **Its figure captions are currently stale** — see `docs/STATE.md`'s open items. |
| `docs/superpowers/plans/2026-08-05-nominal-rigidities.md` | Implementation plan for the sticky-price / nominal-deposit workstream (Tasks 1–16) |
| `docs/superpowers/specs/2026-08-01-policy-experiments-design.md` | Design spec for the `experiments/` package |
| `docs/superpowers/plans/2026-08-03-policy-experiments.md` | Implementation plan for the same |
