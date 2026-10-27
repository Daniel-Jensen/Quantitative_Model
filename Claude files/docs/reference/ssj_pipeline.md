# Sequence-space (SSJ) pipeline: reference

Relocated from CLAUDE.md on 2026-10-26. Load this file when working on `linear_ssj/`,
`experiments/` or `diagnostics/regimes/`. Interpreter: `/opt/anaconda3/envs/ssj/bin/python` only.
Text marked `[corrected 2026-10-26]` was checked against the repository during the move; the rest
is verbatim. The SSJ calibration (`calibration/ssj.py`) is separate from the global one; the SSJ
open-issue register and its dated numbers are in `STATE.md` Part I.

[STATUS 2026-10-26] The viva paper describes the global solver. CLAUDE.md formerly said the paper's
current figures run on this pipeline; that may still hold for the PRE_VIVA draft and
`experiments/paper/`, but it is not true of the viva paper.

## Running and testing

### A. Sequence-space pipeline (`linear_ssj/`)

**Structural regression test** — the full pipeline is the regression test. Run after any
equation change and inspect the printed residuals:
```bash
/opt/anaconda3/envs/ssj/bin/python -m linear_ssj.main     # or: python3 run.py ssj
```

> The former `audit_artifacts/` harness (`run_audit.py` + targeted scripts and JSON logs)
> was removed on 2026-07-30. It carried its own hardcoded copy of the calibration rather
> than importing `get_calibration()`, so it silently tested a *different* model than
> `linear_ssj/main.py` and its results were misleading. Recover from git history if needed.

**Fast unit tests** — the sticky-price/nominal-deposit blocks, ~1 s, no model solve. Run
these first; they catch a wiring mistake in a second rather than in twelve minutes:
```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest "Claude files/linear_ssj/tests/test_nkpc_blocks.py" -v          # 22 test functions (static count 2026-10-26; was "17 tests")
/opt/anaconda3/envs/ssj/bin/python -m pytest "Claude files" -v  # every pytest-style test: 51 pass as recorded (not re-run 2026-10-26); the 9 in test_eba_calibration need the raw EBA files in data/
```

**Acceptance thresholds** (from the former `Claude files/docs/verification_report.md`, deleted 2026-10-14; git history):
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
produced are recorded in `Claude files/docs/STATE.md` (and in the former `audit.md`, deleted 2026-10-14); the scripts themselves are in
git history (last present at `0c99013`).


### A. Sequence-space architecture (`linear_ssj/`)

The model is implemented in the `sequence_jacobian` (SSJ) library. Blocks are defined as `@simple` or `@het` decorated Python functions in three equation files, then assembled and solved by the modular pipeline (`linear_ssj/main.py`).

### Equation files (edit these; the pipeline imports them)

- `linear_ssj/model/equations_D.py` — Country D (Greece): household EGM het block (`hh_D`), deposit return, bank steady-state and intermediation, production, capital, government fiscal, bond pricing/default
- `linear_ssj/model/equations_F.py` — Country F (Germany): symmetric analogues of all D blocks
- `linear_ssj/model/equations_global.py` — global goods market, external account, bond clearing, the cross-border GK portfolio FOC (`gk_cross_border_foc`), trade balance, bond yield formula, terms of trade, union inflation

### Production pipeline (run this)

- `linear_ssj/main.py` — orchestrator in two stages: `compute` (calibration → steady state → IC-δ / depreciation calibration → Jacobian + baseline IRFs → TPI experiment → save to `results/SSJ/data/`) and `plot` (every figure, from the saved data, into `results/SSJ/figures/`). `--stage all` (the default) runs both; `run.py ssj` runs them for you.
- `calibration/ssj.py`, `linear_ssj/solve/steady_state.py`, `linear_ssj/solve/ic_delta_calibration.py`, `linear_ssj/solve/depreciation_calibration.py`, `linear_ssj/model/full_model.py` — the calibration/solve stages `main.py` calls.
- `linear_ssj/solve/tpi.py`, `linear_ssj/reporting/tpi_plots.py`, `linear_ssj/reporting/irf_plots.py` — TPI experiment and figure generation.

**`full_model.build_block_list()` is the single model definition.** `full_model.py`,
`linear_ssj/solve/tpi.py` and `diagnostics/regimes/regime_model.py` all call it; the TPI layer supplies
its four `_tpi` swaps (`budget_residual_D/F`, `external_account_D`, `domestic_bond_clearing`)
through `tpi.tpi_overrides()`. Add a block once, there. The three used to hardcode their own
copies of the `sj.create_model([...])` list.

The legacy `code/model_v12.ipynb` has been removed; the modular pipeline above (added in PR #28) is the source of truth. (The notebook-era `equation_reconstruction.md` was retired 2026-10-01; git history.)

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
grep -rn "\.solve_jacobian(" --include="*.py" linear_ssj experiments diagnostics | grep -v solve_jacobian_padded
```

must stay **empty**. [corrected 2026-10-26]: the search formerly named `code/`, which now holds only `__pycache__`; the SSJ model lives in `linear_ssj/`. A 25×25 rewrite folding the four targets into existing equations was
considered and rejected: it hits the identical defect with smaller numbers.

### Policy experiments (`experiments/`, added 2026-08-03)

The paper's standard results set. **`linear_ssj/` is deliberately untouched by this package** so
`linear_ssj/main.py` stays usable as the regression test.

- `experiments/run_all.py` — runs everything, renders `Claude files/docs/experiments_results.md`.
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
`Claude files/docs/paper_draft_results.md`. Regenerating E1–E3 does **not** regenerate any of that. Run
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
*module attribute*, so a module-level `from calibration.ssj import get_calibration` binds the
original and silently misses the override — import the module and resolve at use time.
And percentages must divide by their own SS level (`common.pct_of_ss`): `n_inter_D_ss=2.138`
and `K_D_ss=10.8` are not ≈1, and a past bug mislabelled exactly those by 2.1× and 10×.

### Routines

- `routines/grids.py` — deposit and income grids; supports both standard Rouwenhorst Markov chains and GMAR discrete-time process (loaded from `Discretisation/Outputs/`)
- `routines/income.py`, `routines/calculate_gini.py` — income process and distributional statistics

### Audit artifacts

Removed 2026-07-30 (see *Running and testing*). `linear_ssj/main.py` is now the only regression
path; findings live in `Claude files/docs/STATE.md` (the former `audit.md` was deleted 2026-10-14; git history).

## Key modelling choices (SSJ)

These are deliberate design decisions — do not "fix" them without checking `Claude files/docs/SPEC.md`:

- **`Y = F(K_t)` (current-period capital):** production uses same-period capital stock; capital producer receives `mpk·(K−K(-1))` to close capital income accounting (W-1 fix). The alternative `K(-1)` timing eliminates this term but is equally valid.
- **Predetermined deposit rate:** the rate is locked at t−1 and deposit contracts are non-contingent. Since 2026-08-06 the contracted rate is **nominal**: `i_dep_D/F` is the solver unknown, `rdep_D/F` keeps its name as the derived **ex-ante** real rate (t → t+1), and `rdep_expost_D/F` is the realised rate carrying the inflation surprise. Using a period-t rate instead was T-2, the critical doom-loop sign inversion. Note `rdep_expost` carries its own `(-1)` internally — writing `rdep_expost_D(-1)` double-lags it.
- **Hatchondo-Martinez perpetuity:** bond coupon decays at rate `1−delta_b`; duration ≈ 1/delta_b quarters. **Duration alone does NOT generate the MTM capital losses — duration *interacting with* `rho_s` (persistence of the latent default-risk factor) does.** Corrected 2026-08-28; the former wording ("This is what generates MTM capital losses") credited duration by itself and is wrong. Holding duration at the calibrated 17 quarters and varying persistence, the perpetuity and a one-quarter bill lose *exactly the same* on the reported impulse when `rho_s = 0` (both −1.01%), diverge mildly at `rho_s = 0.5` (−1.27% vs −1.01%), and separate only at the calibrated `rho_s = 0.95` (−5.58% vs −1.01%, a factor of 5.5). A transitory risk shock marks a long bond down no more than a short one. The loss is large because a persistent `s_t` keeps the hazard elevated across the claim's whole remaining life. Two corollaries: a short bond's *principal* is marked down too, so the distinction is never coupon-vs-principal; and any text describing the MTM channel must carry the persistence, not just the maturity.

**Sovereign pricing is structural — no wedges (2026-08-18, `gk-structural-foc` stages 2–5).**
Read `Claude files/docs/STATE.md` → *GK structural refactor* before touching the bank block.

- **One source of truth for the bond payoff.** `bond_return_D/F` emits three objects:
  `rb_exp_D` (expected payoff over the default distribution — the ONLY return the pricing
  equations read), `rb_actual_D` (the realised return on the branch the IRF traces, gated by
  `writeoff_enabled_D = 0`), and `EL_load_D` (a pure diagnostic, read by `linear_ssj/solve/tpi.py`'s CB
  P&L and by nothing else). Never price off `rb_actual`.
- **`zeta_writeoff_D/F = 1`.** A default writes down the perpetuity's continuation value as
  well as its coupon. Coupon-only pricing under-states the loss on a 12.9-quarter claim by
  ~12.6×. `zeta_writeoff` governs what is PRICED; `writeoff_enabled` governs what is
  REALISED. They are independent and both are needed.
- **`psi_spread_D/F` and `EL_price_D/F` are DELETED, not recalibrated.** So are
  `divert_bond_foc_D/F`, `divert_portfolio_adj`, `domestic_bond_foc_D/F`,
  `portfolio_adj_cost`, `bond_price_ss_D/F` and the `excess_return_*_ss` anchors.
  `Claude files/linear_ssj/tests/test_nkpc_blocks.py::test_no_ad_hoc_sovereign_spread_wedge_anywhere` AST-scans
  `linear_ssj/**/*.py` + `calibration/ssj.py` and fails if any of those names reappears in live code. There must be no
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
  `linear_ssj/` and the paper disagree; the deletion is the fix, not a caveat in the text. The old
  `psi_lambda_B = 3.01` counterfactual is void and must never be described as "the
  non-fundamental share of the spread".
- **Do not write "x% fundamental / y% non-fundamental".** This is a linearised equilibrium
  model; the channels operate jointly. Describe mechanisms: direct expected-default-loss
  pricing, intermediary balance-sheet amplification, and (if enabled) risk-sensitive
  collateral amplification.
- **Walras redundancy:** `ca_res_D` and `goods_mkt_F` are *dropped* from the solver target system (not a bug). Post-fix they hold to machine tolerance; monitoring them is the primary regression check.
- **p-conversion in F-bank returns:** F-bank's D-bond book is denominated in D-goods; returns must be converted via `p(-1)/p` to F-goods before entering the F-goods budget constraint (W-2 fix). Missing this causes `goods_mkt_F` to leak up to 2% of GDP.


**Nominal rigidities (added 2026-08-06 on `add-nkpc`; see `Claude files/docs/SPEC.md` for the full
rationale and `Claude files/docs/STATE.md` for numbers):**

- **Rotemberg price NKPC, subsidy-neutralised.** `pi = beta*pi(+1) + kappa_p*(mu_p*mc − 1)` in both countries; markup wedge `w = mu_p*mc*(1−alpha)*Y/N` in `labor_demand_D/F`. **Wages stay flexible** — `labor_market_D/F` is untouched. `mc_ss = 1/mu_p` neutralises the markup, so `mu_p*mc = 1`, `profit_ss = 0`, `pi_ss = 0` and **the steady state is bit-identical to the flex-price model**. `kappa_p → ∞` recovers flexible prices exactly — the standing equivalence gate.
- **Markup rent proportional to `e`, not lump-sum.** `firm_profit_D/F` route `(1 − mu_p*mc)(1−alpha)Y` through `income_D/F` in proportion to productivity, which makes labour-plus-profit income per unit of `e` exactly `(1−alpha)Y·e` — identical to the flex-price model, so the wedge bites only on hiring, never on household income. A lump-sum rebate was rejected as countercyclical (it would manufacture a progressive incidence result as an artifact of the rebate rule).
- **No policy rate.** The union-inflation normalisation `omega_pi_D*pi_D + (1−omega_pi_D)*pi_F = 0` (`phi_pi → ∞` limit of an ECB rule on union PPI) pins the inflation level; `p/p(-1) = (1+pi_F)/(1+pi_D)` pins the differential off the existing unknown `p`. No contract in the model carries a policy rate, so no Fisher relation is needed. `omega_pi_D = 0.071` is the renormalised capital key, **not** GDP weights — GDP weights would erase the 93/7 Greek-deflation split.
- **Nominal deposits against REAL sovereign bonds** — a deliberate asymmetry that maximises banks' Fisher exposure (nominal debtors, real creditors). **Must be stated as a modelling choice in the paper.** Nominal sovereign bonds are a candidate extension, not a correction.
- **27×27 solver system** (was 23×23): `+mc_D, pi_D, mc_F, pi_F` unknowns (and `rdep_D/F → i_dep_D/F`), `+nkpc_p_res_D/F, tot_res, union_pi_res` targets.


## Typical iteration (SSJ)

1. Edit equation files (`equations_D.py`, `equations_F.py`, `equations_global.py`). New blocks go in `full_model.build_block_list()` — one place, not three.
2. Run the fast tests first: `/opt/anaconda3/envs/ssj/bin/python -m pytest "Claude files/linear_ssj/tests/test_nkpc_blocks.py" -v` (22 test functions, ~1 s). Then re-run the pipeline: `/opt/anaconda3/envs/ssj/bin/python -m linear_ssj.main` (calibration → steady state → Jacobian → IRFs → TPI).
3. Inspect residuals: `goods_mkt_D`, `goods_mkt_F`, `ca_res_D`, `deposit_mkt_D/F` — all ≤ 1e−7 — plus the four sticky-price targets `nkpc_p_res_D`, `nkpc_p_res_F`, `tot_res`, `union_pi_res`, which must be exactly `0.000000e+00` at SS (they are zero by construction at `mu_p*mc = 1`, `pi = 0`; anything else means the subsidy neutralisation broke).
4. Check the default shock: `n_inter_D[0]` (D bank net worth on impact) and `Y_D[0]` (D output on impact) are expected to fall. A rise must be investigated before anything is reported; the first suspect is deposit-rate timing (the T-2 inversion). Values recorded 2026-08-18 under the SSJ calibration then live: `−11.4073%` and `−1.9742%` of SS (not re-verified 2026-10-26). [corrected 2026-10-26]: was "must both fall (positive = timing bug)". Also check the GK FOC table `main.py` now prints — all four legs must read `OK`.
5. Confirm the IC-δ consistency check and Walras residuals printed by `main.py` are unchanged.
6. Update the living docs after any calibration or structural change — **STATE.md and PROGRESS.md (changelog entry)**, plus TODO.md for open work (not just CLAUDE.md). This is **enforced** by two hooks that block the commit otherwise:
   - `.claude/hooks/require-docs-before-commit.sh` — PreToolUse gate, fires when Claude Code runs the commit.
   - `.githooks/pre-commit` — git-native twin, covers terminal commits. **Enable once per clone: `git config core.hooksPath .githooks`.**

   Both fire only when the commit stages any `*.py` (the old `code/**` pattern is now inert); doc-only commits pass. Bypass a false positive with `git commit --no-verify`. Keep the required-doc set in the two files in sync. (`docs/PROCESS.md` was retired 2026-07-30, superseded by PROGRESS.md.)
7. Commit the changed `.py` files, with the doc updates in the same commit.
