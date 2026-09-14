# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Two-country heterogeneous-agent model of a monetary union with Gertler-Karadi
financial intermediaries and sovereign default risk, calibrated to the
2010–2012 Greek sovereign debt crisis. The default mechanism follows
**Bocola (2016, JPE) "The Pass-Through of Sovereign Risk"**: an EXOGENOUS
rise in the *priced* probability of default π_t (his s-shock, eqs. 11–12 —
an input path, never a function of debt) makes bond prices fall, banks take
mark-to-market losses, the single-λ occasionally-binding incentive
constraint tightens, lending spreads rise and output falls — with no default
ever realized. Only D is default-risky; F bonds are safe. Application: ECB
asset purchases (TPI). Primary output is a research paper (Overleaf:
https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c).
(The 2026-07-16 Bocola-faithful rewrite replaced the earlier Cole-Kehoe
crisis-zone wrapper, the always-binding IC, and the patched default branch
*within the global-projection pipeline*; see git history on branch
`bocola-rewrite`.)

**The repository carries TWO solvers of this model.** The sequence-space (SSJ)
pipeline in `code/*.py` is the linearised sticky-price / nominal-deposit model
and is what `experiments/` (E1–E4), `diagnostics/regimes/` and the paper's
current figures run on. The global-projection pipeline in `code/global/` is the
nonlinear occasionally-binding solution. They share no code and no interpreter;
neither supersedes the other. Say which one you mean before quoting a number:
the two are calibrated differently and their impulse magnitudes are not
comparable.

## Environment

Two solvers now live in this repository, at non-overlapping paths, and they
need DIFFERENT interpreters. Pick the one that matches the code you are editing.

**A. Sequence-space (SSJ) pipeline — `code/*.py`.** Always use
`/opt/anaconda3/envs/ssj/bin/python`. The base Anaconda environment has a broken
`liblapack` symlink that causes silent numerical failures.

```bash
conda activate ssj
/opt/anaconda3/envs/ssj/bin/python code/main.py
```

**B. Global projection pipeline — `code/global/`.** Plain `python3`
(numpy/scipy/matplotlib). It does NOT use `sequence_jacobian` and must not be
run under the `ssj` environment.

```bash
cd code/global && python3 main.py
```

## Model code B — global projection (`code/global/`)

**Package layout.** The modules are grouped into subpackages; imports are
absolute from the `code/global/` root (`from blocks.bank import …`). `main.py`
sits at the root (run `python3 main.py`) and is the CHEBYSHEV-SMOLYAK PROJECTION
driver — there is NO perfect-foresight / representative-branch machinery
(`solver_pf/` was deleted 2026-08-11; git history preserves it).
- `main.py` — projection driver: SS → TFP → risk pass-through → OMT/TPI (all recursive)
- `config/` — `calibration.py`, `steady_state.py`
- `blocks/` — economic blocks (solver-agnostic): `bank.py`, `government.py`,
  `household.py`, `distribution.py`, `rouwenhorst.py`, `fast_kernels.py`,
  `firms.py`, `capital.py`, `trade.py`
- `solver_recursive/` — the ONLY solver: recursive global solution by GLOBAL
  CHEBYSHEV COLLOCATION (Bocola's own design): `state_grid.py`,
  `decision_rules.py`, `point_map.py`, `collocation.py` (the Newton),
  `recursive_main.py` (time iteration, now only a warm start),
  `recursive_experiment.py` (risk + TFP), `ltro_experiment.py` (the LTRO backstop)
- `reporting/` — `prints.py` (SS table), `plots.py` (activation-IRF figure)
- `tests/` — regression suite

The model is solved GLOBALLY as recursive decision rules on a Smolyak sparse
grid (Chebyshev interpolation), over the 10-state vector
`[K_D, K_F, P_D, P_F, b_DD, b_DF, b_FD, V_dep, s, Z_D]` — two capital
stocks, two banks' gross deposit obligations, the three carried sovereign
holdings, the cross-border deposit position, the sovereign-risk factor s, and the
TFP state Z_D (deterministic AR(1); the TFP experiment reads the IRF along a
Z-decay path). The CB backstop adds no state. At each grid point
THIRTEEN unknowns are solved (the per-period image of the old stacked system)
with Bocola's closed-form occasionally-binding μ. Expectations are genuine
multi-branch Gauss-Hermite quadrature over the s-innovation × a COMPOUND regime
`(default d′, CB-active m′)` — see `decision_rules.regime_table`.

**Driver: GLOBAL COLLOCATION NEWTON (`solver_recursive/collocation.py`), 2026-08-28.**
The policy VALUES at the collocation points are the unknowns and there is no inner
root find — Bocola's `residual_model.m` + `parsolve.m` exactly. Every stored rule is
an unknown (19 per point per regime: the 13 market-clearing/Euler unknowns plus the
six objects that used to be READ OFF a frozen continuation — alpha, C, r_wc per
country — which now carry Bocola's identity residual `log(guess/implied)`). The whole
coefficient vector goes to one damped Newton with a finite-difference Jacobian
(`parsolve`, dense) or Newton-Krylov (`krylov`, Jacobian-free) on the refined grid.
Solve ladder, also his: coarse μ=1 grid → d=0 at π=0 → d=1 by haircut homotopy
(0.85/0.70/0.55/0.45) → joint → SEED the s-refined grid and re-solve there.
Time iteration (`recursive_main.time_iteration`) survives ONLY as the warm start that
puts the Newton inside its basin — it is not a convergent solver here: its binding
mode is the franchise-value recursion at 0.990 per sweep, so runs reported
`max|F| = 1e-14` and "rule-change tol not reached" simultaneously.

**Grid: μ=1 Smolyak × a DENSE Chebyshev factor in s** (`SmolyakGrid(refine=(dim, m))`).
Raising the Smolyak level instead raises the GLOBAL budget; the tensor factor buys
degree m−1 in the one dimension that carries curvature (the logistic p^d(s)) and full
interaction with the sparse basis. Measured relative RMS error on this model's
curvature profile: μ=1 21pts **1.9e-1**, μ=2 221pts **3.9e-2**, m=5 95pts **2.5e-2**,
m=9 171pts **1.1e-3**. `S_REFINE = 5` ships (95 points, ~70 min); `S_REFINE = 9` is
Bocola's own resolution and the ladder walks 5 → 9, at ~4 h, because the dense Jacobian
is m+1 = 19·2·n+1 residual evaluations and the solve scales as n².
**CONVERGENCE CHECKED 2026-08-29** at the 100 bp calibration: going 5 → 9 moves the
impact output response from −0.1105% to −0.1087% (fitted) and −0.1278% to −0.1234%
(exact) — 1.6% and 3.4%, both well inside the 13%-wide identification bracket — and
every other reported number in the third digit (credit spread +90.6 → +90.8 bp/yr,
Q_bD −9.263 → −9.304%, Euler ALL −4.37 → −4.38). **The solution is converged at 5 for
every reported object**; 9 is the confirmation, not the working setting. NB this was
NOT true before the recalibration, when refining was fighting a Gibbs phenomenon at
the KKT kink. Hot kernels (household EGM backward, distribution
forward) are numba-JITed with an exact pure-numpy fallback (`cal["use_numba"]`).

| File | Contents |
|------|----------|
| `calibration.py` | All parameters. Single λ per bank (Bocola IC); Bocola/Greece anchors documented inline. Credit spread 100 bp/yr (NOT his 8 — see the kink note above); leverage 5, exposure 7.6%, recovery 0.45 are his. f = exit/payout share; Ω = β·[f + (1−f)α′] (Bocola's ψ = 1−f survival weight on the franchise value). |
| `steady_state.py` | Two-stage SS solve: {rk_D, rk_F, p} on capital markets + current account, then {β_D, β_F} on deposit markets. Symmetric SS required (see docstring). |
| `bank.py` | GK/Bocola bank block. `bank_backward` (α, μ, bond prices, cross-border FOC holdings), `bank_forward` (net worth, dividends, deposit supply; portfolio shares on ACTUAL net worth). PRICED (`def_price_D`) vs REALIZED (`def_real_D`) default split; only D is risky, F bonds are safe. |
| `government.py` | HM perpetuity bonds, Bohn rule. `govt_transition` forward-integrates the debt stock in one pass. Default risk is exogenous (no crisis zones). |
| `solver_recursive/point_map.py` | The per-point period map (image of the old stacked system): 13 residuals at one grid point given the frozen continuation rules. Bocola closed-form μ, with the LTRO facility entering it as `(n+m)/(lev-λm)`; quadrature over the s-innovation × the compound regime table; Z_D read from the state. |
| `solver_recursive/state_grid.py` | Smolyak sparse grid + Chebyshev basis, with `refine=(dim, m)` for a dense tensor factor on one dimension; `build_state_box`, `default_prob`, `s_process_params`. |
| `solver_recursive/collocation.py` | THE SOLVER. `make_residual` (the global F(theta), image of `residual_model.m`), `parsolve` (port of his damped FD Newton), `krylov_solve`, `solve_collocation`. |
| `solver_recursive/recursive_main.py`, `recursive_experiment.py`, `ltro_experiment.py` | Time iteration (warm start only) + SS anchors; the risk + TFP experiments and the solve ladder; the LTRO-backstop activation comparison (E1 never-fired path, E2 bond decomposition, E3 franchise-value counter-test). |
| `fast_kernels.py` | numba kernels for EGM backward + distribution forward; exact numpy fallback when numba is absent (`cal["use_numba"]`). |
| `household.py`, `distribution.py` | EGM with GHH utility; stationary distribution and forward iteration. |
| `trade.py` | CES basket and bilateral flows with PER-COUNTRY home bias and the country-mass ratio (`size_ratio`); `omega_home_F` is derived from `omega_home_D` and the sizes so trade balances at p = 1. |
| `firms.py`, `capital.py` | Flexible-price production with the Neumeyer-Perri working-capital wedge (w ÷ (1+ζ·r_wc), the spread→output channel; ζ=0 nests exactly — Bocola §V.C's own open-economy fix), Jermann adjustment costs, CES/Armington trade. |
| `prints.py` | Console reporting: `banner`, `print_ss_table`, and THE UNIT CONVENTION (`bp_ann`, `ann_pct`, `ann_prob`, and the `BOCOLA_IRF_*` benchmarks). Rates are annualised bp; p^d is printed quarterly AND annual; flow responses in level % with a ×4 annualised companion — Bocola's Table 5 unit. |
| `plots.py` | `plot_activation_irf` (the OMT/TPI activation overlay), written to `output/`. |
| `main.py` | Projection driver: SS → TFP → risk pass-through → OMT/TPI, each a full time-iteration solve. Heavy by design (~20–30 min). |
| `tests/` | Regression suite (see below). |

## Running and testing

### A. Sequence-space pipeline (`code/*.py`)

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

### B. Global projection pipeline (`code/global/`)

```bash
cd code/global
python3 main.py                                       # full projection pipeline (SS+TFP+risk+TPI), ~20-30 min
python3 -m solver_recursive.recursive_experiment      # risk pass-through only
python3 -m solver_recursive.ltro_experiment           # LTRO backstop only (phi = 0/50/100%)
python3 tests/test_ss_identities.py          # SS theory identities (fast)
python3 tests/test_bank_block.py             # bank FOC/no-arbitrage identities (fast)
python3 tests/test_fast_kernels.py           # numba/numpy kernel equivalence (fast)
python3 tests/test_state_grid.py             # Smolyak grid exactness (fast)
python3 tests/test_collocation.py            # THE SOLVER: packing, the six identity
                                             #   residuals, the refined grid, and a real
                                             #   d=0 solve to max|F| ~ 1e-9 (~90 s)
python3 tests/test_recursive_nesting.py      # SS rest point (N1) + the pi=0 grid-wide
                                             #   solve (N2, a hard gate since the
                                             #   collocation Newton replaced time iteration)
```

**Comment convention** (enforced across `code/global/`): every module and every
function carries exactly ONE leading ALL-CAPS comment saying what it is; any
further explanation is lowercase `#` comments attached to the specific hard
line. No docstrings, no bold markers, no prose blocks inside function bodies.
Console output lives in `prints.py`, never inside the model blocks.

**Acceptance thresholds** (all enforced in tests):
- Global collocation: Bocola's own test, `sum(F^2) <= m*(1e-9)^2` — the sum a
  uniform `max|F| = 1e-9` (`collocation.TOL_MAXF`) would give — at EVERY stage, over
  the 19 equations × points × regimes. This replaces the old two-part
  time-iteration test (settled rule AND every point clearing), which could pass on
  residuals while the rules were still moving. 1e-9 rather than machine zero because
  the period map's arithmetic floor is ~1e-10: the capital and bond Eulers difference
  O(1) expectations down to O(1e-4), and no Newton step improves on that. It is still
  four orders below any economic signal (the headline shock moves μ by 7.5e-3).
- goods_D (imposed) ≤ 1e−9; goods_F (Walras-redundant diagnostic) ≤ 2e−6 —
  including when the debt stock moves.  (The Newton solver typically lands
  goods_D near 1e−13; acceptance is `tol_transition` = 1e−10 normalized.)
- Zero-shock transition stays at SS to ≤ 1e−5.
- Risk-only shock (exogenous π): Q_bD↓, n_D↓, n_F↓, Y_D[0]↓, C_D[0]↓,
  lending spread↑, b_gov↑, Tax↑ (a positive Y or n response to sovereign
  risk = bug).
- Complementarity on every solved path: μ ≥ 0, slack = αn − λ·assets ≥ 0,
  μ·slack ≈ 0 (`out["mu_D/F"]`, `out["slack_D/F"]` from point_map.py).
  Known open item: risk-on n_D[0] can sit above
  risk-off (M1 deposit-rate channel; test warning, not assert) and
  post-impact Y_D runs mildly positive — both die with the union deposit
  market (docs/sunspot_transition_study.md §8).

## Key modelling choices — do not "fix" without checking docs/SPEC.md

- **Single λ (Bocola 2016 eq. 3):** all three asset classes carry the same
  divertability. Diverging them re-opens the portfolio-substitution margin
  that made sovereign risk *expansionary* pre-rework.
- **Priced vs realized default:** `def_price` enters bond pricing and
  expected-return FOCs; `def_real` enters realized returns and government
  flows. The baseline experiment prices risk but never realizes it
  (Bocola's pass-through design); a realized-default variant just passes
  `def_real ≠ 0`.
- **Endogenous debt in clearing:** the government's end-of-period stock is
  forward-integrated inside every residual evaluation and absorbed by banks.
  Clearing against a fixed `B_gov_ss` instead re-opens a Walras leak of
  ~0.5% of GDP per 5% debt deviation.
- **ASYMMETRIC COUNTRY SIZE, SYMMETRIC PER-CAPITA STEADY STATE (2026-08-28):**
  `size_F/size_D = 8`. Every variable is PER CAPITA of its own country and the
  per-capita SS is UNCHANGED (p_ss = 1, identical n_ss, leverage 5, μ_ss = 0.001,
  identical deposit supply); the mass ratio enters ONLY where D and F quantities
  are aggregated — goods market, union deposit clearing, both sovereign markets,
  the union wealth identity `W_F = P_F − V/(sz·p)`. Sovereign holdings are carried
  in the ISSUER's per-capita units, so `b_DD + b_DF = B_D` still clears the D
  market and the F bank's own book holds `b_DF/sz`. Home bias MUST scale with size
  or trade cannot balance: `(1−ω_F) = (1−ω_D)·size_D/size_F`, so D imports 15% of
  its basket and F imports 1.875% of its (`omega_home_F` is DERIVED in
  calibration.py). WHY: with a symmetric union D is half the union, so D's own
  sovereign shock moved the union real deposit rate 45 bp/yr and cancelled 78% of
  the credit-spread rise before it reached any firm's wage bill — the 2026-08-28
  audit's finding. Bocola's §V.C open economy has no such feedback: his
  `R = 1/β + 0.01·(B_for/gdp)` is a WORLD rate. `size_F = size_D` nests the old
  symmetric model exactly.
- **Symmetric steady state (in per-capita ratios):** other country asymmetries
  enter through shocks only.
  An asymmetric SS (e.g. δ_b_D ≠ δ_b_F) shifts p_ss off 1 and opens an
  O(1e−4) SS goods-market wedge (p is weakly identified by external balance
  at trade elasticity 0.5; see steady_state.py docstring).
- **Occasionally-binding IC (Bocola):** the leverage constraint enters the
  stacked system as the Fischer-Burmeister complementarity between μ (from
  the capital FOC, valid in both regimes) and slack = αn − λ·assets, scaled
  by μ_ss and n_ss (FB's zero set is scaling-invariant). At the SS the
  constraint binds (μ_ss ≈ 0.02, slack = 0), where FB is smooth. Portfolio
  shares and branch initial conditions divide by ACTUAL net worth, not n_IC.
- **Ω-kernel weights (Bocola Prop. 1):** f = exit/payout share, so
  Ω = β·[f + (1−f)·α′] — weight 1−f ≈ 0.95 on the franchise value α′
  (Bocola's survival ψ). beta_inter ≈ β_hh ≈ 0.99 proxies the household SDF;
  values ≪ 1/(1+rdep) drive α_ss below 1 and mute the franchise channel
  (the pre-rewrite code had the weights swapped AND beta_inter = 0.96).
- **Risk channel = genuine multi-branch quadrature (solver_recursive/), NOT a
  representative branch.** The default fork enters `point_map.py`'s banker FOCs
  as a real probability-weighted integral: Gauss-Hermite over the s-innovation ×
  the default realization d′∈{0,1} weighted by π_t (EXOGENOUS input path), where
  the default state is the SAME fitted decision rules evaluated at a reachable
  next-period point — never a frozen stand-in economy. The premium is endogenous
  (Ω^d > Ω^nd on the low default payoffs). `pi ≡ 0` nests the risk-neutral model
  exactly (test_recursive_nesting). The earlier perfect-foresight
  representative-branch pricing got the sign wrong (expansionary); the entire PF
  stack (`solver_pf/`: transition + solvers + risk_branch) was deleted 2026-08-11
  and the Chebyshev-Smolyak projection solver is now the ONLY machinery (TFP is a
  deterministic 7th state Z_D).
- **TPI = A STOCHASTIC LTRO BACKSTOP (2026-08-31), Bocola's own instrument.** With
  per-period probability `cal["phi_ltro"]` (a per-experiment scalar, NOT a state) the
  CB offers collateralised credit of size `cal["ltro_D"]`. It is his
  `residual_model_ltro_firstperiod.m` exactly: CB funding both LEAVES the divertable
  base and COUNTS as equity in the constraint,
  `mu_ratio = N'/(lambda*A')  ->  (N'+m)/(lambda*(A'-m))`. To first order that is
  `(1 + leverage) = 6x` the constraint relief of a bond purchase of the same size, and
  the numerator term is a margin NO quantity of bond-buying can reach.
  **IT IS A ONE-EQUATION CHANGE.** Lent at the deposit rate, the facility changes the
  COMPOSITION of the bank's funding, not its size or its cost: `P'` is algebraically
  unchanged, the household swaps one claim for another at the same rate so union
  clearing and `nfa` are unchanged, and the CB lends at the rate it pays so its carry is
  zero and NO remittance identity is needed. `test_recursive_nesting` N4 asserts exactly
  that — deposit clearing, `dep_D`, `P'`, `V'` and `n_D` bit-identical with the facility
  on, `mu` strictly lower. No new state, no new unknown, no complementarity.
  **FOUR regimes**, `(d,m)` orthogonal: the facility supports BANKS, so it is available
  in the default state too, and it has to be — the default branch carries little
  probability mass but the largest payoff deviation, so it dominates `cov(Om, payD)`,
  which is the term a credible backstop compresses.
  **SIZE IS THE CALIBRATION DECISION AND BOCOLA'S OWN IS A TRAP:** 2.0% of quarterly GDP
  unbinds the constraint at the SS and 3.4% unbinds it in the crisis state, against his
  40%. At his size `mu = 0` with huge margin in every relieved regime, so the whole m=1
  coefficient set sits ON the KKT kink. `ltro_D = 0.012` ships (halves the crisis
  multiplier, keeps `mu > 0` in both regimes).
  **THE HEADLINE READ IS THE NEVER-FIRED PATH** — regime `(0,0)`, announced and not
  drawn, which is the OMT fact. Two channels decide the sign and they oppose: the
  facility lowers `Om'` most where `payD` is lowest, shrinking `cov(Om, payD)` and
  raising the price everywhere (stabilising); but a looser future lowers `alpha'`, hence
  `E[Om]`, which RAISES today's `mu` (the charter-value channel, destabilising and NOT
  second-order). `ltro_experiment.run` reports both. **PREDECESSOR, RETIRED:** a
  one-sided yield peg with real purchases was built, solved and measured — purchases can
  only remove the LIQUIDITY premium (0.2-0.7% of the price here, 0.63% at the crisis
  corner against a 22.2% gap) because they work by pushing `mu` down and `mu` is floored
  at zero. `liquidity_ceiling_report` is that diagnostic, kept; see
  `docs/ltro_backstop_plan.md` and git history for the implementation.
- **Predetermined deposit rate:** the rate paid at t was locked at t−1
  throughout (bank funding legs, household EGM returns, μ timing).
- **Predetermined capital (Bocola eq. 6):** the stock producing at t was
  bought at t−1 (`Kap_prod[t] = Kap[t−1]`); mpk is the marginal product of
  the bank-held vintage, so impact output moves through hours alone. The
  old contemporaneous timing let the sovereign-risk investment boom raise
  Y_0 directly — reverting it re-opens the comovement problem.
- **Union deposit market (deposit-UIP):** deposits are own-good claims at
  national rates; a frictionless union interbank replaces the two national
  clearings with ONE union-wide clearing (D-good units) plus real-rate
  parity (1+rdep_D) = (1+rdep_F)·p′/p — the flexible-price image of one
  nominal union rate + national inflation differentials (BKK/Baxter-Crucini
  single-traded-bond margin). UIP makes the interbank pass-through
  zero-profit → no Walras leak; the cross-border deposit position
  (`out["nfa_dep_D"]`) is the absorption margin that broke the national
  S=I trap (the M1 comovement mechanism). A literal rdep_D=rdep_F with
  own-good legs is WRONG (unassigned RER valuation profit → Walras leak).
  Stage-2 SS imposes β_F = β_D (symmetric-SS doctrine).
- **Hatchondo-Martinez perpetuity:** stock decays at rate 1−δ_b; duration
  ≈ 1/δ_b quarters (0.036 ⇒ ~7y). Long duration is what makes priced risk
  generate large MTM losses — an interlude with δ_b=0.25/recovery=0.80 cut
  the repricing ~6x and made the risk channel expansionary (study §8).
- **Default branch = ONE pure-haircut feared event (Bocola):** the branch
  solves a single deterministic event — a full write-down to recovery
  `recovery_rate_D` = 0.45 (Greek PSI; Bocola's D = 0.55) on the whole
  claim, with the default-state recession arising endogenously through bank
  balance sheets. There are no scarring add-ons: the old Arellano output
  cost, GK ξ_K capital-quality loss and HFSF recap FLAGS were all 0 at the
  Bocola-pure baseline and were deleted in the 2026-07-21 cleanup (see git
  history if a variant needs them back). The recap machinery survives only
  as `_RECAP_LADDER`, a warm-start continuation used when the direct branch
  solve stalls. If the event is infeasible after that, the branch RAISES.
  Bohn taxes respond to the SURVIVING stock (taxing the pre-haircut stock at
  t=0 was a ~31%-of-GDP artifact).
- **Working capital (Neumeyer-Perri):** ζ_wc=1 × wage bill pre-financed at
  r_wc = rdep(−1) + λμ/Ω̃; the wedge is the only channel from spreads into
  impact output (without it Y_D moved −0.2% even at Q_bD −30%, n_D −20%).
  The LOAN is a bank asset inside the divertable base at the same λ, and the
  financing income accrues to BANK net worth through the deposit obligation
  `P' = R(QK'+qB'+L−N') − R_W·L` — Bocola's `residual_model_open.m` exactly.
  (`cal["wc_rebate"] = 1.0` instead hands it to households as a dividend, which
  makes the spread a pure intra-period transfer and, with no GHH wealth effect to
  offset it, turns the risk channel expansionary. Default is 0.)
- **Walras redundancy:** goods_F and the current account are *dropped* from
  the residual system and monitored as diagnostics.
- **Policy rules present:** the Bohn tax, and the LTRO backstop above
  (`phi_ltro`/`ltro_D`). No macroprudential policy, by design.

### A. Sequence-space architecture (`code/*.py`)

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
- **Hatchondo-Martinez perpetuity:** bond coupon decays at rate `1−delta_b`; duration ≈ 1/delta_b quarters. **Duration alone does NOT generate the MTM capital losses — duration *interacting with* `rho_s` (persistence of the latent default-risk factor) does.** Corrected 2026-08-28; the former wording ("This is what generates MTM capital losses") credited duration by itself and is wrong. Holding duration at the calibrated 17 quarters and varying persistence, the perpetuity and a one-quarter bill lose *exactly the same* on the reported impulse when `rho_s = 0` (both −1.01%), diverge mildly at `rho_s = 0.5` (−1.27% vs −1.01%), and separate only at the calibrated `rho_s = 0.95` (−5.58% vs −1.01%, a factor of 5.5). A transitory risk shock marks a long bond down no more than a short one. The loss is large because a persistent `s_t` keeps the hazard elevated across the claim's whole remaining life. Two corollaries: a short bond's *principal* is marked down too, so the distinction is never coupon-vs-principal; and any text describing the MTM channel must carry the persistence, not just the maturity.

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

### B. Global projection — known limitations

## Known limitations (documented, next thesis phases)

- Comovement problem RESOLVED (2026-07-18): predetermined capital + the
  union deposit market restored the impact contraction at the headline
  shock (Y_D[0] and I_D[0] both negative at π = 1%·0.95^t). Remaining
  next-phase dial: NK/union nominal block (needed for the TPI
  application); real interest parity currently plays the role of the
  single policy rate.
- Risk channel (recursive) approximations: Λ^nd ≡ beta_inter, rep-agent
  income-SDF proxy for Λ^d, household-side π-blindness (the deposit Euler never
  weights the default branch — faithful to Bocola, where household deposits are
  riskless too).
- **THE DETERMINISTIC SS IS NOT THE STOCHASTIC REST POINT** (2026-08-29). Verified:
  solved at pi == 0 the model sits on the deterministic SS for 200 quarters to six
  decimals with mu = 0.001001 and max|F| ~ 1e-7, so the SS and the solver are exact.
  But the risk-pricing rules rest elsewhere -- measured Y_D -0.13%, C_D -0.13%,
  I_D +0.23%, n_D +2.1%, K_D -0.58%, b_DD +2.6%, and **mu_D falls to EXACTLY 0**.
  Grid-independent (coarse vs s-refined agree to ~10% of the gap) and a UNIQUE GLOBAL
  ATTRACTOR (eight perturbed starts converge to the same state to 1e-12). Bocola has
  the same gap (his ergodic q = 0.979 against a deterministic 1.000, debt +2.0%) and
  locates it the same way: `recursive_experiment.stochastic_rest_point` is his
  `generate_irf.m` step 1, a zero-shock simulation to convergence. EVERY IRF starts
  there AND is differenced against an unshocked path (his
  `gdp = mean(gdp_s) - mean(gdp_nos)`). Before that fix the IRF charged the walk
  between the two rest points to the shock: 54-63% of the reported post-impact hump
  and of the +4.2% bank-net-worth overshoot was drift, and the "capital grinding
  down" was almost entirely drift.
- **THE ECONOMY RESTS ON THE KKT KINK, AND THAT IS THE BINDING ACCURACY LIMIT.**
  mu = max{.,0} is C0 and mu = 0 exactly at the rest point, so a Chebyshev interpolant
  returns mu > 0 in a neighbourhood where the truth is 0. Reading the fitted rules and
  clearing the period map exactly at the same state therefore disagree by MORE than the
  response: Y_D at p^d = 1.98% is -0.081% fitted against -0.008% cleared. The gap
  shrinks with resolution (0.087 -> 0.073 pp from 21 to 95 points) but slowly, as a
  Gibbs phenomenon does. THE SAME PATHOLOGY IS IN BOCOLA'S OWN SOLUTION: his fitted mu
  policy returns a 28.4 bp liquidity premium on impact where the exact multiplier gives
  2.1 bp, and 28.4 is his published number. `impact_table` and `dynamic_irf` print BOTH
  reads (`read_exact`); the pair is the honest object and the level of the output
  response is NOT identified at the current resolution.
  **CURED 2026-08-29** by moving `credit_spread_target` 8 -> **100 bp/yr**, which is
  where the constraint starts binding at the rest point (mu_rest 0 -> 0.0098) and the
  identification gap collapses 4x (0.087 -> 0.021 pp). f does NOT do this:
  `calibrate_bank_targets` forces alpha_ss = lambda*theta at the SS whatever f is, so
  the binding margin is ~theta*mu_ss and mu_ss is set by the SPREAD (measured: 8 bp ->
  mu_rest 0; 25 -> 0; 100 -> 0.0098; 250 -> 0.0286 with no further identification
  gain). 100 bp is Gertler-Kiyotaki's own target and inside the 100-300 bp periphery
  lending spreads of 2011-12; Bocola's 8 bp is his ESTIMATE, but it belongs to his
  closed model where mu contributes 2 bp to output and the wedge channel does not
  exist -- his SS V.C transmission and his closed-model mu^bg cannot both be imported.
- THE OUTPUT CHANNEL IS THE WORKING-CAPITAL WEDGE ALONE. GHH removes Bocola's
  closed-economy channel (the labour-supply wealth effect: in his benchmark
  `dlog l = −1.25·dlog c` exactly, and the leverage multiplier contributes 2 bp),
  so output moves only through `r_wc = rdep + λμ/E[Ω]`. BOTH legs matter, which
  is why country size is now asymmetric — see the key-choices list.
- THE BENCHMARK. Bocola's Table 5 (−1.05/−1.44/−1.53) is a cumulated quarterly
  GROWTH gap ×400 over an 8-quarter estimated shock sequence — its output LEVEL
  equivalent is −0.26/−0.36/−0.38%. The like-for-like single-shock IRF targets,
  rescaled to p^d = 1.98%/qtr, are **−0.157% (his §V.C open economy, whose GHH +
  working-capital transmission this model shares)** and −0.222% (his closed
  benchmark). `reporting/prints.py` carries these as constants and `dynamic_irf`
  prints them next to the trough.

**Nominal rigidities (added 2026-08-06 on `add-nkpc`; see `docs/SPEC.md` for the full
rationale and `docs/STATE.md` for numbers):**

- **Rotemberg price NKPC, subsidy-neutralised.** `pi = beta*pi(+1) + kappa_p*(mu_p*mc − 1)` in both countries; markup wedge `w = mu_p*mc*(1−alpha)*Y/N` in `labor_demand_D/F`. **Wages stay flexible** — `labor_market_D/F` is untouched. `mc_ss = 1/mu_p` neutralises the markup, so `mu_p*mc = 1`, `profit_ss = 0`, `pi_ss = 0` and **the steady state is bit-identical to the flex-price model**. `kappa_p → ∞` recovers flexible prices exactly — the standing equivalence gate.
- **Markup rent proportional to `e`, not lump-sum.** `firm_profit_D/F` route `(1 − mu_p*mc)(1−alpha)Y` through `income_D/F` in proportion to productivity, which makes labour-plus-profit income per unit of `e` exactly `(1−alpha)Y·e` — identical to the flex-price model, so the wedge bites only on hiring, never on household income. A lump-sum rebate was rejected as countercyclical (it would manufacture a progressive incidence result as an artifact of the rebate rule).
- **No policy rate.** The union-inflation normalisation `omega_pi_D*pi_D + (1−omega_pi_D)*pi_F = 0` (`phi_pi → ∞` limit of an ECB rule on union PPI) pins the inflation level; `p/p(-1) = (1+pi_F)/(1+pi_D)` pins the differential off the existing unknown `p`. No contract in the model carries a policy rate, so no Fisher relation is needed. `omega_pi_D = 0.071` is the renormalised capital key, **not** GDP weights — GDP weights would erase the 93/7 Greek-deflation split.
- **Nominal deposits against REAL sovereign bonds** — a deliberate asymmetry that maximises banks' Fisher exposure (nominal debtors, real creditors). **Must be stated as a modelling choice in the paper.** Nominal sovereign bonds are a candidate extension, not a correction.
- **27×27 solver system** (was 23×23): `+mc_D, pi_D, mc_F, pi_F` unknowns (and `rdep_D/F → i_dep_D/F`), `+nkpc_p_res_D/F, tot_res, union_pi_res` targets.

## Branch convention

- `bocola-rewrite` — the global-projection working branch (Bocola-faithful trim,
  2026-07-16); merged with `main` so both solvers coexist.
- `add-nkpc` — the sequence-space working branch. Sticky prices + nominal deposit
  contracts (Tasks 1–16, 2026-08-05/06).
- `main` — merge target. Contains the six structural fixes (W-1, W-2, W-3, T-2, A-2,
  TPI-1, PR #27), the modular-file reorganisation (PR #28), and the sticky-price /
  nominal-deposit merge.
- `global` — pre-rewrite snapshot (CK zones, always-binding IC).
- `audit` / `AB-audit` — historical audit branches. `AB-audit` was merged into `main`
  (PR #27); `audit` (PR #26) was closed as superseded. Do not reuse.
- `bank-cal` — old calibration branch predating structural fixes. **Do not merge.**
  Port calibration values only (see `docs/bank_cal_review.md`).

## History

The `code/global/` model was built in July 2026 as a standalone Bocola-faithful
replacement for the sequence-space code. **It did not replace it.** Both are live:
the sequence-space pipeline carries the sticky-price / nominal-deposit model, the
`experiments/` results set (E1–E4) and the paper's current figures; the global
pipeline carries the nonlinear occasionally-binding solution. The SSJ-era audit
trail (six structural fixes W-1…TPI-1, Walras forensics) lives in `docs/audit.md`,
`docs/walras_forensics.md` and `docs/verification_report.md`.

## Current model state and open issues (sequence-space pipeline)

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
