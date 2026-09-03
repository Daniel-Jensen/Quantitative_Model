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
crisis-zone wrapper, the always-binding IC, and the patched default branch;
see git history on branch `bocola-rewrite`.)

## Environment

Plain `python3` (numpy/scipy/matplotlib). **Do not use the old
`/opt/anaconda3/envs/ssj` environment or the `sequence_jacobian` library** —
that was the previous implementation (see "History" below); the path no
longer exists.

## Model code (`code/global/`)

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

## Branch convention

- `bocola-rewrite` — current working branch (Bocola-faithful trim, 2026-07-16).
- `global` — pre-rewrite snapshot (CK zones, always-binding IC).
- `main` — merge target.
- `audit`, `bank-cal` — historical SSJ-era branches; do not use for new work.

## History

The previous implementation used the `sequence_jacobian` (SSJ) library
(`code/model_v12.ipynb`, `equations_*.py`, `audit_artifacts/`) — superseded
by the standalone `code/global/` model in July 2026. The SSJ-era audit trail
(six structural fixes W-1…TPI-1, Walras forensics) lives in `docs/audit.md`,
`docs/walras_forensics.md`, `docs/verification_report.md` and git history.
`docs/STATE.md` records the current model state and calibration.
