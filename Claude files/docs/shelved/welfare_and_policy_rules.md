# Shelved: the welfare section and the policy-rule economies

Shelved on 2026-10-27. Nothing here is on `main`. The code is on the local branch
`shelved/policy-rules-welfare` (commit `bc1a83a`, the warranted floor and the policy-economy
driver, and `ca58548`, a work-in-progress snapshot of the welfare section). Solved economies are
in the git-ignored `results/WELFARE/` and `results/POLICY/`. Everything below was computed on the
code calibration (leverage 5, 100 bp wedge, union deposit clearing), not the paper's.

The long-term list in `TODO.md` carries the two ideas worth reviving: the level-versus-change
backstop with its welfare evaluation, and the warranted floor.

## From STATE.md as of 2026-10-27, verbatim

### 4. Welfare and the level-versus-change backstop (VIVA section 6, 2026-10-08)

Code calibration. Rule: holdings target theta = gamma_L (sp - 0.568) + gamma_D (sp - sp_lag),
a share of the steady-state stock, capped at 0.5, no sales, nothing bought in default;
13-state grid (125 points), every economy solved to max|F| <= 1e-9 at every node. Welfare is
the consumption equivalent (% of rest-point consumption) at the no-backstop rest point.
Headline numbers are the Monte Carlo with every quarter cleared exactly (40 x 600 quarters,
common shocks and defaults, collocation value at T with weight 0.17). Union is
population-weighted (D:F = 1:8).

| Rule | Union | D | F | Announcement U / D | Profit PV (% qGDP) | Buys / book |
|---|---|---|---|---|---|---|
| Level 0.025 | +0.0031 | -0.028 | +0.0070 | +0.0032 / -0.010 | 0.9 | 10.7% / 1.2% |
| Level 0.05 | +0.0046 | -0.047 | +0.0110 | +0.0049 / -0.016 | 1.6 | 11.4% / 2.2% |
| Level 0.10 | +0.0067 | -0.075 | +0.0168 | +0.0073 / -0.026 | 2.5 | 11.3% / 3.9% |
| Level 0.15 | +0.0080 | -0.097 | +0.0211 | +0.0090 / -0.033 | 3.1 | 10.9% / 5.3% |
| Level 0.20 | +0.0091 | -0.115 | +0.0246 | +0.0105 / -0.041 | 3.6 | 10.5% / 6.5% |
| Change 0.025 | +0.0026 | -0.029 | +0.0066 | +0.0028 / -0.012 | 0.5 | 8.0% / 0.6% |
| Change 0.05 | +0.0053 | -0.060 | +0.0135 | +0.0056 / -0.028 | 1.0 | 9.7% / 1.3% |
| Change 0.10 | +0.0105 | -0.127 | +0.0276 | +0.0116 / -0.067 | 1.8 | 10.8% / 2.5% |
| Combined (0.05, 0.025) | +0.0063 | -0.068 | +0.0156 | +0.0068 / -0.027 | 1.8 | 10.3% / 2.6% |
| Combined (0.05, 0.05) | +0.0082 | -0.092 | +0.0207 | +0.0090 / -0.041 | 2.0 | 9.8% / 3.0% |
| Combined (0.10, 0.025) | +0.0078 | -0.088 | +0.0199 | +0.0087 / -0.034 | 2.6 | 10.5% / 4.1% |
| Combined (0.10, 0.05) | +0.0095 | -0.106 | +0.0240 | +0.0106 / -0.046 | 2.6 | 9.9% / 4.4% |
| Cap 200 bp | +0.0095 | -0.127 | +0.0266 | +0.0102 / -0.055 | 1.3 | 1.5% / 4.6% |
| kappa_D = 1: level 0.05 | +0.0158 | -0.085 | +0.0284 | +0.0150 / -0.082 | 1.7 | 11.0% / 2.3% |
| kappa_D = 1: level 0.10 | +0.0187 | -0.098 | +0.0332 | +0.0184 / -0.104 | 2.8 | 11.1% / 4.2% |
| kappa_D = 1: level 0.15 | +0.0186 | -0.092 | +0.0324 | +0.0188 / -0.105 | 3.8 | 11.0% / 5.9% |
| kappa_D = 1: level 0.20 | +0.0171 | -0.076 | +0.0287 | +0.0176 / -0.094 | 4.4 | 11.0% / 7.4% |
| kappa_D = 1: change 0.025 | +0.0092 | -0.049 | +0.0165 | +0.0091 / -0.046 | 0.5 | 8.0% / 0.7% |
| kappa_D = 1: change 0.05 | +0.0166 | -0.098 | +0.0309 | +0.0168 / -0.101 | 1.0 | 9.7% / 1.3% |
| kappa_D = 1: change 0.10 | +0.0264 | -0.162 | +0.0500 | +0.0271 / -0.188 | 1.8 | 10.7% / 2.5% |
| kappa_D = 1: cap 200 bp | +0.0214 | -0.123 | +0.0394 | +0.0211 / -0.132 | 2.4 | 1.8% / 6.0% |

Rows with kappa_D = 1 remit the Eurosystem's profit and loss on D bonds to D. Sampling s.e.:
D 0.0006-0.009 (cap 0.012), F <= 0.001 (cap 0.0014); union and between-rule s.e. are
conservative bounds. Results: `results/WELFARE/<economy>/mc_exact.json`; tables refreshed
with `welfare_tables.py --refresh`, prose numbers with `welfare_numbers.py`, figure with
`welfare_figures.py` (all in `Claude files/global_projection/policy_rules/`).

- **Capital key.** The union gains 0.003-0.011% of consumption, all from the announcement
  (after-shock component -0.0014 to -0.0001); D loses 0.03-0.13%, F gains 0.007-0.028%. All
  rule classes lie on one frontier, with a union gain of 0.08-0.11 per point of D's loss,
  falling with strength. The union gain is still rising at the strongest rules solved: no
  interior optimum has been found.
- **kappa_D = 1.** The union gain is 2-3.5 times the capital-key value for the same rule; D
  still loses and F gains more. At equal D loss (0.098) the level rule (0.10) beats the change
  rule (0.05) by 0.0021 (4.0 conservative s.e.), three quarters of it in the announcement.
  The level rule saturates: the rest-point lending spread is 53.5 / 51.7 / 55.8 bp at
  gamma_L = 0.10 / 0.15 / 0.20 while D's crisis component keeps rising, so D's loss falls
  with gamma_L beyond 0.1 and gamma_L = 0.20 weakly Pareto-dominates 0.05. The change rule
  keeps relaxing the constraint (rest lending spread 5.5 bp at 0.10, mu 0.0007, almost slack;
  the least accurate economy).
- **Channels.** Under the capital key the backstop is a transfer problem (terms of trade, net
  exports); under kappa_D = 1 it is a terms-of-trade externality. The announcement works by
  raising the bank's expected marginal value of wealth (E[Omega_D] 1.159 to 1.164 under the
  capital key, 1.173 under kappa_D = 1) with net worth unchanged, so mu falls (0.0120 to
  0.0105 / 0.0066).
- **Households.** D's loss is even across income quintiles; F's gain rises with income.
- **Accuracy (numerical).** Grid-solved D welfare overstates the loss 1.5-2.3x under the
  capital key and understates it 6-24% under kappa_D = 1; fitted-rule simulation understates
  D's loss 1.8-3.7x under kappa_D = 1 (7.7x for the cap). Rank rules only on the exact Monte
  Carlo. Grid present values of profits and purchases are contaminated by unreachable
  negative-book nodes: use simulated values.

### 5. Policy economies in progress

The economies of VIVA section 5 were re-solved on the household-fix model from 2026-10-07 with
`Claude files/global_projection/policy_rules/run_policy.py`. Saved in `results/POLICY/<job>/`
(as of 2026-10-27): the warranted floor at tau = 16 and 0 bp, caps 300 and 360 bp, the 200 bp
cap with profits remitted to D, and kappa_nfa (premium on the cross-border position) = 0.0025
and 0.04. No saved economy for the tightest non-triggering cap or the 9-node refinement. All
on the code calibration, so none of them is a paper number.


## From PROGRESS.md as of 2026-10-27, verbatim

### 2026-10-08: Welfare and the level-versus-change backstop (VIVA section 6), uncommitted on `policy-rules`

- **The rule** (`solver_recursive/gamma_rule.py`): the Eurosystem's D-bond holdings target is
  theta = gamma_L (sp - sbar) + gamma_D (sp - sp_lag), a share of the SS stock (sp in pp/yr, sbar
  = the no-backstop rest spread 0.5678), smoothly capped at 0.5, no sales, nothing bought in
  default. The change term needs a 13th state, last quarter's spread, appended only on grids built
  with `sp_band` (`state_grid.ISP`, `STATE_NAMES_LAG`); 12-state grids, pickles and scripts are
  untouched (W1 below: the embedding is exact). `point_map`, `recursive_experiment.advance`,
  `stochastic_rest_point` and `accuracy` carry the state when present.
- **The solver** (`solver_recursive/rule_homotopy.py`): continuation in (gamma_L, gamma_D) from the
  embedded no-backstop economy with a chain-rule Jacobian through the interpolant (221 period-map
  calls per node instead of 5000 per column; equal to the finite difference to ~1e-6). Every
  economy solved to max|F| <= 1e-9 at every node, ~1-5 min each on a loaded machine. A node on the
  mu = 0 kink makes Newton converge at ~0.5 a step whatever the Jacobian; the chord step is kept
  while the residual falls by 20%.
- **Welfare** (`solver_recursive/welfare.py`, `welfare_households.py`): the linear Bellman operator
  on the collocation grid for any flow, the quiet history, regime-aware simulation (fitted or
  exactly cleared), and the incomplete-markets households re-solved under the solved aggregate law
  of motion. Drivers, all in `Claude files/global_projection/policy_rules/`: `run_welfare.py`
  (solve), `welfare_eval.py` (grid welfare), `quiet_exact.py`, `welfare_mc_exact.py` (headline),
  `welfare_sim.py`, `welfare_channels.py`, `welfare_hh_eval.py`, `welfare_tables.py`,
  `welfare_figures.py`. Results in `results/WELFARE/<economy>/`.
- **Tests:** `Claude files/global_projection/tests/test_welfare_rule.py` (W1 embedding exact, W2
  purchase algebra at every node, W3 the operator preserves constants, W4 chain-rule Jacobian =
  finite difference); the existing fast global tests pass.
- **Accuracy finding (numerical, not economic):** D's rest point is <1% of the bank's deposit
  obligation from the mu = 0 kink; the grid-solved welfare of D is inconsistent with the fitted
  dynamics by ~5e-4 in period utility per quarter, the size of D's welfare effects. Nine s-nodes do
  not help. The grid present values of Eurosystem P&L and purchases are contaminated by unreachable
  negative-book nodes. Headline welfare is therefore the exactly cleared Monte Carlo.
- **All 21 economies measured** (13 capital key, 8 with kappa_D = 1). New economic finding from the
  second batch: when D keeps the P&L the level rule's announcement SATURATES (rest lending spread
  52-56 bp for gamma_L 0.10-0.20) while its crisis insurance for D keeps growing, so D's loss falls
  with gamma_L beyond 0.1; the change rule keeps relaxing (rest lending spread 5.5 bp at gamma_D =
  0.10, constraint almost slack) and alone reaches union gains above ~0.019, at about twice the
  marginal cost to D. Corrected text ranges: union after-shock component -0.0014 to -0.0001 (was
  stated as -0.0011 to -0.0002); D's after-shock component under the capital key down to -0.074 (was
  -0.06); D's sampling s.e. up to 0.009 (was "below 0.006"). `welfare_numbers.py` now prints
  conservative union and between-rule standard errors (per-history values are kept only for the
  no-backstop economy, so the D-F and between-rule covariances are bounded, not estimated).
- **Paper:** in `06-welfare.tex` the profit-allocation subsection now precedes "Level against
  change", which it explains; frontier figure labels separated.
- Results: see STATE.md, "Welfare and the level-versus-change backstop".

### 2026-10-07: Parity kept; the warranted floor and a policy-rule driver committed (`policy-rules`)

- **Decision (user):** real deposit parity stays the global model's baseline. The union-real-rate
  shortcut (D pays F's real rate plus the debt-elastic premium, the expected terms-of-trade change
  dropped) is shelved as `git stash` "SHELVED 2026-10-07 (world-deposit-rate)". It leaves the
  valuation gain on the cross-border deposit position unassigned (the F goods market then fails to
  clear by ~1.4e-4 of GDP in the crisis state, ten times parity's ~1e-5), its full-grid solve
  stopped at max|F| = 4.8e-3, and it overshoots (impact output -1.1% coarse, -2.6% full). Its
  unconverged results are kept in `results/GLOBAL_union_unconverged_2026-10-06`; `results/GLOBAL`
  again holds the converged 6 Oct household-fix run. The nominal block (`nominal_block_scope.md`)
  is shelved too.
- **The warranted floor** (`solver_recursive/warranted.py`): the cap moves with the expected-loss
  spread, cap(s) = EL(s) + tau (+ a uniform shift used only as a homotopy), where EL(s) is the
  spread of the price an unconstrained risk-neutral investor discounting at the riskless F yield
  would pay. Written 5-6 Oct on top of 78f53c5 and never committed until now; `_tpi_floor` now
  takes s. Constant caps are untouched (`cap_bp` returns `tpi_cap_bp` unless
  `tpi_rule == "warranted"`).
- **`Claude files/global_projection/policy_rules/run_policy.py`** solves a policy economy the way
  `run.py` solves the 200 bp cap (coarse baseline, the policy walked in on the coarse grid, the
  s-refined grid 3 -> 5 with the smoothing eased in, the GZ hand-back) and evaluates it on the
  exact-clearing basis VIVA section 5 uses. Results go to `results/POLICY/<job>/`. The pre-fix
  drivers are preserved in `Claude files/global_projection/warranted_floor/`.
- Calibration unchanged: the 6 Oct run uses Daniel's pre-fix calibration (nothing in
  `calibration/`, `run.py` or the steady state changed between 78f53c5 and 583caa2).

