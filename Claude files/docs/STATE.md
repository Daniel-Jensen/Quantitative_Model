# Project State

Current truth only, global model first. Every number names its solver, its calibration and the
run that produced it. What is superseded is moved, verbatim, to
`Claude files/docs/shelved/STATE_archive_2026-10-27.md` (the full history to 27 Oct 2026);
the timeline is `PROGRESS.md`; open work is `TODO.md` at the repository root; symbols are in
`GLOSSARY.md`; data provenance is in `DATA.md`.

Rules for this file: replace a section when it is superseded and move the old text to the
archive; never stack dated notes on top of each other; keep it under about 300 lines.

---

# Part I: Global projection model (the paper's model)

## 1. Calibration status, as of 2026-10-27: the paper is ahead of the code

The viva paper (sections 2 and 4, appendix B) and the code calibration
(`calibration/global_projection.py`) describe different economies. Until they agree, every
model number in the paper stays [TBC] and no code number goes into the paper.

| | Paper (Overleaf, 2026-10-22) | Code (every global result below) |
|---|---|---|
| Deposit market | national clearing | one union market, real deposit-rate parity |
| Leverage | 8 (Greek financial accounts, 7.98) | 5 (Bocola) |
| Intermediation wedge r^k - r | 5 bp/yr (mu_bar = 0.0010) | 100 bp/yr (mu_bar ≈ 0.012) |
| Own sovereign / bank assets | 14.4% (D), 6.5% (F) | 7.6% (Bocola's Italian banks) |
| Foreign share of D debt | 12.7% (German banks) | 19.6% |
| Relative size F/D | 11.3 | 8 |
| sigma_s (risk-factor innovation s.d.) | 0.41 | 0.4455 |
| rho_s (risk-factor persistence) | 0.9 | 0.95 |
| Households, trade, fiscal | asymmetric (sigma 1.35, Frisch 0.5, home bias 0.68, eta 1.5, G 0.204/0.191, delta_b 0.056/0.039, Bohn 0.04) | symmetric per capita (log utility, Frisch 2, home bias 0.85, eta 0.5, G = 0, delta_b 0.056, debt root 0.93) |

Unreconciled wording: the 6 Oct results call the code calibration "Daniel's calibration",
while the 14 Oct calibration merge records Daniel's values as the paper's. The paper's values
are the ones agreed with the coauthor.

## 2. Open numerical risk under the paper's wedge

At 100 bp the incentive constraint binds at the stochastic rest point and the fitted and
exactly cleared output responses agree closely. At 8 bp the rest point sat on the kink
(mu = 0) and the output response was not numerically identified. The paper's 5 bp at
leverage 8 gives mu_bar = 0.0010, the same order, so the problem may return: priced default
risk leaves banks richer at the rest point, and the cushion before the constraint goes slack
is only 0.0125% of assets. Unresolved until the paper calibration is solved and both
responses are compared.

## 3. Reference solution: the 6 Oct 2026 household-fix run

Code calibration (section 1), real deposit parity, 115-point grid, 2.00%/qtr headline shock,
every number cleared exactly each quarter (`results/POLICY/base/eval.json`; the solution is
`results/GLOBAL`). Spreads in bp/yr; output and net worth in % of the economy's own rest
point. Every global impulse-response level from before this run is void, because households
previously received a constant in place of the working-capital flow.

| | No backstop | Cap 200 bp |
|---|---|---|
| Rest point: sovereign / lending spread | 56.8 / 98.1 | 51.1 / 72.6 |
| Rest point: mu_D (multiplier on the D bank constraint) | 0.0120 | 0.0090 |
| Rest point: output, investment vs no backstop | 0 | +0.22%, +0.49% |
| Crisis: spread sustained by the promise / on impact | 299 / 299 | 224 / 200 |
| Crisis: purchases on impact (% of stock) | 0 | 37.3 |
| Crisis: lending spread, bank net worth on impact | 185, -3.9% | 39, -1.3% |
| Crisis: output on impact, trough | -0.07% (q0) | -0.12%, -0.30% (q9) |

- **No backstop.** Output is back at rest by q2 and +0.06% in q4-6. Consumption rises 0.07%
  on impact, because the deposit rate falls 48 bp and cushions about 60% of the
  lending-spread rise, then sits about 0.10% below rest from year three for a decade: F banks
  earn the recovery return on their 20% of D debt along the no-default history (2.4% of
  quarterly GDP over ten years). Output on impact decomposes into credit spread -0.18 pp,
  deposit rate +0.10, terms of trade +0.01.
- **Cap 200 bp: a transfer problem.** The promise alone holds the spread at 224 bp; the
  Eurosystem buys 37% of the stock on impact and nothing after. Bank net worth falls by a
  third as much, yet impact output falls more (-0.115% against -0.072%) and the downturn is
  longer (-0.22% at q40). The Eurosystem's excess return on its 37% goes 92.9% to F by the
  capital key (3.6% of quarterly GDP over ten years), financed by a D trade surplus of
  0.2-0.28% of GDP and a 1.3-1.8% real depreciation, which lowers the consumption wage and
  hours under GHH preferences. The external accounts close to 1.6e-7.
- **Price impact per euro is small.** At fixed rules, buying 30% of the stock lifts the D bond
  price 0.72% at the headline shock (0.28% at rest; -19 and -7 bp of spread), against +3.7%
  in the solved backstop economy, 3.1 points of it through the continuation value. A cap far
  below the market spread is met only at the corner (21 of 23 nodes at p^d 4.8%/qtr); at the
  headline shock the purchase is interior and the D bank keeps about 45% of the stock.
- **Accuracy (numerical).** No-backstop rules: mean log10 Euler error -4.95 (max -2.54).
  Backstop rules: -3.21 (max -1.14), the error concentrated in the D-bond first-order
  conditions and the floor condition because the floor kinks the price and purchase rules.
  Backstop numbers are therefore reported on exact clearing (current quarter exact,
  expectations fitted).
- **Convergence in s (numerical).** Checked on 2026-08-29 at the 100 bp calibration, before the
  household fix: going from 5 to 9 nodes in the risk factor moves impact output by 1.6%
  (fitted) and 3.4% (exact), well inside the 13% identification bracket, and every other
  reported object only in the third digit. The solution is converged at 5 nodes; 9 is the
  confirmation. This was not true before the recalibration, when refinement fought a Gibbs
  phenomenon at the kink.

## 4. Welfare and the level-versus-change backstop (VIVA section 6, 2026-10-08)

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

## 5. Policy economies in progress

The economies of VIVA section 5 were re-solved on the household-fix model from 2026-10-07 with
`Claude files/global_projection/policy_rules/run_policy.py`. Saved in `results/POLICY/<job>/`
(as of 2026-10-27): the warranted floor at tau = 16 and 0 bp, caps 300 and 360 bp, the 200 bp
cap with profits remitted to D, and kappa_nfa (premium on the cross-border position) = 0.0025
and 0.04. No saved economy for the tightest non-triggering cap or the 9-node refinement. All
on the code calibration, so none of them is a paper number.

## 6. Standing facts about the global model

- Default risk is exogenous: the priced probability pi^d = 1/(1 + exp(-s)) depends on the
  latent factor s, never on debt. Only D bonds are risky.
- The default branch is solved on the same grid; regimes d' in {0, 1}. The LTRO backstop was
  deleted on 2026-10-01 (bit-for-bit neutral for the two-regime model).
- The backstop (TPI) adds the Eurosystem book as two states (M_cb, O_cb) and purchases m_cb as
  an unknown under a complementarity condition; 12 states, 14 unknowns. A full run with the
  backstop needs 4-5 hours of awake machine time on 8 cores.
- Known leak: the F treasury is held fixed (a documented simplification), so the F goods
  market clears only up to about 2e-4 of F output on impact of the headline shock.
- Deliberate modelling choices and their reasons: `Claude files/docs/reference/global_modelling_choices.md`.
  Solver, grid and tolerances: `Claude files/docs/reference/global_solver.md`.

---

# Part II: Sequence-space (SSJ) model

First-order, sticky prices, nominal deposits; used for the linear policy experiments
(E1-E4). Not the model the viva paper describes, and its impulse magnitudes are not
comparable with Part I.

## 7. Live calibration (`calibration/ssj.py`, as of 2026-08-18; not re-verified since)

EBA 2011 moments with `BANK_SCOPE="broad"` (banks intermediate the whole capital stock),
sticky prices, nominal deposits and the Gertler-Karadi structural refactor. Sovereign pricing
comes from the state-contingent bond payoff inside the bank's portfolio condition; there is
no ad hoc spread wedge. `psi_lambda_B` (collateral-friction amplification dial) = 0,
`zeta_writeoff` (share of principal written off in default) = 1, `writeoff_enabled` = 0
(responses trace the no-default branch), all four divertibility ratios `Delta` = 0.20,
`recovery_rate` = 0.30, `phi_lamb` (Bohn coefficient) = 0.15, `mv_rule` = 0 (par debt rule).
Solved steady state: both bond yields 80.0 bp/yr, steady-state spread exactly 0,
`lambda_gk_D` (divertible fraction) = 2.1087, `Omega_D` = 10.3462. Derivation of the moments:
`Claude files/docs/eba_calibration.md`.

On a 1 pp default-probability shock, with nothing tuned to a moment: peak spread 205.9 bp,
output -1.97% on impact, consumption -2.51%, D bank net worth -11.4%, German yield -16.2 bp
(flight to quality), German banks shed D bonds (-2.1%) while D banks absorb them (+2.1%).
TPI loading 0.520 / 0.504 / 0.482 at gamma = 2 / 5 / 10, below one: the Eurosystem does not
over-compensate.

## 8. SSJ open issues

- **E1-E4 and every `experiments/paper/fig0*.png` are stale.** The steady state moved with the
  structural refactor (`q_b_D` 0.968941 to 0.974906). Rebuild in order:
  `diagnostics/regimes/regime_model.py --force`, `experiments/run_all.py`,
  `experiments/e4_distribution.py`, `experiments/paper_outputs.py`.
- **`psi_lambda_B` is to be deleted (decided 2026-08-24).** It is zero and will stay zero; the
  paper asserts that no balance-sheet friction depends on the default probability, while
  `collateral_quality_D` still contains the term at coefficient zero.
- **F-1 (fiscal rule).** `mv_rule = 1` with `phi_lamb = 0.15` is a hard break, not a mild
  near-unit root; porting empirical duration is a two-parameter move.
- **Investment amplitude.** The model produces the sign and ordering of the Greek
  investment-output asymmetry but about 60% of its relative amplitude (investment/output fall
  ratio 1.5 in the model against 2.5 in 2007Q4-2015Q3).

Resolved issue codes (C-1, S-1, GK-1, GK-2, EL-1, PT-1) are defined in `GLOSSARY.md` and
documented in the archive.
