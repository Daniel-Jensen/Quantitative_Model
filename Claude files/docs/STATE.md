# Project State

Current truth only, global model first. Every number names its solver, its calibration and the
run that produced it. What is superseded is moved, verbatim, to `Claude files/docs/shelved/`:
`STATE_archive_2026-10-27.md` (the full history to 27 Oct 2026), `STATE_archive_2026-11-05.md`,
`STATE_archive_2026-11-12.md` (the old code calibration's last sections) and
`STATE_calibration_branch_2026-11-09.md` (the cumulative STATE Daniel's PR #38 carried). The
timeline is `PROGRESS.md`; open work is `TODO.md` at the repository root; symbols are in
`GLOSSARY.md`; data provenance is in `DATA.md`.

Rules for this file: replace a section when it is superseded and move the old text to the
archive; never stack dated notes on top of each other; keep it under about 300 lines.

---

# Part I: Global projection model (the paper's model)

## 1. Calibration status, as of 2026-11-12: the code carries appendix B, with four departures

PR #38 (merged 2026-11-12) put the paper's appendix-B calibration into
`calibration/global_projection.py` and calibrates the bank block at the stochastic rest point,
not the deterministic steady state. The paper and the code still differ in the rows marked
below, so paper model numbers stay [TBC] until the two agree.

| | Paper (Overleaf, 2026-10-22) | Code (every global result below) |
|---|---|---|
| Deposit market | national clearing | **differs:** national in the steady state (zero cross-border deposits, p_ss = 0.9893); one union market with deposit-UIP and a debt-elastic premium (kappa_nfa 0.01) in the dynamics |
| Bohn coefficient (annual surplus response to debt, both over GDP) | 0.04 | **differs: 0.25** (gamma_tau 0.0625, debt root 0.937); a faster rule chosen so that impulse responses return within 20 quarters. Greek estimates: 0.06 (2001-19), 0.15 (2010-19) |
| Leverage | 8 | instruments 8.109 (D), 7.969 (F), set at the Bohn 0.06 rest point; **7.96 at today's rest point** |
| Intermediation wedge r^k - r | 5 bp/yr (mu_bar = 0.0010) | 5 bp at the deterministic SS (mu_D 0.00101); **0 at the rest point**, where the D constraint is slack (section 2); the F instrument sits at its 1 bp bound |
| Own sovereign / bank assets | 14.4% (D), 6.5% (F) | **13.0% (D)** and 6.5% (F) at the rest point; the D instrument (12.9% at the SS) was set at the Bohn 0.06 rest point |
| Foreign share of D debt | 12.7% | 12.7% at the SS (12.45% at rest); D banks hold 0.5% of the F stock (the paper is silent) |
| Relative size F/D, sigma_s, rho_s, pd_rest | 11.3, 0.41, 0.9, 0.10%/qtr | same |
| Households, trade, fiscal spending, maturities | sigma 1.35, Frisch 0.5, home bias 0.68, eta 1.5, G 0.204/0.191, delta_b 0.056/0.039 | same |

Code-only choices the paper does not yet describe: the D banks' claim on the Eurosystem is not
divertible (lambda_cb_D, divertible share of that claim, = 0, Gertler-Karadi reserves); a
funded net-worth floor at 50% of n_ss that acts in default only (D's top-up is new debt the
banks absorb, F's a lump-sum tax); no short sales on the cross-border bond books; the risk
grid runs to p^d (quarterly priced default probability) = 12%.

## 2. The D bank is unconstrained at the rest point (economic result)

At the deterministic steady state the D multiplier mu_D (on the D bank's incentive constraint)
is 0.0010 by construction. Evaluated with priced risk at the steady-state state it is 0.017;
at the rest point it is zero, and D bank net worth is 7.6% above its steady-state value. The
bank expects to earn the default premium in every history without default, and that wealth
slackens the constraint. Daniel's sweep (Bohn 0.06) found no wedge target up to 60 bp that
makes it bind at rest. The kink risk recorded here before (rest point on mu = 0) has resolved
in the slack direction: there is no financial friction at rest, and the 5 bp wedge describes
only the deterministic steady state.

The constraint binds from p^d of about 1.2-1.5%/qtr (exactly cleared impact: mu_D 0 at 1.10%,
0.030 at 1.98%). The model's own risk process rarely gets there: simulated p^d has mean 0.147%
and 95th percentile 0.443%, and the 2010-11 Greek average (9.1%) is a 4.9 s.d. event. Under
Bohn 0.06 the threshold was about 5%; how the fiscal rule moves it is not documented.

## 3. Reference solution: the 12 Nov 2026 run of main after PR #38

The run of `main` with PR #38 (model code identical to Daniel's 0c80116, which reports the same
numbers): appendix-B calibration (section 1), 75-point fine grid, risk box to p^d 12%, warm start
`data/global_warm_start.pkl` (the Bohn 0.25 rules), a 1.5%/qtr headline shock from the 0.10%
rest value decaying at rho_s = 0.9, and the backstop as a 200 bp cap held by unanticipated
purchases on the no-backstop rules, bought each quarter as needed and held to maturity.
Solved in 10 minutes from the warm start. Source: `results/GLOBAL/run.log`; the 5 Nov run on
the old code calibration is in `results/GLOBAL_2026-11-05_code_calibration`. Every number is
exactly cleared at the visited state. Spreads in bp/yr; output, consumption and net worth in %
of the no-shock path. "Lending spread" is lambda_K mu_D / alpha_D (lambda_K, divertible share of
capital; alpha_D, franchise value per unit of net worth).

| 1.5% shock | No backstop | Cap 200 bp |
|---|---|---|
| Purchases (% of the SS D debt) | 0 | 16.7 on impact, none later |
| Impact: sovereign / lending spread | 222 / 68 | 200 / 20 |
| Impact: mu_D | 0.0149 | 0.0042 |
| Impact: GDP, investment | -0.076%, -1.49% | -0.045%, -0.80% |
| Impact: D bank net worth, D bond price | -10.0%, -6.2% | -7.3%, -5.4% |
| D consumption: impact, q8 | +0.13%, -0.02% | -0.04%, -0.18% |
| GDP: q8, q24 | -0.041%, -0.028% | -0.042%, -0.036% |
| D bank net worth: q8, q12 | +2.6%, +2.9% | +0.3%, +0.6% |

- **The constraint binds on impact only.** The lending spread is zero from quarter 1; output
  persistence after that runs through capital (-0.11% at q8).
- **Output on impact decomposes** into the credit spread -0.053 pp, the terms of trade
  -0.051 pp and the deposit rate +0.028 pp. The D deposit rate falls 37 bp and cancels 55% of
  the credit spread in the loan rate.
- **D consumption rises on impact** (wrong sign by the project's diagnostic). Daniel classifies
  it as a general-equilibrium offset: households draw down deposits as the deposit rate falls,
  and smooth a temporary income loss. That deposit-rate fall comes from the union deposit
  market of the dynamics, which the paper does not have (section 1). Not yet resolved.
- **D bank net worth overshoots from q4** (economic result): responses condition on no default,
  so the holder of D bonds collects the default premium as the price recovers.
- **The D bond's expected excess return rises 270 bp on impact:** 201 bp risk premium
  (covariance of the banker's discount factor with the bond's return) and 70 bp liquidity
  premium (the constraint). The cap lowers it 88 bp (-50 liquidity, -39 risk).
- **The cap trades impact relief for later D losses.** Under the capital key (0.071) 93% of the
  Eurosystem's no-default profit on the bonds goes to F, and D consumption runs below the
  no-backstop path from q1 (Daniel's report `Claude files/reports/tpi_root_cause_2026-11-05.md`
  section 10).
- **Price impact.** At the rest point a 20% purchase raises the D bond price 0.23% and lowers
  the spread 6 bp. At the shock a 5% purchase tightens the constraint (mu_D 0.0149 to 0.0160)
  before 20% switches it off (bond price +1.1%, spread -29 bp): the franchise-value channel.
- **A 2% shock** gives GDP -0.119%, lending spread +137 bp, sovereign spread 268 bp.
- **Rest-point targets missed after the Bohn change** (calibration): leverage 7.96 against 8,
  D wedge 0 against 5 bp, D exposure 13.0% against 14.4%; the transmission row (overdraft-spread
  rise at p^d 2%) is 142 bp against a 108 bp upper bound. Other validation rows: GR-DE spread at
  rest 61.5 bp (data 25.3), deposits 259% of annual GDP (84%), bank-held D debt 44% of GDP at
  market value (40%), duration of the D book 4.15 years (4.4).
- **Accuracy (numerical).** The collocation solves converge (joint 6.6e-9, hard-max polish
  1.4e-8). Mean log10 Euler error -4.37 on the ergodic set; worst -1.98, in the D
  franchise-value identity. The D banks' deposit-obligation state leaves the box in 1.5% of simulated periods.
  The fitted rules overstate the impact recession by a factor of two (GDP -0.160% against
  -0.076%, mu_D 0.031 against 0.015), so only exactly cleared numbers are quotable, and even
  these take today's franchise value, loan rate and banker consumption from the fitted rules
  (the inconsistency the 4 Nov audit found to halve mu_D at the old headline shock).

## 4. Standing facts about the global model

- Default risk is exogenous: the priced probability pi^d = 1/(1 + exp(-s)) depends on the
  latent factor s, never on debt. Only D bonds are risky. The default branch is solved on the
  same grid; regimes d' in {0, 1}.
- The backstop adds the Eurosystem book as two states (M_cb, O_cb) and purchases m_cb as an
  unknown. Three modes: "cap" (the reference, unanticipated), "relief" (buy at the shock until
  mu_D = 0) and "rule" (the cap anticipated as a rule). Under "cap" and "relief" no one prices
  the promise in advance; the anticipated rule has not solved below about 1100 bp on the 12% box.
- The three equilibrium conditions aligned with the paper on 2026-10-27 (terms of trade in the
  cross-border bond conditions, F budget balance at the current F bond price, the Bohn response
  set at the traded bond price) are in PR #38's code; the Bohn coefficient is now the
  calibrated primitive rather than the debt root.
- Equilibrium selection: at some normal-times grid points (p^d at most 1%) the period
  conditions also admit a self-fulfilling crisis solution (D net worth below half its SS value,
  first met at the capital-box edges). The solver selects the calm one. This bears on the paper's "no bad equilibrium" sentences.
- The multiplier switch mu = max(z, 0) and the net-worth floor are hard maxima in the final
  solve (smoothing walked down, then polished out); the collocation Jacobian uses the complex
  step. The rest point is a fixed point of the exactly cleared period map.
- The banker's discount factor omits the ratio of consumer price indices, as the paper's
  (branch-sdf) does; recorded as a deliberate choice in the modelling-choices reference.
- Deliberate modelling choices and their reasons: `Claude files/docs/reference/global_modelling_choices.md`.
  Solver: `Claude files/docs/reference/global_solver.md` (before PR #38) and Daniel's
  `Claude files/docs/guide_global.md`; how to run, resume and cache: `Claude files/docs/RUN_GLOBAL.md`.
  Daniel's analyses behind sections 1-3: `Claude files/reports/` (appendix-B validation,
  fix A, TPI mechanism, TPI root cause).

---

# Part II: Sequence-space (SSJ) model

First-order, sticky prices, nominal deposits; used for the linear policy experiments
(E1-E4). Not the model the viva paper describes, and its impulse magnitudes are not
comparable with Part I.

## 5. Live calibration (`calibration/ssj.py`, as of 2026-08-18; not re-verified since)

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

## 6. SSJ open issues

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
