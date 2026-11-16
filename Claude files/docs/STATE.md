# Project State

Current truth only, global model first. Every number names its solver, its calibration and the
run that produced it. What is superseded is moved, verbatim, to `Claude files/docs/shelved/`:
`STATE_archive_2026-10-27.md` (the full history to 27 Oct 2026), `STATE_archive_2026-11-05.md`,
`STATE_archive_2026-11-12.md` (the old code calibration's last sections),
`STATE_archive_2026-11-16.md` (the Bohn 0.25 reference and the surprise cap) and
`STATE_calibration_branch_2026-11-09.md` (the cumulative STATE Daniel's PR #38 carried). The
timeline is `PROGRESS.md`; open work is `TODO.md` at the repository root; symbols are in
`GLOSSARY.md`; data provenance is in `DATA.md`.

Rules for this file: replace a section when it is superseded and move the old text to the
archive; never stack dated notes on top of each other; keep it under about 300 lines.

---

# Part I: Global projection model (the paper's model)

## 1. Calibration status, as of 2026-11-16: appendix B with Bohn 0.15 and a 2% headline shock

PR #38 put the paper's appendix-B calibration into `calibration/global_projection.py` and
calibrates the bank block at the stochastic rest point. On 2026-11-13 the authors set the Bohn
coefficient to 0.15 and the headline shock to p^d (quarterly priced default probability) = 2%,
keeping 1.5% as the lower-stress case (`run.py` RISK_SHOCK_PD). The paper and the code still
differ in the rows marked below.

| | Paper (Overleaf, 2026-10-22) | Code (every global result below) |
|---|---|---|
| Deposit market | national clearing | **differs:** national in the steady state (p_ss = 0.9893); one union market with deposit-UIP and a debt-elastic premium (kappa_nfa 0.01) in the dynamics. Left as is for the viva draft (authors, 2026-11-16) |
| Bohn coefficient (annual surplus response to debt, both over GDP) | 0.04 | **differs: 0.15**, the Greek 2010-19 estimate (s.e. 0.02); debt root 0.9635, half-life about 19 quarters |
| Leverage | 8 | instruments 8.109 (D), 7.969 (F); **7.95 at the rest point** |
| Intermediation wedge r^k - r | 5 bp/yr | 5 bp at the deterministic SS; **0 at the rest point** (section 2) |
| Own sovereign / bank assets | 14.4% (D), 6.5% (F) | **13.2% (D)** and 6.5% (F) at the rest point |
| Reserves (the D banks' claim on the Eurosystem) | divertible (in the asset base of the IC) | **differs:** not divertible (lambda_cb_D = 0) |
| Foreign share of D debt, F/D size, sigma_s, rho_s, pd_rest | 12.7%, 11.3, 0.41, 0.9, 0.10%/qtr | same |

Other code-only choices: a funded net-worth floor at 50% of n_ss that acts in default only; no
short sales on the cross-border books; the risk grid runs to p^d = 12%. The previous solution
(Bohn 0.25) is kept in `results/GLOBAL_2026-11-12_pr38` (main repo).

## 2. The D bank is unconstrained at the rest point (economic result)

At the rest point mu_D (the multiplier on the D bank's incentive constraint) is zero and D bank
net worth is 8.0% above its deterministic-steady-state value: the bank earns the default premium
in every history without default. The constraint binds from p^d = 1.31%/qtr (exact one-quarter
clearing; 1.30% with two quarters cleared jointly). The 2% headline shock is 3.2 unconditional
s.d. of the risk factor above its mean (7.4 innovations); the stationary distribution puts 0.07%
of quarters above it (1.5%: 2.9 s.d., 0.19%).

## 3. Reference solution: the canonical run of 13 Nov 2026 and the announced k = 2 rule

Outputs in `results/GLOBAL_2026-11-13_bohn015` (main repo): `no_policy/` (the pipeline run,
1.2e-9 at the nodes), `rules_bft_k2.pkl` (the rule, 1.3e-8), `headline.json`, `forensic/`.
Scripts in `Claude files/diagnostics/bft_rule/bohn015/canonical/`; narrative in
`Claude files/diagnostics/bft_rule/FINDINGS.md`. The policy is the proportional purchase rule
M' q = k q_ss B_ss (ER - ER_ref), never selling, known to all agents, k = 2, ER_ref 15.45 bp/yr
(the no-policy rest point's expected excess return). Each economy is measured from its own rest
point; numbers clear the impact quarter exactly, two-quarter joint clearing in parentheses.

| Impact, 2% shock | No backstop | Rule | Effect |
|---|---|---|---|
| Lending spread (bp) | 115 (119) | 0 (8) | -115 (-110) |
| mu_D | 0.025 | 0 (0.002) | |
| Sovereign spread (bp, rise) | 207 | 194 | -13 (-11) |
| D bank net worth (%) | -12.5 | -9.6 | +2.9 (+2.5) |
| Investment / output (%) | -1.77 / -0.103 | -0.93 / -0.024 | +0.84 / +0.079 |
| Purchases (% of SS D debt) | | 1.42 (1.49) | |

- **1.5% shock:** the rule removes a 35 bp lending spread with purchases of 1.06%; output
  -0.055% to -0.020%.
- **Announcement effect at rest:** sovereign spread -1.7 bp, output +0.04%, investment +0.10%.
- **Mechanism (controlled counterfactuals at the common no-policy state, 2%):** the rule known
  but no purchase in the impact quarter gives 97% of the rise in the price of capital (+0.44% of
  +0.46%) and the whole cut in the lending spread; the same purchase as a surprise gives
  +0.015%. Q_K: expected resale value +0.28, constraint premium +0.29 (the slack constraint
  removing the no-policy spread), deposit rate -0.11. Net worth +3.6%: capital revaluation 3.0,
  sovereign gain 0.6.
- **Without a backstop:** sovereign spread 62 to 269 bp, net worth -12.5%; output decomposes
  into lending spread -0.087 pp, terms of trade -0.060, deposit rate +0.043 (the D deposit rate
  falls 57 bp). Excess return +381 bp: 267 risk premium, 115 constraint premium.
- **Surprise 200 bp cap:** holds at 1.5% with 17.1% purchases (the purchase lowers next
  quarter's spread and with it the fair spread, 207 to about 191 bp); fails at 2% for purchases
  up to the 28% book band. Dropped from the paper (2026-11-16).
- **Accuracy (numerical).** Identities exact at every reported state (net worth 2e-15, bond
  market 1e-11, mu = max(psi, 0) 3e-11). Off-grid, on an exactly cleared 1,200-quarter
  simulation: mean / max log10 Euler error -4.4 / -1.9 (no policy), -3.8 / -1.7 (rule). The
  fitted rules overstate next quarter's franchise value by up to 3.4%. The fitted simulation of
  the solver's own accuracy report turns the Eurosystem book negative, which switches the rule
  off: not an accuracy measure for any BFT solution. The policy rest point must be the exact
  fixed point from a zero book (the fitted iteration lands at a book of -3.4% of the debt).
- **Open:** whether the rule leaves the constraint exactly slack at 2% (threshold 2.39% one
  quarter, 1.92% two); k >= 4 fails at the 27.9% book band; D consumption rises on impact.

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
