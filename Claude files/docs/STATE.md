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

The run of `main` after the household budget fix (PR #36), on the code calibration (section 1):
real deposit parity, 115-point grid, a 2.00%/qtr headline shock decaying at rho_s = 0.95, and
the 200 bp spread cap as the backstop. Source: `results/GLOBAL/run.log`. Spreads in bp/yr;
output and net worth in % of the economy's own rest point. Every global impulse-response level
from before this run is void, because households previously received a constant in place of
the working-capital flow (an accounting error that leaked into the goods market).

| | No backstop | Cap 200 bp |
|---|---|---|
| Rest point: sovereign / lending spread | 56.9 / 99.0 | 51.0 / 63.9 |
| Rest point: mu_D (multiplier on the D bank constraint) | 0.0123 | 0.0080 |
| Rest point: output, investment vs no backstop | 0 | +0.25%, +0.59% |
| Headline shock, impact: sovereign / lending spread | 299 / 180 | 200 / 19 |
| Headline shock, impact: D bank net worth, D bond price | -3.9%, -8.8% | -1.4%, -5.7% |
| Headline shock: peak purchases (% of the stock) | 0 | 35.9 |
| Headline shock: output on impact, trough | -0.069% (impact is the trough) | -0.131%, -0.301% |

- **Announcement.** With nothing bought, the promise raises the bank's expected marginal value
  of wealth (E[Omega_D] 1.159 to 1.164), so mu_D falls by a third, the lending spread by
  35 bp and output rises 0.25%.
- **No backstop, impact.** The lending spread rises 80 bp but the deposit rate falls 48 bp,
  cancelling 60% of it; output on impact decomposes into credit spread -0.18 pp, deposit rate
  +0.10, terms of trade +0.01.
- **Cap 200 bp.** The cap halves the fall in bank net worth and nearly closes the lending
  spread, yet output falls more on impact and much more at the trough. The D bond price gain
  (+3.7% against no backstop) comes almost entirely through the continuation value (+3.1
  points: the floor under the bond's future resale price), not the risk premium. Why output
  falls more is not established on `main`; the shelved policy-rule work attributes it to the
  capital-key transfer of the Eurosystem's profits to F
  (`Claude files/docs/shelved/welfare_and_policy_rules.md`).
- **Price impact per euro is small.** At the rest point with rules fixed, buying 20% of the
  stock raises the D bond price 0.20% and lowers the spread 5 bp.
- **Accuracy (numerical).** No-backstop rules: mean log10 Euler error -4.95 on the ergodic set
  (max -2.54). The floor kinks the backstop economy's price and purchase rules, so backstop
  numbers should be read on exact clearing (current quarter exact, expectations fitted).
- **Convergence in s (numerical).** Checked on 2026-08-29 at the 100 bp calibration, before the
  household fix: going from 5 to 9 nodes in the risk factor moves impact output by 1.6%
  (fitted) and 3.4% (exact), well inside the 13% identification bracket, and every other
  reported object only in the third digit.

## 4. Standing facts about the global model

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
