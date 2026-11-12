# STATE archive: sections 1-4 as they stood on 2026-11-12

> Moved verbatim from STATE.md when PR #38 (the appendix-B calibration in the global solver) was merged and its solution replaced the 5 Nov run of the old code calibration as the reference. Nothing below was edited. The 5 Nov output files are kept in results/GLOBAL_2026-11-05_code_calibration.

## 1. Calibration status, as of 2026-10-27: the paper is ahead of the code

The viva paper (sections 2 and 4, appendix B) and the code calibration
(`calibration/global_projection.py`) describe different economies. Until they agree, every
model number in the paper stays [TBC] and no code number goes into the paper.

| | Paper (Overleaf, 2026-10-22) | Code (every global result below) |
|---|---|---|
| Deposit market | national clearing | one union market, real deposit-rate parity |
| Leverage | 8 (Greek financial accounts, 7.98) | 5 (Bocola) |
| Intermediation wedge r^k - r | 5 bp/yr (mu_bar = 0.0010) | 100 bp/yr (mu_bar ≈ 0.012) |
| Own sovereign / bank assets | 14.4% (D), 6.5% (F) | 7.4% (D), 9.0% (F) (Bocola's Italian banks: 7.6%) |
| Foreign share of D debt | 12.7% (German banks) | 20.0% |
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

## 3. Reference solution: the 5 Nov 2026 run of corrected main

The run of `main` at 118908b on the code calibration (section 1), after the three
equilibrium-condition fixes of 2026-10-27 (section 4): real deposit parity, 115-point grid, a
2.00%/qtr headline shock decaying at rho_s = 0.95, and the 200 bp spread cap as the backstop.
Source: `results/GLOBAL/run.log`; the 6 Oct run it replaces is kept in
`results/GLOBAL_2026-10-06_reference` and its section here in
`Claude files/docs/shelved/STATE_archive_2026-11-05.md`. Spreads in bp/yr; output, consumption
and net worth in % of the economy's own rest point. "Lending spread" is the intermediation
wedge lambda_K mu_D / alpha_D (the excess return on capital the binding constraint requires).

| | No backstop | Cap 200 bp |
|---|---|---|
| Rest point: sovereign / lending spread | 56.4 / 98.6 | 50.8 / 63.5 |
| Rest point: mu_D (multiplier on the D bank constraint) | 0.0123 | 0.0079 |
| Rest point: output vs no backstop | 0 | +0.24% |
| Headline shock, impact: sovereign / lending spread | 298 / 178 | 200 / 20 |
| Headline shock, impact: D bank net worth, D bond price | -3.8%, -8.8% | -1.4%, -5.7% |
| Headline shock, quarter 8: D bank net worth, lending spread | +0.85%, 64 | -0.58%, 100 |
| Headline shock: peak purchases (% of the stock) | 0 | 36.4 |
| Headline shock: output on impact, trough | -0.059% (impact is the trough) | -0.130%, -0.299% (q8) |
| Headline shock: D consumption on impact, quarter 8 | +0.08%, -0.03% | -0.45%, -0.68% |

- **Effect of the 27 Oct fixes.** Impact output without the backstop falls by 0.059% against
  0.073% on 6 Oct (exact clearing: 0.059% against 0.069%); the terms-of-trade contribution
  rises from +0.011 to +0.018 pp. Every financial object and the whole backstop economy move
  only in the third digit.
- **Cap 200 bp: output falls more, cause not established.** The cap stabilises the bond and
  the banks on impact (lending spread 157 bp below the no-backstop value), yet output falls
  more on impact and for longer; exact clearing confirms the impact number. The responses are
  a redistribution from D to F: summed over 25 quarters, D consumption falls short of the
  no-backstop economy by 11.7% of quarterly D GDP while F consumption gains 10.6%. Only about a
  quarter of that is the Eurosystem's profit remitted to F (the bond's excess return over the D
  deposit rate on the 36% bought, 3.1% of quarterly D GDP after the 0.071 key, on the
  no-default trajectory). Two further candidates, neither decomposed: the backstop removes the
  fall in the D deposit rate that cushions the shock without it (61 bp higher on impact), and D
  banks lose the recapitalisation from the bond-price recovery (mu_D above its no-backstop
  value from quarter 4, lending spread +37 bp at quarter 8). Until decomposed, the backstop's
  output response is not classified (economic offset or accounting error) and is not quotable.
- **Price impact per euro is small.** At the rest point with rules fixed, buying 20% of the
  stock raises the D bond price 0.19% and lowers the sovereign spread 5 bp; at the shock,
  0.46% and 12 bp.
- **Accuracy (numerical).** Every reported collocation solve converges (no-backstop 115-point
  max error 3.7e-8; backstop 4.6e-11). The first step of the working-capital homotopy
  (no working capital, coarse grid) stops at 3.5e-3 in the F franchise-value identity; it only
  seeds the next step, which converges to 4.5e-9. No-backstop rules: mean log10 Euler error
  -4.92 on the ergodic set (max -2.50, D franchise value). Euler errors are not computed for
  the backstop rules. Complementarity: floor binds at 43 of 115 nodes, D banks hold no D bonds
  at 21; max |m x gap| 2.5e-9. Exact clearing still takes today's consumption and loan rate
  from the fitted rules.
- **Convergence in s (numerical).** Checked on 2026-08-29 at the 100 bp calibration, before the
  household fix: going from 5 to 9 nodes in the risk factor moves impact output by 1.6%
  (fitted) and 3.4% (exact), well inside the 13% identification bracket, and every other
  reported object only in the third digit. Not repeated on this run.

## 4. Standing facts about the global model

- Default risk is exogenous: the priced probability pi^d = 1/(1 + exp(-s)) depends on the
  latent factor s, never on debt. Only D bonds are risky.
- The default branch is solved on the same grid; regimes d' in {0, 1}. The LTRO backstop was
  deleted on 2026-10-01 (bit-for-bit neutral for the two-regime model).
- The backstop (TPI) adds the Eurosystem book as two states (M_cb, O_cb) and purchases m_cb as
  an unknown under a complementarity condition; 12 states, 14 unknowns. A full run with the
  backstop needs 4-5 hours of awake machine time on 8 cores.
- Three equilibrium conditions were brought into line with the paper on 2026-10-27 and are
  in the 5 Nov reference run: (i) each bank prices the foreign sovereign at its return in
  home goods, including the expected change in the terms of trade p_{t+1}/p_t (net worth
  already valued the realised legs that way, so the portfolio choice and the accounting were
  inconsistent); (ii) F taxes balance F's budget each period at the current F bond price with
  the stock held at B_F, closing the F budget leak that left the F goods market clearing only
  to about 2e-4 of F output on impact; (iii) the Bohn response gamma_tau is inverted from the
  debt root at the price banks pay rather than the risk-free price (a 0.8% fall in gamma at
  the Bocola calibration: the debt root moves from 0.9294 to its 0.93 target). None changes the
  deterministic steady state; all three move the dynamics, (i) through the terms of trade and
  (ii) through F household income. Their measured effect is in section 3.
- The banker's discount factor omits the ratio of consumer price indices, as the paper's
  (branch-sdf) does; recorded as a deliberate choice in the modelling-choices reference.
- Deliberate modelling choices and their reasons: `Claude files/docs/reference/global_modelling_choices.md`.
  Solver, grid and tolerances: `Claude files/docs/reference/global_solver.md`.

