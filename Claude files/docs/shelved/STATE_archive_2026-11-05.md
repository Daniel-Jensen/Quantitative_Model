# STATE archive: section 3 as it stood on 2026-11-05

> Moved verbatim from STATE.md when the 5 Nov 2026 run of corrected main replaced the 6 Oct run as the reference global solution. Nothing below was edited. The 6 Oct output files are kept in results/GLOBAL_2026-10-06_reference.

## 3. Reference solution: the 6 Oct 2026 household-fix run

The run of `main` after the household budget fix (PR #36), on the code calibration (section 1)
and before the three equilibrium-condition fixes of 2026-10-27 (section 4):
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

