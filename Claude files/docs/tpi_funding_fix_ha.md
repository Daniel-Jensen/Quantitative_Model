# The TPI with TARGET2 funding and non-divertible Eurosystem claims, deposits cleared by heterogeneous households

*Results note, 14 October 2026. Global-projection model: branch `domestic-funding-fix` (`0dd0732`) plus the
switches in `Claude files/diagnostics/tpi_national/tpi_switches_v3.patch`, all off by default. This note is
complete: the full-grid solve of the TPI economy (run F2, 115 points, 7.9 hours) finished, and §2.4 and §3.2
report it.*

---

## 0. Results in brief

**Full grid (F2, §2.4 and §3.2) — the headline numbers.**
- The announcement is modest and robust: the sovereign spread falls from 56.5 to 50.7 bp, the lending spread
  from 98 to 65 bp, and output rises 0.27% (exact clearing; the fitted rules agree to 0.3 bp).
- On the headline shock the TPI holds the cap: the spread is 200 bp instead of 298, and the Eurosystem buys
  38% of the stock on impact. The banks are protected: D bank net worth falls 0.6% instead of 3.8%.
- But D output and consumption do worse with the TPI, from both starting points: output −0.14% on impact
  against −0.06%, and the 25-quarter consumption loss is −15 against −0.9.
- Heterogeneous households: with the TPI they cut consumption slightly less than the representative household
  (−0.57% on impact against −0.59%; 25-quarter loss −13.8 against −15.2), the reverse of the no-TPI economy.

**The quick grid (Q8) and the no-TPI household results:**

1. **The configuration solves and holds its 200 bp cap** (Q8).
2. **The announcement does most of the work, with nothing bought.** Under exact clearing at the rest point:
   - the sovereign spread falls from 62 to 24 bp;
   - the lending spread falls from 95 bp to 0, because the D banks' constraint goes slack;
   - output rises 0.9%, investment 1.8% and D bank net worth 1.9%.
3. **Non-divertibility changes quantities, not the spread.** Against the shipped TPI with divertible claims:
   - the rest-point spread is 24 against 26 bp;
   - output rises 0.87% against 0.53%, and bank net worth 1.9% against 0.3%;
   - purchases move the price 1.4–2.3 times more per euro;
   - the risk premium at the shock is 92 against 132 bp.
4. **The price impact per euro stays small**: −0.55 bp per 1% of the stock at the shock, against −15 to −30
   bp in the evidence.
5. **On the headline shock the spread on impact is 132 bp instead of 414, with nothing bought.** The quick
   grid's impulse responses are pathological (output −1.3% on impact), so the dynamics wait for F2.
6. **Heterogeneous households (no-TPI union economy, full grid)**, with D's deposit market cleared through them:
   - D consumption falls 0.15% on impact, where the representative household's rises 0.08%;
   - over 25 quarters the consumption loss is 2.8 times the representative household's (−2.5 against −0.9);
   - the poorest fifth cuts four times as much as the richest;
   - clearing moves D's deposit rate by less than 3 bp/yr.
7. **Clearing every price through the households along the path is not determined in the global solution.**
   The hybrid is unstable through the terms of trade. The fix inside the global framework is an aggregate
   consumption rule fitted to the heterogeneous households, Krusell–Smith style (§4).

---

## 1. What was run

The experiment changes two things in the shipped TPI and nothing else:

| Element | This run | Shipped TPI | First fix (previous note) |
|---|---|---|---|
| D banks' deposit funding | **union market**: a TARGET2 margin at the parity rate | union market | national |
| The Eurosystem's claim Z (what a D bank holds after selling bonds to the Eurosystem) | **non-divertible**, λ_Z = 0 (λ_Z: the share of Z a banker can divert) | divertible, λ_Z = λ | non-divertible |
| Spread cap | fixed 200 bp | fixed 200 bp | fixed 200 bp |
| Profit and loss | capital key, 7.1% to D | capital key | all to D |
| Default response to purchases | none | none | logit, φ = 1 |
| Seniority | pari passu | pari passu | pari passu |

**Funding fix.** Under the union market a D bank funds itself at the union deposit rate, with D's rate tied to
F's by deposit-UIP plus a small premium on D's cross-border position (κ_nfa = 0.01, the premium's slope). The
previous note reads this as the model's reduced form of TARGET2: the Eurosystem is the marginal funder of a
stressed periphery bank, as it was in Greece in 2011–12.

**Two layers of solution.**
1. **The global solution.** The Chebyshev–Smolyak rules for banks, firms, the government and the Eurosystem
   carry everyone's expectations, including the default risk. In this layer D has a representative household.
   - Quick grid, run Q8: 25 points, 2.5 hours, done.
   - Full grid, run F2: 115 points, 5 nodes in s (the sovereign-risk factor), still running.
2. **Heterogeneous households along the impulse path.** D's representative household is replaced by the
   steady state's heterogeneous households: idiosyncratic productivity, a borrowing limit, solved by EGM with
   a Young lottery for the distribution. Their discount factor is recalibrated so that, at the rest point,
   they hold exactly the representative household's deposits (β = 0.99435 against 0.99444 at the SS). Three
   versions, from least to most general equilibrium:
   - **Prices fixed (PE).** The households face the global solution's paths of every price and income.
   - **Deposit market cleared.** D's deposit rate is the price that clears. Under deposit-UIP, deposits the
     households hold beyond the global path's are lent abroad, and each unit of that position, relative to
     GDP, lowers D's rate by κ_nfa. Newton on the rate path; all other prices stay on the global path.
   - **Every price cleared along the path.** Each quarter the period map is re-cleared with D's deposit rate
     and the households' deposits given; an outer Newton matches both paths. This one is not determined (§4).

   F (eight times D's size) keeps its representative household, and the global rules' expectations stay the
   representative household's.

---

## 2. The global solution

### 2.1 The announcement: the rest point, nothing bought (quick grid)

Exact clearing; the fitted rules give the same spreads to within 2–13 bp.

| Rest point, against no TPI | No TPI | Shipped TPI (divertible Z) | **This run (non-divertible Z)** |
|---|---:|---:|---:|
| Sovereign spread | 62.4 bp | 26.2 bp | **24.4 bp** |
| Lending spread | 94.6 bp | 0 | **0** |
| μ_D (the D bank's constraint multiplier) | 0.012 | 0 (slack) | **0 (slack)** |
| D output | — | +0.53% | **+0.87%** |
| D investment | — | +0.90% | **+1.84%** |
| D consumption | — | +0.29% | **+0.23%** |
| D bank net worth | — | +0.30% | **+1.94%** |
| D deposit rate | — | −1.3 bp/yr | **−14.4 bp/yr** |
| p (D's terms of trade) | — | +1.4% | **+3.0%** |
| D taxes | — | −12.3% | **−13.1%** |
| D debt | — | −1.1% | **−1.2%** |

**Reading.**
- **The announcement works.** With nothing bought, the promise of a backstop at 200 bp lowers the sovereign
  spread by 38 bp and relaxes the D banks' constraint completely, so the lending spread falls to zero.
- **Non-divertibility moves quantities, not the spread.** It lowers the rest-point spread by only 2 bp more
  than the shipped TPI. But it doubles the output gain and multiplies the bank-capital gain by six. It comes
  with a D deposit rate 14 bp/yr lower and a terms-of-trade appreciation twice as large. Nothing is bought at
  rest, so all of this works through expectations of what the Eurosystem would do in bad states. I have not
  decomposed it further.

### 2.2 Purchases at fixed rules (quick grid)

Change in the sovereign spread, bp/yr, when the Eurosystem buys a share of the stock once, rules held fixed:

| Purchase | Shipped TPI: rest | **This run: rest** | Shipped TPI: shock | **This run: shock** |
|---:|---:|---:|---:|---:|
| 5% | −1.2 | **−2.8** | −2.4 | **−3.4** |
| 10% | −2.4 | **−5.4** | −4.6 | **−6.5** |
| 20% | −4.5 | **−10.3** | −8.4 | **−12.0** |
| 30% | −6.4 | **−14.6** | −11.5 | **−16.5** |

- **Non-divertibility raises the price impact per euro** 2.3 times at rest and 1.4 times at the shock: the
  Gertler–Karadi capacity channel.
- **It is still small.** About −0.55 bp per 1% of the stock at the shock, against −15 to −30 bp in the PEPP
  and SMP evidence. On the full grid the faithful design's slope was about twice the quick grid's, so F2 should
  be read against that.

### 2.3 The headline shock (quick grid)

The headline shock raises the quarterly default probability p^d to 2%. The D-bond price legs at the shock are
one-period premia in bp/yr, from exact reads:

| At the shock | Expected loss | Risk premium | Liquidity premium | Sovereign spread on impact | Bought |
|---|---:|---:|---:|---:|---:|
| No TPI | 428 | +119 | 172 | 414 bp | — |
| Shipped TPI | 423 | +132 | 0 | 127 bp | 0% |
| **This run** | **425** | **+92** | **0** | **132 bp** | **0%** |

- **The spread on impact falls from 414 to 132 bp with nothing bought.** The cap is never reached: the
  promise alone keeps the spread below 200 bp.
- **Non-divertibility removes 40 bp of risk premium at the shock** (+132 → +92 bp). With divertible claims the
  risk premium is actually higher than without a TPI.
- **The impulse responses on this grid are not usable.** Output falls 1.3% on impact and consumption 3.7%,
  against −0.08% without a TPI. D's terms of trade jump about 11% and D's deposit rate falls 64 bp/yr. This is the
  coarse-grid pathology the shipped TPI shows too (output −1.8% on impact). The full grid removed it for the
  shipped TPI, so the dynamics of this run wait for F2.

### 2.4 Full grid (F2)

Solved in 7.9 hours; every stage met the collocation tolerance, and the cap reached 200 bp. The purchase pair
binds at 43 of the 115 nodes, 21 of them at the corner where the Eurosystem holds the D bank's whole book.

**The announcement** (rest point, nothing bought):

| Rest point | No TPI | TPI, fitted read | TPI, exact clearing |
|---|---:|---:|---:|
| Sovereign spread | 56.5 bp | 50.4 bp | **50.7 bp** |
| Lending spread | 98 bp | 48 bp | **65 bp** |
| μ_D (the D bank's constraint multiplier) | 0.012 | — | 0.008 |
| D output, against no TPI | — | +0.32% | **+0.27%** |
| D investment | — | +0.76% | +0.59% |
| D bank net worth | — | — | +0.25% |
| p (D's terms of trade) | — | — | +1.1% |

Unlike the faithful design (F1c), the fitted and exact reads agree: the constraint stays binding at rest, away
from the kink. The quick grid overstated the announcement (there: spread −38 bp, constraint slack, output
+0.9%).

**The headline shock** (exact clearing each quarter; % or bp/yr from the unshocked path):

| | q0 | q1 | q2 | q4 | q8 | q12 | q24 | Sum q0–q24 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Sovereign spread, no TPI (level, bp) | 298 | 269 | 246 | 209 | 161 | 132 | 88 | |
| Sovereign spread, TPI (level, bp) | **200** | 193 | 185 | 172 | 148 | 129 | 90 | |
| Eurosystem book (% of the stock) | **38** | 36 | 34 | 32 | — | — | — | |
| D output, no TPI | −0.06 | +0.00 | +0.04 | +0.07 | +0.06 | +0.04 | +0.00 | +0.74 |
| D output, TPI | **−0.14** | −0.14 | −0.14 | −0.14 | −0.18 | −0.20 | −0.19 | **−4.45** |
| D consumption, no TPI | +0.08 | +0.07 | +0.06 | +0.03 | −0.02 | −0.06 | −0.09 | −0.89 |
| D consumption, TPI | **−0.59** | −0.61 | −0.62 | −0.62 | −0.64 | −0.64 | −0.54 | **−15.23** |
| D investment, no TPI | −0.49 | −0.13 | +0.09 | +0.29 | +0.31 | +0.24 | +0.09 | +4.14 |
| D investment, TPI | +0.62 | +0.57 | +0.52 | +0.40 | +0.15 | +0.04 | −0.07 | +3.59 |
| D bank net worth, no TPI | −3.81 | −1.89 | −0.68 | +0.47 | +0.85 | +0.69 | +0.34 | |
| D bank net worth, TPI | **−0.57** | −0.73 | −0.92 | −1.27 | −1.66 | −1.31 | −0.36 | |
| Eurosystem P&L (% of D's quarterly GDP) | 0 | +0.30 | +0.27 | +0.23 | +0.17 | +0.13 | +0.05 | |

Started from the exact rest point instead, the TPI path barely changes: output −0.06% on impact falling to
−0.11% by q5, consumption −0.52%, 25-quarter sums −4.0 (output) and −14.8 (consumption), 37.6% bought.

**Reading.**
- **The backstop works on prices and on banks.** It holds the spread at the 200 bp cap by buying 38% of the
  stock, removes the liquidity premium (182 → 0 bp) and most of the risk premium (+122 → +20 bp), and D bank
  net worth falls 0.6% instead of 3.8%.
- **It does not work on D's economy.** Output and consumption fall more, and for longer, than without it. This is
  robust to the starting point, unlike F1c's version of the same pattern. The model notes trace it to the
  conditioning on no default: in those quarters the bonds' excess return goes to the Eurosystem (about 0.3% of
  D's quarterly GDP a quarter) and, by the capital key, 92.9% of it on to F. I have not decomposed how much of
  the 15-point consumption loss that transfer explains.
- **Purchases still move the price little per euro.** At fixed rules 30% of the stock lowers the spread 27 bp at
  the shock (about −0.9 bp per 1%) and 10 bp at rest: still 15–30 times below the PEPP and SMP evidence.

---

## 3. Deposits cleared by heterogeneous households

### 3.1 The no-TPI union economy (full grid)

This is the no-TPI baseline every TPI result is measured against. The numbers use F1c's no-TPI rules. F2's
no-TPI economy is the same model, because the TPI switches act only when the TPI is on, and the chain re-runs it
on F2's own rules. The headline shock is differenced against the unshocked path, over a 60-quarter horizon. The
deposit-market clearing converges to 6×10⁻¹³ in three Newton steps (condition number 5).

| D consumption (% from the unshocked path) | q0 | q1 | q2 | q4 | q8 | q12 | q24 | Sum q0–q24 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Representative household (global solution) | +0.083 | +0.071 | +0.059 | +0.029 | −0.024 | −0.057 | −0.093 | −0.89 |
| Heterogeneous households, prices fixed | −0.190 | −0.164 | −0.147 | −0.130 | −0.117 | −0.109 | −0.088 | −2.87 |
| **Heterogeneous households, deposit market cleared** | **−0.146** | **−0.122** | **−0.106** | **−0.093** | **−0.092** | **−0.096** | **−0.102** | **−2.51** |

| D deposits (%) and the clearing rate | q0 | q1 | q2 | q4 | q8 | q12 | q24 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Representative household | −0.077 | −0.230 | −0.303 | −0.322 | −0.217 | −0.126 | −0.067 |
| Heterogeneous, prices fixed | −0.049 | −0.178 | −0.229 | −0.213 | −0.059 | +0.062 | +0.149 |
| Heterogeneous, deposit market cleared | −0.054 | −0.188 | −0.247 | −0.249 | −0.132 | −0.045 | −0.020 |
| Clearing shift in D's deposit rate (bp/yr) | −0.7 | −1.3 | −1.7 | −2.3 | −2.6 | −2.5 | −1.4 |

**Who cuts consumption on impact** (% change, by wealth quintile of the rest-point distribution):

| Quintile | 1 (poorest) | 2 | 3 | 4 | 5 (richest) |
|---|---:|---:|---:|---:|---:|
| Prices fixed | −0.37 | −0.19 | −0.16 | −0.15 | −0.14 |
| Deposit market cleared | −0.34 | −0.14 | −0.12 | −0.10 | −0.08 |

**Reading.**
- **The sign of the consumption response flips.** The representative household raises consumption on impact;
  the heterogeneous households cut it. Their resources fall:
  - on impact, non-labour income other than deposit interest falls by 0.77% of consumption, and labour income by
    0.09%;
  - a quarter later the deposit return falls by 1.1% of consumption, because the contracted rate fell 47 bp/yr
    on a deposit stock of about ten quarters of consumption.

  The representative household smooths through these and consumes more because the deposit rate fell. The
  heterogeneous households cut, most of all those with little to smooth with: the poorest fifth, who hold less
  than a third of average deposits, cut four times as much as the richest.
- **Clearing barely moves the deposit rate, and moves consumption only modestly.** The union market lends the
  households' extra deposits abroad at a rate at most 2.6 bp/yr lower. Even so, consumption falls 23% less on
  impact than at fixed prices, because these households are very sensitive to a persistent rate path.
- **For the TPI question:** heterogeneous households roughly triple the 25-quarter consumption cost of the
  shock relative to the representative household. Whether the TPI recovers it is the F2 run (§3.2).

### 3.2 The TPI economy (F2, full grid)

The deposit-market clearing converges to 4×10⁻¹³ in five Newton steps (condition number 5).

| D consumption (% from the unshocked path) | q0 | q1 | q2 | q4 | q8 | q12 | q24 | Sum q0–q24 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Representative household (global solution) | −0.593 | −0.611 | −0.620 | −0.624 | −0.640 | −0.641 | −0.539 | −15.23 |
| Heterogeneous households, prices fixed | −0.506 | −0.511 | −0.511 | −0.504 | −0.512 | −0.523 | −0.453 | −12.51 |
| **Heterogeneous households, deposit market cleared** | **−0.568** | **−0.574** | **−0.574** | **−0.567** | **−0.572** | **−0.578** | **−0.486** | **−13.82** |
| Clearing shift in D's deposit rate (bp/yr) | +0.1 | +0.2 | +0.3 | +0.6 | +1.2 | +1.6 | +2.0 | |

**Who cuts consumption on impact** (cleared): poorest fifth −0.70%, then −0.63, −0.56, −0.52, richest −0.48%.

**Reading.**
- **With the TPI the households' loss is broad and persistent.** Every quintile cuts about half a percent, and
  keeps cutting: a lasting income loss, which the representative household and the heterogeneous households
  both price in. The heterogeneous households cut slightly less (−13.8 against −15.2 over 25 quarters),
  because part of the loss lands on deposit income, which the poorest do not have.
- **Against no TPI, heterogeneous households make the TPI look somewhat less costly:** their 25-quarter
  consumption loss rises from −2.5 to −13.8, against −0.9 to −15.2 for the representative household. The
  sign of the TPI's consumption effect is the same in both.
- **Clearing again moves D's deposit rate by at most 2 bp/yr.**

---

## 4. Accuracy

**The global solution.**
- Q8 (quick grid): the fitted and exact-clearing rest points agree to 2 bp in the spread, so the announcement
  result is robust. The TPI economy's impulse responses are not: output −1.3% on impact against −0.08%
  without a TPI. That is the known coarse-grid pathology, and the full grid removed it for the shipped TPI.
- The faithful design on the full grid (F1c) showed that a TPI economy can rest next to the banks'
  constraint kink, where the fitted and exact reads disagree (`tpi_fix_explained.md` §3.4). F2 gets the same
  check (exact rest point and an impulse response from it) in the chain.

**Heterogeneous households.**
- Deposit-market clearing: Newton to 6×10⁻¹³ on both paths. The boundary conditions are the rest point's own
  distribution, policy and discount factor. Horizon 60 quarters; the deposit contracted in the last quarter
  pays its own rate.
- Every price cleared along the path: not determined. After 120 Levenberg–Marquardt steps the best iterate
  still leaves a deposit gap of 3×10⁻⁴ of deposits in the first decade, about as large as the gap at fixed
  prices. A hundredfold larger finite-difference step barely helps (4.5 → 3.3 × 10⁻⁴ over its 60 steps), so it
  is not numerical noise. The Jacobian shows why:
  - Inside the period map, a change in household deposits is a change in D's consumption flow (deposits are
    7.7 quarters of GDP). It moves the terms of trade so much that the deposit-UIP residual responds 100 times
    more than through the κ_nfa premium, and with the opposite sign.
  - The households' own deposits then respond 113% to the deposit path the map uses, so the clearing condition
    has a slope of 0.13 instead of 1. Over 60 quarters the feedback compounds into an unstable system.
  - The cause is the hybrid: the households see the path, but every other agent forms expectations with the
    representative household's global rules, which do not react to the households' saving.
- **The fix inside the global framework is Krusell–Smith approximate aggregation.** The steps:
  1. replace D's representative-household Euler in the period map with an aggregate consumption rule over the
     model's existing states;
  2. solve the global model with it;
  3. simulate the heterogeneous households on the solved rules;
  4. refit the rule, and repeat until it is consistent.

  Every agent's expectations then reflect the heterogeneous households. It is a multi-day project: an outer
  loop around a global solve that takes 2–11 hours per pass.

---

## 5. Files

| File | Contents |
|---|---|
| `Claude files/docs/tpi_funding_fix_ha.md` | this note |
| `Claude files/docs/tpi_fix_explained.md` (and `.pdf`) | the faithful-TPI note; §3.4 now has the full-grid run F1c |
| `Claude files/diagnostics/tpi_national/ha_dmkt.py` | D's deposit market cleared through the heterogeneous households (§3) |
| `Claude files/diagnostics/tpi_national/ha_pe.py`, `ha_dist.py` | the households at fixed prices; consumption by wealth quintile |
| `Claude files/diagnostics/tpi_national/ha_path.py`, `ha_homotopy.py`, `ha_ge2.py` | the every-price clearing attempt and its diagnostics (§4) |
| `Claude files/diagnostics/tpi_national/exact_rest.py`, `exact_irf.py` | rest points under exact clearing; impulse responses from the exact rest point |
| `results/TPI_NATIONAL/` | runs Q8, F2 (when done), F1c, `uni`: JSON, logs and HA outputs (git-ignored) |
