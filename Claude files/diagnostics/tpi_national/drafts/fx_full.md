### 3.4 Full grid: the faithful baseline (F1c)

The 115-point run (union / TARGET2 funding, CNT θ_q = 0.0046, capital key, cap = fundamental + 50 bp,
divertible reserves) solved in 10.8 hours. Every stage met the collocation tolerance (max|F| ≤ 10⁻⁹), and the
cap reached its +50 bp target. The purchase pair binds at 39 of the 115 nodes, 3 of them at the corner where the
Eurosystem holds the D bank's whole book.

The full grid confirms one result of the quick grid and overturns two others.

**What holds: the price impact per euro is far below the evidence.** Purchases at fixed rules (change in the
sovereign spread, bp/yr):

| Purchase (% of the stock) | Rest point | Headline shock |
|---:|---:|---:|
| 5% | −2.0 | −5.5 |
| 10% | −3.9 | −11.1 |
| 20% | −7.8 | −22.9 |
| 30% | −11.2 | **−36.0** |

At the shock that is about −1.2 bp per 1% of the stock, twice the quick grid's −0.6 bp but still **12–25 times
below** the −15 to −30 bp of the PEPP and SMP evidence. The elasticity gap of §4 is real, not a grid artefact.

**What does not hold (1): the announcement effect.** On this grid the TPI economy's rules are not accurate
near its own rest point:

| Rest point | No TPI | Faithful TPI, fitted read | Faithful TPI, exact clearing |
|---|---:|---:|---:|
| Sovereign spread | 56.5 bp | −23.8 bp | **72.5 bp** |
| μ_D (the D bank's constraint multiplier) | 0.012 | 0 (slack) | 0.010 (binds) |
| Lending spread | 98 bp | 0 bp | 85 bp |
| D output, against no TPI | — | −0.60% | −0.24% |

The "fitted read" evaluates the solved rules at the rest point. "Exact clearing" solves the period map exactly
each quarter, given the rules for the continuation, and runs the economy without shocks until it settles (400
quarters). Without the TPI the two agree to 0.1 bp. With the TPI they disagree in sign:
- the fitted rules put the rest point on the constraint's kink (μ_D = 0), where the Chebyshev solution is least
  accurate (the Gibbs problem in the model notes);
- cleared exactly, the economy settles on the other side of the kink: the constraint binds (μ_D = 0.010), the
  market spread is 72.5 bp, and purchases are negligible. The warranted cap there is 91 bp (fundamental 41 bp +
  50 bp), so the floor does not bind at rest.

**The likely cause is how little the warranted cap promises in bad states.** At the headline shock it is 293 bp
(fundamental 243 bp + 50 bp), against 200 bp for the fixed cap. So the announcement only partly relaxes the D
banks' constraint and leaves the rest point next to the kink. With the fixed 200 bp cap (run Q8, reported
separately) the constraint goes slack at rest in both reads, and they agree to 2 bp. **The announcement effect
of the faithful design is therefore not identified at 115 points, not even in sign.**

**What does not hold (2): the impact of the headline shock.** Each economy is shocked from its rest point and
cleared exactly each quarter. The TPI economy is run from both rest points:

| Headline shock | No TPI | Faithful TPI, from the fitted rest | Faithful TPI, from the exact rest |
|---|---:|---:|---:|
| Sovereign spread on impact | 297.9 bp | 289.5 bp | 293.1 bp |
| Bought on impact | — | 0% | 6.3% of the stock |
| Lending spread on impact (change) | +79 bp | +10 bp | +55 bp |
| D output on impact | −0.059% | −0.138% | −0.151% |
| D output, 25-quarter sum (% × quarters) | +0.74 | −1.90 | −1.97 |
| D consumption on impact | +0.08% | −0.34% | −0.21% |
| D consumption, 25-quarter sum | −0.89 | −7.11 | −6.41 |
| D bank net worth on impact | −3.81% | −2.70% | −3.52% |

- **The shock barely reaches the cap.** On the full grid the no-TPI impact spread is 298 bp, not the quick
  grid's 414 bp. The warranted cap at the shock is 293 bp, so the TPI buys between 0% and 6% of the stock. From
  the exact rest point it holds the spread exactly at the cap (293.1 bp). The quick grid's "414 → 293 bp, 83%
  bought" (§3.2) was a coarse-grid result.
- **The bond-price legs still move a lot.** At the shock the TPI lowers the risk premium from +122 to +91 bp and
  the liquidity premium from 182 to 10 bp. The spread falls only 4–8 bp because the deposit-rate and continuation
  legs offset most of it.
- **Output and consumption do worse with the TPI, from both starting points.** On the 25-point grid the same
  pattern was about ten times larger (Q6: output −1.5%, consumption −3.2% at q1). It shrinks with resolution, so
  most of it is probably approximation error. I do not report it as a result. Settling it needs Bocola's 9-node
  resolution in s, which is about 4 hours of no-TPI solve plus the TPI walk.

**What this means for the note's conclusions.**
- The elasticity gap (§0, §4) stands at the full grid, narrowed from 10–50× to 12–25×.
- The quick-grid magnitudes for the faithful design are indicative only: the 414 → 293 bp compression, the 80%
  book and the §3.3 insurance numbers.
- The warranted cap is harder to solve accurately than a fixed cap. Its weaker promise in bad states leaves the
  rest point next to the banks' constraint kink. A paper version should either run it at 9 nodes in s or report
  the fitted/exact pair as the honest object, as the model notes already do for μ.
