### 3.2 Results on the quick grid

**Can the backstop be run, and what does it buy?**

| Run | Funding | CNT θ_q | κ_D | Claim | Lowest tolerance held (target +50 bp) | Bought on impact at the shock | Sovereign spread on impact: no TPI → TPI |
|---|---|---:|---:|---|---:|---:|---|
| **Q6 faithful baseline** | union | 0.0046 | 0.071 | divertible | **+50 bp** ✓ | **83%** of the stock | 414 → **293 bp** |
| Q7c | union | **0** | 0.071 | divertible | +68.75 bp | 81% | 414 → 312 bp |
| Q4c | union | 0.0046 | 0.071 | **non-div.** | +68.75 bp | 82% | 414 → 312 bp |
| **Q5c** | **national** | 0.0046 | 0.071 | divertible | **+344 bp**, never binds | 0 | 339.5 → 339.5 bp |
| Q3–Q5, Q3b–Q5b | union / national | 0.0046 | **0.814** (income) | either | no-TPI baseline unsolvable (§2.3) | — | — |
| `uni` (shipped TPI) | union | 0 | 0.071 | divertible | fixed 200 bp ✓ | coarse-grid pathology | 414 → 127 bp |

- **Funding is the precondition.** With the TARGET2 margin the warranted cap holds; under national deposits it
  never binds, and a purchase *raises* the spread (+8.6 bp per 30% of the stock at the shock).
- **The CNT channel is what lets the cap reach the 50 bp target.** Without it the walk stalls at +68.75 bp.
- **Every union variant has to buy about 80% of the stock**, which is the edge of the state box, to bring the
  spread down to the fundamental level plus the tolerance. That is the elasticity gap of §4.

**Purchases at fixed rules** (change in the sovereign spread, bp/yr):

| Purchase | Q6 rest | Q6 shock | Q7c shock (no CNT) | Q4c shock (non-div.) | Q5c shock (national) | `uni` shock (shipped) |
|---:|---:|---:|---:|---:|---:|---:|
| 5% | −1.4 | −4.2 | −2.4 | −5.4 | +1.7 | −2.4 |
| 10% | −2.9 | −7.3 | −4.6 | −10.7 | +3.4 | −4.6 |
| 20% | −6.3 | −13.1 | −8.6 | −19.7 | +6.4 | −8.4 |
| 30% | −9.6 | **−19.3** | −12.0 | **−24.9** | **+8.6** | −11.5 |

Channel attribution at the shock, for 30% of the stock (Q6 = −19.3 bp):
- **The CNT rollover channel** adds −7.3 bp: Q6 against Q7c, about 38% of Q6.
- **The Gertler–Karadi capacity channel**, as an upper bound, would add a further −5.6 bp: Q4c against Q6,
  about 29% on top.
- **The risk-transfer remainder** is −12.0 bp (Q7c).

For comparison, Krishnamurthy, Nagel & Vissing-Jorgensen attribute 37% of SMP/OMT to default risk and 50% to
segmentation.

**The D-bond price legs** (one-period premia, bp/yr; union rules):

| | expected loss | risk premium | liquidity premium | Ω^d/Ω^nd |
|---|---:|---:|---:|---:|
| no TPI, rest | 21.2 | +8.4 | 93.5 | 1.29 |
| no TPI, shock | 428.1 | **+118.9** | 171.5 | 1.25 |
| faithful TPI (Q6), rest | 21.5 | +8.6 | 0.6 | 1.25 |
| faithful TPI (Q6), shock | 415.7 | **−14.7** | 0.0 | 0.92 |

At the shock the faithful TPI **extracts the whole risk premium** (+119 → −15 bp), and the constraint stops
binding (liquidity premium 0). The expected loss falls slightly through the CNT channel (428 → 416 bp).

**The announcement** (rest point, nothing bought): in Q6 the sovereign spread moves from 62.5 to 66.1 bp, and
the lending spread from 94.8 to 5.9 bp. Output at rest is 0.34% *lower*. The union rest points with the TPI are
the least reliable objects on this grid: Q7c even shows a negative spread. The full-grid figures are in §3.4.
