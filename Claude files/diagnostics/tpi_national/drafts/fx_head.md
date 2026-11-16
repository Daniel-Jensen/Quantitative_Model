# The TPI fix, explained and made faithful to the literature and the euro-area record

*Working note, 14 October 2026. Global-projection model, code of branch `domestic-funding-fix` (`0dd0732`) plus
the experimental switches listed in §6, all off by default and all nesting the shipped model exactly.
Simulations are on the 25-point quick grid unless marked "full grid". The full-grid run of the faithful design
(F1c, 115 points) is in §3.4, and it revises several quick-grid magnitudes.*

---

## 0. In one page

**What the first fix was.** To make the TPI work under national deposit funding, the first fix:
- made the Eurosystem's claim non-divertible;
- let the default probability fall with the private float, via a logit with elasticity 1;
- returned all of the Eurosystem's profit and loss to D.

It worked mechanically, holding a 200 bp cap and halving the spread at the headline shock. Measured against
the literature and the record, though, three of its four pieces are wrong.

| Piece of the first fix | Verdict | Faithful replacement |
|------------|--------------------------|-----------------|
| national deposits + profits returned to D | **wrong root fix.** D's deposit rate collapsing to −710 bp/yr in default has no counterpart in Greece, where deposit rates were about 3.9% in 2011 and the Eurosystem was the marginal funder (ELA about €107–124 bn, TARGET2 about −€97 bn) | **a TARGET2 funding margin.** A Eurosystem facility at the parity rate *is* the union deposit market (§2.1) |
| default response, logit φ = 1 | **about 4.5× too strong** against the only estimate | **Costain–Nuño–Thomas linear rollover channel**, θ_q = 0.0046 from θ = 2.02×10⁻⁵ per € bn (§2.2) |
| κ_D = 1 (all P&L to D) | not the TPI's documented rule; PSPP is 80/20; SMP income shared, then returned to Greece | **capital key** baseline; PSPP retention as a sensitivity, numerically fragile here (§2.3) |
| non-divertible reserves | **against the evidence**: reserves use balance-sheet space (Diamond–Jiang–Ma) | **divertible**, as shipped; non-divertible as the Gertler–Karadi bound (§2.6) |
| fixed 200 bp cap | not the TPI's rule ("unwarranted, disorderly") | **cap = fundamental spread + 50 bp** (§2.5) |
| pari passu | ✓ matches the TPI | unchanged (§2.4) |

**What the faithful TPI does** (quick grid, §3; full-grid revisions in brackets, §3.4):
- **It holds the spread at the fundamental level + 50 bp**: 414 → 293 bp on impact of the headline shock.
  [Full grid: the shock lifts the no-TPI spread only to 298 bp, against a 293 bp cap, so the TPI buys 0–6% of
  the stock.]
- **It extracts the whole risk premium**: from +119 bp to −15 bp. [Full grid: from +122 to +91 bp; the
  liquidity premium falls from 182 to 10 bp.]
- **It is an insurance contract.** In no-default quarters the Eurosystem earns about 1.5% of D's quarterly GDP,
  93% of it going to F. In a realized default it loses about 31% of D's quarterly GDP, 93% borne by F; D's bank
  net worth then rises 2% instead of falling 32%.
- **It works only with the TARGET2 margin.** Under national deposits it never binds, and purchases raise
  spreads.

**What it cannot do faithfully yet.** Its price impact per euro is **12–25 times smaller** than the evidence
on the full grid: about −1.2 bp per 1% of the stock at the shock (−0.6 bp on the quick grid), against −15 to −30
bp (PEPP, SMP). So on the quick grid holding the spread near fundamentals takes about 80% of the stock. The
literature's slope comes from risk-averse holders alongside inelastic ones (Costain–Nuño–Thomas), and the same
ingredient would lift the model's stress spreads toward the 1,000+ bp seen in Greece. That is the next model
change (§4).

**What the full grid could not settle (§3.4).**
- The announcement effect of the faithful design: the fitted rules say the spread falls 80 bp, while exact
  clearing says it rises 16 bp. The rest point lies next to the banks' constraint kink.
- Its impulse responses show weaker output and consumption with the TPI than without. The same pattern was ten
  times larger on the quick grid, so most of it is approximation error.

Both need Bocola's 9-node resolution in s.

---
