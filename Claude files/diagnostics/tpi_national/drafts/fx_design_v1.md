## 3. Simulations

### 3.1 Design

Every run solves the economy without and with the TPI (Bocola's ladder, then the cap homotopy) on the
25-point quick grid. It then compares the two:
- rest points (the announcement);
- the headline shock (the quarterly default probability p^d rises to 2%), cleared exactly each quarter;
- the price impact of forced purchases at fixed rules;
- the D-bond price legs;
- a realized default at quarter 1 (§3.3).

The faithful baseline (`Q6`) combines every element of §2; each other run changes one element:

| Run | Funding | Default response (θ_q) | P&L: D's share κ_D | Activation | Eurosystem claim | What it isolates |
|---|---|---:|---:|---|---|---|
| **`Q6` faithful baseline** | union (TARGET2) | CNT, 0.0046 | 0.071 (capital key) | fundamental + 50 bp | divertible (λ_Z = λ) | — |
| `Q7c` | union | **0** | 0.071 | fundamental + 50 bp | divertible | the CNT rollover channel |
| `Q4c` | union | 0.0046 | 0.071 | fundamental + 50 bp | **non-divertible** (λ_Z = 0) | Gertler–Karadi capacity |
| `Q5c` | **national** | 0.0046 | 0.071 | fundamental + 50 bp | divertible | the TARGET2 margin |
| `Q3`–`Q5`, `Q3b`–`Q5b` | union / national | 0.0046 | **0.814** on income (losses 0.814 or 0.071) | fundamental + 50 bp | either | PSPP-style retention (no-TPI baseline unsolvable, §2.3) |
| `F1c` (full grid) | union | 0.0046 | 0.071 | fundamental + 50 bp | divertible | the baseline at 115 points (§3.4) |
| `uni` (previous note) | union | 0 | 0.071 | **fixed 200 bp** | divertible | the shipped TPI |
| `nat_phi1_lz0_key1` (previous note) | national | logit φ = 1 | 1 | fixed 200 bp | non-divertible | the first, mechanical fix |

**Two cautions about the quick grid.**
- **Union impulse responses with the TPI are not usable there.** They show a known coarse-grid pathology
  (`tpi_backstop_plan.md` §9.2: impact output −1.8% and no purchase on impact). The union rest points,
  price impacts and price legs are comparable across runs; the union impulse responses need the full grid.
  The national impulse responses are clean.
- **The full-grid run of the faithful baseline (F1c) is the reference for dynamics** (§3.4). It needs 4–5 hours
  of awake machine time.
