## 4. The gap that remains: private demand for the D bond is far too elastic

Even the faithful design leaves one large mismatch with the evidence: **the price impact of a purchase**.

| Measure | Price impact per 1% of the debt stock |
|---|---:|
| Eser & Schwaab (2016, *JFE*), SMP, five-year yields, stress period | about **−30 bp** (−3 bp per 0.1%) |
| Costain–Nuño–Thomas (2025), PEPP announcement, Italy, five-year (≈90 bp for ≈6% of the debt) | about **−15 bp** |
| Altavilla–Giannone–Lenza (2016), OMT announcement, Italy/Spain two-year | −200 bp in total, with no purchases at all |
| **Model, faithful baseline Q6, fixed rules, at the headline shock** | **−0.6 bp** |
| Model, faithful baseline Q6, equilibrium (414 → 293 bp for an 83% book) | −1.5 bp |

The model is **10–50 times too small**, and the reason is structural. The D bank values the bond almost
linearly: its first-order condition (P) has no quantity in it, so in equilibrium it absorbs whatever the
Eurosystem does not hold, at almost the same price. The only curvature is the response of its kernel Ω to
net worth (formula (I) of the companion note). In the literature the slope comes from elsewhere:
- **Preferred-habitat and limited-risk-bearing investors** (Vayanos & Vila 2021; Gabaix & Maggiori 2015).
  Demand is downward-sloping, and a purchase removes risk the marginal holder is unwilling to bear.
- **Costain, Nuño & Thomas specifically** combine inelastic preferred-habitat investors, holding 44% of the debt
  in their calibration, with risk-averse arbitrageurs holding 56%. The arbitrageurs' risk aversion γ is
  identified from the German term premium. The default sensitivity θ is identified from the PEPP
  announcement.

**What a faithful version would add.** A mean-variance term in the D bank's (and the F bank's) bond condition:

$$
q_t\big(\mathbb{E}_t[\Omega]R_t+\lambda\mu_t\big)\big(1+\gamma\,\mathcal V_t\,x_t-\gamma\,\mathcal V_{ss}x_{ss}\big)=\mathbb{E}_t[\Omega\,\Xi],
\qquad x_t=\frac{q_tb^D_{D,t+1}}{n_{D,t}},\quad \mathcal V_t=\mathrm{Var}_t\Big(\frac{\Xi_{t+1}}{q_t}\Big).
$$

Here x is the bank's D-bond exposure relative to its net worth and 𝒱 the one-quarter variance of the bond's
return, computed on the existing quadrature. The term is anchored at the steady state, so the steady state does
not change. A purchase lowers x, and with it the premium, which is the credit-risk extraction of
Costain–Nuño–Thomas.

**Calibration arithmetic.** At the headline shock 𝒱 ≈ π(1−π)·0.55² ≈ 0.006 per quarter. Buying 1% of the stock
lowers x by about 0.004. To reach −15 to −30 bp per 1%, γ must be about 15–30.

At that γ the same term adds several hundred basis points to the stress spread, because 𝒱 is about 20 times
its rest value. In the shipped model that would be a problem. In the data it is not: Greek spreads were well
above 1,000 bp through 2011–12, against the model's 300–400 bp at the headline shock. **So the missing slope and
the missing level of stress spreads are the same missing ingredient.** The two-type structure of
Costain–Nuño–Thomas reconciles a large slope with small normal-time premia, and it is the natural next
extension. It is not implemented here: it changes the no-TPI model too, and needs its own calibration (§5).

**Why this matters for the TPI's design.** With elastic private demand, holding the spread at
"fundamental + 50 bp" in stress requires buying most of the stock (83% on the quick grid). With CNT-sized
slopes the same compression would take a few percent of the stock, the scale of actual SMP purchases (Greek
holdings of €40–50 bn, about 12–14% of Greek debt). The model's TPI book at the headline shock, about 20% of
annual GDP, is close to the Greek SMP's €40–50 bn (about 21% of Greek 2011 GDP) only because the model's debt
stock is small: 25% of annual GDP, calibrated to bank exposure, against Greece's 172%.
