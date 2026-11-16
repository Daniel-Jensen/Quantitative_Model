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
## 1. Starting point: why the shipped TPI fails under national deposits

The companion note (`tpi_national_analysis.md`) established five facts, all on the quick grid:

- **The shipped TPI does nothing under national deposits.** The lowest cap the solver can hold is 494 bp;
  the announcement moves the sovereign spread by +0.1 bp; a forced purchase *raises* the spread.
- **A purchase is a swap, not a balance-sheet relief.** The Eurosystem pays with a claim Z of equal value
  and the same divertibility, so nothing moves today. The price can move only through tomorrow's banker
  discount factor Ω.
- **Under national deposits Ω falls in both regimes.** D's deposit market must clear at home, so a
  purchase reaches D households' marginal utility through D's own deposit rate.
- **There is no risk premium to remove.** The D bond is a hedge for the D banker: its risk premium is
  −16 bp/yr at the shock, because captive depositors accept a −710 bp/yr real deposit rate in default.
- **Most of the spread is untouchable.** It is expected default loss, and in the shipped model that is
  exogenous.

The first, mechanical fix combined three changes:
- non-divertible Eurosystem claims;
- a default probability that falls with the private float, as logit π = s + φ·ln((B−M)/B) with φ = 1;
- the Eurosystem's profit and loss returned to D.

It made the cap feasible at 200 bp and halved the spread at the headline shock. This note asks whether each
piece of that fix survives contact with the literature and the euro-area record, replaces the pieces that do
not, and re-runs the model with the faithful version.

---

## 2. The fix, element by element

Every element acts on one term of the D bank's pricing condition (`tpi_national_analysis.md` §3):

$$
q_t=\frac{\mathbb{E}_t[\Omega_{t+1}\Xi_{t+1}]}{\mathbb{E}_t[\Omega_{t+1}]R_t+\lambda\mu_t}
\quad\Longleftrightarrow\quad
\underbrace{\tfrac{\mathbb{E}_t[\Xi_{t+1}]}{q_t}-R_t}_{\text{expected excess return}}
=\underbrace{\tfrac{\lambda\mu_t}{\mathbb{E}_t[\Omega_{t+1}]}}_{\text{liquidity premium}}
\;\underbrace{-\tfrac{\mathrm{Cov}_t(\Omega_{t+1},\Xi_{t+1})}{q_t\,\mathbb{E}_t[\Omega_{t+1}]}}_{\text{risk premium}} .
$$

Here q is the D bond's price, Ξ its payoff (coupon plus continuation value, times the 45% recovery in
default), Ω the banker's stochastic discount factor, R the gross D deposit rate, λ the divertible fraction of
assets and μ the multiplier on the bank's incentive constraint. The expected default loss sits inside E[Ξ].

### 2.1 Funding: give D a TARGET2 margin (replaces "remit the P&L to D" as the root fix)

**What goes wrong without it.** Under national deposits, each bank is funded only by its own households.
Every change in the D bank's funding needs, including the extra saving the Eurosystem's profits drain from D,
must clear in D's own deposit market at D's own real rate. In the default branch that rate falls to about
−710 bp/yr. That is what turns the D bond into a hedge (Ω^d/Ω^nd = 0.95) and makes purchases tighten the
constraint.

**What the euro area actually did.** In 2010–12 Greek banks lost about a third of their deposits, and the
Eurosystem became their marginal funder:
- Emergency Liquidity Assistance peaked at roughly €107–124 bn in 2012, and about a third of Greek bank
  balance sheets were Eurosystem-financed;
- the Bank of Greece's TARGET2 liability was about −€97 bn by mid-2011.

The Eurosystem intermediated the cross-border flow: deposits that left Greece landed in core banks as reserves,
and Greek banks borrowed the same amount from their national central bank. And the Greek deposit rate did not
collapse. Household time-deposit rates (up to one year) averaged 3.88% in April 2011 and were rising (Bank of
Greece data), and banks were paying about 2 points above normal to keep deposits. With inflation of about 3%,
that is a real rate near zero, not −7%.

**What it is in the model.** A Eurosystem facility ℓ that lends to the D bank at the rate at which the
Eurosystem itself borrows (from F banks' reserves, at F's rate converted at the expected terms of trade p) is
*equivalent to the union deposit market*. The D bank arbitrages deposits against the facility, so

$$
R_{D,t}=R^\ell_t=R_{F,t}\frac{\mathbb{E}_t[p_{t+1}]}{p_t}-\kappa\,\frac{\ell_t-\bar\ell}{\bar Y_D}.
$$

- The facility's profit and loss per unit is R_F·(E[p′] − p′)/p: zero in expectation, terms-of-trade noise only.
- The two national clearings plus the facility add up to the union clearing, with V = −ℓ.
- Real-rate parity is the code's UIP condition, with the penalty slope κ in the role of κ_nfa.

So the union deposit market in `point_map.py` is the reduced form of TARGET2 intermediation, and running the
TPI on it, as the shipped code does, is the faithful choice. What is *not* faithful is the national market
without a lender of last resort. That describes a currency board or a system under capital controls, not the
euro area. (A fully explicit facility with a collateral constraint is option D3 of `deposit_portfolio_review.md`;
in the slack-collateral case it collapses to the equation above.)

**Calibration.** Union market as shipped (κ_nfa = 0.01). D3 would add the collateral constraint and the
ELA–MRO penalty spread.

### 2.2 Default risk that responds to the backstop: Costain–Nuño–Thomas, not φ = 1

**The literature's version.** Costain, Nuño & Thomas (2025, *JF*; working paper November 2021) derive the
default probability from rollover crises à la Calvo (1988) and Cole & Kehoe (2000):

- investors sometimes coordinate on not rolling over, with Poisson rate η;
- the government then defaults if the stochastic cost of default χ (c.d.f. Φ) is below the emergency taxes it
  would need to raise.

So

$$
\psi_t=\eta\,\Phi\!\Big(\int_t^\infty e^{-(\hat r+\varphi)(s-t)}\big(d_s+S_s(0)+\bar\Gamma\big)\,ds\Big)
\;\;\overset{\Phi\ \text{uniform}}{\Longrightarrow}\;\;
\psi_t=\psi_{ss}+\theta\,(F_t-F_{ss}),\qquad \theta\equiv\eta\Phi' .
$$

- d_s are future deficits, and S_s(0) **the redemptions of bonds in private hands**, net of central-bank
  holdings; F is their discounted sum ("fiscal pressure", at r̂ + φ = 0.5 a year).
- By buying, the Eurosystem lowers future private redemptions, hence F, hence the default probability.
- They estimate **θ = 2.02×10⁻⁵ per € billion** from the PEPP announcement, which lowered Italian 5-year yields
  by almost 90 bp, with ψ_ss = 14 bp a year and a 25% haircut.

**Their key finding.** In their calibration the expected default loss falls by only about **1 bp a year**.
Most of the PEPP effect is the **credit risk premium**: risk-averse arbitrageurs demand less compensation when
less default risk is outstanding ("default risk extraction"). The default-probability channel is real but
small; risk extraction does the heavy lifting.

**Mapping into this model.** Convert θ to GDP units with Italian 2019 GDP (€1.79 tn): θ_GDP = 0.036 per year
per unit of F/GDP. One unit of debt held by the Eurosystem lowers the present value of private redemptions by
c_PV = δe^{−0.125}/(1−(1−δ)e^{−0.125}) = 0.296, given the perpetuity's 5.6% quarterly redemption rate and their
discounting. Scaling to Greek debt (about 172% of GDP) and to quarterly frequency gives

$$
\pi_t=\max\Big(p^d(s_t)-\theta_q\,\frac{M_{t+1}}{B_{t+1}},\;10^{-5}\Big),\qquad
\theta_q=\tfrac14\,\theta_{GDP}\,c_{PV}\,(B/Y)_{GR}=0.0046 .
$$

- **Size:** holding 43% of the stock (the full-grid union purchase at the headline shock) lowers the quarterly
  default probability by **0.20 pp** (from 2.0% to 1.8%).
- **Against the first fix:** φ = 1 gave 0.86 pp. **The first fix overstated the default channel about 4.5-fold
  relative to the only quantitative estimate in the literature.** The linear form is also theirs; the logit made
  the channel proportional to π.
- **Nesting:** θ_q = 0 is the shipped model, and the rest point is untouched (M = 0 there).

**Cross-checks.** Krishnamurthy, Nagel & Vissing-Jorgensen (2018) attribute 37% of the SMP/OMT yield decline
to *default risk*, measured from CDS. That includes the default risk *premium*, so it is consistent with a
small direct probability effect plus a large premium effect. Bocola & Dovis (2019, *AER*) find that
non-fundamental (rollover) risk played a limited role in Italy in 2008–12.

### 2.3 Who keeps the Eurosystem's profit and loss: capital key as the baseline, PSPP as the sensitivity

**Real-world arrangements.**

| Programme | Loss and income sharing | D's share κ_D in the model |
|---|---|---:|
| SMP (2010–12) | Income accrued to all euro-area NCBs (shared), but the Eurogroup (27 Nov 2012) decided to pass Greece an amount equal to the SMP income accruing to their NCBs from budget year 2013 (€7.8 bn of net interest income in 2012–16, ECB); transfers suspended 2015, restored 2018 | ≈ 0.071 by design; ≈ 1 for income in practice |
| OMT (2012, never used) | Pari passu; sharing not specified in the sources found | — |
| PSPP / PEPP | **80% of government-bond purchases not loss-shared:** NCBs buy their own sovereign at their own risk, and the income on those holdings is computed at the MRO rate (the deposit facility rate since Nov 2024), so the excess yield stays national; 20% shared by capital key | **0.8 + 0.2 × 0.071 = 0.814** |
| TPI (2022, never used) | Pari passu; sharing not specified in the sources found | open |

**Why it matters in the model.** In every no-default quarter the Eurosystem's book earns the bond's excess
return. Under the capital key 92.9% of it goes to F. Under national funding that drain caused the deposit-rate
explosion (1,158–1,215 bp/yr in quarter 2). Costain–Nuño–Thomas *assume* national retention for their own
mechanism to work: during a crisis "the central bank transfers income ... back to the national government that
issued those bonds", which they justify by current Eurosystem practice.

**Baseline used: κ_D = 0.071 (capital key), with PSPP-style national retention as the sensitivity.**
- **Why the capital key.** The TPI's sharing rule is not public, the OMT's is not either in the sources found,
  and the SMP's income accrued to all NCBs.
- **Why not PSPP retention.** National retention of the Eurosystem's *income* (κ_D = 0.814 on profits),
  combined with the CNT default channel, makes the no-TPI baseline unsolvable on the current state box. That
  holds whether losses are kept nationally or shared. The trouble sits at box corners where the obligation O
  is mispriced relative to the book M and the P&L reaches ±20% of quarterly GDP. Retaining 81% of that
  nationally moves D's issuance far outside the ±12% band for bank bond holdings. These are states the economy
  never visits, so this is a numerical limit, not an economic result. Narrowing the obligation band (`o_band`,
  0.25 now) is the likely cure (§5).

With a TARGET2 margin (§2.1), the drain of profits to F matters much less for feasibility. The choice is about
distribution: who bears losses in default (§3.3).

### 2.4 Creditor treatment: pari passu

The TPI press release states that the Eurosystem accepts the same (pari passu) treatment as private creditors,
and the model's Eurosystem book takes the same 55% haircut. The Greek precedent was different: the Eurosystem's
SMP holdings (about €40–50 bn) were swapped out of the 2012 PSI, so the Eurosystem was de facto senior. Fitch
noted that this meant private bondholders "must take on more losses" for a given debt reduction, and the OMT,
then the TPI, committed to pari passu precisely to avoid that subordination. The model follows the TPI. A senior
Eurosystem would *raise* private recovery risk and is not simulated.

### 2.5 Activation: only "unwarranted" spreads (replaces the fixed 200 bp cap)

**Real-world rule.** The TPI counters "unwarranted, disorderly market dynamics". Its other features:
- eligibility rests on fiscal and macro criteria;
- purchases are in the secondary market, at 1–10 year maturities, with no ex-ante limit;
- purchases end "upon a durable improvement in transmission, or ... if persistent tensions are due to country
  fundamentals".

A fixed 200 bp cap does not do this. At the headline shock it binds even though most of that spread is the
fundamental expected loss.

**The model's version.** Compute the *fundamental spread* as the yield gap of a risk-neutral perpetuity that
faces only the exogenous default process p^d(s), at the SS deposit rate. Then let the Eurosystem cap the spread
at the fundamental spread plus a tolerance c:

$$
y^D_t-y^F_t\;\le\;\underbrace{\Delta y^{\rm fund}(s_t)}_{\text{warranted}}+c,\qquad c=50\ \text{bp} .
$$

| p^d per quarter | 0.10% (rest) | 0.50% | 1.00% | 1.58% | 2.00% (headline) | 4.80% (box edge) |
|---|---:|---:|---:|---:|---:|---:|
| fundamental spread, bp/yr | 41 | 98 | 152 | 206 | **243** | 468 |
| cap = fundamental + 50 bp | 91 | 148 | 202 | 256 | **293** | 518 |

The market spread at the headline shock is 302 bp on the full union grid, so only about 60 bp of it is
"unwarranted" in this sense. The 50 bp tolerance keeps the TPI inactive at rest: the full-grid rest spread is
56 bp against a cap of 91 bp. This is the operational version of the draft's "disorderly market" design, which
targets the default-warranted price.

### 2.6 Reserves and bank capacity: keep λ_Z = λ (Gertler–Karadi λ_Z = 0 as the bound)

**Theory.** In Gertler–Karadi the central bank's liabilities do not use private-bank capacity, so purchases
relax the incentive constraint. In the model, λ_Z = 0 does that.

**Evidence.** QE-created reserves *do* use bank balance-sheet capacity. Diamond, Jiang & Ma (2024, *JFE*,
"The reserve supply channel of unconventional monetary policy") find that each dollar of reserves crowds out
about 7.7 cents of bank lending and raises loan rates by 7.4 bp, through balance-sheet costs such as the
leverage ratio. In the euro area, reserves and euro-area sovereign bonds both carry a zero risk weight and both
count in the leverage-ratio exposure (apart from the temporary exclusion of central-bank reserves during the
pandemic). A swap of one for the other relieves no regulatory constraint.

**Faithful baseline:** λ_Z = λ, as shipped. λ_Z = 0 is the Gertler–Karadi upper bound, reported as
sensitivity.

### 2.7 Summary: the first fix against the faithful fix

| Element | First fix | Faithful fix | Anchor |
|---|---|---|---|
| funding | national + profits returned to D (κ_D = 1) | **union / TARGET2 margin** | ELA, TARGET2; a facility at parity = the union market (§2.1) |
| default response | logit, φ = 1 | **CNT linear, θ_q = 0.0046** (4.5× weaker) | Costain–Nuño–Thomas θ = 2.02e-5/€bn |
| P&L sharing | κ_D = 1 | **κ_D = 0.071 (capital key)**; PSPP 0.814 as a sensitivity (numerically fragile) | SMP income shared; TPI/OMT rule not public; PSPP 80% national |
| creditor treatment | pari passu | pari passu | TPI press release |
| activation | fixed 200 bp cap | **fundamental spread + 50 bp** | "unwarranted, disorderly"; fundamentals test |
| reserves | non-divertible (λ_Z = 0) | **divertible (λ_Z = λ)**; λ_Z = 0 as the bound | Diamond–Jiang–Ma; leverage ratio |

---

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
- **The full-grid run of the faithful baseline (F1c) is the reference** (§3.4). It took 10.8 hours, and it
  revises the quick-grid magnitudes.

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

### 3.3 A realized default: what the insurance pays (quick grid, faithful baseline Q6)

Every impulse response so far conditions on no default, so it shows only the premium side of the backstop. Here
the headline shock hits at q0 and the default (the 55% haircut, pari passu for the Eurosystem) is realized at
q1. Each path is differenced against the same economy's unshocked path. Debt and fiscal flows are in % of D's
quarterly GDP.

| q1 response | No TPI: no default | No TPI: **default** | Faithful TPI: no default | Faithful TPI: **default** |
|---|---:|---:|---:|---:|
| D output | −0.00% | **−2.51%** | −1.53% | **+0.29%** |
| D consumption | +0.05% | −1.36% | −3.21% | +0.08% |
| D bank net worth | −2.93% | **−32.5%** | −7.2% | **+2.1%** |
| D debt stock | +1.4 | −54.3 | +0.7 | −52.7 |
| Eurosystem P&L | 0 | 0 | **+1.48** | **−30.9** |
| F's taxes (in D-GDP units) | 0 | 0 | **−1.37** (F gains) | **+28.7** (F pays) |
| Eurosystem book (% of the stock) | 0 | 0 | 84 | 34 (after the haircut) |

**Reading.**
- **The bargain.** In no-default quarters D pays a premium: the bonds' excess return accrues to the
  Eurosystem, 93% to F by the capital key, and D bank net worth and output are lower. In the default quarter D
  receives the claim: the Eurosystem takes the haircut on its book, F's taxpayers cover 93% of a loss worth 31%
  of D's quarterly GDP, and D's banks come through with net worth up instead of down a third.
- **The real-world counterpart.** This is the euro area's own bargain, with one difference. In Greece the
  Eurosystem avoided the loss by seniority, being exempted from the PSI, and returned its income to Greece
  afterwards. A pari-passu TPI with a capital key turns that into an explicit cross-border transfer in default.
  Which one a TPI implies is the political-economy question the ECB has not answered publicly.
- **A caution on dynamics.** After q0 these quick-grid paths are clipped at the edge of the TPI book's box (80%
  of the stock), so q1 onwards is indicative. On the full grid the headline shock barely reaches the warranted
  cap, and the TPI economy's dynamics are not yet identified (§3.4), so these magnitudes remain quick-grid ones.

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

---

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

---

## 5. What to take into the paper, and what to run next

**The faithful TPI, as a specification:**

1. **Funding through a TARGET2 margin.** Use the union deposit market, interpreted as Eurosystem
   intermediation of cross-border funding (§2.1). The TPI's results should not be reported on the national
   market: there the deposit rate insures the banker and the backstop has nothing to act on.
2. **A rollover component of default risk** à la Costain–Nuño–Thomas: π_t = p^d(s_t) − θ_q·M_{t+1}/B_{t+1}, with
   θ_q = 0.0046 mapped from their θ = 2.02×10⁻⁵ per € billion (§2.2). State in the text that p^d(s) is the
   fundamental part and θ_q the backstop-sensitive part.
3. **Pari passu**, as the TPI and OMT texts commit (§2.4).
4. **A fundamentals-based activation rule**: cap = fundamental (expected-loss) spread + 50 bp, the model
   counterpart of "unwarranted, disorderly" (§2.5). Drop the fixed 200 bp cap: below the expected-loss spread
   it is a transfer, not a backstop.
5. **P&L by capital key** as the baseline, PSPP-style national income retention as the sensitivity, and a
   statement that the TPI's actual sharing rule is not public (§2.3).
6. **Reserves as divertible as bonds** (λ_Z = λ), as the regulatory treatment and Diamond–Jiang–Ma suggest,
   with λ_Z = 0 as the Gertler–Karadi bound (§2.6).

**What the faithful TPI does in the model (§3):** it caps the spread at the fundamental level plus 50 bp,
lowers the risk premium (on the quick grid it extracts all of it), and in a realized default becomes a large
cross-border transfer from F to D.

**What it costs:** it needs a very large book, about 80% of the stock on the quick grid, because private
demand is too elastic. On the full grid the headline shock barely reaches the warranted cap, so the book stays
small; the slope, −1.2 bp per 1% of the stock, is still 12–25 times below the evidence.

**Next steps, in order:**

1. **The elasticity gap (§4).** Add CNT's two-type demand: inelastic preferred-habitat holders plus
   mean-variance arbitrageurs, with γ calibrated to the Eser–Schwaab slope, and the 2011–12 Greek spread level
   as a joint target. This is the single most important remaining step toward the real-world magnitudes.
2. **PSPP-style retention.** Narrow the TPI obligation band (`o_band`, currently ±25% of the book's half-width)
   so that the corners with ±20%-of-GDP P&L disappear, then re-run Q3b.
3. **The faithful design at 9 nodes in s.** At 115 points its announcement effect and dynamics are not
   identified (§3.4). Only the baseline (F1c) has been run on the full grid; the other variants have not.
4. **Implement D3** (the explicit Eurosystem facility with a collateral constraint and the ELA penalty) as the
   structural version of the TARGET2 margin.
5. **Welfare.** Weigh the no-default path (D pays the insurance premium) and the default path (F pays the
   claim) with the default probabilities. The two paths in §3.3 are the inputs.

---

## 6. Files

| File | Contents |
|---|---|
| `Claude files/docs/tpi_fix_explained.md` | this report |
| `Claude files/docs/tpi_national_analysis.md` | the previous note: why the shipped TPI fails under national deposits |
| `Claude files/diagnostics/tpi_national/tpi_switches_v3.patch` | every experimental switch as one patch against `domestic-funding-fix` (`0dd0732`): the national-TPI unlock, `lambda_Z_D`, `phi_float`, `theta_cnt`, the warranted cap, and the income/loss key split |
| `Claude files/diagnostics/tpi_national/tpi_national.py` | solve + compare driver (`--theta`, `--warranted`, `--keyD`, `--keyLoss`, `--lamZ`, `--market`, `--s-refine`) |
| `Claude files/diagnostics/tpi_national/default_path.py` | the realized-default paths of §3.3 |
| `Claude files/diagnostics/tpi_national/summarise.py` | the table numbers |
| `results/TPI_NATIONAL/` | every run's JSON and log (git-ignored) |

---

## 7. References

- Altavilla, C., D. Giannone and M. Lenza (2016). The financial and macroeconomic effects of OMT announcements. *International Journal of Central Banking* 12(3).
- Bocola, L. (2016). The pass-through of sovereign risk. *Journal of Political Economy* 124(4).
- Bocola, L. and A. Dovis (2019). Self-fulfilling debt crises: A quantitative analysis. *American Economic Review* 109(12).
- Calvo, G. (1988). Servicing the public debt: The role of expectations. *American Economic Review* 78(4).
- Cole, H. and T. Kehoe (2000). Self-fulfilling debt crises. *Review of Economic Studies* 67(1).
- Corsetti, G. and L. Dedola (2016). The mystery of the printing press. *Journal of the European Economic Association* 14(6).
- Costain, J., G. Nuño and C. Thomas (2025). The term structure of interest rates in a heterogeneous monetary union. *Journal of Finance* 80(4). Working-paper version, November 2021, used for the equations and θ.
- Diamond, W., Z. Jiang and Y. Ma (2024). The reserve supply channel of unconventional monetary policy. *Journal of Financial Economics* 159.
- Eser, F. and B. Schwaab (2016). Evaluating the impact of unconventional monetary policy measures: Empirical evidence from the ECB's Securities Markets Programme. *Journal of Financial Economics*.
- Gabaix, X. and M. Maggiori (2015). International liquidity and exchange rate dynamics. *Quarterly Journal of Economics* 130(3).
- Gertler, M. and P. Karadi (2011). A model of unconventional monetary policy. *Journal of Monetary Economics* 58(1).
- Krishnamurthy, A., S. Nagel and A. Vissing-Jorgensen (2018). ECB policies involving government bond purchases: Impact and channels. *Review of Finance* 22(1).
- Lorenzoni, G. and I. Werning (2019). Slow moving debt crises. *American Economic Review* 109(9).
- Vayanos, D. and J.-L. Vila (2021). A preferred-habitat model of the term structure of interest rates. *Econometrica* 89(1).

**Institutional and data sources consulted:**
- [ECB press release on the TPI, 21 July 2022](https://www.ecb.europa.eu/press/pr/date/2022/html/ecb.pr220721%7e973e6e7273.en.html) · [Bundesbank TPI page](https://bundesbank.de/en/tasks/monetary-policy/outright-transactions/transmission-protection-instrument-tpi--896158)
- [ECB, technical features of OMT, 6 Sep 2012](https://www.ecb.europa.eu/press/pr/date/2012/html/pr120906_1.en.html)
- [ECB FAQ on the PSPP](https://www.ecb.europa.eu/mopo/implement/app/html/ecb.faq_pspp.en.html) · [ECB FAQ on the PEPP](https://www.ecb.europa.eu/mopo/implement/pepp/html/ecb.faq_pepp.en.html) · [Eurosystem monetary income (Springer, 2025)](https://link.springer.com/article/10.1007/s11079-025-09804-w)
- [Eurogroup and the SMP income transfers to Greece: ECB correspondence](https://www.ecb.europa.eu/ecb/access_to_documents/document/correspondence/shared/data/ecb.dr.cor20171010Chountis.fi.pdf) · [Private sector involvement (Greek PSI, SMP exemption)](https://en.wikipedia.org/wiki/Private_sector_involvement) · [Xafa, "Life after debt"](https://mirandaxafa.com/mediaupload/articles/life_after_debt.pdf)
- [Greek ELA and Eurosystem funding: Yale case study](https://elischolar.library.yale.edu/journal-of-financial-crises/vol4/iss2/39/) · [TARGET2 balances 2011: CESifo Forum](https://www.cesifo.org/DocDL/forum-0112-special-5.pdf) · [Greek deposit rates, April 2011 (Bank of Greece data)](https://neoskosmos.com/en/2011/06/15/news/greece/greek-banks-deposits-down-2-4-billion-euros-april/)
- [Costain–Nuño–Thomas, BIS WP 1165](https://www.bis.org/publ/work1165.pdf) · [Diamond–Jiang–Ma, NBER w31693](https://www.nber.org/papers/w31693) · [Altavilla–Giannone–Lenza, ECB WP 1707](https://ECB.eu/pub/pdf/scpwps/ecbwp1707.pdf) · [Krishnamurthy–Nagel–Vissing-Jorgensen, NBER w23985](https://www.nber.org/papers/w23985) · [Eser–Schwaab, EconStor](https://econstor.eu/handle/10419/154075)
