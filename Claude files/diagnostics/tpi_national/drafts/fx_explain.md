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
