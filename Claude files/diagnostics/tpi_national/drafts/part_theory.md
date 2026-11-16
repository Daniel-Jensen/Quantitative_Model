## 3. Theory: what a TPI purchase can and cannot move in this model

### 3.1 The pricing equation and its three wedges

Everything the TPI can do to the D spread has to pass through the D bank's first-order condition
for the D bond, which sets the price q (D-bond price) whenever the bank holds a positive amount:

$$
q_t\,\big(\mathbb{E}_t[\Omega_{t+1}]\,R_t+\lambda\mu_t\big)=\mathbb{E}_t\big[\Omega_{t+1}\,\Xi_{t+1}\big],
\qquad
\Xi_{t+1}=h_{t+1}\big(\delta+(1-\delta)\,q_{t+1}\big).
\tag{P}
$$

Here:
- Ω is the banker's stochastic discount factor, Ω = Λ·(f + (1−f)α′), with Λ the household-type discount
  factor, f the banker exit share and α′ tomorrow's franchise value per unit of net worth;
- R is the gross D deposit rate, λ the divertible fraction of assets and μ the multiplier on the
  incentive constraint;
- Ξ is the bond's payoff, with h = 1 without default and h = 0.45 (the recovery rate) in default, and δ the
  coupon/decay rate of the perpetuity.

Splitting E[ΩΞ] = E[Ω]E[Ξ] + Cov(Ω, Ξ) turns (P) into a statement about the one-period expected excess return:

$$
\underbrace{\frac{\mathbb{E}_t[\Xi_{t+1}]}{q_t}-R_t}_{\text{expected excess return}}
=\underbrace{\frac{\lambda\mu_t}{\mathbb{E}_t[\Omega_{t+1}]}}_{\text{liquidity premium}}
\;\underbrace{-\;\frac{\mathrm{Cov}_t(\Omega_{t+1},\Xi_{t+1})}{q_t\,\mathbb{E}_t[\Omega_{t+1}]}}_{\text{risk premium}} .
$$

The expected default loss sits inside E[Ξ]. With a default probability π and recovery h it is
π(1−h)·(δ + (1−δ)q′_nd)/q per quarter. Against the F bond, which the F bank prices with the same
equation and never defaults, the D–F yield spread is approximately

$$
y^D-y^F\;\approx\;\underbrace{\pi(1-h)}_{\text{expected loss (exogenous)}}
+\underbrace{(\mathrm{rp}^D-\mathrm{rp}^F)}_{\text{risk premia}}
+\underbrace{\Big(\tfrac{\lambda\mu_D}{\mathbb{E}[\Omega_D]}-\tfrac{\lambda\mu_F}{\mathbb{E}[\Omega_F]}\Big)}_{\text{liquidity premia}}
+\underbrace{(r_D-r_F)}_{\text{funding-rate gap}},
$$

accumulated over the bond's life through q′. That gives a full list of what a central-bank purchase *could*
move:

| Spread component | What a TPI would have to change | Can the current TPI change it? |
|---|---|---|
| expected loss π(1−h) | the default probability π or the recovery h | **No.** π is Bocola's exogenous s-process and h is fixed. The TPI buys nothing in default and does not alter the haircut. |
| risk premium | the covariance of Ω with the payoff, i.e. who bears the loss in default | **Only through tomorrow's net worth** (the risk transfer, §3.3). |
| liquidity premium | μ, the tightness of the bank's constraint | **No, within the period**: the purchase is paid with a claim Z at the bond's own λ (§3.2). Only indirectly, through E[Ω]. |
| funding-rate gap | the D deposit rate relative to F | Only through general equilibrium; under national deposits r_D is set by D's own market. |
| continuation price q′ | all of the above in future states | **Yes, through the floor:** future no-default states where the cap binds get a price floor (§3.5). |

### 3.2 Proposition 1: a TPI purchase is neutral within the period

A purchase dm at price q is paid with Z = q·dm, a one-period claim on the Eurosystem held by the
D bank, paying R next period and divertible at λ_Z = λ (`point_map.py:311`). Within the period, the
bank's assets 𝒜, its divertible base λ𝒜, its deposits, its net worth n and therefore

$$
\mu_t=1-\frac{\mathbb{E}_t[\Omega_{t+1}]R_t\,n_t}{\lambda\mathcal{A}_t}
$$

are all unchanged. Holding tomorrow's kernel Ω′ and price q′ fixed, (P) is unchanged too, so ∂q/∂m = 0. Any
price effect must therefore come from the response of Ω′ and q′ to the next-period state the purchase leaves
behind (the book M′, the obligation O′, the bank's bond book b′_DD). That response is Proposition 2.

This is the open-market-operation irrelevance of Wallace (1981) and Eggertsson & Woodford (2003), restated for an
intermediary. If the central bank swaps one asset for another of equal market value, and the swap changes no
agent's constraint set, prices do not move. The project's TPI was designed this way on purpose: it is
"a swap today and a risk transfer tomorrow" (`tpi_backstop_plan.md` §4), so as not to repeat the yield peg,
which could only act on the liquidity premium. The cost of that choice is that **the channel through which
central-bank purchases are powerful in Gertler–Karadi models is switched off by construction.**

### 3.3 Proposition 2: the risk-transfer price impact

Tomorrow, in regime j (a node of the s′ × d′ quadrature, with probability w_j), the bank is owed R·q·dm
instead of Ξ_j·dm. Its net worth changes by

$$
dn'_j=(1-f+\omega\kappa_D)\,(Rq-\Xi_j)\,dm ,
$$

where ω is the entrant transfer and κ_D the D share of the Eurosystem's profit and loss: positive in default, slightly negative
without default. That moves tomorrow's kernel through the franchise value. Where the constraint binds tomorrow,
Bocola's closed form gives

$$
\alpha=\frac{\mathbb{E}[\Omega'']R}{1-\mu},\quad \mu=1-\frac{\mathbb{E}[\Omega'']Rn}{\lambda\mathcal{A}}
\;\;\Longrightarrow\;\;
\frac{\partial\alpha}{\partial n}=-\frac{\alpha^2}{\lambda\mathcal{A}}\quad(\mu>0),\qquad
\frac{\partial\alpha}{\partial n}=0\quad(\mu=0).
$$

So dΩ_j = k_j·q·e_j·dm, where e_j ≡ Ξ_j/q − R is the bond's excess return in regime j, and

$$
k_j\;\equiv\;\Lambda_j(1-f)(1-f+\omega\kappa_D)\,\frac{\alpha_j'^2}{\lambda\mathcal{A}'_j}\,\mathbf{1}\{\mu'_j>0\}\;\ge 0 .
$$

Differentiating (P), and using dμ = −(1−μ)·dE[Ω]/E[Ω] (n and 𝒜 do not move), gives at fixed rules:

$$
\boxed{\;
\frac{d\ln q}{dm}
=\frac{q}{\mathbb{E}[\Omega]R+\lambda\mu}\Big[\underbrace{\textstyle\sum_j w_j k_j e_j^2}_{\text{risk transfer}\;\ge 0}
+\underbrace{\tfrac{\lambda(1-\mu)}{\mathbb{E}[\Omega]}\textstyle\sum_j w_j k_j e_j}_{\text{franchise channel}}\Big]
+\underbrace{\frac{\sum_j w_j\Omega_j\,d\Xi_j/dm}{\mathbb{E}[\Omega\Xi]}}_{\text{continuation (smaller float)}}
\;}
\tag{I}
$$

Formula (I) says three things.

1. **The power of a purchase is a kernel-weighted variance.** The leading term is the variance of the bond's
   excess return, weighted by k_j, the sensitivity of the banker's marginal value of wealth to net worth in
   each regime. A purchase moves the price only to the extent that the default state is a state in which
   banker net worth is scarce. This is the GK-model analogue of the "risk extraction" channel in
   preferred-habitat models (Vayanos & Vila 2021; Costain, Nuño & Thomas 2025), where the impact is risk
   aversion × variance.
2. **The franchise term can offset it.** If k is much larger in default than without default, then Σw k e < 0:
   a safer future lowers E[Ω], which raises today's μ and *lowers* the price (`tpi_backstop_plan.md` §4
   anticipated this).
3. **The impact is zero where tomorrow's constraint is slack (k_j = 0) and small where Ω_j does not respond.**

The closed form for k_j keeps only the net-worth effect. In the solved model the purchase also moves every other
next-period state it touches (the book M′, the obligation O′, the bank's bond book b′_DD), and with them Λ_j and
α′_j. So the runs measure dΩ_j directly from the solved rules and evaluate (I) with the measured dΩ_j. The
formula reproduces the directly implied price change to within 0.005 pp on the no-TPI rules of all three economies, at rest and at the shock (§2.4).

### 3.4 Why national deposit funding reverses the sign: measured

Formula (I) leaves the sign of the impact to two equilibrium objects: how Ω_j responds to the purchase in each
regime, and whether the default state is a high-Ω state. Both differ between the two deposit markets, and the
runs in §2 measure them.

**(a) The default state is not a bad state for the D banker under national funding.** The kernel is
Ω_j = Λ_j(f + (1−f)α′_j), with Λ_j = β(x_t/x_{t+1,j})^σ built from D's GHH consumption composite x. Under
national deposits the D bank cannot fund itself abroad, so in the default branch D's deposit market must clear
at home while the bank deleverages. The next-period D deposit rate in default is **−860 bp/yr at rest and
−710 at the shock** (union: −417 and −361). D's captive depositors insure the D banker: they accept a deeply
negative real return exactly when the bank's net worth collapses. That cheapens the banker's future funding, and
with it α′ (tomorrow's franchise value). In addition, D's consumption composite rises in that branch (by about 9%
on the coarse grid, Note 2), which lowers Λ^d. Both push Ω^d down. The result:

| | Ω^d/Ω^nd, rest | Ω^d/Ω^nd, shock | risk premium, rest | risk premium, shock |
|---|---:|---:|---:|---:|
| national | **0.97** | **0.95** | +0.7 bp/yr | **−16.4 bp/yr** |
| union | 1.29 | 1.25 | +8.4 bp/yr | **+118.9 bp/yr** |

Under national funding Cov(Ω, Ξ) > 0 at the shock, so the D bond is a *hedge* for the D banker. The risk
premium that the TPI's risk transfer is designed to remove does not exist; it is negative.

**(b) The purchase reaches the D household's marginal utility.** Per 5% of the stock bought (fixed rules, at
rest):

| | dΛ/Λ (no default / default) | dα′/α′ (no default / default) | dr′_D (no default / default) |
|---|---|---|---|
| national | **−0.51% / −0.86%** | −0.88% / −4.96% | +7 / **−66 bp** |
| union | +0.04% / −0.06% | −0.00% / −5.00% | +1 / +28 bp |

Under union funding the purchase does what §3.3 intends: the default-state franchise value falls (the bank is
insured), and nothing else moves. Under national funding the change in the D bank's future funding needs must
be absorbed by D households through D's own deposit rate, so their consumption composite and Λ move. The
banker's kernel falls in *both* regimes; E[Ω] falls; today's μ = 1 − E[Ω]Rn/(λ𝒜) rises; and the
**franchise term of (I) dominates with the wrong sign**: −0.267% of the price per 5% of the stock at rest,
against a risk-transfer term of −0.001% and a continuation term of +0.173%.

**(c) Making the Eurosystem's claim non-divertible is not enough under national funding.** With λ_Z = 0
(reserves outside the incentive constraint, the Gertler–Karadi channel) the direct relaxation of μ is
switched on. The fixed-rule effect of a 30% purchase at the shock falls from +10.2 bp to +4.3 bp on the
sovereign spread, but it is still the wrong sign: the franchise and household channels are larger than the
balance-sheet relief.

### 3.5 Proposition 3: the cap cannot be held where the price impact is negative

The cap is a price floor, q̲ = δ/(δ + y_F + c), where y_F is the F bond's yield and c the cap. At a node where
the no-TPI price is below q̲, the Eurosystem must buy an m* such that the bank's valuation
v(m) ≡ E[ΩΞ]/(E[Ω]R + λμ) reaches the floor: v(m*) = q̲. Under national funding v′(m) < 0 (§3.4), so buying
moves the bank's valuation *away* from the floor and **no interior m* exists**. In the national runs the walk
holds 493.75 bp and fails at 487.5 bp on the TPI equation itself (max|F| ≈ 6e-4 to 2e-2), at the one node where
the floor first binds; the D bank there still holds about 5% of the stock (min b_DD/B = 0.048), so this is not
a corner.

The only remaining equilibrium is the corner in which the bank sells everything and the Eurosystem holds the
float. Under national deposits the F banks' D-bond books are also no-short-sale pairs, so they exit too, and
the corner is M ≈ B. That needs a book larger than the state box allows (±b_DD,ss = ±0.784, 80% of the stock);
this is consistent with the commit note that the refined national grid stalled at full-book corners.

Even where the corner exists, the instrument changes nature below the expected-loss spread. At the headline
shock the one-period expected loss alone is π(1−h) = 0.02 × 0.55 = **440 bp/yr** (452 bp in the national
run). A cap below the risk-neutral spread means the Eurosystem pays more than the bond is worth to a
risk-neutral holder. It is then financing the D government at a subsidised yield, with expected losses that
the capital key sends 92.9% to F. `tpi_backstop_plan.md` §8 says it directly: "Below the expected-loss spread
the instrument is a cross-border transfer, not a backstop."

### 3.6 Proposition 4: the announcement is a put struck in remote states

At the rest point nothing is bought. The TPI raises today's price only because it floors *future* prices in
states where the cap would bind:

$$
q_t^{\rm TPI}-q_t^{\rm noTPI}\;\approx\;\mathbb{E}_t\Big[\sum_{k\ge1}\tilde\Omega_{t,t+k}\,(1-\delta)^k
\,\mathbf{1}\{\text{no default to }t+k\}\,\big(\underline{q}_{t+k}-q^{\rm noTPI}_{t+k}\big)^+\Big]
+\text{(liquidity and risk-premium feedback)} .
$$

Here Ω̃_{t,t+k} is the bank's kernel compounded from t to t+k, and the last term collects the feedback on μ and
on the risk premium. The put has two features that make it small:
- **It is struck in remote states.** On the full union grid the cap binds from p^d ≈ 1.6%/qtr (20 of
  the 23 nodes at 1.58%/qtr, every node at 4.8%/qtr). From the rest point (0.10%/qtr), reaching that in one
  quarter takes Δs = logit(0.0158) − logit(0.001) ≈ 2.8, which is **6.2 standard deviations** of the
  s-innovation (σ_s = 0.4455).
- **It pays nothing in default.** The Eurosystem does not buy in default and recovery is fixed, so the floor
  never insures the 55% haircut, which is the component that dominates the spread.

That is why, even under union funding on the full grid, the announcement moves the *sovereign* spread by only
**−4.3 bp** (55.8 → 51.6) while it moves the *lending* spread by −36 bp: what the announcement mostly does is
lower μ through the franchise channel, not insure bondholders (`tpi_backstop_plan.md` §9.2).

The draft's alternative design, the "disorderly market" TPI (`model_tracked_revision.tex` eq. q-warranted),
targets the default-warranted price q⋆. That is the price a risk-neutral bank would pay, given the same deposit
rate and liquidity premium. It is always feasible and never a transfer in expectation, but by construction it
removes **only** the risk premium. Under national deposits, where that premium is about zero or negative (§2.4), it would do
nothing either.
