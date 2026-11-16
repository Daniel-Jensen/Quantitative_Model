# Why the TPI cannot lower the D spread under national deposit funding, and how to make it

*Working note, 12 October 2026. Global-projection model, code of branch `domestic-funding-fix` (`0dd0732`) in a
detached worktree, with three experimental switches that are off by default (§1.2). Every run is on the
25-point quick grid: good for signs and for comparisons between economies solved the same way, not for levels.
Full-grid union numbers are quoted from `tpi_backstop_plan.md` §9.2.*

---

## 0. Bottom line

1. **As shipped, the TPI does nothing under national deposits.**
   - The lowest cap the solver can hold is **493.75 bp**.
   - The announcement moves the sovereign spread by +0.1 bp and the lending spread by 0.
   - A forced purchase of 30% of the stock **raises** the sovereign spread by 10 bp at the headline shock.
   - This is not a box artefact: the run uses the TPI book's full box, not the 0.1% band the branch ships for
     national mode.

2. **The cause is a sign, and the formula in §3.3 makes it measurable.** A TPI purchase swaps the bond for a
   Eurosystem claim Z of equal value and the same divertibility. It is neutral within the period
   (Proposition 1) and moves the price only through tomorrow's banker kernel Ω (formula (I)).
   - **Union funding:** Ω falls in default and rises without default. That is the intended risk transfer, and
     the price rises.
   - **National funding:** Ω falls in *both* regimes. D's deposit market clears at home, so a purchase reaches
     D households' marginal utility through D's own deposit rate (−66 bp in default per 5% of the stock).
   - **No risk premium to remove:** the default state is not a bad state for the D banker. His kernel is lower
     there (ratio 0.95), because captive depositors take a −710 bp/yr real deposit rate. The D bond is a hedge,
     with a **−16 bp risk premium** against **+119 bp** under union.
   - **Net effect:** the constraint (franchise) term of (I) dominates, a purchase tightens today's constraint,
     and the price falls.

3. **The announcement cannot reach the dominant term.** Most of the spread is expected default loss (452 bp/yr
   one-period at the headline shock), and the Bocola s-process makes it exogenous. The cap is a put struck in
   states six standard deviations away, and it pays nothing in default. Even under union on the full grid the
   announcement moves the sovereign spread by only −4.3 bp.

4. **The literature gets large effects from two channels this model switches off:**
   - **balance-sheet capacity (Gertler–Karadi):** Krishnamurthy–Nagel–Vissing-Jorgensen attribute 50% of the
     SMP/OMT yield decline to it ("segmentation");
   - **default risk that responds to the backstop:** rollover and self-fulfilling crises (Calvo; Cole–Kehoe;
     Lorenzoni–Werning; Corsetti–Dedola; Camous–Cooper) and default risk extraction in Costain–Nuño–Thomas
     (2025, *JF*). This accounts for 37% of the decline.

5. **The concrete fix has three parts, and each was tested.**
   - **(i) Non-divertible Eurosystem claims** (λ_Z = 0).
   - **(ii) Default risk extraction:** logit π = s + φ·ln((B − M)/B), at φ = 1.
   - **(iii) Remit the Eurosystem's P&L to the issuer** (κ_D = 1). This part is specific to national funding:
     under the capital key, the no-default profits on a large book go 92.9% to F, and with no external funding
     margin D's deposit rate explodes (1,158–1,215 bp/yr two quarters after the shock).

   | Economy (national deposits) | Lowest cap held | Announcement effect on the sovereign spread |
   |---|---:|---:|
   | (i) alone | 494 bp | 0 |
   | (ii) alone | 406 bp | wrong-signed |
   | (i) + (ii) | 375 bp | −12 bp |
   | **(i) + (ii) + (iii)** | **200 bp** (the target) | **−23 bp** |

   With all three, the backstop **halves the spread at the headline shock (339 → 166 bp) without buying
   anything**, and turns the national recession (cumulative output −6.0 over 25 quarters) into no recession
   (+1.9).

   A structural alternative to (iii) is the Eurosystem-refinancing funding closure of the companion note (D3).
   It gives D an external funding margin and makes default a bad state again. It is not tested here.

---

## 1. What was run

### 1.1 Setup

| Run | Deposit market | Eurosystem claim Z | Default probability | P&L remittance | Purpose |
|---|---|---|---|---|---|
| `nat` | national | divertible at λ (as shipped) | exogenous p^d(s) | capital key (7.1% D) | the question asked |
| `uni` | union | divertible at λ | exogenous | capital key | control: the economy the TPI was built for |
| `nat_lz0` | national | **non-divertible** (λ_Z = 0) | exogenous | capital key | (i) capacity |
| `nat_phi1` | national | divertible | **logit π = s + ln((B−M)/B)** | capital key | (ii) default risk extraction |
| `nat_phi1_lz0` | national | non-divertible | logit π = s + ln((B−M)/B) | capital key | (i) + (ii) |
| `nat_phi1_lz0_key1` | national | non-divertible | logit π = s + ln((B−M)/B) | **all to D** (κ_D = 1) | (i) + (ii) + (iii) |
| `uni_phi1_lz0` | union | non-divertible | logit π = s + ln((B−M)/B) | capital key | (i) + (ii) with a funding margin |

Each run solves the no-TPI economy (Bocola's ladder), then walks the cap from 700 bp towards 200 bp in 50 bp
rungs (halving a failed rung down to 6.25 bp), then compares the two economies:
- the rest points (the announcement);
- the headline-shock (p^d → 2%/qtr) impulse responses, cleared exactly each quarter;
- the price impact of a forced purchase at fixed rules (5–30% of the stock);
- the complementarity conditions;
- the D-bond price split into expected loss, risk premium and liquidity premium, with the banker kernel by regime.

### 1.2 Code changes (worktree only; nothing on any repo branch)

| File | Change | At its default |
|---|---|---|
| `recursive_experiment.py` | lifted the `assert` that forbids the TPI under national deposits | — |
| `recursive_experiment.py` | the 0.1% `NATIONAL_M_BAND` applies only if `cal["national_m_band"]` (set False here) | shipped behaviour |
| `point_map.py:311` | Z enters the divertible base at `cal.get("lambda_Z_D", lambda_bD_D)` | identical |
| `point_map._continuation` | `p^d(s + phi_float·ln(max((B′−M′)/B′, 0.05)))`, with a smooth floor | identical (φ = 0) |
| calibration key | `tpi_key_D` (D's share of the Eurosystem P&L) set to 1 in `nat_phi1_lz0_key1` | 0.071 |

**Files.**
- **The four switches as a patch** against `domestic-funding-fix`:
  `Claude files/diagnostics/tpi_national/tpi_switches.patch`. Apply it with `git apply` on a checkout of that
  branch.
- **The drivers:** `Claude files/diagnostics/tpi_national/tpi_national.py` (solve, compare, price legs) and
  `tpi_mechanism.py` (formula (I), measured).
  - Point them at the patched checkout with `WT=/abs/path/to/checkout`.
  - Example: `WT=… python3 tpi_national.py --market national --phi 1 --lamZ 0 --keyD 1 --tag mine`.
  - Add `--s-refine 5` for the full grid.
- **Every run's data and log:** `results/TPI_NATIONAL/` (git-ignored). This holds `<tag>.json` (the full
  comparison) and `mechanism*.json` (the formula-(I) decompositions).

---

## 2. Results

### 2.1 Can the cap be held, and does the announcement work?

| Run | Lowest cap held | Announcement, sovereign spread (bp/yr) | Announcement, lending spread (bp/yr) | Output at the rest point, TPI vs none |
|---|---:|---|---|---:|
| `nat` | 493.75 bp | 53.4 → 53.5 | 90.0 → 90.0 | 0 |
| `nat_lz0` | 493.75 bp | 53.4 → 53.5 | 90.0 → 90.0 | 0 |
| `nat_phi1` | 406.25 bp | 53.4 → 57.7 | 90.0 → 9.6 | −0.44% |
| `nat_phi1_lz0` | 375 bp | 53.4 → 41.4 | 90.0 → 16.2 | −0.32% |
| **`nat_phi1_lz0_key1`** | **200 bp** ✓ | **53.4 → 30.2** | **90.0 → 42.7** | **+0.46%** |
| `uni` | 200 bp ✓ | 62.5 → 13.6 | 94.8 → 0.0 | +1.96% |
| `uni_phi1_lz0` | 200 bp ✓ | 62.5 → 29.2 | 94.8 → 0.0 | +0.59% |

- **The national failures.** In `nat` and `nat_lz0` the walk fails on the TPI equation itself (max|F| 6e-4 to
  2e-2) at the one node where the floor first binds. The D bank still holds about 5% of the stock there
  (min b_DD/B = 0.048), so it is not a corner.
- **The φ runs under the capital key.** In `nat_phi1` and `nat_phi1_lz0` the announcement crushes the lending
  spread but *lowers* output at rest. Relaxing the D bank's constraint under national funding raises D's
  deposit rate (132 → 165 bp at rest in `nat_phi1`) and crowds out investment, the national S = I constraint
  again.
- **The coarse union numbers are inflated.** They overstate the full-grid union announcement several-fold
  (full grid: −4.3 bp sovereign, −36 bp lending), so compare within the quick grid only.

### 2.2 Purchases: price impact at fixed rules (change in the sovereign spread)

| Purchase (share of the stock) | `nat` rest | `nat` shock | `nat_lz0` shock | `nat_phi1` shock | `nat_phi1_lz0` shock | **`…_key1` rest** | **`…_key1` shock** | `uni` shock | `uni_phi1_lz0` shock |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5% | +1.8 | +2.1 | +1.3 | −1.7 | −2.7 | −2.4 | −5.5 | −2.4 | −7.7 |
| 10% | +3.6 | +4.1 | +2.3 | −2.9 | −5.0 | −4.6 | −10.9 | −4.6 | −15.2 |
| 20% | +7.0 | +7.6 | +3.6 | −4.3 | −8.1 | −8.8 | −21.3 | −8.4 | −29.5 |
| 30% | **+9.8** | **+10.2** | +4.3 | −4.3 | −9.4 | **−12.4** | **−29.1** | −11.5 | **−43.0** |

All values are in bp/yr and use each run's solved TPI-economy rules. For scale, Eser & Schwaab (2016) estimate
about −3 bp on SMP five-year yields per purchase of 0.1% of the outstanding stock, about −30 bp per 1%. The best
model economy here gives about −1.0 to −1.4 bp per 1% at the shock: still one order of magnitude short at fixed
rules. The equilibrium effect is much larger, because the solved TPI economy prices the bond differently in
every state (§2.5). The model spread is the HM perpetuity yield and the SMP figure a five-year yield, so read
the comparison as an order of magnitude.

### 2.3 Why: the price legs and the kernel

**D-bond price legs (one-period premia, bp/yr), no-TPI rules:**

| | expected loss | risk premium | liquidity premium | Ω^d/Ω^nd | r_D today | r′_D, no default | r′_D, default |
|---|---:|---:|---:|---:|---:|---:|---:|
| national, rest | 22.6 | +0.7 | 81.5 | 0.971 | 132 | 124 | **−860** |
| national, shock | 451.7 | **−16.4** | 154.8 | **0.952** | 96 | 105 | **−710** |
| union, rest | 21.2 | +8.4 | 93.5 | 1.290 | 124 | 122 | −417 |
| union, shock | 428.1 | **+118.9** | 171.5 | **1.254** | 82 | 89 | −361 |

**Default risk extraction at work.** In `nat_phi1`, at the shock, the TPI economy buys 79% of the stock on
impact, and the one-period expected loss falls from **451.7 to 115.3 bp/yr**: the default odds fall with the
float.

### 2.4 Formula (I), measured

Fixed rules, 5% of the stock, no-TPI rules; changes in % of the price:

| | risk transfer | constraint (franchise) | continuation | (I) total | implied | fully re-solved |
|---|---:|---:|---:|---:|---:|---:|
| national, rest | −0.001 | **−0.267** | +0.173 | −0.094 | −0.097 | −0.086 |
| national, shock | +0.036 | −0.052 | −0.027 | −0.043 | −0.043 | −0.082 |
| national λ_Z=0, rest | −0.001 | −0.315 | +0.177 | −0.139 | −0.144 | −0.047 |
| national λ_Z=0, shock | +0.034 | −0.098 | −0.024 | −0.088 | −0.089 | −0.053 |
| union, rest | +0.004 | +0.006 | +0.040 | +0.049 | +0.049 | +0.054 |
| union, shock | **+0.061** | **+0.081** | −0.017 | **+0.125** | +0.124 | +0.109 |

"Implied" is the price at the perturbed continuation with today's other unknowns held; "fully re-solved" also
re-clears the other 13 unknowns. (I) matches "implied" to within 0.005 pp everywhere.

**The kernel response, split** (5% of the stock, fixed rules, rest point):

| | dΩ/Ω: no default / default | dΛ/Λ (household part) | dα′/α′ (franchise part) | dr′_D (bp) |
|---|---|---|---|---|
| national | −1.33% / −5.45% | **−0.51% / −0.86%** | −0.88% / −4.96% | +7 / **−66** |
| union | +0.03% / −4.80% | +0.04% / −0.06% | −0.00% / −5.00% | +1 / +28 |

At the shock: national −0.20% / −3.84% (household part −0.35% / −0.69%); union +0.53% / −3.76% (household part
+0.08% / −0.02%).

### 2.5 The headline shock: impulse responses (quick grid)

| Economy | Bought on impact | Sovereign spread on impact | Impact output | Trough output | Σ output, 25 q | Σ investment, 25 q | Bank net worth on impact | Box escapes |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| national, no TPI | — | 339.5 bp | −0.030% | −0.306% (q10) | −5.99 | −12.00 | −4.48% | none |
| `nat` / `nat_lz0` with TPI | 0 (cap 494 never binds) | 339.5 bp | identical | | | | | none |
| `nat_phi1` with TPI | **79% of the stock** | 406 bp (at cap) | −0.71% | −3.93% (q2) | −76.2 | −269 | −30.5% | **yes, M at the box edge** |
| `nat_phi1_lz0` with TPI | 79% | 375 bp (at cap) | −3.13% | −4.60% (q2) | −56.7 | −211 | −36.4% | **yes** |
| **`nat_phi1_lz0_key1` with TPI** | **0** | **166.0 bp** | **−0.011%** | **−0.011% (q0)** | **+1.91** | **+20.87** | **+0.66%** | none |
| union, no TPI | — | 414.4 bp | −0.077% | −0.077% (q0) | +1.22 | +5.71 | −5.03% | none |
| `uni` / `uni_phi1_lz0` with TPI | 0 | 126.7 / 164.7 bp | −1.78% / −1.41% | | | | | none (coarse pathology) |

- **The capital-key φ runs.** In `nat_phi1` and `nat_phi1_lz0` the Eurosystem must buy 79% of the stock to hold
  a cap of 375–406 bp. It then remits 92.9% of the book's no-default profit to F, and with no external margin
  D's deposit rate reaches **1,158–1,215 bp/yr in quarter 2**. Investment collapses on impact (−14% to −17%).
  After quarter 0 those paths are clipped at the box edge, so only the impact quarter is clean.
- **The full fix.** With the P&L remitted to D (`key1`) the same economy never needs to buy. Along the
  shocked path the spread stays below the cap, and the announcement does all the work.
- **One counter-intuitive feature.** In the full-fix economy the lending spread rises *more* on impact
  (+267 bp on exact reads, against +71 bp without the TPI), because a backstopped future lowers the franchise
  value and tightens the constraint in stress. Bank net worth and investment still rise, because the bond price
  falls only half as much (−5.3% against −10.3%).
- **The union TPI impulse responses on the coarse grid are pathological** in both union economies, as already
  recorded for the shipped union TPI.

---
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
removes **only** the risk premium. Under national deposits, where that premium is about zero or negative (§2.3), it would do
nothing either.

---

## 4. How the literature gets central-bank purchases to lower spreads

Five mechanisms appear in the literature. Each maps onto one term of the spread decomposition in §3.1, so
each can be checked against this model.

### M1. Balance-sheet capacity: the central bank is not constrained (Gertler–Karadi)

**Papers.** Gertler & Karadi (2011, *JME*) and Gertler & Kiyotaki (2010, *Handbook of Monetary Economics*);
Gertler & Karadi (2013, *IJCB*) for large-scale asset purchases; Dedola, Karadi & Lombardo (2013, *JME*) for
the two-country case.

**Mechanism.** Private intermediaries face the same incentive constraint as here. The central bank is less
efficient at intermediating (a resource cost τ per unit), but it can raise funds elastically by issuing
riskless government debt or reserves that *households* hold, or that do not count against bank capacity. A
purchase therefore *removes* the asset from the constrained balance sheet:

$$
\lambda\mathcal{A}^{\rm private}_t=\lambda\,(Q_tK_t+q_tb_t-q_tM_t)\quad\Longrightarrow\quad
\frac{\partial\mu_t}{\partial M_t}=-(1-\mu_t)\frac{q_t}{\mathcal{A}_t}<0 .
$$

That lowers the excess return on *every* asset the bank holds: the lending spread and the bond's
liquidity premium. This is what Krishnamurthy, Nagel & Vissing-Jorgensen (2018) call the **market-segmentation**
channel, and it accounts for **50%** of the yield declines they attribute to SMP and OMT for Italy, Spain and
Portugal.

**In this model.** It is switched off because Z enters the divertible base at λ_Z = λ. The concrete
implementation is one parameter: **λ_Z = 0**, reserves that cannot be diverted. They sit at the central bank;
Gertler & Kiyotaki's interbank parameter ω plays the same role for claims that are easier to enforce. It does
not reopen the portfolio-substitution margin that the single-λ rule guards against, because the bank does not
choose Z: the Eurosystem does.

### M2. Limited risk-bearing capacity: risk extraction (Vayanos–Vila; Costain–Nuño–Thomas)

**Papers.** Vayanos & Vila (2021, *Econometrica*) for preferred-habitat arbitrage; Gabaix & Maggiori (2015,
*QJE*) for limited intermediary risk capacity; and, for exactly this question, **Costain, Nuño & Thomas
(2025, *Journal of Finance*), "The Term Structure of Interest Rates in a Heterogeneous Monetary Union."**

**Mechanism.** Risk-averse arbitrageurs price bonds with a market price of risk proportional to the risk they
hold. Central-bank purchases reduce the duration and default risk left in private hands, so they compress the
term premium and the **credit risk premium**. Costain, Nuño & Thomas build a two-country monetary union with one
defaultable issuer, solve it in affine form, and find that both the level and the shifts of the Italian spread
are mainly the credit risk premium. They report that the PEPP announcement lowered Italian yields by about
80 bp, with its flexibility adding about 15 bp. Their extension makes the default probability fall as the
central bank holds more of the debt, which they call **default risk extraction** and justify by rollover crises.

**In this model.** Formula (I) is the GK analogue: the "price of risk" is the k-weighted variance. Under
national deposits Ω^d ≲ Ω^nd, so the price of risk is about zero and there is nothing to extract (§2.3). The
concrete implementations are either (a) a funding closure that keeps Ω^d > Ω^nd (D3, the Eurosystem refinancing
of the companion note, or union parity), or (b) a risk budget on bank sovereign exposure (P2 of the companion
note), which adds a genuine risk-bearing limit.

### M3. Default risk that the central bank can change: multiple equilibria and rollover

**Papers.** Calvo (1988); Cole & Kehoe (2000, *REStud*); Lorenzoni & Werning (2019, *AER*), "Slow Moving Debt
Crises"; Corsetti & Dedola (2016, *JEEA*), "The Mystery of the Printing Press"; Camous & Cooper (2019, *AEJ:
Macro*), "'Whatever It Takes' Is All You Need"; Aguiar, Amador, Farhi & Gopinath (2015, *QJE*); Bacchetta,
Perazzi & van Wincoop (2018, *JIE*). The quantitative counterpart is Bocola & Dovis (2019, *AER*).

**Mechanism.** Default is partly self-fulfilling. High yields raise the debt burden and make default more likely
(slow-moving crises), or investors refuse to roll over (rollover crises). A credible backstop removes the bad
equilibrium:
- Corsetti & Dedola show that above a critical scale of central-bank purchases the government strictly prefers
  fiscal adjustment to default when fundamentals are sound;
- Camous & Cooper show that a state-contingent commitment overcomes the strategic uncertainty, while simple rules
  do not;
- Bacchetta et al. show that monetary policy can preclude a slow-moving crisis only through sustained inflation.

**This is the only mechanism in which an announcement lowers the spread at zero purchases.** The OMT was never
used. Krishnamurthy et al. attribute **37%** of the SMP/OMT yield declines to lower default risk and **13%** to
lower redenomination risk. Bocola & Dovis find that non-fundamental risk played a *limited* role in Italy in
2008–12, so the size of this channel is an empirical question, not a given.

**In this model.** It is absent by construction: π_t = p^d(s_t) is exogenous (Bocola's s-shock), so neither
purchases nor the cap can change the expected-loss term, which is most of the spread. The concrete
implementation is to make the default probability depend on the private float, as in Costain–Nuño–Thomas:

$$
\pi_t=\frac{1}{1+\exp\!\big(-[\,s_t+\phi\,\ln((B_{t+1}-M_{t+1})/B_{t+1})\,]\big)},\qquad \phi\ge0 ,
$$

or on the market-value debt burden, as in Lorenzoni–Werning:

$$
\pi_t=p^d\big(s_t+\phi_b\,(\delta B_{t+1}/q_t-\overline{\delta B/q})\big).
$$

φ = 0 nests the current model. φ would be calibrated to the default-risk share of the OMT effect (KNVJ: 37%), or
bounded by Bocola–Dovis.

### M4. Fiscal backstop: loss-sharing and transfers

**Papers.** Roch & Uhlig (2018, *JIE*), "The dynamics of sovereign debt crises and bailouts" (a bailout guarantee, explicitly motivated by OMT, removes the self-fulfilling run); Corsetti & Dedola
(2016) on central-bank losses and inflation.

**Mechanism.** A third party commits resources in default states, which raises recovery for private holders or
removes the incentive to default.

**In this model.** It is present but mis-targeted. The Eurosystem's book takes the haircut in default and the
capital key sends 92.9% of the loss to F, but private holders' recovery is unchanged, the Eurosystem buys nothing
in default, and π does not respond. So the transfer insures the D banks *on the bonds the Eurosystem has already
bought*: it relieves bank net worth, not the market price. Lowering the spread through this channel would require
raising recovery for the bonds still in private hands. Under national funding the capital key also works
*against* the backstop in no-default states: the book's excess return is a transfer from D to F that D cannot
smooth (§2.5, §5.2 (iii)).

### M5. Redenomination risk

**Paper.** Krishnamurthy, Nagel & Vissing-Jorgensen (2018), who identify redenomination from the gap between
local-law euro bonds and foreign-law or dollar bonds: 13% of the effect on average, and none for Italy.

**In this model.** Absent: there are no currencies and no exit option. It belongs to a nominal extension.

### What the evidence says about the size of a purchase's price impact

| Source | Measured impact |
|---|---|
| Eser & Schwaab (2016, *JFE*), SMP | about **−3 bp** on 5-year yields per purchase of 0.1% of the outstanding stock, i.e. about −30 bp per 1% of the stock, plus large announcement effects |
| Ghysels, Idier, Manganelli & Vergote (2017, *JEEA*) | SMP lowered program-country yields at high frequency (daily regressions understate it, because purchases responded to price falls) |
| Krishnamurthy, Nagel & Vissing-Jorgensen (2018) | SMP + OMT: 2-year yields −400 bp (Italy, Spain), −500 bp (Ireland, Portugal), −1000 bp (Greece); 37% default risk / 13% redenomination / 50% segmentation |
| Costain, Nuño & Thomas (2025) | PEPP announcement about −80 bp on Italian yields, +15 bp from flexibility; mostly credit risk premium |

These set the target for the model. A credible TPI should deliver, at the stress node, a price impact of the
order of **tens of bp per 1% of the stock**, through segmentation and default-risk extraction.

---

## 5. The most concrete way to make the TPI lower spreads

### 5.1 The design rule

By §3.1 a backstop can lower the spread only by acting on a term of (P) that is both **large** and
**responsive to purchases**:

| Term of (P) | Size at the shock (national) | Responds to purchases in the shipped model? | What makes it respond | Literature |
|---|---:|---|---|---|
| expected loss | **452 bp/yr** | no: π exogenous | default risk that falls with the private float | M3: Costain–Nuño–Thomas; Cole–Kehoe; Lorenzoni–Werning; Corsetti–Dedola |
| liquidity premium | 155 bp/yr | no: Z divertible at λ | Eurosystem liabilities outside bank capacity | M1: Gertler–Karadi; KNVJ "segmentation" |
| risk premium | −16 bp/yr | yes, but it is about zero or negative | a funding closure in which default is a bad state for the banker | M2: Vayanos–Vila; Costain–Nuño–Thomas |

The largest term, and the only one that can deliver an announcement effect without purchases, is the expected
loss. That requires default risk to depend on the backstop. The second is capacity, which is how purchases
compress premia in every intermediary model with credit policy. Under **national** funding there is a third
requirement, which the runs uncovered: the backstop's profit and loss must not drain the issuer's resources,
because D cannot borrow abroad to smooth them.

### 5.2 The changes

**(i) Non-divertible Eurosystem claims, λ_Z = 0.** One line in the incentive constraint (`point_map.py:311`):

$$
\lambda\mathcal{A}_t=\lambda_K Q_tK_{t+1}+\lambda_b\,q_tb^D_{D,t+1}+\lambda_b^F p_tq^F_tb^F_{D,t+1}+\lambda_K L^{wc}_t+\underbrace{\lambda_Z}_{=0}\,q_tM_{t+1}.
$$

A purchase then relaxes μ within the period: ∂μ/∂M = −(1−μ)·λ_b·q/(λ𝒜) < 0. The justification is
institutional: reserves sit at the central bank and cannot be absconded with, which is Gertler & Kiyotaki's
(2010) device for claims that are easier to enforce. It reopens no portfolio-substitution margin, because the
bank does not choose Z. Alone it halves the wrong-signed effect but does not reverse it (§2.2).

**(ii) Default risk extraction: π depends on the private float.** In `_continuation`:

$$
\pi_t=\Lambda\Big(s_t+\phi\,\ln\frac{B_{t+1}-M_{t+1}}{B_{t+1}}\Big),\qquad \Lambda(z)=\frac{1}{1+e^{-z}},\qquad \phi\ge0 .
$$

Here π_t is the probability of default at t+1, B_{t+1} − M_{t+1} the end-of-period stock left in private hands,
and Λ the logistic function. It is the model analogue of Costain–Nuño–Thomas's extension, in which the default
probability falls as the central bank holds more of the debt, which they rationalise by rollover crises.

- **Nesting:** φ = 0 is the shipped model. At the rest point M = 0, so the rest point and the deterministic SS
  are untouched for any φ.
- **Why it produces an announcement effect:** markets anticipate purchases at the floor in stress states, so
  the expected float there is smaller, π along the expected path is lower, and the continuation price rises
  today.
- **Calibration:** φ is the semi-elasticity of the default odds with respect to the private float. Three
  targets are available: the default-risk share of the SMP/OMT yield decline (KNVJ: 37%); the stock effect per
  1% of the stock (Eser & Schwaab: about −30 bp on five-year yields); or an upper bound from Bocola & Dovis's
  finding that non-fundamental risk was limited in Italy. φ = 1 (used here) means that buying 43% of the stock,
  the full-grid union purchase at the headline shock, cuts the default odds by 43%.
- **Alone, under the capital key,** it is destructive under national funding (§2.5).

**(iii) Remit the Eurosystem's P&L to the issuer, κ_D = 1.** In `_government`, the profit
Π = Ξ·M_lag − O_lag (the book's payoff minus its obligation) goes to D's treasury, which lowers D's issuance,
instead of 92.9% to F's taxpayers. Why it is needed under national funding:
- to hold a low cap at the stress nodes the Eurosystem must hold most of the stock (79% at the headline shock
  with φ = 1);
- in every no-default quarter that book earns the bond's excess return, about 5.8%/yr at the shock;
- under the capital key 92.9% of it leaves D;
- D cannot borrow abroad, so D's saving falls, its deposit market clears at a deposit rate above 1,100 bp/yr, and
  the bond price, investment and bank net worth collapse with it.

Remitting to D removes the drain. With (i) + (ii) + (iii) the national cap walk reaches the 200 bp target, and
the announcement alone halves the spread at the headline shock (§2.5).

The cost is that the backstop no longer shares losses across countries: in default D's treasury bears the
haircut on the Eurosystem's book. The Eurosystem made exactly this choice for most PSPP holdings, which national
central banks hold without loss-sharing. It is a design choice the paper should state, and it needs the default
branch to evaluate (the IRF conditions on no default).

**(iv) The structural alternative: a funding closure that makes default a bad state (untested here).** The
Eurosystem refinancing closure (D3 in `deposit_portfolio_review.md`) pins D's marginal funding rate to the
policy rate outside collateral-binding states. That should do three things:
- remove the captive-depositor insurance, making the default state again a high-Ω state with a positive risk
  premium;
- stop purchases from reaching the household discount factor (§2.4);
- give D the external margin that the capital-key remittance needs.

`uni_phi1_lz0` is the test with union funding as that margin: it holds the 200 bp cap under the capital key,
with the largest price impact of all runs (−43 bp per 30% at the shock).

### 5.3 Test results: the fix ladder (national deposits, quick grid)

| Economy | Lowest cap | Announcement: sovereign / lending spread | 30% purchase at the shock | Headline shock with TPI: spread on impact / Σ output 25 q |
|---|---:|---|---:|---|
| shipped TPI | 493.75 bp | +0.1 / 0 bp | **+10.2 bp** | 339.5 bp / −5.99 (nothing bought) |
| (i) λ_Z = 0 | 493.75 bp | +0.1 / 0 bp | +4.3 bp | identical to no TPI |
| (ii) φ = 1 | 406.25 bp | +4.3 / −80 bp | −4.3 bp | at the cap, 79% bought; D deposit rate 1,215 bp in q2 |
| (i) + (ii) | 375 bp | −12.0 / −74 bp | −9.4 bp | at the cap, 79% bought; D deposit rate 1,158 bp in q2 |
| **(i) + (ii) + (iii)** | **200 bp** ✓ | **−23.2 / −47 bp** | **−29.1 bp** | **166 bp / +1.91 (nothing bought)** |
| *union with (i) + (ii)* | 200 bp ✓ | −33.3 / −95 bp | −43.0 bp | coarse-grid pathology |

The no-TPI benchmark at the headline shock is 339.5 bp on impact, with Σ output −5.99 over 25 quarters.

### 5.4 Implementation order and checks

1. **Ship (i), (ii) and (iii) behind switches** (`lambda_Z_D`, `phi_float`, and the existing `tpi_key_D`), with
   defaults that nest the shipped model. Checks:
   - bit-identical period-map outputs at the defaults (the 227-output harness);
   - N1: the rest point is unchanged at M = 0 for any φ;
   - N3: the within-period invariance of the swap now holds for everything except μ when λ_Z = 0, so it needs a
     λ_Z = 0 variant asserting dμ = −(1−μ)·λ_b·q·dm/(λ𝒜);
   - N4: λ_Z enters no budget and κ_D only reallocates Π between the two treasuries, so the union budget
     identity holds unchanged.
2. **Lift the national-TPI `assert`.** Apply the 0.1% M band only when the TPI is off; the TPI needs the
   full box.
3. **Make the cap target coherent with fundamentals.** Cap at the default-warranted price q⋆ (the draft's
   "disorderly market" design, eq. q-warranted) or at the rest-point spread (`TPI_CAP = "rest"`). A fixed
   200 bp cap below the expected-loss spread is a transfer, not a backstop (Proposition 3).
4. **Full grid** (`S_REFINE = 5`) for every reported number. The coarse grid overstates the union announcement
   several-fold and produces the pathological union TPI impulse responses.
5. **Calibrate φ** to one of the three targets and report the TPI under both deposit markets. The paper's claim
   then becomes:

   > With national bank funding, an asset-swap backstop protects transmission only through bank capacity and
   > the non-fundamental part of default risk, and only if its profits stay in the issuing country. The
   > risk-transfer channel is closed, because captive depositors already provide the insurance the backstop
   > sells.

6. **Then (iv),** which the companion note recommends for its own reasons. Re-run the ladder on it.

---

## 6. Caveats and what to run next

- **Quick grid.** Every run here uses the 25-point coarse grid. Signs and cross-economy comparisons are
  robust to that, since every economy is solved the same way. Levels are not: the coarse union announcement
  (−49 bp sovereign, −95 bp lending) is several times the full-grid one (−4.3 / −36 bp). The next step is the
  full grid (`S_REFINE = 5`) for `nat`, `nat_phi1_lz0_key1` and `uni_phi1_lz0`.
- **The box.** In the φ economies the Eurosystem's impact purchase (79% of the stock) sits at the edge of
  the TPI book's box (±80%), and every later quarter of those impulse responses is clipped back into the box
  (the run logs print the escapes). Only their impact quarter is clean. Widen `m_band` above 1 before
  reading the dynamics.
- **φ is illustrative.** φ = 1 was chosen to show the channel, not calibrated. §5.4 lists three calibration
  targets.
- **κ_D = 1 is a policy choice, not a free lunch.** It removes the no-default drain to F, but also the
  cross-border loss-sharing in default. The impulse responses condition on no default, so they show only the
  benefit. A welfare comparison needs the default branch.
- **The lending-spread jump in the full-fix economy** (+267 bp on impact on exact reads) comes with a smaller
  fall in the bond price and a rise in bank net worth and investment. The franchise-value reading in §2.5 is an
  interpretation, not a decomposition. Run the output decomposition on that economy before writing it up.
- **Pre-existing open items** that these runs inherit, unrelated to the TPI:
  - impact D consumption rises under national deposits (+0.19% here), the C_D gate;
  - the coarse union TPI impulse response is pathological (impact output −1.8%), as already recorded on
    2026-10-02 to 05.
- **Nothing is committed.** The worktree patches (§1.2) live in the session scratchpad. Porting them means
  adding `lambda_Z_D` and `phi_float` to `calibration/global_projection.py` with their nesting defaults, and
  N-tests as in §5.4.

## 7. References

- Bacchetta, P., E. Perazzi and E. van Wincoop (2018). Self-fulfilling debt crises: What can monetary policy do? *Journal of International Economics* 110.
- Bocola, L. (2016). The pass-through of sovereign risk. *Journal of Political Economy* 124(4).
- Bocola, L. and A. Dovis (2019). Self-fulfilling debt crises: A quantitative analysis. *American Economic Review* 109(12).
- Calvo, G. (1988). Servicing the public debt: The role of expectations. *American Economic Review* 78(4).
- Camous, A. and R. Cooper (2019). "Whatever it takes" is all you need: Monetary policy and debt fragility. *AEJ: Macroeconomics* 11(4).
- Cole, H. and T. Kehoe (2000). Self-fulfilling debt crises. *Review of Economic Studies* 67(1).
- Corsetti, G. and L. Dedola (2016). The mystery of the printing press: Monetary policy and self-fulfilling debt crises. *Journal of the European Economic Association* 14(6).
- Costain, J., G. Nuño and C. Thomas (2025). The term structure of interest rates in a heterogeneous monetary union. *Journal of Finance* 80(4).
- Dedola, L., P. Karadi and G. Lombardo (2013). Global implications of national unconventional policies. *Journal of Monetary Economics* 60(1).
- Eggertsson, G. and M. Woodford (2003). The zero bound on interest rates and optimal monetary policy. *Brookings Papers on Economic Activity* 2003(1).
- Eser, F. and B. Schwaab (2016). Evaluating the impact of unconventional monetary policy measures: Empirical evidence from the ECB's Securities Markets Programme. *Journal of Financial Economics*.
- Gabaix, X. and M. Maggiori (2015). International liquidity and exchange rate dynamics. *Quarterly Journal of Economics* 130(3).
- Gertler, M. and P. Karadi (2011). A model of unconventional monetary policy. *Journal of Monetary Economics* 58(1).
- Gertler, M. and P. Karadi (2013). QE 1 vs. 2 vs. 3...: A framework for analyzing large-scale asset purchases as a monetary policy tool. *International Journal of Central Banking*.
- Gertler, M. and N. Kiyotaki (2010). Financial intermediation and credit policy in business cycle analysis. *Handbook of Monetary Economics* 3A.
- Ghysels, E., J. Idier, S. Manganelli and O. Vergote (2017). A high-frequency assessment of the ECB Securities Markets Programme. *Journal of the European Economic Association*.
- Krishnamurthy, A., S. Nagel and A. Vissing-Jorgensen (2018). ECB policies involving government bond purchases: Impact and channels. *Review of Finance* 22(1).
- Lorenzoni, G. and I. Werning (2019). Slow moving debt crises. *American Economic Review* 109(9).
- Roch, F. and H. Uhlig (2018). The dynamics of sovereign debt crises and bailouts. *Journal of International Economics* 114.
- Vayanos, D. and J.-L. Vila (2021). A preferred-habitat model of the term structure of interest rates. *Econometrica* 89(1).
- Wallace, N. (1981). A Modigliani–Miller theorem for open-market operations. *American Economic Review* 71(3).
- Aguiar, M., M. Amador, E. Farhi and G. Gopinath (2015). Coordination and crisis in monetary unions. *Quarterly Journal of Economics* 130(4).

**Web sources consulted** (publication details and headline findings):
[Costain–Nuño–Thomas, Banco de España summary](https://www.bde.es/wbe/en/areas-actuacion/analisis-e-investigacion/investigacion-destacada/research-features/term-structure-interest-rates-heterogeneous-monetary-union.html) ·
[Camous–Cooper, AEJ: Macro](https://www.aeaweb.org/doi/10.1257/mac.20170167) ·
[Krishnamurthy–Nagel–Vissing-Jorgensen, NBER w23985](https://www.nber.org/papers/w23985) ·
[NBER Digest summary of KNVJ](https://www.nber.org/digest/feb18/quantifying-impact-ecb-policies-during-debt-crisis) ·
[Bocola–Dovis, AER](https://www.aeaweb.org/doi/10.1257/aer.20161471) ·
[Bacchetta–Perazzi–van Wincoop, NBER w21158](https://www.nber.org/papers/w21158) ·
[Eser–Schwaab, EconStor](https://econstor.eu/handle/10419/154075) ·
[Corsetti–Dedola, IDEAS](https://ideas.repec.org/a/bla/jeurec/v14y2016i6p1329-1371.html) ·
[Lorenzoni–Werning, AER](https://www.aeaweb.org/doi/10.1257/aer.20141766) ·
[Aguiar–Amador–Farhi–Gopinath, NBER w20277](https://www.nber.org/papers/20277) ·
[Roch–Uhlig, IDEAS](https://ideas.repec.org:443/p/imf/imfwpa/2016-136.html) ·
[Bocola (2016) and the LTRO, Minneapolis Fed](https://minneapolisfed.org/article/2015/risky-business) ·
[Gertler–Karadi (2011), IDEAS](https://ideas.repec.org/a/eee/moneco/v58y2011i1p35-38.html) ·
[DNB WP 778, fragmentation](https://www.dnb.nl/en/publications/research-publications/working-paper-2023/778-identifying-financial-fragmentation-do-sovereign-spreads-in-the-emu-reflect-differences-in-fundamentals/)
