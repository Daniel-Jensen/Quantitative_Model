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
