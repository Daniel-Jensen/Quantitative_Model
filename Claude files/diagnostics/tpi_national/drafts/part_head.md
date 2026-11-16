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

The drivers are `tpi_national.py` (solve, compare, price legs) and `tpi_mechanism.py` (formula (I), measured),
in the session scratchpad. The two worktrees are `dff` (λ_Z switch) and `dff2` (λ_Z and φ switches).

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
