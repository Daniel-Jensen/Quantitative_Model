# A TPI backstop in the global model: design, derivations, obstacles

**Date:** 2026-10-02 · **Branch:** `OMT-fix` · **Replaces:** the stochastic LTRO backstop
(`ltro_backstop_plan.md`, deleted from the code 2026-10-01, kept here as history).
**Pipeline:** global projection only (`global_projection/`). The sequence-space model's TPI
(`linear_ssj/solve/tpi.py`) is a different object and is not touched.

**Status:** implemented and gated (steps 1–4 of the plan); the full runs are in §9.

---

## 1. The instrument

In the no-default regime the Eurosystem buys the D sovereign whenever its spread over the F
bond would exceed a cap, and pays for the bonds with a **safe claim on itself** held by the D
banks. Four design choices, each the user's:

| choice | decision |
|---|---|
| the floor | a cap `tpi_cap_bp` (annualised bp) on the D–F spread in HM flow yields |
| eligibility | the no-default regime only; nothing is bought in the default regime |
| holdings | held to maturity: purchases `m ≥ 0`, never sales |
| P&L | remitted every period and shared by the euro-area capital key, `tpi_key_D = 0.071` to D |

The HM flow yield of a perpetuity priced `Q` with coupon decay `δ` is `y(Q) = δ(1−Q)/Q`, so the
cap is a floor on the D price:

```
y_D − y_F ≤ c̄    ⟺    Q_bD ≥ Q̲ = δ_D / (δ_D + y_F(Q_bF) + c̄),      c̄ = tpi_cap_bp / 4e4
```

## 2. New objects (deliverable 1)

| object | where | meaning |
|---|---|---|
| `x_cb` | unknown and stored rule (14th) | the TPI variable: purchases `m = B·max(x,0)`, or the price's slack over the floor `max(−x,0)` |
| `M_cb` | state (11th) | Eurosystem holdings carried in, issuer per-capita units |
| `O_cb` | state (12th) | the Eurosystem's gross obligation to the D banks, `(1+R)·Z` |
| `M = (1−δ_D)·surv·M_lag + m` | derived | the book, held to maturity (haircut in default) |
| `Z = Q_bD·M` | derived | the safe claim, held by the D banks, paying `1 + rdep_D` |
| `Π = Ξ_D·M_lag − O_lag` | derived | the Eurosystem's P&L, `Ξ_D = surv·(δ_D + (1−δ_D)·Q_bD)` |
| `Q̲` | derived | the floor |

`O` is a state because it carries last period's price, which nothing else remembers.

## 3. Equations that change (deliverables 2–6), all in `point_map.py`

1. **Gross debt and clearing.** `B_D = b_DD_lag + b_DF_lag + M_lag`; `b_DD = B' − b_DF − M`.
   Coupons and the Bohn tax apply to the gross stock (the Eurosystem is pari passu).
2. **D-bank book.** `assets_D += Z`, payoff `X_D += O_lag`; net worth, deposits, `P'`,
   the entrant transfer and dividends keep their formulas.
3. **Divertable base.** `lev_D += λ_bD·Z`. Under the single λ the swap leaves it unchanged.
4. **Treasuries.** D: `new_D = (G_D + coupon_D − Tax_D − κ_D·Π)/Q_bD`. F (whose `B_F` is fixed):
   `Tax_F = Tax_F,ss − (1−κ_D)·Π/(sz·p)`.
5. **The purchase rule — Garcia–Zangwill.** On, in `d = 0`, one smooth variable carries both
   sides of the complementarity: `m = B_ss·max(x, 0)`, `Q_bD = Q̲ + Q_bD,ss·max(−x, 0)`, so
   `m ≥ 0`, `Q_bD ≥ Q̲` and `m·(Q_bD − Q̲) = 0` hold by construction. Residual 14 ties the stored
   `Q_bD` rule to that price; in `d = 1` (and with the TPI off) it holds `x = 0`, `m = B·x`.
6. **The D-bank bond FOC in KT form** (residual 8, when on): `FB(−foc_D, b_DD/B_ss) = 0`, with
   `foc_D = (E[Ω·payD] − Q_bD·(E[Ω]R + λ_bD·μ))/E[Ω]`. The banker is linear, so a price above his
   valuation clears at the corner `b_DD = 0` with the Eurosystem holding the rest. In the
   no-default regime this residual pins `x`; in default it pins `Q_bD`.
7. **Laws of motion.** `M' = M`, `O' = (1 + rdep_D)·Z`.

`FB(a, b) = a + b − √(a² + b² + ε²)` (`tpi_eps = 1e-5`): zero iff `a, b ≥ 0` and `ab = ε²/2`.

**Three formulations were tried; two failed, and why they failed fixed the third.**
- A **smooth min** for the purchase pair has no slope in `m` while the price sits below the
  floor, so the Newton could only move the price: stalled at max|F| = 1.7e-2.
- **Fischer–Burmeister on a stored `m`** converges at every cap down to 200 bp, but `m` is a
  kinked function of the state and its Chebyshev fit is not: on the coarse grid the quadratic
  in `s` through (0, 0, 0.29) **bought 19.7% of the stock at the headline shock, where the
  floor does not bind**, and a degree-4 fit turns the same kink into sales between nodes.
- **Garcia–Zangwill with next period's price computed as `Q̲' + max(−x', 0)`** puts the kink
  inside every expectation (at the quadrature points), and the finite-difference Newton
  stalled at the 500 bp rung (max|F| 5e-3, steps cut to 1/64).
- **Garcia–Zangwill with the continuation reading the stored `Q_bD` rule** removes that kink,
  but the Newton still stalled walking the cap down (500 bp rung, steps cut to 1/64) where the
  Fischer–Burmeister form takes full steps to 200 bp.
- **The shipped version roots in one form and reads in the other.** The two forms describe
  the same allocation at every node (measured: the GZ residual at the converged FB solution is
  2.2e-9 and the GZ Newton accepts it in one step), so the conversion is exact there:
  `x_GZ = x_FB` where the floor binds, `−gap` where it does not (`recursive_experiment._tpi_form`).
  The cap homotopy and the refined-grid solve run in the FB form (`cal["tpi_gz"] = False`, the
  14th slot the purchase share, residual 14 `FB(x, (Q_bD − Q̲)/Q_bD,ss)`); each hands back in
  the GZ form with a one-to-two-step polish (`_tpi_polish`), and every reader sees GZ rules:
  kink-exact, never buying where the floor is slack. The continuation reads next period's price
  off the stored `Q_bD` rule in both forms.

**The equilibrium** at each of 12 states and in each of two regimes: 14 unknowns, 14
residuals, the 6 identity rules for the derived objects, 20 stored rules per regime.

## 4. The risk transfer (deliverable 7)

A marginal purchase `dm` at price `Q`: the bank's bonds fall by `dm`, its claim `Z` rises by
`Q·dm`, so assets, the divertable base, deposits, `P'` and `μ` are unchanged **today**.
**Tomorrow**, in regime `d'`, the bank is owed `R·Q·dm` instead of the bonds' payoff `Ξ'(d')·dm`:

```
dn'(d') = [(1−f) + ω_ent·κ_D] · (R·Q − Ξ'(d')) · dm          dΠ'(d') = −(R·Q − Ξ'(d')) · dm
```

(the `ω_ent·κ_D` term: D's share of the P&L is new D debt the D bank absorbs, and the entrant
transfer is a share of assets). Without default `Ξ' ≈ R·Q` plus the premium, so the bank gives
up a little; in default `Ξ' = rec·(δ + (1−δ)Q')` and it is insured. `test_recursive_nesting` N3
asserts both to 1e-12 (measured: −2.1e-5 and +4.6e-3 per `dm = 0.01`).

The chain is then `n'_d ↑ → μ'_d ↓ → α'_d ↓ → Ω_d ↓`: `cov(Ω, payD)` shrinks, the risk premium
falls and `Q` rises. Against it runs the franchise channel: a safer future lowers `α'`, hence
`E[Ω]`, which raises today's `μ`. The sign is measured, not assumed (§9).

## 5. Why there is no ψ_D^H (deliverable 8)

The banker is linear, so his bond FOC is a price condition with no quantity in it. A purchase
moves `Q` only through equilibrium objects: `μ` (through the balance sheet, which the swap
neutralises) and `Ω'` (through the continuation states `b_DD'`, `M'`, `O'`). So
`Q = Q(B − M; n, μ, Ω, p^d)` is a general-equilibrium inverse demand and `∂Q/∂M` is
endogenous — measured in §9, small at rest and larger in stress. A home-bias cost `ψ_D^H`
would put a quantity into the FOC by assumption; the TPI's price effect would then be that
assumption.

## 6. What does not change (deliverable 9)

The Ω kernel, the α recursion, the IC and Bocola's closed-form μ; both capital Eulers, labour
FOCs and deposit Eulers; deposit-UIP and union clearing; goods_D; the D-bond FOC itself (only
its KT form); the F bank's D-bond FOC with `ψ_bD_F` and both F-bond FOCs; working capital,
loans, dividends; firms, capital and trade; the `s` process, `p^d(s)` and the haircut; the
Bohn rule's formula; **the deterministic steady state** (`M = O = Z = Π = 0`).

## 7. The state box: a shear, not the natural axes

The book's two states are zero at the SS and their bands are symmetric round 0, so the SS is a
node and **every node with `M = O = 0` is a node of the old 10-state grid**: on that slice the
12-state interpolant is the 10-state one, and with the TPI off the model nests the 10-state
solution to solver tolerance (measured: rules 5.5e-10, risk IRF 4.5e-9).

The box is drawn in `z_bDD = b_DD + M` and `z_O = O − ρ·M` (`ρ = (1 + rdep_ss)·Q_bD,ss`), not in
`b_DD` and `O`. On the natural axes the nodes are states no purchase reaches (bonds sold with
nothing owed; a claim owed with nothing held), a purchase walks `b_DD` out of its own band, and
the measured price impact of a forced purchase came out **non-monotone** (−0.5% at 5% of the
stock, +0.3% at 30%) with `μ_D` jumping 4×. On the shear it is monotone and small (§9). The
corner `b_DD = 0` sits at the edge of the `M` band (`m_band = 1`). `SmolyakGrid`'s existing
rotation does the work, with centre 0, so every other coordinate and the slice are untouched.

## 8. Obstacles (deliverable 10) and how each was met

1. **The price impact per euro is small, so the cap decides the size.** At fixed rules a
   purchase of 30% of the stock lifts `Q_bD` 0.24% at rest and 0.61% at the headline shock
   (115-point grid; 0.46% for 31% on the coarse grid). The solved TPI economy prices the bond
   3.8% above the no-TPI one at the shock, most of it through the continuation (§9.2). A cap
   far below the market spread is met only at the corner, the Eurosystem holding the D bank's
   whole book: at 21 of the 23 nodes at p^d = 4.8%/qtr. At the headline shock (302 bp against
   the 200 bp cap) the purchase is interior.
   Below the expected-loss spread the instrument is a cross-border transfer, not a backstop.
2. **Two kinks** (the purchase switch and the `b_DD` corner). Met with Garcia–Zangwill for
   the purchase pair, Fischer–Burmeister for the corner, a **smoothing continuation** on the
   refined grid (ε = 1e-2 → 1e-3 → `tpi_eps` = 1e-4; at the final 1e-5 first chosen, the
   finite-difference step of 1e-6 saw a near-kink and the 5-node solve stalled), and a **cap homotopy**
   (`recursive_experiment._solve_tpi`: 700 bp, where the cap binds nowhere on the box, down to
   the target in 50 bp rungs in the FB form, a failed rung retried at half the step, then
   handed back in the GZ form). §3 records the formulations that failed and why.
3. **The franchise channel** may offset the risk-premium channel: measured in §9.
4. **Cost.** 12 states: 25 coarse / 115 refined points, 20 rules × 2 regimes; a quick run went
   528 s → 824 s, and each coarse cap rung takes 4–7 Newton steps (~65 s each).
5. **Box.** The shear in §7.
6. **Walras.** The global model had no goods_F check. Adding one found a **pre-existing** leak
   (STATE Part II, PROGRESS 2026-10-02): the rep-agent anchor `hh_T` (0.583) is a constant
   standing in for the working-capital flow. The identity behind it holds at any point, so
   `test_recursive_nesting` N4 checks that the TPI's flows cancel exactly: remainder 2e-10 with
   the TPI on, 3.5e-10 off, at random states with a live book.
7. **Cross-border loss sharing.** In default F households pay `(1−κ_D)·|Π|` through taxes; it
   closes through the existing budgets (N4).
8. **Flight to safety.** `Q̲` moves with `y_F`, so when the F yield falls in stress the floor
   rises and the Eurosystem buys more exactly when buying is costly.

## 9. Results

### 9.1 Preview: the coarse grid (`run.py global --quick`, NOT converged), cap 200 bp

Grid μ=1, 25 points × 2 regimes, three nodes in `s`. **The full run (§9.2) supersedes every
number here**; they are recorded because they shaped the design (exact reads, §8).

**Solve.** The cap homotopy 700 → 200 bp roots every rung in 4–6 Newton steps (FB form), then
hands back in GZ form. On this grid the floor binds only at the riskiest node (p^d 4.8%/qtr):
interior down to 400 bp — purchases 3 / 15 / 29 / 46 / 67% of the SS stock at
600 / 550 / 500 / 450 / 400 bp — and at the corner (`b_DD = 0`, the Eurosystem holding 81% of
the stock, `foc_D < 0`) from 350 bp down. KT conditions on the solved rules hold to 1e-9.

**The announcement (rest point, nothing bought).**

| object | no TPI | TPI | diff |
|---|---|---|---|
| sovereign spread, bp/yr | 61.8 | 34.7 | −27.0 |
| credit spread, bp/yr | 79.0 | 71.0 | −8.0 |
| `mu_D` (IC multiplier) | 0.0098 | 0.0089 | −0.0010 |
| `alpha_D` (franchise value) | 1.1769 | 1.1801 | +0.0031 |
| `Y_D` (D output) | 0.99996 | 1.00060 | +0.064% |

The franchise channel does not dominate here: `alpha_D` rises with the TPI.

**The headline risk shock (p^d → 1.98%/qtr), both paths cleared at every quarter.** The TPI
never fires: the impact spread is 164 bp, under the cap, against 419 bp without it — the whole
compression is the announcement. The credit spread stays at zero (151 bp without) and bank
net worth falls 3.9% (5.2%). **Output falls more on impact, −0.31% against −0.12%**: the exact
output decomposition gives the credit-spread leg +0.13 pp and the relative-price leg −0.40 pp
(p +3.0% against +0.25%), because D demand falls harder (consumption −0.57%, investment −1.28%)
as the D deposit rate rises 19 bp instead of falling 44 bp. Whether this survives the refined
grid is the open question §9.2 answers.

**Fitted reads are not usable for the TPI path.** Read off the fitted rules the same shock
bought 53% of the stock on impact and then asked for more bonds than the D bank held
(`b_DD` −120%), running into the box wall for 24 of 25 quarters; three nodes in `s` cannot
place the floor's boundary, even for the smooth `x`. `tpi_experiment` therefore clears both
economies' IRFs at every quarter.

### 9.2 The full resolution (`S_REFINE = 5`, 115 points), cap 200 bp — 2026-10-05

**Solve.** The no-TPI model on the 115-point grid reproduces the published numbers (impact
output −0.1105% read off the fitted rules and −0.1278% cleared exactly, bond price −9.2633%). The
cap ladder roots every rung 700 → 200 bp; the TPI's refined solve walks 3 → 5 nodes in `s` with the
smoothing continuation ε = 1e-2 → 1e-3 → 1e-4 (§8) and ends at max|F| 2.2e-9. Without the
continuation the 5-node solve stalled: the refined grid puts nodes right at the floor's boundary —
the floor binds at **43 of the 115 no-default nodes**, the whole p^d = 4.8%/qtr level and 20 of the
23 nodes at p^d = 1.58%/qtr — and with the final smoothing from the start each of those is a
near-kink. 21 nodes, all at p^d = 4.8%/qtr, sit at the corner (the D bank holds none of the bond).
On the solved rules: nothing is bought in default (x ≡ 0); the spread never exceeds the cap (max
199.999995 bp; at nodes on the boundary it sits up to 1 bp below it, the ε = 1e-4 smoothing); the
complementarity products sit at the smoothing's bias, |m·gap| = ε²/4 = 2.5e-9 at every node (the
GZ read) and |b_DD·foc| ≤ 6.6e-9 (the FB pair's ε²/2 = 5e-9 plus the solve's tolerance).

**The announcement (rest point, nothing bought).**

| object | no TPI | TPI | diff |
|---|---|---|---|
| sovereign spread, bp/yr | 55.8 | 51.6 | −4.3 |
| lending spread, bp/yr | 88.7 | 52.6 | −36.1 |
| `mu_D` (IC multiplier) | 0.0110 | 0.0066 | −0.0045 |
| `n_D` (D-bank net worth) | 1.9395 | 1.9480 | +0.44% |
| `Y_D` (D output) | 0.99941 | 1.00116 | +0.18% |
| `I_D` (D investment) | 0.20460 | 0.20574 | +0.56% |

**The headline shock (p^d → 1.98%/qtr, both paths cleared exactly every quarter).** The TPI
fires once: the Eurosystem buys **43.3% of the SS stock on impact**, the spread is held at the cap
(200 bp, against 302 without), and from the next quarter the spread sits below the cap, nothing
more is bought and the book runs off at the bonds' maturity rate (11.5% of the stock still held
entering q24). The purchase is interior: the D bank keeps 39% of the stock and its FOC holds (the
corner is reached only at the p^d = 4.8%/qtr nodes).

| | no TPI | TPI |
|---|---|---|
| output, impact | −0.128% | **−0.080%** |
| output, trough | −0.128% (q0); above 0 from q3 | **−0.147% (q6)**; −0.071% at q24 |
| consumption, impact / trough | −0.02% / −0.04% | −0.24% / −0.34% |
| bank net worth | −4.17% on impact, above its no-shock path from q3 | −1.47% on impact, −1.23% at q4, −0.63% at q8 |
| lending spread, impact | 180 bp | 34 bp, rising to 83 bp at q6 |
| D bond price vs no-shock path, impact | −8.91% | −5.68% |

Impact output, exact decomposition (pp): credit spread −0.218 → −0.024, deposit rate +0.116 →
+0.001 (without the TPI the deposit rate falls 54 bp and cushions the wage bill; with it, 0.6 bp),
relative price −0.026 → −0.057.

**Why the bond is worth more on impact** (+3.82%, the FOC legs, residual 0.0000%): the
continuation price +3.25 (the floor under future resale values), the liquidity premium +0.37 (the
bank constraint relaxed), the risk premium +0.30 (the banks de-risked), the expected loss +0.05,
discounting −0.16. A one-off purchase at fixed rules moves the price +0.24% (rest) / +0.61%
(shock) per 30% of the stock.

**Why the downturn is longer with the TPI.** Along this path no default occurs, so the risk
transfer of §4 runs in its no-default direction: the D banks sold 43% of the stock to the
Eurosystem at the floor and hold the safe claim instead, so the bonds' excess return over the
deposit rate — the default premium and the price recovery, which in the no-TPI economy rebuild
bank net worth within three quarters — accrues to the Eurosystem (its P&L is +0.34% of quarterly
GDP at q1, declining to 0.05% at q24), and the capital key sends 92.9% of it to F households. D
bank net worth stays down, the lending spread rises after impact, and D consumption falls. **The
IRF conditions on no default, so it shows what the insurance costs and none of what it pays**:
the welfare comparison needs the default branch, where the Eurosystem — and, by the capital key,
mostly F — takes the haircut instead of the D banks.

**The coarse preview (§9.1) is superseded.** On the coarse grid the TPI never fired at the
headline shock and output fell more with it; both were coarse-grid artefacts, and its +9.2% price
effect fell to +3.8%.
