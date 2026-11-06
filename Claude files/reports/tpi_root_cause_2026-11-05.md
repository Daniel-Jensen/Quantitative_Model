# Why the TPI cap walk cannot reach 200 bp: economics first, computation second (2026-11-05)

Global solver only (`global_projection/`), branch `calibration`, Bocola's risk box (`PD_BOX_TOP` 0.062), the
`python3 run.py` started 2026-11-04 19:54. Every number comes from the scripts in
`Claude files/diagnostics/tpi_root_cause/`, read-only on the saved rungs in `results/GLOBAL_tpi_walk/`.

**Verdict.** The 200 bp hard cap asks for an equilibrium that this model can only deliver as a corner: the Eurosystem
holds the D banks' whole book in every crisis state. The state box cannot represent that corner. The cap walk also
wastes hours on box nodes that a TPI economy never visits. More compute will not reach 200 bp.

## 1. What the run did

| time | rung | result |
|---|---|---|
| 20:41–21:23 | 644 → 602.6 bp | solved; one node buying (node 42, 5% of the debt) |
| 23:36 | 595.4 bp | accepted at max\|F\| 9.8e-5 (path tolerance 1e-4): **not a solution**, the standard is 1e-9 |
| 23:36 → 09:35 | 587.9, 589.2, 589.1, 589.3, 589.3 bp | each stalls at max\|F\| ≈ 2e-4 on `bondD_D` (the D bank's D-bond FOC), ~2 h each |

- **Halving the step does nothing.** The walk keeps each rung 3 bp (`margin`) from the next switching cap, 592.4 bp.
  From 595.4 bp every candidate (steps 8.96, 4.48, 2.24, 1.12, 0.56) is therefore pushed to ≤ 589.4 bp: the same ~6 bp
  jump, five times.
- **The 595.4 rung landed on a switch.** Node 69's switch was predicted at 601.4 bp from the 602.6 bp rules. It moved
  to ~595.5 bp as node 42 bought, so the rung sat on it (both arguments of the purchase complementarity ≈ 0).
- **The nodes being crossed are unphysical.** Node 69 is the `M_cb`− edge (the Eurosystem's D-bond book at −26% of
  the debt, the D banks holding 113% of it). Node 73, the next switch, is the `O_cb`− edge (the Eurosystem's
  debt to the D banks below zero). Purchases are ≥ 0 and the book starts at 0, so a TPI economy never visits either.
- **How it ends.** The last step falls below `min_step` (0.5 bp) and the walk stops. The run then polishes and
  reports the TPI at 595.4 bp, saved as the 200 bp experiment (`cap_bp` 200, `cap_solved` 595.4).

## 2. Economics: the cap lies below the fair spread at 24 of the 75 nodes

The D–F spread has three parts:
- expected default loss, set by the exogenous p^d (priced default probability per quarter);
- a risk premium, −Cov(Ω, R_b)/E[Ω] (Ω: the banker's kernel; R_b: the D bond's return);
- a liquidity premium, λμ/E[Ω] (μ: the multiplier on the D bank's incentive constraint).

Purchases reach only the last two.

**Fair (expected-loss-only) spread** of the D bond, risk-neutral, on the s process alone (`fair_spread.py`). The last
column gives the excess return of a risk-neutral holder who buys at the 200 bp price, with the floor defended next
quarter.

| p^d per quarter | 0.6% | 2.1% | 3.0% | 4.7% | 6.2% |
|---|---|---|---|---|---|
| fair spread | 68 bp | 145 bp | 186 bp | 255 bp | 310 bp |
| one-quarter excess return at the 200 bp price | — | +899 | +289 | **−548** | **−1077** bp/yr |

- **Who will hold the bond.** Above p^d ≈ 3.5% nobody private holds the D bond at the 200 bp price, not even a
  risk-neutral bank. The grid has 24 nodes there (23 at 6.2%, 1 at 4.7%).
- **The only equilibrium at those nodes is a corner:** the Eurosystem holds the D banks' whole book.
- **One node's purchases as the cap falls** (`node_demand.py`; node 42, p^d 6.2%; rules held at the 602.6 bp rung, so
  this is a one-off read, not the full equilibrium):

| cap (bp) | 600 | 550 | 500 | 450 | 400 | 200 |
|---|---|---|---|---|---|---|
| bought, % of D debt | 5.8 | 15.9 | 23.5 | 30.5 | 89.1 | 92.1 |
| D banks' book, % of D debt | 80.4 | 70.9 | 63.9 | 57.6 | 0 | 0 |

- **Little lever for purchases.** The D bank's constraint is slack at central states of this box, so there is no
  liquidity premium to remove. Purchases move the spread only by taking risk off the banks, about 3 bp per 1% of the
  debt bought. Nodes 58 and 49 behave the same way: at the corner by 400 bp and 300 bp.
- **So a 200 bp hard cap here is monetisation plus a fiscal transfer, not the ECB's TPI.** The Eurosystem holds ~90%
  of D debt whenever p^d is above ~3.5%, at a price above fair value. Its expected loss is 3.4% of the book per quarter
  at p^d 6.2%. Of that, 92.9% falls on F, through `tpi_key_D` (D's share of the Eurosystem's P&L, 0.071).
- **The ECB's own design excludes this.** The real TPI's eligibility criteria rule out spreads driven by fundamentals.

## 3. Computation: why the walk fails and why each attempt costs hours

**3.1 The box cannot hold the target equilibrium.**
- The `M_cb` band is ±27.9% of the D debt (`box_specs.json`, `M_cb` 0.32 × the D banks' SS book). The 200 bp corner
  needs 70–92% at 24 nodes.
- The box is sheared: it bounds b_DD + M (b_DD: the D banks' book). So whatever is bought beyond the edge is read back
  as the bank's own bonds (`clip_readback.py`, node 42):

| bought, % of D debt | 20 | 30 | 40 | 60 | 86 |
|---|---|---|---|---|---|
| true D-bank book / Eurosystem book | 66 / 20 | 56 / 30 | 46 / 40 | 26 / 60 | 0 / 86 |
| the state the continuation reads | 66 / 20 | 57 / 29 | 56.7 / 29.3 | 56.7 / 29.3 | 56.7 / 29.3 |

- **Consequence.** Past the edge the continuation stops seeing purchases, the bank's valuation stops moving, and only
  the corner solves. That is the 30% → 89% jump between 450 and 400 bp above.
- **The walk would fail near 450–500 bp even with perfect switch handling.**

**3.2 Half of the `M_cb` and `O_cb` dimensions is unphysical** (negative books), and the walk is stuck on exactly those
nodes now (§1).

**3.3 Walk mechanics.**
- The switching caps are predicted from the previous rung's rules, and they move with the solution.
- Rungs are accepted at 1e-4.
- The margin rule makes step halving useless next to a switch.

**3.4 Cost of a Newton iteration** (`time_residual.py`):
- **A Jacobian costs 48 core-minutes.** It needs 3,000 complex-step columns, and each column is a full residual: 150
  point maps in Python, 0.9 s. That is ~6 min on 8 cores.
- **Most of that work is wasted.** Each column changes only ~1/3 of the residuals (1,061 of 3,000). The coupling
  between points runs through a linear refit (basis × Φ⁻¹), so a structured Jacobian would need ~20× fewer point maps.
  That block is: each point's own unknowns, plus the chain rule through the continuation values at its quadrature
  points.
- **For comparison:** Bocola's readme gives 30 minutes for his whole solution (12 Newton solves, 3,112 unknowns,
  vectorised MATLAB with precomputed integrals, 16 cores).

## 4. Against Bocola's replication (`replication/bocola2016/Model`)

- **Same baseline method.**
  - Smolyak/Chebyshev collocation.
  - Newton on all residuals with a numerical Jacobian (`parsolve.m`).
  - Haircut ladder 0 → 0.55.
  - Next-period states clipped to the box.
  - Closed-form multiplier max(·, 0) and a net-worth floor.
  - Our no-TPI solve converges at 2.4e-10, and nothing here questions it.
- **Differences in the baseline method** (`residual_model.m`, `precomp_integral_def.m`):
  - He parameterises the expectation integrands, with the integrals precomputed: one basis evaluation per point, not 5
    quadrature nodes × 2 regimes.
  - He softens the haircut at 71 extreme points ("tapering").
  - His residual is vectorised across points.
- **His only policy experiment is a much easier object** (LTRO, `ltro_policies.m`):
  - a fixed quantity for 12 quarters;
  - solved backward from the no-policy rules, so each period's continuation is known;
  - no price rule, no complementarity, no state for the central bank's book.
- **Nothing in the package solves an infinite-horizon spread cap,** so Bocola offers no template for this TPI solve.

## 5. Options (author decisions)

1. **Stop the current run.** It cannot reach 200 bp, and the no-TPI results are already saved (step 5).
2. **What the TPI caps:**
   - **(a)** the premium over the fair spread, the "unwarranted" part;
   - **(b)** the 200 bp cap only below an eligibility threshold in p^d (the ECB's debt-sustainability test);
   - **(c)** a smooth purchase rule (SPEC.md's closed-loop rule), with no switches, so a standard Newton solves it;
   - **(d)** keep the 200 bp hard cap and accept Eurosystem holdings of ~90% in crises.
   - Whatever the rule, check the cap against what purchases can deliver at each node (`node_demand.py`) before
     solving.
3. **Computational fixes for any of them:**
   - an `M_cb` band of [0, the largest book the rule implies], not ±28%;
   - rungs accepted only at full tolerance, and switches re-predicted by node-local solves (seconds);
   - the structured Jacobian.

## Reproduce (each runs in under a minute, alongside a solve)

| script | what |
|---|---|
| `fair_spread.py` | §2 fair spread and the risk-neutral excess return at the 200 bp price |
| `rung_anatomy.py` | §1 every saved rung: purchases, books, box position, residuals by node |
| `switch_nodes.py` | §1 the next switching nodes and their box edges |
| `node_demand.py [nodes]` | §2 one node's purchase as the cap falls to 200 bp |
| `clip_readback.py` | §3.1 the true next state against the one the continuation reads |
| `time_residual.py` | §3.4 cost of one residual and one Jacobian column |

## 6. The illustrative TPI run (2026-11-05/06): a smaller shock, the cap solved

**The user's request:** make the shock smaller, so that spreads exceed the cap by little, keep it illustrative, and run.

**Design (run.py):**
- `PD_BOX_TOP` 0.02: the model's risk range is cut from 6.2% to the 2% headline shock. The highest grid spread falls
  from 618 to 300 bp.
- The headline shock stays at 2%. Without the TPI it opens the spread at 223 bp, 23 bp above the cap; a smaller shock
  stays below 200 bp.
- `TPI_CAP` stays at 200 bp; the user chose the hard cap over a smooth purchase rule.

**Why the box had to change too** (each tried, in order):

| Eurosystem-book band (`TPI_M_BAND`) | edge, % of D debt | result |
|---|---|---|
| 0.32 (old) | 27.9% | the cap needs 27% at the shock and up to ~90% at the riskiest nodes |
| 0.65 | 56.7% | 5 nodes need more; the direct solve stalled at 6.1e-4 |
| 1.15 | 100.4% | the D bank's book is negative at the M+ edge, a no-TPI corner that corrupted the rules along M |
| 0.95 | 82.9% | the riskiest nodes sell their whole book (~88%), just past the edge |
| **1.0** | **87.3%** | the whole book fits in the soft edge; the default-regime corner at M+ is −0.5 to −0.8%, mild |

- Purchases move the spread only ~1 bp per 1.2% of the debt, because the D IC is slack (μ_D = 0) on this box.
- At the M− edge the Eurosystem is short 87% and the D bank holds 175%. A TPI path never reaches it, and the cap had no
  root there, so the TPI does not buy where the Eurosystem's book is negative.

**The solver (`recursive_experiment._tpi_direct`, `TPI_WALK` direct):**
- **Cap ladder.** The cap goes 270 → 245 → 225 → 210 → 200 bp. Each rung starts from its nodes above the cap solved
  alone, with the complementarities smoothed at 1e-2. Each rung converged in 8–10 iterations (56 min in all).
  - One jump straight to 200 bp stalled at 1.2e-2.
- **Smoothing walk at 200 bp.** Factor-2 steps, halved on a failure: 1e-2 → 7.1e-3 → 4.2e-3 → 1.9e-3 → 6e-4 → 1e-4,
  all converged (5e-3 failed first and was halved).
  - A walk that stalls keeps the smallest smoothing solved (checkpoints carry it). This run did not need that.
- **Then** the GZ polish and the hard polish (μ, book and floor smoothings to 0) converged at 1.5e-8.

**Two bugs found and fixed on the way:**
1. **`calibrate_household_anchors` read the SS with the TPI on.** In the GZ form, x = 0 puts the bond price on the
   floor (0.875 instead of 0.947). The anchor then gave D households a lump-sum transfer of 0.79% of quarterly GDP,
   baked into the hard-polished TPI rules (rest-point C_D +1.8%).
   - It now reads with the TPI off. The final polish re-ran from the pre-polish checkpoint.
2. **`tpi_experiment.rest_point_row` read the fitted rules.** Off the nodes, the GZ purchase rule blurs the kink, and
   the read showed a spurious 23 bp fall in the rest spread. It now reads exactly.

**Results** (2% shock, exact reads each quarter; TPI economy solved at smoothing 1e-4 with every other condition hard):

| | no TPI | TPI |
|---|---|---|
| rest-point spread | 68.5 bp | 75.3 bp |
| rest-point n_D (D bank net worth), α_D (franchise value) | 1.498, 1.116 | 1.479, 1.105 |
| impact: sovereign spread | 223 bp | 200 bp (cap holds) |
| impact: purchases | — | 41.4% of the D debt |
| impact: D bond price vs no-shock path | −6.0% | −4.9% |
| impact: D bank net worth | −5.8% | −7.9% |
| impact: investment, GDP | +0.10%, +0.009% | −1.11%, −0.056% |
| quarter 4 / 24: spread | 133 / 73 bp | 141 / 81 bp |
| quarter 24: GDP | −0.007% | −0.114% |
| Bocola split of the excess-return change at impact | — | risk premium −142 bp, liquidity premium 0 |

**What the run shows:**
- **The TPI works where it can.** It holds the spread at 200 bp and supports the bond price by 1.1 pp. It does so only
  by taking 41% of the debt, and only through the risk premium: the D bank is unconstrained, so there is no liquidity
  premium to remove and no lending channel (lending spread 0 throughout).
- **It does not relieve the economy.** D bank net worth falls more on impact, and investment, GDP and consumption fall
  instead of rising.
  - Spreads run above the no-TPI path after quarter 2. The rest-point spread is ~7 bp higher, and the bank's
    franchise value lower.
  - Along the no-default path the Eurosystem earns the bonds' excess return (P&L +0.72% of quarterly GDP at q1),
    93% of it remitted to F. At 200 bp on this box the cap lies above the fair (expected-loss) spread everywhere
    (129 bp at 2%), so the Eurosystem buys cheap relative to risk-neutral value.
- **Not robust yet.** The post-impact dynamics run through the Eurosystem-book state, which has 3 grid points, and the
  one-off purchase at a node moved with the box width (44% → 71% → 76% of the debt at the centre node). Read them as
  directions to check, not results.
  - The natural next step is more nodes in M (Smolyak level 2 in that dimension), and a decomposition of the
    investment response.

**Reproduce:** `python3 run.py` (with the cache and `results/GLOBAL_tpi_walk` it reuses everything). Console:
`Claude files/diagnostics/tpi_root_cause/runs/run_box2pct_final.out`. Results: `results/GLOBAL/data/tpi.*`;
figures: `results/GLOBAL/figures/tpi_irf.png`, `tpi_mechanism.png`. The 6.2%-box results are kept in
`results/GLOBAL_box6.2pct_2026-11-04`.

## 7. Why the TPI economy does worse after impact (2026-11-06)

Three of the reported gaps are computational. The fall in D output after impact is the model's economics, and it does
not depend on the cap's own rules.

### 7.1 The rest point: a bug in how it is found

`stochastic_rest_point` iterates the fitted reads (`read_at`). In the TPI economy's GZ form, that read takes `x_cb` (the
purchase, or the price gap above the floor) from its fitted rule, and `x_cb` *sets* the bond price. The rule
interpolates a kinked function: purchases of 0.2–0.9 of the debt at the 23 top-risk nodes, price gaps near −0.05
elsewhere.

| At the fitted-read rest point | Fitted read | Exact read |
|---|---|---|
| TPI: largest residual | 1.2e-2 | 1e-14 |
| TPI: sovereign spread | 46.4 bp | 75.3 bp |
| no TPI: sovereign spread | 69.8 bp | 68.5 bp |

With the price 1.2% too high, the Bohn rule (γ_τ 0.015, the tax response to debt) settles the debt where that price
implies: B/B_ss = (γ − δ(1−Q_ss))/(γ − δ(1−Q)) = 1.081, so 2.008, the 2.010 found. At that state the exact read has
the debt still growing.

Solved instead as S = advance(read_exact(S)), Newton on the 10 endogenous states, |S − g(S)| ≤ 4e-13:

| Exact rest point | no TPI | TPI | Gap (was) |
|---|---|---|---|
| Sovereign spread | 69.1 bp | 75.3 bp | +6.3 bp (+6.7) |
| n_D, D bank net worth | 1.5123 | 1.5141 | +0.12% (−1.3%) |
| D debt | 2.105 | 2.130 | +1.2% (−4.7%) |
| C_D, Y_D, I_D | | | −0.06%, −0.09%, −0.23% |

The impulse responses from the exact rest points barely move: impact net worth −5.8% / −7.7%, investment +0.09% / −0.99%,
GDP +0.009% / −0.050%.

### 7.2 The spread gaps: interpolation leakage

**What economics allows.** The cap binds only at the top-risk nodes (pd ≥ 1.6% at the box centre). From the rest point
(pd 0.1%), next quarter's risk spans pd 0.03–0.33%. The chance of reaching the cap region within the bonds' life
(~16 quarters) is about 1%. So the TPI's true effect at rest is a fraction of a bp, and downward: it is a put on
future prices.

**What the rules show:**
- At the rest point's other states, the TPI rules price D bonds 0.25% below the no-TPI rules at the rest pd (0.10%),
  0.46% below at 0.28% and 0.54% below at 0.75%. Even the lowest-risk nodes (pd 0.003–0.08%) differ: α_D (franchise
  value) +0.23 to +0.34%.
- At the 23 top-risk nodes the TPI's effect varies wildly by direction: α_D from −18% (at P_D = +1) to +39% (at
  Z_D = −1), Q_bD from −1.3% to +3.9%.
- The risk × state cross terms are quadratic in s through three points (pd 0.003%, 0.08%, 2%). They carry each
  top-node difference to pd 0.28% with weight 0.26 and to pd 0.75% with weight 0.60.

**Not a one-step error.** Solving next quarter exactly leaves the rest gap at +5.9 bp: the leakage is in the whole
fixed point.

**After impact,** measured as changes from each economy's own rest point: the TPI economy's spread runs about 3 bp
above no-TPI at quarters 4–12. A one-off 41% purchase with the no-TPI rules (run C below) runs 2–6 bp *below* no-TPI,
as it should.

### 7.3 The real effects after impact: the model's economics

**The cap's rules are not needed.** A one-off purchase of the same 41% at impact, then the no-TPI rules (no cap in
them), reproduces the TPI economy. Values are % vs each economy's unshocked path; spreads are bp vs its rest point:

| Impact (quarter 8) | A: no TPI | B: TPI | C: no TPI + one-off 41% | D: one-off 41%, no shock |
|---|---|---|---|---|
| Spread | +156 (+32) | +125 (+35) | +138 (+30) | −5 (−2) |
| D bank net worth | −5.8 (+1.0) | −7.7 (−0.4) | −7.5 (−0.7) | +0.4 (−0.4) |
| Investment | +0.09 (+0.04) | −0.99 (−0.53) | −0.75 (−0.41) | +0.06 (−0.02) |
| GDP | +0.009 (+0.005) | −0.050 (−0.083) | −0.038 (−0.075) | −0.007 (−0.007) |
| D consumption | +0.06 (−0.02) | −0.16 (−0.29) | −0.18 (−0.36) | −0.15 (−0.11) |

**Robust to next quarter's rule error.** The impact was re-solved with next quarter cleared exactly at all 5 risk nodes
× 2 regimes, using the fitted rules only from t+2. B moves from GDP −0.050 to −0.048%, investment −0.99 to −1.06% and
net worth −7.7 to −7.9%. C keeps its sign: investment −0.75 → −0.53%, net worth −7.5 → −6.6%. A is unchanged.

**The chain:**
1. **The Eurosystem holds bonds bought at a crisis price.** Along the no-default path it earns the default premium: a
   P&L (`Pi_cb`) of 0.71% of quarterly D GDP in quarter 1 and 4.87% over 25 quarters. The capital key (key_D 0.071)
   sends 4.53% of it to F.
   - In the default branch next quarter the Eurosystem loses 36% of D GDP, of which F bears 33.5%.
   - The contract is actuarially fair (0.98 × 0.66 ≈ 0.02 × 33.5). The no-default IRF shows only the premium.
2. **D bank net worth falls through its capital book.** The bond support adds 1.1 pp; the capital price costs 3.0 pp.

   | Impact net worth, pp | Capital book | D bonds | Total |
   |---|---|---|---|
   | A: no TPI | +0.24 | −5.97 | −5.79 |
   | B: TPI | −2.78 | −4.86 | −7.71 |

3. **The capital price falls because capital's required return rises.** The linearised no-default capital Euler splits
   the impact −0.50% (B) exactly:

   | Impact capital price (%) | Deposit-rate path | Capital risk premium | Marginal product | After q24 | Total |
   |---|---|---|---|---|---|
   | A: no TPI | −0.05 | +0.12 | −0.00 | −0.02 | +0.04 |
   | B: TPI | −0.06 | −0.43 | +0.03 | −0.03 | −0.50 |
   | C: one-off 41% | −0.11 | −0.24 | +0.02 | −0.04 | −0.38 |

   - The banker's kernel becomes steep in risk. Next quarter's α_D (marginal value of bank net worth) is 1.38–1.41 at
     the two highest-risk nodes against 1.21 at the others; without the TPI it is 1.16 against 1.12–1.14. Capital pays
     poorly exactly in those states.
   - Without the cap in the rules (C) the premium is −0.24 pp. B's extra sits on the top-risk states where the TPI
     rules are least accurate (§7.2).
   - The D deposit rate, including the premium on D's net foreign debt (`kappa_nfa`), explains only 0.06–0.11 pp.
4. **Investment follows the capital price.** It falls 1% (Jermann ξ 0.5 gives investment an elasticity of 2 to the
   capital price).
5. **Impact GDP moves only with hours, since capital is predetermined.** Lower D demand cheapens D goods: `p` (the
   price of F goods in D units) rises 0.54%, which raises D's consumer price index and cuts the real consumer wage.
   With GHH preferences hours respond only to that wage: −0.07%, so GDP −0.05%. Later capital falls (−0.22% by q24)
   and GDP reaches −0.11%.
6. **Nothing offsets it.** The D incentive constraint is slack at the rest point and along both paths (`mu_D`, the
   constraint multiplier, is 0). It binds at 52 of 75 box nodes, all off the path, so the purchases free no lending.
   The 2% shock alone causes no recession (A: GDP +0.009%), so the TPI has nothing to repair and only its side
   effects show.
7. **D households gain from default in this calibration.** Their consumption in the default branch is 0.5% above the
   no-default branch without the TPI, and 2.5% above with it. Insurance against default moves resources towards the
   states where D's marginal utility is lowest.

### 7.4 What remains uncertain

- **The Eurosystem-book dimension.** It has 3 nodes, and the lowest is unphysical (the Eurosystem short, D banks
  holding 175% of the debt, the constraint binding in default). At the post-impact states the default-branch rules
  miss the exact solve by +20% (α_D), −0.86% (C_D), +1.1% (p) and −160 bp/yr (rdep_D, the D deposit rate). One-step
  exact solves move the impact responses by 0–30%; errors further out are not tested.
- **The impact state sits on the box's top edge** (pd 2% = `PD_BOX_TOP`). The 23% of next quarter's risk mass that
  lies above the box is read at its edge.

### 7.5 Options (author decisions)

1. **Fix the rest point:** solve it as the fixed point of the exact period map. Cheap; changes the rest table.
2. **Reduce the leakage and the book error,** then re-solve both economies (~2–4 h). Put the M box on [0, M+] so all
   three M nodes are physical, add risk nodes to the cross terms, and move the box top above the shock.
3. **Economics.** On the 2% box the TPI can only shift risk. For it to help D, it must act where the D constraint binds:
   a larger shock or a tighter calibration, which brings back the corners of §2. The remittance is also a choice:
   key_D 0.071 makes the TPI cross-border default insurance, while key_D 1 isolates its domestic effect.

**Reproduce** (all read-only, on the saved rules; logs in `Claude files/diagnostics/tpi_root_cause/runs/`):
`post_impact.py` (rest states, same-state reads, rule error along both paths), `exact_rest.py` (exact rest points and
both IRFs from them), `one_off_purchase.py` (runs A–D), `regime_anatomy.py` (the banker's kernel by regime, the
default-branch rules), `one_step_exact.py` (impact with next quarter solved exactly), `leakage_profile.py` (the TPI
rules against the no-TPI rules along s and at the top nodes), `capital_price_decomp.py` (the capital-price and
net-worth splits).

## 8. A shock that causes a recession, and what the TPI can then do (2026-11-06)

The 2% experiment has no recession to repair (§7.3). This section moves to crisis-sized shocks. All reads are exact,
with the rest point solved as a fixed point of the exact period map (now the library default, §8.5).

### 8.1 Where the model has a recession

- **2% box and 6.2% box: nowhere.** On the 6.2% box the D constraint stays slack on impact up to its top (μ_D, the D
  constraint multiplier, is 0 at 6.2%). GDP moves −0.02% to +0.02%, and consumption and investment rise.
- **Appendix-B 12% box (`data/global_warm_start.pkl`, fix A): from about 5%.** No-TPI impact against the unshocked path:

| Quarterly default probability | Spread rise | D bond price | D bank net worth | GDP | Investment | D consumption | μ_D | Lending spread |
|---|---|---|---|---|---|---|---|---|
| 4% | +356 bp | −12.8% | −16.8% | −0.036% | −1.44% | +0.22% | 0 | 0 |
| 5% | +430 bp | −15.1% | −19.9% | −0.080% | −1.72% | +0.19% | 0.019 | +81 bp |
| 6% | +505 bp | −17.2% | −22.9% | −0.142% | −2.03% | +0.13% | 0.046 | +197 bp |
| 8% | +646 bp | −21.0% | −28.0% | −0.232% | −2.45% | +0.05% | 0.088 | +368 bp |
| 9% | +715 bp | −22.8% | −30.1% | −0.266% | −2.61% | +0.02% | 0.105 | +434 bp |

- **The constraint binds only on impact.** From quarter 1 it is slack again, GDP sits near −0.05% and investment
  recovers slowly (−0.26% at q12). The 2010–11 Greek average, 9.1%/qtr, lies in this region.
- **All acceptance signs hold except D consumption.** The D household's budget at 9% shows why
  (`consumption_sign.py`): on impact its income falls 4.9% of quarterly GDP (bank dividends −3.4%) and the real value of
  its deposits 2.4% (consumer prices +0.25%). It runs its deposits down by 7.3% of GDP (12–13% from quarter 1), helped by
  a 180 bp fall in the deposit rate, and keeps consumption flat. The loss is temporary, so a representative household
  with a deposit Euler smooths it; only hand-to-mouth households would cut consumption.

### 8.2 What purchases can do there

**A 200 bp cap is out of reach.** At 9% the spread is still ~700 bp after the purchase that switches the constraint
off. That remainder is expected default loss, which no purchase can remove.

**Purchases relieve the constraint.** One-off purchase at impact, rules held fixed (`crisis_shock.py`):

| Default probability | Purchase that switches the constraint off | Impact GDP: no purchase → with it |
|---|---|---|
| 5% | 7.7% of the debt | −0.080% → −0.040% |
| 6% | 11.8% | −0.142% → −0.045% |
| 8% | 14.9% | −0.232% → −0.050% |
| 9% | 15.6% | −0.266% → −0.052% |

**A relief rule cannot be unconditional.** On the 12% box the D constraint binds at 54 of 75 grid points
(`node_relief.py`), including the low-risk box centre (μ_D 0.01–0.03 at default probabilities of 0.003–0.2%).
- The box is centred on the deterministic steady state, where the constraint binds. The stochastic rest point lies far
  from it (b_DD +15.6%, K_D −2.6%) and is about 10% slack.
- Relieving the constraint at every binding point would buy a median 25% and up to 57% of the debt, in calm states too.
- So any relief rule must be switched on by the sovereign crisis, not by the constraint alone.

### 8.3 After impact: the remittance decides the sign

- **Capital key (0.071).** After the 15.6% relief purchase at 9%, the Eurosystem's no-default profit is 1.7% of
  quarterly D GDP in quarter 1, of which 93% goes to F.
  - D GDP is −0.144% at q12 against −0.069% without any purchase.
  - D consumption is −0.46% at q8 against −0.08%, and D bank net worth +2.2% against +7.7%.
- **Profit and loss kept in D (tpi_key_D 1):** see §8.4.
- **The common thread.** Impulse responses condition on no default, so whoever holds the D bonds collects the whole
  default premium as the price recovers.
  - Without the TPI that is the D banks: net worth overshoots, and dividends lift D consumption (§8.1).
  - With the TPI it is the Eurosystem, which under the capital key passes it to F.
  - In the default branch the capital key works the other way: F would bear 93% of the Eurosystem's haircut.

### 8.4 The remittance decides how purchases work, not only who gains

The no-TPI rules were re-solved with the Eurosystem's profit and loss kept in D (`tpi_key_D` 1). The direct jump
stalled at max|F| 0.14 (F's default-branch Euler at the book edges); the continuation 0.071 → 0.25 → 0.45 → 0.7 → 1
converged at every rung (last 1.3e-9; `key_homotopy.py`, rules in `runs/rules_12pct_key1.pkl`).

**On impact (9%, rules held fixed):**

| Purchase | Capital key: μ_D | Bond price | Bank net worth | Kept in D: μ_D | Bond price | Bank net worth |
|---|---|---|---|---|---|---|
| 0 | 0.105 | | | 0.105 | | |
| 15.6% | 0.000 | +2.3% | +7.2% | 0.095 | +0.4% | +2.5% |
| 40% | 0 | +5.1% | +9.9% | 0.036 | +1.8% | +8.8% |

- **Under the capital key, purchases relieve the D banks mainly by moving D's default risk to F,** which bears 93% of
  the Eurosystem's haircut in default. That lifts the bond price, bank net worth and the banks' expected discount factor
  (+2.6%) enough to switch the constraint off.
- **Kept in D, the risk never leaves D.** The price barely moves, the discount factor falls 3%, and the swap of bonds
  for non-divertable reserves alone frees little room.

**After impact** (GDP % vs the unshocked path at impact / q1 / q4 / q8 / q12):

| Path | GDP |
|---|---|
| No TPI | −0.266 / −0.047 / −0.059 / −0.064 / −0.069 |
| Capital key, buy 15.6% | −0.052 / −0.072 / −0.106 / −0.128 / −0.144 |
| Kept in D, buy 15.6% | −0.237 / −0.031 / −0.037 / −0.032 / −0.028 |
| Kept in D, buy 40% (past the 28% book edge: continuation clipped, indicative) | −0.082 / −0.047 / −0.010 / +0.002 / +0.014 |

- Kept in D, investment, the capital price and the sovereign spread are also better than without the TPI at every
  horizon.
- D consumption is 0.05–0.08 pp lower from q4: the income now reaches households through lower taxes (Bohn rule)
  rather than bank dividends.
- **The trade-off.** The capital key buys strong impact relief with cross-border risk sharing, and pays for it with the
  premium F collects afterwards. Kept in D, there is no transfer, but relief needs purchases ~2.5× larger, beyond the
  current book box.
- **A third design, not yet solved:** income to D, losses by key (the 2012 Greek arrangement, under which SMP income
  was returned to Greece). It should combine the capital key's impact relief with no premium flowing to F, at the
  cost of an expected subsidy from F.

### 8.5 Rest-point fix (in the code, 2026-11-06)

- `recursive_experiment.stochastic_rest_point` now polishes the fitted iteration's end on the exact period map
  (`_exact_rest`, Newton on the 10 endogenous states) and caches it separately.
- The solver's box re-centring and the calibration moments keep the fitted point (`exact=False`), so neither the
  solver path nor the calibration targets move.
- It reproduces the diagnostic's exact rest points to 0.0 on both 2%-box economies.

## 9. The implemented experiment: a 6% crisis shock and the TPI as a relief intervention (2026-11-06)

**The author's choices:** shock 6%/qtr, capital key 0.071, and the TPI as a relief intervention: at the shock,
unanticipated, the Eurosystem buys just enough D bonds to switch the D banks' constraint off and holds them to
maturity, and the no-TPI rules continue. `run.py`: `RISK_SHOCK_PD` 0.06, `PD_BOX_TOP` 0.12, `TPI_M_BAND` 0.32,
`TPI_MODE` "relief"; `tpi_experiment.intervention`. The full run takes 17 min plus figures (console in
`diagnostics/tpi_root_cause/runs/run_crisis6pct.out`; results in `results/GLOBAL`). Every path below is cleared
exactly each quarter.

| % vs the unshocked path (spreads bp/yr) | Impact | q1 | q4 | q8 | q12 | q24 |
|---|---|---|---|---|---|---|
| GDP, no TPI | −0.142 | −0.045 | −0.053 | −0.055 | −0.058 | −0.073 |
| GDP, TPI | **−0.045** | −0.057 | −0.077 | −0.090 | −0.099 | −0.121 |
| Lending spread, no TPI / TPI | 188 / **0** | 0 / 0 | | | | |
| Sovereign spread (level), no TPI / TPI | 567 / 532 | 421 / 403 | 219 / 212 | 128 / 127 | 95 / 96 | 68 / 70 |
| D bank net worth, no TPI / TPI | −22.9 / −20.0 | −11.2 / −12.9 | +1.4 / −1.8 | +5.7 / +2.4 | +5.7 / +3.0 | +3.0 / +1.8 |
| Investment, no TPI / TPI | −2.03 / −1.37 | −1.06 / −1.11 | −0.45 / −0.64 | −0.20 / −0.46 | −0.20 / −0.46 | −0.34 / −0.51 |
| D consumption, no TPI / TPI | +0.13 / +0.07 | +0.15 / −0.04 | +0.01 / −0.22 | −0.05 / −0.29 | −0.06 / −0.29 | −0.05 / −0.22 |

- **Impact: the TPI works as it should.** 11.8% of the debt switches the D constraint off (μ_D 0.046 → 0), which removes
  the 188 bp lending spread. It cuts the GDP fall by two thirds and lifts the bond price, bank net worth and investment.
  - The bond's expected excess return falls 291 bp: −197 bp liquidity premium, −94 bp risk premium.
- **After impact: worse under the capital key.** The Eurosystem's no-default profit on the crisis-priced bonds is
  0.98% of quarterly D GDP in quarter 1 and 3.66% over 25 quarters, of which 3.40% goes to F.
  - The D banks no longer hold the recovering bonds, so their net worth recovers more slowly.
  - GDP is below the no-TPI path from quarter 1: −0.81 pp-quarters over 25 quarters, roughly break-even over the
    first four.
  - Kept in D (§8.4), GDP after impact is above the no-TPI path at every horizon, but the same purchase barely
    relieves the constraint. "Income to D, losses by key" is untested.

**Fitted and exact reads disagree at crisis states.** Step 4 of `run.py` (the risk IRF figure, the persistence IRF,
the impact table, the output decomposition) reads the fitted rules, by an earlier documented choice
(`output_decomposition._refine`: exact re-solves slipped onto the slack branch at a barely binding SS). At the 6% shock
the fitted read has residuals of 4–6e-3, and it doubles the impact recession:

| No TPI, 6% shock | GDP impact / q1 / q4 / q24 | Credit spread on impact | Deposit rate on impact |
|---|---|---|---|
| Fitted reads (step 4) | −0.309 / −0.145 / −0.061 / −0.098 | 311 bp | −44 bp/yr (level) |
| Exact reads (step 6) | −0.142 / −0.045 / −0.053 / −0.073 | 188 bp | +36 bp/yr |

- Here the exact read stays on the binding branch (no slip). The risk IRF figure and the TPI figure therefore show
  different no-TPI recessions for the same shock.
- The output decomposition's "labour-FOC residual" (−0.12%) is this gap.
- The cure is a more accurate solution: a box re-centred on the rest point (the 12% box is centred on the
  deterministic SS, 15.6% off in the bank's bond book), or more nodes in the bank states. Until then the exact read is
  the one consistent with the TPI comparison.

**Resolved the same day (author: "solve the model exactly"):** every readout of `run.py` now clears the period map
exactly at each visited state.
- What switched: the rest-point and impact tables, the persistence and risk IRFs, the TFP IRF, the calibration table,
  and the output decomposition (`simulate(exact=True)`).
- What stays fitted: the impact table keeps one fitted column as the accuracy check, and the Euler-error report still
  measures the fitted rules.
- The calibration routine keeps its fitted reads, so nothing is re-calibrated.
- **Results:**
  - The risk IRF and the TPI comparison now show the same no-TPI path (impact GDP −0.142%).
  - The output decomposition closes without its residual: lending spread −0.149, deposit rate +0.070, terms of trade
    −0.063, labour-FOC residual 0.000, total −0.142%.
  - The TFP IRF barely moves (impact GDP +1.281% vs +1.285% fitted).
  - The calibration table's transmission rows are now read at p^d 2% as labelled: 0.7 bp against the 108 bp target.
    The 337.5 bp before was the 6% shock's fitted response under a 2% label.
- Run: `diagnostics/tpi_root_cause/runs/run_crisis6pct_exact.out` (76 s on the cached solves). The fitted-read run is
  kept in `results/GLOBAL_crisis6pct_fitted_2026-11-06`.
