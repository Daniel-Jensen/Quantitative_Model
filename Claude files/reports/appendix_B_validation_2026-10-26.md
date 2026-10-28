# Appendix B (`B-calibration.tex`) against the global solver — validation, solvability, fix plan

**Summary.**
- **Phase 1:** the code was an older calibration vintage, with ~15 parameters off; all are now ported.
  The model equations match the text, apart from one Walras-leaking bug (F taxes, fixed) and four documented
  departures the text should state (smoothed KKT switch, now hard; union deposit market; representative-agent
  dynamics; the 15% floor).
- **Phase 2:** the appendix calibration's steady state solves and reproduces every closed form. The global
  model does not solve there. Its first stage is near-singular on the coarse grid: cold starts fail, and a
  continuation path breaks short of the appendix (at the F coupon decay, with G, Bohn and the wedge still to go).
- **Phase 3:** Tier-1 remedies are tested and ranked. The remaining fixes are calibration choices (Tier 2) that
  need your decision.

2026-10-26, branch `calibration`. Solver: the global projection model in `global_projection/`; the SSJ
pipeline is not involved. Rates are annualised unless marked "/qtr"; default probabilities are quarterly.

## Phase 1 — calibration and equations

### 1.1 Parameters (appendix → code before today → code now)

| Parameter (code key) | Appendix B | Code before 2026-10-26 | Code now |
|---|---|---|---|
| Relative size κ (`size_F`) | 11.3 | 8 | 11.3 |
| Risk aversion σ (`sigma_X`) | 1.35 | 1.0 | 1.35 |
| Frisch ν (`frisch_X`) | 0.5 | 2.0 | 0.5 |
| ρ_e, σ_e, a_min | 0.9, 0.2, 0 | same | same |
| Trade elasticity η (`epsilon_trade`) | 1.5 | 0.5 | 1.5 |
| Home bias ϖ_D (`omega_home_D`) | 0.68 | 0.85 | 0.68 (ϖ_F = 0.9717 derived) |
| α, δ, ξ, ζ | 0.30, 0.025, 0.50, 1 | same | same |
| Retail elasticity ε (`epsilon_X`) | 11 | 6 | 11 |
| Real deposit rate r̄ | 0.003 | same | same |
| Payout share f | 0.08 | same | same |
| Banker discount β_int | 1/(1+r̄) = 0.997009 | 0.997 | 1/(1+r̄), derived |
| Foreign share of D debt (`foreign_share_D`) | 0.127 | 0.20 | 0.127 |
| Cross-border cost ψ | 2 | same | same |
| ρ_s | 0.9 | 0.95, **hard-coded in `state_grid.s_process_params`** | 0.9, in the calibration |
| σ_s | 0.41 | 0.4455 | 0.41 |
| Recovery ϱ_D | 0.45 | same | same |
| π^d at rest (`pd_rest`) | 0.10% | 0.02% | 0.10% |
| Gauss–Hermite nodes O (`N_GH`) | 5 | same | same |
| G/GDP (`G_X`) | 0.204; 0.191 | 0; 0 | 0.204; 0.191 |
| Coupon decay δ_b (`delta_b_X`) | 0.056; 0.039 | 0.056; 0.056 | 0.056; 0.039 |
| Bohn coefficient, annual | 4γ_τ = 0.04 | root 0.93 targeted (4γ ≈ 0.28) | `bohn_annual` 0.04 (the primitive; root 0.9926 reported) |
| Leverage target | 8 | 8 | 8 |
| Wedge r̄^k − r̄ | 0.05% p.a. (5 bp) | 1% p.a. (100 bp) | 5 bp (rest-point target and starting instrument) |
| Own sovereign / assets | 14.4%; 6.5% | same | same |
| KKT switch | hard max (eq. mu-closed) | smoothed, mu_eps 1e-2 | hard (`mu_eps` 0; solved smoothed first, hard polish last) |

**Code parameters the appendix does not list** (kept by your decision 2026-10-26, each needs a line in the
appendix): the D banks' share of the F stock (`foreign_share_F` 0.5%); the funded net-worth floor in default
(`nw_floor_frac` 15%); the debt-elastic premium on the cross-border deposit position (`kappa_nfa` 0.01);
the no-short-sale complementarity on cross-border books (`book_eps`, now hard); the floor's max
(`floor_eps`, now hard); Z_D's deterministic AR(1) (`rho_z` 0.9, TFP experiment only); the guards
`alpha_cap` 40 and the μ cap 0.95 (never active at a solution).

**Steady state at the ported calibration** (`solve_steady_state`, 9 s, every `_check_targets` assertion
passes). The appendix's closed forms reproduce exactly:

| Object | Appendix | Code |
|---|---|---|
| ς̄θ̄ | 0.0010 | 0.00100 |
| μ̄ | 0.0010 | 0.00100 |
| φ̄ | 1.013 | 1.01261 |
| λ | 0.127 | 0.12658 |
| ω | 0.0095 | 0.00954 |
| r̄^k | 0.0031 (1.25%) | 0.003125 (1.250%) |
| q̄^D | 0.947 | 0.9471 |
| Duration | 4.2 years | 4.24 years |
| Debt root ρ_B | 0.993 | 0.99257 |

The appendix's TBC entries at these deterministic-SS instruments (they move once the rest-point loop
changes the instruments, see Phase 2):

| Object | Value | Data |
|---|---|---|
| β_D; β_F | 0.992412; 0.992219 | |
| Bank-held D debt / annual output | 0.50 | EBA-held 0.40 |
| D bank assets / annual output | 3.03 | 1.7 |
| Liquid deposits / GDP, D | 2.65 | 0.84 |
| p_ss | 0.98868 | |

### 1.2 Equations

Checked line by line against `02-model.tex` (eq. numbers by label) and the appendix:

| Equation | Code | Verdict |
|---|---|---|
| Jermann coefficients (jermann-coeffs-app), Q and ι (tobins-q), capital law | `blocks/capital.py:5-44` | match |
| Labour demand and mpk with working capital (labour-demand), loan demand | `blocks/firms.py:20`, `point_map.py:218` | match |
| Return on capital (rk) | `point_map.py` `_regime_continuation` | match |
| Balance sheet, asset payoff, bond payoffs with the haircut h on coupon and continuation | `point_map._balance_sheets`, `_bond_demands` | match |
| P-state, gross wealth n^g | `point_map._balance_sheets` | match |
| Ω (omega, branch-sdf) on the GHH composite | `point_map._discount_kernels` | match |
| IC, μ closed form (mu-closed), φ recursion (phi-recursion) | `point_map._incentive_constraint` | match; the max was smoothed until today |
| Capital Euler, bond pricing (bond-pricing) | `point_map._incentive_constraint`, `_tpi_rule` | match |
| Cross-border FOCs with p′/p and the friction (bondD/F-foreign) | `point_map._bond_demands`, `_residual_vector` | match; the friction's reference equals b̄ when spreads are equal (the appendix case) |
| r^wc (rwc) | `point_map._wages_and_dividends` | match |
| Net-worth law and bank dividends (nw-lom, bank-div), household dividends (dividends) | same | match (working-capital receipts pass through the household budget and its claim) |
| s process, logistic π^d, the expectation (s-process, pd, expectation) | `point_map._next_states`, `_continuation`, `gh_nodes` | match (probabilists' rule) |
| Government budget, debt law, Bohn rule with the default anchor ϱB̄ (govt-budget, debt-lom, bohn, debt-root, gamma-root) | `point_map._government`, `blocks/government.py` | match for D; **F was wrong (E1)** |
| SS closed forms (ss-closedform-app) | `blocks/bank.calibrate_bank_targets` | match (numbers above) |
| CES baskets, import demand, NX, home bias and size (ces-D/F, nx, home-bias-size) | `blocks/trade.py` | match |

**Discrepancies**

- **E1 — a Walras leak, fixed today.** The text fixes the F stock and lets F taxes close the F budget.
  The code held F taxes at their SS value (`point_map._government`), so F's budget leaked
  δ_F·B̄_F·(q_F,t − q̄_F) every quarter. On the saved 2026-10-24 run, goods_F (the F goods-market residual)
  was −1.2e-5 at rest and −3.7e-5 at the 2% shock, against the guide's 2e-6. Both equal the leak to 3 digits.
  The nesting test N4 had subtracted this term as a "known" one. Now `Tax_F = G_F + δ_F·B̄_F·(1 − q_F,t)` minus
  the Eurosystem remittance, plus the floor top-up. Re-clearing the saved states gives goods_F −2.5e-7, flat
  across states (that remainder is the saved run's household anchor). N4 now passes without the term.
- **E2 — the KKT switch.** The code smoothed max(z, 0) with eps 1e-2. At the appendix's 5 bp wedge, μ̄ is
  0.0010, so this smoothing is ten times the multiplier. It gave λ 0.1236 instead of 0.1266, and dμ/dz = 0.04 at
  the SS where the hard max gives 1. Now a genuine max, by your decision.
- **E3 — deposit market.** `02-model.tex` eq. dep-national clears deposits country by country. The code has a
  union-wide market with deposit-UIP, the V_dep state and kappa_nfa: a deliberate choice since 2026-07 (it is what
  restored the impact contraction). The text should describe the union market.
- **E4 — households in the dynamics.** The text and panel B describe heterogeneous households with their own
  β_X. The global solution's household block is a representative-agent closure with β_eff = 1/(1 + r̄)
  (`point_map._households`). The HA block and β_X shape only the SS distribution. The text should say so.
- **E5 — the net-worth floor.** The 15% funded floor (`point_map._funded_floor`: D tops up with new debt, F with
  a lump-sum tax) is not in the text. It bound in default at ~22 of 25 nodes on the 2026-10-20 solve (old
  calibration); at appendix B the default regime has not been reached (Phase 2).
- **E6 — productivity.** The text has a constant Z_X; the code carries Z_D as a deterministic AR(1) state for the
  TFP experiment. It is inert in the risk experiment.
- **E7 — the debt-stability statement.** "Debt is stable only if γ_τ > δ_b(1 − q̄)" is exact at q̄, but the root
  (1 − δ_b) + (δ_b − γ_τ)/q exceeds 1 whenever q < 1 − γ_τ/δ_b = 0.821. In the saved 2026-10-24 run (old
  calibration) every shock above p^d ≈ 0.55% priced the D bond below that, so at 4γ = 0.04 a comparable crisis
  makes the debt rule locally explosive.

### 1.3 Code changes made today

| Change | Where |
|---|---|
| F taxes close the F budget at today's price (E1) | `point_map._government` |
| ρ_s read from the calibration | `state_grid.s_process_params` |
| s box 2.75 → 3.65 unconditional sd, so it reaches p^d 3%/qtr (2.75 stopped at 1.31% under ρ_s 0.9, σ_s 0.41) | `state_grid.S_COVER_SD` |
| Floor smoothing as a calibration key (`floor_eps`; 0 = hard) | `point_map._funded_floor` |
| Genuine conditions: every stage runs at working smoothings (`HARD_WORK`: mu_eps 1e-4 as the walk's floor, book_eps 1e-4, floor_eps 1e-3), then one hard polish; if that fails, hardened one at a time (IC, books, floor), and the rules record what they solve (`rules.smoothing`) | `recursive_experiment._soften`, `_hard_polish`, `solve_recursive`, `solve_no_default` |
| Bohn coefficient as the primitive | `blocks/government.py`, calibration `bohn_annual` |
| Appendix calibration ported | `calibration/global_projection.py` |
| Rest-point loop: spread instruments from 1 bp/yr (they stopped at 40 bp) | `calibrate_stochastic.PAIRS` |
| Tests tied to the calibration rather than to Bocola's numbers; N4 drops the old F-treasury term | `test_state_grid`, `test_ss_identities`, `test_recursive_nesting` |

**Tests at the ported calibration:**
- Pass: `test_ss_identities`, `test_fast_kernels`, `test_state_grid`, `test_no_unbound_names`,
  `test_calibrate_stochastic` (C1–C3), `test_cross_border`, `test_post_default`, and nesting gates N1, N3, N4
  (3.8e-10) and N5.
- `test_collocation`'s real solve fails at max|F| 2.75e-6. It failed before today as well (7.7e-4, worklog 2026-10-23/24).
  The test solves d0 straight, without the solver's ladders.
- **N2 fails** at max|F| 4.1e-5 (69 min). This is the grid-wide π = 0 solve through the solver's own ladder. The
  smoothing walk could not leave eps = 1e-2 (it halved down to a step of 1.02 and failed). The hard polish from
  1e-2 stalls on the IC. N2 also failed before today (6.7e-6). See Phase 2.

## Phase 2 — solvability

### 2.1 The machinery

- **Method:** time-iteration warm starts, then Newton on the whole collocation system: 20 residuals per node
  (14 equilibrium conditions plus 6 rule identities), on a level-1 Smolyak grid (25 nodes) over 12 states
  and 2 regimes. A finite-difference Jacobian runs on every core, with a line search and an LM fallback.
- **Acceptance:** sum F² ≤ m·(1e-9)².
- **The ladders:**
  1. the working-capital flow (0 → 1) and the cross-border valuation p′/p (0 → 1);
  2. the haircut (0.85 → 0.45);
  3. the joint solve;
  4. the KKT smoothing walk (1e-2 down);
  5. the instrument walk from a 100 bp economy;
  6. the s refinement (axis, then cross).
- **The genuine max:** a hard polish at the end, added today.

### 2.2 What was measured today

**Deterministic SS:** solvable at the appendix calibration (Phase 1).

**The coarse no-default stage d0 (π = 0) does not solve at the appendix calibration.** It is the first rung
of every ladder, so neither the risk solve nor the rest-point loop can start. Five setups, all on the
production V5 box:

| Setup | Result |
|---|---|
| SS with the genuine max (λ 0.1266), 12 warm sweeps, 15% floor (production) | the warm sweeps diverge (max\|F\| 2.5 → 0.25 after 2 sweeps → 20 after 8); Newton/LM stall at 12.6 on the first rung |
| same, no warm sweeps | Newton reaches 1.7e-2 in 19 steps, LM stalls at 5.3e-3 (goods_D, the D goods-market residual) |
| SS solved consistent with eps 1e-2 (λ 0.1236), no warm sweeps | stalls at 4.4e-3 (dep_clear, the union deposit-market residual) |
| no floor, 12 warm sweeps (the N2 setup) | first rung stalls at 6.6e-3 (id_alpha_F, the F franchise-value identity); N2's later rungs and the hard IC polish end at 4.1e-5 |
| hard max from any of these | inherits the stall |

The sweeps diverge for a concrete reason. The franchise recursion α = E[Λ(f + (1−f)α′)]R/(1−μ) is a contraction
only while μ < f = 0.08. The V5 box swings bank net worth by ±13–16% at its edge nodes, which puts μ at about 0.15
there. So iteration on α is locally expansive, and the 15% floor makes it worse.

**Which appendix values break it** (each one alone, cold-started from the pre-port calibration through the
production ladder; that calibration itself solves in 5 steps at 1.0e-8):

| Appendix value alone | Full move | Half move |
|---|---|---|
| G = 0.204/0.191 (from 0) | fails, 3.9e-4 (C_D 0.79 → 0.59) | converges, 6 steps |
| Bohn 4γ = 0.04 (from 0.275; debt root ≈ 0.93 → 0.99) | fails, 3.1e-5 | converges, 6 steps |
| δ_b^F = 0.039 | converges, 5 steps | – |
| κ = 11.3 (from 8) | converges, 6 steps | – |
| ϖ_D = 0.68 (from 0.85) | fails, 3.3e-4 | fails, 6.6e-5 |
| η = 1.5 (from 0.5) | fails, 1.6e-5 | crashes: η = 1 exactly divides by zero in `trade.ces_price` (a latent singularity at the Cobb–Douglas point) |

The sweep stopped at that crash. Foreign share, markup, σ, Frisch, banker β and the wedge were not tested alone.
The continuation path below covers them.

So no single parameter is the culprit. **The coarse d0 collocation tolerates only moderate parameter moves
from a solved point.** Cold starts at the appendix calibration fail because it is several such moves away from
any economy the ladder can solve from scratch.

**Tier 1 that works so far — a continuation path with SS-relative seeding** (`diag_path.py`): each rung's
rules are the new economy's SS plus the last solution's deviations from its own SS, in the rules' fit space.
A failed rung is halved. Results:

| Group (in path order) | Result |
|---|---|
| households: σ 1 → 1.35, Frisch 2 → 0.5 | one rung, 8 steps, 8.0e-10 |
| markup: ε 6 → 11 | one rung, 4 steps, 8.0e-9 |
| size 8 → 11.3 and η 0.5 → 1.5 (both broke the cold start alone) | one rung, 19 steps, 1.0e-10 |
| foreign share 0.20 → 0.127 | one rung, 6 steps, 1.0e-10 |
| home bias 0.85 → 0.68 | 1/16 rungs at first, then adaptive. **The endpoint floors at 6.03e-8 from two different rungs**; kept, flagged |
| δ_b^F 0.056 → 0.039 (converged alone from the old calibration) | from the home-bias endpoint: full 4.7e-7, ½ 1.7e-7, ¼ 1.0e-7, ⅛ 8.1e-8 (kept), 3/16 9.2e-8 (kept), then 1.2e-7 and 1.0e-7. **The path breaks at 0.19 of this move**: every further rung lands above 1e-7 |
| G, Bohn, banker β, wedge 100 → 5 bp | never reached |

**Is the floor an artifact?** No. From the 6.03e-8 point, a Newton step with a central-difference Jacobian
(`probe_cfd.py`, h = 1e-5) finds condition 1.5e10 and no descent. The same economy on a box with half the
bank-stock bands (`probe_box.py`) is worse, 8.4e-5. The coarse collocation system for d0 has a near-fold at the
appendix's open-economy values. Each further appendix move raises the attainable precision (1e-10 → 6e-8 → ~1e-7);
past a fifth of the F coupon-decay move, no rung gets below 1e-7.

### 2.3 Verdict

- **Code as of 2026-10-25** (old calibration): solvable on the production grid. The no-TPI cross solve reached
  2.0e-9, but the acceptance signs fail (Y_D rises on impact, the lending spread falls 106 bp). The cause is the
  franchise value, read off the saved rules (exact reads at the rest point, other states held):

  | p^d | n_D | E[Ω] = (1−p^d)·Ω^nd + p^d·Ω^d | α′ no default | α′ default | μ_D |
  |---|---|---|---|---|---|
  | 0.02% (rest) | 1.305 | 1.270 = 1.269 / 6.58 | 1.295 | 7.31 | 0.0247 |
  | 0.5% | 1.193 | 1.396 = 1.370 / 6.53 | 1.411 | 7.26 | 0.0093 |
  | 1% | 1.127 | 1.491 = 1.440 / 6.46 | 1.487 | 7.16 | 0.0022 |
  | 2% | 1.025 | 1.665 = 1.570 / 6.33 | 1.623 | 6.99 | 0.00065 |

  Net worth falls 21%, but E[Ω] rises 31%. The no-default kernel rises 24% because next period's franchise value
  anticipates future tightness. The default leg adds p^d × 6.3, because the D bank loses ≥ 63% of net worth in
  default (exposure/NW 1.15 × 55%) and sits on the 15% floor, where α ≈ λθ/0.15. In μ = 1 − E[Ω]Rn/(λA) the
  kernel wins, so priced risk relaxes the IC. That is the opposite of the appendix's narrative, and at 5 bp
  (slack at rest) it leaves the impact credit-spread channel shut.
- **Appendix B:**
  - The deterministic SS solves.
  - **The global solution does not.** The first stage (no default, π = 0, coarse) has no cold-start
    solution. Continuation reaches the household, markup, size, η, foreign-share and home-bias values, the last
    only at a precision floor of 6e-8. It then breaks a fifth of the way through the F coupon-decay move, above 1e-7,
    at a near-singular Jacobian. G, Bohn, banker β and the 5 bp wedge are never reached.
  - Still ahead of that stage: the 5 bp wedge (μ̄ = 0.001, the kink at rest), the genuine hard max, the default
    regime with an insolvent D bank, and the rest-point loop. Each was harder than d0 at the old calibration.
- **Failure points, in order of appearance:**
  1. warm time-iteration sweeps diverge where μ > f;
  2. the cold start lands on a non-root (LM stalls);
  3. a near-fold at the appendix's trade openness (ϖ_D 0.68 with η 1.5 and κ 11.3) and fiscal values (G, Bohn);
  4. (ahead) the kink at μ̄ = 0.001 under a genuine max;
  5. (ahead) the default regime;
  6. (ahead) TPI box coverage.

**Why the steps must be small** (`probe_svd.py`, at the stalled 1/8 home-bias rung): the column-scaled
collocation Jacobian has condition ≈ 2.5e5. Its weakest direction (σ/σ_max = 4e-6) carries 65% of the stalled
residual. That direction is one node, K_F+ (F capital at the upper edge), where α_D (D franchise value), K′_D
(D capital choice) and A_D (D household deposits) move together. It maps onto the C identities and UIP at the
centre. Removing the residual needs a large move along that direction: Newton's line search rejects it and LM
regularises it away. The coarse system is near a fold there, as the 2026-10-21 record found ("solvability flips
with ~1e-3 parameter changes"). The same-sized rung converges from a nearby solved point (1/16 → 3/16) where it
stalls from 0.

## Phase 3 — fix plan (strict escalation; each tier only after the one before fails)

### Tier 1 — computational

| # | Fix | Status / evidence | Trigger for the next tier |
|---|---|---|---|
| T1.1 | Genuine max via a final hard polish after the smoothing walk (IC, book pairs, floor) | built (`_soften`, `_hard_polish`) | the hard polish cannot root d0 at μ̄ = 0.001 |
| T1.2 | s box 3.65 sd, so the 2% headline shock lies inside | built | – |
| T1.3 | No warm time-iteration sweeps before the first Newton at this calibration (they diverge: the α recursion is expansive where μ > f) | measured: cuts the first-rung residual from 12.6 to 5.3e-3, not enough alone | – |
| T1.4 | Continuation path from a solved economy, SS-relative seeding, adaptive rungs | built (`diag_path.py`); reaches size, η, markup, households, foreign share and home bias, then breaks at 0.19 of the δ_b^F move (floor > 1e-7) | the precision floor rises toward the appendix (1e-10 → 6e-8 → 1.2e-7) |
| T1.5 | Better Newton accuracy near the fold: central-difference Jacobian columns, or an analytic Jacobian for the IC and α identities; row/column equilibration before the LU | central differences tested 2026-10-26: cond 1.5e10, no descent, so ruled out; the analytic Jacobian is untested | floor stays above 1e-8 |
| T1.6 | Box: narrow the edges where net worth swings ±13–16% (μ reaches 0.15 > f there), or the aligned two-box design | half bands tested 2026-10-26: worse (8.4e-5; the kink at μ_ss is sharper in unit coordinates); wider bands push μ further above f; the two-box design is untested at this calibration | – |
| T1.7 | Denser D-bank block (Smolyak level 2 on K_D, P_D, b_DD) so the kink at μ̄ = 0.001 is resolved | historically folds at leverage 8 (2026-10-22) | – |
| T1.8 | Integrate T1.4 into `solve_recursive` (a parameter walk with SS-centred boxes, replacing the raw-value `_seed_from` in `_walk_instruments`) | to build once d0 reaches the appendix | – |

### Tier 2 — calibration (small, economically defensible; each needs your decision)

| # | Change | Why | Evidence / magnitude |
|---|---|---|---|
| T2.1 | A recapitalisation floor in default, 15% → ~50% of n_ss | caps the default-branch franchise value α ≈ λθ/floor (6.8 → ~2), the leg that makes priced risk relax the IC | saved run: the default leg adds 0.127 of the +0.395 rise in E[Ω] at p^d 2%; the 2012–13 HFSF recapitalisation came with the PSI |
| T2.2 | Bohn 4γ_τ 0.04 → ~0.06 (the 2001–19 estimate in the appendix's own table) | moves the debt rule away from its explosive region (q* = 1 − γ/δ_b: 0.821 → 0.732) and the near-unit root (0.993 → 0.988) that breaks d0 alone | Bohn alone fails at 0.04, solves at 0.16 |
| T2.3 | A rest-point wedge above 5 bp if the IC cannot bind at rest | the guide's 2026-08-29 record: at ≤ 25 bp μ_rest = 0 and the economy rests on the kink (fitted vs exact gap 4×) | the rest-point loop will show it |
| T2.4 | Banker exit f 0.08 → ~0.10 | widens the region where α's recursion contracts (μ < f) | not tested |

### Tier 3 — theory (last resort)

| # | Change | Addresses |
|---|---|---|
| T3.1 | Recapitalisation as an explicit default-state policy rule rather than a numerical floor | the sign problem: risk relaxes the IC through the default-branch franchise value |
| T3.2 | Price the working-capital loan off the default-correlated kernel (a risky intra-period loan), so the risk channel reaches impact output even when μ falls | Y_D rises on impact when μ falls; the loan spread is the only impact channel |
| T3.3 | National deposit clearing as in 02-model, or rewrite the text for the union market | E3 |
| T3.4 | A KS-style aggregate rule for the HA block, or state the RA closure | E4 |

## Follow-up (2026-10-27/28) — Tier 2 run: 50% recapitalisation floor, Bohn 4γ_τ = 0.06

Your decision, run with every fix from this session. `run.NW_FLOOR` is 0.5 and `bohn_annual` is 0.06 (debt
root 0.987, explosive only below q = 0.732). Drivers: `diag_path.py` (d0 path), `m1_coarse.py` / `m1_kband.py`
(the risk solve). Logs: `runs/diag_path3..6.log`, `runs/m1_*.log`.

**No-default stage d0.**
- **Cold start:** still fails (first rung 1.4e-3 on `dep_clear`; it was 12.6 at Bohn 0.04).
- **Continuation**, from the household/markup/size/foreign-share checkpoint, reordered (fiscal first, home bias last):

  | Group | Result |
  |---|---|
  | δ_b^F | one rung, 5.7e-9 |
  | Bohn 0.275 → 0.06 | one rung, 3.7e-8 (loose) |
  | G | half floors at 2.7e-7 from two different starts; reached at 1.4e-6 (loose) |
  | banker β | one rung, 1.4e-6 (loose) |
  | wedge 100 → 5 bp | converged properly at 5.1e-9; the low wedge eases the fold |
  | home bias 0.85 → 0.68 | in eighths, endpoint 1.46e-6 (loose) |

  **The full appendix-B economy (Bohn 0.06) is reached at d0**, with rungs kept at ≤ 2e-6 (flagged). The repo's
  acceptance is 1e-9.
- **Diagnosis of the floors:** not FD accuracy, not the box size, not `kappa_nfa` (0.05 makes it worse). At the
  stalled points the solution's K_F ± nodes show a spurious symmetric fall in D bank net worth (0.93 and 0.84 of
  SS), which looks like a second branch: F is 11.3 × D, so a ±2% band in F capital is ~23% of D capital in union
  units.

**Risk solve with the 50% floor** (coarse, from the d0 endpoint):
- **Haircut ladder:** converges at every rung to the calibrated 45% recovery (≤ 1e-8).
- **Joint solve on the SS-centred V5 box:** stalls at 2.3e-5. The rest point is outside it: K_D +2.07% against
  the ±2% band, b_DD +9.1% (the debt drift of a near-unit-root Bohn rule).
- **Re-centring the box on the unconverged rest point diverges** (4.1e-4, then 4.2, with the rest point running
  away). Rest points read off unconverged rules are not reliable.
- **Wider capital band, SS-centred:**

  | | k_band 3% | k_band 4% |
  |---|---|---|
  | Joint residual (60 steps) | 7.8e-7 | 1.7e-5 |
  | Rest point inside the box | no, K_D at the edge | yes, b_DD at its edge |
  | μ_D / loan spread at rest | 0.0001 / 0.4 bp (slack) | 0.013 / 72 bp (binding) |
  | Leverage / sovereign spread at rest | 6.6 / 9 bp | 7.55 / 59 bp |

  Two qualitatively different rest points. The coarse solution is not unique (as recorded on 2026-10-22).

**The mechanism at the 4% box** (exact reads; fitted agree in sign except C_D):

| p^d | ΔQ_bD | Δn_D | Δn_F | ΔY_D | ΔC_D | Δ lending spread | μ_D | E[Ω] | α′ default |
|---|---|---|---|---|---|---|---|---|---|
| 0.1% (rest) | | | | | | | 0.013 | 0.921 | |
| 0.5% | −4.4% | −6.9% | −0.14% | −0.056% | +0.12% | +73 bp | 0.028 | 0.964 | 1.57 |
| 1% | −7.0% | −11.0% | −0.23% | −0.089% | +0.19% | +115 bp | 0.037 | 0.992 | 1.57 |
| 2% | −10.1% | −16.1% | −0.35% | −0.133% | +0.28% | +175 bp | 0.051 | 1.030 | 1.57 |

**The 50% floor fixes the sign.** The default-state franchise value falls from ~7 to ~1.57, so the default leg
no longer swamps the net-worth loss. The IC tightens with risk, the lending spread rises and output falls, as
appendix B describes. Every acceptance sign holds except C_D, which rises in the exact read (≈ 0 fitted). These
numbers come from an unconverged joint solve (1.7e-5), before the genuine-max polish, so they are indicative.

**Converging the 4% solution:**
- 60 more Newton steps on the same box: no progress (1.71e-5 again, an identical point). It is a hard stall.
- Widening the debt band to 18%: the solve runs off to the box corners (rest point at every edge, leverage 14,
  sovereign spread 1,000 bp, E[Ω] 31). Discarded.

### Where this leaves the global model (2026-10-28)

1. **Economics:** the Tier-2 pair (50% floor, Bohn 0.06) turns the risk channel the right way round: the IC
   tightens, the lending spread rises and output falls. The 15% floor made it relax.
2. **Numerics:** at appendix B the coarse collocation system is near-singular. The no-default stage reaches
   the appendix only by continuation, at residuals ~1e-6. The joint stage stalls at ~1e-5. Its rest point depends
   on the box (slack at k_band 3%, binding at 4%), so the solution is not unique.
3. **Not yet done:** the genuine-max polish, the smoothing walk, the rest-point calibration loop, the fine check.
   All need a converged joint stage first.
4. **Next escalation (Tier 1, method level):**
   - an analytic or automatic-differentiation Jacobian (FD at cond 1e10 cannot root it);
   - a grid that resolves the K_F dimension's effect on D (or drops K_F's quadratic term, since F is 11.3 × D);
   - the aligned two-box design, built on 2026-10-21 and switched off.
