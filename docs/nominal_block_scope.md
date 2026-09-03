# Scope: a nominal block for the monetary union

**Date:** 2026-08-29
**Status:** design, not implemented.
**Why now:** with `size_F/size_D = 8` the *union-wide* deposit-rate leak is closed
(rdep_F moves −1.4 bp/yr against −13.2 bp/yr at equal country size). What remains is
**D-specific and flexible-price**, and it is now the single largest gap to Bocola's
−0.157%.

---

## 1. The problem, measured

Under GHH the only channel from the financial block into output is the
working-capital rate

```
r_wc = rdep + λ_K·μ / E[Ω]
```

and the two legs move in opposite directions. At the s-refined solve, impact of a
p^d → 1.98%/qtr shock, read against the model's own rest point:

| leg | bp/yr |
|---|---|
| credit spread λμ/E[Ω] | **+42.9** |
| deposit rate rdep_D | **−39.4** |
| **net r_wc** | **+3.8** |

**92% of the credit-spread rise is cancelled before it reaches a firm's wage bill.**
Holding rdep fixed instead — which is exactly what Bocola's §V.C small open economy
does, since his `R = 1/β + 0.01·(B_for/gdp)` is a *world* rate — takes the output
response from −0.081% to roughly −0.20%, i.e. onto his −0.157%.

### Where the 39 bp comes from

Two components, and only the first was fixed by country size:

1. **Union-wide** (≈1.4 bp): banks delever, households want to save, the union real
   rate falls. This is a correct general-equilibrium response and Bocola's *closed*
   model has it too (his R falls ~16 bp/yr).
2. **D-specific** (≈38 bp): D's terms of trade `p` jump **+0.15% on impact and revert**.
   Real deposit-UIP (residual 6) prices that reversal as a low D real rate:

   ```
   (1 + rdep_D) = (1 + rdep_F)·E[p′]/p + κ_nfa·nfa/Y
   ```

   With `E[p′] < p`, rdep_D must sit below rdep_F. Under **flexible prices the terms of
   trade are a jump variable**, so the whole adjustment happens in one quarter and the
   implied real-rate differential is large.

### Why this is the wrong sign empirically

In 2011–12 periphery **bank funding costs rose** — deposit flight, closed wholesale
markets, TARGET2 balances. This model has D's real funding cost *falling* 39 bp exactly
when its sovereign is under stress. That is the counterfactual signature of a missing
nominal block, and `CLAUDE.md` has flagged it since the rework: *"real interest parity
currently plays the role of the single policy rate."*

---

## 2. What a nominal block changes

The union has **one** nominal policy rate. National *real* rates then differ only by
expected inflation differentials, and with sticky prices the terms of trade move
**slowly** instead of jumping — which is precisely the mechanism generating the 38 bp.

### 2.1 Price setting: Rotemberg, not Calvo

Recommend **Rotemberg** quadratic adjustment costs.

- The monopolistic structure is already in place: `epsilon_D = epsilon_F = 6` and
  `markup_ss = (ε−1)/ε`. Today `mc` is *fixed* at that value (`firms.solve_firm_path`
  hard-codes `mc = markup_ss(cal, country)`). Making `mc` endogenous is the change.
- Rotemberg has **no price-dispersion state**. Calvo would add one per country, i.e.
  +2 states on a grid whose cost already scales as n². At a first-order-equivalent
  calibration the two are observationally close; the state saving is decisive here.

New equation per country (producer-price inflation π):

```
π_t (1 + π_t) = (ε/φ_p)·(mc_t − (ε−1)/ε) + β·E_t[ Λ' · π_{t+1}(1 + π_{t+1})·(Y_{t+1}/Y_t) ]
```

with `φ_p` the Rotemberg cost parameter, calibrated to a target slope (equivalently a
Calvo duration of ~4 quarters).

### 2.2 Monetary rule

```
1 + i_t = (1 + i*)·(1 + π^union_t)^{φ_π} · (Y^union_t / Y^union)^{φ_y}
```

with union inflation the **mass- and price-weighted** average — `size_ratio` already
exists in `trade.py` for exactly this kind of aggregation. `φ_π = 1.5`, `φ_y = 0.125/4`
is the standard starting point. `φ_π → ∞` (a strict inflation target) is the clean
limiting case worth reporting: it pins π^union = 0 and makes the union real rate
constant, which is the cleanest test of the mechanism.

### 2.3 Deposits become nominal

This is the substantive change to the financial block. Today deposits are *own-good*
real claims at *national* real rates, tied by real UIP. Under the union they are
**nominal** claims at the **common** rate `i_t`, and the realised real return differs
across countries by realised CPI inflation:

```
1 + rdep_D,t+1 = (1 + i_t) / (1 + π^CPI_D,t+1)
1 + rdep_F,t+1 = (1 + i_t) / (1 + π^CPI_F,t+1)
```

Consequences for `point_map.py`:

- `rdep_D`, `rdep_F` **leave the unknown vector** (−2). They become functions of `i`
  and next-period inflation, so they enter the Eulers *inside the expectation* rather
  than as deterministic returns — this is a real retiming, not a substitution.
- **Residual 6 (real deposit-UIP) is deleted** (−1). It becomes an identity: one
  nominal claim, two realised real returns.
- The bank's deposit obligation `P' = (1+rdep)·dep − (1+r_wc)·L_wc` is now a **nominal**
  obligation deflated by realised inflation. The working-capital loan is intra-period
  and unaffected.
- `r_wc = rdep + λμ/E[Ω]` becomes `r_wc = i + λμ/E[Ω]` in nominal terms, deflated in the
  labour FOC. **This is the whole point:** the firm's financing cost is now anchored to
  the ECB rate rather than to D's own real rate.

### 2.4 The terms of trade become a state

Under sticky prices `p` (D-goods per F-good, in producer prices) is no longer a jump
variable:

```
p_{t+1} = p_t · (1 + π_F,t+1) / (1 + π_D,t+1)
```

- `p` **leaves the unknown vector and joins the state vector** (−1 unknown, +1 state).
- **Residual 7 (goods_D) no longer pins `p`.** Under sticky prices output is
  demand-determined at the posted price, so goods-market clearing pins *quantities*.
  **This is the one place the mapping is not mechanical and needs deriving before any
  code is written** — the labour FOC (residuals 3–4) and goods clearing (residual 7)
  have to be re-sorted into (labour supply, labour demand, market clearing) with `mc`
  endogenous.

### 2.5 Net accounting

| | now | with the nominal block |
|---|---|---|
| states | 10 | **11** (+p) |
| solved unknowns / point / regime | 13 | **12** (−rdep_D, −rdep_F, −p, +π_D, +π_F) |
| stored rules (collocation unknowns) | 19 | ~18–20 |
| coarse grid points (μ=1) | 21 | 23 |
| s-refined points (m=5) | 95 | 115 |
| dense Jacobian cost | — | **≈ +18%** (scales as n²) |

Nothing in `collocation.py`, `state_grid.py` or the solve ladder changes. The bank,
household, government and trade blocks are untouched except for the nominal retiming
of the deposit contract.

---

## 3. Why the collocation rework makes this affordable

The Phillips curve introduces a **new forward-looking recursion** — π depends on E[π′] —
which is a slow mode of exactly the kind that made time iteration unusable (the
franchise-value recursion contracts at 0.990/sweep, 235 sweeps per decade). A damped
fixed-point iteration would inherit a second such mode and compound the problem.

The global Newton has no contraction rate to leak through: it drives the residual to
its arithmetic floor regardless of how many slow forward-looking blocks the system
contains. **This is the change that makes the nominal block practical**, and it is worth
saying so explicitly when the two pieces of work are written up together.

---

## 4. Expected payoff, and how to falsify it

Predicted, from the measured decomposition: the D-specific 38 bp offset largely
disappears, the net wedge goes from +3.8 bp/yr to something near the full +42.9 bp, and
the impact output response moves from **−0.081% to roughly −0.15%/−0.20%** — at or just
past Bocola's open-economy −0.157%.

### The falsification test — RUN, and it passes decisively

Implemented as `cal["union_nominal_rate"] = True` in `point_map.py` (residual 6 becomes
a literal `rdep_D = rdep_F`). **This is not an equilibrium** — with own-good deposit
legs the real-exchange-rate valuation profit is unassigned — and the solve does not
reach the acceptance floor (max|F| = 9.6e−3). Diagnostic only. Measured on the coarse
grid at the *8 bp* calibration:

| | real UIP | `rdep_D = rdep_F` |
|---|---|---|
| deposit-rate offset | −45.0 bp/yr | **−7.3 bp/yr** |
| credit spread | +54.4 bp/yr | **+187.5 bp/yr** |
| net wedge r_wc | +9.9 bp/yr | **+185.7 bp/yr** |
| Y_D fitted / exact | −0.094% / −0.007% | **−1.51% / −1.57%** |

Three things follow. **(1)** Pinning the rate removes essentially the whole offset, so
the terms-of-trade/UIP channel *is* the dominant remaining gap — the block is worth
building. **(2)** It also cures the KKT-kink identification problem as a side effect
(fitted and exact agree to 3.5% instead of differing by 10×), because μ moves far off
the kink. **(3)** The crude version **overshoots badly**: −1.5% is ten times Bocola's
−0.157%, because with the rates literally equal and no inflation to share the
adjustment, the constraint absorbs all of it. So −1.5% is an **upper bound** and
−0.11% (the current calibration) a lower one; his −0.157% sits comfortably inside, and
a proper sticky-price block should land near it rather than at either end.

---

## 5. Sequencing

1. The one-line falsification test in §4. One coarse solve.
2. Derive the §2.4 re-sorting of labour supply / labour demand / goods clearing under
   endogenous `mc`. Paper, not code.
3. `firms.py`: endogenous `mc`, Rotemberg PC, per-country π.
4. `point_map.py`: nominal deposit contract, `i` from the rule, delete residual 6,
   add two PC residuals, move `p` to the state.
5. `state_grid.py`: 11-state box (mechanical — `refine=` and the box builder already
   take the dimension from `STATE_NAMES`).
6. `steady_state.py`: zero-inflation SS. Every existing SS object is unchanged at
   π = 0; the new content is `φ_p` and the rule's coefficients, neither of which binds
   at the SS.
7. Re-run the ladder; check `test_recursive_nesting` N1/N2 still hold at π = 0 — the
   nominal block **must** nest the current model exactly when prices are flexible
   (`φ_p → 0`) and the rule is a real-rate peg.

**Nesting is the acceptance test.** `φ_p = 0` has to reproduce today's solution to the
solver's floor, the same way `size_F = size_D` reproduces the symmetric model and
`π = 0` reproduces the risk-free one.

---

## 6. Known open questions

- §2.4's re-sorting of residuals 3, 4 and 7 is the one genuine derivation.
- Whether the **household** deposit Euler needs the inflation risk priced separately, or
  whether the existing GHH kernel handles it once the return is inside the expectation.
- Whether `κ_nfa` (the SGU stationarity premium) is still needed once nominal rates are
  common, or whether it double-counts. It is currently load-bearing for stationarity —
  the rest point is a unique global attractor with it on.
- Bocola himself has **no** nominal block; his §V.C fix is to make the rate exogenous.
  Making D small enough that the union rate is effectively exogenous to it is the
  cheaper approximation, and §4's test is also a test of *that* route.
