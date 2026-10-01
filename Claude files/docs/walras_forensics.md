# HHBANK v12 — Walras Law Forensics

**Status:** analytical derivation complete; numerical confirmation pending (run_audit.py in progress).

## 1. Setup

Transition system: 23 unknowns / 23 targets. Targeted: all D/F internal market clearing + FOCs + `goods_mkt_D`. **Untargeted residuals computed by the model: `ca_res_D`, `goods_mkt_F`, `global_goods_res`.**

- Theoretically redundant equations: in a 2-country economy with two goods markets, two deposit markets, two bond markets and an external account, summing **all** agents' budget constraints country-by-country yields exactly **two** aggregate identities; given all targeted equations, `ca_res_D` and `goods_mkt_F` should hold automatically.
- Implemented redundant equations: the authors drop `ca_res_D` and `goods_mkt_F` (and `global_goods_res`, which is the linear combination `goods_mkt_D + p·goods_mkt_F`, not independent).
- Actually redundant equations: **derived below — they are NOT redundant in this implementation.** The implementation leaks.

## 2. Derivation of country-D aggregate identity (from the implemented equations, not from theory)

All in D-goods. P ≡ P_CES_D. Household aggregation of the EGM budget (exact, from `hh_D` code):

```
P·C + Dep_t = (1+rdep)·Dep_{t−1} + w·N + div − P·TAX          (HH)    [Dep_t ≡ P_t·DEP_t]
```

Bank definitions (`bank_return_D`, `intermediation_P2_D`, `banker_div_res_D`, `intermediation_P3_D`):

```
(1+rn_t)n_{t−1} = (1+rk_t)Q(−1)K(−1) + (1+rb_D,t)q_bD(−1)b_DD(−1) + (1+rb_F,t)q_bF(−1)b_FD(−1) − (1+rdep_t)Dep_{t−1}
n_t  = (1−f)[(1+rn)n_{t−1} + cap_profit] + m
div  = f[(1+rn)n_{t−1} + cap_profit] − m − Phi − T
⇒ n_t + div = (1+rn)n_{t−1} + cap_profit − Phi − T
Q·K_t + q_bD·b_DD,t + q_bF·b_FD,t = n_t + Dep_t                (balance sheet / P3 + deposit_mkt)
```

Substituting HH + bank + `capital_producer_profit` ((1+rk)Q(−1)K(−1) + cap_profit = mpk·K_{t−1} + Q_t·K_t − I_t):

```
P·C + I + Phi + T + q_bD·b_DD,t + q_bF·b_FD,t
   = w·N + mpk_t·K_{t−1} − P·TAX
   + (1+rb_D)q_bD(−1)b_DD(−1) + (1+rb_F)q_bF(−1)b_FD(−1)
```

Government budget (`budget_residual_D` = 0, targeted) + bond clearing (b_DD = b_gov − b_DF) + the bond cash-flow identity `(1+rb_D)q_bD(−1)b = [δ_b(1−def·h·mult) + (1−δ_b)q_bD·surv]·b` collapse the domestic-bond terms into cross-border terms only:

```
P·C + I + G + Phi + T + nfa_D,t − (1+rb_F)q_bF(−1)b_FD(−1) + (1+rb_D)q_bD(−1)b_DF(−1)
   = w·N + mpk_t·K_{t−1}
```

Comparing with the definitions of `goods_mkt_D` and `ca_res_D`:

```
┌─────────────────────────────────────────────────────────────────────┐
│  goods_mkt_D + ca_res_D  =  Y_t − w_t·N_t − mpk_t·K_{t−1}           │  (★)
└─────────────────────────────────────────────────────────────────────┘
```

**The aggregate identity does NOT say `goods_mkt_D + ca_res_D = 0`.** It equals the factor-income gap. With the implemented production function Y = Z·K_t^α·N^{1−α} and w = (1−α)Y/N, mpk = αY/K_t:

```
Y − wN − mpk·K_{t−1} = αY·(K_t − K_{t−1})/K_t   ≠ 0 off steady state (first order in dK)
```

### Consequence
`goods_mkt_D` is targeted to 0 ⇒ **ca_res_D,t = αY·(K_t − K_{t−1})/K_t ≈ (αY/K)·(dK_t − dK_{t−1})** along any transition with investment dynamics. The current account identity CA = ΔNFA **fails at first order**, and the leak is *not* in the external block — `external_account_D` is correct — it is caused by the **capital-timing mismatch** between `labor_D`/`capital_adj_D` (production uses K_t) and `bank_return_D` (banks are paid mpk_t on K_{t−1}).

The income on newly installed capital, mpk_t·(K_t − K_{t−1}), is paid to **no agent**. It is produced (it's in Y) and absorbed (goods market is cleared by construction) but never enters any budget constraint.

### Two fixes (either restores Walras)
1. Standard SSJ timing: `Y_t = Z·K_{t−1}^α·N^{1−α}`, `mpk_t = α·Y_t/K_{t−1}` in both `labor_D` and `capital_adj_D` (and F analogues). Then Y = wN + mpk·K(−1) exactly.
2. Keep Y(K_t) but pay the capital producer the within-period product of new capital: cap_profit += mpk_t·(K_t − K_{t−1}).

**Adopted: option 2 (author decision).** Commit 4c810e1 on `audit`. Verified numerically equivalent on the Walras criterion (ca_res_D ≤ 5.8e−8, goods_mkt_F ≤ 8e−10, stationary, correct default-shock signs).

## 3. Country-F aggregate identity — second, independent leak

Same derivation in F-goods, with the extra complication that F-bank bond positions are D-good claims. The **true** F-good payoff on bonds held from t−1 is `(1+rb,t)·q(−1)·b(−1)/p_t`. The implemented `bank_return_F` pays `(1+rb,t)·q(−1)·b(−1)/p_{t−1}` (shares φ converted at p(−1), returns left in D-good terms). Hence F-bank's measured income misses the terms-of-trade revaluation:

```
goods_mkt_F + [F external-flow residual] = (Y_F − w_F N_F − mpk_F K_F,{t−1})        (factor gap, F)
                                          + Payoff_t^{D-goods}·(1/p_{t−1} − 1/p_t)  (ToT revaluation gap)
```

where Payoff^{D-goods} = (1+rb_F)q_bF(−1)b_FF(−1) + (1+rb_D)q_bD(−1)b_DF(−1).

Since `ca_res_D` is already forced non-zero by (★) and there is no F external residual in the model, **`goods_mkt_F` absorbs: its own factor gap + the ToT revaluation gap + p·(D-leaks) via the global identity**. Numerical decomposition in `run1_results.json` (`decomp_*` keys) tests:
- `gapD` ≈ `ca_res_D` path (prediction 1)
- `goods_mkt_F` ≈ `gapF + rev_gapF` (prediction 2)

## 4. Hypothesis status (pre-numerics)

| Hypothesis | Status | Evidence |
|---|---|---|
| H1: Walras leak in external block | **REVISED** | External block itself is internally correct; leak originates in production/capital timing (★) and F-bank return units. ca_res_D is the *symptom*. |
| H2: D/F conversion inconsistency | **CONFIRMED (analytical)** | `bank_return_F` and `divert_bond_foc_F` use D-good returns vs F-good rdep without p-conversion; `divert_portfolio_adj` does convert. |
| H3: sovereign return timing | REJECTED so far | bond_return/budget/default timing internally consistent (def_t from b_gov_{t−1}, priced via def(+1)). |
| H8: convergence masks structural errors | **CONFIRMED (analytical)** | 23×23 system is solvable regardless of the leak; residual-free targets coexist with violated untargeted identities. |

## 5. What evidence would prove this wrong?

- If `ca_res_D` and `goods_mkt_F` IRF paths are numerically ~0 (≤1e-10) at first order, derivation (★) is wrong somewhere → re-derive.
- If `ca_res_D` ≠ gapD, some other budget element (m, Phi, T, dividend flow) was mis-transcribed → re-audit P2/div.
- Numerical verdict: see §6 below once run completes.

## 6. Numerical verdict (run_audit.py, T=500, baseline notebook pipeline replicated exactly)

Targeted residuals at machine zero as expected: `goods_mkt_D` ≤ 4e−16, `deposit_mkt_D/F` ≤ 2e−15.

**Untargeted residuals leak at first order, exactly as derived:**

| Shock (1%, ρ=0.8) | max\|ca_res_D\| | max\|goods_mkt_F\| |
|---|---|---|
| TFP Z_D | 1.49e−4 | **2.02e−2** |
| Default shock_def_D | 6.25e−5 | 1.26e−3 |

Attribution (predictions from §2–§3 vs computed paths):

| Test | Result |
|---|---|
| `ca_res_D` vs D factor-timing gap `dY − w·dN − N·dw − mpk·dK(−1) − K·dmpk` | max diff 2.3e−7 (0.16% of leak), corr 0.99999983 → **W-1 proven: leak ≡ capital-timing gap** |
| `goods_mkt_F` vs F factor gap + ToT revaluation gap | corr 0.9998; residual after subtraction = 1.49e−4 = exactly the D-gap re-entering through the global identity → **W-2 proven and dominant** |
| Relative size | ToT revaluation gap (2.1e−2) ≈ **140×** the factor-timing gap → W-2 is the first-order problem for cross-country results |

Scale interpretation: on a standard 1% TFP shock, country F's goods market fails to clear by **2% of F GDP** on impact. Every F-country IRF (Y_F, C_F, p, NX) and every spillover/contagion statement from this model is contaminated at first order.

### Answers to the three Walras questions
- **Theoretically redundant:** `ca_res_D` and `goods_mkt_F` (given the 21 other targets + goods_mkt_D).
- **Implemented as redundant:** the same two (dropped from targets) — the *choice* of dropped equations is correct.
- **Actually redundant:** **neither.** Both fail off-SS because two budget-constraint elements (capital income timing; F-bank bond revaluation) are inconsistent with the rest of the system.

## 7. Fix verification (fix_test.py — patches applied to block *copies*, model files untouched)

Patches: W-1 (Y and mpk use K(−1), both countries), W-2 (`bank_return_F` converts bond returns to F-goods with `p(−1)/p`), optionally W-3 (`divert_bond_foc_F` converts with `p/p(+1)`). Same SS (unchanged by construction), same 23×23 system.

| max abs residual | baseline | fix12 (W-1+W-2) | fix123 |
|---|---|---|---|
| ca_res_D (TFP) | 1.49e−4 | **4.1e−7** | 3.4e−7 |
| goods_mkt_F (TFP) | 2.02e−2 | **4.2e−8** | 4.5e−8 |
| ca_res_D (def) | 6.25e−5 | 2.9e−8 | 2.4e−8 |
| goods_mkt_F (def) | 1.26e−3 | 5.2e−9 | 5.0e−9 |

Leaks collapse by 3–6 orders of magnitude, to the floor set by the SS Broyden tolerance (~4e−7). **Root-cause attribution is proven by repair.** W-3 does not affect Walras (it is an optimality condition, not a budget element), as predicted.

### Contamination of headline results (baseline vs fix123, peak responses)

| Variable, 1% TFP-D shock | baseline | fixed | error |
|---|---|---|---|
| C_F | 0.0102 | 0.0063 | **+63% overstated** |
| U_F (welfare) | 0.0139 | 0.0076 | **+81% overstated** |
| n_inter_F | −0.041 | −0.188 | **understated 4.6×** |
| q_b_F | 0.0394 | 0.0367 | +7% |
| p (ToT) | 0.0197 | 0.0264 | −25% understated |
| Y_F | 0.00152 | 0.00110 | +39% overstated |
| Y_D | 0.0288 | 0.0283 | +2% (D-side barely touched) |

| Variable, 1pp default-D shock | baseline | fixed | error |
|---|---|---|---|
| Y_D | +1.47e−4 | +8.8e−5 | +67% (still **positive** — perverse sign not caused by leaks; see audit.md T-2) |
| C_F | 5.7e−4 | 3.1e−4 | +87% |
| U_F | 8.2e−4 | 3.7e−4 | +124% |
| n_inter_F | −0.045 | −0.034 | +34% |
| spread_rb | 0.00260 | 0.00244 | +7% (D-side spread robust) |

Conclusion: **all cross-country spillover and welfare-spillover results are first-order wrong in the baseline**; domestic D-side responses (Y_D, spread) are only mildly affected. Any contagion/spillover claims must be re-run after the fix.

## 8. TPI extension — third, *injected* leak (no CB budget constraint)

`domestic_bond_clearing_tpi` sets `b_D_D = b_gov_D − b_D_F − cb_buy_D` with no block recording who pays `q_D·cb_buy` or receives the CB's coupon income. Extending §2's algebra, the D identity becomes:

```
goods_mkt_D + ca_res_D = factor_gap_D + [ q_D·cb_t − (1+rb_D,t)·q_D(−1)·cb_{t−1} ]
```

i.e. the private sector receives/loses unbacked resources whenever the CB position changes or earns returns. The TPI welfare figures (Figures TPI-6..9) are computed on top of this flow.

**Numerical verdict (tpi_test.py):** at γ=10, max |ca_res_D| = 2.6e−2 (2.6% of quarterly GDP per period), 417× the no-policy leak; corr with the analytical hole formula = −0.99999 (sign = bookkeeping convention). Discounted over 100q, the unbacked flow is ≈ 40% of the measured welfare gain ΔW_D (+0.103% of SS consumption). TPI welfare conclusions are not interpretable until a CB budget constraint is added (remit purchases and coupon income through `budget_residual_D`).

**FIXED (2026-06-11):** `budget_residual_D_tpi` (notebook cells TPI-1/TPI-2) consolidates the CB with the government: remittance `rem_cb_D = coupon-on-cb(−1) + q_b·surv·(1−δ_b)·cb(−1) − q_b·cb`, which exactly cancels the hole term above. Verified post-fix with all other fixes applied (tpi_log2.txt): max |ca_res_D| = 5.3e−8 **at every γ ∈ {0,2,5,10}**, max |goods_mkt_F| ≤ 9e−10. With accounting closed and the doom loop live (post T-2), the TPI results change qualitatively: the default shock now costs D households W_D = −4.30 (was −0.86), TPI at γ=10 recovers ΔW_D = +1.88 while ΔW_F = −1.90 — i.e. **TPI is approximately a zero-sum transfer of the crisis burden from D to F**, and the spread itself is *not* compressed (peak rises slightly with γ; def_rate is debt-driven). These regenerated numbers supersede every TPI figure in the notebook.

## 9. Final hypothesis ledger

| Question | Answer |
|---|---|
| Theoretically redundant eqs | ca_res_D, goods_mkt_F |
| Implemented as redundant | same two (correct choice) |
| Actually redundant | neither, pre-fix; **both, post-fix** (proven: residuals → SS-tolerance floor) |
| Leak sources | W-1 capital timing (D & F), W-2 F-bank ToT revaluation; TPI adds missing CB budget |
| Leak NOT caused by | external account block, trade identities, sovereign cash flows, household aggregation, deposit clearing — all verified exact |

