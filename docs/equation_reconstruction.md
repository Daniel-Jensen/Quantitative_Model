# HHBANK v12 — Mathematical Reconstruction (from code, not docs)

**Status:** Phase 1 complete (reconstructed from `equations_D.py`, `equations_F.py`, `equations_global.py`, `model_v12.ipynb` cells 2–21). Discrepancies flagged inline with ⚠.

Notation: country D is numeraire (P_D-good ≡ 1). `p` = price of F-good in D-goods (terms of trade). All F-country real quantities are in F-goods unless noted. Time-`t` variables unsubscripted; `(-1)`/`(+1)` denote lags/leads as in SSJ.

---

## 0. Key denomination conventions (reconstructed, nowhere documented in one place)

| Object | Units | Source |
|---|---|---|
| `Y_D, I_D, G_D, w_D, K_D (×Q_D), n_inter_D, D_supply_D, div_D, Phi_D, T_D` | D-goods | `market_clearing_D` comment |
| `C_D, DEP_D, TAX_D, X_D, z_D` | D-bundle | `income_D` divides by `P_CES_D` |
| `Y_F, I_F, G_F, w_F, n_inter_F, D_supply_F, div_F, Phi_F, T_F` | F-goods | `market_clearing_F` |
| `C_F, DEP_F, TAX_F, X_F` | F-bundle | `income_F` |
| **Both** sovereign bonds `b_gov_D`, `b_gov_F` (face), prices `q_b_D`, `q_b_F` | **D-goods** | `smart_steady_F` comment: "Bonds are D-good (numeraire) claims"; `budget_residual_F` divides bond flows by `p` |
| `rb_actual_D, rb_actual_F` | D-good returns | `bond_return_*` |
| `rdep_D` | D-good rate | `deposit_return_D` |
| `rdep_F` | **F-good** rate | `deposit_return_F`: household real return = `(1+rdep_F)·P_CES_F(-1)/P_CES_F` ⇒ deposit is an F-good claim |
| `rk_F, rn_F` | F-good returns | `capital_adj_F`, `bank_return_F` |

⚠ Consequence: any equation comparing `rb_actual_F` (D-good return) to `rdep_F` (F-good return) without a `p`-conversion mixes units. See findings W-2/W-3 in `audit.md`.

Both sovereigns issue D-good ("union currency") denominated debt. This is a monetary-union design choice, consistent in the government budget blocks, but **not consistently carried through F-bank return accounting**.

---

## 1. Households (per country, shown for D)

Heterogeneous agents, GHH preferences, single asset (bank deposit).

- Idiosyncratic productivity: 15-state Rouwenhorst, ρ=0.9, σ=0.3. (`nZ=19` branch loading `Discretisation/Outputs/Px_GMAR.txt` is dead in current calibration.)
- Disutility: v(N) = φ·N^{1+1/ν}/(1+1/ν), ν = frisch = 0.5. N is **uniform across households** (no idiosyncratic labor choice; N_D is an aggregate unknown).
- GHH composite: x = c − v(N). Utility u(x) = x^{1−1/σ}/(1−1/σ), σ = eis = 0.5.

Budget (bundle units), from EGM code (`hh_D`):
```
c_t + dep_t = Rgross_t·dep_{t-1} + z_t
Rgross_t = (1+rdep_t)·P_CES(-1)/P_CES        (deposit is D-good claim)
z_t   = λ·(y_pre)^{1−τ} − T_ls,   y_pre = (w·N·e + div)/P_CES
t_paid = y_pre − z                            (HSV progressive tax + lump-sum Bohn component)
```
Euler (EGM): u'(x_t) = β·E[Rgross_{t+1}·u'(x_{t+1})], borrowing limit dep ≥ 0 (grid min).
Aggregates: `C_D = E[c]`, `DEP_D = E[dep']`, `TAX_D = E[t_paid]`, `UCE_D = E[x^{-1/σ}]` (computed but unused downstream ⚠ dead output).

Labor supply: intratemporal GHH FOC imposed at aggregate level (`labor_market_D`):
`w/P_CES = φ·N^{1/ν}` (UCE cancels under GHH). `vphi` calibrated from SS via `labor_ss_D` with μ_w = 1.

⚠ Dividends `div` enter `y_pre` **uniformly** (not proportional to e), then are taxed progressively.
⚠ Welfare object: `U = X/C_ss` — a *composite-consumption* index, not utility u(x); valid only as a linear consumption-equivalent proxy.

## 2. Production & capital

- `labor_D`: **Y_t = Z_t·K_t^α·N_t^{1−α}** — uses *current* K_t.
- `labor_demand_D`: w = (1−α)Y/N.
- `capital_adj_D`:
  - ι_t = I_t/K_{t−1}
  - **mpk_t = α·Z_t·K_t^{α−1}·N^{1−α}** (current K)
  - **rk_t = (mpk_t + (1−δ)Q_t)/Q_{t−1} − 1** ⚠ pays current-K marginal product on capital bought at t−1
  - Q_t = 1/(γ0(1−ξ)ι^{−ξ}) (Hayashi-style, ξ=0.5)
  - K_t = (1−δ)K_{t−1} + (γ0·ι^{1−ξ} + γ1)·K_{t−1} (investment installed and **productive within the same period**)
- `capital_producer_profit_D`: cap_profit = Q(K − (1−δ)K(−1)) − I → accrues to bankers via gross_income.

⚠ **Timing clash (Finding T-1):** with Y_t = F(K_t, N_t), factor payments are w·N + mpk·K_t = Y_t. But capital owners (banks) are paid mpk_t per unit of K_{t−1} (via rk_t), and the capital producer receives only Q·ΔK − I. Income paid out = wN + mpk·K_{t−1} + cap_profit ⇒ a first-order income gap αY·(K_t−K_{t−1})/K_t exists off-SS. See `walras_forensics.md` for the exact identity.

## 3. Banking sector (Gertler–Karadi multi-asset, per country)

Balance sheet (D-bank, D-goods): Q·K + q_bD·b_DD + q_bF·b_FD = θ·n = n + D_supply.
F-bank (F-goods): Q_F·K_F + (q_bF·b_FF + q_bD·b_DF)/p = θ_F·n_F.

Portfolio return (`bank_return_D`):
```
rn_t = κ_{t−1}(rk_t − rdep_t) + φ_bD,t−1(rb_D,t − rdep_t) + φ_bF,t−1(rb_F,t − rdep_t) + rdep_t
κ_{t−1} = θ_{t−1} − φ_bD,t−1 − φ_bF,t−1,  φ_i,t−1 = q_i(−1)b_i(−1)/n(−1)
```
⚠ F-bank version converts shares with `p(−1)` but applies **D-good returns rb_actual without p(−1)/p_t revaluation** (Finding W-2).

Net-worth evolution (`intermediation_P2`): n_t = (1−f)·[(1+rn_t)n_{t−1} + cap_profit_t] + m, with m a **constant** (SS-computed). Exit prob f = 0.12.

Dividends (`banker_div_res`): div = f·[(1+rn)n(−1) + cap_profit] − m − Phi − T.
⚠ SS block `smart_steady_*` sets m = n(1−(1−f)(1+rn)) + Phi + T; dynamic P2 requires m = n(1−(1−f)(1+rn)) (no Phi+T). Consistent **only because Phi=T=0 in baseline** (Finding A-2, latent).

Franchise-value Bellman (`intermediation_P1`), per unit net worth:
```
ν_i,t = SDF_banker,t · Ω_{t+1} · (r_i,t+1 − rdep_{t+1}),   i ∈ {K, bD, bF}
η_t  = SDF_banker,t · Ω_{t+1} · (1 + rdep_{t+1})
Ω_{t+1} = f + (1−f)·λ_gk·θ_{t+1}
SDF_banker,t = β_inter·(X_{t+1}/X_t)^{−1/σ}   (aggregate GHH composite of households)
```
Incentive constraint (`intermediation_IC`), multi-asset divertability:
```
value = ν_K·κ + ν_bD·φ_bD + ν_bF·φ_bF + η
θ = value/λ_gk + (1−Δ_bD^eff)·φ_bD + (1−Δ_bF^eff)·φ_bF
Δ_i^eff = Δ_i + ψ_λB·def_rate_{i,t+1}    (default risk worsens bond collateral, ψ_λB = 3)
```
θ, ν, η solved each period by an inner `.solved()` block (5 unknowns / 5 targets).
λ_gk is a **constant** in transition (SS value from `steady_auxilliary`: λ = f/(θ(1/(β_inter(1+rn)) − (1−f))), the **single-asset** formula). Δ's are then **back-solved** in notebook cell 12 to make the multi-asset IC bind at SS. ⚠ This forces Δ_cross = 2·Δ_own with Δ_cross possibly > 1 (more-than-fully divertable — economically dubious; check calibration_review.md).

Bond FOCs (replacing market clearing for portfolio composition):
- Own bonds (`divert_bond_foc_D`): E_t[rb_D,t+1] − rdep_{t+1} = xs_ss + ψ_spread·def_rate_{t+1} + ψ_bD(φ_bD − φ_ss) + τ_mp, where ψ_spread = λ_gk·ψ_λB/(β_inter·Ω) (anchored post-SS). ψ_bD_D = 0 ⇒ pure pricing equation for q_b_D.
- ⚠ F-bank own-bond FOC (`divert_bond_foc_F`) compares D-good return `rb_actual_F(+1)` with F-good `rdep_F(+1)` with **no p-conversion**, although the dead block `domestic_bond_foc_F` did convert (`(1+rb)·p/p(+1)−1`) and the cross FOCs in `divert_portfolio_adj` do convert (Finding W-3).
- Cross bonds (`divert_portfolio_adj`): converted expected returns, level-penalty ψ·(b − b_ss), ψ_bF_D = ψ_bD_F = 0.5 on face-value levels.

## 4. Sovereign sector (per country, D shown; F divides flows by p)

Default rate (`government_default`):
```
def_rate_t = shock_def_t + def_scale·[(b_gov_{t−1}/Y_ss + c0)^{0.5} − (b_ss/Y_ss + c0)^{0.5}],  c0 = 0.05, def_scale = 0.25
```
Bond cash flows (Woodford geometric perpetuity, decay δ_b = 0.1; coupon = δ_b per face unit; (1−δ_b) of face survives):
```
haircut = 1 − recovery_rate = 1;  multiplier = writeoff_enabled = 0  ← BASELINE
rb_actual_t = [δ_b(1 − def_t·h·mult) + (1−δ_b)q_b,t(1 − ζ·def_t·h·mult)]/q_b,t−1 − 1
```
⚠ **With writeoff_enabled = 0 (baseline), default NEVER touches cash flows**: no haircut on coupons, no face write-off, recovery_rate irrelevant. "Default" is purely a *risk-premium/divertability* shock transmitted through ψ_spread·def_rate(+1) in the bond FOCs and Δ^eff in the IC. There are no realized default losses anywhere in the baseline model (Finding S-1, interpretation-critical).

Government budget (`budget_residual_D`):
```
coupon_t = δ_b(1 − def_t·h·mult)·b_gov_{t−1}
net_issuance_t = q_b,t·(b_gov_t − surv·(1−δ_b)·b_gov_{t−1}),  surv = 1 − ζ·def·h·mult
0 = coupon_t + G_t − P_CES·TAX_t − net_issuance_t      (G fixed at SS level; b_gov is the unknown)
```
Fiscal rule (`tax_rule`): T_ls,t = φ_λ·(b_gov_{t−1} − b_ss), φ_λ = 0.02, entering household income via z.
SS closure (`government_ss`): G = P_CES·TAX + net_iss − coupon (G residual at SS).

Consistency: SS and dynamic government blocks agree at SS. Bond-return and budget blocks use the same haircut conventions (no double-count of losses *in the budget*; with mult=0 nothing to double-count — must re-verify if writeoff_enabled is ever set to 1).

## 5. Open economy

- CES bundles, home bias ω = 0.85, trade elasticity ε = 1.5:
  - `P_CES_D = (ω + (1−ω)p^{1−ε})^{1/(1−ε)}`, `P_CES_F = (ω + (1−ω)(1/p)^{1−ε})^{1/(1−ε)}` ✓ mutually consistent.
  - Import demand: `IM_D = (1−ω)(P_CES_D/p)^ε C_D` (quantity of F-goods), `IM_F = (1−ω)(P_CES_F·p)^ε C_F` (quantity of D-goods) ✓.
- Trade balance: `NX_D = IM_F − p·IM_D` (D-goods), `NX_F = IM_D − IM_F/p` (F-goods). Identity NX_D + p·NX_F ≡ 0 ✓.
- External account (D, all D-goods):
```
nfa_D,t = q_bF,t·b_FD,t − q_bD,t·b_DF,t                  (market value, end of period)
ca_res_D = NX_D + (1+rb_F,t)q_bF(−1)b_FD(−1) − (1+rb_D,t)q_bD(−1)b_DF(−1) − nfa_D,t
```
  This is the correct valuation-consistent external flow identity (gross payoff incl. market value of surviving bonds finances the new position + NX). Targeted **only at SS**; in transition it must hold via Walras — see `walras_forensics.md`.
- Bond market clearing (`domestic_bond_clearing`): b_DD = b_gov_D − b_DF, b_FF = b_gov_F − b_FD (domestic bank is residual holder).
- Goods market (D): `goods_mkt_D = Y_D − (P_CES_D·C_D + I_D + G_D + Phi_D + T_D) − NX_D` ✓ algebraically = D-goods market clearing after substituting C_DD + p·IM_D = P_CES_D·C_D.

## 6. System assembly (transition)

23 unknowns / 23 targets (notebook cell 20). Targeted residuals: deposit markets (both), bank balance sheet & net worth & dividends (both), capital accumulation & Q (both), government budgets (both), cross-bond FOCs (both), labor FOCs (both ×2), `goods_mkt_D`, own-bond FOCs (both).
**Untargeted (must hold by Walras): `goods_mkt_F`, `ca_res_D`** (`global_goods_res` = goods_mkt_D + p·goods_mkt_F is a linear combination, not independent).

SS solve: 3 unknowns (β_D, β_F, p), 3 targets (deposit_mkt_D, deposit_mkt_F, ca_res_D). `rdep` fixed at 0 at SS; `goods_mkt_D/F` *not* SS targets (rely on SS Walras).

## 7. Documentation vs code discrepancies

| # | Doc claim | Code reality |
|---|---|---|
| 1 | `model_notes`: "portfolio adjustment cost from deviating from steady state level of foreign bonds holding (simplified)" | Implemented as **level** penalty on face value in `divert_portfolio_adj` (consistent), but own-bond FOCs additionally carry share-based penalty with ψ=0 (inert) |
| 2 | calibration comment "beta_inter → 1.7%/yr govt bond yield" | check numerically; SS rb_actual = δ_b(1/q−1) with q from β_inter pricing ⇒ ≈1.0%/yr (β_inter=0.99752, δ_b=0.1) — verify in run |
| 3 | calibration comment `Delta_bD_D: 0.2 (preferred)` etc. | Values **overwritten** by cell-12 back-solve; printed values materially different (incl. >1) |
| 4 | `lambda_gk_D: 0.2` calibration entry | Overwritten by `steady_auxilliary` SS output in `ha`; in `ha_full` it is a constant at the *solved* value, not 0.2 |
| 5 | `bond_yield` comment: "old formula overstated yield ~20×" | current formula correct for this perpetuity ✓ |
