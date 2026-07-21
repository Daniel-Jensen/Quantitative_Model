# EBA-anchored calibration

Source: EBA **2011 EU-wide stress test** disclosure (`data/DATA_DISCLOSURE.CSV`),
base **31 Dec 2010 actual**; nominal GDP from Eurostat `nama_10_gdp` (2010).
Country map: **D = Greece** (banks GR030–GR035, issuer GR), **F = Germany**
(banks DE0xx, issuer DE). Moments produced by `code/eba_calibration.py` →
`data/eba_moments.json`. Code map validated: `DE017 = Deutsche Bank`
(2010 RWA €346.6 bn, CT1 €30.4 bn, 8.76% — matches published).

## Parameter → moment map

| Parameter | Old | EBA value | Moment (EUR m, 2010) |
|---|---|---|---|
| `n_inter_D` | 3.0 | **0.408** | CT1_GR 22,778 / quarterly GDP 55,898 |
| `n_inter_F` | 3.0 | **0.175** | CT1_DE 114,317 / quarterly GDP 653,815 |
| `phi_bD_D_ss` (GR banks' GR book / capital) | 0.25 | **2.39** | 54,447 / 22,778 |
| `phi_bF_F_ss` (DE banks' Bund / capital) | 0.25 | **2.76** | 315,313 / 114,317 |
| `phi_bD_F_ss` (DE banks' GR / capital, contagion) | 0.15 | **0.069** | 7,934 / 114,317 |
| `phi_bF_D_ss` (GR banks' Bund / capital) | 0.15 | **0.018** | 411 / 22,778 |
| `B_supply_D`, `b_gov_D` (bank-held GR debt) | 2.4 | **1.19** | 62,381 / qGDP (27.9% ann. GDP) |
| `B_supply_F`, `b_gov_F` (bank-held Bund) | 2.4 | **0.591** | 323,247 / qGDP (12.1% ann. GDP) |

`phi` = q·(sovereign book)/(bank net worth). The own-holding ratios (2.39×, 2.76×)
are the doom-loop vulnerability; cross-holdings (0.069, 0.018) are the thin direct
contagion channel — a **10–30× reduction** from the previous 0.25 placeholder.

**Design note (b1):** `B_supply` is the *bank-held* sovereign stock (~28%/12% of GDP),
NOT headline debt/GDP (~150% GR). The non-bank/official/ECB residual is intended to
enter via the SMP conduit (phase 2). PSI recovery (Stage 3) is deferred.

## Structural change: `omega_K` capital-intermediation fund

EBA over-identifies the bank block: matching the 2.39× sovereign/capital ratio needs
thin net worth (`n_inter_D=0.408`), but the model ties net worth to intermediating the
*entire* capital stock. Resolution (user-approved option A): banks hold `omega_K` of
capital, a **passive deposit-funded capital fund** holds `(1-omega_K)K` and rebates its
spread `(rk-rdep)` to households as `div_fund` (`capital_fund_D/F`, into `income_D/F`).
`omega_K` set to preserve `K≈10.8` (K/annualY≈2.7): `omega_K = N·(θ-φ_own-φ_cross)/(Q·K_target)`
→ **0.060 (D), 0.019 (F)**.

**Verification:** the change is Walras-clean to machine precision at SS **and** along
IRFs (`ca_res_D`, `goods_mkt_F` ~1e-8; solved targets ~1e-16). Verified as a no-op at
`omega_K=1` (reproduces the pre-EBA SS exactly).

## Stability: what works, what's open

With `psi_lambda_B` rescaled **3.0 → 0.31** (doom-loop amplification ~ `psi·phi_own`;
`phi_own` rose ~10×, so `psi` falls ~10× to keep comparable amplification and keep
`Delta_eff = Delta + psi·def_rate` below 1), plus `mv_rule=1` and `phi_lamb=0.60`:

- ✅ **Doom-loop signs correct**: a default shock gives spread ↑, `q_b` ↓, `rb` ↓,
  bank net worth ↓ (previously inverted by the C-1 `Delta→1` collateral flip).
- ✅ **Walras clean** on all IRFs.
- ⚠️ **Residual ~1.2%/period explosive root** (`K_D`, `b_gov_D` drift up to t=499).
  `phi_lamb` reduces the amplitude (`K_D[499]` 9.5→3.4 as `phi_lamb` 0.3→0.6) but not
  the root.

### Open issue: fixed-`omega_K` over-couples K to bank net worth

With `omega_K` a **fixed parameter**, `k_balance_sheet` pins total capital
`K = (θ·N − bonds)/(omega_K·Q)`, so `K` is levered `1/omega_K ≈ 16–67×` on the bank's
thin net worth — any `N`/bond fluctuation is amplified into `K` and mildly explodes.
Economically, the passive fund holds a *fixed share* with no downward-sloping capital
demand; an unconstrained fund earning `rk−rdep` would arbitrage without a friction
(the Gertler-Karadi rationale for why only constrained banks hold capital).

**Candidate fixes (decision pending):**
1. **Microfound the fund** — households/fund hold `(1-omega_K)K` subject to a management
   cost (Gertler-Karadi-Prestipino), giving a downward-sloping capital demand so `K` is
   pinned by the real block, not by `1/omega_K`. Cleanest; real equation work.
2. **Endogenous `omega_K`** — make `omega_K` an unknown pinned by `k_balance_sheet`, with
   `K` pinned by capital accumulation + Q Euler. Decouples `K` from the amplification but
   needs the SSJ unknown/target rebalance to be well-posed.
3. **Reduced-form + accept** — keep fixed `omega_K`, treat the slow root as a truncation
   artifact and report bounded-horizon IRFs only. Weakest.
4. **Relax the 2.39× target** — a smaller sovereign/capital ratio needs less-thin capital
   and less-extreme `omega_K` (revisits the user's cross-holding-faithfulness choice).

## Reproduce

```bash
/opt/anaconda3/envs/ssj/bin/python code/eba_calibration.py     # moments
/opt/anaconda3/envs/ssj/bin/python code/main.py                # full pipeline
```
