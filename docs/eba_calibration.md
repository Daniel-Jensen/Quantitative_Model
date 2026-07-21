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

### Why C-1 forces `Delta→1` (analytical; explains the collateral flip)

C-1 is **not** a calibration input — it is an artifact of an inconsistency between the SS
and dynamic blocks. `steady_auxilliary` builds `lambda_gk` from the *single-asset* GK
formula (`lambda_gk = f/(θ·(1/(β_inter(1+rn)) − (1−f)))`), which forces the identity
`value/lambda_gk = θ`. The *dynamic* IC is multi-asset with divertability weights `Delta`.
Substituting the identity into `ic_delta_calibration._ic_delta` collapses the back-solve to
a function of the portfolio ratio alone:

```
Delta_own = (phi_own + phi_cross) / (phi_own + ratio·phi_cross)        [ratio = 2.0, hardcoded]
```

Verified exactly against all three observed cases:

| case | phi_own | phi_cross | predicted `Delta_cross` | observed |
|---|---|---|---|---|
| baseline (pre-EBA) | 0.250 | 0.150 | 1.4545 | **1.45** (C-1 in CLAUDE.md) |
| EBA D-bank | 2.390 | 0.018 | 1.9852 | **1.9852** |
| EBA F-bank | 2.760 | 0.069 | 1.9524 | **1.9524** |

So `Delta_cross > 1` whenever `Delta_own > 0.5` — i.e. essentially always. The pre-EBA
calibration merely sat far enough from the boundary (`Delta_own=0.73`) to limp; EBA's
realistic concentration (`phi_own=2.39`) drives `Delta_own→0.99`, and `Delta_eff = Delta +
psi·def_rate` then crosses 1, flipping the sign of the collateral channel. **This is exactly
why `psi_lambda_B` had to be rescaled 3.0 → 0.31** to restore correct signs.

*Permanent fix (not yet done):* solve `lambda_gk` from the multi-asset IC jointly —
`lambda_gk = value/(θ − (1−Delta_own)phi_own − (1−Delta_cross)phi_cross)`, a fixed point since
`value` depends on `Omega(lambda_gk, θ)`. That would make `Delta` a genuine free parameter and
let `Delta_D=0.2 / Delta_F=0.4` be hardcoded as CLAUDE.md prefers.

*Ruled out:* using EBA **total-asset** leverage as `θ` (GR 16.56, DE 42.62). `θ` multiplies
only the GK book (capital+sovereign), whereas EBA total assets include low-yield
loans/reserves; `θ=16.56` on the model's fixed `rk−rdep=0.74%` spread implies ~52% annual
banker ROE and the SS **does not converge** (verified). Do not retry.

### Open issue: the EBA sovereign doom loop is intrinsically explosive here

Full stationarity is **not** achieved. `phi_lamb` (fiscal feedback) damps the amplitude
(`K_D[499]` 9.5→3.4 as `phi_lamb` 0.3→0.6) but leaves a residual ~1.2%/period explosive
root. Two capital-structure variants were tested; **neither** removes the root:

- **Fixed `omega_K`** (bank holds `omega_K·K`): `K = (θN−bonds)/(omega_K·Q)` levers `K`
  `1/omega_K≈16–67×` on thin net worth. *Least unstable of the two* (default-shock
  `b_gov[499]=0.17` at `phi_lamb=0.6`). **This is the current committed state.**
- **Endogenous `omega_K` = fixed fund quantity** (bank holds `Kbank=K−Kfund_ss`): tested
  (commit history) and made it **worse** (`b_gov[499]=79.5`). `Kbank=0.65` is a small
  difference of large numbers (`10.8−10.15`), so it becomes hyper-sensitive to `K` and
  `kappa=Kbank·Q/N` inherits it. Reverted.

**Conclusion:** the instability is not the capital plumbing — it is the **strong EBA
sovereign doom loop itself** (`phi_bD_D=2.39`: Greek banks held 2.39× their capital in
Greek debt). Both variants are correctly signed (spread↑, `q_b`↓, net worth↓ on default)
and Walras-clean; both mildly explode.

**This is arguably economically meaningful, not just a nuisance:** a self-reinforcing,
(locally) explosive doom loop is *precisely the condition a backstop is meant to
neutralise*. The natural next step — feeding the measured **SMP purchase path through the
TPI conduit (phase 2)** — is the candidate stabiliser, and the instability strengthens
that narrative rather than undercutting it.

**Options for the author:**
1. **Proceed to phase 2 (SMP conduit)** and test whether the backstop restores
   stationarity — the intended paper mechanism. *Recommended.*
2. **Microfound the fund** (Gertler-Karadi-Prestipino management cost) for a downward-
   sloping capital demand — the rigorous fix, real equation work.
3. **Relax the 2.39× target** toward a value where the baseline is stationary — revisits
   the "faithful to EBA" cross-holding choice.
4. **Report bounded-horizon IRFs** where signs are correct and Walras holds (weakest).

## Reproduce

```bash
/opt/anaconda3/envs/ssj/bin/python code/eba_calibration.py     # moments
/opt/anaconda3/envs/ssj/bin/python code/main.py                # full pipeline
```
