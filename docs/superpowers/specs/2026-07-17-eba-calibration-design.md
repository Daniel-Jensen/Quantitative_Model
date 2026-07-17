# EBA-anchored recalibration + measured SMP deployment — design

**Date:** 2026-07-17
**Branch:** `eba-calibration` (from `ecb-balance-sheet`)
**Status:** Approved (design), building

## Amendments (during build)

- **Exercise identified:** EBA **2011 EU-wide stress test** (published 15 Jul 2011). Base
  year **31 Dec 2010 actual** (scenario blank); 2011/2012 are baseline(100)/adverse(105)
  *projections*. Sovereign exposures are as-of **31 Dec 2010** per methodology (the
  worksheet-5 `INFO_DATE_CODE=20111231` is a reference-year label, not the exposure date).
  → **Calibrate to the 2010 actual base with 2010 GDP.**
- **Code map validated:** `DE017 = Deutsche Bank` (2010 RWA €346.6 bn, CT1 €30.4 bn,
  8.76% — matches published). Confirms `30010=RWA`, `30014=Core Tier 1`,
  `34010=gross direct long sovereign`, `30029=total assets`.
- **Bank identification fork (resolved: option A).** EBA over-identifies the bank block:
  sovereign-book/capital = 2.39× (GR) needs thin net worth, but the model ties net worth to
  intermediating the *entire* capital stock `K`. Resolution: **add a bank
  capital-intermediation share `omega_K_D/F`** — banks hold `omega_K·K`, the rest is held by
  a passive unlevered capital fund whose return is rebated to households (no two-asset HANK
  portfolio choice). Lets us match `n_inter=CT1/GDP`, the 2.39× ratio, K/Y, and stability
  simultaneously. Touches `equations_D/F.py`, `steady_state.py`, household income accounting.

## Locked 2010-base moments (EUR m; GDP: GR 223,590.5 / DE 2,615,260 annual)

| | Greek banks (D) | German banks (F) |
|---|---|---|
| CT1 (30014) | 22,778 | 114,317 |
| RWA (30010) | 222,466 | 1,222,402 |
| Total assets (30029) | 377,200 | 4,872,189 |
| own-sovereign book (34010) | 54,447 (GR) | 315,313 (DE) |
| cross-sovereign book (34010) | 411 (DE) | 7,934 (GR) |

Derived model targets (symmetric Y=1): `n_inter_D=0.408`, `n_inter_F=0.175`,
`phi_bD_F_ss=0.069`, `phi_bF_D_ss=0.018`, `phi_bD_D=2.39`, `phi_bF_F=2.76`,
leverage TA/CT1 = 16.6×/42.6×, bank-held GR debt = 27.9% ann. GDP → `B_supply_D≈1.12`.

## Objective

Use `data/DATA_DISCLOSURE.CSV` (EBA 2011 EU-wide stress-test disclosure) to
(1) **correctly calibrate** the two-country model to measured 2011 bank/sovereign
data, and (2) *(phase 2)* feed the **actual ECB SMP purchase path** for Greek debt
through the TPI conduit as the *measured* deployment case, replacing the synthetic
γ-sweep with a data-anchored point.

**Priority: calibration first (phase 1).** SMP feed is phase 2.

Country map: **D = Greece** (banks `GR030–GR035`, issuer `GR`),
**F = Germany** (banks `DE017–DE029`, issuer `DE`).

## Data source

`data/DATA_DISCLOSURE.CSV` — EBA 2011 disclosure, `latin-1` encoded, 3456 rows × 97 cols.
Worksheets: `1 - Aggregate information` (capital/RWA/assets), `2 - Capital composition`,
`3 - Mitigating measures`, `4 - EADs`, `5 - Sovereign exposures`.
Row key = WORKSHEET + INFORMATION_CODE + INFO_DATE_CODE + SCENARIO_CODE + MATURITY_CODE +
MEASURE_CODE; columns are banks. Sovereign worksheet is a **2011-12-31 snapshot**
(post-first-PSI-provisioning); aggregate info covers 2010/2011/2012 × baseline(100)/adverse(105).

**Inferred code map (to be confirmed against the official EBA 2011 data dictionary — download step):**
- `34010` = gross direct long sovereign exposure (worksheet 5)
- `30010` = total RWA, `30014` = Core Tier 1 capital, `30029` = total assets (worksheet 1, MEASURE 120)

## Measured moments (EBA 2011, €m, baseline)

Cross-holding matrix (worksheet 5, `34010`, maturity 999, summed over bank groups):

| Holder \ Issuer | Greek sovereign | German sovereign |
|---|---|---|
| Greek banks (D) | 54,447 | 411 |
| German banks (F) | 7,934 | 315,313 |

Bank aggregates (worksheet 1, MEASURE 120):

| Moment | Greek banks (D) | German banks (F) |
|---|---|---|
| Core Tier 1 (30014) | 22,590 | 116,390 |
| Total RWA (30010) | 228,551 | 1,353,393 |
| Total assets (30029) | 371,100 | 4,813,251 |

## Scope (approved)

In: **cross-border sovereign holdings**, **bank capital & leverage**, **debt/GDP & PSI recovery**.
Out: bond duration from maturity buckets (deferred; also F-1 explosive).

Cross-holdings calibrated **faithful to EBA 2011** (risk-premium framing).

## Decisions (approved recommendations)

**(a) Country-size normalization.** Keep symmetric `Y_D = Y_F = 1`; match each country's
*ratios to its own nominal GDP*. No asymmetric sizing.

**(b) Debt/GDP reconciliation.** The model forces the two banks to hold *all* government
debt (`b_D_D + b_D_F = B_supply_D`), but GR+DE banks held only ~€62 bn of Greece's ~€356 bn
(~170% GDP). **Now:** set `B_supply`/`b_gov` to the *bank-held* stock (~30% GDP) — realistic
bank channel; "debt/GDP" in-model ≠ headline. **Later:** revisit full 170% with a non-bank
residual holder = the SMP/official sector (couples to phase 2).

**(c) PSI recovery vs writeoff.** Calibrate PSI as an **SS provisioning write-down of the
sovereign book** (consistent with `writeoff_enabled=0`), NOT by enabling the dynamic writeoff
channel (which per S-1/F-1 makes the default response perverse under the market-value rule).
Writeoff stays OFF.

## Parameter mapping

Cross-holding shares (GDP-free; `phi = q·book / net-worth`):
- `phi_bD_F_ss` = q_b_D·(DE banks' GR book)/(DE CT1) = 7,934/116,390 ≈ **0.068** (was 0.25)
- `phi_bF_D_ss` = q_b_F·(GR banks' Bund)/(GR CT1) = 411/22,590 ≈ **0.018** (was 0.25)

Level moments (need nominal GDP → download):
- `n_inter_D` = CT1_GR / quarterly-GDP_GR ≈ 22,590/51,000 ≈ **0.44** (was 3.0)
- `n_inter_F` = CT1_DE / quarterly-GDP_DE (per-GDP ratio, symmetric-normalized)
- `theta` / leverage from total-assets / CT1 (GR ≈ 16×, DE ≈ 41×)
- `B_supply_D`, `b_gov_D` from bank-held Greek stock / GDP (decision b1)
- `recovery_rate_D` / SS book write-down from 2011 PSI provisioning (decision c)

## Architecture

- **`code/eba_calibration.py`** — reads the CSV, decodes via an explicit `CODE_MAP`
  dict (confirmed against the data dictionary), computes the moment set, writes
  `data/eba_moments.json` with values **and provenance** (worksheet/code/date/scenario
  for every number). Pure, testable, no side effects beyond the JSON.
- **`calibration.py`** — recalibrated in **stages**, each changed value commented with
  its EBA moment + `eba_moments.json` key.
- **`docs/eba_calibration.md`** — maps every changed parameter to its EBA moment and
  records the GDP sources.

## Downloads (into `data/`, for user review)

1. EBA 2011 EU-wide stress-test **data dictionary** (confirm code meanings).
2. **GR & DE nominal GDP 2010–2012** (Eurostat/AMECO) for level normalization.
3. *(phase 2)* ECB weekly **SMP** settlement series (ECB Data Portal / weekly financial statements).

## Staged execution + stability gates

After **each** stage: run `code/main.py`; assert `goods_mkt_D ≤ 1e-14`, `goods_mkt_F ≤ 1e-7`,
`ca_res_D ≤ 1e-7`, `deposit_mkt_D/F ≤ 1e-13`; default-shock sign check (`n_inter_D[0]`,
`Y_D[0]` both negative).

- **Stage 1 — clean wins:** cross-holdings (`phi_bD_F_ss`, `phi_bF_D_ss`, `b_*` book) +
  bank capital/leverage (`n_inter_D/F`, `theta`). No stability risk expected.
- **Stage 2 — debt/GDP (gated):** rescale `B_supply`/`b_gov` to bank-held stock. Check
  stationarity; if explosive, `mv_rule=1` is the documented fix.
- **Stage 3 — PSI recovery:** SS book write-down; `writeoff_enabled` stays 0.
- **Phase 2 — measured SMP:** reconstruct quarterly Greek SMP cumulative path from ECB
  data, feed as exogenous `cb_buy_D` (open-loop), report where it lands on the loading
  schedule ℓ(γ). Expected: far right of the current γ = {0,2,5,10} grid.

## Deliverables

`code/eba_calibration.py`, `data/eba_moments.json` + downloaded sources, recalibrated
`code/calibration.py` (staged, provenance-commented), `docs/eba_calibration.md`,
updated `docs/STATE.md`.
