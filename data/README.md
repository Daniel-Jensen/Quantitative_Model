# Calibration data sources

Raw sources are **not tracked** (large; see `.gitignore`). Derived
`eba_moments.json` **is** tracked. Regenerate with:

```bash
/opt/anaconda3/envs/ssj/bin/python code/eba_calibration.py
```

## Files

| File | Source | Notes |
|------|--------|-------|
| `DATA_DISCLOSURE.CSV (1).csv` | EBA 2011 EU-wide stress test disclosure database (published 15 Jul 2011) | latin-1 encoded. Worksheets 1–5; base 31 Dec 2010 actual + 2011/2012 baseline(100)/adverse(105) projections. Sovereign exposures as-of 31 Dec 2010. |
| `EBA_2011_stress_test_summary_report.pdf` | eba.europa.eu | Methodology + aggregate results. |
| `EBA_2011_disclosure_templates.pdf` | eba.europa.eu | Disclosure template definitions (code map reference). |
| `eurostat_gdp_raw.json` | Eurostat `nama_10_gdp` (B1GQ, CP_MEUR), geo EL+DE, 2009–2013 | Nominal GDP for GDP-normalization. |
| `../Empirics/outputs/spreads_fred.csv` | FRED long-term (10y) government bond yields, monthly | Sovereign yields for the maturity-ladder repricing (duration → `delta_b`) and the Acharya–Steffen MTM moment. Already tracked under `Empirics/`. |
| `eba_moments.json` | **derived** by `code/eba_calibration.py` | Model calibration targets, maturity ladder, MTM moments, and the **identification ledger**. **Tracked** — and read directly by `code/calibration.py` / `code/steady_state.py`, which must never carry their own copy of these numbers. |

## Code map (validated)

`DE017 = Deutsche Bank` (2010 RWA €346.6 bn, CT1 €30.4 bn, ratio 8.76% — matches published):

- Worksheet 1, MEASURE 120: `30010`=RWA, `30014`=Core Tier 1, `30029`=total assets
- Worksheet 4 (credit-risk EAD by counterparty country, `TO`=total):
  `33010`=Institutions, `33011`=Corporate ex-CRE, `33012`=Retail ex-CRE,
  `33018`=Commercial Real Estate, `33021`=Total exposures
- Worksheet 5 (sovereign, by issuer × residual maturity): `34010`=gross direct
  long, `34011`=of which loans & advances, `34012`=net direct positions,
  `34013`=AFS banking book, `34014`=FVO banking book, `34015`=trading book,
  `34016`=derivatives, `34017`=indirect trading-book
- `MATURITY_CODE`: `125`..`155` are the seven residual-maturity buckets (upper
  bounds 3M/1Y/2Y/3Y/5Y/10Y/15Y); `999` is the total. The seven buckets sum
  exactly to the total — asserted by `test_ladder_exhausts_total`.

**Not used, deliberately:** `SCENARIO_CODE=105` (adverse) CT1 depletion. The 2011
exercise excluded sovereign default in the banking book, so its capital
depletion understates sovereign pass-through by construction and must not
identify the model's amplification. Guarded by `test_adverse_scenario_not_used`.

## Country map

D = Greece (banks GR030–GR035, issuer GR) · F = Germany (banks DE0xx, issuer DE)

## Re-downloading raw sources

```bash
# Eurostat GDP
curl -sSL "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/nama_10_gdp?format=JSON&na_item=B1GQ&unit=CP_MEUR&geo=EL&geo=DE&time=2009&time=2010&time=2011&time=2012&time=2013" -o data/eurostat_gdp_raw.json
# EBA PDFs
curl -sSL "https://www.eba.europa.eu/documents/10180/15935/54a9ec8e-3a44-449f-9a5f-e820cc2c2f0a/EBA_ST_2011_Summary_Report_v6.pdf" -o data/EBA_2011_stress_test_summary_report.pdf
```
The EBA disclosure CSV is provided by the user (EBA 2011 stress-test results database).
