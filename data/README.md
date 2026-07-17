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
| `eba_moments.json` | **derived** by `code/eba_calibration.py` | Model calibration targets + provenance. **Tracked.** |

## Code map (validated)

`DE017 = Deutsche Bank` (2010 RWA €346.6 bn, CT1 €30.4 bn, ratio 8.76% — matches published):

- Worksheet 1, MEASURE 120: `30010`=RWA, `30014`=Core Tier 1, `30029`=total assets
- Worksheet 5, `34010`, MATURITY 999: gross direct long sovereign exposure (total)

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
