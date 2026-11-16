# Greek (D) / German (F) calibration data, 2008–2013: macro, trade, sovereign risk, household heterogeneity, and literature values

Conventions in these notes. "pp" = percentage points; 1 pp = 100 bp. "Hazard" = per-period default probability. Eurostat figures were pulled live from the Eurostat dissemination API on 2026-10-23 (current vintage, ESA 2010, chain-linked 2010 euros where real). They will differ slightly from 2013-vintage figures in IMF and EC programme documents. World Bank figures are from the WDI API, same date. Values marked "computed" are my own arithmetic on the cited series. The repository file `Empirics/outputs/spreads_fred.csv` holds monthly averages of the 10-year benchmark yields (OECD/FRED "IRLTLT01" series) and spreads over the Bund in pp.

## 1. Sovereign risk: yields, spreads, CDS, default probabilities, persistence, recovery

### Takeaway
On a 10-year basis the Greek spread over Bunds averaged 6.35 pp in 2010, 13.1 pp in 2011 and 21.0 pp in 2012, and peaked at 27.4 pp (monthly average) in February 2012, just before the PSI. With Bocola's haircut D = 0.55, the 2010 average implies roughly a 2.9% quarterly (11.5% annual) default hazard, against Bocola's own Italian ergodic mean of about 0.09% per quarter. The PSI wiped out 53.5% of face value. The Zettelmeyer–Trebesch–Gulati present-value haircut was 59–65% (headline 64.6%), and the CDS auction set recovery at 21.5%.

### Cited Findings
**10-year yields and spreads (monthly averages, from the repository's `Empirics/outputs/spreads_fred.csv`, OECD long-term rates via FRED; computed)** — [FRED IRLTLT01GRM156N](https://fred.stlouisfed.org/series/IRLTLT01GRM156N), [FRED IRLTLT01DEM156N](https://fred.stlouisfed.org/series/IRLTLT01DEM156N)

| Year | GR 10y (%) | DE 10y (%) | Spread (pp) |
|---|---|---|---|
| 2008 | 4.80 | 3.98 | 0.82 |
| 2009 | 5.17 | 3.22 | 1.95 |
| 2010 | 9.09 | 2.74 | 6.35 |
| 2011 | 15.75 | 2.61 | 13.14 |
| 2012 | 22.50 | 1.50 | 21.00 |
| 2013 | 10.05 | 1.57 | 8.48 |

- Monthly spread path in pp (same file): 2008-01 0.37, 2008-09 0.79, 2008-12 2.03, 2009-03 2.85, 2009-08 1.21 (2009 trough), 2009-12 2.35, 2010-01 2.76, 2010-04 4.77, 2010-05 5.24 (first programme), 2010-06 6.56, 2010-09 9.04, 2010-12 9.10, 2011-03 9.23, 2011-06 13.80, 2011-09 15.95, 2011-12 19.21, 2012-01 24.09, **2012-02 27.39 (peak)**, 2012-03 17.23 (after the exchange), 2012-06 26.52 (election-period second peak), 2012-09 19.42, 2012-12 12.03, 2013-06 8.54, 2013-12 6.86. Peak Greek 10-year yield: 29.24% (Feb 2012). Caveat: after March 2012 the series refers to the new post-PSI bonds.
- Pre-crisis calm (2001–2007 monthly): mean spread 0.297 pp, s.d. 0.107 pp. Computed from the same file.

**AR(1) fits on the spread (computed from the same file; OLS x_t = c + ρ x_{t−1} + e, σ = s.d. of e, in pp)**

| Sample, frequency | ρ | σ_e (pp) | mean (pp) | n |
|---|---|---|---|---|
| 2001–2019, monthly | 0.978 | 1.17 | 4.68 | 226 |
| 2001–2019, quarterly averages | 0.955 | 1.64 | 4.70 | 75 |
| 2008–2013, monthly | 0.956 | 1.98 | 8.62 | 71 |
| 2008–2013, quarterly | 0.907 | 2.72 | 8.62 | 23 |
| 2009–2012, quarterly | 0.905 | 3.01 | 10.61 | 15 |
| 2001–2007, quarterly | 0.759 | 0.044 | 0.30 | 27 |

**In Bocola's logit units.** Bocola models the default probability as pd = e^s/(1+e^s) with s an AR(1). I converted the quarterly-average spread to a quarterly hazard, spread/(4·D) with D = 0.55, and took the logit (computed):
- 2008–2012: ρ_s = 0.904, σ_s = 0.318, mean logit −3.82 (pd ≈ 2.1%/qtr), n = 19.
- 2008–2013: ρ_s = 0.876, σ_s = 0.314, mean −3.72, n = 23.
- 2001–2012: ρ_s = 1.018 (unit root, because the sample is a regime switch), σ_s = 0.284, mean −5.48.
- Bocola's own Italian estimates for comparison: ρ_s = 0.95, σ_s = 0.63, s* = −7.06 (see Section 6).

**Spread-implied default hazards (computed, hazard ≈ spread / D, quarterly ≈ annual/4)**

| Loss-given-default D | 2010 annual / qtr | 2011 annual / qtr | 2012 annual / qtr | Feb-2012 peak annual / qtr |
|---|---|---|---|---|
| 0.55 (Bocola) | 11.5% / 2.9% | 23.9% / 6.0% | 38.2% / 9.6% | 49.8% / 12.4% |
| 0.646 (ZTG PV haircut) | 9.8% / 2.5% | 20.3% / 5.1% | 32.5% / 8.1% | 42.4% / 10.6% |
| 0.785 (1 − 21.5% CDS recovery) | 8.1% / 2.0% | 16.7% / 4.2% | 26.8% / 6.7% | 34.9% / 8.7% |

**CDS**
- 5-year Greek CDS was about 420 bp in late January 2010, about 698 bp in late April 2010 and about 1,003 bp in mid-June 2010. These are point-in-time figures from search-result snippets across market blogs and news ([ritholtz.com](https://ritholtz.com/?p=50377), [hedgeco.net](https://www.hedgeco.net/blogs/?p=1744), [Central Banking](https://www.centralbanking.com/node/1720902)); the individual pages were not opened, so treat them as indicative.
- 5-year CDS was 1,375 bp in May 2011 after the S&P downgrade. — [Calculated Risk, May 2011](https://calculatedriskblog.com/2011/05/another-downgrade-for-greece.html)
- A March 2012 Fitch poll put the 6-month contract above 20,000 bp, effectively pricing an imminent default. — [Institutional Investor](https://institutionalinvestor.com/article/b14zplwg719r9q/greek-cds-holders-cautiously-optimistic-before-auction)
- ISDA declared a restructuring credit event on 9 March 2012, triggered by Greece's use of collective action clauses. The 19 March 2012 auction set a **final price (recovery) of 21.5 cents**, so 78.5% of notional was paid out. Net payouts were about €2.89bn, and gross notional was about €80bn, of which roughly US$3bn was net. — [Orrick derivatives blog](https://blogs.orrick.com/derivatives/?p=8); [Institutional Investor](https://inv.institutionalinvestor.com/article/2bsvlexhhb6xk5hg0d98g/portfolio/greek-cds-holders-cautiously-optimistic-before-auction). ZTG give the same 21.5-cent auction price, which matched the price of the cheapest new (2042) bond. — [ZTG, CESifo WP 4333](https://www.cesifo.org/DocDL/cesifo1_wp4333.pdf)

**2012 PSI haircut (Zettelmeyer, Trebesch & Gulati 2013, Economic Policy; CESifo WP 4333)** — [ZTG CESifo WP 4333](https://www.cesifo.org/DocDL/cesifo1_wp4333.pdf); [Economic Policy draft](https://www.economic-policy.org/wp-content/uploads/2013/04/zettelmeyer.pdf)
- Exchange of about €200bn in March/April 2012. Participation was €199.2bn, or 96.9% of eligible principal. Greek-law bonds were €177.3bn, of which 82.5% consented before the CACs were applied. Foreign-law participation was 71% in the end, and holdouts held €6.4bn. Per 100 of old face value, investors received 15 in short-term EFSF notes plus 31.5 in new long bonds (plus GDP warrants worth about 0.23). **Face value haircut = 100 − 46.5 = 53.5%**; ZTG report a face-value decline of about €107bn, or 52% of eligible debt.
- **ZTG Table 3**, present-value haircut 100·(1 − PV_new/PV_old) at the stated discount rate:
  - At the 15.3% average exit yield: 64.6%
  - At 16.3%: 64.4%
  - At 18.7%: 64.0%
  - Imputed yield curve: 59.6%, with sensitivity range of about 55–61%.
- **ZTG Table 3**, "market haircut" (100 − PV_new): 76.9 / 77.5 / 78.8 / 77.2, respectively. PV_new is 21.2–23.1 per 100 of old face. ZTG give "59–65%" as the range and use 64.6% for historical comparison.
- The weighted average exit yield was 15.3%, and new bonds traded at roughly 14–17.5% on 12 March 2012.
- Losses varied by maturity: 75% or more on bonds maturing within a year, under 50% on bonds maturing after 2025.
- The exchange transferred about €100bn in PV to Greece, roughly 50% of 2012 GDP net of bank recapitalisation costs.
- The December 2012 buyback was done at about 34 cents per unit of new face value. Combining the exchange and the buyback, a participating creditor received 15 + 0.315·34 = 25.7 cents per original euro, against about 73.9 for holding out, a **combined haircut of about 65%**. Private bond holdings fell from about €206bn (Feb 2012) to €35bn.
- Hardouvelis & Vayanos count a face-value fall from €205.6bn to €98.5bn (−52.1%). — [Hardouvelis & Vayanos, LSE GreeSE 180](https://www.lse.ac.uk/asset-library/hellenic-observatory/greese-no180.pdf)

### Inferences
- The 10-year spread is not the one-period hazard. With prices far below par and short-dated bonds suffering larger PV losses (ZTG), the spread/D rule overstates term-structure-consistent near-term risk in 2012. A 1-year CDS series, Bocola's choice, would be cleaner.
- In the model's units, the 2010 episode is roughly a 2.5–3% quarterly priced hazard. That is about 30× Bocola's Italian ergodic mean (0.09%) and similar to his "high s" conditioning value (2.5% per quarter).
- The logit AR(1) on 2008–2013 (ρ ≈ 0.88–0.90, σ ≈ 0.31–0.32) has lower persistence and half the volatility of Bocola's Italian estimate (ρ = 0.95, σ = 0.63). Longer samples give ρ near or above 1 because the pre-2008 and post-2008 regimes stack. Either value can be defended, so the choice should be flagged to the user rather than defaulted.
- For the "pure haircut" default branch, the defensible range for D is 0.535 (face) to 0.646 (PV at exit yield). Bocola's 0.55 sits at the low end, and the CDS recovery of 21.5% implies an upper value of 0.785.

### Gaps
- I did not obtain a full daily or monthly 5-year or 1-year Greek CDS series (Markit/Bloomberg are proprietary). Only press point values are cited above. No annual-average CDS figure from an official source was found.
- No published AR(1) estimate of Greek default probability was found. The fits above are my own, on 10-year spreads.

## 2. Fiscal: debt, primary balance, spending, revenue, fiscal reaction

### Takeaway
Greek gross debt rose from 110.9% of GDP (2008) to 175.1% (2011). It dipped to 164.1% after the PSI (2012) and rose again to 180.4% (2013). The primary deficit peaked at about 10.3% of GDP in 2009. Total revenue was about 41–50% of GDP and total expenditure about 51–64%. I found no published Greece-specific Bohn coefficient. Euro-area panel estimates are 0.03–0.05.

### Cited Findings
**Eurostat general government (S13), % of GDP, current vintage** — [Eurostat gov_10dd_edpt1](https://ec.europa.eu/eurostat/databrowser/view/gov_10dd_edpt1/default/table), [gov_10a_main](https://ec.europa.eu/eurostat/databrowser/view/gov_10a_main/default/table), [gov_10a_taxag](https://ec.europa.eu/eurostat/databrowser/view/gov_10a_taxag/default/table)

| % GDP | 2007 | 2008 | 2009 | 2010 | 2011 | 2012 | 2013 |
|---|---|---|---|---|---|---|---|
| EL gross debt | 104.6 | 110.9 | 128.5 | 147.8 | 175.1 | 164.1 | 180.4 |
| DE gross debt | 63.7 | 65.2 | 72.4 | 81.0 | 78.5 | 79.8 | 77.4 |
| EL net lending (B9) | −6.8 | −10.3 | −15.4 | −11.4 | −10.5 | −9.3 | −13.6 |
| DE net lending | 0.2 | −0.3 | −3.2 | −4.4 | −0.8 | −0.1 | 0.1 |
| EL interest paid (D41) | 4.6 | 4.9 | 5.1 | 6.1 | 7.7 | 5.4 | 4.2 |
| **EL primary balance (B9+D41, computed)** | −2.2 | −5.4 | −10.3 | −5.3 | −2.8 | −3.9 | −9.4 |
| DE primary balance (computed) | 2.9 | 2.4 | −0.6 | −2.0 | 1.7 | 2.2 | 1.9 |
| EL total revenue (TR) | 41.0 | 41.2 | 39.4 | 41.7 | 44.6 | 48.6 | 50.3 |
| EL total expenditure (TE) | 47.8 | 51.5 | 54.8 | 53.1 | 55.1 | 57.9 | 63.9 |
| DE TR / TE | 43.7/43.5 | 44.2/44.4 | 45.1/48.3 | 43.7/48.1 | 44.4/45.3 | 45.0/45.1 | 45.2/45.1 |
| EL taxes + net SSC | 34.0 | 34.2 | 33.4 | 34.6 | 36.8 | 40.3 | 39.5 |
| DE taxes + net SSC | 39.0 | 39.4 | 39.8 | 38.4 | 38.8 | 39.5 | 39.8 |

- The 2013 Greek deficit and primary balance include one-off bank recapitalisation support, so −9.4 is not the underlying stance. The 2012 debt dip reflects the PSI.
- Government final consumption, % of GDP (WDI): EL 21.0 (2008), 23.7 (2009), 22.5 (2010), 22.1 (2011), 22.6 (2012), 21.0 (2013); DE 18.4–20.0; euro area 20.2–21.8. — [World Bank WDI NE.CON.GOVT.ZS](https://data.worldbank.org/indicator/NE.CON.GOVT.ZS?locations=GR-DE)
- Hardouvelis & Vayanos: "a primary deficit of 10.3% of GDP in 2009 turned into a primary surplus of 0.4% in 2014". — [Hardouvelis & Vayanos, LSE GreeSE 180](https://www.lse.ac.uk/asset-library/hellenic-observatory/greese-no180.pdf)
- IMF (July 2013 staff report, 2013-vintage data): gross debt path 130 / 148 / 170 / 157 / 176 / 176 / 174% of GDP, primary balance from −10.5 to +1.5. The year alignment is unverified in the snippet. — [IMF CR 13/241](https://www.imf.org/external/pubs/ft/scr/2013/cr13241.pdf)
- **Fiscal reaction:** an ECB working paper on fiscal fatigue finds the euro-area primary balance rises about 0.03–0.05 pp per 1 pp of debt/GDP, in a sample that includes Greece. — [ECB WP 2036](https://www.ecb.europa.eu/pub/pdf/scpwps/ecbwp2036.en.pdf)
- A panel-spline study finds strong evidence of debt responsiveness only when Greece is included, implying a high Greek response in pooled data. — [PMC9374487](https://pmc.ncbi.nlm.nih.gov/articles/PMC9374487)
- An NBS (Slovakia) paper finds Greece and Italy react less than the EU average. These two findings conflict. — [Hajnovič et al., NBS WP 1/2012](https://www.nbs.sk/_img/documents/_publikacie/vyskum/konferencia-14092012/zeman-hajnovic-zilinsky-fiscal_space_wp%205-2012-rev.pdf)
- Everaert & Jansen find Greek fiscal policy "sustainable" under their most general specification but stress that Bohn's condition is weak. — [NBB WP 320](https://aei.pitt.edu/97460/1/wp320en.pdf)

### Inferences
- In calibration units, a debt target near 110–150% of annual GDP (440–600% of quarterly GDP) corresponds to 2008–2010.
- Greece's primary-balance swing (−10.3 in 2009 to about 0 by 2013–14 underlying) was programme-imposed. A backward-looking Bohn coefficient estimated on 2008–2013 would be contaminated by the programme conditionality.

### Gaps
- I found no published point estimate of the Greece-specific Bohn coefficient with a standard error. The Stoica–Leonte study exists, but its coefficient was not accessible.
- The IMF 13/241 table columns were not verified year by year.

## 3. Macro: GDP collapse, expenditure shares, capital, labour share, relative size

### Takeaway
Greek real GDP fell **28.2% from peak (2007Q2) to trough (2013Q1)** on quarterly chain-linked data, and 27.0% on annual 2008→2013 data. Investment (GFCF) fell 66% in real terms over 2007–2013. In 2010, Greece was 8.5% of German nominal GDP and 2.3% of the euro area, but 13.6% of Germany's population. Its size share in a GR+DE union is therefore about 0.08 by GDP and 0.12 by population.

### Cited Findings
**Real GDP (Eurostat namq_10_gdp, chain-linked 2010 €m, seasonally and calendar adjusted; computed)** — [Eurostat namq_10_gdp](https://ec.europa.eu/eurostat/databrowser/view/namq_10_gdp/default/table)
- EL: 2007Q2 62,614 (peak), 2008Q1 62,467, 2009Q1 58,400, 2010Q1 58,921, 2011Q1 52,378, 2012Q1 47,582, 2012Q4 45,305, 2013Q1 44,980 (trough), 2013Q4 45,535.
- Changes: peak→trough −28.2%; 2008Q1→2013Q1 −28.0%; 2010Q1→2012Q4 −23.1%.
- DE: 2008Q1 670,016 → 2009Q1 622,914 (−7.0%), back above the pre-crisis peak by 2011Q1.

**Annual real aggregates (Eurostat nama_10_gdp, CLV10 €m; computed)** — [Eurostat nama_10_gdp](https://ec.europa.eu/eurostat/databrowser/view/nama_10_gdp/default/table)
- EL GDP 2008 247,276 → 2013 180,522 (−27.0%).
- EL household consumption 161,056 → 118,128 (−26.7%).
- EL GFCF 2007 60,127 → 2013 20,680 (−65.6%).
- IMF real GDP growth rows, 2013 vintage: −4.9 (2009), −7.1 (2010), −6.0 (2011), −6.4 (2012), −4.2 (2013), with GFCF growth of −15.0, −19.6, −14.4, −19.2, −4.0 (year alignment inferred). — [IMF CR 13/156](https://www.imf.org/external/pubs/ft/scr/2013/cr13156.pdf)
- IMF summary: "output contracted by 22 percent during 2008–12", with unemployment about 27%. — [IMF PIN on 2013 Article IV](https://www.imf.org/en/news/articles/2015/09/28/04/53/pn1364)
- Real GDP per capita fell 25.2% from 2007 to 2014 (€22,500 → €16,830). Investment/GDP fell from 25.7% to 11.6%, and unemployment rose from 8.4% to 26.5%. — [Hardouvelis & Vayanos, LSE GreeSE 180](https://www.lse.ac.uk/asset-library/hellenic-observatory/greese-no180.pdf)

**Expenditure shares, % of GDP (World Bank WDI)** — [WDI](https://data.worldbank.org/indicator/NE.GDI.FTOT.ZS?locations=GR-DE)

| % GDP | 2008 | 2009 | 2010 | 2011 | 2012 | 2013 |
|---|---|---|---|---|---|---|
| EL GFCF | 23.3 | 20.4 | 16.8 | 13.9 | 12.3 | 11.5 |
| DE GFCF | 20.0 | 18.8 | 19.2 | 20.0 | 20.0 | 19.6 |
| EL household consumption | 65.7 | 66.3 | 66.7 | 69.5 | 69.6 | 69.2 |
| DE household consumption | 55.0 | 57.2 | 55.8 | 54.8 | 55.2 | 54.8 |
| EL gov. consumption | 21.0 | 23.7 | 22.5 | 22.1 | 22.6 | 21.0 |

- EL gross capital formation incl. inventories (Eurostat P5G / GDP, computed): 25.7% (2007), 24.8% (2008), 18.9% (2009), 18.2% (2010), 14.1% (2011), 12.3% (2012), 12.1% (2013).

**Capital–output (Eurostat nama_10_nfa_bs, net fixed assets N11N, total economy, current prices, all assets including dwellings; computed)** — [Eurostat nama_10_nfa_bs](https://ec.europa.eu/eurostat/databrowser/view/nama_10_nfa_bs/default/table)
- EL K/Y (annual): 2.77 (2007), 2.83 (2008), 2.91 (2009), 2.98 (2010), 3.14 (2011), 3.35 (2012), 3.38 (2013). The rise reflects the fall in Y, not capital accumulation.
- DE K/Y: 3.10 (2007), 3.16 (2008), 3.23 (2010), 3.22 (2012).

**Labour share (Eurostat D1, employees and self-employed; computed)**
- Raw D1/GDP: EL 34–37%, DE 48–51%. Greece's raw share is low because about a third of employment is self-employed.
- Adjusted share, D1 × (total employment / employees) / GDP, the AMECO-style self-employed correction:
  - EL: 51.9% (2008), 54.4% (2010), 55.8% (2012).
  - DE: 55.2% (2008), 56.6% (2010), 57.2% (2012).
- On GDP net of production taxes (GDP − D2+D3): EL 58.5% / 61.3% / 63.4%; DE 61.2% / 62.6% / 63.5%.
- Sources: [Eurostat nama_10_gdp](https://ec.europa.eu/eurostat/databrowser/view/nama_10_gdp/default/table), [nama_10_pe](https://ec.europa.eu/eurostat/databrowser/view/nama_10_pe/default/table). Employment in thousands: EL 4,867 total / 3,212 employees (2008); DE 40,851 / 36,372.

**Labour market and prices (Eurostat)**
- Unemployment rate: EL 9.8 (2009), 12.9, 18.1, 24.8, 27.8 (2013); DE 7.3 → 4.9. — [une_rt_a](https://ec.europa.eu/eurostat/databrowser/view/une_rt_a/default/table)
- HICP inflation: EL 4.2 (2008), 1.3, 4.7, 3.1, 1.0, −0.9 (2013); DE 2.8, 0.2, 1.1, 2.5, 2.2, 1.6. — [prc_hicp_aind](https://ec.europa.eu/eurostat/databrowser/view/prc_hicp_aind/default/table)
- Employment: EL −11.6% over 2008→2013 (4,867k → 4,301k).

**Relative size (computed)**
- Nominal GDP, €m (Eurostat): EL 238,735 / 223,590 / 185,899 and DE 2,589,280 / 2,615,260 / 2,800,680 for 2008 / 2010 / 2012. EA20 was 9.72 / 9.64 / 9.95 trillion.
- GR/DE nominal GDP: 0.092 (2008), 0.086 (2010), 0.066 (2012).
- GR/EA20: 0.025, 0.023, 0.019.
- GR/(GR+DE): 0.084, 0.079, 0.062.
- Population (WDI), 2010: GR 11.12m, DE 81.78m, euro area 348.4m. GR/DE = 0.136; GR/(GR+DE) = 0.120; GR/EA = 0.032. — [WDI SP.POP.TOTL](https://data.worldbank.org/indicator/SP.POP.TOTL?locations=GR-DE)
- Per-capita GDP 2010, GR relative to DE: about 0.63 (0.0855 / 0.136).

### Inferences
- Under symmetric per-capita SS, the size parameter is a choice between population (≈0.12) and GDP weights (≈0.08) in a GR–DE union. Against the whole euro area, Greece is about 2–3%. Because Greek per-capita GDP is about 63% of Germany's, symmetric per-capita SS cannot match both weights at once.
- The adjusted labour shares are close in the two countries (≈0.55–0.57 of GDP), so a common α ≈ 0.33–0.45 is defensible depending on the tax/mixed-income treatment. Bocola uses α = 0.30.

### Gaps
- AMECO's official adjusted wage share (ALCD0) and net capital stock (OKND) were not pulled directly. My adjusted share is the standard approximation, not the AMECO series.
- Productive (non-residential) capital stock was not separated from dwellings.

## 4. Trade and external position

### Takeaway
Greece is a relatively closed economy for its size. Imports were about 27–32% of absorption over 2008–2012 and exports 19–29% of GDP. The current-account deficit was about 15% of GDP in 2008 and still about 10% in 2010–2011, and NIIP was about −100% of GDP in 2010, falling to −119% in 2012. Germany ran a surplus of about +5.5% of GDP with NIIP of +22–25%.

### Cited Findings
**Eurostat nama_10_gdp at current prices (computed)** — [Eurostat nama_10_gdp](https://ec.europa.eu/eurostat/databrowser/view/nama_10_gdp/default/table)

| | 2008 | 2009 | 2010 | 2011 | 2012 | 2013 |
|---|---|---|---|---|---|---|
| EL imports / absorption (P7/(P3+P5G)) | 30.8 | 25.4 | 27.2 | 29.3 | 32.0 | 31.9 |
| EL exports / GDP | 22.8 | 18.7 | 21.7 | 25.3 | 29.0 | 30.3 |
| EL imports / GDP | 34.3 | 27.6 | 29.2 | 31.0 | 33.4 | 32.6 |
| EL net exports / GDP | −11.5 | −8.9 | −7.5 | −5.7 | −4.5 | −2.3 |
| DE imports / absorption | 37.3 | 32.2 | 36.2 | 38.7 | 39.0 | 38.1 |
| DE exports / GDP | 40.9 | 35.5 | 39.5 | 41.6 | 42.7 | 41.6 |
| DE net exports / GDP | 5.8 | 4.8 | 5.1 | 4.9 | 6.0 | 5.8 |

- The World Bank WDI agrees: GR imports 29.19% of GDP in 2010, trade 50.9% of GDP. — [WITS Greece 2010](https://wits.worldbank.org/CountryProfile/en/Country/GRC/Year/2010)
- **Current account, % of GDP (WDI BN.CAB.XOKA.GD.ZS):** EL −14.2 (2007), −14.6 (2008), −11.0 (2009), −10.2 (2010), −10.1 (2011), −2.6 (2012), −2.1 (2013); DE +6.5, +5.4, +5.7, +5.6, +6.1, +7.0, +6.6. — [WDI](https://data.worldbank.org/indicator/BN.CAB.XOKA.GD.ZS?locations=GR-DE)
- Hardouvelis & Vayanos give an older-vintage EL CA deficit of 15.9% of GDP in 2007 (versus 5.3% in 1998). — [LSE GreeSE 180](https://www.lse.ac.uk/asset-library/hellenic-observatory/greese-no180.pdf)
- The OECD survey CA row of −9.9 / −3.4 / −0.4 / 1.3 refers to later or projection years. — [OECD Economic Survey Greece 2013](https://www.oecd.org/content/dam/oecd/en/publications/reports/2013/11/oecd-economic-surveys-greece-2013_g1g38894/eco_surveys-grc-2013-en.pdf)
- **NIIP, % of GDP (Eurostat tipsii10):** EL −95.1 (2007), −76.8 (2008), −88.7 (2009), −100.0 (2010), −90.4 (2011), −119.2 (2012), −132.3 (2013); DE +18.5, +17.9, +24.6, +25.3, +22.8, +22.0, +24.3. — [Eurostat tipsii10](https://ec.europa.eu/eurostat/databrowser/view/tipsii10/default/table)
- Hardouvelis & Vayanos quote EL NIIP of −19.3% (1998) and −99.9% (2007), an older vintage. — [LSE GreeSE 180](https://www.lse.ac.uk/asset-library/hellenic-observatory/greese-no180.pdf)
- Greek exports are concentrated in oil products, shipping and tourism, and more than 55% go to the EU. — [Piraeus Bank, "A Closer Look to External Trade"](https://www.piraeusholdings.gr/-/jssmedia/Com/2016/Files/oikonomiki-analisi-ependitiki-stratigiki/oikonomiki-analisi-ellinikis-oikonomias-kladikes-meletes/A-Closer-Look-to-External-Trade-Greeces-Exports-Performance.pdf)

### Inferences
- Measured against all partners, home bias in Greek absorption is about 0.70–0.73, and about 0.62–0.64 for Germany. In a two-country model where F is "Germany / rest of the union", the import share of absorption is roughly 0.27–0.31. A strict bilateral GR–DE import share would be far smaller, since Germany supplies only about one-tenth of Greek imports (not verified here).
- The 2008–2012 current-account adjustment, from −14.6 to −2.6, came mostly through import compression, with imports falling with absorption. The export ratio rose mostly because GDP fell.

### Gaps
- No Greece-specific trade (Armington) elasticity estimate was found in this pass. Macro models conventionally use about 1–1.5, and the trade literature about 4–6, but neither is sourced here.
- Bilateral GR–DE import shares (Comtrade/OECD TiVA) were not pulled.
- The Lane–Milesi-Ferretti External Wealth of Nations data were not accessed directly. The Eurostat NIIP is used instead.

## 5. Household heterogeneity: wealth, deposits, income risk, consumption by group

### Takeaway
In HFCS wave 1 (Greek fieldwork in 2009), median Greek household net wealth was €101.9k and the mean €147.8k, while Germany's was €51.4k median and €195.2k mean. Only 73.4% of Greek households held a deposit account (99.0% in Germany), and the conditional median deposit was €3.6k. Deposits are 80.7% of Greek household financial assets. Median net liquid assets were just 4.9% of income (22.3% in Germany), which points to a large hand-to-mouth share. Household deposits were about 85–93% of GDP over 2008–2012, and domestic deposits at Greek banks fell 30% from end-2009 to April 2012.

### Cited Findings
**ECB HFCS first wave, Statistics Paper No. 2 (April 2013). Greek reference year 2009; German 2010.** — [ECB Statistics Paper 2](https://www.ecb.europa.eu/pub/pdf/scpsps/ecbsp2.en.pdf)

| Item (table) | Greece (2009) | Germany (2010) | Euro area |
|---|---|---|---|
| Median net wealth, €k (T4.1) | 101.9 (s.e. 2.5) | 51.4 (3.2) | 109.2 |
| Mean net wealth, €k (T4.1) | 147.8 (5.0) | 195.2 (11.9) | 230.8 |
| Share of EA net wealth / of EA households, % (T4.1) | 1.9 / 3.0 | 24.3 / 28.7 | 100 / 100 |
| Participation: any financial asset / deposits, % (T2.4) | 74.5 / 73.4 | 99.3 / 99.0 | 96.8 / 96.4 |
| Participation: mutual funds / bonds / shares, % (T2.4) | 1.2 / 0.5 / 2.7 | 16.9 / 5.2 / 10.6 | 11.4 / 5.3 / 10.1 |
| Conditional median: financial assets / deposits, €k (T2.5) | 4.4 / 3.6 | 17.1 / 7.9 | 11.4 / 6.1 |
| Deposits as share of financial assets, % (T2.6) | 80.7 | 44.4 | 42.9 |
| Participation: any debt / mortgage / non-mortgage, % (T3.1) | 36.6 / 17.5 / 26.1 | 47.4 / 21.5 / 34.6 | — |
| Median debt-to-asset / debt-to-income / debt-service-to-income, % (T3.4) | 14.8 / 47.2 / 14.7 | 28.4 / 37.3 / 10.9 | 21.8 / 62.0 / 13.9 |
| **Median net liquid assets / income, % (T3.4)** | **4.9** | 22.3 | 18.6 |
| Household size / members 16+ / employed (T1.1) | 2.64 / 2.22 / 1.05 | 2.04 / 1.72 / 0.95 | 2.32 / 1.93 / 0.97 |

- Greek participation (T2.1): any real asset 92.2%, HMR ownership 72.4%, other real estate 37.9%, self-employment business 9.8%. 18.9% of Greek household reference persons are self-employed (euro area 9.0%). — [ECB Statistics Paper 2](https://www.ecb.europa.eu/pub/pdf/scpsps/ecbsp2.en.pdf)
- Greek wave 1 was fielded in 2009 with 2,971 households, and wave 2 in 2014. The change in median net wealth between the waves is significant at 1%. — [Bank of Greece paper on HFCS](https://www.bankofgreece.gr/RelatedDocuments/4520170702_en.pdf)
- **Household deposits (Eurostat nasa_10_f_bs, S14_S15, F2 assets; computed /GDP):** EL €201.9bn (2008, 84.6% of GDP), €213.1bn (2009, 90.9%), €201.5bn (2010, 90.1%), €187.3bn (2011, 92.1%), €173.4bn (2012, 93.3%), €170.4bn (2013, 95.7%). DE: 64.5% (2008) to 68.8% (2012). — [Eurostat nasa_10_f_bs](https://ec.europa.eu/eurostat/databrowser/view/nasa_10_f_bs/default/table)
- Domestic deposits at Greek banks fell from €237.5bn (end-2009) to €166bn (April 2012), −30.1%. Their share of bank liabilities fell from 58.1% to 49.8% (end-2012), while ECB funding rose from 10.3% to 31.4% of liabilities. — [Hardouvelis & Vayanos, LSE GreeSE 180](https://www.lse.ac.uk/asset-library/hellenic-observatory/greese-no180.pdf)
- **Household gross saving rate (Eurostat nasa_10_ki):** EL 9.2% (2008), 11.1 (2009), 5.0 (2010), 1.0 (2011), −3.4 (2012), −5.3 (2013); DE about 16.4–17.2 throughout. — [Eurostat nasa_10_ki](https://ec.europa.eu/eurostat/databrowser/view/nasa_10_ki/default/table)
- **Consumption by group:** household spending fell almost 30% between 2008 and 2013 on the Household Expenditure Survey, and consumption inequality rose on all measures. — [Kaplanoglou, "Consumption inequality and poverty in Greece"](https://eclass.uoa.gr/modules/document/file.php/ECON208/Some%20Interesting%20papers%20(not%20for%20the%20exams)/2022%20Kaplanoglou_Inequality.pdf)
- ELSTAT HBS average monthly household expenditure fell from €2,401.44 (2008) to €1,637.10 (2012), −22.7% nominal and −31.8% at constant prices. — [ResearchGate, Kaplanoglou & Rapanos](https://www.researchgate.net/publication/311627332_Evolutions_in_Consumption_Inequality_and_Poverty_in_Greece_The_Impact_of_the_Crisis_and_Austerity_Policies); [ELSTAT HBS](https://www.statistics.gr/documents/20181/985029/Results+of+the+Household+Budget+Survey+(HBS)/2a63d7b8-666a-4daf-816c-6afcc38c215a?version=1.0)
- The bottom 20% roughly kept its expenditure share while the middle lost in both absolute and relative terms. — [Kaplanoglou](https://eclass.uoa.gr/modules/document/file.php/ECON208/Some%20Interesting%20papers%20(not%20for%20the%20exams)/2022%20Kaplanoglou_Inequality.pdf); see also [Bank of Greece WP 174](https://www.bankofgreece.gr/Publications/Paper2014174.pdf)

### Inferences
- Greek wealth is housing-heavy and liquid wealth is very thin. A median net-liquid-assets-to-income ratio of 4.9% and 26.6% of households with no deposit account imply a hand-to-mouth share well above Germany's. A conditional median deposit of €3.6k means a substantial fraction of households holds almost no liquid buffer.
- Aggregate household deposits (about 90% of GDP) are concentrated at the top, since the median household's deposits are small. A model matching deposits/GDP and the HFCS median at the same time needs a skewed liquid-wealth distribution.

### Gaps
- No Kaplan–Violante–Weidner-style hand-to-mouth share for Greece was found. KVW (2014, BPEA) cover DE, FR, IT, ES but not GR.
- No Greek earnings-risk estimate (persistent/transitory variances from EU-SILC panels) was found in this pass.
- Greek median gross household income from HFCS Table 5.1 was not extracted.

## 6. Literature calibration values

### Takeaway
Bocola (2016) is quarterly and Italian:
- Leverage 5, sovereign exposure 7.6% of bank assets, average bond maturity 18 months (fraction maturing 0.056 per quarter), haircut D = 0.55 (from ZTG).
- Logit default process with s* = −7.06 (about 0.09% per quarter), ρ_s = 0.95, σ_s = 0.63.
- Frisch elasticity 2, capital share 0.30.
- Interbank-spread sample average of 8 bp (annualised).

Costain–Nuño–Thomas set the haircut at δ = 0.25, steady-state default intensity ψ_ss = 14 bp per year, and risk aversion γ = 0.16. Neumeyer–Perri set the working-capital share θ = 1 (the entire wage bill paid in advance).

### Cited Findings
**Bocola (2016, JPE 124(4), 879–926), Table 1, parameters set externally (quarterly, Italy)** — [Bocola 2016 (MFM copy)](https://mfm.uchicago.edu/wp-content/uploads/2020/06/Bocola-The-Pass-Through-of-Sovereign-Risk.pdf)
- exp_bg (bond holdings / bank total assets) = 0.076
- π (fraction of bonds maturing per quarter) = 0.056 (average maturity 18 months, EBA 2011 stress test plus Bankscope, end-2010, five largest Italian banks)
- exp{g*} (spending/output) = 0.198; ρ_g = 0.894; σ_g = 0.012 (AR(1) over 1999Q1–2009Q4)
- i/y = 0.213; hours l = 0.318; R_bg = 1.003 (gross quarterly)
- Capital share α = 0.300; inverse Frisch ν = 0.500 (Frisch = 2)
- Leverage = 5.000 ("in the range of values used in the literature")
- Fiscal rule γ_τ = 1; adjustment cost at BGP = 0; bond price normalised to 1

**Bocola (2016), Table 2, posteriors**
- μ_bg×100 = 0.10 [0.02, 0.21], i.e. a Lagrange multiplier of about 0.001, roughly one-tenth of the Gertler–Karadi value.
- ω = 0.96 [0.94, 0.98] (banker survival; prior centred on GK)
- ψ = 0.42 [0.31, 0.52] (prior mean 0.5, "conventional")
- g×400 = 0.43 [0.17, 0.72]
- ρ_z = 0.51; σ_z×100 = 0.59
- **s* = −7.06 [−9.73, −4.31]; ρ_s = 0.95 [0.88, 0.99]; σ_s = 0.63 [0.52, 0.76].** The minus signs are lost in the PDF text extraction; −7.06 is consistent with the 0.09% quarterly default probability quoted in the text.

**Bocola, default process and measurement**
- pd = e^s/(1+e^s), with s′ = (1−ρ_s) log(s*) + ρ_s s + σ_s ε. This is eq. 12 as printed; the "log" is Bocola's notation.
- The default probability is backed out of **1-year** Italian CDS (Markit, CR clause, quarterly premium, in bp) using a risk-adjusted formula (eq. 19), with **haircut D = 0.55**, "consistent with the Greek debt restructuring of 2012 (Zettelmeyer, Trebesch, and Gulati 2013)".
- Table 3 data moments of pd (% per quarter): median 0.02, mean 0.19, s.d. 0.36, autocorrelation 0.82, skewness 3.02, kurtosis 12.74.
- Table 4 conditioning values: "low s" = 0.09% and "high s" = 2.5% default probability next quarter. Risk premia on firms' claims rise by about 30 bp annualised between them.
- The sample average interbank spread is **8 basis points** (annualised), and the model-implied average liquidity premium is 36 bp annualised.
- Data: 2002Q2–2012Q4 (multiplier and GDP growth); monthly rates 2002:M4–2012:M12.

**Costain, Nuño & Thomas** — [CNT, Banco de España record, published J. Finance 80(4), Aug 2025, pp. 2389–2434, DOI 10.1111/jofi.13463](https://www.bde.es/wbe/en/publicaciones/analisis-economico-investigacion/documentos-trabajo/term-structure-interest-rates-heterogeneous-monetary-union.html); values read from a locally held November 2021 draft. Later drafts: [BdE WP 2223](https://repositorio.bde.es/bitstream/123456789/22545/1/dt2223e.pdf), [BIS WP 1165](https://www.bis.org/publ/work1165.pdf)
- Two-country preferred-habitat/arbitrageur term-structure model, calibrated to Germany (core) and Italy (periphery) over 2013–2019.
- Short rate: r̄ = −49 bp, σ = 32 bp, mean reversion κ = 0.12 (monthly autocorrelation 0.99). Maturities run from 1 month to 10 years.
- Preferred-habitat elasticity α = 500. Arbitrageurs hold 56% of net debt.
- Arbitrageur risk aversion γ = 0.16, fitted to the German 10-year term premium.
- Haircut δ = 0.25 ("reasonable … given the international evidence", Cruces–Trebesch 2013).
- Steady-state default intensity ψ_ss = 14 bp per year, an expected default loss δψ_ss under 4 bp per year. Default sensitivity to fiscal pressure θ = 2.02×10⁻⁵ (fitted to the PEPP announcement).
- Rollover-crisis discounting: r̂ + φ = 0.5 (expected crisis duration 2 years).
- Robustness: γ = 0.078 with ψ_ss = 32 bp and θ = 4.32×10⁻⁵ (longer sample); γ = 0.242 with ψ_ss = 13 bp and θ = 9.50×10⁻⁵ (α = 2000).
- Caveat: these values come from the 2021 working paper held locally. The 2025 JF version may differ.

**Neumeyer & Perri (2005, JME 52, 345–380)**, quarterly, Argentina 1983–2001 — [Neumeyer & Perri, JME 2005 (IGC copy)](http://www.theigc.org/sites/default/files/2016/06/Neumeyer_Perri_JME.pdf)
- "We assume that all the wage bill is paid in advance and set θ = 1." With θ = 0 the model reduces to the standard neoclassical one.
- Labour share = (1−α)/(1+(R̄−1)θ). They target a labour share of 0.6.
- GHH utility curvature σ = 5 (Reinhart–Végh). Labour curvature ν = 1.6, so the labour-supply elasticity is 1/(ν−1) ≈ 1.67.
- β = 0.93 (GHH) or 0.98 (Cobb–Douglas). Real rate 14.8% per year. Hours 20%. I/Y = 0.21. NFA/Y = −42%.
- Table 2 shock processes: TFP ρ_A = 0.95; international rate ρ₁ = 0.81, σ = 0.63%; independent country risk ρ₂ = 0.78, σ_D = 2.59%; induced country risk η̄ = 1.04, σ_I = 1.7%.
- A Turkish replication finds results "critically depend on the magnitude of the working capital parameter", with country-spread shocks explaining under 9% of output volatility. — [BCB WP 218 / related](https://www.bcb.gov.br/pec/wps/ingl/wps218.pdf)
- A Brazilian paper restates the 100% wage-bill-in-advance assumption. — [BCB WP 218](https://www.bcb.gov.br/pec/wps/ingl/wps218.pdf)

**Typical preference values**
- Bocola: Frisch 2. Neumeyer–Perri: GHH labour elasticity about 1.67 and risk aversion 5.
- Risk aversion of 2 is "standard in the sovereign default literature". — [Latin American Economic Review calibration table citing Arellano](https://latinaer.springeropen.com/articles/10.1186/s40503-018-0060-2/tables/1); [arXiv 1512.06960](https://arxiv.org/pdf/1512.06960)

### Inferences
- Bocola's leverage of 5 and exposure of 7.6% are Italian. Greek banks' domestic-sovereign exposure was much larger relative to assets, and Hardouvelis–Vayanos document that the PSI made all Greek banks insolvent. Bocola's values are therefore a lower bound for D-country exposure; the parallel bank-calibration audit should hold the Greek EBA numbers.
- In units: Bocola's 8 bp is an **annualised** interbank-minus-risk-free spread, used to measure the IC multiplier. It is not a lending spread. His π = 0.056 per quarter implies an average maturity of 1/0.056 ≈ 18 quarters as a perpetuity-decay parameter, but the text states 18 months. The two readings conflict and should be checked against the model equation before porting.
- CNT's δ = 0.25 (cross-country average from Cruces–Trebesch) and Bocola's D = 0.55 (Greek PSI) bracket the haircut choice. For a Greek-calibrated D country, 0.535–0.65 is the empirically grounded range.

### Gaps
- **Arellano (2008, AER)** parameters were not verified from the primary source in this pass. My recollection is β = 0.953, σ = 2, r = 1.7% quarterly, re-entry probability 0.282, output-cost threshold ŷ = 0.969·E[y]; check against AER Table 1 before use. A secondary source cites a different set (β = 0.81, re-entry 0.10, likely annual), and these conflict.
- The CNT 2025 JF published table was not accessed. WP 2021 values are reported.
- No Greece-specific working-capital share estimate was found. Neumeyer–Perri's θ = 1 is an assumption, not an estimate. Mendoza (2010) and Uribe–Yue (2006) estimates were not retrieved.
- No source was found for the Greek bank lending spread or interbank spread in 2008–2013 (needed for a Greek analogue of Bocola's 8 bp).
