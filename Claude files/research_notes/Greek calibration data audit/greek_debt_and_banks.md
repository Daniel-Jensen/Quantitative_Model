# Greek sovereign debt and Greek/German banking-system data, 2009-2012 (calibration fact sheet)

Scope: Greece (D) and Germany (F), 2009-2012. Units are EUR unless noted. "bn" means EUR billion.
The repository file `data/eba_moments.json` was read (not modified) and is cross-checked below. In this file, "EBA file" means the numbers in that file: EBA 2011 stress test, base date 31-Dec-2010, 6 Greek banks (GR030-035) and 12 German banks.

## Q1. Greek general government gross debt and nominal GDP, 2009-2012

### Takeaway
Under the current Eurostat vintage (pulled 2026-10-23), Maastricht debt was 301.1 / 330.6 / 356.2 / 305.1 bn in 2009-2012. That is 128.5 / 147.8 / 175.1 / 164.1 % of GDP. Nominal GDP fell from 234.3 bn (2009) to 185.9 bn (2012). The 2012 fall in debt is the March 2012 PSI and the December 2012 buyback, net of new official loans.

### Cited Findings
Source for every number in this list: [Eurostat gov_10dd_edpt1 and nama_10_gdp via the dissemination API, last updated 2026-10-21 and 2026-10-22](https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/gov_10dd_edpt1?geo=EL&sector=S13&na_item=GD&unit=MIO_EUR). The series are: S13 general government; GD = consolidated gross debt at nominal (face) value; B1GQ = GDP at current prices.

| | 2009 | 2010 | 2011 | 2012 |
|---|---|---|---|---|
| EL gross debt, EUR m | 301,062 | 330,570 | 356,235 | 305,085 |
| EL gross debt, % GDP | 128.5 | 147.8 | 175.1 | 164.1 |
| EL nominal GDP, EUR m | 234,322.6 | 223,590.5 | 203,475.0 | 185,899.2 |
| DE gross debt, EUR m | 1,804,696 | 2,119,241.8 | 2,155,875.7 | 2,234,802.5 |
| DE gross debt, % GDP | 72.4 | 81.0 | 78.5 | 79.8 |
| DE nominal GDP, EUR m | 2,493,970 | 2,615,260 | 2,747,210 | 2,800,680 |

- The 2010 GDP figures (EL 223,590.5 m, DE 2,615,260 m) are exactly the `GDP_ann_D/F` in the EBA file. The repository and the current Eurostat vintage therefore agree. — [Eurostat nama_10_gdp](https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/nama_10_gdp?geo=EL&geo=DE&na_item=B1GQ&unit=CP_MEUR)
- Vintage risk: Eurostat's April 2010 notification carried a reservation on the Greek data and warned that the 2009 debt ratio could be revised by 5-7 pp of GDP. The 2009 deficit was later revised from 13.6% to 15.4% of GDP (Nov 2010). — [Eurostat April 2010 release (mirror)](https://www.nuevatribuna.es/media/nuevatribuna/editmaker/pdf/eurostat.pdf); [MNI/investinglive, 15 Nov 2010](https://investinglive.com/news/!/eurostat-revises-greece-2009-deficit-to-15-4-of-gdp-vs-13-6-20101115/)
- Over 2010-2011, total assets of the Greek bank sector (MFIs ex-BoG) were roughly 2.9x GDP (654.6 bn at Dec-2010 against 223.6 bn GDP). — [IMF CR 11/351, Table 2](https://www.imf.org/external/pubs/ft/scr/2011/cr11351.pdf) combined with Eurostat

### Inferences
- Debt-to-GDP rose 2010 to 2011 mostly because GDP fell 9.0% in nominal terms, not only because debt grew 7.8%.
- Quarterly nominal GDP for the model (annual/4) is 55.9 bn in 2010 and 50.9 bn in 2011.
- Debt ratios published at the time (ESA95, pre-2014 GDP benchmark) differ from the current vintage. The EUR amounts are essentially unchanged; the ratios move with GDP revisions. The ratio vintage should be stated whenever a % of GDP target is quoted.

### Gaps
- The exact ratios as first published in the April 2012 and April 2013 EDP notifications (ESA95 vintage) were not retrieved. The search did not surface the Greek tables, so I did not quote them.

## Q2. Debt composition (securities vs loans) and the March 2012 PSI

### Takeaway
At end-2011, 70.6% of Greek debt was securities (251.6 bn) and 29.1% was loans (103.5 bn, mostly GLF and IMF). By end-2012 this had flipped to 30.6% securities and 69.0% loans. The PSI targeted 205.6 bn of privately held bonds. 199.2 bn (96.9%) was exchanged at a 53.5% nominal haircut, and the present-value haircut was 59-65%. The PV recovery was about 21-23 cents per euro of face, or 35-40% of the pre-exchange present value.

### Cited Findings
Composition. Source for this table: [Eurostat gov_10dd_edpt1, EL, S13](https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/gov_10dd_edpt1?geo=EL&sector=S13&na_item=GD_F3&na_item=GD_F4&unit=MIO_EUR). Values are EUR m, nominal value.

| | 2009 | 2010 | 2011 | 2012 |
|---|---|---|---|---|
| F2 currency & deposits | 1,663 | 1,203 | 1,106 | 1,055 |
| F3 debt securities | 262,773 | 262,555 | 251,622 | 93,432 |
| - F31 short-term (T-bills) | 10,024 | 8,273 | 10,850 | 15,309 |
| - F32 long-term (bonds) | 252,749 | 254,282 | 240,772 | 78,123 |
| F4 loans | 36,626 | 66,812 | 103,507 | 210,598 |
| - F42 long-term loans | 33,213 | 63,643 | 100,589 | 208,107 |

Official lending and programme facts:
- First programme, May 2010: 80 bn of EU bilateral loans (GLF) plus 30 bn of IMF credit. Only about 65 bn of the 110 bn had been disbursed when the July 2011 plan was drawn up. — [Zettelmeyer, Trebesch & Gulati (ZTG) 2013, CESifo WP 4333, Sec. 2](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
- PSI-related official lending: Greece borrowed 29.7 bn from the EFSF to pay the short-term EFSF-note "sweetener". In April 2012 it borrowed another 25 bn from the EFSF for bank recapitalisation. — [ZTG 2013, Sec. 3.3 and fn. 34](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)

PSI eligible stock and timeline:
- The offer covered all privately held sovereign bonds issued before 2012 (face 195.7 bn) plus 36 sovereign-guaranteed bonds (face just under 10 bn). The total is 205.6 bn. — [ZTG 2013, p. 10-11](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf); [Hardouvelis & Vayanos 2023, LSE GreeSE 180](https://www.lse.ac.uk/Hellenic-Observatory/Assets/Documents/Publications/GreeSE-Papers/GreeSE-No180.pdf) also gives "from €205.6bn to €98.5bn".
- 177.3 bn of the sovereign bonds were Greek-law. The rest were mostly English-law. One LSE blog puts Greek-law bonds at 93.1% of the total. That share was not checked against ZTG, and the base may differ. — [ZTG 2013](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf); [LSE Greece@LSE blog, 23 Jan 2012](https://blogs.lse.ac.uk/greeceatlse/2012/01/23/the-psis-enigma-and-a-possible-solution/)
- Timeline (all from [ZTG 2013](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)):
  - 21 Feb 2012: deal announced.
  - 23 Feb 2012: Greek Bondholder Act retrofitted collective action clauses (CACs) to Greek-law bonds.
  - 24 Feb 2012: formal offer.
  - 9 Mar 2012: 82.5% of the 177.3 bn Greek-law bonds accepted, and the CACs were activated.
  - 26 Apr 2012: last foreign-law bonds settled.
- Final participation was 199.2 bn, or 96.9% of eligible principal. Holdouts held 6.4 bn of face value across 25 bonds. — [ZTG 2013](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
- The consideration per 100 of old face had four parts. The new face totals 46.5, which is the 53.5% nominal haircut. — [ZTG 2013, p. 11](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
  - 15 in 1-2 year EFSF notes.
  - 31.5 in 20 new English-law GGBs maturing 2023-2042, with coupons of 2-4.3%.
  - A GDP warrant.
  - Accrued interest paid in 6-month EFSF notes.
- Payout: 29.7 bn of EFSF notes plus 62.4 bn of new bonds. Face value fell by about 107 bn, or 52% of eligible debt. — [ZTG 2013](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf); [Hardouvelis & Vayanos 2023](https://www.lse.ac.uk/Hellenic-Observatory/Assets/Documents/Publications/GreeSE-Papers/GreeSE-No180.pdf)
- The principal-weighted average exit yield on the new bonds was 15.3%. — [ZTG 2013](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
- [ZTG 2013, Table 3](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf). All values are in % of old principal:

| Discount rate | PVnew | Market-convention haircut (100-PVnew) | PVold | PV haircut |
|---|---|---|---|---|
| 15.3% | 23.1 | 76.9 | 65.3 | 64.6 |
| 16.3% | 22.5 | 77.5 | 63.3 | 64.4 |
| 18.7% | 21.2 | 78.8 | 59.0 | 64.0 |
| Imputed yield curve | 22.8 | 77.2 | 56.5 | 59.6 (sensitivity range 55-61) |

  ZTG's preferred estimate is 64.6%.
- The earlier July 2011 IIF offer implied a PV haircut of only 11.5% at a 9% discount rate. The range across discount rates was -2.3 to 21.4%. — [ZTG 2013, Table 1](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
- PV debt relief net of bank recapitalisation costs was 98-106 bn, or about 51-55% of GDP. — [ZTG 2013, Sec. 3.3](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
- December 2012 buyback, 12 Dec 2012: 11.3 bn of EFSF financing retired 31.9 bn of the new GGBs at an average price of about 34 cents. Face debt fell by a net 20.6 bn. Greek banks were among the sellers. — [ZTG 2013, Sec. 4.1](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf); [Hardouvelis & Vayanos 2023, fn. 14](https://www.lse.ac.uk/Hellenic-Observatory/Assets/Documents/Publications/GreeSE-Papers/GreeSE-No180.pdf)

### Inferences
- Recovery-rate choices for the model:
  - (a) Nominal recovery: 46.5% of face.
  - (b) Market-value recovery per unit of face: about 21-23% (PVnew).
  - (c) Recovery relative to the no-default value of the same claim: 1 - PV haircut, so 35.4% (at 15.3%) to 40.4% (yield-curve case).
  - The repository's `recovery_rate = 0.30` (EBA file, "from the PSI NPV literature") lies between (b) and (c). Which concept the Hatchondo-Martinez pure-haircut branch needs depends on whether the haircut applies to face or to the market price. That is a modelling question to raise with the user, not something to settle silently.
- The Eurostat F3 fall from 251.6 bn to 93.4 bn over 2011-2012 is consistent with about 199 bn exchanged into 62.4 bn of new bonds, less 31.9 bn bought back, plus new T-bills. EFSF notes and loans moved into F4.

### Gaps
- A bond-by-bond maturity profile of the 205.6 bn eligible stock was not extracted. ZTG's Figure 5 (bond-by-bond haircuts by remaining duration) and the PDMA invitation memorandum would supply it.

## Q3. Who held GGBs before the PSI

### Takeaway
The best-documented blocks before the PSI are:
- Eurosystem: ECB SMP 42.7 bn (16.3% of bonds, Feb 2012) and national central banks (NCBs) 13.5 bn. These were swapped out and exempted, together with 0.3 bn at the EIB, for 56.7 bn in all.
- Greek banks: about 45 bn nominal (IMF, late 2011).
- Greek social security funds: 23.7 bn.
- Foreign private holders: inferred as the large remainder of the 205.6 bn eligible.

German banks held about 8 bn (EBA, Dec-2010, 12 banks) to 10-17 bn (Bundesbank-based estimates, broader German sector). I found no single official holder table.

### Cited Findings
- ECB (SMP) held 42.7 bn of Greek bonds in Feb 2012, 16.3% of holdings and the largest single bondholder. NCBs held 13.5 bn (about 5%) and the EIB 0.315 bn. All were swapped into new series with identical terms and excluded from the PSI. — [ZTG 2013, p. 11](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
- Debts to non-participating creditors were 6.4 bn to holdouts and 56.7 bn to the ECB and NCBs. — [ZTG 2013](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
- The RBA gives ECB face value of about 43 bn and NCBs 13.5 bn. — [RBA Statement on Monetary Policy, May 2012, Box B](https://www.rba.gov.au/publications/smp/2012/may/pdf/box-b.pdf)
- SMP disclosure of 21 Feb 2013, as of 31 Dec 2012: Greece 33.9 bn nominal, 30.8 bn book value, average remaining maturity 3.6 years. The Eurosystem SMP total was 218.0 bn nominal (208.7 bn book), with average remaining maturity 4.3 years. — [ECB press release 21 Feb 2013](https://www.ecb.europa.eu/press/pr/date/2013/html/pr130221_1.en.html)
- Greek banks held GGBs of about 45 bn nominal (39 bn book value after the June 2011 impairments), against aggregated core capital of 22 bn. The six largest banks held 97% of the domestic GGB holdings. — [IMF Country Report 11/351 (5th review, Dec 2011), Box 3](https://www.imf.org/external/pubs/ft/scr/2011/cr11351.pdf)
- Greek social security funds held about 23.7 bn of government bonds. — [IMF CR 11/351, para. on social security](https://www.imf.org/external/pubs/ft/scr/2011/cr11351.pdf)
- Estimated holdings of creditor-committee members, June 2011, from Barclays (2011) and the IIF. — [ZTG 2013, Table 2](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
  - Greek: NBG 13.7, Piraeus 9.4, Alpha Eurobank 3.7, Marfin 2.3.
  - French: BNP Paribas 5.0, Société Générale 2.9, CNP 2.0, Groupama 2.0, Axa 1.9, BPCE 1.2, Crédit Agricole 0.6.
  - German: Commerzbank 2.9, Deutsche Bank 1.6, LBBW 1.4, Allianz 1.3.
  - Other: Dexia 3.5, Generali 3.0, Bank of Cyprus 1.8, Ageas 1.2, RBS 1.1, Unicredit 0.9, Intesa 0.8, HSBC 0.8.
  - The 32-member committee held an estimated 30-40% of privately held debt.
- German institutions cut their Greek public bond holdings from 16.6 bn to about 10 bn within a year (to mid-2011). This is Bundesbank data as reported by OSW, a secondary source. — [OSW, 15 Jun 2011](https://www.osw.waw.pl/en/publikacje/analyses/2011-06-15/german-financial-institutions-are-withdrawing-greek-bonds)
- German banks' foreign claims on Greece on a consolidated ultimate-risk basis were about 32 bn in March 2010 and 3.9 bn at end-2012. The unit (USD in BIS) and the original BIS table were not verified. — [World Economic Forum, Jan 2015](https://www.weforum.org/stories/geographies-in-depth/who-is-invested-in-greece/)
- At end-2010, French banks held USD 56 bn and German banks USD 40 bn of total foreign claims on Greece, all sectors. — [BIS QR coverage via Calculated Risk](https://www.calculatedriskblog.com/2010/06/bis-reports-bank-exposure-to-euro-area.html); secondary.
- EBA file, Dec-2010, gross direct long exposure (item 34010) to the Greek sovereign:
  - 12 German banks: 7,933.6 m (`b_D_F`).
  - 6 Greek banks: 54,447.1 m (`b_D_D`).
  - Greek banks' exposure to Germany: 410.7 m (`b_F_D`).
  - Source: `data/eba_moments.json`, raw_EURm.
- Bruegel's argument is that, from autumn 2011, non-residents cut their exposure to Southern European sovereigns while domestic banks raised theirs. In Germany the non-resident share rose. — [Merler & Pisani-Ferry 2012, Bruegel PC 2012/02 (landing page)](https://ideas.repec.org/p/bre/polcon/695.html)

### Inferences
- Rough decomposition of the 205.6 bn eligible: about 45 bn Greek banks (plus their foreign and insurance subsidiaries), about 23.7 bn social security funds, and about 135 bn other. The "other" block is mostly foreign banks and insurers, plus Greek non-bank private holders. Adding the 56.7 bn Eurosystem/EIB block gives about 262 bn of bonds, which matches the Eurostat F3 figure of 251.6 bn at end-2011 plus the guaranteed bonds. This is an inference from mixed dates; treat it as a sanity check, not data.
- On disagreements over German exposure:
  - The EBA 7.9 bn figure covers only 12 stress-tested banking groups at Dec-2010.
  - The OSW/Bundesbank 16.6 to 10 bn figure covers all German financial institutions, probably including FMS Wertmanagement/HRE and insurers.
  - Use EBA for the F-bank `phi_bD_F` (7.9 / 114.3 = 6.9% of CT1). Do not use the sector-wide figure.
- The drop in SMP holdings from 42.7 bn (Feb 2012, ZTG) to 33.9 bn (Dec 2012, ECB) is consistent with redemptions of SMP-held bonds during 2012. ZTG note that Greece repaid the ECB in August 2012 with T-bills sold to Greek banks. The 42.7 bn may also include holdings outside the SMP. Not verified.
- Face vs book: SMP Greek bonds were bought below par. The Dec-2012 book/nominal ratio was 30.8/33.9 = 0.909.

### Gaps
- Bruegel sovereign-holdings dataset numbers for Greece (domestic banks' share, non-resident share by quarter): the PDF fetch failed (connection refused). The dataset is at bruegel.org, "Sovereign bond holdings in the euro area".
- Bank of Greece / PDMA official holder breakdown of GGBs for 2010-2011: not found in this session.
- BIS CBS primary tables (immediate vs ultimate risk, public-sector claims on Greece by German and French banks): not retrieved. Only secondary reports were found.

## Q4. Greek banks: assets, capital, sovereign exposure, PSI losses, recapitalisation

### Takeaway
Greek banks' own-sovereign exposure was about 2.0-2.4x core capital:
- IMF: 45 bn GGBs against 22 bn core capital (late 2011).
- EBA: 54.4 bn gross sovereign exposure against 22.8 bn CT1 (Dec-2010, 6 banks). That is 14.4% of the 6 banks' total assets.

The PSI loss on state exposures was 37.7 bn gross (31.9 bn net of provisions), or 1.7x (1.4x net) the 22.1 bn CT1 of Dec-2011. Capital was wiped out. The Bank of Greece set total capital needs at 40.5 bn (27.5 bn for the four core banks). The HFSF put 25.5 bn of public money into the four core banks, completed in July 2013, plus 3.1 bn of private capital. In total, 38.9 bn of the 50 bn facility was used for recapitalisation and resolution.

### Cited Findings
- EBA file, Dec-2010, 6 Greek banks:
  - CT1 22,778 m; RWA 222,466 m; total assets 377,200 m.
  - Own-sovereign gross direct long 54,447 m.
  - Corporate EAD 90,304 m; CRE EAD 19,034 m.
  - Fair-valued share of the own-sovereign book: 11.1%.
  - Source: `data/eba_moments.json`.
- Ratios computed from the EBA file:
  - CT1/RWA = 10.2%.
  - Own-sovereign/CT1 = 2.39 (`phi_bD_D_ss`).
  - Own-sovereign/TA = 14.4%.
  - TA/CT1 = 16.56.
- NBG's own disclosure: aggregate direct and indirect exposure to the Greek sovereign (bonds, loans and guarantees) was 23.3 bn at end-2010. — [NBG announcement, ATHEX, 2011](https://www.athexgroup.gr/documents/10180/43440/57_17431_2011_English_2.pdf/ab4f0af4-1952-44a4-9e98-85a3646ec17b)
- After the July 2011 PSI agreement, banks impaired GGBs by 16% on average, against a 21% benchmark. Accounting losses were about 5 bn, about 2 pp off CT1. The system's aggregated CT1 ratio fell to 8.6%. — [IMF CR 11/351](https://www.imf.org/external/pubs/ft/scr/2011/cr11351.pdf)
- IMF estimate: under a 50% haircut, the six largest banks would need about 17 bn to reach a 10% CT1 ratio. Total HFSF needs were preliminarily put at about 40 bn (30 bn of additional support). — [IMF CR 11/351, Box 3 and DSA annex](https://www.imf.org/external/pubs/ft/scr/2011/cr11351.pdf)
- An IMF report of 16 Mar 2012 said the PSI "will trigger impairments of about €22 billion". — cited in [ZTG 2013, fn. 34](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf). This disagrees with the BoG figure of 37.7 bn: the IMF number is net of prior impairments, has a narrower perimeter and was made at an earlier date.
- Financial soundness indicators. — [IMF CR 11/351, Table 3](https://www.imf.org/external/pubs/ft/scr/2011/cr11351.pdf)
  - Regulatory capital/RWA: 11.7% (2009), 12.2% (2010Q4), 10.6% (2011Q2).
  - Tier 1/RWA: 10.7% (2009), 10.9% (2010Q4), 9.6% (2011Q2).
  - NPL/gross loans: 7.7% (2009), 10.4% (2010Q4), 12.8% (2011Q2).
- Aggregated MFI balance sheet (excluding the BoG), in bn. — [IMF CR 11/351, Table 2](https://www.imf.org/external/pubs/ft/scr/2011/cr11351.pdf)
  - Total assets: 579.2 (Dec-09), 654.6 (Dec-10), 620.3 (Dec-11 proj.).
  - Securities: 116.0 (Dec-09), 101.2 (Dec-10).
  - Credit to government: 51.0 (Dec-09), 63.0 (Dec-10).
  - Capital and reserves: 41.6 (Dec-09), 46.7 (Dec-10).
  - Capital/assets: 7.2% (Dec-09), 7.1% (Dec-10), 7.8% (Dec-11 proj.).
  - Liabilities to the BoG (Eurosystem funding): 49.7 (Dec-09), 97.8 (Dec-10).
  - Column alignment note: the PDF text extraction dropped the Dec-09 header. Columns are assigned on the assumption that the 12 values run Dec-06 to Dec-11. Check against the PDF.
- [Bank of Greece "Report on the Recapitalisation and Restructuring of the Greek Banking Sector" (Dec 2012), Chart I.1, reproduced in Hardouvelis & Vayanos 2023, Table 1](https://www.lse.ac.uk/Hellenic-Observatory/Assets/Documents/Publications/GreeSE-Papers/GreeSE-No180.pdf). Values in bn. PSI loss = losses on GGBs and loans to the state.

| Bank | CT1 Dec-2011 | PSI loss | PSI provisions (06/2011) | Credit-loss projection | Loan-loss reserves | Capital generation 2012-14 | Target CT1 Dec-2014 | Capital need |
|---|---|---|---|---|---|---|---|---|
| NBG | 7.29 | -11.74 | 1.65 | -8.37 | 5.39 | 4.68 | 8.66 | 9.76 |
| Eurobank | 3.52 | -5.78 | 0.83 | -8.23 | 3.51 | 2.90 | 2.60 | 5.84 |
| Alpha | 4.53 | -4.79 | 0.67 | -8.49 | 3.12 | 2.43 | 2.03 | 4.57 |
| Piraeus | 2.62 | -5.91 | 1.01 | -6.28 | 2.57 | 1.08 | 2.41 | 7.34 |
| ATEbank | 0.38 | -4.33 | 0.84 | -3.38 | 2.34 | 0.47 | 1.23 | 4.92 |
| Postbank | 0.56 | -3.44 | 0.57 | -1.48 | 1.28 | -0.32 | 0.90 | 3.74 |
| Total (14 banks) | 22.12 | -37.73 | 5.86 | -46.83 | 24.73 | 11.38 | 20.06 | 40.54 |

  Net-of-provision sovereign losses were 31.87 bn, which "wiped out completely the capital of the banks". Net credit losses were 22.10 bn. The four core banks' needs sum to 27.51 bn.
- The BoG confirmed that the 50 bn earmarked under the second programme was adequate. The report was completed in March 2012 and published in late December 2012. — [Eurobank Greece Macro Focus, 25 Jan 2013](https://www.eurobank.gr/Uploads/Reports/GREECE%20MACRFOCUS%20%20January%2025%202013.pdf)
- Bridge recapitalisation, 28 May 2012: the HFSF deposited EFSF floating-rate notes counted as equity. The split was NBG 6.9, Piraeus 5.0, Eurobank 4.2 and Alpha 1.9 bn, 18 bn in total. — [thecorner.eu](https://thecorner.eu/?p=12463); [Vietnamnet/Xinhua](https://vietnamnet.vn/en/greek-leading-banks-recapitalized-with-18-bln-euros-E22690.html). These are secondary sources, and the search snippet said "May 28" without the year. The 2012 dating comes from the subsequent EFSF disbursement and is my inference.
- Final recapitalisation, completed July 2013. — [Hardouvelis & Vayanos 2023, Sec. 4.1](https://www.lse.ac.uk/Hellenic-Observatory/Assets/Documents/Publications/GreeSE-Papers/GreeSE-No180.pdf)
  - Public funds to the four core banks: 25.5 bn (Alpha 4.0, Eurobank 5.8, NBG 8.7, Piraeus 7.0).
  - Private capital: 3.1 bn (Alpha 0.6, NBG 1.1, Piraeus 1.4). Eurobank failed to raise the required 10% privately and became fully state-controlled.
  - Resolutions used another 13.4 bn, so 38.9 bn of the 50 bn EFSF loan was used in total.
  - Contradicted by an academic summary surfaced in search, which says 48.2 bn was injected by June 2013. That figure probably includes the 2012 bridge notes later converted, plus resolution funding; treat it as a different definition. — [Emerald JFRC paper (search snippet)](https://www.emerald.com/insight/content/doi/10.1108/JFRC-12-2014-0050/full/html)
- Funding stress (from [Hardouvelis & Vayanos 2023](https://www.lse.ac.uk/Hellenic-Observatory/Assets/Documents/Publications/GreeSE-Papers/GreeSE-No180.pdf)):
  - Between the pre-crisis peak and April 2012, interbank borrowing fell from 78.6 to 40.7 bn and domestic deposits from 237.5 to 166 bn (-30.1%).
  - ECB direct loans rose from 6.3 bn (Mar 2008) to 48.1 bn (Mar 2010).
  - Direct loans plus ELA peaked at 159.2 bn.
- Market concentration: in 2013 the five largest banks held 94% of banking-system assets. — [Hardouvelis & Vayanos 2023](https://www.lse.ac.uk/Hellenic-Observatory/Assets/Documents/Publications/GreeSE-Papers/GreeSE-No180.pdf)
- Greek banks were excluded from the EBA December 2011 capital exercise, which started with 71 banks and ended with 61. No new benchmarks were set because the programme backstop exceeded the EBA exercise, so the exercise has no Greek data. — [EBA, Dec 2011 recommendation and results](https://www.eba.europa.eu/the-eba-publishes-recommendation-and-final-results-of-bank-recapitalisation-plan-as-part-of-co-ordinated-measures-to-restore-confidence-in-the-banking)

### Inferences
- Where the 45 bn (IMF) and 54.4 bn (EBA) exposure figures disagree, the cause is definitional. EBA item 34010 is gross direct long exposure, which includes loans to the state, booked gross of provisions, at Dec-2010. The IMF figure is GGBs only, at late 2011. Either gives an exposure/capital ratio of 2.0-2.4.
- The BoG PSI loss of 37.73 bn also includes "other loans to the Greek state" and subsidiaries' holdings, so it exceeds 53.5% x 45 bn (about 24 bn). For a mark-to-market moment, the consistent pair is PSI loss/CT1 = 1.71 gross or 1.44 net.
- The EBA file's MTM loss on the own-sovereign ladder over 2011 (9.05 bn, 39.7% of CT1) is a pure price-revaluation number, with no default. It sits well below the realised PSI loss, which is the expected ordering for a priced-not-realised calibration.

### Gaps
- The Bank of Greece Financial Stability Review figures (system total assets, CT1 and GGB holdings split banking/trading book, 2010-2011) were not fetched directly. The banking-vs-trading split is available only via the EBA file's 11.1% fair-valued share.
- The BoG December 2012 report itself was not opened; the numbers above come via Hardouvelis & Vayanos and the EC decision snippet. The EC snippet's NBG 11.7, Piraeus 5.9, Eurobank 5.8 and Alpha 4.8 bn agree with it.

## Q5. German banks (F bank): GGB and Bund holdings, assets, capital

### Takeaway
At Dec-2010 the 12 German EBA banks had:
- CT1 114.3 bn, RWA 1,222 bn and total assets 4,872 bn.
- Bund exposure of 315.3 bn (2.76x CT1, 6.5% of assets).
- Greek sovereign exposure of 7.9 bn (6.9% of CT1).

German Maastricht debt was 2,119 bn (81.0% of GDP) in 2010.

### Cited Findings
- EBA file, DE017-DE029 (12 banks), Dec-2010:
  - CT1 114,317 m; RWA 1,222,402 m; TA 4,872,189 m.
  - Bund exposure (`b_F_F`) 315,313 m; GGB exposure (`b_D_F`) 7,934 m.
  - Corporate EAD 960,328 m; CRE EAD 237,207 m.
  - Fair-valued share of the own-sovereign book: 45.4%.
  - Source: `data/eba_moments.json`.
- Computed: CT1/RWA = 9.35%; TA/CT1 = 42.6; Bund/TA = 6.47%; Bund/CT1 = 2.76; GGB/CT1 = 0.069.
- On the Bund ladder, the WAM is 4.86 y and the modified duration 4.22. The 2011 MTM was a gain of 14.1 bn as Bund yields fell from 2.91% to 1.93%. — `data/eba_moments.json` ladder.b_F_F
- German government debt and GDP: see the Q1 table. — [Eurostat](https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/gov_10dd_edpt1?geo=DE&sector=S13&na_item=GD&unit=MIO_EUR)

### Inferences
- The F bank's exposure to D-sovereign risk is small relative to capital (about 7% of CT1). Contagion to F through direct GGB losses is therefore weak in the data. A 64.6% PV haircut on 7.9 bn is about 5.1 bn, roughly 4.5% of German CT1.

### Gaps
- Bundesbank Monthly Report or Financial Stability Review series on German banks' GGB holdings by quarter, 2010-2012: not retrieved directly. Only OSW's secondary report was found.
- EBA 2011 published German bank-level Greek exposures as of Dec-2010 only. The Dec-2011 capital-exercise sovereign disclosures (Sep-2011 data) for German banks were not pulled.

## Q6. Maturity and duration of Greek debt and of bank-held GGBs

### Takeaway
Greek banks' own-sovereign book had a weighted-average residual maturity of 5.1 years. At the Dec-2010 yield of 12.0% its Macaulay duration was 3.5 years and its modified duration 3.1 years. Eurosystem SMP Greek holdings had 3.6 years remaining at end-2012. Post-PSI, the new bonds amortise over 2023-2042. No official pre-PSI average maturity for total Greek debt was found.

### Cited Findings
- Greek banks' own-sovereign maturity ladder, EBA file, 7 buckets, Dec-2010:
  - Face 54.4 bn; WAM 5.13 y.
  - Yield 12.01% (Dec-2010), 21.14% (Dec-2011).
  - Market value 41.9 bn.
  - Macaulay duration 3.49 y; modified duration 3.12 y.
  - 2011 MTM loss: 21.6% of MV (9.05 bn).
  - Implied HM perpetuity decay `delta_b_D` = 0.0777 per quarter, at an assumed 4.7% coupon.
  - Source: `data/eba_moments.json`.
- Duration sensitivity to the coupon assumption is small: modified duration is 3.10 at a 3% coupon and 3.13 at 6%. — same file
- German banks' GGB ladder: WAM 7.41 y and modified duration 4.05, longer than Greek banks' own book. — same file
- SMP Greek holdings had an average remaining maturity of 3.6 y at 31 Dec 2012, against 4.3 y for the whole SMP. — [ECB, 21 Feb 2013](https://www.ecb.europa.eu/press/pr/date/2013/html/pr130221_1.en.html)
- Bond amortisations due mid-2011 to 2020, excluding ECB holdings: 203 bn, of which 89 bn fell due mid-2011 to mid-2014. The ECB and other central banks held about 53 bn maturing 2012-2020 (26 bn in 2012-2014). — [ZTG 2013, fn. on the IIF offer](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf)
- New PSI bonds mature 2023-2042, 20 series amortising. Second-programme EFSF loans have a weighted-average maturity of 42.5 years, and ESM loans 32.5 years. — [ZTG 2013](https://www.ifo.de/DocDL/cesifo1_wp4333.pdf); EFSF/ESM maturity from a [secondary FES study](https://collections.fes.de/publikationen/download/ftpack/plain/450116)
- The initial (July 2011) PSI covered paper maturing up to 2020. PSI-plus covered about 206 bn with maturities up to 2035. — [Eurobank, Greece Debt Sustainability Outlook, 31 Oct 2011](https://www.eurobank.gr/Uploads/Reports/31%2010%2011%20Greece%20Debt%20Sustainability%20Outlook.pdf)

### Inferences
- If 203 bn of the roughly 260 bn bonded stock fell due by 2020, the bonded debt was short. A WAM of about 5-7 years is plausible and consistent with the EBA bank ladder (5.1 y), but this is an inference and not an official figure.
- `delta_b_D` = 0.0777 per quarter implies an HM perpetuity with an expected life of about 1/0.0777 = 12.9 quarters, or 3.2 years. That matches the measured duration.

### Gaps
- Official average residual maturity of Greek central government debt for 2009-2012, before and after the PSI, from the PDMA Public Debt Bulletin: not found in search. The PDMA bulletins (pdma.gr) need a direct fetch.

## Q7. Bank leverage measures and literature targets

### Takeaway
Leverage depends heavily on the definition:
- Greek banks (Dec-2010): about 14x on accounting assets/equity (IMF MFI data), 16.6x on TA/CT1 (EBA) and 5.5x on credit-EAD/CT1 (the EBA file's `theta_D`).
- German banks: 42.6x on TA/CT1 (inflated by derivatives and interbank assets) and 6.9x on EAD/CT1 (`theta_F`).

The model's choice of leverage 8 (per the calibration-branch memory) sits between these definitions.

### Cited Findings
- Greek MFIs (ex-BoG): capital and reserves / total assets was 7.2% at Dec-09 and 7.1% at Dec-10, which is assets/equity of about 14. — [IMF CR 11/351, Table 2](https://www.imf.org/external/pubs/ft/scr/2011/cr11351.pdf)
- From the EBA file, Dec-2010, `leverage_alternatives` and `model_targets`:
  - Greek banks: TA/CT1 16.56; total-book/CT1 7.21; (corporate + CRE + sovereign) EAD/CT1 5.51 (`theta_D`).
  - German banks: TA/CT1 42.62; total-book/CT1 13.30; EAD/CT1 6.94 (`theta_F`).
  - Source: `data/eba_moments.json`.
- Risk-based ratios: Greek Tier 1/RWA was 10.9% at 2010Q4 and 9.6% at 2011Q2 (IMF FSI). The aggregated CT1 ratio was 8.6% after the July 2011 impairments. EBA CT1/RWA was 10.2% for Greece and 9.35% for Germany. — [IMF CR 11/351](https://www.imf.org/external/pubs/ft/scr/2011/cr11351.pdf); `data/eba_moments.json`
- The EBA capital exercise target was a 9% CT1 ratio after a sovereign buffer. The 2011 stress-test benchmark was 5% CT1. — [EBA Dec 2011](https://www.eba.europa.eu/the-eba-publishes-recommendation-and-final-results-of-bank-recapitalisation-plan-as-part-of-co-ordinated-measures-to-restore-confidence-in-the-banking)

### Inferences
- German TA/CT1 (42.6) is not comparable with Greek TA/CT1 (16.6). German total assets carry large derivative and interbank books that GK-style lending does not capture. The EAD-based `theta` (5.5 vs 6.9) or the total-book measure (7.2 vs 13.3) are the more model-consistent choices.

### Gaps
- Literature leverage targets (Gertler-Karadi 2011 JME, Gertler-Kiyotaki, Bocola 2016 JPE's Italian bank leverage, and euro-area GK applications) were not verified in this session. I recall a GK (2011) steady-state leverage of 4, but that is unverified and should be checked against the papers before it is cited.
- Basel III-style leverage ratios (Tier 1/exposure) for Greek and German banks in 2010-2012 were not found.
