# Greek bank funding and the Eurosystem backstop, 2009–2013 (plus PSPP/PEPP/TPI rules)

Notes for calibrating the funding side and the central-bank backstop of the two-country model.
Numbers marked "ECB Data Portal" were downloaded on 2026-10-23 straight from the ECB SDMX API
(`https://data-api.ecb.europa.eu/service/data/<dataset>/<key>`). The annual means and the picked
months were computed from those monthly series in this session. Units are stated on each line.

## 1. Greek bank funding: deposits, Eurosystem refinancing, ELA, TARGET2

### Takeaway
Private-sector deposits at Greek MFIs fell by about EUR 88bn (−36%), from a peak of EUR 244.2bn
(Dec 2009) to a trough of EUR 156.2bn (Jun 2012). Eurosystem lending replaced them. Direct
monetary-policy loans and ELA together peaked at about EUR 159bn (Feb 2012), and ELA alone
peaked at about EUR 124bn (May 2012). The Bank of Greece TARGET2 liability peaked at
EUR −109.3bn (Nov 2011) and stayed above EUR 100bn in absolute value through mid-2012.

### Cited findings
**Deposits (ECB BSI, stocks at end of month, EUR million, Greek MFIs excl. ESCB, deposits from euro-area non-MFIs excl. central government, all maturities and currencies)**
- Series `BSI.M.GR.N.A.L20.A.1.U2.2300.Z01.E`. Peak 244,237 in Dec 2009 (244,175 in Sep 2009). Then 228,350 (Apr 2010), 215,639 (Dec 2010), 195,069 (Jun 2011), 180,085 (Dec 2011), 170,176 (Feb 2012), 163,090 (May 2012), trough 156,202 (Jun 2012), 167,930 (Dec 2012), 170,874 (Dec 2013). In 2015 it fell further, to 125,904 (Jul 2015). Annual means: 2009 240.2bn, 2010 224.5bn, 2011 195.3bn, 2012 164.6bn, 2013 170.0bn. Dec 2009 to Jun 2012 is −88.0bn (−36.0%) — [ECB Data Portal, BSI](https://data-api.ecb.europa.eu/service/data/BSI/M.GR.N.A.L20.A.1.U2.2300.Z01.E?format=csvdata)
- Households only (`...U2.2250...`): 197,419 (Dec 2009) to 127,819 (Jun 2012), −69.6bn (−35.3%). Then 135,443 (Dec 2012) and 135,401 (Dec 2013) — [ECB Data Portal, BSI households](https://data-api.ecb.europa.eu/service/data/BSI/M.GR.N.A.L20.A.1.U2.2250.Z01.E?format=csvdata)
- Non-financial corporations only (`...U2.2240...`): 37,129 (Dec 2009) to 18,181 (Jun 2012), −18.9bn (−51%). Then 20,754 (Dec 2012) — [ECB Data Portal, BSI NFC](https://data-api.ecb.europa.eu/service/data/BSI/M.GR.N.A.L20.A.1.U2.2240.Z01.E?format=csvdata)
- Vives/Haliassos-style case study on Bank of Greece data: "between the end of 2009 and April 2012 … interbank borrowing dropped from €78.6bn to €40.7bn (−48.2%) and domestic deposits dropped from €237.5bn to €166bn (−30.1%)". Deposits fell from 58.1% of liabilities to 49.8% by end-2012 — [DNB conference paper "The Greek Economic Crisis and the Banks"](https://www.dnb.nl/media/ec2dqfwd/presentation-session-iii.pdf) (same paper: [LSE GreeSE No. 180](https://www.lse.ac.uk/Hellenic-Observatory/Assets/Documents/Publications/GreeSE-Papers/GreeSE-No180.pdf))
- A Bank of Greece working paper puts the deposit loss at roughly 40%, about EUR 90bn — [BoG WP 262, Louri & Migiakis](https://www.bankofgreece.gr/Publications/Paper2019262.pdf). The definition differs: the ~30% figure counts domestic private deposits, the 36–40% figures count a wider aggregate or a later trough.

**Eurosystem refinancing and ELA**
- "Direct loans [ECB monetary-policy operations] … averaged €44.2bn in 2009, rose to an average of €95.4bn in the second half of 2010 and €93.4bn in the first half of 2011. During that period, ELA was essentially zero … The sum of direct loans and ELA reached a peak of €159.2bn in February 2012. ELA was wound down during 2013." Direct loans had already risen from EUR 6.3bn (Mar 2008) to EUR 48.1bn (Mar 2009) — [DNB / GreeSE paper](https://www.dnb.nl/media/ec2dqfwd/presentation-session-iii.pdf)
- ELA outstanding peaked at EUR 124bn in May 2012 and again at EUR 90bn in May 2015. ELA borrowing peaked at about 65% of Greek nominal GDP in 2012 — [Yale Journal of Financial Crises 4(2), 2022, Greek ELA case study](https://elischolar.library.yale.edu/cgi/viewcontent.cgi?article=1328&context=journal-of-financial-crises) (search snippet; the PDF returned 403)
- A March 2012 IMF report says Eurosystem exposure peaked at EUR 103bn in mid-2011, before ELA took over from the ECB window — [IMF Country Report 12/57](https://www.imf.org/external/pubs/ft/scr/2012/cr1257.pdf) (search snippet)
- In July 2015 the ELA ceiling was held at EUR 89bn — [Business Standard / PTI, 8 Jul 2015](https://www.business-standard.com/amp/article/pti-stories/ecb-leaves-emergency-funding-to-greece-unchanged-at-89bn-euros-115070801425_1.html)
- The sources disagree on the peak of total Eurosystem support. The DNB/GreeSE paper gives EUR 159.2bn (Feb 2012, direct loans plus ELA). A GreeSE snippet also cites EUR 126.6bn for a combined figure without a date, and the IMF gives EUR 103bn (mid-2011).

**TARGET2 (ECB TGB, end-of-month position, EUR million; negative = liability to the Eurosystem)**
- Bank of Greece: −49,036 (Dec 2009), −82,604 (Apr 2010), −87,088 (Dec 2010), −96,802 (Jun 2011), **minimum −109,315 (Nov 2011)**, −104,750 (Dec 2011), −107,267 (Feb 2012), −101,555 (May 2012), −105,987 (Jun 2012), −98,355 (Dec 2012), −51,116 (Dec 2013), −107,702 (Jun 2015). Annual means: 2010 −79.5bn, 2011 −94.2bn, 2012 −105.0bn, 2013 −62.0bn — [ECB Data Portal, TGB GR](https://data-api.ecb.europa.eu/service/data/TGB/M.GR.N.A094T.U2.EUR.E?format=csvdata)
- Bundesbank, for comparison: +177,723 (Dec 2009), +463,134 (Dec 2011), **peak +751,449 (Aug 2012)**, +655,670 (Dec 2012) — [ECB Data Portal, TGB DE](https://data-api.ecb.europa.eu/service/data/TGB/M.DE.N.A094T.U2.EUR.E?format=csvdata)
- The Nov 2011 peak of EUR 109.3bn is confirmed independently — [ResearchGate, "TARGET2 Imbalances and the ECB as Lender of Last Resort"](https://www.researchgate.net/publication/283664597_TARGET2_Imbalances_and_the_ECB_as_Lender_of_Last_Resort). Whelan gives EUR 87bn for Jan 2013 — [Whelan, "TARGET2 and Central Bank Balance Sheets", 2013](https://www.karlwhelan.com/Papers/T2Paper-March2013.pdf)

**ELA pricing and collateral**
- The Bank of Greece provided ELA at a 100–150 bp premium over the ECB refinancing (MRO) rate — [Yale JFC case study](https://elischolar.library.yale.edu/cgi/viewcontent.cgi?article=1328&context=journal-of-financial-crises) (snippet). A Fed note lists euro-area ELA generically at "+25 bp over MRO", with terms "determined by the NCB" — [Federal Reserve FEDS Note, 26 Feb 2025](https://www.federalreserve.gov/econres/notes/feds-notes/central-bank-liquidity-facilities-around-the-world-20250226.html). The two figures conflict, and the 100–150 bp figure is the Greece-specific one. The ECB's 2024 ELA Agreement requires only "a penalty interest rate" and sets no number — [ECB ELA Agreement, 27 Sep 2024](https://www.ecb.europa.eu/pub/pdf/other/ecb.agreementemergencyliquidityassistance202410~b8b78cd4f5.en.pdf)
- ELA is granted by the NCB at its own risk and against less stringent collateral, at "a significantly higher interest rate" than ECB funding — [DNB / GreeSE paper](https://www.dnb.nl/media/ec2dqfwd/presentation-session-iii.pdf)
- Waiver chronology:
  - The ECB revoked the first waiver on 27 Feb 2012, after the PSI exchange — [Yale JFC](https://elischolar.library.yale.edu/cgi/viewcontent.cgi?article=1328&context=journal-of-financial-crises)
  - On 8 Mar 2012 the ECB again accepted Greek-issued or Greek-guaranteed instruments "without the minimum credit rating threshold … until further notice", backed by an EFSF buy-back scheme of EUR 35bn nominal — [ECB press release 8 Mar 2012](https://www.ecb.europa.eu/press/pr/date/2012/html/pr120308_1.en.html)
  - Yale JFC lists a second waiver revoked on 18 Jul 2012 — [Yale JFC](https://elischolar.library.yale.edu/cgi/viewcontent.cgi?article=1328&context=journal-of-financial-crises)
  - On 4 Feb 2015 the ECB lifted the waiver, effective at the maturity of the then-current MRO (11 Feb 2015) — [ECB press release 4 Feb 2015](https://www.ecb.europa.eu/press/pr/date/2015/html/pr150204.en.html); see also [Bruegel, "ECB collateral damages on Greece"](https://www.bruegel.org/blog-post/ecb-collateral-damages-greece)
- PSI context: GGBs with face value EUR 199.2bn were exchanged for new GGBs of EUR 62.4bn plus EFSF notes of EUR 29.7bn. Greek banks lost EUR 37.73bn on GGBs (EUR 31.87bn net of EUR 5.86bn provisions), against aggregate CT1 capital of EUR 22.12bn (Dec 2011). The Dec 2012 buyback retired EUR 31.9bn face value for EUR 11.3bn — [DNB / GreeSE paper](https://www.dnb.nl/media/ec2dqfwd/presentation-session-iii.pdf)

### Inferences
- For the model, take the deposit-flight shock as about −36% of private deposits in 30 months (Dec 2009 to Jun 2012, ECB BSI). Use −30% (Bank of Greece "domestic deposits" to Apr 2012) as the lower bound.
- The TARGET2 liability (about EUR 100–109bn) is roughly the same size as the deposit outflow (about EUR 88bn) plus the lost interbank funding (about EUR 38bn) less Greek repayments. This fits the view that deposit flight was financed by Eurosystem credit and settled through TARGET2.
- The ELA share of Eurosystem funding rose sharply after GGB collateral became ineligible (Feb 2012, Jul 2012). The cost to banks of this backstop is therefore MRO + 100–150 bp, not the MRO rate.

### Gaps
- I found no primary Bank of Greece monthly series splitting Eurosystem lending into MRO/LTRO and ELA. ECB ILM keys for GR returned no data, and the BoG Annual Report 2012 PDF returned 403. The EUR 124bn ELA peak (May 2012) comes from a secondary source.
- Official ELA rates have never been published. The 100–150 bp figure comes from a secondary source (press reports and the Yale case study).
- I could not confirm the dates on which the waiver was reinstated after Jul 2012 (reportedly Dec 2012). Haircut schedules for GGB collateral were not retrieved.

## 2. Interest rates: deposit and lending rates in Greece vs Germany, and the MRO path

### Takeaway
In 2012 Greek banks paid households about 4.8% on new term deposits, against 1.2% in Germany (a spread of about 3.6 pp). They charged non-financial corporations about 6.6% (composite cost of borrowing), against 2.9% in Germany. The policy rate (MRO) was 0.75–1.00%. The Greek loan-to-policy spread is therefore about 5.7 pp and the loan-to-deposit spread about 1.8 pp. In Germany the loan-to-policy spread is about 2.0 pp and the loan-to-deposit spread about 1.8 pp.

### Cited findings
**MFI interest rates, new business, % p.a. Each line gives monthly values (Sep 2009, Dec 2009, Dec 2010, Dec 2011, Jun 2012, Dec 2012, Dec 2013) and then annual means for 2010 / 2011 / 2012 / 2013.**
- Household deposits with agreed maturity up to 1 year:
  - Greece: 2.08, 2.10, 3.66, 4.86, **5.01**, 4.68, 2.80. Means 3.25 / 4.17 / **4.77** / 3.66 — [ECB MIR GR](https://data-api.ecb.europa.eu/service/data/MIR/M.GR.B.L22.F.R.A.2250.EUR.N?format=csvdata)
  - Germany: 1.01, 0.87, 1.08, 1.38, 1.15, 0.75, 0.61. Means 1.10 / 1.40 / **1.16** / 0.72 — [ECB MIR DE](https://data-api.ecb.europa.eu/service/data/MIR/M.DE.B.L22.F.R.A.2250.EUR.N?format=csvdata)
- NFC deposits with agreed maturity up to 1 year:
  - Greece: means 2.53 / 3.55 / **4.19** / 3.35 (Jun 2012: 4.17) — [ECB MIR GR NFC](https://data-api.ecb.europa.eu/service/data/MIR/M.GR.B.L22.F.R.A.2240.EUR.N?format=csvdata)
  - Germany: means 0.57 / 1.11 / **0.44** / 0.21 — [ECB MIR DE NFC](https://data-api.ecb.europa.eu/service/data/MIR/M.DE.B.L22.F.R.A.2240.EUR.N?format=csvdata)
- Composite cost-of-borrowing indicator for NFC loans (`A2I.AM`):
  - Greece: 4.60, 4.61, 5.83, 7.10, **6.94**, 6.53, 5.41. Means 5.21 / 6.46 / **6.57** / 6.25 — [ECB MIR GR CCBI](https://data-api.ecb.europa.eu/service/data/MIR/M.GR.B.A2I.AM.R.A.2240.EUR.N?format=csvdata)
  - Germany: 3.13, 3.20, 3.35, 3.50, 3.06, 2.61, 2.60. Means 3.16 / 3.53 / **2.94** / 2.57 — [ECB MIR DE CCBI](https://data-api.ecb.europa.eu/service/data/MIR/M.DE.B.A2I.AM.R.A.2240.EUR.N?format=csvdata)
- The Greek series for NFC loans with floating rate or fixation up to 1 year, all amounts (`A2A.A`), has gaps for 2011–2012 in the API. German values: means 2.76 / 3.15 / 2.51 / 2.11 — [ECB MIR DE A2A](https://data-api.ecb.europa.eu/service/data/MIR/M.DE.B.A2A.A.R.A.2240.EUR.N?format=csvdata)
- The BoG Annual Report 2012 shows charts of Greek and euro-area new household deposit rates and the Greek loan–deposit spread. It also reports that the Greek–German 10y spread fell from about 3,000 bp (May 2012) to under 1,000 bp — [Bank of Greece Annual Report 2012](https://www.bankofgreece.gr/Publications/Annrep2012.pdf) (search snippet; the PDF returned 403)

**ECB key rates. Each line gives the effective date, then MRO / deposit facility / marginal lending facility in %:**
- 13 May 2009: 1.00 / 0.25 / 1.75
- 13 Apr 2011: 1.25 / 0.50 / 2.00
- 13 Jul 2011: 1.50 / 0.75 / 2.25
- 9 Nov 2011: 1.25 / 0.50 / 2.00
- 14 Dec 2011: 1.00 / 0.25 / 1.75
- 11 Jul 2012: 0.75 / 0.00 / 1.50
- 8 May 2013: 0.50 / 0.00 / 1.00
- 13 Nov 2013: 0.25 / 0.00 / 0.75
- 11 Jun 2014: 0.15 / −0.10 / 0.40
- 10 Sep 2014: 0.05 / −0.20 / 0.30
- Before May 2009: 21 Jan 2009 2.00, 11 Mar 2009 1.50, 8 Apr 2009 1.25 (MRO)
- Source: [ECB key interest rates](https://www.ecb.europa.eu/stats/policy_and_exchange_rates/key_ecb_interest_rates/html/index.en.html)

### Inferences
- The 2012 average MRO was about 0.88% (1.00% to 11 Jul, 0.75% after). Computed against it:
  - Greek household term-deposit rate: about +3.9 pp over MRO (Germany about +0.3 pp).
  - Greek NFC CCBI: about +5.7 pp over MRO (Germany about +2.1 pp).
  - Greek lending spread over the deposit rate: 6.57 − 4.77 ≈ 1.8 pp, almost the same as Germany's 2.94 − 1.16 ≈ 1.8 pp.
- So in the data, almost all of the Greek–German gap in lending rates is in the funding cost (the deposit rate), not in the loan–deposit margin. If the model's lending spread is measured over the deposit rate, it should be roughly equal across countries. If it is measured over the policy rate, the Greek spread should be about 3.6 pp higher.
- Greek deposit rates rose by about 2.9 pp from Sep 2009 to Jun 2012 while the MRO was flat or falling. This is the clearest measure of the rise in Greek bank funding costs.

### Gaps
- Monthly Greek NFC loan rates by fixation period for 2011–2012 were missing from the API pull. Use the CCBI instead.
- I did not retrieve rates on outstanding amounts (stocks), which may matter for a predetermined deposit rate.

## 3. SMP: start, size, Greek holdings, price impact, PSI exemption

### Takeaway
The SMP began on 10 May 2010 and ended on 6 Sep 2012 when OMT was announced. At end-2012 it held EUR 218.0bn nominal (EUR 208.7bn book). Greek bonds were EUR 33.9bn nominal and EUR 30.8bn book, with 3.6 years average remaining maturity. Eser & Schwaab estimate that each EUR 1bn of purchases lowered Greek 5y yields by about 16–21 bp. Eurosystem holdings were swapped for new bonds with the same terms in Feb 2012 and so were exempt from the PSI haircut.

### Cited findings
- Table of SMP holdings at 31 Dec 2012 (nominal EUR bn / book EUR bn / average remaining maturity in years), held to maturity at amortised cost — [ECB press release 21 Feb 2013](https://www.ecb.europa.eu/press/pr/date/2013/html/pr130221_1.en.html):
  - Ireland: 14.2 / 13.6 / 4.6
  - **Greece: 33.9 / 30.8 / 3.6**
  - Spain: 44.3 / 43.7 / 4.1
  - Italy: 102.8 / 99.0 / 4.5
  - Portugal: 22.8 / 21.6 / 3.9
  - **Total: 218.0 / 208.7 / 4.3**
- OMT press release: "the Securities Markets Programme (SMP) is herewith terminated. The liquidity injected through the SMP will continue to be absorbed as in the past" — [ECB press release 6 Sep 2012](https://www.ecb.europa.eu/press/pr/date/2012/html/pr120906_1.en.html)
- Greek SMP purchases in May–June 2010 were about EUR 40bn, reportedly bought at a discount — [Positive Money Europe](https://positivemoney.org/eu/archive/ecb-smp-profits-billions-lost-for-greece/) (advocacy source, secondary)
- Eser & Schwaab (2016, JFE 119(1), 147–167): purchases of 1/1000 of outstanding debt lowered 5y yields by about 3 bp (pooled). Per EUR 1bn, the effect ranged from −1 to −2 bp (Italy) to −16 to −21 bp (Greece) at the 5y maturity. The working-paper version reports up to −17 to −21 bp for Greece. Sample: Greece, Ireland, Italy, Portugal and Spain, 2010–11 — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0304405X15000963); [RePEc](https://ideas.repec.org/a/eee/jfinec/v119y2016i1p147-167.html); [ECB WP 1587](https://www.ecb.europa.eu/pub/pdf/scpwps/ecbwp1587.pdf)
- Ghysels, Idier, Manganelli & Vergote ("A high frequency assessment of the ECB SMP", ECB WP 1642, 2014; published 2017): daily regressions show no significant yield effect of SMP interventions, because purchases are triggered by price drops. At 15-minute frequency the expected negative sign appears — [ECB WP 1642](https://www.ecb.europa.eu/pub/pdf/scpwps/ecbwp1642.pdf); [CEPR DP 9778](https://ideas.repec.org/p/cpr/ceprdp/9778.html)
- Krishnamurthy, Nagel & Vissing-Jorgensen (2018, Review of Finance 22(1)): event studies of SMP and OMT find large yield declines, concentrated at 2y maturities. They work mainly through lower default-risk premia and market segmentation, with a smaller role for redenomination risk — [NBER WP 23985](https://www.nber.org/system/files/working_papers/w23985/w23985.pdf); [RoF](https://academic.oup.com/rof/article-abstract/22/1/1/4665133). The NBER Digest summary gives combined policy-event declines in 2y yields of about 400 bp (Italy, Spain), about 500 bp (Ireland, Portugal) and **about 1,000 bp (Greece)** — [NBER Digest Feb 2018](https://www.nber.org/digest/feb18/quantifying-impact-ecb-policies-during-debt-crisis). This is a secondary summary that sums several announcements, so check the paper's tables before using a single-event number.
- PSI exemption: "the Bank, like the other Eurosystem central banks, has exchanged all the Greek government bonds that it held in its various portfolios for new securities issued by the Hellenic Republic … same nominal value, coupon and redemption dates, but … not … eligible … under the PSI" — [National Bank of Belgium](https://nbb.be/en/node/2715). The swap was agreed in Feb 2012 between the ECB and the Hellenic Republic and covered SMP and ANFA holdings — [Wikipedia, Private sector involvement](https://en.wikipedia.org/wiki/Private_sector_involvement) (secondary). The Eurogroup of 21 Feb 2012 set a 53.5% nominal haircut for private holders and noted that Eurosystem holdings were held for public-policy purposes — [Eurogroup statement 21 Feb 2012](https://www.esm.europa.eu/sites/default/files/2012-02-21_eurogroup_statement_bailout_for_greece.pdf)

### Inferences
- For the model's purchase-elasticity of the D bond price, use about 16–21 bp per EUR 1bn of 5y-equivalent purchases for Greece. In debt-share terms, use about 3 bp per 0.1% of the outstanding stock (Eser & Schwaab). Both are flow (impact) effects. Eser & Schwaab also report a stock effect, which was not retrieved here.
- SMP holdings of GGBs were about 33.9/(GGB stock) ≈ 12–17% of the Greek bond stock. This is my approximation: the post-PSI marketable stock was much smaller than the EUR 205.6bn pre-PSI face value cited by the DNB paper.
- Because of the PSI exemption, the Eurosystem was de facto senior in 2012. OMT reversed this by accepting pari passu treatment. For the model's "pure-haircut default branch", this means that central-bank-held D bonds took no haircut in the Greek case.

### Gaps
- I did not retrieve the SMP's peak holdings (often quoted as about EUR 219–220bn in early 2012) or the weekly purchase series.
- Ghysels et al.'s Greece-specific magnitudes were not retrieved. Only the method and sign were.

## 4. Profit sharing, capital keys, PSPP/PEPP risk-sharing rules, TPI

### Takeaway
On 27 Nov 2012 the Eurogroup committed member states to pass to Greece an amount equal to the SMP income their NCBs earned on GGBs from budget year 2013 onward, with ANFA income added. EUR 1.5bn (Jul 2013) and EUR 0.5bn (Dec 2013) were paid. Transfers were suspended in 2015 and resumed from 2017 in tranches of about EUR 600m. From 1 Jan 2009 the ECB capital key was Greece 1.9649% and Germany 18.9373%, which is about 2.81% and 27.1% of euro-area paid-up capital. The PSPP and PEPP public-sector purchases are 80% non-risk-shared (each NCB buys its own sovereign at its own risk) and 20% risk-shared. For the TPI, the ECB has published no loss-sharing rule.

### Cited findings
**Eurogroup Nov 2012 and transfers**
- The Eurogroup statement of 27 Nov 2012 contains "a commitment by Member States to pass on to Greece's segregated account, an amount equivalent to the income on the SMP portfolio accruing to their national central bank as from budget year 2013" — [Eurogroup statement 27 Nov 2012](https://www.consilium.europa.eu/uedocs/cms_Data/docs/pressdata/en/ecofin/133857.pdf)
- A release of EUR 2bn was mandated. EUR 1.5bn of SMP income from euro-area NCBs ("not counted in the financial envelope") went to Greece in Jul 2013 and EUR 0.5bn in Dec 2013 — [European Commission, Financial assistance to Greece](https://economy-finance.ec.europa.eu/eu-financial-assistance/euro-area-countries/financial-assistance-greece_en)
- Transfers were suspended during the 2015 Eurogroup–Varoufakis standoff, and funds already in the ESM account were blocked. They were reactivated in 2017, excluding the 2015–16 profits. Under the 22 Jun 2018 Eurogroup, future refunds were planned as EUR 600m semi-annual tranches for 2017–2022. Positive Money estimates Eurosystem SMP profits on GGBs at over EUR 13bn and ANFA profits at about EUR 4bn, and quotes an FT figure of EUR 7.8bn — [Positive Money Europe](https://positivemoney.org/eu/archive/ecb-smp-profits-billions-lost-for-greece/) (advocacy source, secondary; treat the profit totals as estimates)
- Later payments: in Dec 2019 the Eurogroup approved EUR 767m (SMP/ANFA) — [Keep Talking Greece, 4 Dec 2019](https://www.keeptalkinggreece.com/2019/12/04/eurogroup-greece-smp-anfa-768millioneuros/). The ESM transferred EUR 644m from SMP/ANFA income — [Hellenic News](https://hellenicnews.com/esm-transfers-to-greece-644-million-euros-from-income-earned-on-smp-anfa-holdings/). For the policy framework see the [ESM explainer on the 2018 final disbursement and debt relief](https://www.esm.europa.eu/content/explainer-final-disbursement-esm-financial-assistance-greece-and-medium-term-debt-relief)

**Capital key**
- Capital key weights from 1 Jan 2009, as % of the subscribed capital of all EU NCBs — [ECB press release 1 Jan 2009](https://www.ecb.europa.eu/press/pr/date/2009/html/pr090101_1.en.html):
  - **Bank of Greece 1.9649%** (EUR 113,191,059.06 paid up)
  - **Deutsche Bundesbank 18.9373%** (EUR 1,090,912,027.43)
  - Euro-area NCB subtotal 69.7915% (EUR 4,020,445,721.56)
  - Banque de France 14.2212%, Banca d'Italia 12.4966%, Banco de España 8.3040%
- Computed euro-area shares (key divided by the euro-area subtotal):
  - 2009–2010 (16 members): **Greece 2.8154%, Germany 27.1341%**
  - From 1 Jan 2011, when Estonia joined at 0.1790% and the subtotal became 69.9705%: **Greece 2.8082%, Germany 27.0647%**
  - The ECB capital increase of 29 Dec 2010 did not change the key (Bank of Greece subscribed EUR 211,436,059.06 of EUR 10,760,652,402.58 = 1.9649%) — [ECB Decision 2010/27, OJ L11](https://www.ecb.europa.eu/ecb/access_to_documents/document/legal_act/shared/data/ecb.dr.leg2011_L11_54.ga.pdf)
- After the 1 Jan 2014 recalculation the key is Greece 2.0332% and Germany 17.9973%. Their shares of euro-area paid-up capital are 2.9055% and 25.7184% (Bundesbank table, euro area including Lithuania from 2015) — [Bundesbank, Understanding the capital key](https://www.bundesbank.de/en/tasks/topics/understanding-the-capital-key-666388)

**PSPP and PEPP risk-sharing**
- PSPP (announced 22 Jan 2015): 20% of purchases are subject to loss sharing. This share was kept at 20% when the supranational share was cut and ECB purchases rose from 8% to 10%. NCBs buy only their own sovereign's bonds, and losses on those are not shared — [ECB PSPP FAQ](https://www.ecb.europa.eu/mopo/implement/app/html/ecb.faq_pspp.en.html); [ECB PSPP presentation, OMCG 17 Mar 2015](https://www.ecb.europa.eu/paym/groups/pdf/omg/2015/20150317/Presentation_PSPP_OMCG.pdf), which labels the 80% "non-loss-sharing"
- PEPP (Mar 2020): "the same risk-sharing principles apply as for the APP". Public-sector purchases follow the 20% shared / 80% NCB-own-risk split, and private-sector PEPP assets are fully risk-shared — [ECB PEPP FAQ](https://www.ecb.europa.eu/mopo/implement/pepp/html/ecb.faq_pepp.en.html)

**TPI**
- The TPI was approved on 21 Jul 2022. It requires a cumulative list of four criteria:
  - compliance with the EU fiscal framework
  - absence of severe macroeconomic imbalances
  - fiscal sustainability
  - sound and sustainable macroeconomic policies
- Purchases are secondary-market buys of public-sector securities with 1–10y remaining maturity, "not restricted ex ante" — [ECB press release 21 Jul 2022](https://www.ecb.europa.eu/press/pr/date/2022/html/ecb.pr220721~973e6e7273.en.html); for analysis see [European Parliament, "The ECB's TPI: a legal & economic analysis"](https://www.europarl.europa.eu/cmsdata/253891/QA-07-22-986-EN-N.pdf)
- Loss sharing is not specified: "it is not yet specified whether the risk (and returns) of TPI purchases will be borne individually by the NCBs … or will be shared according to the ECB capital key" — [OMFIF, Aug 2022](https://www.omfif.org/2022/08/ecb-transmission-protection-instrument-needs-a-user-manual/)
- OMT was intended to be fully risk-shared, and OMT "is incompatible with national risk-sharing" — [OMFIF](https://www.omfif.org/2022/08/ecb-transmission-protection-instrument-needs-a-user-manual/); [VoxEU, "ECB's bond purchase programmes and limits of national risk sharing"](https://cepr.org/voxeu/columns/ecbs-bond-purchase-programmes-and-limits-national-risk-sharing)

### Inferences
- In the model, the share of union-wide central-bank P&L that falls on country D should be about 2.8% under capital-key sharing, with Germany at about 27%. The model's D/F size weights may differ. In practice SMP income on GGBs went back to Greece only from 2013 and only in part, so Greece effectively bore the 2010–2012 losses and gains in proportion to its key.
- For the post-2015 regime, the cleanest parameterisation is: share of the D-bond purchase risk retained by the D national central bank = 0.8 (PSPP/PEPP). For SMP/OMT/TPI the retained share is 0 (fully mutualised by capital key), or unspecified for TPI.

### Gaps
- I found no official ECB figure for SMP income on GGBs by year. The total amount withheld in 2015 was not quantified by any primary source found.
- I did not confirm whether SMP income was shared by capital key. This is the standard Eurosystem monetary-income rule and is implied by the Eurogroup wording, but no SMP-specific primary text was retrieved.
- Whether TPI losses are shared remains unresolved in public sources (as of 2022 analysis). I found no later ECB clarification.

## 5. OMT: announcement effects and why Greece was not eligible

### Takeaway
OMT was announced in Draghi's 26 Jul 2012 speech and through the 2 Aug and 6 Sep 2012 decisions. It lowered Italian and Spanish 2y yields by about 2 pp and left German and French yields unchanged. Greece was ineligible because OMT requires strict conditionality under an EFSF/ESM programme. For countries already under an adjustment programme, OMT applies only "when they will be regaining bond market access", and Greece had no market access.

### Cited findings
- OMT requires "strict and effective conditionality" attached to an EFSF/ESM programme. For programme countries it would apply only "for Member States currently under a macroeconomic adjustment programme when they will be regaining bond market access". It targets maturities of 1–3 years, the Eurosystem accepts pari passu treatment, and no ex ante quantitative limits are set — [ECB press release 6 Sep 2012](https://www.ecb.europa.eu/press/pr/date/2012/html/pr120906_1.en.html); [Bundesbank OMT page](https://www.bundesbank.de/en/tasks/monetary-policy/outright-transactions/outright-monetary-transactions-omt--896164) ("to date, the Eurosystem has not conducted any OMTs")
- Altavilla, Giannone & Lenza: "OMT announcements decreased the Italian and Spanish 2-year government bond yields by about 2 percentage points, while leaving unchanged the bond yields of the same maturity in Germany and France". In a multi-country VAR, this is associated with significant increases in activity, credit and prices in Italy and Spain — [ECB WP 1707 (2014)](https://www.ecb.europa.eu/pub/pdf/scpwps/ecbwp1707.pdf); published in [IJCB 12(3), Sep 2016](https://www.ijcb.org/journal/ijcb16q3a1.htm); [CEPR DP 10025](https://cepr.org/publications/DP10025)
- The CJEU found OMT compatible with EU law (Gauweiler, 16 Jun 2015) — [CJEU press release 70/15](https://curia.europa.eu/site/upload/docs/application/pdf/2015-06/cp150070en.pdf)
- The Greek–German 10y spread fell from about 3,000 bp (May 2012) to under 1,000 bp by early 2013 — [BoG Annual Report 2012](https://www.bankofgreece.gr/Publications/Annrep2012.pdf) (snippet). This happened even though Greece was ineligible, which points to spillovers from the OMT, the Nov 2012 Eurogroup package and the Dec 2012 buyback acting together.

### Inferences
- In the model's terms, OMT is a conditional, unlimited backstop: it caps the priced default probability, and its effect appears without any purchases. Greece is the counterfactual case where the backstop was absent, because it had no market access.

### Gaps
- I did not retrieve exact announcement-day changes for Greek yields around 26 Jul, 2 Aug and 6 Sep 2012.

## 6. Cross-border funding-cost proxies for a debt-elastic premium on the national / TARGET2 position

### Takeaway
No direct market price exists for the TARGET2 position itself. Within the Eurosystem it is settled at the policy rate, which is general Eurosystem practice and not verified in this session. Usable proxies are:
- (a) the ELA penalty of 100–150 bp over MRO
- (b) the Greek–German gap in deposit rates (about 3.6 pp in 2012, household term deposits) against a TARGET2 liability of about EUR 105bn
- (c) the Greek–German sovereign spread (about 3,000 bp peak)

### Cited findings
- Greek minus German household term-deposit rates (annual mean, pp, computed from the ECB MIR series in section 2): 2010 2.15; 2011 2.77; 2012 3.61; 2013 2.94 — [ECB MIR GR](https://data-api.ecb.europa.eu/service/data/MIR/M.GR.B.L22.F.R.A.2250.EUR.N?format=csvdata), [ECB MIR DE](https://data-api.ecb.europa.eu/service/data/MIR/M.DE.B.L22.F.R.A.2250.EUR.N?format=csvdata)
- Greek TARGET2 liability (annual mean, EUR bn): 2010 79.5; 2011 94.2; 2012 105.0; 2013 62.0 — [ECB TGB](https://data-api.ecb.europa.eu/service/data/TGB/M.GR.N.A094T.U2.EUR.E?format=csvdata)
- ELA premium of 100–150 bp over MRO — [Yale JFC](https://elischolar.library.yale.edu/cgi/viewcontent.cgi?article=1328&context=journal-of-financial-crises)
- Before the crisis, Greece's NIIP went from −19.3% of GDP (1998) to −99.9% (2007), and its current-account deficit from 5.3% to 15.9% of GDP — [DNB / GreeSE paper](https://www.dnb.nl/media/ec2dqfwd/presentation-session-iii.pdf)

### Inferences
- A rough elasticity: between 2010 and 2012 the TARGET2 liability rose by about EUR 25bn (annual means) and the deposit-rate spread rose by about 1.5 pp. That is about 6 bp per EUR 1bn, or about 0.11 pp per 1% of Greek GDP if 2012 GDP is about EUR 190bn. The GDP figure is not sourced here and should be checked. This mixes sovereign-risk pass-through with the funding position, so treat it as an upper bound on a pure debt-elastic premium.
- A conservative alternative is to set the premium on the marginal unit of central-bank funding equal to the ELA penalty (1.0–1.5 pp), switched on only when collateral eligibility is lost.

### Gaps
- The remuneration of TARGET2 intra-Eurosystem balances at the MRO rate is not sourced here (ECB Guideline on TARGET2 / annual accounts).
- I found no academic estimate of a debt-elastic premium on TARGET2 positions specifically. Bruegel and CESifo TARGET2 trackers were not fetched.
