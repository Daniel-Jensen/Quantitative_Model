# Data

What was computed from data, from which source, by which script, and where the result is
stored: the record needed to replicate the calibration and, later, to estimate the model.
Checked against the files on 2026-10-27. Values themselves are reported in the paper
(appendix B) and in the outputs named here; this file records provenance, not results.

Status legend used below:

- **In tree**: the script and its output are on the working branch and the number reproduces.
- **Recovered**: the script was shelved in commit `c2d2337` on branch `cross-border-empirics`
  (a work-in-progress snapshot of 1 Oct 2026) and restored to `main` on 2026-10-27. The output
  file is in `Empirics/outputs/`, which is git-ignored.
- **Untracked**: the script is in the working tree but not committed.
- **No script**: the number was computed in a chat session and never saved as code.

## 1. The replication risk in brief

1. `Empirics/outputs/` is git-ignored (only the Markov-switching outputs and `spreads_fred.csv`
   are committed). The files most of the paper's calibration rests on (EBA exposures, Greek
   leverage, the Euribor-Eurepo wedge, CDS series, risk-process and Bohn estimates, lending
   rates) exist only on this laptop.
2. The code that produced the EBA, leverage, wedge and CDS numbers was recovered from the
   shelved branch on 2026-10-27; only `bohn.py` and `risk_process.py` have been re-run since.
3. Some calibration inputs (government purchases, coupon decay, import share, relative size,
   and several validation moments) have no script anywhere (section 4).

## 2. Raw sources

| Source | Identifier | Coverage | Local file | In git | Fetched by |
|---|---|---|---|---|---|
| EBA 2011 EU-wide stress test, disclosure database | Bank-level sovereign exposures by country and residual maturity; CT1, RWA, total assets | 90 banks, 31 Dec 2010 | `data/DATA_DISCLOSURE.CSV (1).csv`; methodology PDFs beside it | No (raw data ignored) | Manual download, EBA website |
| EBA 2010 stress test, 2011 and 2012 capital exercises | Sovereign exposure by country and maturity | 2010-2012 | expected in `Empirics/cross_border/data_raw/`, absent | No | `fetch_eba.py`; see the README in that folder |
| Bocola (2016) replication archive | `Data.xls`: Italian bank leverage, one-week Euribor and Eurepo, his multiplier series; `model_param.m` | 1999-2012 | `data/bocola2016/` | Yes | Journal archive |
| Eurostat quarterly financial accounts | `nasq_10_f_bs`, sectors S121+S122+S123, total financial assets and equity liabilities (F51) | GR, IT, quarterly | cached in `bocola_greece_series.csv` | No | `bocola_greece/moments.py` |
| Eurostat government debt | `gov_10q_ggdebt` (S13, GD) quarterly; `gov_10dd_edpt1` annual (Greek debt 2010: EUR 330.57bn, hard-coded) | GR | `greece_macro.csv` | No | `motivation_figures.py`, `fiscal_limit.py` |
| Eurostat national accounts | `namq_10_gdp` B1GQ, P51G, chain-linked 2010, SCA | GR quarterly | `greece_macro.csv` | No | `motivation_figures.py` |
| Eurostat annual GDP | nominal GDP | EBA countries | `data/eurostat_gdp_raw.json` | No | `calibration/eba_moments.py` |
| FRED 10-year yields | `IRLTLT01{IT,ES,GR,PT,IE,FR,DE}M156N`; recession dummies `{ITA,ESP,GRC,PRT}RECM` | monthly, 1995-2026 | `Empirics/outputs/spreads_fred.csv` | Yes | `graph_spreads.py`, `ms_regime.py` |
| ECB MFI interest rates | `MIR/M.{GR,DE}.B.A2I.AM.R.A.2240.EUR.N` (composite cost of borrowing, NFCs); `MIR/M.{GR,DE}.B.A2Z.A.R.A.2240.EUR.N` (overdrafts) | monthly, 2003-2015 | `lending_rates.csv` | No | `eba_descriptives.py` |
| ECB financial market data | `FM/M.U2.EUR.RT.MM.EURIBOR3MD_.HSTA` (3-month Euribor) | monthly | `lending_rates.csv` | No | `eba_descriptives.py` |
| LSEG Workspace sovereign CDS | Composite par spreads, `<CTY>GV<TENOR><CCY>AB=R` (2003 ISDA cum-restructuring), GR and IT, 6M-10Y | daily, Jul 2008 to the March 2012 credit event | `Empirics/greek_cds/data/cds_panel.csv` | No, and never: licensed | `fetch_lseg.py`; needs the Workspace desktop app and an app key kept outside the repository |
| IMF DataMapper | `GGXONLB_G01_GDP_PT` (primary balance), `GGXWDG_NGDP` (gross debt), `NGDP_RPCH` (real growth) | GR annual; downloaded 2026-10-13 | `Empirics/greek_fiscal/data/imf_greece.json` | No | written during a session; `bohn.py` reads it |

## 3. Numbers in the paper and where they come from

Paper calibration (section 4, appendix B). "Output" is in `Empirics/outputs/` unless stated.

| Number | Construction | Script, output | Status |
|---|---|---|---|
| Own sovereign / bank assets 14.4% (GR), 6.5% (DE) | Gross direct long exposure to own sovereign over total assets; 6 Greek and 12 German EBA banks, 31 Dec 2010 | `eba_descriptives.py`, `eba_bank_panel.csv` (`sov_own / ta`) | Recovered; reproduces (0.1443, 0.0647). Net of shorts, Greek share 12.8% (`bocola_greece_moments.json`) |
| Foreign share of D debt 12.7% | German banks' Greek holdings over Greek plus German banks' Greek holdings (7.93 / 62.38 EUR bn) | `eba_descriptives.py`, `eba_holder_groups.csv` | Recovered; reproduces (0.1272) |
| D banks' F bonds / assets 0.11% (code only; not in the paper) | Six Greek banks' gross direct long Bund exposure over their total assets (EUR 0.41bn / 377.2bn) | `calibration/eba_moments.py` reader (`sovereign_book(df, "DE", GR_BANKS)`), computed in session 2026-10-27; no stored output | In tree reader; no script for the ratio |
| All-foreign share 39.6%; EBA sample holds 27% of Greek debt (EUR 90.1bn of 330.6bn); holder shares GR 60.4, FR 11.2, DE 8.8, CY 6.4, BE 4.3, other 8.9% | Same panel | `eba_discipline.json` (`stocks`), `eba_holder_groups.csv` | Recovered; reproduces |
| Greek sovereign / CT1 2.39 (EUR 54.4bn / 22.8bn) | Six Greek banks | `eba_discipline.json`, `eba_bank_panel.csv` | Recovered; reproduces |
| Repricing loss 40% of Greek CT1; price fall 21.6% | Frozen 31 Dec 2010 maturity ladder revalued from the Dec 2010 to the Dec 2011 ten-year yield (12.01% to 21.14%), flat curve | `eba_discipline.json` (`shock`) | Recovered; reproduces (0.397, 0.216) |
| Leverage 7.98 | Greek MFIs' total financial assets over equity, quarterly financial accounts, mean 2002Q2-2007Q2 (21 quarters), on Bocola's definition | `bocola_greece/moments.py`, `bocola_greece_moments.json` (`leverage`) | Recovered; reproduces |
| Wedge 5 bp (and 8 bp over 2002-2012) | Mean one-week Euribor minus Eurepo, annualised, 2002Q2-2009Q4 (Bocola's estimation sample) | same, from Bocola's `Data.xls` (`wedge_bp_pa`) | Recovered; reproduces (5.10, 8.23). Greek-bank variants (2-12 bp, 29 bp) in the same file |
| Pre-crisis Greek-German ten-year spread 25.3 bp | Mean 2002Q2-2007Q2 | `greek_cds_bocola.json`; also `spreads_fred.csv` (25.2 on monthly averages) | Recovered / in tree |
| Ten-year spread 306 bp in 2010Q1 | Quarterly mean of the monthly spread | `spreads_fred.csv` | In tree; reproduces |
| Risk process: ρ_s and σ_s table (full sample 0.985 / 0.24, tranquil 0.75 / 0.19, CDS 0.75 / 0.65, ...), rescaling ratio 1.74, σ_s = 0.41 | OLS of the logit of the quarterly default probability implied by spreads at LGD 0.55; quarterly averages of one-year CDS and the ten-year spread | `greek_cds/risk_process.py`, `risk_process_estimates.json` | In tree; reproduces (re-run 2026-10-27 at ρ_s = 0.9, σ_s = 0.41). Bias correction to 0.89 (and ρ_s = 0.9) is hand algebra in appendix B |
| Default probability at rest, 0.092% (Greek one-year CDS, 2008Q3) | CDS-implied quarterly hazard at LGD 0.55 | not located | No script. The stored quarterly series (`greek_cds_series.csv`) gives 0.102% |
| Greek CDS 1.5-2.3%/qtr between Oct 2009 and 23 Apr 2010; above 4% from May 2010 | Monthly averages of the CDS-implied hazard | `cds_panel.csv` | No script |
| Bohn coefficients (1995-2008 -0.10, 1995-2019 0.01, 2001-2019 0.06, 2010-2019 0.15) | Primary balance on lagged debt and the HP(100) output gap, annual, Newey-West 2 lags | `greek_fiscal/bohn.py`, `bohn_estimates.json` | In tree; coefficients and appendix B's implied roots reproduce (re-run 2026-10-27 at q̄ = 0.947, Bohn 0.04) |
| Lending spreads 290 bp (GR), 180 bp (DE), 2003-2007 | Composite cost of borrowing minus three-month Euribor | `lending_rates.csv` | Not reproduced: the obvious window gives 280 / 176 (section 6) |
| Overdraft-spread rise 108 bp | Greek minus German overdraft rate, crisis (2010-2012H1) minus pre-crisis | `lending_rates.csv`; the `eba_descriptives.py` caption reports 84 to 193 bp | Pre-crisis window not recorded; 2005-2007 vs 2010-2012 gives 124 |
| Investment fell 69.7% and GDP 27.4%, 2007Q4 to 2015Q3; PSI cut EUR 75bn and 33.7 points of debt ratio | Eurostat series above | `motivation_figures.py`, `greece_macro.csv` | In tree (output ignored) |
| Duration of Greek banks' Greek book 4.4 years | Macaulay duration at a 4.5% discount rate | not located | No script. The stored ladder gives 3.5 years at the end-2010 Greek yield (`eba_discipline.json`) and an average maturity of 25.6 quarters (`bocola_greece_moments.json`) |
| PSI recovery 46.5% face value, 35-41% NPV | Zettelmeyer, Trebesch and Gulati (2013) | literature | Cited |

## 4. Numbers with no script

Computed or looked up during chat sessions; the source series is named in the paper but no
code or cached file records the calculation. Each needs a script before it can be replicated.

- Government purchases / GDP, 0.204 (GR) and 0.191 (DE): Eurostat `nama_10_gdp`.
- Coupon decay 0.056 (GR) and 0.039 (DE): Eurostat government-securities statistics
  (`eurostat_irt` in the bibliography). 0.056 coincides with Bocola's Italian value.
- Home bias 0.68: Greek import share of consumption 32% in 2012.
- Relative size 11.3: German over Greek nominal GDP; the year is still a placeholder.
- Real deposit rate 1.11%, 1999-2007 (validation of Bocola's 1.2%).
- Bond-implied hazards 0.23% (1999-2007) and 0.02-0.05% (basis-adjusted).
- Liquid deposits 84% of Greek GDP.
- Greek bank assets 1.7 times GDP; EBA-held Greek debt 40% of Greek GDP.
- Debt outstanding 162% (GR) and 80% (DE) of GDP in 2012: Eurostat `gov_10dd`.
- Greek primary surpluses of 3.5-4.2% of GDP in 2016-2018.
- Default probability 9.1% averaged over 2010-2011 (the exceedance benchmark).

## 5. Data behind the SSJ model (documented elsewhere)

- EBA moments for the SSJ calibration: `calibration/eba_moments.py` reads the EBA disclosure
  database, Eurostat GDP and FRED yields and writes `data/eba_moments.json` (tracked), with
  an identification ledger. Derivation: `Claude files/docs/eba_calibration.md`.
- Default persistence `rho_def` = 0.9408: three-state Markov-switching model on peripheral
  spreads, `Empirics/ms_regime.py`, outputs `ms_regime_*` (tracked).
- Fiscal limit `def_scale_D`: `Empirics/fiscal_limit.py` (Eurostat quarterly debt and FRED
  spreads, credit triangle at the SSJ recovery 0.30).

## 6. Known discrepancies

- Lending spreads, overdraft rise, the 2008Q3 default probability and the 4.4-year duration do
  not reproduce from the stored files (section 3).

## 7. Rebuilding

Run from the repository root with `/opt/anaconda3/envs/ssj/bin/python` (the base environment's
linear algebra library is broken).

```bash
python Empirics/graph_spreads.py                   # FRED yields -> spreads_fred.csv
python Empirics/motivation_figures.py --refresh    # Eurostat debt and activity
python Empirics/eba_descriptives.py                # EBA panel, exposures, repricing, MIR rates
python Empirics/bocola_greece/moments.py           # leverage, Euribor-Eurepo, multiplier
python Empirics/greek_cds/risk_process.py          # risk-process table
python Empirics/greek_fiscal/bohn.py               # Bohn regressions
python calibration/eba_moments.py                  # SSJ EBA moments
```

The CDS pull (`fetch_lseg.py`) needs LSEG Workspace running locally and
`LSEG_APP_KEY=$(cat ~/.config/lseg/appkey)`; the extract is licensed and stays out of git.

## 8. Toward estimation

Quarterly observables already assembled, with their files:

| Observable | Coverage | File |
|---|---|---|
| Greek-German and other periphery ten-year spreads | 1995-2026 | `spreads_fred.csv` |
| Greek one-year (and 6M-10Y) CDS | 2008Q3-2012Q1 | `cds_panel.csv`, `greek_cds_series.csv` |
| Greek and Italian MFI leverage, Euribor-Eurepo, Bocola's multiplier | 1995Q1-2026Q1 (Euribor-Eurepo and the multiplier only to 2012, from Bocola's archive) | `bocola_greece_series.csv` |
| Greek real GDP, investment, government debt | 1995Q1-2026Q1 | `greece_macro.csv` |
| Greek and German lending and overdraft rates, 3-month Euribor | 2003-2015 | `lending_rates.csv` |
| Greek primary balance, debt, growth (annual) | 1980 (primary balance 1990) to 2031; IMF projections after 2025 | `imf_greece.json` |
| Bank-level sovereign exposures | end-2010 snapshot | `eba_bank_panel.csv` |

The shelved branch also holds an estimation design, including a data-requirements note:
`diagnostics/bayesian_estimation/DATA_REQUIREMENTS.md` in commit `c2d2337`.
