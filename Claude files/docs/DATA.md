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

1. Since 2026-10-27 the derived tables in `Empirics/outputs/` (json, csv) are tracked in git;
   figures are regenerated, and the quarterly CDS series (`greek_cds_series.csv`) stays out
   because it derives from the licensed LSEG extract. Raw downloads under `data/` and the CDS
   extract remain untracked.
2. The code that produced the EBA, leverage, wedge and CDS numbers was recovered from the
   shelved branch and re-run on 2026-10-27: every output is identical to the stored one.
3. The chat-computed moments are now scripted (`Empirics/validation_moments.py`, section 4),
   except coupon decay, the import share and the basis-adjusted hazard.
4. Two scripts in the tree still hard-code the old model calibration in their model-side
   columns (section 6).

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
| All-foreign share 39.6%; EBA sample holds 27% of Greek debt (EUR 90.1bn of 330.6bn); holder shares GR 60.4, FR 11.2, DE 8.8, CY 6.4, BE 4.3, other 8.9% | Same panel | `eba_discipline.json` (`stocks`), `eba_holder_groups.csv` | Recovered; reproduces |
| Greek sovereign / CT1 2.39 (EUR 54.4bn / 22.8bn) | Six Greek banks | `eba_discipline.json`, `eba_bank_panel.csv` | Recovered; reproduces |
| Repricing loss 40% of Greek CT1; price fall 21.6% | Frozen 31 Dec 2010 maturity ladder revalued from the Dec 2010 to the Dec 2011 ten-year yield (12.01% to 21.14%), flat curve | `eba_discipline.json` (`shock`) | Recovered; reproduces (0.397, 0.216) |
| Leverage 7.98 | Greek MFIs' total financial assets over equity, quarterly financial accounts, mean 2002Q2-2007Q2 (21 quarters), on Bocola's definition | `bocola_greece/moments.py`, `bocola_greece_moments.json` (`leverage`) | Recovered; reproduces |
| Wedge 5 bp (and 8 bp over 2002-2012) | Mean one-week Euribor minus Eurepo, annualised, 2002Q2-2009Q4 (Bocola's estimation sample) | same, from Bocola's `Data.xls` (`wedge_bp_pa`) | Recovered; reproduces (5.10, 8.23). Greek-bank variants (2-12 bp, 29 bp) in the same file |
| Pre-crisis Greek-German ten-year spread 25.3 bp | Mean 2002Q2-2007Q2 | `greek_cds_bocola.json`; also `spreads_fred.csv` (25.2 on monthly averages) | Recovered / in tree |
| Ten-year spread 306 bp in 2010Q1 | Quarterly mean of the monthly spread | `spreads_fred.csv` | In tree; reproduces |
| Risk process: ρ_s and σ_s table (full sample 0.985 / 0.24, tranquil 0.75 / 0.19, CDS 0.75 / 0.65, ...), rescaling ratio 1.74, σ_s = 0.41 | OLS of the logit of the quarterly default probability implied by spreads at LGD 0.55; quarterly averages of one-year CDS and the ten-year spread | `greek_cds/risk_process.py`, `risk_process_estimates.json` | In tree; reproduces. Bias correction to 0.89 (and ρ_s = 0.9) is hand algebra in appendix B |
| Default probability at rest, 0.102% (Greek one-year CDS, 2008Q3) | Mean of daily one-year CDS from the first quote (21 Jul 2008) to end-September, 22.4 bp, as a quarterly hazard at LGD 0.55 | `cds_panel.csv`; `greek_cds_series.csv` | Reproduces. The paper's former 0.092% was not reproducible and was corrected on Overleaf on 2026-10-27 |
| Greek CDS 1.5-2.3%/qtr between Oct 2009 and 23 Apr 2010; above 4% from May 2010 | Monthly averages of the CDS-implied hazard | `cds_panel.csv` | No script |
| Bohn coefficients (1995-2008 -0.10, 1995-2019 0.01, 2001-2019 0.06, 2010-2019 0.15) | Primary balance on lagged debt and the HP(100) output gap, annual, Newey-West 2 lags | `greek_fiscal/bohn.py`, `bohn_estimates.json` | In tree; coefficients reproduce. Implied roots in the json use the old calibration (section 6) |
| Lending spreads 280 bp (GR), 176 bp (DE), 2003-2007 | Composite cost of borrowing minus three-month Euribor, monthly mean Jan 2003 to Dec 2007 | `lending_rates.csv` | Reproduces (279.5, 175.9). The paper's former 290 / 180 corrected on Overleaf, 2026-10-27 |
| Overdraft-spread rise 108 bp | Greek minus German overdraft rate, crisis (2010-2012H1) minus pre-crisis | `eba_descriptives.py` (section 3 of its console report), `lending_rates.csv` | Reproduces: 84.6 bp (2005-2007) to 192.2 bp (2010-2012H1), a rise of 107.5 |
| Investment fell 69.7% and GDP 27.4%, 2007Q4 to 2015Q3; PSI cut EUR 75bn and 33.7 points of debt ratio | Eurostat series above | `motivation_figures.py`, `greece_macro.csv` | In tree (output ignored) |
| Duration of Greek banks' Greek book 4.4 years | Macaulay duration of the six Greek banks' end-2010 maturity ladder at a 4.5% yield and 4.7% coupon | `calibration/eba_moments.py`: `ladder_value_and_duration(sovereign_ladder(df, "GR", GR_BANKS), 0.045, COUPON_GGB)` | Reproduces (4.36). At the end-2010 Greek yield of 12.0% the same ladder gives 3.5 years (`eba_discipline.json`); the gap is the discount rate |
| PSI recovery 46.5% face value, 35-41% NPV | Zettelmeyer, Trebesch and Gulati (2013) | literature | Cited |

## 4. Numbers once computed only in chat

`Empirics/validation_moments.py` (2026-10-27) now computes these from Eurostat, the ECB, FRED and
the stored EBA, IMF and CDS outputs, caching the downloads in `validation_inputs.json` and
writing `validation_moments.json` (paper value, reproduced value, construction).

| Number | Construction | Status |
|---|---|---|
| Government purchases / GDP 0.204 (GR), 0.191 (DE) | General-government final consumption / GDP, mean 1999-2010 | Reproduces (0.2041, 0.1915). Window recovered by search; the original was not recorded |
| Relative size 11.3 | German / Greek nominal GDP, 2006 | Reproduces (11.32). Year recovered; the paper's placeholder can read 2006 |
| Liquid deposits 84% of Greek GDP | Households' and NPISHs' currency and deposits (F2) / GDP, 2008 | Reproduces (84.6%) |
| Greek bank assets 1.7 times GDP | Total assets of the six Greek EBA banks / GDP, 2010 | Reproduces (1.69) |
| EBA-held Greek debt 40% of GDP | EBA sample's Greek sovereign / GDP, 2010 | Reproduces (40.3%) |
| Debt 164% (GR), 80% (DE) of GDP, 2012 | Eurostat `gov_10dd_edpt1` | Reproduces on the current vintage (164.1, 79.8); the paper's former 162% corrected 2026-10-27 |
| Greek primary surpluses 3.5-4.2% of GDP, 2016-2018 | IMF DataMapper | Reproduces (3.50, 4.19, 4.14) |
| Bond-implied hazard 0.23%, 1999-2007 | Mean GR-DE ten-year spread (49.8 bp) as a quarterly hazard at LGD 0.55 | Reproduces (0.226%) |
| Default probability 9.1%, 2010-2011 | Greek one-year CDS through Bocola's conversion, LGD 0.55 (`greek_cds/greek.py`) | Reproduces (9.09%) |
| Real deposit rate 1.15%, 1999-2007 | Euro-area three-month Euribor minus HICP inflation | Reproduces the corrected value; the paper's former 1.11% had no recorded construction |
| Home bias 0.68 (import share 32%, 2012) | Greek imports / GDP is 33.4% in 2012 and 31.0% in 2011; imports / final consumption 36.3% | Not reproduced |
| Bond-implied hazard 0.02-0.05%, basis-adjusted | Construction not recovered | Not reproduced |
| Coupon decay 0.056 (GR), 0.039 (DE) | Eurostat government-securities statistics (`eurostat_irt`); 0.056 is Bocola's Italian value | No script |
| Greek CDS 1.5-2.3%/qtr, Oct 2009 to 23 Apr 2010; above 4% from May 2010 | Monthly averages of the CDS-implied hazard | No script |

## 5. Data behind the SSJ model (documented elsewhere)

- EBA moments for the SSJ calibration: `calibration/eba_moments.py` reads the EBA disclosure
  database, Eurostat GDP and FRED yields and writes `data/eba_moments.json` (tracked), with
  an identification ledger. Derivation: `Claude files/docs/eba_calibration.md`.
- Default persistence `rho_def` = 0.9408: three-state Markov-switching model on peripheral
  spreads, `Empirics/ms_regime.py`, outputs `ms_regime_*` (tracked).
- Fiscal limit `def_scale_D`: `Empirics/fiscal_limit.py` (Eurostat quarterly debt and FRED
  spreads, credit triangle at the SSJ recovery 0.30).

## 6. Known discrepancies

- `risk_process.py` uses ρ_s = 0.95 and σ_s = 0.4455 for its model-side columns (share of
  quarters above 2%, 5%, 9.1%), so those columns describe the old process, not the paper's.
- `bohn.py` uses q̄ = 0.91057 and a debt root of 0.93, so the json's implied roots
  (1.034, 1.002, 0.989, 0.964) and its "model" block differ from appendix B's (1.030, 0.9997,
  0.988, 0.963 at q̄ = 0.947, Bohn 0.04).
- Re-running `greek_cds/greek.py` on top of its own cached series appends duplicate columns
  (`*.1`) to `greek_cds_series.csv`, because the script concatenates fresh frames onto the file
  it reads. Values are unaffected; restore the file or run with a clean cache.

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
python Empirics/greek_cds/greek.py                 # CDS -> default probability (Bocola's conversion)
python Empirics/validation_moments.py              # chat-era moments (section 4)
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
