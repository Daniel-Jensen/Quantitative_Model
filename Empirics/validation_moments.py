# THE CALIBRATION APPENDIX'S VALIDATION MOMENTS AND SET-DIRECTLY INPUTS THAT WERE COMPUTED ONLY IN CHAT.
# Eurostat national and financial accounts, ECB rates, FRED spreads and the stored EBA, IMF and CDS
# outputs. Each entry records the paper's value, the value reproduced here and its construction.
# Fetched series are cached in Empirics/outputs/validation_inputs.json; --refresh re-downloads them.
# Writes Empirics/outputs/validation_moments.json.
import argparse
import io
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Empirics/outputs"
CACHE = OUT / "validation_inputs.json"
LGD = 0.55  # loss given default, the model's 1 - recovery
EUROSTAT = "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/"
ECB = "https://data-api.ecb.europa.eu/service/data/"

# every fetched series: name -> (source, dataset or flow, filters)
SERIES = {
    "gdp_EL": ("es", "nama_10_gdp", dict(geo="EL", unit="CP_MEUR", na_item="B1GQ")),
    "gdp_DE": ("es", "nama_10_gdp", dict(geo="DE", unit="CP_MEUR", na_item="B1GQ")),
    "gcons_EL": ("es", "nama_10_gdp", dict(geo="EL", unit="CP_MEUR", na_item="P3_S13")),
    "gcons_DE": ("es", "nama_10_gdp", dict(geo="DE", unit="CP_MEUR", na_item="P3_S13")),
    "imports_EL": ("es", "nama_10_gdp", dict(geo="EL", unit="CP_MEUR", na_item="P7")),
    "fincons_EL": ("es", "nama_10_gdp", dict(geo="EL", unit="CP_MEUR", na_item="P3")),
    "debt_pc_EL": ("es", "gov_10dd_edpt1", dict(geo="EL", unit="PC_GDP", sector="S13", na_item="GD")),
    "debt_pc_DE": ("es", "gov_10dd_edpt1", dict(geo="DE", unit="PC_GDP", sector="S13", na_item="GD")),
    "deposits_hh_EL": ("es", "nasa_10_f_bs", dict(geo="EL", unit="MIO_EUR", sector="S14_S15",
                                                  finpos="ASS", na_item="F2", co_nco="NCO")),
    "euribor3m": ("ecb", "FM", "M.U2.EUR.RT.MM.EURIBOR3MD_.HSTA"),
    "hicp_U2": ("ecb", "ICP", "M.U2.N.000000.4.ANR"),
}


def eurostat(dataset, filters):
    url = EUROSTAT + dataset + "?format=JSON&lang=EN&" + "&".join(f"{k}={v}" for k, v in filters.items())
    with urllib.request.urlopen(url, timeout=90) as r:
        d = json.load(r)
    inv = {v: k for k, v in d["dimension"]["time"]["category"]["index"].items()}
    return {inv[int(i)]: float(v) for i, v in d["value"].items()}


def ecb(flow, key):
    url = f"{ECB}{flow}/{key}?format=csvdata&startPeriod=1999-01&endPeriod=2008-12"
    with urllib.request.urlopen(url, timeout=90) as r:
        d = pd.read_csv(io.StringIO(r.read().decode()))
    return {str(k): float(v) for k, v in zip(d["TIME_PERIOD"], d["OBS_VALUE"])}


def inputs(refresh):
    if CACHE.exists() and not refresh:
        raw = json.loads(CACHE.read_text())
    else:
        raw = {k: eurostat(a, b) if src == "es" else ecb(a, b) for k, (src, a, b) in SERIES.items()}
        CACHE.write_text(json.dumps(raw, indent=1))
    return {k: pd.Series(v).sort_index() for k, v in raw.items()}


def annual(s):
    s = s.copy()
    s.index = s.index.astype(int)
    return s


def entry(paper, value, construction, reproduces):
    return dict(paper=paper, reproduced=value, construction=construction, reproduces=reproduces)


def main(refresh=False):
    x = inputs(refresh)
    gdp_EL, gdp_DE = annual(x["gdp_EL"]), annual(x["gdp_DE"])
    g_EL = (annual(x["gcons_EL"]) / gdp_EL).loc[1999:2010].mean()
    g_DE = (annual(x["gcons_DE"]) / gdp_DE).loc[1999:2010].mean()
    imp = annual(x["imports_EL"])
    eba = json.loads((OUT / "eba_discipline.json").read_text())["stocks"]
    panel = pd.read_csv(OUT / "eba_bank_panel.csv", index_col=0)
    ta_GR = panel.loc[panel.index.str.startswith("GR"), "ta"].sum() / 1e3  # EUR bn
    spread = pd.read_csv(OUT / "spreads_fred.csv", index_col=0, parse_dates=True)["GRC_spread"] * 100.0
    s_9907 = spread.loc["1999":"2007"].mean()
    rr = pd.concat([x["euribor3m"], x["hicp_U2"]], axis=1).dropna().loc["1999-01":"2007-12"]
    imf = json.loads((ROOT / "Empirics/greek_fiscal/data/imf_greece.json").read_text())
    pb = imf["series"]["GGXONLB_G01_GDP_PT"]["values"]
    cds = json.loads((OUT / "greek_cds_bocola.json").read_text())["greece_answer"]

    m = {
        # set directly (calibration inputs)
        "G_over_GDP": entry([0.204, 0.191], [g_EL, g_DE],
                            "general-government final consumption / GDP, current prices, mean 1999-2010 "
                            "(window recovered: the original was not recorded)", True),
        "relative_size_DE_EL": entry(11.3, float(gdp_DE[2006] / gdp_EL[2006]),
                                     "German / Greek nominal GDP, 2006 (year recovered)", True),
        "import_share_EL_2012": entry(0.32, dict(imports_over_gdp=float(imp[2012] / gdp_EL[2012]),
                                                 imports_over_final_consumption=float(
                                                     imp[2012] / annual(x["fincons_EL"])[2012]),
                                                 imports_over_gdp_2011=float(imp[2011] / gdp_EL[2011])),
                                      "Greek imports of goods and services, 2012; no national-accounts "
                                      "ratio gives 32%", False),
        # validation moments
        "real_deposit_rate_1999_2007": entry(1.11, float((rr.iloc[:, 0] - rr.iloc[:, 1]).mean()),
                                             "euro-area three-month Euribor minus HICP inflation (yoy), "
                                             "monthly mean, % a year", False),
        "bond_hazard_1999_2007": entry(0.23, float(100 * (1 - np.exp(-s_9907 / 1e4 / 4 / LGD))),
                                       "mean GR-DE ten-year spread 1999-2007 as a quarterly hazard at "
                                       "LGD 0.55, %", True),
        "bond_hazard_basis_adjusted": entry([0.02, 0.05], None,
                                            "construction not recovered", False),
        "liquid_deposits_over_gdp_EL": entry(0.84, float(annual(x["deposits_hh_EL"])[2008] / gdp_EL[2008]),
                                             "households' and NPISHs' currency and deposits (F2) / GDP, "
                                             "2008", True),
        "greek_bank_assets_over_gdp": entry(1.7, float(ta_GR / (gdp_EL[2010] / 1e3)),
                                            "total assets of the six Greek EBA banks / Greek GDP, 2010", True),
        "eba_held_greek_debt_over_gdp": entry(0.40, float(eba["GR_sov_held_by_sample_EURbn"]
                                                          / (gdp_EL[2010] / 1e3)),
                                              "Greek sovereign held by the EBA sample / Greek GDP, 2010", True),
        "debt_over_gdp_2012": entry([162, 80], [annual(x["debt_pc_EL"])[2012], annual(x["debt_pc_DE"])[2012]],
                                    "general-government gross debt, % of GDP, 2012 (current vintage)", "DE only"),
        "greek_primary_surplus_2016_2018": entry([3.5, 4.2], [pb[str(y)] for y in (2016, 2017, 2018)],
                                                 "IMF DataMapper primary balance, % of GDP", True),
        "default_probability_2010_2011": entry(9.1, cds["sovereign_crisis_2010Q1_2011Q4"]["p_phys_pct_q_D55"],
                                               "Greek one-year CDS, 2010Q1-2011Q4, Bocola's conversion at "
                                               "LGD 0.55, % a quarter (from greek_cds/greek.py)", True),
    }
    path = OUT / "validation_moments.json"
    path.write_text(json.dumps(m, indent=2, default=float))
    for k, e in m.items():
        print(f"{k:<34} paper {e['paper']!s:<14} reproduced {e['reproduced']!s:<60.60} {e['reproduces']}")
    print(f"wrote {path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true", help="re-download the Eurostat and ECB series")
    main(ap.parse_args().refresh)
