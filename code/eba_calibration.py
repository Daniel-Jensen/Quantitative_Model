"""EBA 2011 EU-wide stress-test -> model calibration moments.

Reads the EBA 2011 disclosure database (``data/DATA_DISCLOSURE.CSV``) and Eurostat
nominal GDP, decodes them via an explicit, validated code map, and emits
``data/eba_moments.json`` with every moment plus its provenance.

Exercise:  EBA 2011 EU-wide stress test (published 15 Jul 2011).
Base:      31 Dec 2010 ACTUAL (aggregate rows carry a blank SCENARIO_CODE; the
           2011/2012 rows are baseline=100 / adverse=105 *projections*).
Sovereign: gross direct long exposures as-of 31 Dec 2010 (worksheet 5; the
           INFO_DATE_CODE=20111231 field is a reference-year label per EBA method).
Country map: D = Greece (banks GR030-GR035, issuer GR),
             F = Germany (banks DE017-DE029, issuer DE).

Code map (validated: DE017 = Deutsche Bank, 2010 RWA 346.6bn / CT1 30.4bn / 8.76%):
    worksheet "1 - Aggregate information", MEASURE 120:
        30010 = total risk-weighted assets (RWA)
        30014 = Core Tier 1 capital
        30029 = total assets
    worksheet "5 - Sovereign exposures", INFORMATION_CODE 34010, MATURITY 999:
        gross direct long sovereign exposure (total across maturities)

Run:  /opt/anaconda3/envs/ssj/bin/python code/eba_calibration.py
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

# ── Paths ────────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
DATA_DIR = os.path.join(_ROOT, "data")
CSV_PATH = os.path.join(DATA_DIR, "DATA_DISCLOSURE.CSV (1).csv")
GDP_PATH = os.path.join(DATA_DIR, "eurostat_gdp_raw.json")
OUT_PATH = os.path.join(DATA_DIR, "eba_moments.json")

# ── Country / bank map ───────────────────────────────────────────────────────
GR_BANKS = ["GR030", "GR031", "GR032", "GR033", "GR034", "GR035"]
# German banks: all DE-prefixed disclosure columns.
DATE_2010 = 20101231.0

WS_AGG = "1 - Aggregate information"
WS_SOV = "5 - Sovereign exposures"
CODE_RWA, CODE_CT1, CODE_TA = 30010, 30014, 30029
CODE_SOV_GROSS_LONG, MATURITY_TOTAL = 34010, 999.0


def load_eba(csv_path: str = CSV_PATH) -> pd.DataFrame:
    """Load the EBA disclosure CSV (latin-1: it contains non-UTF-8 bytes)."""
    return pd.read_csv(csv_path, low_memory=False, encoding="latin-1")


def _to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _de_banks(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith("DE") and c[2:].isdigit()]


def group_sum(df, worksheet, info_code, date, banks, *, scenario="actual",
              maturity=None):
    """Sum one disclosure line item across a group of bank columns.

    scenario: "actual" -> SCENARIO_CODE is blank (the 2010 base year);
              100/105  -> baseline/adverse projection.
    """
    w = df[(df["WORKSHEET"] == worksheet)
           & (df["INFORMATION_CODE"] == info_code)
           & (df["INFO_DATE_CODE"] == date)]
    if scenario == "actual":
        w = w[w["SCENARIO_CODE"].isna()]
    else:
        w = w[w["SCENARIO_CODE"] == scenario]
    if maturity is not None:
        w = w[w["MATURITY_CODE"] == maturity]
    if len(w) == 0:
        return float("nan")
    return float(_to_num(w[banks].iloc[0]).sum())


def sovereign_book(df, issuer, banks):
    """Gross direct long sovereign exposure (as-of 31 Dec 2010), EUR m."""
    w = df[(df["WORKSHEET"] == WS_SOV)
           & (df["C_COUNTRY_CODE"] == issuer)
           & (df["INFORMATION_CODE"] == CODE_SOV_GROSS_LONG)
           & (df["MATURITY_CODE"] == MATURITY_TOTAL)]
    return float(_to_num(w[banks].iloc[0]).sum())


def load_gdp(gdp_path: str = GDP_PATH) -> dict:
    """Parse Eurostat JSON-stat -> {geo: {year: annual_nominal_GDP_EURm}}."""
    with open(gdp_path) as fh:
        js = json.load(fh)
    geo_idx = js["dimension"]["geo"]["category"]["index"]   # e.g. {"DE":0,"EL":1}
    time_idx = js["dimension"]["time"]["category"]["index"]  # {"2009":0,...}
    n_time = len(time_idx)
    values = js["value"]
    out: dict = {}
    for geo, gi in geo_idx.items():
        out[geo] = {}
        for yr, ti in time_idx.items():
            flat = gi * n_time + ti           # size [...,geo=2,time=n]; row-major
            v = values.get(str(flat))
            if v is not None:
                out[geo][yr] = float(v)
    return out


def compute_moments(df=None, gdp=None) -> dict:
    """Compute the full moment set with provenance. Returns a JSON-able dict."""
    if df is None:
        df = load_eba()
    if gdp is None:
        gdp = load_gdp()
    de_banks = _de_banks(df)

    # Bank aggregates (2010 actual base).
    ct1_D = group_sum(df, WS_AGG, CODE_CT1, DATE_2010, GR_BANKS)
    ct1_F = group_sum(df, WS_AGG, CODE_CT1, DATE_2010, de_banks)
    rwa_D = group_sum(df, WS_AGG, CODE_RWA, DATE_2010, GR_BANKS)
    rwa_F = group_sum(df, WS_AGG, CODE_RWA, DATE_2010, de_banks)
    ta_D = group_sum(df, WS_AGG, CODE_TA, DATE_2010, GR_BANKS)
    ta_F = group_sum(df, WS_AGG, CODE_TA, DATE_2010, de_banks)

    # Cross-holding matrix (gross direct long, 31 Dec 2010).
    b_D_D = sovereign_book(df, "GR", GR_BANKS)   # GR banks' Greek book
    b_D_F = sovereign_book(df, "GR", de_banks)   # DE banks' Greek book (contagion)
    b_F_D = sovereign_book(df, "DE", GR_BANKS)   # GR banks' Bund book
    b_F_F = sovereign_book(df, "DE", de_banks)   # DE banks' Bund book

    # Nominal GDP (2010, annual); quarterly = /4. Eurostat geo: EL=Greece.
    gdp_ann_D = gdp["EL"]["2010"]
    gdp_ann_F = gdp["DE"]["2010"]
    qgdp_D = gdp_ann_D / 4.0
    qgdp_F = gdp_ann_F / 4.0

    # Model targets (symmetric Y=1 normalization; ratios to own quarterly GDP).
    targets = {
        "n_inter_D": ct1_D / qgdp_D,
        "n_inter_F": ct1_F / qgdp_F,
        # sovereign book / capital (unit-free; the doom-loop vulnerability)
        "phi_bD_D_ss": b_D_D / ct1_D,   # GR banks' own-sovereign / capital
        "phi_bF_F_ss": b_F_F / ct1_F,   # DE banks' own-sovereign / capital
        "phi_bD_F_ss": b_D_F / ct1_F,   # DE banks' GR / capital (cross)
        "phi_bF_D_ss": b_F_D / ct1_D,   # GR banks' Bund / capital (cross)
        "leverage_D": ta_D / ct1_D,
        "leverage_F": ta_F / ct1_F,
        # bank-held gov debt as fraction of annual own GDP (decision b1: B_supply)
        "bankheld_gov_D_pct_annGDP": (b_D_D + b_D_F) / gdp_ann_D,
        "bankheld_gov_F_pct_annGDP": (b_F_F + b_F_D) / gdp_ann_F,
        "B_supply_D_qgdp": (b_D_D + b_D_F) / qgdp_D,
        "B_supply_F_qgdp": (b_F_F + b_F_D) / qgdp_F,
    }

    return {
        "meta": {
            "source_exercise": "EBA 2011 EU-wide stress test",
            "base_date": "2010-12-31 (actual; scenario blank)",
            "sovereign_asof": "2010-12-31 gross direct long (34010, maturity total)",
            "gdp_source": "Eurostat nama_10_gdp B1GQ CP_MEUR 2010",
            "country_map": {"D": "Greece (GR030-035, issuer GR)",
                            "F": "Germany (DE0xx, issuer DE)"},
            "de_banks": de_banks,
            "units": "EUR million unless noted",
        },
        "raw_EURm": {
            "CT1_D": ct1_D, "CT1_F": ct1_F,
            "RWA_D": rwa_D, "RWA_F": rwa_F,
            "TA_D": ta_D, "TA_F": ta_F,
            "b_D_D": b_D_D, "b_D_F": b_D_F, "b_F_D": b_F_D, "b_F_F": b_F_F,
            "GDP_ann_D": gdp_ann_D, "GDP_ann_F": gdp_ann_F,
        },
        "model_targets": targets,
    }


def validate(df=None) -> None:
    """Self-check: DE017 must reproduce Deutsche Bank's 2010 figures."""
    if df is None:
        df = load_eba()
    rwa = group_sum(df, WS_AGG, CODE_RWA, DATE_2010, ["DE017"])
    ct1 = group_sum(df, WS_AGG, CODE_CT1, DATE_2010, ["DE017"])
    ratio = ct1 / rwa
    assert abs(rwa - 346_608) < 5, f"DE017 RWA {rwa:,.0f} != Deutsche Bank 346,608"
    assert abs(ratio - 0.0876) < 0.0010, f"DE017 CT1 ratio {ratio:.4f} != 0.0876"


def main() -> dict:
    df = load_eba()
    validate(df)
    moments = compute_moments(df)
    with open(OUT_PATH, "w") as fh:
        json.dump(moments, fh, indent=2)
    t = moments["model_targets"]
    print(f"wrote {OUT_PATH}")
    print("── EBA 2010-base model targets ──")
    for k, v in t.items():
        print(f"  {k:28} {v:12.4f}")
    return moments


if __name__ == "__main__":
    main()
