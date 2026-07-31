"""EBA 2011 EU-wide stress test -> model calibration moments (rebuilt 2026-07-31).

Reads the EBA 2011 disclosure database (``data/DATA_DISCLOSURE.CSV``), Eurostat
nominal GDP and FRED long-term sovereign yields, and emits ``data/eba_moments.json``
with every moment, its provenance, and an explicit **identification ledger**
recording which model parameters the moments pin down and which are left free.

Exercise:  EBA 2011 EU-wide stress test (published 15 Jul 2011).
Base:      31 Dec 2010 ACTUAL (aggregate rows carry a blank SCENARIO_CODE; the
           2011/2012 rows are baseline=100 / adverse=105 *projections*).
Sovereign: gross direct long exposures as-of 31 Dec 2010 (worksheet 5; the
           INFO_DATE_CODE=20111231 field is a reference-year label per EBA method).
Country map: D = Greece (banks GR030-GR035, issuer GR),
             F = Germany (banks DE017-DE029, issuer DE).

What changed vs. the 2026-07-22 build (see docs/eba_calibration.md "Rebuild"):

  * The **maturity ladder** (MATURITY_CODE 125..155) is now read. Previously only
    the 999 "total" row was used, so `delta_b` had no EBA counterpart at all and
    was carried over from an unrelated sovereign-average-maturity target.
  * `theta` is measured as CT1 / **GK-eligible** assets (corporate + commercial
    real estate + sovereign), not CT1/total assets. The total-assets version
    (14.9 / 32.9) was correctly rejected before: `theta` multiplies only the GK
    book, while total assets include interbank, reserves and retail.
  * `omega_K` is **measured** from the corporate+CRE EAD book against a K/Y
    target, instead of being back-solved as the residual that made a *assumed*
    theta=4.0 consistent with a capital-stock target. It is no longer a free plug.
  * Amplification is disciplined by an **exposure-ladder x observed-yield-move
    mark-to-market loss** (Acharya-Steffen "greatest carry trade" construction),
    NOT by the stress test's adverse-scenario CT1 depletion. The 2011 adverse
    scenario deliberately excluded sovereign default in the banking book, so its
    capital depletion understates sovereign pass-through by construction.

Run:  /opt/anaconda3/envs/ssj/bin/python code/eba_calibration.py
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from scipy.optimize import brentq

# ── Paths ────────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
DATA_DIR = os.path.join(_ROOT, "data")
CSV_PATH = os.path.join(DATA_DIR, "DATA_DISCLOSURE.CSV (1).csv")
GDP_PATH = os.path.join(DATA_DIR, "eurostat_gdp_raw.json")
YLD_PATH = os.path.join(_ROOT, "Empirics", "outputs", "spreads_fred.csv")
OUT_PATH = os.path.join(DATA_DIR, "eba_moments.json")

# ── Country / bank map ───────────────────────────────────────────────────────
GR_BANKS = ["GR030", "GR031", "GR032", "GR033", "GR034", "GR035"]
DATE_2010 = 20101231.0

WS_AGG = "1 - Aggregate information"
WS_EAD = "4 - EADs"
WS_SOV = "5 - Sovereign exposures"

# Worksheet 1 (MEASURE 120)
CODE_RWA, CODE_CT1, CODE_TA = 30010, 30014, 30029
# Worksheet 4 - credit-risk EAD by counterparty country
CODE_EAD_INST, CODE_EAD_CORP, CODE_EAD_CRE, CODE_EAD_TOT = 33010, 33011, 33018, 33021
# Worksheet 5 - sovereign exposure columns
CODE_SOV_GROSS_LONG, CODE_SOV_LOANS, CODE_SOV_NET = 34010, 34011, 34012
CODE_SOV_AFS, CODE_SOV_FVO, CODE_SOV_TRADING = 34013, 34014, 34015
MATURITY_TOTAL = 999.0

# Residual-maturity buckets. The EBA 2011 sovereign template labels the seven
# columns 3M / 1Y / 2Y / 3Y / 5Y / 10Y / 15Y; these are the bucket UPPER bounds.
# The seven buckets sum exactly to the MATURITY_TOTAL row (asserted in validate),
# so the last one is effectively "10Y and over" with a 15Y nominal label.
BUCKETS: dict[int, tuple[float, float]] = {
    125: (0.0, 0.25), 130: (0.25, 1.0), 135: (1.0, 2.0), 140: (2.0, 3.0),
    145: (3.0, 5.0),  150: (5.0, 10.0), 155: (10.0, 15.0),
}

# Representative coupons on the legacy stock outstanding at end-2010. GGB: the
# PSI-eligible stock averaged ~4.7%; Bunds ~3.5%. Used only for the repricing
# geometry; `coupon_sensitivity` in the output reports how much this matters.
COUPON_GGB, COUPON_BUND = 0.047, 0.035

# Model-side constants needed to translate a duration into `delta_b`.
BETA_INTER = 0.9975155088          # code/calibration.py beta_inter_D/F
K_OVER_Y_ANNUAL = 2.7              # model K=10.8 at quarterly Y=1


# ── Loaders ──────────────────────────────────────────────────────────────────
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


def sovereign_item(df, issuer, banks, code=CODE_SOV_GROSS_LONG,
                   maturity=MATURITY_TOTAL):
    """One sovereign-exposure cell (as-of 31 Dec 2010), EUR m."""
    w = df[(df["WORKSHEET"] == WS_SOV)
           & (df["C_COUNTRY_CODE"] == issuer)
           & (df["INFORMATION_CODE"] == code)
           & (df["MATURITY_CODE"] == maturity)]
    return float(_to_num(w[banks].iloc[0]).sum()) if len(w) else float("nan")


def sovereign_book(df, issuer, banks):
    """Gross direct long sovereign exposure, total across maturities, EUR m."""
    return sovereign_item(df, issuer, banks)


def sovereign_ladder(df, issuer, banks) -> dict[int, float]:
    """Gross direct long exposure by residual-maturity bucket, EUR m."""
    return {m: sovereign_item(df, issuer, banks, maturity=float(m))
            for m in BUCKETS}


def ead_item(df, code, country, banks):
    """One credit-risk EAD cell (31 Dec 2010), EUR m. country='TO' is the total."""
    w = df[(df["WORKSHEET"] == WS_EAD)
           & (df["INFORMATION_CODE"] == code)
           & (df["C_COUNTRY_CODE"] == country)]
    return float(_to_num(w[banks].iloc[0]).sum()) if len(w) else float("nan")


def load_gdp(gdp_path: str = GDP_PATH) -> dict:
    """Parse Eurostat JSON-stat -> {geo: {year: annual_nominal_GDP_EURm}}."""
    with open(gdp_path) as fh:
        js = json.load(fh)
    geo_idx = js["dimension"]["geo"]["category"]["index"]
    time_idx = js["dimension"]["time"]["category"]["index"]
    n_time = len(time_idx)
    values = js["value"]
    out: dict = {}
    for geo, gi in geo_idx.items():
        out[geo] = {}
        for yr, ti in time_idx.items():
            v = values.get(str(gi * n_time + ti))
            if v is not None:
                out[geo][yr] = float(v)
    return out


def load_yields(path: str = YLD_PATH) -> pd.DataFrame:
    """FRED long-term (10y) government bond yields, monthly, percent."""
    return pd.read_csv(path, parse_dates=["date"]).set_index("date")


# ── Bond maths ───────────────────────────────────────────────────────────────
def _bond_cashflows(t: float, coupon: float):
    """(times, cashflows) per 1 of face for an annual-coupon bond of maturity t."""
    if t < 1.0:                                   # sub-annual bucket: single payment
        return np.array([t]), np.array([1.0 + coupon * t])
    n = int(round(t))
    ts = np.arange(1, n + 1, dtype=float)
    cf = np.full(n, coupon)
    cf[-1] += 1.0
    return ts, cf


def bond_price(t: float, y: float, coupon: float) -> float:
    ts, cf = _bond_cashflows(t, coupon)
    return float((cf / (1.0 + y) ** ts).sum())


def ladder_value_and_duration(lad: dict[int, float], y: float, coupon: float):
    """Market value (EUR m) and MV-weighted Macaulay duration (years) of a ladder."""
    mv = 0.0
    num = 0.0
    for m, (lo, hi) in BUCKETS.items():
        t = 0.5 * (lo + hi)
        ts, cf = _bond_cashflows(t, coupon)
        disc = cf / (1.0 + y) ** ts
        p = float(disc.sum())
        mv += lad[m] * p
        num += lad[m] * float((ts * disc).sum())
    return mv, num / mv


def ladder_reprice(lad: dict[int, float], y0: float, y1: float, coupon: float):
    """(loss_fraction_of_initial_MV, initial_MV, terminal_MV) under a parallel shift."""
    mv0, _ = ladder_value_and_duration(lad, y0, coupon)
    mv1, _ = ladder_value_and_duration(lad, y1, coupon)
    return (mv0 - mv1) / mv0, mv0, mv1


def ladder_modified_duration(lad, y, coupon, bp=1e-4) -> float:
    """Effective modified duration (years): -dP/P per unit yield, by repricing."""
    loss, _, _ = ladder_reprice(lad, y, y + bp, coupon)
    return loss / bp


# ── Model-side translation: modified duration <-> delta_b ────────────────────
def model_modified_duration(delta_b: float, r_q: float = 1.0 / BETA_INTER - 1.0,
                            annual: bool = True) -> float:
    """Modified duration of the Hatchondo-Martinez perpetuity used in the model.

    The bond pays ``delta_b*(1-delta_b)**j`` at t+1+j, so
    ``q_b = delta_b / (r_q + delta_b)`` at a per-quarter discount rate ``r_q``,
    and ``-dlog q_b / d r_q = 1/(r_q + delta_b)`` quarters. Divide by 4 for years
    (a 100bp *annual* yield move is 25bp per quarter).
    """
    d_quarters = 1.0 / (r_q + delta_b)
    return d_quarters / 4.0 if annual else d_quarters


def delta_b_for_duration(mod_dur_years: float,
                         r_q: float = 1.0 / BETA_INTER - 1.0) -> float:
    """Invert `model_modified_duration`: the delta_b reproducing a target duration."""
    return brentq(lambda d: model_modified_duration(d, r_q) - mod_dur_years,
                  1e-8, 0.999)


# ── Moment construction ──────────────────────────────────────────────────────
def compute_moments(df=None, gdp=None, yld=None) -> dict:
    """Compute the full moment set with provenance. Returns a JSON-able dict."""
    if df is None:
        df = load_eba()
    if gdp is None:
        gdp = load_gdp()
    if yld is None:
        yld = load_yields()
    de_banks = _de_banks(df)

    # ---- bank aggregates (2010 actual base) --------------------------------
    ct1_D = group_sum(df, WS_AGG, CODE_CT1, DATE_2010, GR_BANKS)
    ct1_F = group_sum(df, WS_AGG, CODE_CT1, DATE_2010, de_banks)
    rwa_D = group_sum(df, WS_AGG, CODE_RWA, DATE_2010, GR_BANKS)
    rwa_F = group_sum(df, WS_AGG, CODE_RWA, DATE_2010, de_banks)
    ta_D = group_sum(df, WS_AGG, CODE_TA, DATE_2010, GR_BANKS)
    ta_F = group_sum(df, WS_AGG, CODE_TA, DATE_2010, de_banks)

    # ---- sovereign cross-holding matrix ------------------------------------
    b_D_D = sovereign_book(df, "GR", GR_BANKS)   # GR banks' Greek book
    b_D_F = sovereign_book(df, "GR", de_banks)   # DE banks' Greek book (contagion)
    b_F_D = sovereign_book(df, "DE", GR_BANKS)   # GR banks' Bund book
    b_F_F = sovereign_book(df, "DE", de_banks)   # DE banks' Bund book

    # ---- credit-risk EAD: the bank claim on productive capital -------------
    # Corporate (ex-CRE) + Commercial Real Estate = claims on *business* capital.
    # Residential mortgages are excluded: the model has no housing stock.
    # Own-country is the right concept for a two-country model in which each
    # bank funds its own country's capital; the total book is reported as a
    # sensitivity (it embeds DE banks' large non-euro-area corporate lending).
    corp_D, cre_D = ead_item(df, CODE_EAD_CORP, "GR", GR_BANKS), ead_item(df, CODE_EAD_CRE, "GR", GR_BANKS)
    corp_F, cre_F = ead_item(df, CODE_EAD_CORP, "DE", de_banks), ead_item(df, CODE_EAD_CRE, "DE", de_banks)
    corp_D_tot, cre_D_tot = ead_item(df, CODE_EAD_CORP, "TO", GR_BANKS), ead_item(df, CODE_EAD_CRE, "TO", GR_BANKS)
    corp_F_tot, cre_F_tot = ead_item(df, CODE_EAD_CORP, "TO", de_banks), ead_item(df, CODE_EAD_CRE, "TO", de_banks)
    kbank_D, kbank_F = corp_D + cre_D, corp_F + cre_F

    # ---- GDP ----------------------------------------------------------------
    gdp_ann_D, gdp_ann_F = gdp["EL"]["2010"], gdp["DE"]["2010"]
    qgdp_D, qgdp_F = gdp_ann_D / 4.0, gdp_ann_F / 4.0

    # ---- theta: GK leverage on the GK-eligible book ------------------------
    gk_assets_D = kbank_D + b_D_D + b_F_D
    gk_assets_F = kbank_F + b_F_F + b_D_F
    theta_D, theta_F = gk_assets_D / ct1_D, gk_assets_F / ct1_F

    # ---- omega_K: measured bank share of the capital stock -----------------
    K_D, K_F = K_OVER_Y_ANNUAL * gdp_ann_D, K_OVER_Y_ANNUAL * gdp_ann_F
    omega_K_D, omega_K_F = kbank_D / K_D, kbank_F / K_F

    # ---- yields and the maturity ladder ------------------------------------
    y_GR_0 = float(yld.loc["2010-12-01", "GRC"]) / 100.0
    y_GR_1 = float(yld.loc["2011-12-01", "GRC"]) / 100.0
    y_DE_0 = float(yld.loc["2010-12-01", "DEU"]) / 100.0
    y_DE_1 = float(yld.loc["2011-12-01", "DEU"]) / 100.0

    ladders = {
        "b_D_D": (sovereign_ladder(df, "GR", GR_BANKS), y_GR_0, y_GR_1, COUPON_GGB),
        "b_D_F": (sovereign_ladder(df, "GR", de_banks), y_GR_0, y_GR_1, COUPON_GGB),
        "b_F_F": (sovereign_ladder(df, "DE", de_banks), y_DE_0, y_DE_1, COUPON_BUND),
        "b_F_D": (sovereign_ladder(df, "DE", GR_BANKS), y_DE_0, y_DE_1, COUPON_BUND),
    }
    ladder_stats = {}
    for key, (lad, y0, y1, cpn) in ladders.items():
        face = sum(lad.values())
        mv0, mac = ladder_value_and_duration(lad, y0, cpn)
        mod = ladder_modified_duration(lad, y0, cpn)
        loss11, _, _ = ladder_reprice(lad, y0, y1, cpn)
        ladder_stats[key] = {
            "face_EURm": face,
            "shares_by_bucket": {str(m): lad[m] / face for m in sorted(BUCKETS)},
            "wavg_residual_maturity_y": sum(
                lad[m] * 0.5 * (lo + hi) for m, (lo, hi) in BUCKETS.items()) / face,
            "yield_2010_12": y0, "yield_2011_12": y1,
            "market_value_2010_EURm": mv0,
            "macaulay_duration_y": mac,
            "modified_duration_y": mod,
            "mtm_loss_2011_frac_of_MV": loss11,
            "mtm_loss_2011_EURm": loss11 * mv0,
            "delta_b_implied": delta_b_for_duration(mod),
        }

    # ---- the amplification moment (Acharya-Steffen construction) -----------
    # Mechanical mark-to-market hit to Core Tier 1 per 100bp on the OWN sovereign,
    # and the realised 2011 hit. Full book (not just AFS/trading): the model's
    # n_inter is economic net worth and marks the whole book through q_b.
    def mtm_block(lad_key, ct1):
        s = ladder_stats[lad_key]
        return {
            "dNW_per_100bp": -s["modified_duration_y"] * 0.01
                             * s["market_value_2010_EURm"] / ct1,
            "mtm_loss_2011_over_CT1": -s["mtm_loss_2011_EURm"] / ct1,
        }

    mtm_D, mtm_F = mtm_block("b_D_D", ct1_D), mtm_block("b_F_F", ct1_F)

    # Accounting vs economic: what share of the book was actually fair-valued.
    fv_share_D = ((sovereign_item(df, "GR", GR_BANKS, CODE_SOV_AFS)
                   + sovereign_item(df, "GR", GR_BANKS, CODE_SOV_FVO)
                   + sovereign_item(df, "GR", GR_BANKS, CODE_SOV_TRADING)) / b_D_D)
    fv_share_F = ((sovereign_item(df, "DE", de_banks, CODE_SOV_AFS)
                   + sovereign_item(df, "DE", de_banks, CODE_SOV_FVO)
                   + sovereign_item(df, "DE", de_banks, CODE_SOV_TRADING)) / b_F_F)

    # ---- coupon sensitivity (the one free convention in the bond maths) ----
    coupon_sens = {}
    for cpn in (0.030, COUPON_GGB, 0.060):
        lad, y0, _, _ = ladders["b_D_D"]
        coupon_sens[f"{cpn:.3f}"] = ladder_modified_duration(lad, y0, cpn)

    # ---- BROAD-SECTOR scope (2026-07-31) -----------------------------------
    # The CT1 scope below takes n_inter = Core Tier 1 of the EBA stress-test
    # sample. That is the wrong object for this model: n_inter is the net worth
    # of the agent intermediating the WHOLE capital stock, whereas CT1 is the
    # capital of the sovereign-exposed sub-sample. Using CT1 forces a tiny
    # omega_K (banks fund ~12% of K), and holding K/Y at target then makes the
    # accelerator gain dK/dN = theta*K/(N*(theta-phi)) blow up ~ 1/n_inter.
    #
    # Broad scope instead defines the intermediary as the entire capital-funding
    # sector and lets its net worth follow from the measured leverage and the
    # balance sheet:
    #       N_broad = (Q*K + sovereign_book) / theta
    # with K = K/Y_annual * annual GDP and theta the MEASURED GK leverage.
    # Then omega_K = 1 by construction (the sector funds all capital, so the
    # passive-fund device disappears) and phi_own = sovereign / N_broad.
    #
    # What is kept measured: theta (a ratio), the sovereign book (a level), K/Y.
    # What is given up: n_inter as directly-observed CT1, and phi_own = 2.39 as a
    # MODEL parameter. 2.39 remains a true statement about the stress-tested
    # slice; it is concentration *within* that slice, not within the whole
    # capital-funding sector, and only the latter is what this model's phi_own
    # means.
    # KEY ASSUMPTION, flagged in the ledger: applying the EBA sample's theta to
    # the whole sector assumes non-stress-tested capital funding levers the same
    # way. If its true leverage is lower, N_broad is larger and phi_own smaller.
    K_lvl_D, K_lvl_F = K_D, K_F                       # EURm, = K/Y_ann * GDP_ann
    sov_book_D, sov_book_F = b_D_D + b_F_D, b_F_F + b_D_F
    N_broad_D = (K_lvl_D + sov_book_D) / theta_D
    N_broad_F = (K_lvl_F + sov_book_F) / theta_F
    targets_broad = {
        "n_inter_D": N_broad_D / qgdp_D,
        "n_inter_F": N_broad_F / qgdp_F,
        "phi_bD_D_ss": b_D_D / N_broad_D,
        "phi_bF_F_ss": b_F_F / N_broad_F,
        "phi_bD_F_ss": b_D_F / N_broad_F,
        "phi_bF_D_ss": b_F_D / N_broad_D,
        "theta_D": theta_D, "theta_F": theta_F,
        "omega_K_D": 1.0, "omega_K_F": 1.0,       # sector funds all capital
        "delta_b_D": ladder_stats["b_D_D"]["delta_b_implied"],
        "delta_b_F": ladder_stats["b_F_F"]["delta_b_implied"],
        "B_supply_D_qgdp": (b_D_D + b_D_F) / qgdp_D,
        "B_supply_F_qgdp": (b_F_F + b_F_D) / qgdp_F,
    }

    targets = {
        "n_inter_D": ct1_D / qgdp_D,
        "n_inter_F": ct1_F / qgdp_F,
        "phi_bD_D_ss": b_D_D / ct1_D,
        "phi_bF_F_ss": b_F_F / ct1_F,
        "phi_bD_F_ss": b_D_F / ct1_F,
        "phi_bF_D_ss": b_F_D / ct1_D,
        "theta_D": theta_D,
        "theta_F": theta_F,
        "omega_K_D": omega_K_D,
        "omega_K_F": omega_K_F,
        "delta_b_D": ladder_stats["b_D_D"]["delta_b_implied"],
        "delta_b_F": ladder_stats["b_F_F"]["delta_b_implied"],
        "B_supply_D_qgdp": (b_D_D + b_D_F) / qgdp_D,
        "B_supply_F_qgdp": (b_F_F + b_F_D) / qgdp_F,
    }

    return {
        "meta": {
            "source_exercise": "EBA 2011 EU-wide stress test",
            "rebuilt": "2026-07-31",
            "base_date": "2010-12-31 (actual; scenario blank)",
            "sovereign_asof": "2010-12-31 gross direct long (34010), 7 maturity buckets",
            "gdp_source": "Eurostat nama_10_gdp B1GQ CP_MEUR 2010",
            "yield_source": "FRED long-term (10y) govt bond yields via Empirics/outputs/spreads_fred.csv",
            "country_map": {"D": "Greece (GR030-035, issuer GR)",
                            "F": "Germany (DE0xx, issuer DE)"},
            "de_banks": de_banks,
            "units": "EUR million unless noted",
            "coupon_assumption": {"GGB": COUPON_GGB, "Bund": COUPON_BUND},
            "K_over_Y_annual": K_OVER_Y_ANNUAL,
            "beta_inter": BETA_INTER,
            "NOT_used": ("adverse-scenario (105) CT1 depletion -- the 2011 exercise "
                         "excluded sovereign default in the banking book, so its "
                         "capital depletion understates sovereign pass-through by "
                         "construction and must not identify amplification"),
        },
        "raw_EURm": {
            "CT1_D": ct1_D, "CT1_F": ct1_F,
            "RWA_D": rwa_D, "RWA_F": rwa_F,
            "TA_D": ta_D, "TA_F": ta_F,
            "b_D_D": b_D_D, "b_D_F": b_D_F, "b_F_D": b_F_D, "b_F_F": b_F_F,
            "corp_own_D": corp_D, "cre_own_D": cre_D,
            "corp_own_F": corp_F, "cre_own_F": cre_F,
            "corp_tot_D": corp_D_tot, "cre_tot_D": cre_D_tot,
            "corp_tot_F": corp_F_tot, "cre_tot_F": cre_F_tot,
            "K_bank_D": kbank_D, "K_bank_F": kbank_F,
            "GK_assets_D": gk_assets_D, "GK_assets_F": gk_assets_F,
            "GDP_ann_D": gdp_ann_D, "GDP_ann_F": gdp_ann_F,
        },
        "ladder": ladder_stats,
        "mtm": {
            "D": mtm_D, "F": mtm_F,
            "fair_valued_share_of_own_book": {"D": fv_share_D, "F": fv_share_F},
            "note": ("dNW_per_100bp is the MECHANICAL revaluation of the whole "
                     "own-sovereign book against Core Tier 1. Any excess of the "
                     "model's general-equilibrium net-worth response over this "
                     "number is amplification and is what psi_lambda_B governs."),
        },
        "leverage_alternatives": {
            "theta_D_total_assets": ta_D / ct1_D,
            "theta_F_total_assets": ta_F / ct1_F,
            "theta_D_total_book": (corp_D_tot + cre_D_tot + b_D_D + b_F_D) / ct1_D,
            "theta_F_total_book": (corp_F_tot + cre_F_tot + b_F_F + b_D_F) / ct1_F,
            "omega_K_D_total_book": (corp_D_tot + cre_D_tot) / K_D,
            "omega_K_F_total_book": (corp_F_tot + cre_F_tot) / K_F,
        },
        "coupon_sensitivity_modified_duration_b_D_D": coupon_sens,
        "model_targets": targets,            # CT1 scope (historical)
        "model_targets_broad": targets_broad,  # broad capital-funding sector
        "identification": _identification_ledger(),
    }


def _identification_ledger() -> dict:
    """Which parameters this moment set pins down, and which remain free."""
    return {
        "identified_by_EBA_2011": {
            "n_inter_D/F": "CT1 / own quarterly nominal GDP",
            "phi_bD_D_ss, phi_bF_F_ss": "own-sovereign gross long / CT1",
            "phi_bD_F_ss, phi_bF_D_ss": "cross-border sovereign gross long / CT1",
            "theta_D/F": "(corporate ex-CRE + CRE + sovereign) EAD / CT1",
            "delta_b_D/F": ("modified duration of the own-sovereign maturity "
                            "ladder at the 31-Dec-2010 market yield, inverted "
                            "through the HM perpetuity"),
            "B_supply_D/F": "bank-held sovereign stock / own quarterly GDP",
        },
        "identified_jointly_with_a_standard_macro_target": {
            "omega_K_D/F": ("measured corporate+CRE EAD divided by K, where K "
                            "comes from the conventional K/Y_annual = 2.7. Not a "
                            "residual: the EAD book is observed, and the implied "
                            "theta - phi_own - phi_cross balance-sheet identity "
                            "reproduces it (asserted in validate())."),
        },
        "still_NOT_identified": {
            "Delta_bD_D, Delta_bF_D, Delta_bF_F, Delta_bD_F": (
                "divertability weights in the multi-asset IC. Hardcoded 0.2/0.4. "
                "No EBA counterpart and no moment attached. These set the "
                "collateral value of sovereign debt and are the single largest "
                "unidentified input to the doom loop."),
            "psi_lambda_B_D/F": (
                "state-dependent divertability / amplification. The EBA ladder "
                "now pins the MECHANICAL MTM channel, so psi_lambda_B no longer "
                "absorbs a mis-calibrated mechanical loss -- but its own level "
                "still has to be matched to something. Matching it to the 150bp "
                "GR-DE spread is a single-moment tune; identifying it properly "
                "needs bank equity returns regressed on the EBA exposure "
                "cross-section (Acharya-Steffen), which is NOT in this repo."),
            "def_scale_D/F": "shock scaling; hand-set at 0.25, no moment.",
            "f_D/F": "banker exit rate; 0.12, from the GK literature, not EBA.",
            "theta ratio-vs-level": (
                "theta is measured on the OWN-country book. The total-book "
                "alternative (7.21 / 13.30) is equally defensible for F and is "
                "reported under leverage_alternatives."),
            "phi_lamb_D/F": "fiscal feedback; ~Bohn, not an EBA object.",
            "recovery_rate_D/F": "0.30 from the PSI NPV literature (EL-1), not EBA.",
        },
        "deliberately_rejected": {
            "adverse_scenario_CT1_depletion": (
                "EBA 2011's adverse scenario excluded sovereign default in the "
                "banking book. Its capital depletion is not a measure of "
                "sovereign-stress pass-through and must not identify "
                "amplification."),
            "theta_from_total_assets": (
                "CT1/total assets (14.9 D, 32.9 F) includes interbank, reserves "
                "and retail, which theta does not multiply. Previously verified "
                "not to converge in steady state."),
        },
    }


# ── Validation ───────────────────────────────────────────────────────────────
def validate(df=None, moments=None) -> None:
    """Self-checks on the decode and on the internal consistency of the moments."""
    if df is None:
        df = load_eba()
    # 1. Code map: DE017 must reproduce Deutsche Bank's published 2010 figures.
    rwa = group_sum(df, WS_AGG, CODE_RWA, DATE_2010, ["DE017"])
    ct1 = group_sum(df, WS_AGG, CODE_CT1, DATE_2010, ["DE017"])
    assert abs(rwa - 346_608) < 5, f"DE017 RWA {rwa:,.0f} != Deutsche Bank 346,608"
    assert abs(ct1 / rwa - 0.0876) < 0.0010, f"DE017 CT1 ratio {ct1/rwa:.4f} != 0.0876"

    # 2. The seven maturity buckets must exhaust the reported total.
    de_banks = _de_banks(df)
    for issuer, banks, lab in [("GR", GR_BANKS, "GR/GR"), ("DE", de_banks, "DE/DE"),
                               ("GR", de_banks, "DE/GR"), ("DE", GR_BANKS, "GR/DE")]:
        lad = sovereign_ladder(df, issuer, banks)
        tot = sovereign_book(df, issuer, banks)
        assert abs(sum(lad.values()) - tot) < max(1.0, 1e-4 * tot), \
            f"{lab}: ladder {sum(lad.values()):,.0f} != total {tot:,.0f}"

    # 3. delta_b round-trip: duration -> delta_b -> duration.
    for d in (0.02, 0.07, 0.15):
        assert abs(delta_b_for_duration(model_modified_duration(d)) - d) < 1e-9

    if moments is None:
        return

    # 4. Balance-sheet identity: theta - phi_own - phi_cross must reproduce the
    #    measured bank capital book, i.e. omega_K is not a free residual.
    t, raw = moments["model_targets"], moments["raw_EURm"]
    for c, ct1k, kbk in [("D", "CT1_D", "K_bank_D"), ("F", "CT1_F", "K_bank_F")]:
        own = t[f"phi_bD_D_ss"] if c == "D" else t[f"phi_bF_F_ss"]
        cross = t[f"phi_bF_D_ss"] if c == "D" else t[f"phi_bD_F_ss"]
        implied = (t[f"theta_{c}"] - own - cross) * raw[ct1k]
        assert abs(implied - raw[kbk]) < max(1.0, 1e-6 * raw[kbk]), \
            f"{c}: balance sheet implies K_bank {implied:,.0f} != measured {raw[kbk]:,.0f}"

    # 5. The implied delta_b must reproduce the ladder's modified duration.
    for key, cty in [("b_D_D", "D"), ("b_F_F", "F")]:
        s = moments["ladder"][key]
        assert abs(model_modified_duration(t[f"delta_b_{cty}"])
                   - s["modified_duration_y"]) < 1e-8


def main() -> dict:
    df = load_eba()
    moments = compute_moments(df)
    validate(df, moments)
    with open(OUT_PATH, "w") as fh:
        json.dump(moments, fh, indent=2)

    t = moments["model_targets"]
    print(f"wrote {OUT_PATH}")
    print("\n── EBA 2011 model targets (rebuilt) ──")
    for k, v in t.items():
        print(f"  {k:22} {v:12.4f}")

    print("\n── maturity ladder ──")
    for k, s in moments["ladder"].items():
        print(f"  {k}: face {s['face_EURm']:10,.0f}  wavg-mat {s['wavg_residual_maturity_y']:5.2f}y "
              f" ModDur@{s['yield_2010_12']:.2%} {s['modified_duration_y']:5.2f}y "
              f" -> delta_b {s['delta_b_implied']:.4f}")

    print("\n── mechanical MTM (Acharya-Steffen construction) ──")
    for c in ("D", "F"):
        m = moments["mtm"][c]
        print(f"  {c}: dNW/CT1 per 100bp = {m['dNW_per_100bp']:7.2%}   "
              f"realised 2011 = {m['mtm_loss_2011_over_CT1']:7.2%} of CT1")
    fv = moments["mtm"]["fair_valued_share_of_own_book"]
    print(f"  fair-valued share of own book: D {fv['D']:.1%}, F {fv['F']:.1%} "
          f"(rest at amortised cost -- economic, not accounting, loss)")

    print("\n── BROAD-SECTOR scope (the live one) ──")
    for k, v in moments["model_targets_broad"].items():
        print(f"  {k:22} {v:12.4f}")

    print("\n── still unidentified ──")
    for k in moments["identification"]["still_NOT_identified"]:
        print(f"  * {k}")
    return moments


if __name__ == "__main__":
    main()
