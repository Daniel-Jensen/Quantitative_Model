"""Greek counterparts to Bocola (2016)'s intermediary calibration moments.

Purpose (paper role)
--------------------
Bocola (2016, JPE 124(4)) disciplines the Gertler-Karadi block of his Italian
model with two data objects and one estimate. This module (a) reproduces his
construction exactly from his own replication archive, and (b) rebuilds the same
objects for Greece from official sources, so the paper's intermediary block is
not disciplined entirely from Italian evidence.

It is an EMPIRICAL AUDIT. It changes no model parameter and is imported by
nothing in `code/`.

What Bocola actually does (paper + archive; see docs/greek_bocola_calibration.md)
---------------------------------------------------------------------------------
  * `lev_bg = 5.000` is Table 1, "chosen so that the model produces a leverage
    ratio of 5 in a balanced growth path, in the range of values used in the
    literature". It is NOT an Italian moment: his own measured Italian flow-of-
    funds leverage averages 8.3.
  * The Euribor-Eurepo spread averages 8bp ANNUALISED over 2002Q2-2012Q4. It
    enters only through the constructed multiplier series, as an observable.
  * `mu_bg` is ESTIMATED (Table 2), posterior mean 0.001, prior centred on
    Gertler-Karadi (2011).

So (lev, spread) -> mu is a MEASUREMENT chain, not a calibration chain. The only
reason the three numbers cohere is that eq. (B2) at lev=5 and 8bp/yr returns
0.001, which is where the likelihood also lands.

Run:
    /opt/anaconda3/envs/ssj/bin/python Empirics/bocola_greece/moments.py

Outputs (Empirics/outputs/):
    bocola_greece_series.csv    cached quarterly series (offline rebuild)
    bocola_greece_moments.json  the moment table, for the memo
"""
from __future__ import annotations

import io
import json
import os
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "Empirics", "outputs")
XLS = os.path.join(ROOT, "data", "bocola2016", "Data.xls")
SERIES_CSV = os.path.join(OUT, "bocola_greece_series.csv")

# ── Bocola's Table 1 / Table 2, transcribed with provenance ────────────────────
# Values are quoted from the published paper. `source` records how each is
# obtained, in his words where they matter.
BOCOLA_TABLE1 = {
    "exp_bg":  (0.076, "EBA 2011 net direct holdings of Italian government debt, "
                       "5 largest Italian banks, over Bankscope consolidated total "
                       "assets at end-2010. Archive: 160091/2098237.6 = 0.076298."),
    "pi":      (0.056, "Fraction of bonds maturing per quarter. Archive computes a "
                       "17.90-quarter average maturity from the EBA maturity ladder; "
                       "1/0.056 = 17.9 quarters. The paper's text says '18 months', "
                       "which is inconsistent with its own Table 1 and archive - read "
                       "as 18 QUARTERS."),
    "g_ratio": (0.198, "Average government-spending-output ratio; AR(1) on the "
                       "spending-output ratio, 1999Q1-2009Q4."),
    "rho_g":   (0.894, "Same AR(1)."),
    "sigma_g": (0.012, "Same AR(1)."),
    "i_over_y": (0.213, "Sample average of the investment-output ratio."),
    "hours":   (0.318, "Sample average of hours worked."),
    "R_bg":    (1.003, "Sample average of the risk-free rate, gross quarterly "
                       "(0.3%/qtr = 1.2%/yr)."),
    "alpha":   (0.300, "Sample average of the labour income share."),
    "nu_inv":  (0.500, "Inverse Frisch elasticity; chosen for a Frisch elasticity of 2."),
    "lev_bg":  (5.000, "NOT a data moment. 'The parameter q is chosen so that the "
                       "model produces a leverage ratio of 5 in a balanced growth "
                       "path, in the range of values used in the literature.'"),
    "gamma_t": (1.000, "Fiscal rule; lump-sum taxes so the choice is near-inert."),
    "adj_bg":  (0.000, "Capital adjustment costs zero along the BGP, by convention."),
    "q_b_bg":  (1.000, "Bond price normalised to one."),
}
BOCOLA_MU_POSTERIOR = {
    "mean": 0.0010, "ci90": (0.0002, 0.0021),
    "source": "Table 2, mu_bg x 100: prior Normal(0.96, 25) - i.e. centred on the "
              "Gertler-Karadi (2011) calibrated value - posterior mean 0.10, "
              "90% credible set [0.02, 0.21]. Observables: GDP growth and the "
              "constructed multiplier series, 2002Q2-2009Q4 (31 quarters).",
}
# Windows. Bocola's own two are the only ones with a documented rationale.
WINDOWS = {
    "bocola_estimation_2002Q2_2009Q4": ("2002-Q2", "2009-Q4"),
    "bocola_figure1_2002Q2_2012Q4":    ("2002-Q2", "2012-Q4"),
    "greek_precrisis_1999Q1_2007Q2":   ("1999-Q1", "2007-Q2"),
    "greek_precrisis_2002Q2_2007Q2":   ("2002-Q2", "2007-Q2"),
    "greek_crisis_2010Q1_2012Q4":      ("2010-Q1", "2012-Q4"),
}


# ── eq. (16)/(B2) ─────────────────────────────────────────────────────────────
def mu_of(spread_q, lev, R_gross=1.003, divide_by_R=True):
    """Bocola eq. (B2). `spread_q` is R_interbank - R in GROSS QUARTERLY units.

    The published equation divides the spread by R_t; his spreadsheet does not.
    The two differ by a factor 1/R ~ 0.3%, immaterial at these magnitudes, but
    `divide_by_R=False` reproduces the archive bit-for-bit.
    """
    x = (spread_q / R_gross if divide_by_R else spread_q) * lev
    return x / (1.0 + x)


def wedge_of(mu, lev, R_gross=1.003):
    """Invert (B2): the per-quarter excess return implied by a multiplier."""
    return (mu / (1.0 - mu)) * R_gross / lev


# ── Bocola's Italian data, from his own archive ───────────────────────────────
def _qidx(dates):
    return pd.PeriodIndex(pd.to_datetime(dates, format="%Y/%m/%d", errors="coerce"),
                          freq="Q").astype(str).str.replace("Q", "-Q")


def bocola_italy():
    """Leverage, Euribor-Eurepo spread and the multiplier, from Data.xls."""
    x = pd.ExcelFile(XLS)

    nfa = x.parse("National Financial Accounts", header=None).iloc[11:]
    lev = pd.DataFrame({
        "TA": pd.to_numeric(nfa[2], errors="coerce"),
        "TL": pd.to_numeric(nfa[5], errors="coerce"),
        "EQ": pd.to_numeric(nfa[8], errors="coerce"),
        "MF": pd.to_numeric(nfa[11], errors="coerce"),
        "lev_xls": pd.to_numeric(nfa[14], errors="coerce"),
    }).set_index(_qidx(nfa[0])).dropna(subset=["TA"])
    # Appendix B.C: equity = total financial assets - total DEBT liabilities,
    # debt liabilities = total liabilities - shares and other equity - mutual fund shares.
    lev["equity"] = lev.TA - lev.TL + lev.EQ + lev.MF
    lev["lev"] = lev.TA / lev.equity

    ib = x.parse("Interbank Spread", header=None).iloc[11:]
    sp = pd.DataFrame({
        "euribor_1w_q": pd.to_numeric(ib[2], errors="coerce"),
        "eurepo_1w_q": pd.to_numeric(ib[4], errors="coerce"),
        "spread_xls": pd.to_numeric(ib[6], errors="coerce"),
    }).set_index(_qidx(ib[0]))
    sp["spread"] = sp.euribor_1w_q - sp.eurepo_1w_q

    lm = x.parse("Lagrange_Multiplier", header=None).iloc[11:]
    mu = pd.DataFrame({
        "mu_xls": pd.to_numeric(lm[6], errors="coerce"),
        "gdp_growth": pd.to_numeric(lm[8], errors="coerce"),
    }).set_index(_qidx(lm[0]))

    df = pd.concat([lev[["TA", "equity", "lev", "lev_xls"]],
                    sp[["euribor_1w_q", "eurepo_1w_q", "spread", "spread_xls"]],
                    mu], axis=1)
    df["mu"] = mu_of(df.spread, df.lev, divide_by_R=False)
    return df


# ── Greek and Italian leverage on Bocola's definition, from Eurostat ──────────
EUROSTAT = ("https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/"
            "data/nasq_10_f_bs")
# S122 alone is not published for EL or IT; S121_S122_S123 (central bank +
# deposit-takers + MMFs) is the closest available match to Bocola's "monetary
# and financial institutions". The definitional cost is documented in the memo.
EUROSTAT_SECTOR = "S121_S122_S123"


def _eurostat(geo, item, finpos, sector=EUROSTAT_SECTOR):
    qs = (f"?format=JSON&lang=EN&unit=MIO_EUR&geo={geo}&sector={sector}"
          f"&na_item={item}&finpos={finpos}")
    with urllib.request.urlopen(EUROSTAT + qs, timeout=180) as r:
        d = json.load(r)
    idx = d["dimension"]["time"]["category"]["index"]
    inv = {v: k for k, v in idx.items()}
    return pd.Series({inv[int(k)]: v for k, v in d["value"].items()}).sort_index()


def eurostat_leverage(geo):
    """lev = total financial assets / (assets - liabilities + F51_liab + F52_liab)."""
    A = _eurostat(geo, "F", "ASS")
    L = _eurostat(geo, "F", "LIAB")
    E = _eurostat(geo, "F51", "LIAB")
    M = _eurostat(geo, "F52", "LIAB").reindex(A.index).fillna(0.0)
    df = pd.concat({"A": A, "L": L, "E": E, "M": M}, axis=1).dropna(subset=["A", "L", "E"])
    df["equity"] = df.A - df.L + df.E + df.M
    df["lev"] = df.A / df.equity
    return df


# ── ECB series: the BSI leverage cross-check and the funding-wedge candidates ──
ECB = {
    # MFI balance sheet: total assets, and capital and reserves. Monthly.
    "bsi_ta_GR":  "BSI/M.GR.N.A.T00.A.1.Z5.0000.Z01.E",
    "bsi_cap_GR": "BSI/M.GR.N.A.L60.X.1.Z5.0000.Z01.E",
    "bsi_ta_IT":  "BSI/M.IT.N.A.T00.A.1.Z5.0000.Z01.E",
    "bsi_cap_IT": "BSI/M.IT.N.A.L60.X.1.Z5.0000.Z01.E",
    # Money market. EUREPO is an EBF index, discontinued in 2015 and not on the
    # ECB API; Bocola's own 1-week series (Data.xls) is the only copy we have, so
    # EONIA stands in only for the post-2012 extension and is flagged as such.
    "euribor_3m": "FM/M.U2.EUR.RT.MM.EURIBOR3MD_.HSTA",
    "eonia":      "FM/M.U2.EUR.4F.MM.EONIA.HSTA",
    # Bank funding premia: MIR rates on NEW household time deposits up to 1 year.
    # Liability-side, so an imperfect stand-in for Bocola's asset-side wedge.
    "dep_hh_GR":  "MIR/M.GR.B.L22.F.R.A.2250.EUR.N",
    "dep_hh_DE":  "MIR/M.DE.B.L22.F.R.A.2250.EUR.N",
    "dep_hh_IT":  "MIR/M.IT.B.L22.F.R.A.2250.EUR.N",
}


def _ecb(key, start="1999-01", end="2013-12"):
    u = (f"https://data-api.ecb.europa.eu/service/data/{key}"
         f"?format=csvdata&startPeriod={start}&endPeriod={end}")
    with urllib.request.urlopen(u, timeout=120) as r:
        d = pd.read_csv(io.BytesIO(r.read()))
    d = d[["TIME_PERIOD", "OBS_VALUE"]].dropna()
    s = pd.Series(d.OBS_VALUE.values, index=pd.PeriodIndex(d.TIME_PERIOD, freq="M"))
    return s.groupby(s.index.asfreq("Q")).mean()      # quarterly average


def build_series(refresh=False):
    """Everything on one quarterly index. Cached: the memo must rebuild offline."""
    if os.path.exists(SERIES_CSV) and not refresh:
        return pd.read_csv(SERIES_CSV, index_col=0)

    it = bocola_italy()
    frames = {
        "lev_bocola_IT": it.lev,
        "spread_1w_q": it.spread,
        "mu_bocola_IT": it.mu,
        "gdp_growth_IT": it.gdp_growth,
    }
    for geo, tag in (("IT", "IT"), ("EL", "GR")):
        frames[f"lev_eurostat_{tag}"] = eurostat_leverage(geo).lev
    ecb = {}
    for name, key in ECB.items():
        s = _ecb(key)
        s.index = s.index.astype(str).str.replace("Q", "-Q")
        ecb[name] = s
        print(f"  pulled {name}: {len(s)} quarters")
    for tag in ("GR", "IT"):
        frames[f"lev_bsi_{tag}"] = ecb[f"bsi_ta_{tag}"] / ecb[f"bsi_cap_{tag}"]
    # Funding premia, percent per annum -> gross quarterly decimal.
    for tag in ("GR", "IT", "DE"):
        frames[f"deprate_{tag}"] = ecb[f"dep_hh_{tag}"]
    frames["euribor_3m"] = ecb["euribor_3m"]
    frames["eonia"] = ecb["eonia"]

    df = pd.concat(frames, axis=1).sort_index()
    df = df.loc[[i for i in df.index if isinstance(i, str) and "-Q" in i]]
    # Candidate Greek wedges, all in GROSS QUARTERLY units to match `spread_1w_q`.
    pa_to_q = lambda pa: pa / 400.0
    df["wedge_depGR_over_eonia"] = pa_to_q(df.deprate_GR - df.eonia)
    df["wedge_depGR_over_depDE"] = pa_to_q(df.deprate_GR - df.deprate_DE)
    df["wedge_depIT_over_eonia"] = pa_to_q(df.deprate_IT - df.eonia)
    df.to_csv(SERIES_CSV)
    return df


# ── Bocola's exp_bg and pi, and their Greek counterparts, from the EBA panel ──
# Bocola: net direct position in the OWN sovereign for the 5 largest domestic
# banks, over consolidated total assets at end-2010 (Bankscope); average maturity
# from the same maturity ladder, bucket upper bounds in quarters.
BUCKET_Q = {125: 1, 130: 4, 135: 8, 140: 12, 145: 20, 150: 40, 155: 60}
BOCOLA_TOP5 = {"gr_net": 160091.0, "ta": 2098237.6, "avg_maturity_q": 17.902924}


def eba_exposure_moments():
    """Greek counterpart to exp_bg and pi, plus Bocola's own Italian numbers."""
    panel = os.path.join(OUT, "eba_bank_panel.csv")
    out = {"Italy, Bocola top-5 (archive)": {
        "net_over_ta": BOCOLA_TOP5["gr_net"] / BOCOLA_TOP5["ta"],
        "gross_over_ta": float("nan"),
        "avg_maturity_q": BOCOLA_TOP5["avg_maturity_q"], "n_banks": 5}}
    if not os.path.exists(panel):
        print(f"  (skipped Greece: run Empirics/eba_descriptives.py first)")
        return out
    p = pd.read_csv(panel, index_col=0)
    p["nat"] = [b[:2] for b in p.index]
    gr = p[p.nat == "GR"]
    lad = {m: float(gr[f"gr_{m}"].sum()) for m in BUCKET_Q}
    face = sum(lad.values())
    out["Greece, 6 EBA banks"] = {
        "net_over_ta": float(gr.gr_net.sum() / gr.ta.sum()),
        "gross_over_ta": float(gr.gr_gross.sum() / gr.ta.sum()),
        "avg_maturity_q": sum(w * BUCKET_Q[m] for m, w in lad.items()) / face,
        "n_banks": int(len(gr))}
    return out


# ── reporting ─────────────────────────────────────────────────────────────────
def window_stats(df, col, win):
    lo, hi = win
    s = pd.to_numeric(df.loc[(df.index >= lo) & (df.index <= hi), col],
                      errors="coerce").dropna()
    if not len(s):
        return None
    return {"n": int(len(s)), "mean": float(s.mean()), "median": float(s.median()),
            "min": float(s.min()), "max": float(s.max())}


def main():
    os.makedirs(OUT, exist_ok=True)
    it = bocola_italy()
    rep = {}

    # ── 1. verification: can the published construction be recreated? ─────────
    lev_err = float(np.abs(it.lev - it.lev_xls).max())
    sp_err = float(np.abs(it.spread - it.spread_xls).dropna().max())
    mu_err = float(np.abs(it.mu - it.mu_xls).dropna().max())
    mu_errR = float(np.abs(mu_of(it.spread, it.lev, divide_by_R=True)
                           - it.mu_xls).dropna().max())
    # Taking the spread as the archive stores it, (B2) is exact. The residual
    # above is the archive's own rounding of the rate LEVELS to six decimals.
    mu_err0 = float(np.abs(mu_of(it.spread_xls, it.lev, divide_by_R=False)
                           - it.mu_xls).dropna().max())
    print("=== 1. BOCOLA RECONSTRUCTION (his archive, his definitions) ===")
    print(f"  leverage identity  TA/(TA-TL+F51+F52)   max abs err {lev_err:.2e}")
    print(f"  Euribor-Eurepo 1w  difference of levels max abs err {sp_err:.2e}")
    print(f"  multiplier eq.(B2) archive convention   max abs err {mu_err:.2e}")
    print(f"                     published /R_t form  max abs err {mu_errR:.2e}")
    print(f"  multiplier from the STORED spread column max abs err {mu_err0:.2e}"
          "   <- (B2) is exact; the rows above carry the archive's 6-dp rounding "
          "of the rate levels")
    rep["verification"] = {"lev_max_err": lev_err, "spread_max_err": sp_err,
                           "mu_max_err_archive": mu_err, "mu_max_err_paper": mu_errR,
                           "mu_max_err_given_stored_spread": mu_err0}

    print("\n  Italian moments by window (spread annualised bp; lev; mu):")
    rep["italy"] = {}
    for name, win in WINDOWS.items():
        if name.startswith("greek_precrisis_1999"):
            continue
        sl = it.loc[(it.index >= win[0]) & (it.index <= win[1])]
        sp = float(sl.spread.mean() * 4e4)
        lv = float(sl.lev.mean())
        mu = float(sl.mu.mean())
        rep["italy"][name] = {"spread_bp_pa": sp, "lev_mean": lv, "mu_mean": mu,
                              "mu_median": float(sl.mu.median()), "n": int(len(sl))}
        print(f"    {name:34s} spread {sp:6.2f}bp/yr   lev {lv:6.3f}   "
              f"mu mean {mu:.6f} median {sl.mu.median():.6f}")
    print(f"\n  Table 1 lev_bg = {BOCOLA_TABLE1['lev_bg'][0]:.3f} vs measured Italian "
          f"flow-of-funds leverage {it.lev.mean():.3f} (1995Q1-2012Q4). "
          "Table 1 is a literature convention, not this series.")
    print(f"  eq.(B2) at lev=5 and the 8bp/yr sample spread: "
          f"mu = {mu_of(8.19e-4/4, 5.0):.6f}  "
          f"vs estimated posterior mean {BOCOLA_MU_POSTERIOR['mean']:.4f}")

    # ── 2. the Greek counterparts ────────────────────────────────────────────
    df = build_series()
    print("\n=== 2. LEVERAGE: Bocola's definition applied to Greece ===")
    cmp = pd.concat({"bocola": it.lev,
                     "eurostat": pd.to_numeric(df.lev_eurostat_IT, errors="coerce")},
                    axis=1).dropna()
    corr = float(cmp.bocola.corr(cmp.eurostat))
    print(f"  comparability test, Italy: Bank of Italy (Bocola) vs Eurostat "
          f"{EUROSTAT_SECTOR}\n    corr {corr:.4f}, means {cmp.bocola.mean():.3f} vs "
          f"{cmp.eurostat.mean():.3f}, max abs gap {np.abs(cmp.bocola-cmp.eurostat).max():.3f}")
    rep["comparability_italy"] = {"corr": corr, "mean_bocola": float(cmp.bocola.mean()),
                                  "mean_eurostat": float(cmp.eurostat.mean())}

    rep["leverage"] = {}
    for col in ("lev_eurostat_GR", "lev_eurostat_IT", "lev_bsi_GR", "lev_bsi_IT"):
        rep["leverage"][col] = {k: window_stats(df, col, w) for k, w in WINDOWS.items()}
    print("\n  leverage by window:")
    hdr = f"    {'window':34s}" + "".join(f"{c.replace('lev_',''):>18s}"
                                          for c in rep["leverage"])
    print(hdr)
    for name in WINDOWS:
        row = f"    {name:34s}"
        for col in rep["leverage"]:
            st = rep["leverage"][col][name]
            row += f"{st['mean']:10.2f}/{st['median']:<7.2f}" if st else f"{'-':>18s}"
        print(row)
    print("    (cells are mean/median)")
    # Where measured equity crosses zero the ratio is not a leverage reading at
    # all, and neither mean nor median of the window means anything. Say so.
    rep["negative_equity_quarters"] = {}
    for geo, tag in (("EL", "GR"), ("IT", "IT")):
        lv = eurostat_leverage(geo)
        bad = list(lv.index[lv.equity <= 0])
        rep["negative_equity_quarters"][tag] = bad
        thin = lv.index[(lv.equity / lv.A) < 0.03]
        print(f"    {tag}: measured equity <= 0 in {bad if bad else 'no quarter'}; "
              f"below 3% of assets in {len(thin)} quarters"
              + (f" ({thin[0]}-{thin[-1]})" if len(thin) else ""))

    # ── 3. the wedge ─────────────────────────────────────────────────────────
    print("\n=== 3. FUNDING WEDGE: annualised bp ===")
    rep["wedge_bp_pa"] = {}
    cols = ["spread_1w_q", "wedge_depGR_over_eonia", "wedge_depGR_over_depDE",
            "wedge_depIT_over_eonia"]
    for col in cols:
        rep["wedge_bp_pa"][col] = {}
        for name, w in WINDOWS.items():
            st = window_stats(df, col, w)
            rep["wedge_bp_pa"][col][name] = ({k: v * 4e4 for k, v in st.items()
                                              if k != "n"} | {"n": st["n"]}) if st else None
        line = f"  {col:26s}"
        for name in WINDOWS:
            st = rep["wedge_bp_pa"][col][name]
            line += f"{st['mean']:9.1f}" if st else f"{'-':>9s}"
        print(line)
    print("  " + " " * 26 + "".join(f"{n.split('_')[-1][:8]:>9s}" for n in WINDOWS))

    # ── 4. the implied multiplier, Greece ────────────────────────────────────
    print("\n=== 4. IMPLIED MULTIPLIER (eq. B2, published /R form) ===")
    rep["mu_greece"] = {}
    for lev_col, sp_col, tag in [
            ("lev_eurostat_GR", "spread_1w_q", "GR lev x euro-area 1w spread"),
            ("lev_eurostat_GR", "wedge_depGR_over_eonia", "GR lev x GR deposit premium"),
            ("lev_eurostat_GR", "wedge_depGR_over_depDE", "GR lev x GR-DE deposit gap"),
            ("lev_eurostat_IT", "spread_1w_q", "IT lev x euro-area 1w spread")]:
        s = pd.DataFrame({"lev": pd.to_numeric(df[lev_col], errors="coerce"),
                          "sp": pd.to_numeric(df[sp_col], errors="coerce")}).dropna()
        s["mu"] = mu_of(s.sp, s.lev)
        rep["mu_greece"][tag] = {k: window_stats(
            s.assign(**{"mu": s.mu}), "mu", w) for k, w in WINDOWS.items()}
        line = f"  {tag:32s}"
        for name in WINDOWS:
            st = rep["mu_greece"][tag][name]
            line += f"{st['mean']:11.6f}" if st else f"{'-':>11s}"
        print(line)
    print("  " + " " * 32 + "".join(f"{n.split('_')[-1][:9]:>11s}" for n in WINDOWS))

    # ── 5. the model mapping ─────────────────────────────────────────────────
    print("\n=== 5. MODEL MAPPING: mu-bar, 1/mu, (1-mu)/mu ===")
    rows = []
    for tag, theta, sigma_pa in [
            ("paper as written (GK wedge)", 5.0, 100.0),
            ("Bocola Table 1 + 8bp/yr", 5.0, 8.19),
            ("Bocola estimated mu_bg", 5.0, None),
            ("Greek pre-crisis lev + 8bp/yr", 7.98, 8.19),
            ("Greek pre-crisis lev + GR dep premium", 7.98, None),
    ]:
        if tag == "Bocola estimated mu_bg":
            mu = BOCOLA_MU_POSTERIOR["mean"]
            sigma_pa = wedge_of(mu, theta) * 4e4
        elif sigma_pa is None:
            st = rep["wedge_bp_pa"]["wedge_depGR_over_eonia"]["greek_precrisis_2002Q2_2007Q2"]
            sigma_pa = st["mean"]
            mu = mu_of(sigma_pa / 4e4, theta)
        else:
            mu = mu_of(sigma_pa / 4e4, theta)
        rows.append({"calibration": tag, "theta": theta, "sigma_bp_pa": sigma_pa,
                     "mu": mu, "inv_mu": 1 / mu, "odds": (1 - mu) / mu,
                     "sigma_theta_bp_pa": sigma_pa * theta})
    t = pd.DataFrame(rows).set_index("calibration")
    print(t.round(6).to_string())
    rep["mapping"] = t.reset_index().to_dict("records")

    # ── 6. the sovereign-exposure and maturity moments ───────────────────────
    print("\n=== 6. SOVEREIGN EXPOSURE AND MATURITY (Bocola's exp_bg and pi) ===")
    rep["exposure"] = eba_exposure_moments()
    for k, v in rep["exposure"].items():
        print(f"  {k:28s} net/TA {v['net_over_ta']:7.4f}  gross/TA {v['gross_over_ta']:7.4f}"
              f"  avg maturity {v['avg_maturity_q']:6.2f}q  implied pi {1/v['avg_maturity_q']:.4f}")
    print(f"  Bocola Table 1: exp_bg = {BOCOLA_TABLE1['exp_bg'][0]:.3f}, "
          f"pi = {BOCOLA_TABLE1['pi'][0]:.3f} (1/pi = {1/BOCOLA_TABLE1['pi'][0]:.1f}q)")

    with open(os.path.join(OUT, "bocola_greece_moments.json"), "w") as f:
        json.dump({"bocola_table1": {k: {"value": v[0], "source": v[1]}
                                     for k, v in BOCOLA_TABLE1.items()},
                   "bocola_mu_posterior": BOCOLA_MU_POSTERIOR,
                   "windows": WINDOWS, **rep}, f, indent=2)
    print(f"\n  wrote {os.path.join(OUT, 'bocola_greece_moments.json')}")
    return rep


if __name__ == "__main__":
    main()
