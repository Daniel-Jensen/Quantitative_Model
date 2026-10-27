"""Greek counterpart to Bocola (2016)'s CDS -> physical default probability.

Run:  /opt/anaconda3/envs/ssj/bin/python Empirics/greek_cds/greek.py

What this does and does not do
------------------------------
Bocola's conversion (his eq. 19) is

    p^d_t = CDS_t * CORRECTION_t / Delta ,

so it splits cleanly into three ingredients:

  * CDS_t          the 1-year sovereign CDS spread as a quarterly rate.  For
                   Greece this is NOT obtainable from any open source; see
                   ``cds_input.py``.  Everything else here is computable.
  * CORRECTION_t   the change of measure from the banker's risk-adjusted
                   probability to the physical one.  Country-specific through
                   bank leverage and consumption growth only.  Built here for
                   Greece from Eurostat on Bocola's own definitions.
  * Delta          the haircut.  Bocola sets 0.55 off the 2012 Greek PSI.  The
                   TPI model's ``recovery_rate = 0.30`` implies 0.70, so both
                   are reported: the conversion must use the haircut that the
                   model it feeds actually assumes.

Because the CDS leg is missing, the script inverts the relation instead and
reports the BREAK-EVEN CDS spread: the 1-year Greek sovereign CDS level that
Bocola's transformation would map into the model's ``def_rate_ss`` of 0.10% per
quarter.  That is a sufficient statistic - once the Greek CDS series arrives,
comparing it with this number answers the calibration question directly.

A bond-based cross-check is computed in full, including the Italian
bond-versus-CDS wedge, which is measurable because Italy has both objects.
"""

import io
import json
import os
import urllib.request

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from greek_cds import bocola as B  # noqa: E402
from greek_cds import cds_input    # noqa: E402

OUT = os.path.join(B.ROOT, "Empirics", "outputs")
SERIES_CSV = os.path.join(OUT, "greek_cds_series.csv")
JSON_OUT = os.path.join(OUT, "greek_cds_bocola.json")

# pi_d = 0.10% per quarter is the priced one-quarter default probability at the
# global model's steady state (branch eba-global-calibration; its own stochastic
# rest point sits at 0.001001). That model already carries Bocola's haircut,
# recovery_rate_D = 0.45, so Delta = 0.55 is the internally consistent conversion
# and is reported first. The TPI model on the working branch instead sets
# recovery_rate = 0.30, i.e. Delta = 0.70, which is reported alongside because a
# CDS spread maps into a probability only jointly with a haircut.
MODEL_DEF_RATE_SS = 0.0010
MODEL_HAIRCUT = 0.70

# Bocola's window, and the Greek windows the data actually support.
WINDOWS = {
    "pre_crisis_2002Q2_2007Q2": ("2002Q2", "2007Q2"),
    "cds_era_2006Q1_2007Q2": ("2006Q1", "2007Q2"),
    "bocola_window_2002Q2_2009Q4": ("2002Q2", "2009Q4"),
    "crisis_2010Q1_2011Q4": ("2010Q1", "2011Q4"),
    "full_2002Q2_2011Q4": ("2002Q2", "2011Q4"),
}

# Windows on which a GREEK CDS observation exists at all. The series opens
# 2008-07-21 and ends at the PSI credit event, so there is no pre-crisis
# quarter, and the only quarter in which the leverage-based correction is also
# still meaningful is 2008Q3.
CDS_WINDOWS = {
    "gfc_2008Q3_2009Q4": ("2008Q3", "2009Q4"),
    "sovereign_crisis_2010Q1_2011Q4": ("2010Q1", "2011Q4"),
    "full_cds_2008Q3_2012Q1": ("2008Q3", "2012Q1"),
}

# Greek leverage pre-crisis, and the correction it delivers. Used wherever the
# contemporaneous Greek leverage series has stopped being a solvency reading.
PRE_CRISIS_CORRECTION = 0.9039

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
EUROSTAT = ("https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/")
# S122 alone is published for neither EL nor IT; this is the closest match to
# Bocola's "monetary and financial institutions".  Same choice as
# Empirics/bocola_greece/moments.py, and the same documented cost.
SECTOR = "S121_S122_S123"


def _eurostat(dataset, **q):
    qs = "?format=JSON&lang=EN&" + "&".join(f"{k}={v}" for k, v in q.items())
    with urllib.request.urlopen(EUROSTAT + dataset + qs, timeout=180) as r:
        d = json.load(r)
    inv = {v: k for k, v in d["dimension"]["time"]["category"]["index"].items()}
    s = pd.Series({inv[int(k)]: v for k, v in d["value"].items()}).sort_index()
    s.index = pd.PeriodIndex([i.replace("-Q", "Q") for i in s.index], freq="Q")
    return s


def eurostat_leverage(geo):
    """Bocola's flow-of-funds identity: TA / (TA - TL + F51_liab + F52_liab)."""
    f = lambda item, pos: _eurostat(  # noqa: E731
        "nasq_10_f_bs", unit="MIO_EUR", geo=geo, sector=SECTOR,
        na_item=item, finpos=pos)
    A, L, E = f("F", "ASS"), f("F", "LIAB"), f("F51", "LIAB")
    M = f("F52", "LIAB").reindex(A.index).fillna(0.0)
    return A / (A - L + E + M)


def eurostat_consumption_growth(geo):
    """Log growth of real private final consumption, chain-linked, SCA."""
    c = _eurostat("namq_10_gdp", geo=geo, na_item="P31_S14_S15",
                  unit="CLV10_MEUR", s_adj="SCA")
    return np.log(c).diff()


def ecb_long_rate(geo):
    """Harmonised long-term (10y) government bond yield, monthly, per cent."""
    key = f"IRS/M.{geo}.L.L40.CI.0000.EUR.N.Z"
    u = (f"https://data-api.ecb.europa.eu/service/data/{key}"
         "?format=csvdata&startPeriod=1999-01&endPeriod=2013-12")
    with urllib.request.urlopen(u, timeout=180) as r:
        d = pd.read_csv(io.BytesIO(r.read()))
    d = d[["TIME_PERIOD", "OBS_VALUE"]].dropna()
    s = pd.Series(d.OBS_VALUE.values, index=pd.PeriodIndex(d.TIME_PERIOD, freq="M"))
    return s.groupby(s.index.asfreq("Q")).mean()


def build(refresh=False):
    if os.path.exists(SERIES_CSV) and not refresh:
        df = pd.read_csv(SERIES_CSV, index_col=0)
        df.index = pd.PeriodIndex(df.index, freq="Q")
        return df
    lm = B.archive_lagrange_multiplier()
    df = pd.DataFrame(index=lm.index)
    df["spread_1w_q"] = lm.spread          # euro-area EURIBOR-EUREPO, 1 week
    df["lev_IT_bocola"] = lm.lev
    df["dc_IT_bocola"] = lm.dc
    df["lev_GR"] = eurostat_leverage("EL").reindex(df.index)
    df["lev_IT_eurostat"] = eurostat_leverage("IT").reindex(df.index)
    df["dc_GR"] = eurostat_consumption_growth("EL").reindex(df.index)
    bund = ecb_long_rate("DE")
    for g, tag in [("GR", "GR"), ("IT", "IT")]:
        df[f"bond_spread_{tag}_bp"] = ((ecb_long_rate(g) - bund) * 100).reindex(df.index)
    df["cds_1y_IT_q"] = B.archive_default_probability().cds.reindex(df.index)
    for tag, tenor in [("1y", "1Y"), ("5y", "5Y")]:
        gr = cds_input.load_greek_cds(tenor)
        it = cds_input.load_italian_cds(tenor)
        df[f"cds_{tag}_GR_q"] = gr.reindex(df.index) if gr is not None else np.nan
        df[f"cds_{tag}_IT_lseg_q"] = it.reindex(df.index) if it is not None else np.nan
    df.to_csv(SERIES_CSV)
    return df


# ---------------------------------------------------------------------------
# The conversion
# ---------------------------------------------------------------------------
def greek_conversion(df):
    coef = B.fit_esdf_projection()
    g = pd.DataFrame(index=df.index)
    g["lev"] = df.lev_GR
    g["dc"] = df.dc_GR
    g["mu"] = B.multiplier(df.spread_1w_q, df.lev_GR)
    g["sdf"] = B.sdf(g.dc, g.lev)
    g["esdf"] = B.esdf_projection(g.dc, g.lev, coef)
    g["correction"] = B.correction(g.lev, g.mu, g.esdf)
    # p^d per quarter = cds_q * correction / haircut
    g["mult_bocola"] = g.correction / B.HAIRCUT
    g["mult_model"] = g.correction / MODEL_HAIRCUT
    # break-even 1-year CDS, annualised basis points, for p^d = 0.10% per quarter
    g["breakeven_cds_bp_bocola"] = 1e4 * 4 * MODEL_DEF_RATE_SS / g.mult_bocola
    g["breakeven_cds_bp_model"] = 1e4 * 4 * MODEL_DEF_RATE_SS / g.mult_model
    return g


def risk_adjusted_prob(cds_q, haircut, survival=True):
    """The banker's risk-adjusted one-quarter default probability from the CDS.

    Bocola's eq. (18) is linear, ``phat = CDS_q / Delta``: the premium is paid
    regardless and the contract is a one-period bill.  That is a small-spread
    approximation.  Requiring instead that the premium be paid only if the
    sovereign survives the quarter,

        CDS_q (1 - phat) = Delta * phat   =>   phat = CDS_q / (Delta + CDS_q),

    which agrees with the linear form to better than 1% at Italian spreads and
    is the only usable version at Greek crisis spreads: the linear formula
    returns a probability above one for any par spread beyond about 2,400bp per
    quarter, and Greek quotes pass that in 2011Q4.
    """
    cds_q = np.asarray(cds_q, float)
    return cds_q / (haircut + cds_q) if survival else cds_q / haircut


def leverage_is_a_steady_state(lev):
    """Greek MFI equity turns negative in 2012Q1, so the ratio changes sign.

    A correction above one would say the physical probability exceeds the
    risk-adjusted one, which cannot happen when default is a high marginal-value
    state.  Any window where that occurs is not a usable reading.
    """
    return bool(np.isfinite(lev).all() and (lev > 0).all() and (lev < 25).all())


def greek_answer(df, g):
    """The physical default probability implied by the Greek CDS series."""
    cds = df.cds_1y_GR_q
    corr_live = g.correction.where(
        g.lev.apply(lambda x: np.isfinite(x) and 0 < x < 25), np.nan)
    out = pd.DataFrame(index=df.index)
    out["cds_1y_bp"] = cds * 4 * 1e4
    out["cds_5y_bp"] = df.get("cds_5y_GR_q", np.nan) * 4 * 1e4
    out["correction_live"] = corr_live
    for tag, haircut in [("D55", B.HAIRCUT), ("D60", 0.60), ("D70", MODEL_HAIRCUT)]:
        phat = risk_adjusted_prob(cds, haircut, survival=True)
        out[f"p_{tag}"] = phat * PRE_CRISIS_CORRECTION
        out[f"phat_{tag}"] = phat
    out["p_D55_linear"] = risk_adjusted_prob(cds, B.HAIRCUT, survival=False) * PRE_CRISIS_CORRECTION
    return out


def window_table(g, df):
    rows = {}
    for name, (a, b) in WINDOWS.items():
        w = g.loc[a:b]
        d = df.loc[a:b]
        rows[name] = {
            "quarters": int(w.correction.notna().sum()),
            "lev_GR_mean": float(w.lev.mean()),
            "leverage_usable": leverage_is_a_steady_state(w.lev.dropna()),
            "correction_mean": float(w.correction.mean()),
            "mult_bocola": float(w.mult_bocola.mean()),
            "mult_model": float(w.mult_model.mean()),
            "breakeven_cds_bp_bocola": float(w.breakeven_cds_bp_bocola.mean()),
            "breakeven_cds_bp_model": float(w.breakeven_cds_bp_model.mean()),
            "bond_spread_GR_bp": float(d.bond_spread_GR_bp.mean()),
            "bond_spread_IT_bp": float(d.bond_spread_IT_bp.mean()),
            "cds_1y_IT_bp": float(1e4 * 4 * d.cds_1y_IT_q.mean()),
        }
    return rows


def italian_benchmark(df):
    """Italy's own answer, and the bond-versus-CDS wedge it identifies."""
    dp = B.archive_default_probability().dropna(subset=["p"])
    out = {}
    for name, (a, b) in WINDOWS.items():
        d = dp.loc[a:b]
        s = df.loc[a:b]
        if len(d) == 0:
            continue
        cds_bp = 1e4 * 4 * d.cds.mean()
        bond_bp = s.bond_spread_IT_bp.mean()
        p_cds = float(d.p.mean())
        p_bond = float((bond_bp / 1e4 / 4 / B.HAIRCUT * d.correction).mean())
        out[name] = {
            "cds_1y_bp": float(cds_bp),
            "bond_10y_bp": float(bond_bp),
            "bond_over_cds": float(bond_bp / cds_bp),
            "p_phys_from_cds_pct_q": 100 * p_cds,
            "p_phys_from_bond_pct_q": 100 * p_bond,
            "ratio_to_model_0p10": p_cds / MODEL_DEF_RATE_SS,
        }
    return out


def main():
    refresh = "--refresh" in sys.argv
    df = build(refresh=refresh)
    g = greek_conversion(df)
    a0, a1, err = B.fit_esdf_projection()
    have_cds = bool(df.cds_1y_GR_q.notna().any())

    res = {
        "parameters": {"lambda": B.LAMBDA, "haircut_bocola": B.HAIRCUT,
                       "haircut_model_tpi": MODEL_HAIRCUT,
                       "haircut_market_convention": 0.60,
                       "sdf_b": B.SDF_B, "sdf_omega": B.SDF_OMEGA,
                       "pre_crisis_correction": PRE_CRISIS_CORRECTION,
                       "esdf_projection": {"a0": a0, "a1": a1, "max_abs_err": err}},
        "greece_conversion": window_table(g, df),
        "italy_benchmark": italian_benchmark(df),
        "greek_cds_available": have_cds,
    }

    print(f"E[Lhat|d=1] projection: a0={a0:.6f} a1={a1:.6f}  max abs err on his "
          f"own series {err:.2e} ({100*err/2.0:.2f}% of level)\n")

    if not have_cds:
        print(cds_input.MISSING_DATA_NOTE)
    else:
        res["vendor_validation"] = cds_input.validate_against_bocola(
            B.archive_default_probability().cds)
        res["market_recovery"] = cds_input.market_recovery()
        res["staleness"] = cds_input.staleness_report()
        ans = greek_answer(df, g)

        v = res["vendor_validation"]
        print("VENDOR VALIDATION - Italian 1-year, this composite vs Bocola's "
              f"Markit/WRDS series over {v['quarters']} overlapping quarters:")
        print(f"  correlation {v['corr']:.4f}   mean abs gap {v['mean_abs_pct_gap']:.1f}%"
              f"   level ratio {v['level_ratio']:.3f}\n")

        st = res["staleness"]["GR_1Y"]
        print(f"GREEK 1-YEAR CDS COVERAGE: {st['quotes']} quoted days "
              f"{st['first']} -> {st['last']}, {100*st['share_unchanged']:.1f}% unchanged "
              f"from the previous day. No pre-crisis observation exists.\n")

        print("GREECE - physical one-quarter default probability, per cent")
        print(f"{'quarter':9s}{'1y CDS bp':>11s}{'5y CDS bp':>11s}{'corr':>7s}"
              f"{'p D=.55':>9s}{'p D=.60':>9s}{'p D=.70':>9s}{'p linear':>10s}")
        a = ans.loc["2008Q3":"2012Q1"]
        for q, r in a.iterrows():
            c = f"{r.correction_live:.3f}" if np.isfinite(r.correction_live) else " n/a "
            print(f"{str(q):9s}{r.cds_1y_bp:11.1f}{r.cds_5y_bp:11.1f}{c:>7s}"
                  f"{100*r.p_D55:9.3f}{100*r.p_D60:9.3f}{100*r.p_D70:9.3f}"
                  f"{100*r.p_D55_linear:10.3f}")
        res["greece_answer"] = {}
        print()
        for n, (x, y) in CDS_WINDOWS.items():
            w = ans.loc[x:y]
            m = float(w.p_D55.mean())
            res["greece_answer"][n] = {
                "cds_1y_bp": float(w.cds_1y_bp.mean()),
                "p_phys_pct_q_D55": 100 * m,
                "ratio_to_model_0p10": m / MODEL_DEF_RATE_SS,
            }
            print(f"  {n:32s} mean 1y CDS {w.cds_1y_bp.mean():8.1f}bp   "
                  f"p = {100*m:7.3f}%/q   = {m/MODEL_DEF_RATE_SS:7.1f} x the "
                  f"model's 0.10%")

    print("\nGREECE - the CDS level that WOULD deliver 0.10%/quarter")
    print(f"{'window':30s}{'n':>3s}{'lev':>7s}{'ok':>4s}{'corr':>7s}"
          f"{'bp D=.55':>10s}{'bp D=.70':>10s}{'10y GR':>8s}")
    for n, r in res["greece_conversion"].items():
        print(f"{n:30s}{r['quarters']:3d}{r['lev_GR_mean']:7.2f}"
              f"{'y' if r['leverage_usable'] else 'NO':>4s}{r['correction_mean']:7.3f}"
              f"{r['breakeven_cds_bp_bocola']:10.1f}{r['breakeven_cds_bp_model']:10.1f}"
              f"{r['bond_spread_GR_bp']:8.1f}")

    print("\nITALY - his own answer, and how far a bond spread overstates it")
    print(f"{'window':30s}{'1y CDS':>9s}{'10y':>9s}{'ratio':>7s}"
          f"{'p(CDS)%/q':>11s}{'p(bond)%/q':>12s}{'xmodel':>8s}")
    for n, r in res["italy_benchmark"].items():
        print(f"{n:30s}{r['cds_1y_bp']:9.1f}{r['bond_10y_bp']:9.1f}"
              f"{r['bond_over_cds']:7.2f}{r['p_phys_from_cds_pct_q']:11.4f}"
              f"{r['p_phys_from_bond_pct_q']:12.4f}{r['ratio_to_model_0p10']:8.2f}")

    os.makedirs(OUT, exist_ok=True)
    with open(JSON_OUT, "w") as f:
        json.dump(res, f, indent=2, default=float)
    frames = [df, g.add_prefix("gr_")]
    if have_cds:
        frames.append(greek_answer(df, g).add_prefix("ans_"))
    pd.concat(frames, axis=1).to_csv(SERIES_CSV)
    print(f"\nwrote {JSON_OUT}\nwrote {SERIES_CSV}")


if __name__ == "__main__":
    main()
