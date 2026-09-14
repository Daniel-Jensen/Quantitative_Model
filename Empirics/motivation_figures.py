"""Motivating evidence for the introduction: the Greek episode, and the
disappearance of euro-area sovereign yield decoupling.

Two figures, both built from published data and cached to Empirics/outputs/:

  fig_greece_motivation      2x2 panel -- debt level, debt/GDP, yields, real
                             activity -- the BFT-style motivation panel, Greece only.
  fig_euro_yield_decoupling  10-year yields for four periphery sovereigns against
                             the Bund, 1995-2026, with the ECB/OMT/TPI markers.

Sources
  Eurostat gov_10q_ggdebt   general government consolidated gross debt (Maastricht),
                            sector S13, na_item GD, quarterly, geo EL
  Eurostat namq_10_gdp      B1GQ and P51G, chain-linked volumes (2010), seasonally
                            and calendar adjusted, quarterly, geo EL
  FRED IRLTLT01xxM156N      10-year benchmark government bond yields, monthly,
                            cached by Empirics/graph_spreads.py to spreads_fred.csv

Run:  /opt/anaconda3/envs/ssj/bin/python Empirics/motivation_figures.py [--refresh]
"""

import json
import os
import sys
import textwrap
import urllib.request

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
MACRO = os.path.join(OUT, "greece_macro.csv")
SPREADS = os.path.join(OUT, "spreads_fred.csv")

# Okabe-Ito, assigned in fixed order and never cycled. Validated for CVD
# separation (worst adjacent pair dE 11.0 deutan) by the dataviz validator.
# Germany is deliberately NOT a categorical hue: it is the benchmark, so it wears
# ink and a dashed stroke.
COL = {"ITA": "#0072B2", "ESP": "#E69F00", "GRC": "#D55E00", "PRT": "#009E73"}
INK = "#1a1a1a"
GRID = "#d8d8d8"
SHADE = "#c8c8c8"

# The calibration window: EBA base date 2010-12-31 sits inside it, and it closes
# with the March 2012 PSI exchange.
CAL0, CAL1 = pd.Timestamp("2010-01-01"), pd.Timestamp("2012-06-30")

CAPTIONS = {
    "fig_greece_motivation":
        "The Greek episode in four series: the March 2012 PSI cut the debt stock by "
        "EUR 75bn and the debt ratio by 34 points, yet the ratio was back at its "
        "pre-exchange level within five quarters and the stock not until 2021, while "
        "the sovereign yield reached 29% and real investment fell 70% against a 27% "
        "fall in GDP -- the asymmetry a model of intermediated sovereign risk has to "
        "reproduce.",
    "fig_euro_yield_decoupling":
        "Ten-year sovereign yields converged to within a few tens of basis points of "
        "the Bund after 1999, decoupled violently in 2010-12, and have been "
        "re-compressed ever since OMT: the phenomenon TPI exists to prevent is one "
        "the post-2012 data no longer contains, which is why the instrument's "
        "effects cannot be estimated and have to be modelled.",
}


# ---------------------------------------------------------------- data ------
def _eurostat(dataset, **kw):
    url = ("https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/"
           + dataset + "?format=JSON&lang=EN&"
           + "&".join(f"{k}={v}" for k, v in kw.items()))
    with urllib.request.urlopen(url, timeout=90) as r:
        d = json.load(r)
    inv = {v: k for k, v in d["dimension"]["time"]["category"]["index"].items()}
    s = pd.Series({inv[int(i)]: v for i, v in d["value"].items()}).sort_index()
    s.index = pd.PeriodIndex(s.index.str.replace("-Q", "Q"), freq="Q")
    return s.astype(float)


def greek_macro(refresh=False):
    if os.path.exists(MACRO) and not refresh:
        df = pd.read_csv(MACRO, index_col=0)
        df.index = pd.PeriodIndex(df.index, freq="Q")
        return df
    debt = dict(geo="EL", sector="S13", na_item="GD")
    nat = dict(geo="EL", s_adj="SCA", unit="CLV10_MEUR")
    df = pd.DataFrame({
        "debt_eurbn":     _eurostat("gov_10q_ggdebt", unit="MIO_EUR", **debt) / 1e3,
        "debt_pct_gdp":   _eurostat("gov_10q_ggdebt", unit="PC_GDP", **debt),
        "inv_real_eurbn": _eurostat("namq_10_gdp", na_item="P51G", **nat) / 1e3,
        "gdp_real_eurbn": _eurostat("namq_10_gdp", na_item="B1GQ", **nat) / 1e3,
    })
    df.index.name = "quarter"
    df.to_csv(MACRO)
    return df


def yields():
    y = pd.read_csv(SPREADS, parse_dates=["date"]).set_index("date")
    return y[["ITA", "ESP", "GRC", "PRT", "DEU"]]


# ---------------------------------------------------------------- style -----
def style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "axes.labelsize": 10, "axes.titlesize": 11,
        "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
        "axes.edgecolor": "#999999", "axes.linewidth": 0.8,
        "figure.dpi": 300, "savefig.dpi": 300,
    })


def tidy(ax, ylabel):
    """Recessive grid and axes; the data is the only thing with weight."""
    ax.set_ylabel(ylabel)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def window(ax, label=False):
    ax.axvspan(CAL0, CAL1, color=SHADE, alpha=0.45, lw=0, zorder=0)
    if label:
        ax.annotate("calibration\nwindow", xy=(CAL0 + (CAL1 - CAL0) / 2, 0.045),
                    xycoords=("data", "axes fraction"), ha="center", va="bottom",
                    fontsize=7.5, color="#5a5a5a")


def save(fig, name):
    chars = int(fig.get_size_inches()[0] * 14)
    fig.text(0.5, -0.015, textwrap.fill(CAPTIONS[name], width=chars),
             ha="center", va="top", fontsize=8, style="italic", color="0.35")
    path = os.path.join(OUT, f"{name}.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


# ---------------------------------------------------------------- fig 1 -----
def fig_greece_motivation(df, y):
    q = df.loc["2000Q1":"2019Q4"]
    t = q.index.to_timestamp(how="end")
    ym = y.loc["2000":"2019"]

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.2))
    (a, b), (c, d) = axes

    # (a) debt level -- the PSI write-down is a level event, and it is visible
    a.plot(t, q["debt_eurbn"], color=COL["GRC"], lw=1.8)
    window(a, label=True)
    psi = pd.Timestamp("2012-03-31")
    a.annotate("PSI exchange\n$-$\\texteuro75bn".replace("\\texteuro", "€"),
               xy=(psi, 281.5), xytext=(pd.Timestamp("2013-06-30"), 215),
               fontsize=8, color=INK,
               arrowprops=dict(arrowstyle="->", color="#5a5a5a", lw=0.9))
    a.set_title("(a) General government debt", loc="left")
    tidy(a, "€ bn")

    # (b) debt ratio -- and its full retracement
    b.plot(t, q["debt_pct_gdp"], color=COL["GRC"], lw=1.8)
    window(b)
    b.axhline(175.1, color="#5a5a5a", ls=":", lw=0.9)
    b.annotate("pre-exchange level regained\nwithin five quarters",
               xy=(pd.Timestamp("2013-06-30"), 175.5),
               xytext=(pd.Timestamp("2014-03-31"), 152), fontsize=8, color=INK,
               arrowprops=dict(arrowstyle="->", color="#5a5a5a", lw=0.9))
    b.set_title("(b) General government debt, per cent of GDP", loc="left")
    tidy(b, "per cent")

    # (c) yields -- two series, one axis, gap shaded: the gap IS the spread
    c.fill_between(ym.index, ym["DEU"], ym["GRC"], where=ym["GRC"] >= ym["DEU"],
                   color=COL["GRC"], alpha=0.13, lw=0)
    c.plot(ym.index, ym["GRC"], color=COL["GRC"], lw=1.8)
    c.plot(ym.index, ym["DEU"], color=INK, lw=1.3, ls="--")
    window(c)
    c.annotate("Greece", xy=(pd.Timestamp("2014-01-01"), 22.5), fontsize=8.5,
               color=COL["GRC"], ha="center")
    c.annotate("Germany", xy=(pd.Timestamp("2006-06-01"), 1.4), fontsize=8.5,
               color=INK, ha="center")
    c.set_title("(c) Ten-year sovereign yield; shaded area is the spread", loc="left")
    tidy(c, "per cent")

    # (d) real activity, indexed -- investment against output, same units, one axis
    base = q.loc["2007Q4"]
    d.plot(t, 100 * q["gdp_real_eurbn"] / base["gdp_real_eurbn"],
           color=INK, lw=1.3, ls="--")
    d.plot(t, 100 * q["inv_real_eurbn"] / base["inv_real_eurbn"],
           color=COL["ITA"], lw=1.8)
    window(d)
    d.axhline(100, color="#bbbbbb", lw=0.8, zorder=0)
    d.annotate("real investment", xy=(pd.Timestamp("2015-06-01"), 27), fontsize=8.5,
               color=COL["ITA"], ha="center")
    d.annotate("real GDP", xy=(pd.Timestamp("2015-06-01"), 84), fontsize=8.5,
               color=INK, ha="center")
    d.set_title("(d) Real GDP and investment, 2007Q4 = 100", loc="left")
    tidy(d, "index")

    for ax in (a, b, c, d):
        ax.set_xlim(pd.Timestamp("2000-01-01"), pd.Timestamp("2019-12-31"))
    fig.tight_layout(h_pad=2.0, w_pad=2.4)
    save(fig, "fig_greece_motivation")


# ---------------------------------------------------------------- fig 2 -----
def fig_euro_yield_decoupling(y):
    ym = y.loc["1995":]
    fig, ax = plt.subplots(figsize=(9.2, 4.3))

    for c in ("ITA", "ESP", "PRT", "GRC"):
        ax.plot(ym.index, ym[c], color=COL[c], lw=1.5, label=c)
    ax.plot(ym.index, ym["DEU"], color=INK, lw=1.3, ls="--", label="DEU (benchmark)")

    for name, when, ha in (("ECB created", "1999-01-01", "left"),
                           ("OMT", "2012-09-06", "left"),
                           ("TPI", "2022-07-21", "right")):
        x = pd.Timestamp(when)
        ax.axvline(x, color="#8a8a8a", lw=0.9, ls=(0, (4, 3)), zorder=0)
        ax.annotate(name, xy=(x, 0.97), xycoords=("data", "axes fraction"),
                    ha=ha, va="top", fontsize=8.5, color="#4a4a4a",
                    xytext=(3 if ha == "left" else -3, 0), textcoords="offset points")

    ax.annotate("convergence", xy=(pd.Timestamp("2003-06-01"), 6.1), fontsize=9,
                color="#4a4a4a", ha="center")
    ax.annotate("decoupling", xy=(pd.Timestamp("2010-01-01"), 22.0), fontsize=9,
                color="#4a4a4a", ha="center")
    ax.annotate("re-compression", xy=(pd.Timestamp("2018-06-01"), 8.0), fontsize=9,
                color="#4a4a4a", ha="center")

    tidy(ax, "per cent")
    ax.set_xlim(pd.Timestamp("1995-01-01"), ym.index.max())
    ax.set_ylim(0, 30)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}"))
    ax.legend(loc="upper left", bbox_to_anchor=(0.008, 0.86), frameon=False, ncol=1)
    fig.tight_layout()
    save(fig, "fig_euro_yield_decoupling")


def main():
    style()
    df = greek_macro(refresh="--refresh" in sys.argv)
    y = yields()
    print("building motivation figures ...")
    fig_greece_motivation(df, y)
    fig_euro_yield_decoupling(y)

    emitted = {f[:-4] for f in os.listdir(OUT)
               if f.startswith("fig_") and f.endswith(".png")}
    assert emitted == set(CAPTIONS), (emitted - set(CAPTIONS), set(CAPTIONS) - emitted)
    print(f"caption coverage OK: {len(emitted)} figures, {len(CAPTIONS)} captions")


if __name__ == "__main__":
    main()
