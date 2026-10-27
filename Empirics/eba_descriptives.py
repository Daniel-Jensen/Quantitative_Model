"""Empirical discipline for the TPI paper: stocks, shock scale, transmission.

Purpose (paper role)
--------------------
Three calibration/validation blocks, none of them a causal claim:

  1. STOCKS.   Who held Greek sovereign paper at 31 Dec 2010 and how much of it
               they held against their own capital. Pins the model's
               sovereign-to-net-worth ratios and the domestic/foreign split of
               the D sovereign stock.
  2. SHOCK.    What the observed 2011 repricing did to a FROZEN 31-Dec-2010
               portfolio, as a share of Core Tier 1. Answers whether the
               net-worth impulse the model feeds itself is of the right scale.
  3. TRANSMISSION. Greek versus German corporate borrowing costs over the
               risk-free rate: the observable counterpart to the model's credit
               spread r_wc - r. Order of magnitude and cross-country sign only -
               observed lending rates also carry borrower risk, NPLs and
               composition, so they bound the intermediation wedge from above.

Quantities are frozen throughout. Nothing here identifies a sale, and nothing
here is a regression.

Run:
    /opt/anaconda3/envs/ssj/bin/python Empirics/eba_descriptives.py

Outputs (Empirics/outputs/):
    eba_bank_panel.csv          one row per EBA bank
    eba_holder_groups.csv       the holder-group table
    eba_discipline.json         the headline moments, for the paper
    lending_rates.csv           cached ECB MIR / Euribor pull
    fig_eba_ownership.png       who held it, and against how much capital
    fig_eba_ladder.png          residual-maturity ladders and implied duration
    fig_eba_mtm.png             frozen-portfolio loss over Core Tier 1
    fig_lending_spreads.png     GR vs DE corporate borrowing spreads
"""
from __future__ import annotations

import json
import os
import sys
import textwrap
import urllib.request

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "outputs")
sys.path.insert(0, ROOT)
from calibration import eba_moments as E  # noqa: E402  (path set above)

# Okabe-Ito, the assignment already used in motivation_figures.py. Germany is
# the benchmark and wears ink, not a categorical hue; "other" wears grey.
COL = {"GR": "#D55E00", "FR": "#56B4E9", "CY": "#CC79A7", "BE": "#009E73",
       "DE": "#1a1a1a", "OTHER": "#9a9a9a"}
INK, GRID = "#1a1a1a", "#d8d8d8"
NAME = {"GR": "Greek", "FR": "French", "DE": "German", "CY": "Cypriot",
        "BE": "Belgian", "OTHER": "Other", "CORE": "Core bloc"}

# The model's two blocs. F = Germany is the draft's choice; the core aggregate is
# reported alongside it because the two give materially different cross-border
# exposure, and that is a calibration decision rather than a detail.
CORE = ("DE", "FR", "BE", "NL", "AT", "LU", "FI")

# 10-year benchmark yields, decimal. Dec-2010 is the EBA base date.
Y_2010_12, Y_2011_12, Y_2012_02 = 0.1201, 0.2114, 0.2924
# Zettelmeyer-Trebesch-Gulati NPV haircut range on the March 2012 exchange.
PSI_HAIRCUT = (0.59, 0.65)

CAPTIONS = {
    "fig_eba_ownership":
        "Greek sovereign paper was concentrated where capital was thinnest. "
        "Greek and Cypriot banks held 67 per cent of the EBA sample's Greek "
        "book, and held it against capital the position exceeded outright. The "
        "French and German banks held a fifth of it at under a tenth of their "
        "capital. Gross direct long exposure, 31 December 2010; 58 of the 90 "
        "banks report a Greek position.",
    "fig_eba_ladder":
        "Holder groups sat on different parts of one maturity structure. "
        "Nineteen per cent of the Greek banks' Greek book matured within a "
        "year, against 47 per cent of the German banks' beyond ten, giving "
        "Macaulay durations of 3.5 and 4.5 years. That is holder selection, "
        "not two securities: the ladders value the observed portfolios and the "
        "model keeps one Greek sovereign duration. Nothing here is a sale.",
    "fig_eba_mtm":
        "The 2011 repricing alone, applied to a frozen 31-December-2010 "
        "portfolio, removed two-fifths of Greek bank capital and a quarter of "
        "Cypriot; marking at the restructuring's present value removes almost "
        "all of Greek capital. The same arithmetic costs the German and French "
        "banks one to two per cent. No sale, redemption, hedge or accounting "
        "choice enters: the book is held fixed and only the price moves.",
    "fig_eba_structure":
        "Greek sovereign paper sat where capital was thinnest, and holders sat "
        "on different parts of one maturity structure. Greek banks held 60 per "
        "cent of the EBA sample's Greek book at 2.4 times their Core Tier 1; "
        "the German and French banks held a fifth of it at under a tenth of "
        "theirs. Nineteen per cent of the Greek book matured within a year "
        "against 47 per cent of the German beyond ten, so the ladders value "
        "the observed portfolios - they are not two securities, and the model "
        "keeps one Greek sovereign duration.",
    "fig_lending_spreads":
        "The model's credit spread has an observable counterpart, and it moved. "
        "The Greek-German gap in the cost of corporate borrowing over the "
        "risk-free rate widened from 110 basis points before the crisis to 271 "
        "on average over 2010-12H1, and peaked at 388 inside that window; on "
        "overdrafts, the "
        "closest analogue to working capital, from 84 to 193. Read the change, "
        "not the level, and read it as an upper bound: the gap also carries "
        "borrower risk, non-performing loans and loan composition.",
}


# ─────────────────────────────────────────────────────────── panel ──────────
def bank_panel() -> pd.DataFrame:
    """One row per EBA 2011 bank: capital, sovereign exposure, Greek ladder."""
    df = E.load_eba()
    meta = ["WORKSHEET", "INFORMATION_CODE", "INFO_DATE_CODE", "C_COUNTRY_CODE",
            "SCENARIO_CODE", "MATURITY_CODE", "MEASURE_CODE"]
    banks = [c for c in df.columns if c not in meta]

    def line(ws, code, issuer=None, maturity="ANY", date=None, measure=None):
        w = df[(df.WORKSHEET == ws) & (df.INFORMATION_CODE == code)]
        if issuer is not None:
            w = w[w.C_COUNTRY_CODE == issuer]
        if maturity == "NAN":
            w = w[w.MATURITY_CODE.isna()]
        elif maturity != "ANY":
            w = w[w.MATURITY_CODE == maturity]
        if date is not None:
            w = w[w.INFO_DATE_CODE == date]
        if measure is not None:
            w = w[w.MEASURE_CODE == measure]
        w = w[w.SCENARIO_CODE.isna()]
        if len(w) != 1:
            raise ValueError(f"{code}/{issuer}/{maturity}: {len(w)} rows, want 1")
        return pd.to_numeric(w[banks].iloc[0], errors="coerce")

    p = pd.DataFrame({
        "ct1": line(E.WS_AGG, E.CODE_CT1, date=E.DATE_2010, measure=120),
        "rwa": line(E.WS_AGG, E.CODE_RWA, date=E.DATE_2010, measure=120),
        "ta": line(E.WS_AGG, E.CODE_TA, date=E.DATE_2010, measure=120),
        "gr_gross": line(E.WS_SOV, 34010, issuer="GR", maturity=999.0),
        "gr_net": line(E.WS_SOV, 34012, issuer="GR", maturity=999.0),
        "sov_all": line(E.WS_SOV, 34010, issuer="TO", maturity="NAN"),
    })
    p["giips"] = sum(line(E.WS_SOV, 34010, issuer=c, maturity=999.0)
                     for c in ("GR", "IE", "IT", "PT", "ES"))
    # Own-sovereign book, for the home-bias column: issuer == bank nationality.
    p["nat"] = [b[:2] for b in p.index]
    own = {}
    for iso in p.nat.unique():
        try:
            own[iso] = line(E.WS_SOV, 34010, issuer=iso, maturity=999.0)
        except ValueError:
            own[iso] = pd.Series(np.nan, index=p.index)
    p["sov_own"] = [own[n][b] for b, n in zip(p.index, p.nat)]

    for m in E.BUCKETS:
        p[f"gr_{m}"] = line(E.WS_SOV, 34010, issuer="GR", maturity=float(m))

    p["gr_over_ct1"] = p.gr_gross / p.ct1
    p["ct1_ratio"] = p.ct1 / p.rwa
    return p


def group_key(nat: str) -> str:
    return nat if nat in ("GR", "CY", "FR", "DE", "BE") else "OTHER"


# ──────────────────────────────────────────────── ladders and repricing ─────
def ladder_of(p: pd.DataFrame, rows) -> dict[int, float]:
    return {m: float(p.loc[rows, f"gr_{m}"].sum()) for m in E.BUCKETS}


def reprice(lad: dict[int, float], y0: float, y1: float,
            coupon: float = E.COUPON_GGB) -> dict:
    """Frozen-portfolio revaluation: same bonds, two yields."""
    mv0, mac0 = E.ladder_value_and_duration(lad, y0, coupon)
    mv1, _ = E.ladder_value_and_duration(lad, y1, coupon)
    face = sum(lad.values())
    return {"face": face, "mv0": mv0, "mv1": mv1, "loss": mv0 - mv1,
            "loss_frac_mv": (mv0 - mv1) / mv0 if mv0 else np.nan,
            "macaulay_y": mac0,
            "modified_y": E.ladder_modified_duration(lad, y0, coupon),
            "wavg_resid_y": sum(w * np.mean(E.BUCKETS[m])
                                for m, w in lad.items()) / face if face else np.nan}


# ───────────────────────────────────────────── the model-side impulse ──────
# Paper calibration, sections 2-3 of the VIVA draft.
R_Q, RHO_S, REC_D = 0.003, 0.95, 0.45
PI_SS, PI_SHOCK = 0.001, 0.0198
S_BAR = np.log(PI_SS / (1 - PI_SS))


def hm_price(pi0: float, delta_b: float, horizon: int = 400) -> float:
    """Risk-neutral price of the Hatchondo-Martinez perpetuity at hazard pi0.

    The latent factor mean-reverts at RHO_S in log-odds, exactly as in the
    paper's s-process, and a default writes the claim down to REC_D. There is
    no risk premium and no constraint discount here, so this is a LOWER bound
    on the model's price fall: the covariance and lambda*mu terms add to it.
    """
    s, q, surv, w = np.log(pi0 / (1 - pi0)), 0.0, 1.0, 1.0
    for j in range(1, horizon + 1):
        s = (1 - RHO_S) * S_BAR + RHO_S * s
        surv *= 1 - (1 / (1 + np.exp(-s))) * (1 - REC_D)
        q += delta_b * w * surv / (1 + R_Q) ** j
        w *= 1 - delta_b
    return q


def impulse_check(observed_fall: float, delta_b: float = 0.056) -> dict:
    """How far does the paper's risk impulse move the bond price?"""
    q0, q1 = hm_price(PI_SS, delta_b), hm_price(PI_SHOCK, delta_b)
    lo, hi, target = PI_SS * 1.1, 0.30, observed_fall
    for _ in range(80):                       # bisection; monotone in pi
        mid = 0.5 * (lo + hi)
        if 1 - hm_price(mid, delta_b) / q0 < target:
            lo = mid
        else:
            hi = mid
    return {"delta_b": delta_b, "q_at_ss": q0, "q_at_shock": q1,
            "actuarial_price_fall": 1 - q1 / q0,
            "observed_price_fall": observed_fall,
            "pi_matching_observed": 0.5 * (lo + hi)}


# ─────────────────────────────────────────────────────── lending rates ──────
MIR = {  # label -> ECB Data Portal series key
    "GR_cob": "MIR/M.GR.B.A2I.AM.R.A.2240.EUR.N",
    "DE_cob": "MIR/M.DE.B.A2I.AM.R.A.2240.EUR.N",
    "GR_ovd": "MIR/M.GR.B.A2Z.A.R.A.2240.EUR.N",
    "DE_ovd": "MIR/M.DE.B.A2Z.A.R.A.2240.EUR.N",
    "eur3m": "FM/M.U2.EUR.RT.MM.EURIBOR3MD_.HSTA",
}
RATES_CSV = os.path.join(OUT, "lending_rates.csv")


def lending_rates(refresh: bool = False) -> pd.DataFrame:
    """ECB composite cost of borrowing, overdraft rates and 3M Euribor.

    Cached to disk: the figure must rebuild without a network call.
    """
    if os.path.exists(RATES_CSV) and not refresh:
        return pd.read_csv(RATES_CSV, parse_dates=["date"]).set_index("date")
    out = {}
    for name, key in MIR.items():
        url = (f"https://data-api.ecb.europa.eu/service/data/{key}"
               "?format=csvdata&startPeriod=2003-01&endPeriod=2015-12")
        with urllib.request.urlopen(url, timeout=60) as r:
            d = pd.read_csv(r)
        d = d[["TIME_PERIOD", "OBS_VALUE"]].dropna()
        out[name] = pd.Series(d.OBS_VALUE.values,
                              index=pd.PeriodIndex(d.TIME_PERIOD, freq="M"))
        print(f"  pulled {name}: {len(out[name])} obs")
    r = pd.DataFrame(out)
    r.index = r.index.to_timestamp()
    r.index.name = "date"
    r.to_csv(RATES_CSV)
    return r


# ──────────────────────────────────────────────────────────── plotting ──────
def style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "axes.labelsize": 10, "axes.titlesize": 10.5,
        "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
        "axes.edgecolor": "#999999", "axes.linewidth": 0.8,
        "figure.dpi": 300, "savefig.dpi": 300,
    })


def tidy(ax, ylabel=None, xgrid=False):
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(True, axis="x" if xgrid else "y", color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def save(fig, name):
    chars = int(fig.get_size_inches()[0] * 15)
    fig.text(0.5, -0.02, textwrap.fill(CAPTIONS[name], width=chars),
             ha="center", va="top", fontsize=8, style="italic", color="0.35")
    path = os.path.join(OUT, f"{name}.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


def fig_ownership(p, grp):
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.2, 4.1),
                               gridspec_kw={"width_ratios": [1, 1.25]})
    order = [k for k in grp.sort_values("gr_gross", ascending=False).index
             if k != "OTHER"] + ["OTHER"]
    vals = grp.loc[order, "gr_gross"] / 1e3
    a.bar(range(len(order)), vals, width=0.62,
          color=[COL[k] for k in order], zorder=3)
    for i, (k, v) in enumerate(zip(order, vals)):
        a.text(i, v + 1.2, f"{v:.1f}\n{grp.loc[k,'share']*100:.0f}%",
               ha="center", va="bottom", fontsize=8.5, color=INK, linespacing=1.15)
    a.set_xticks(range(len(order)))
    a.set_xticklabels([NAME[k] if k != "OTHER" else "Other\n(13 countries)"
                       for k in order], fontsize=9)
    a.set_ylim(0, vals.max() * 1.30)
    tidy(a, "€ bn, gross direct long")
    a.set_title("Greek sovereign holdings, 31 Dec 2010", loc="left")

    # Bank-level intensity. Log scale: the range is 0.02% to 991% of capital.
    q = p[p.gr_gross > 1].copy()
    q["g"] = q.nat.map(group_key)
    ypos = {k: i for i, k in enumerate(["GR", "CY", "BE", "FR", "DE", "OTHER"])}
    rng = np.random.default_rng(7)
    for k, i in ypos.items():
        s = q[q.g == k]
        b.scatter(s.gr_over_ct1 * 100, i + rng.uniform(-.17, .17, len(s)),
                  s=26, color=COL[k], alpha=.85, lw=.5, edgecolor="white", zorder=3)
    b.axvline(100, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=2)
    b.text(88, -0.8, "exposure = capital", fontsize=8.5, color=INK,
           va="bottom", ha="right")
    b.set_xscale("log")
    b.set_yticks(list(ypos.values()))
    b.set_yticklabels([NAME[k] for k in ypos])
    b.set_xlabel("Greek sovereign exposure, % of Core Tier 1 (log scale)")
    b.grid(True, axis="x", color=GRID, lw=0.6, alpha=0.7)
    b.set_axisbelow(True)
    for side in ("top", "right"):
        b.spines[side].set_visible(False)
    b.set_ylim(5.7, -1.3)
    b.set_title("Each bank with a Greek position", loc="left", pad=12)
    fig.tight_layout()
    save(fig, "fig_eba_ownership")


def fig_ladder(ladders, reps):
    labels = ["0–3m", "3m–1y", "1–2y", "2–3y", "3–5y", "5–10y", "10y+"]
    keys = ["GR", "CY", "FR", "DE"]
    fig, ax = plt.subplots(figsize=(9.2, 3.9))
    w, x = 0.2, np.arange(len(labels))
    for j, k in enumerate(keys):
        lad = ladders[k]
        tot = sum(lad.values())
        sh = [lad[m] / tot * 100 for m in E.BUCKETS]
        ax.bar(x + (j - 1.5) * w, sh, width=w * 0.88, color=COL[k], zorder=3,
               label=f"{NAME[k]} banks  (duration {reps[k]['macaulay_y']:.1f}y)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    tidy(ax, "% of the holder's Greek book")
    ax.set_ylim(0, 54)
    ax.legend(frameon=False, ncol=2, loc="upper left", handlelength=1.1,
              columnspacing=1.4)
    ax.set_title("Residual maturity of the Greek sovereign book, 31 Dec 2010",
                 loc="left")
    fig.tight_layout()
    save(fig, "fig_eba_ladder")


def fig_mtm(mtm):
    keys = ["GR", "CY", "BE", "FR", "DE"]
    marks = [("y2011", "Marked at the Dec-2011 yield (21.1%)", "#f0b799"),
             ("y2012", "Marked at the Feb-2012 peak (29.2%)", "#dd8a5c"),
             ("psi", "Marked at the PSI present value", "#D55E00")]
    fig, ax = plt.subplots(figsize=(9.2, 4.0))
    x, w = np.arange(len(keys)), 0.26
    for j, (mk, lab, c) in enumerate(marks):
        v = [mtm[k][mk] * 100 for k in keys]
        ax.bar(x + (j - 1) * w, v, width=w * 0.9, color=c, zorder=3, label=lab)
        for xi, vi in zip(x + (j - 1) * w, v):
            if vi > 4:
                ax.text(xi, vi + 1.5, f"{vi:.0f}", ha="center", va="bottom",
                        fontsize=8, color="0.3")
    ax.axhline(100, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=2)
    ax.text(-0.45, 102, "capital exhausted", fontsize=8.5, color=INK, ha="left")
    ax.set_ylim(0, 138)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{NAME[k]} banks" for k in keys])
    tidy(ax, "loss on a frozen book, % of Core Tier 1")
    ax.legend(frameon=False, loc="upper right", handlelength=1.2)
    ax.set_title("What the repricing alone did to capital", loc="left")
    fig.tight_layout()
    save(fig, "fig_eba_mtm")


def fig_spreads(r):
    r = r.loc["2004":"2015"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.2, 3.9), sharey=True)
    panels = [(a, "cob", "Composite cost of borrowing"),
              (b, "ovd", "Revolving loans and overdrafts")]
    for ax, suf, title in panels:
        ax.axvspan(pd.Timestamp("2010-01-01"), pd.Timestamp("2012-06-30"),
                   color="#c8c8c8", alpha=0.45, lw=0, zorder=0)
        ax.plot(r.index, r[f"GR_{suf}"] - r.eur3m, color=COL["GR"], lw=1.7,
                label="Greece", zorder=3)
        ax.plot(r.index, r[f"DE_{suf}"] - r.eur3m, color=INK, lw=1.3, ls="--",
                label="Germany", zorder=3)
        ax.axhline(0, color="#bbbbbb", lw=0.8, zorder=1)
        tidy(ax)
        ax.set_title(title, loc="left")
    a.set_ylabel("rate to non-financial firms\nless 3-month Euribor, pp")
    a.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    save(fig, "fig_lending_spreads")


# ──────────────────────────────────────────────────────────────── main ──────
def fig_structure(p, grp, ladders, reps):
    """The paper exhibit: holdings by nationality (A) and maturity ladder (B)."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.4, 3.8),
                               gridspec_kw={"width_ratios": [1, 1.5]})
    order = [k for k in grp.sort_values("gr_gross", ascending=False).index
             if k != "OTHER"] + ["OTHER"]
    vals = grp.loc[order, "gr_gross"] / 1e3
    a.bar(range(len(order)), vals, width=0.62,
          color=[COL[k] for k in order], zorder=3)
    for i, (k, v) in enumerate(zip(order, vals)):
        a.text(i, v + 1.2, f"{v:.1f}\n{grp.loc[k,'share']*100:.0f}%",
               ha="center", va="bottom", fontsize=8.5, color=INK, linespacing=1.15)
    a.set_xticks(range(len(order)))
    a.set_xticklabels([NAME[k] if k != "OTHER" else "Other" for k in order],
                      fontsize=8.5, rotation=30, ha="right")
    a.set_ylim(0, vals.max() * 1.32)
    tidy(a, "\u20ac bn, gross direct long")
    a.set_title("A. Who held it", loc="left")

    labels = ["0\u20133m", "3m\u20131y", "1\u20132y", "2\u20133y",
              "3\u20135y", "5\u201310y", "10y+"]
    keys = ["GR", "CY", "FR", "DE"]
    w, x = 0.2, np.arange(len(labels))
    for j, k in enumerate(keys):
        lad = ladders[k]
        tot = sum(lad.values())
        sh = [lad[m] / tot * 100 for m in E.BUCKETS]
        b.bar(x + (j - 1.5) * w, sh, width=w * 0.88, color=COL[k], zorder=3,
              label=f"{NAME[k]}  ({reps[k]['macaulay_y']:.1f}y)")
    b.set_xticks(x)
    b.set_xticklabels(labels, fontsize=8.5)
    tidy(b, "% of the holder's Greek book")
    b.set_ylim(0, 56)
    b.legend(frameon=False, ncol=2, loc="upper left", handlelength=1.0,
             columnspacing=1.2, fontsize=8.5,
             title="Macaulay duration in brackets", title_fontsize=8.5)
    b.get_legend().get_title().set_color("0.35")
    b.set_title("B. Residual maturity of the Greek book", loc="left")
    fig.tight_layout()
    save(fig, "fig_eba_structure")


def main():
    os.makedirs(OUT, exist_ok=True)
    style()
    p = bank_panel()
    p.to_csv(os.path.join(OUT, "eba_bank_panel.csv"))
    p["g"] = p.nat.map(group_key)

    total_gr = p.gr_gross.sum()
    grp = p.groupby("g").agg(n_banks=("gr_gross", "size"),
                             n_holding=("gr_gross", lambda s: int((s > 1).sum())),
                             gr_gross=("gr_gross", "sum"),
                             gr_net=("gr_net", "sum"),
                             ct1=("ct1", "sum"), ta=("ta", "sum"),
                             giips=("giips", "sum"))
    grp["share"] = grp.gr_gross / total_gr
    grp["gr_over_ct1"] = grp.gr_gross / grp.ct1
    grp["giips_over_ct1"] = grp.giips / grp.ct1
    grp.to_csv(os.path.join(OUT, "eba_holder_groups.csv"))

    ladders = {k: ladder_of(p, p.index[p.g == k]) for k in grp.index}
    ladders["CORE"] = ladder_of(p, p.index[p.nat.isin(CORE)])
    reps = {k: reprice(lad, Y_2010_12, Y_2011_12) for k, lad in ladders.items()}

    # Frozen-portfolio loss over capital, three markings.
    ct1_of = {k: float(grp.loc[k, "ct1"]) for k in grp.index}
    ct1_of["CORE"] = float(p.loc[p.nat.isin(CORE), "ct1"].sum())
    mtm = {}
    for k, lad in ladders.items():
        r11 = reprice(lad, Y_2010_12, Y_2011_12)
        r12 = reprice(lad, Y_2010_12, Y_2012_02)
        mv0, face = r11["mv0"], r11["face"]
        psi_lo = (mv0 - face * (1 - PSI_HAIRCUT[1])) / ct1_of[k]
        psi_hi = (mv0 - face * (1 - PSI_HAIRCUT[0])) / ct1_of[k]
        mtm[k] = {"mv0": mv0, "face": face,
                  "y2011": r11["loss"] / ct1_of[k],
                  "y2012": r12["loss"] / ct1_of[k],
                  "psi": (psi_lo + psi_hi) / 2,
                  "psi_lo": psi_lo, "psi_hi": psi_hi}

    # ── the model-side translation ──────────────────────────────────────────
    phi_D = float(grp.loc["GR", "gr_over_ct1"])          # D-bank sovereign / n
    phi_F_de = float(grp.loc["DE", "gr_over_ct1"])
    phi_F_core = float(p.loc[p.nat.isin(CORE), "gr_gross"].sum() / ct1_of["CORE"])
    gr_debt_2010 = 330570.0   # Eurostat gov_10dd_edpt1, GD, S13, EL, 2010, EUR m
    moments = {
        "source": "EBA 2011 EU-wide stress test, sovereign exposures 31 Dec 2010",
        "stocks": {
            "n_banks": int(len(p)), "n_holding_GR": int((p.gr_gross > 1).sum()),
            "GR_sov_held_by_sample_EURbn": total_gr / 1e3,
            "GR_general_govt_debt_2010_EURbn": gr_debt_2010 / 1e3,
            "sample_share_of_GR_debt": total_gr / gr_debt_2010,
            "domestic_share_of_sample_GR_book": float(grp.loc["GR", "share"]),
            "foreign_share_of_sample_GR_book": 1 - float(grp.loc["GR", "share"]),
            "net_over_gross_GR": float(p.gr_net.sum() / total_gr),
            "phi_bD_D_sovereign_over_CT1": phi_D,
            "phi_bD_F_Germany_only": phi_F_de,
            "phi_bD_F_core_bloc": phi_F_core,
        },
        "shock": {
            "yields": {"2010-12": Y_2010_12, "2011-12": Y_2011_12,
                       "2012-02": Y_2012_02},
            "loss_over_CT1": {k: {m: mtm[k][m] for m in
                                  ("y2011", "y2012", "psi_lo", "psi_hi")}
                              for k in mtm},
            "macaulay_duration_y": {k: reps[k]["macaulay_y"] for k in reps},
            "price_fall_2011_frac": reps["GR"]["loss_frac_mv"],
            "note": ("A frozen 31-Dec-2010 book revalued at the 10-year "
                     "benchmark yield. The flat-curve assumption UNDERSTATES "
                     "the 2011H2 loss: the Greek curve inverted and short paper "
                     "traded on recovery, not on yield."),
        },
    }
    with open(os.path.join(OUT, "eba_discipline.json"), "w") as f:
        json.dump(moments, f, indent=2)

    # ── report ──────────────────────────────────────────────────────────────
    pd.set_option("display.width", 200)
    print("\n=== 1. STOCKS: who held Greek sovereign paper, 31 Dec 2010 ===")
    t = grp.loc[["GR", "CY", "FR", "DE", "BE", "OTHER"],
                ["n_banks", "n_holding", "gr_gross", "share", "ct1",
                 "gr_over_ct1", "giips_over_ct1"]].copy()
    t["gr_gross"] /= 1e3
    t["ct1"] /= 1e3
    t.columns = ["banks", "holding", "GGB €bn", "share", "CT1 €bn",
                 "GGB/CT1", "GIIPS/CT1"]
    print(t.round(3).to_string())
    print(f"\n  sample holds €{total_gr/1e3:.1f}bn = "
          f"{total_gr/gr_debt_2010*100:.1f}% of Greek general government debt")
    print(f"  domestic share of the sample's Greek book : "
          f"{grp.loc['GR','share']*100:.1f}%")
    print(f"  net / gross direct long                   : "
          f"{p.gr_net.sum()/total_gr*100:.1f}%")
    print(f"  sovereign/CT1  D = Greek banks            : {phi_D:.2f}")
    print(f"                 F = German banks           : {phi_F_de:.3f}")
    print(f"                 F = core bloc (7 countries): {phi_F_core:.3f}")

    print("\n=== 2. SHOCK: frozen 31-Dec-2010 book, loss as % of Core Tier 1 ===")
    rows = []
    for k in ["GR", "CY", "BE", "FR", "DE", "CORE"]:
        rows.append({"group": k, "face €bn": mtm[k]["face"] / 1e3,
                     "MV/face": mtm[k]["mv0"] / mtm[k]["face"],
                     "Macaulay y": reps[k]["macaulay_y"],
                     "at Dec-11": mtm[k]["y2011"] * 100,
                     "at Feb-12": mtm[k]["y2012"] * 100,
                     "at PSI NPV": mtm[k]["psi"] * 100})
    print(pd.DataFrame(rows).set_index("group").round(2).to_string())
    print(f"\n  Greek book lost {reps['GR']['loss_frac_mv']*100:.1f}% of market "
          f"value over 2011 on the flat-curve marking.")
    phi_D_mv = mtm["GR"]["mv0"] / ct1_of["GR"]
    moments["stocks"]["phi_bD_D_market_value_over_CT1"] = phi_D_mv
    print(f"  Model translation: sovereign book at MARKET value is {phi_D_mv:.2f}x "
          f"Greek bank capital ({phi_D:.2f}x at face), so a fall of x% in q_b_D "
          f"costs {phi_D_mv:.2f}x% of capital.")
    print(f"  The observed {reps['GR']['loss_frac_mv']*100:.1f}% price fall is "
          f"therefore the impulse the model has to generate.")

    print("\n--- is the model's risk impulse of that scale? ---")
    imp = {f"delta_b={db}": impulse_check(reps["GR"]["loss_frac_mv"], db)
           for db in (0.056, 0.0777)}
    for k, v in imp.items():
        print(f"  {k:14s} pi 0.10% -> 1.98% moves q_b_D {v['actuarial_price_fall']*100:5.1f}%"
              f"   (actuarial only; observed {v['observed_price_fall']*100:.1f}%)"
              f"   pi matching observed = {v['pi_matching_observed']*100:.2f}%/qtr")
    moments["shock"]["impulse_check"] = imp
    with open(os.path.join(OUT, "eba_discipline.json"), "w") as f:
        json.dump(moments, f, indent=2)

    print("\n=== 3. TRANSMISSION: corporate borrowing spreads ===")
    r = lending_rates()
    for suf, lab in [("cob", "composite cost of borrowing"),
                     ("ovd", "overdrafts / revolving")]:
        pre = r.loc["2005":"2007"]
        cri = r.loc["2010":"2012-06"]
        gr_pre = (pre[f"GR_{suf}"] - pre.eur3m).mean()
        de_pre = (pre[f"DE_{suf}"] - pre.eur3m).mean()
        gr_cri = (cri[f"GR_{suf}"] - cri.eur3m).mean()
        de_cri = (cri[f"DE_{suf}"] - cri.eur3m).mean()
        print(f"  {lab:32s}  GR {gr_pre:5.2f} -> {gr_cri:5.2f}   "
              f"DE {de_pre:5.2f} -> {de_cri:5.2f}   "
              f"diff-in-diff {(gr_cri-gr_pre)-(de_cri-de_pre):+.2f} pp")

    print("\n=== figures ===")
    emitted = {"fig_eba_ownership", "fig_eba_ladder", "fig_eba_structure",
               "fig_eba_mtm", "fig_lending_spreads"}
    assert emitted == set(CAPTIONS), (emitted ^ set(CAPTIONS))
    fig_ownership(p, grp)
    fig_ladder(ladders, reps)
    fig_structure(p, grp, ladders, reps)
    fig_mtm(mtm)
    fig_spreads(r)
    return moments


if __name__ == "__main__":
    main()
