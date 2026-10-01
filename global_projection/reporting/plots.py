# FIGURES FOR THE RECURSIVE PROJECTION EXPERIMENTS.
# Every figure is written to OUTDIR -- results/GLOBAL/figures by default; run.py points
# it at its own results directory before drawing. The projection experiments print their
# own IRF tables; these are the figures they produce.
import os

import numpy as np
import matplotlib.pyplot as plt

# repo root is two levels up from reporting/
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                      "results", "GLOBAL", "figures")

# CATEGORICAL palette for the output-decomposition channels (Okabe-Ito derived;
# passes the lightness/chroma/CVD-separation/contrast checks in this fixed order --
# the order is load-bearing, hues are assigned by channel and never cycled).
# The residual is a deliberate NEUTRAL: it is an accuracy diagnostic, not a channel.
CHANNEL_COLORS = {"credit_spread": "#D55E00", "deposit_rate": "#0072B2",
                  "capital": "#009E73", "rel_price": "#7570B3",
                  "tfp": "#A6761D", "residual": "#9E9E9E"}
# SEQUENTIAL ramps (ORDERED variables, so one hue light->dark, not categorical hues):
# policy scenarios in increasing strength, and the income quintiles.
ACTIVATION_RAMP = ("#6BAED6", "#2171B5", "#08306B")
QUINTILE_RAMP = ("#C7E0B4", "#8FC98A", "#4DA65B", "#1F7A3D", "#0B4526")
SURFACE = "#FFFFFF"
INK, INK_MUTED = "#1A1A1A", "#5C5C5C"


def _save(fig, filename):
    # LAY OUT AND WRITE THE FIGURE TO OUTDIR.
    # tight_layout is skipped when the figure already has a layout engine. Figures
    # with a secondary_yaxis MUST use the constrained engine: tight_layout does not
    # see secondary axes at all, so it packs the panels as if they were absent and the
    # right-hand annualised labels land on top of the next panel's y-label.
    if fig.get_layout_engine() is None:
        fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _style(ax, title, ylabel, xlabel="quarter", zero=True):
    # RECESSIVE AXES: horizontal grid behind the marks, no top/right spines.
    ax.set_title(title, fontsize=10, color=INK)
    ax.set_ylabel(ylabel, fontsize=8, color=INK_MUTED)
    ax.set_xlabel(xlabel, fontsize=8, color=INK_MUTED)
    ax.tick_params(labelsize=8, colors=INK_MUTED)
    ax.grid(axis="y", color="#E6E6E6", lw=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#CCCCCC")
    if zero:
        # only for DEVIATION panels: on a level series (a bond price near 0.8, a spread
        # near 300bp) forcing zero into view squashes the variation being read
        ax.axhline(0, color=INK, lw=0.8)


def _stacked_channels(ax, dec, chans, title, ylab=None, ann=True,
                      total_label="Total output response"):
    # SIGNED STACKED BARS: POSITIVE AND NEGATIVE CONTRIBUTIONS STACK SEPARATELY, SO
    # THE VISIBLE TOP/BOTTOM OF THE STACK IS THE NET RESPONSE (matplotlib's
    # stackplot cannot do this -- it assumes one sign).
    n = len(dec["total"])
    t = np.arange(n)
    pos = np.zeros(n)
    neg = np.zeros(n)
    for key, label in chans:
        v = dec[key]
        base = np.where(v >= 0, pos, neg)
        ax.bar(t, v, bottom=base, width=0.82, color=CHANNEL_COLORS[key],
               label=label, edgecolor=SURFACE, linewidth=0.9)   # surface gap
        pos = pos + np.maximum(v, 0.0)
        neg = neg + np.minimum(v, 0.0)
    ax.plot(t, dec["total"], color=INK, lw=2.0, marker="o", ms=3.2,
            label=total_label, zorder=5)
    _style(ax, title, ylab or "% deviation from the no-shock path (level)")
    if ann:
        # ONLY for a quarterly FLOW. A bond PRICE gap has no annualised reading, so the
        # bond decomposition passes ann=False rather than inviting the ratio to be read.
        _annual_axis(ax, "annualised, % (Bocola Table 5 unit)")


def plot_output_decomposition(cases, chans, note="",
                              filename="output_decomposition.png"):
    # WHICH FACTORS PRODUCE THE OUTPUT RESPONSE, AND WHAT THE BACKSTOP CHANGES.
    # cases = list of (label, decomposition dict); the first is the no-backstop
    # reference, the last the strongest backstop. Panel 3 is the DIFFERENCE
    # between them, channel by channel -- the transmission OMT actually operates on.
    # A SINGLE case draws the one panel: main.py reads this decomposition off the
    # baseline risk solve, where there is no backstop to difference against, and a
    # three-panel layout with two of them duplicated would misrepresent that.
    if len(cases) == 1:
        # ONE title only, as in plot_bond_decomposition: a suptitle, a fig.text note and
        # a panel title need three panels' worth of width, and on a single panel the
        # three lines land on top of each other.
        lab, dec = cases[0]
        fig, ax = plt.subplots(1, 1, figsize=(8.2, 5.1), layout="constrained")
        fig.suptitle(f"What produces the output response to sovereign risk — {lab}"
                     + (f"\n{note}" if note else ""), fontsize=12, color=INK)
        _stacked_channels(ax, dec, chans, "")
        ax.legend(fontsize=7.5, frameon=False, loc="best")
        _save(fig, filename)
        return os.path.join(OUTDIR, filename)
    lo_lab, lo = cases[0]
    hi_lab, hi = cases[-1]
    n = len(lo["total"])
    t = np.arange(n)
    fig, axes = plt.subplots(1, 3, figsize=(16.8, 5.1), layout="constrained")
    fig.suptitle("What produces the output response to sovereign risk, and what the "
                 "OMT/TPI backstop changes", fontsize=12, y=1.04, color=INK)
    if note:
        fig.text(0.5, 0.985, note, ha="center", fontsize=8, color=INK_MUTED)

    _stacked_channels(axes[0], lo, chans, f"Output decomposition — {lo_lab}")
    _stacked_channels(axes[1], hi, chans, f"Output decomposition — {hi_lab}")
    for key, label in chans:
        axes[2].plot(t, hi[key] - lo[key], color=CHANNEL_COLORS[key], lw=2.0,
                     label=label)
    axes[2].plot(t, hi["total"] - lo["total"], color=INK, lw=2.0, marker="o",
                 ms=3.2, label="Total output response")
    _style(axes[2], f"Backstop effect by channel ({hi_lab} minus {lo_lab})",
           "pp of output, difference")
    ymin = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
    ymax = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
    axes[0].set_ylim(ymin, ymax)                     # one shared scale, never two
    axes[1].set_ylim(ymin, ymax)
    axes[0].legend(fontsize=7.5, frameon=False, loc="best")
    axes[2].legend(fontsize=7.5, frameon=False, loc="best")
    _save(fig, filename)
    return os.path.join(OUTDIR, filename)


def plot_welfare_quintiles(labels, cost, gain, cons, quintile_income, note="",
                           filename="omt_welfare_quintiles.png"):
    # OMT/TPI WELFARE INCIDENCE ACROSS THE INCOME DISTRIBUTION.
    # cost[i, q] = CEV cost of the risk shock under scenario i, quintile q (%);
    # gain[i, q] = cost[i] - cost[0], the backstop's welfare improvement (pp);
    # cons[q, t] = per-quintile consumption path deviation under NO backstop (%).
    n_q = cost.shape[1]
    xs = np.arange(n_q)
    names = [f"Q{k + 1}" for k in range(n_q)]
    ramp = [ACTIVATION_RAMP[min(i, len(ACTIVATION_RAMP) - 1)]
            for i in range(len(labels))]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.9))
    fig.suptitle("Welfare gain from the OMT/TPI backstop by income quintile "
                 "(consumption-equivalent, incomplete-markets overlay)",
                 fontsize=12, y=1.04, color=INK)
    if note:
        fig.text(0.5, 0.985, note, ha="center", fontsize=8, color=INK_MUTED)

    w = 0.8 / len(labels)
    for i, lab in enumerate(labels):
        axes[0].bar(xs + (i - (len(labels) - 1) / 2) * w, cost[i], width=w * 0.9,
                    color=ramp[i], label=lab, edgecolor=SURFACE, linewidth=0.9)
    _style(axes[0], "Welfare effect of the sovereign-risk shock",
           "% permanent consumption (negative = cost)", "income quintile")
    axes[0].set_xticks(xs, names)
    axes[0].legend(fontsize=7.5, frameon=False)

    gl = labels[1:]
    # the gains can be small in absolute terms; label them at a precision that
    # actually resolves them rather than printing a column of zeros
    digits = max(2, min(6, int(np.ceil(-np.log10(max(np.max(np.abs(gain)), 1e-9)))) + 2))
    for i, lab in enumerate(gl):
        off = (i - (len(gl) - 1) / 2) * (0.8 / max(len(gl), 1))
        bars = axes[1].bar(xs + off, gain[i + 1], width=0.8 / max(len(gl), 1) * 0.9,
                           color=ramp[i + 1], label=lab, edgecolor=SURFACE,
                           linewidth=0.5)
        axes[1].bar_label(bars, fmt=f"%.{digits}f", fontsize=6.5, padding=1.5,
                          color=INK_MUTED)
    _style(axes[1], "Welfare improvement from the backstop",
           "pp of permanent consumption", "income quintile")
    axes[1].set_xticks(xs, names)
    axes[1].legend(fontsize=7.5, frameon=False)

    t = np.arange(cons.shape[1])
    for q in range(n_q):
        axes[2].plot(t, cons[q], color=QUINTILE_RAMP[q % len(QUINTILE_RAMP)], lw=2.0)
        axes[2].annotate(names[q], (t[-1], cons[q, -1]), fontsize=7.5,
                         color=QUINTILE_RAMP[q % len(QUINTILE_RAMP)],
                         xytext=(3, 0), textcoords="offset points", va="center")
    _style(axes[2], "Consumption incidence of the shock, no backstop", "% deviation")
    axes[2].set_xlim(0, t[-1] * 1.06)

    sub = "  ".join(f"{names[q]} inc {quintile_income[q]:.3f}" for q in range(n_q))
    axes[0].text(0.0, -0.24, sub, transform=axes[0].transAxes, fontsize=6.5,
                 color=INK_MUTED)
    _save(fig, filename)


# SERIES COLOURS FOR THE IRF PANELS, taken from the SAME Okabe-Ito set the
# decomposition uses so the figures read as one system: the shock and the price it
# moves are neutral ink, the credit spread keeps its channel hue, and the real
# quantities keep theirs. Never introduce a hue that is not already in the palette.
# BOND-PRICE LEGS. Two of them are the SAME OBJECT as a leg of the output
# decomposition and keep its hue: the deposit rate (#0072B2) and the bank constraint,
# which is lambda*mu on both figures and so takes credit_spread's vermilion. The other
# three get hues unused elsewhere here, all CVD-safe, and the residual stays the
# documented neutral grey -- it is an accuracy diagnostic, not a channel.
CHANNEL_COLORS["liquidity_premium"] = CHANNEL_COLORS["credit_spread"]   # lambda*mu
CHANNEL_COLORS["continuation"] = "#E69F00"      # amber -- the dominant leg
CHANNEL_COLORS["expected_loss"] = "#882255"     # dark magenta
CHANNEL_COLORS["risk_premium"] = "#44AA99"      # teal

# ANNUALISATION FACTOR FOR A QUARTERLY-FLOW LEVEL GAP. This is Bocola's Table 5 unit:
# his output losses are cumsum(g_s - g_ns)*400 where g is a quarterly log growth rate,
# so the cumulated object is the log LEVEL gap and the 400 is 100 (to %) x 4 (to an
# annual rate). Reporting both on one axis is what makes his -1.05/-1.44/-1.53 and this
# model's level IRFs readable against each other without a conversion in the reader's head.
ANN = 4.0


def _annual_axis(ax, label="ann. %"):
    # RIGHT-HAND TWIN SHOWING THE SAME SERIES AT AN ANNUAL RATE (x4).
    # A secondary_yaxis, not a twinx: it is a relabelling of the SAME data, so it must
    # not be able to drift out of registration with the left axis.
    sec = ax.secondary_yaxis("right", functions=(lambda v: ANN * v, lambda v: v / ANN))
    sec.set_ylabel(label, fontsize=7, color=INK_MUTED, labelpad=1)
    sec.tick_params(labelsize=7, colors=INK_MUTED, pad=1)
    sec.spines["right"].set_color("#CCCCCC")
    sec.yaxis.set_major_formatter(lambda v, _: f"{v:+.2f}")
    return sec


def plot_bond_decomposition(dec, chans, note="",
                            filename="bond_decomposition.png"):
    # WHY THE D SOVEREIGN REPRICES: the bank's own FOC, split leg by leg.
    # Same stacked-signed-bar treatment as the output decomposition, because it is the
    # same kind of object -- an identity, not an attribution. Bocola's Table 4 splits
    # the EXCESS RETURN into a risk premium and a liquidity premium; this splits the
    # PRICE, and adds the two legs his table takes as given (the discount rate and the
    # continuation price), so the bars sum to the observed repricing.
    # ONE title only. The suptitle-plus-note-plus-panel-title stack the output
    # decomposition uses needs three panels' worth of width; on a single panel the
    # three lines land on top of each other.
    fig, ax = plt.subplots(1, 1, figsize=(8.2, 5.1), layout="constrained")
    fig.suptitle("What reprices the D sovereign  $Q_{b,D}$"
                 + (f"\n{note}" if note else ""), fontsize=12, color=INK)
    _stacked_channels(ax, dec, chans, "", ann=False,
                      ylab="% deviation from the no-shock path",
                      total_label="Total bond-price response")
    ax.legend(fontsize=7.5, frameon=False, loc="best")
    _save(fig, filename)
    return os.path.join(OUTDIR, filename)


# THE PAPER FIGURE SCHEME (risk + TFP IRFs). Serif type, BOLD LETTERED panel titles,
# no gridlines, only the left and bottom rules, the legend inside the first panel and NO
# figure title -- the paper's caption carries it. Every series is shown at its ANNUALISED
# reading, so there is no secondary axis on these figures.
PAPER_RC = {"font.family": "serif",
            "font.serif": ["Palatino", "Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 0.8}
# DOMESTIC IS STANFORD RED, FOREIGN IS OXFORD BLUE. They differ in dash pattern too, so
# the pair survives greyscale printing and colour-vision deficiency.
COUNTRY_STYLE = (("#8C1515", "-", "Domestic"), ("#002147", "--", "Foreign"))

# (D key, F key, title, y label, annualise). ANNUALISE ONLY THE QUARTERLY FLOWS: x4 on a
# level gap is Bocola's Table 5 unit, and it is what an annual rate MEANS. The spread is
# already an annualised rate, and a stock (net worth) or a price (the bond) has no annual
# reading at all -- multiplying those by four would invent one.
PAPER_PANELS = (
    ("Y", "Y_F", "GDP", "annualised % deviation", True),
    ("spread", "spread_F", "Lending spread", "basis points per year", False),
    ("n", "n_F", "Bank net worth", "% deviation", False),
    ("Q_bD", "Q_bF", "Government bond price", "level", False),
    ("C", "C_F", "Consumption", "annualised % deviation", True),
    ("I", "I_F", "Investment", "annualised % deviation", True),
)


def _paper_axes(ax, title, ylabel, zero=True):
    # ONE PANEL OF THE PAPER SCHEME: bold lettered title, bare left/bottom rules.
    # The bold is a request, not a guarantee: the system Palatino ships one weight, so
    # the title renders regular here and bold wherever a bold serif face is installed.
    ax.set_title(title, fontsize=11, color=INK, fontweight="bold", pad=8)
    ax.set_ylabel(ylabel, fontsize=9, color=INK)
    ax.set_xlabel("quarter", fontsize=9, color=INK)
    ax.tick_params(labelsize=8.5, colors=INK, direction="out", length=3)
    ax.grid(False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK)
    if zero:
        # only on a DEVIATION panel: on a level series (a bond price near 0.8, a spread
        # near 100 bp) forcing zero into view squashes the variation being read
        ax.axhline(0.0, color="#B0B0B0", lw=0.7, zorder=1)


def _paper_irf(path, filename):
    # THE FIVE-PANEL COUNTRY-PAIR IRF FIGURE, shared by the risk and TFP experiments.
    # Both experiments record the same five series per country, so one panel spec draws
    # both figures and the two are read against each other panel by panel.
    q = np.arange(len(path["Y"]))
    with plt.rc_context(PAPER_RC):
        fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.0), layout="constrained")
        flat = axes.ravel()
        for i, (kd, kf, title, ylab, ann) in enumerate(PAPER_PANELS):
            ax = flat[i]
            scale = ANN if ann else 1.0
            for key, (colour, ls, lab) in zip((kd, kf), COUNTRY_STYLE):
                if key not in path:
                    continue
                ax.plot(q, scale * np.asarray(path[key], dtype=float), color=colour,
                        ls=ls, lw=1.5, label=lab, zorder=3)
            _paper_axes(ax, f"({chr(97 + i)}) {title}", ylab,
                        zero="deviation" in ylab)
            ax.set_xlim(q[0], q[-1])
        flat[0].legend(fontsize=9, frameon=False, loc="best")
        for ax in flat[len(PAPER_PANELS):]:
            ax.set_visible(False)
        _save(fig, filename)
    return os.path.join(OUTDIR, filename)


def plot_risk_irf(path, filename="risk_irf_recursive.png", note=""):
    # SOVEREIGN-RISK IRF, BOTH COUNTRIES, PAPER SCHEME.
    # note is accepted and NOT drawn: these figures carry no title by design.
    return _paper_irf(path, filename)


def plot_tfp_irf(path, filename="tfp_irf_recursive.png", note=""):
    # TFP IRF ALONG THE Z-DECAY PATH (the no-default rules), BOTH COUNTRIES.
    return _paper_irf(path, filename)
