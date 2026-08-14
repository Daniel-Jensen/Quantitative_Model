# FIGURES FOR THE RECURSIVE PROJECTION EXPERIMENTS.
# Every figure is written to output/. The projection experiments print their
# own IRF tables; these are the figures they produce.
import os

import numpy as np
import matplotlib.pyplot as plt

# output/ lives at the package root (one level up from reporting/)
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")

# CATEGORICAL palette for the output-decomposition channels (Okabe-Ito derived;
# passes the lightness/chroma/CVD-separation/contrast checks in this fixed order --
# the order is load-bearing, hues are assigned by channel and never cycled).
# The residual is a deliberate NEUTRAL: it is an accuracy diagnostic, not a channel.
CHANNEL_COLORS = {"credit_spread": "#D55E00", "deposit_rate": "#0072B2",
                  "capital": "#009E73", "rel_price": "#7570B3",
                  "tfp": "#A6761D", "residual": "#9E9E9E"}
# SEQUENTIAL ramp for backstop strength (an ORDERED variable, so one hue
# light->dark, not categorical hues) and for the income quintiles.
ACTIVATION_RAMP = ("#6BAED6", "#2171B5", "#08306B")
QUINTILE_RAMP = ("#C7E0B4", "#8FC98A", "#4DA65B", "#1F7A3D", "#0B4526")
SURFACE = "#FFFFFF"
INK, INK_MUTED = "#1A1A1A", "#5C5C5C"


def _save(fig, filename):
    # TIGHT-LAYOUT AND WRITE THE FIGURE TO output/.
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _style(ax, title, ylabel, xlabel="quarter"):
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
    ax.axhline(0, color=INK, lw=0.8)


def plot_activation_irf(scenarios, filename="tpi_activation_recursive.png"):
    # OVERLAY THE PROJECTION-SOLVER IRFs UNDER OMT/TPI ACTIVATION SCENARIOS.
    # scenarios = list of (label, paths, color); paths = dict of pre-computed %/bp
    # series (Y_D, C_D, I_D, Q_bD, spread, pd) over the shock-decay horizon.
    n = len(scenarios[0][1]["Y_D"])
    t = np.arange(n)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle("Sovereign-risk shock under OMT/TPI activation (recursive "
                 "Chebyshev projection, three-branch quadrature)", fontsize=12, y=1.01)

    def panel(ax, key, title, ylabel):
        for label, paths, color in scenarios:
            ax.plot(t, paths[key], color=color, lw=1.8, label=label)
        ax.axhline(0, color="k", lw=0.7, ls=":")
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=8); ax.set_xlabel("quarter", fontsize=8)
        ax.legend(fontsize=8)

    panel(axes[0, 0], "Y_D", "Output", "% dev.")
    panel(axes[0, 1], "C_D", "Consumption", "% dev.")
    panel(axes[0, 2], "I_D", "Investment", "% dev.")
    panel(axes[1, 0], "Q_bD", "Sovereign bond price Q_bD", "% dev.")
    panel(axes[1, 1], "spread", "Lending spread", "bps ann. dev.")
    panel(axes[1, 2], "pd", "Priced default probability", "% per quarter")
    _save(fig, filename)


def _stacked_channels(ax, dec, chans, title):
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
            label="Total output response", zorder=5)
    _style(ax, title, "% deviation from the no-shock path")


def plot_output_decomposition(cases, chans, note="",
                              filename="output_decomposition.png"):
    # WHICH FACTORS PRODUCE THE OUTPUT RESPONSE, AND WHAT THE BACKSTOP CHANGES.
    # cases = list of (label, decomposition dict); the first is the no-backstop
    # reference, the last the strongest backstop. Panel 3 is the DIFFERENCE
    # between them, channel by channel -- the transmission OMT actually operates on.
    lo_lab, lo = cases[0]
    hi_lab, hi = cases[-1]
    n = len(lo["total"])
    t = np.arange(n)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.9))
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
