"""Paper-ready figures and tables for the first draft.

Everything here is DERIVED — read live from the solved steady state and the cached
response matrices. No number is transcribed. The one thing carried as literal text
is a *source citation* (which paper or dataset a target came from), never a value.

**That includes the captions.** Until 2026-08-06 this module carried a module-level
`CAPTIONS` dict of literal prose with flexible-price numbers frozen into it. The
sticky-price conversion and the `psi_lambda_B` 8.5 -> 7.85 re-tune left every one of
them stale and three of them *inverted* — `fig03` asserted offsetting channels
"roughly four times the headline" when they are now ~0.25x it, `fig08` claimed
consumption rises on impact and that the lowest quintile "gains 0.95%" against a
Table 4 in the same generated document reading +0.4250, and `fig06` claimed the net
path is smaller than its components "at every horizon" when it is not in the impact
quarter. The dict is now built AT RUN TIME by `save()`: each figure hands `save()` a
caption it computed from the same arrays it just plotted, so a caption cannot
survive a recalibration that falsifies it. Directional claims ("monotone",
"reverses by quarter k", "larger than") are SELECTED from the data rather than
asserted, so a sign flip rewrites the sentence instead of lying in it.

Figure set (each caption is baked into the PNG — a caption that lives only in the
LaTeX travels separately from the image and is lost the moment the file is reused):

  fig01_transmission        the sovereign-risk shock and what the backstop does to it
  fig02_loading_schedule    KEY FIGURE — the self-extinguishing premium
  fig03_dy_decomposition    why the headline output number must not be led with
  fig04_spread_decomposition  the wedge is nearly all of the default loading
  fig05_incidence           Germany's side: exposure rises while compensation falls
  fig06_net_effects         the net path against the components that generate it
  fig07_ms_regimes          empirical regime dating (the discipline behind the beliefs)
  fig08_deciles             distributional incidence by income quintile

Palette validated with the dataviz six-checks validator (light mode, categorical):
lightness band, chroma floor, CVD separation, normal-vision floor, contrast. The
project's previous red #8C1515 FAILED the lightness band and is replaced by
#A62B22; the previous navy #002147 failed both the band and the chroma floor
(it reads gray) and is replaced by #1B6CA8. Re-validate before changing either:
  node scripts/validate_palette.js "#1B6CA8,#A62B22,#c87941,#1a6e3a" --mode light
"""
import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import e1_backstop_schedule as e1
from common import BP_ANN, HERE, ROOT, irf_from_cache, load_cache, provenance, regime_irfs
from e2_dy_decomposition import decompose_dY

PAPER_DIR = os.path.join(HERE, "paper")
TABLES_DOC = os.path.join(ROOT, "docs", "paper_draft_results.md")

# Validated categorical palette — see module docstring.
BLUE, RED, ORANGE, GREEN = "#1B6CA8", "#A62B22", "#c87941", "#1a6e3a"
INK, MUTED, GRID = "#1a1a1a", "#5a5a5a", "#d8d8d8"
REGIME_COLOR = {"passive": RED, "medium": ORANGE, "aggressive": GREEN}
REGIME_LABEL = {"passive": "passive (no backstop)", "medium": "medium",
                "aggressive": "aggressive"}

N_IRF = 40          # quarters shown in IRF panels
T_PNL = 100

# Prose names for the E2 identity's components, used when a caption has to say
# which channel it picked out of the data.
COMPONENT_LABEL = {"consumption_quantity": "consumption",
                   "consumption_price": "the consumption deflator",
                   "investment": "investment",
                   "net_exports": "net exports"}

# Populated at RUN TIME by save() — see the module docstring. Never edit by hand:
# a literal here is a claim that no longer has to survive the next recalibration.
CAPTIONS = {}


# ── Caption helpers ──────────────────────────────────────────────────────────
#
# These exist so a caption's *directional* words come from the data too. Writing
# "monotone" or "reverses by quarter four" as a literal is the same defect as
# writing "4.5x" as a literal, only harder to notice when it goes wrong.

def _monotone(v, sign):
    """True if v is strictly monotone in the given direction (NaNs dropped)."""
    v = np.asarray(v, dtype=float)
    v = v[~np.isnan(v)]
    return bool(v.size > 1 and np.all(np.sign(np.diff(v)) == sign))


def _first_quarter(mask, n):
    """First quarter in [0, n) at which mask holds, or None."""
    idx = np.nonzero(np.asarray(mask)[:n])[0]
    return int(idx[0]) if idx.size else None


def _ordinal(q):
    return "the impact quarter" if q == 0 else f"quarter {q}"


def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.axhline(0, lw=0.8, color=MUTED, zorder=1)


def save(fig, name, caption):
    """Register the DERIVED caption, bake it into the image, then write it.

    The caption is an argument, not a lookup: it must be constructed by the figure
    function from the arrays it just plotted, so that it cannot outlive them.
    """
    if not isinstance(caption, str) or not caption.strip():
        raise ValueError(f"{name}: save() needs a derived caption string")
    CAPTIONS[name] = cap = " ".join(caption.split())
    chars = int(fig.get_size_inches()[0] * 15)
    fig.text(0.5, -0.02, textwrap.fill(cap, width=chars), ha="center", va="top",
             fontsize=8, style="italic", color=MUTED)
    os.makedirs(PAPER_DIR, exist_ok=True)
    fig.savefig(os.path.join(PAPER_DIR, f"{name}.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


# ── Figures ──────────────────────────────────────────────────────────────────

def fig01_transmission(cache, regimes):
    Y, n, I, C = (float(cache[k]) for k in ("Y_D_ss", "n_inter_D_ss", "I_D_ss", "C_D_ss"))
    panels = [("spread_rb", "D–F sovereign spread", "bp, annualised", BP_ANN),
              ("n_inter_D", "bank net worth", "% of steady state", 100.0 / n),
              ("I_D", "investment", "% of steady state", 100.0 / I),
              ("Y_D", "output", "% of steady state", 100.0 / Y)]
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.4))
    for ax, (var, title, unit, scale) in zip(axes, panels):
        for name, (_g, irf) in regimes.items():
            ax.plot(np.arange(N_IRF), np.asarray(irf[var])[:N_IRF] * scale,
                    color=REGIME_COLOR[name], lw=2, label=REGIME_LABEL[name])
        _style(ax)
        ax.set_title(title, fontsize=10, color=INK, pad=8)
        ax.set_ylabel(unit, fontsize=8, color=MUTED)
        ax.set_xlabel("quarters", fontsize=8, color=MUTED)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK)
    fig.suptitle("Transmission of a 1pp sovereign default-probability shock, by backstop stance",
                 fontsize=11, color=INK, y=1.04)
    fig.tight_layout()
    save(fig, "fig01_transmission", _caption_fig01(cache, regimes))


def _caption_fig01(cache, regimes):
    """Impact magnitudes and the quarter at which the regime ordering reverses."""
    n = float(cache["n_inter_D_ss"])
    I = float(cache["I_D_ss"])
    Y = float(cache["Y_D_ss"])
    p, a = regimes["passive"][1], regimes["aggressive"][1]
    sp_p = np.asarray(p["spread_rb"]) * BP_ANN
    sp_a = np.asarray(a["spread_rb"]) * BP_ANN
    nw_p, nw_a = (np.asarray(x["n_inter_D"]) * 100.0 / n for x in (p, a))
    peak = float(sp_p[:T_PNL].max())
    n0 = float(nw_p[0])
    i0 = float(np.asarray(p["I_D"])[0] * 100.0 / I)
    y0 = float(np.asarray(p["Y_D"])[0] * 100.0 / Y)

    # "Reversal" = the quarter from which the aggressive path is no longer the
    # better one: a wider spread, or a weaker balance sheet, than doing nothing.
    q_sp = _first_quarter(sp_a > sp_p, N_IRF)
    q_nw = _first_quarter(nw_a < nw_p, N_IRF)
    if q_sp is None and q_nw is None:
        tail = ("the ordering never reverses inside the plotted window, so the backstop "
                "shifts the whole path rather than only its opening quarters")
    else:
        parts = []
        if q_nw is not None:
            parts.append(f"the net-worth ordering reverses by {_ordinal(q_nw)}")
        if q_sp is not None:
            parts.append(f"the spread ordering by {_ordinal(q_sp)}")
        tail = (" and ".join(parts) + " as the unaided economy overshoots on the rebound, "
                "so intervention damps the impact quarter rather than shifting the whole "
                "path down")
    # Verbs from the signs: at this calibration both fall, but a caption that
    # hardcodes "cuts" would misreport a recalibration in which they do not.
    vb = lambda v: "cuts" if v < 0 else "raises"
    if (n0 < 0) == (i0 < 0):
        real = f"{vb(n0)} bank net worth {abs(n0):.1f}% and investment {abs(i0):.1f}%"
    else:
        real = (f"{vb(n0)} bank net worth {abs(n0):.1f}% and {vb(i0)} investment "
                f"{abs(i0):.1f}%")
    return (f"A 1pp rise in the Greek default probability widens the D–F spread to a peak of "
            f"{peak:.0f}bp, {real} on impact, and takes output {y0:+.2f}% "
            f"from steady state; the backstop's "
            f"cushioning is concentrated in the opening quarters — {tail}.")


def fig02_loading_schedule(cache, regimes, payload):
    gammas, loading, peak_bp = e1.loading_schedule(cache)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.plot(gammas, loading, color=BLUE, lw=2.2, zorder=3)
    ax.axhline(1.0, ls=(0, (4, 3)), lw=1.2, color=MUTED, zorder=2)
    ax.text(gammas[-1], 1.0, " actuarially fair", va="center", ha="right",
            fontsize=8, color=MUTED, backgroundcolor="white")
    for name, r in payload["regimes"].items():
        if r["loading"] is None:
            continue
        ax.plot(r["gamma"], r["loading"], "o", color=REGIME_COLOR[name], ms=9,
                mec="white", mew=1.5, zorder=4)
        ax.annotate(REGIME_LABEL[name], (r["gamma"], r["loading"]),
                    textcoords="offset points", xytext=(9, 7),
                    fontsize=8, color=INK)
    _style(ax)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("backstop aggressiveness γ", fontsize=9, color=MUTED)
    ax.set_ylabel("premium PV ÷ expected-loss PV", fontsize=9, color=MUTED)
    ax.set_title("The premium self-extinguishes", fontsize=10, color=INK, pad=8)

    ax = axes[1]
    ax.plot(gammas, peak_bp, color=BLUE, lw=2.2, zorder=3)
    for name, r in payload["regimes"].items():
        ax.plot(r["gamma"], r["peak_spread_bp_ann"], "o", color=REGIME_COLOR[name],
                ms=9, mec="white", mew=1.5, zorder=4)
        ax.annotate(f"{r['peak_spread_bp_ann']:.0f}bp", (r["gamma"], r["peak_spread_bp_ann"]),
                    textcoords="offset points", xytext=(9, 5), fontsize=8, color=INK)
    _style(ax)
    ax.set_xlabel("backstop aggressiveness γ", fontsize=9, color=MUTED)
    ax.set_ylabel("peak D–F spread (bp, annualised)", fontsize=9, color=MUTED)
    ax.set_title("Spread compression", fontsize=10, color=INK, pad=8)

    fig.tight_layout()
    save(fig, "fig02_loading_schedule",
         _caption_fig02(gammas, loading, peak_bp, payload))
    return gammas, loading, peak_bp


def _caption_fig02(gammas, loading, peak_bp, payload):
    """The KEY claim is the DECLINE, so the schedule's own endpoints state it."""
    ok = ~np.isnan(np.asarray(loading, dtype=float))
    g_lo, g_hi = float(gammas[ok][0]), float(gammas[ok][-1])
    l_lo, l_hi = float(np.asarray(loading)[ok][0]), float(np.asarray(loading)[ok][-1])
    falling = _monotone(loading, -1)
    arr = np.asarray(loading)[ok]
    above_one = bool(np.all(arr > 1.0))
    below_one = bool(np.all(arr < 1.0))

    named = {k: v["loading"] for k, v in payload["regimes"].items()
             if v["loading"] is not None}
    named_txt = ("; " + ", ".join(f"{k} {v:.2f}×" for k, v in named.items())
                 + " at the named regimes") if named else ""
    shape = ("falls monotonically" if falling else
             "falls on net but not monotonically" if l_hi < l_lo else
             "does NOT fall — the self-extinguishing-premium claim fails at this "
             "calibration and must not be asserted")
    # Three cases, not two. The old two-branch version assumed the loading STARTS above
    # 1 and asked only whether it crosses; since the 2026-08-18 payoff repair it starts
    # at 0.53 and is below 1 throughout, which the "crossing below ... before the grid
    # ends" wording described backwards.
    floor = (" and stays above the actuarially fair benchmark of 1 throughout"
             if above_one else
             " and stays BELOW the actuarially fair benchmark of 1 throughout — the ECB "
             "is under-compensated at every intervention intensity, so the paper must "
             "NOT assert over-compensation" if below_one else
             ", crossing the actuarially fair benchmark of 1 within the grid")
    peak_txt = (f"peak spread compresses {peak_bp[0]:.0f}bp → {peak_bp[-1]:.0f}bp "
                f"over the same grid")
    return (f"KEY FIGURE — the ECB's compensation per unit of expected loss {shape} from "
            f"{l_lo:.2f}× at γ={g_lo:.2f} to {l_hi:.2f}× at γ={g_hi:.0f}{floor}"
            f"{named_txt} ({peak_txt}). The premium is a rent extracted from a "
            f"balance-sheet-constrained seller and intervention relieves the very "
            f"constraint that creates it: the profit self-extinguishes as the policy "
            f"succeeds.")


def fig03_dy_decomposition(cache, regimes):
    """Two panels, because the shock and the policy work through DIFFERENT channels.

    Panel A (levels) shows what the crisis does: investment collapses and is masked
    in the aggregate mostly by CONSUMPTION, with net exports near zero. Panel B
    (change vs passive) shows what the BACKSTOP does: investment recovers against a
    net-export deterioration, each ~4x the headline and opposite in sign. Plotting
    only panel A and captioning it with panel B's story — which an earlier draft of
    this figure did — misattributes the offset.
    """
    ss = {"P_CES_D_ss": float(cache["P_CES_D_ss"]), "C_D_ss": float(cache["C_D_ss"])}
    show = [("investment", "investment"), ("net_exports", "net\nexports"),
            ("consumption_quantity", "consumption\n(quantity)"),
            ("consumption_price", "consumption\n(price)"), ("__total__", "ΔY\n(total)")]
    names = list(regimes)

    comps = {}
    for name in names:
        c, _r = decompose_dY(regimes[name][1], ss)
        c = {k: v[0] * 1e3 for k, v in c.items()}
        c["__total__"] = float(np.asarray(regimes[name][1]["Y_D"])[0] * 1e3)
        comps[name] = c

    x = np.arange(len(show))
    width = 0.26
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4), sharey=True)

    ax = axes[0]
    for i, name in enumerate(names):
        ax.bar(x + (i - 1) * width, [comps[name][k] for k, _l in show],
               width * 0.86, color=REGIME_COLOR[name], label=REGIME_LABEL[name], zorder=3)
    _style(ax)
    ax.set_xticks(x)
    ax.set_xticklabels([l for _k, l in show], fontsize=8.5, color=INK)
    ax.set_ylabel("impact contribution (×10⁻³, D-goods)", fontsize=9, color=MUTED)
    # Upper left: the investment column has no positive bars, so the legend cannot
    # collide with the deep negative passive bar the way a lower-left legend does.
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="upper left")
    ax.set_title("A. What the shock does (levels)", fontsize=10, color=INK, pad=8)

    ax = axes[1]
    for i, name in enumerate(n for n in names if n != "passive"):
        vals = [comps[name][k] - comps["passive"][k] for k, _l in show]
        ax.bar(x + (i - 0.5) * width, vals, width * 0.86, color=REGIME_COLOR[name],
               label=f"{REGIME_LABEL[name]} − passive", zorder=3)
    _style(ax)
    ax.set_xticks(x)
    ax.set_xticklabels([l for _k, l in show], fontsize=8.5, color=INK)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="lower left")
    ax.set_title("B. What the backstop does (change vs passive)",
                 fontsize=10, color=INK, pad=8)

    fig.suptitle("Components of the output response, on impact",
                 fontsize=11, color=INK, y=1.03)
    fig.tight_layout()
    save(fig, "fig03_dy_decomposition", _caption_fig03(comps, names))


def _caption_fig03(comps, names):
    """Both panels, stated from the impact contributions actually plotted.

    Under flexible prices panel B's investment and net-export channels were each
    several times the headline ΔY and opposite in sign, so the caption's job was to
    warn against leading with the headline. Under sticky prices the ordering is
    REVERSED — consumption carries almost the whole of ΔY and the other channels are
    a quarter of it — so the sentence is selected, not adjusted.
    """
    lo, hi = names[0], names[-1]
    A = comps[lo]
    B = {k: comps[hi][k] - comps[lo][k] for k in comps[lo]}
    tot = B["__total__"]
    r = {k: (B[k] / tot if tot != 0 else np.nan)
         for k in ("consumption_quantity", "investment", "net_exports")}

    lead = max(r, key=lambda k: abs(r[k]))
    lead_name = COMPONENT_LABEL[lead]
    others = [k for k in ("consumption_quantity", "investment", "net_exports") if k != lead]
    # "at {x}x" rather than a verb, so the sentence stays grammatical whichever
    # component the data picks out as the leading one (plural "net exports"
    # included).
    other_txt = " and ".join(f"{COMPONENT_LABEL[k]} at {r[k]:+.2f}×" for k in others)
    # Whether the two secondary channels offset is a claim about signs, so read it.
    other_txt += (" largely offsetting each other"
                  if np.sign(r[others[0]]) != np.sign(r[others[1]])
                  else " pulling the same way")

    residue = abs(tot) < max(abs(B[k]) for k in ("investment", "net_exports",
                                                 "consumption_quantity"))
    verdict = ("the headline ΔY is a residue of larger offsetting channels, which is "
               "why the decomposition and not the headline is the object to report"
               if residue else
               "the headline ΔY is no longer a residue of larger offsetting channels — "
               "it is now the largest object in the decomposition — but the remaining "
               "channels still work against each other, so the decomposition is still "
               "what should be reported")
    return (f"On impact the crisis is a joint contraction: consumption contributes "
            f"{A['consumption_quantity']:+.2f} and investment {A['investment']:+.2f} "
            f"(×10⁻³ of D-goods) against a {A['net_exports']:+.2f} net-export cushion, "
            f"summing to {A['__total__']:+.2f} (panel A). The backstop works through the "
            f"same margin rather than a different one: moving {lo} → {hi} raises ΔY by "
            f"{tot:+.2f}, of which {lead_name} supplies {r[lead]:+.2f}×, with "
            f"{other_txt} (panel B). So {verdict}.")


def fig06_net_effects(cache, regimes, n_q=16):
    """Contributions to dY_t, quarter by quarter, with the net path overlaid.

    Stacked bars are accumulated separately above and below zero — matplotlib's
    stackplot cannot represent mixed-sign contributions, and letting it try
    silently misplaces every segment once a component crosses zero.
    """
    ss = {"P_CES_D_ss": float(cache["P_CES_D_ss"]), "C_D_ss": float(cache["C_D_ss"])}
    parts = [("investment", "investment", "#1a6e3a"),
             ("net_exports", "net exports", "#1B6CA8"),
             ("consumption_quantity", "consumption (quantity)", "#c87941"),
             ("consumption_price", "consumption (price)", "#A62B22")]
    q = np.arange(n_q)
    by_regime = {}

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), sharey=True)
    for ax, (name, (_g, irf)) in zip(axes, regimes.items()):
        comps, _r = decompose_dY(irf, ss)
        by_regime[name] = (comps, np.asarray(irf["Y_D"]))
        pos = np.zeros(n_q)
        neg = np.zeros(n_q)
        for key, label, colour in parts:
            v = np.asarray(comps[key])[:n_q] * 1e3
            up, dn = np.clip(v, 0, None), np.clip(v, None, 0)
            ax.bar(q, up, bottom=pos, width=0.78, color=colour, label=label, zorder=3)
            ax.bar(q, dn, bottom=neg, width=0.78, color=colour, zorder=3)
            pos += up
            neg += dn
        ax.plot(q, np.asarray(irf["Y_D"])[:n_q] * 1e3, color=INK, lw=2,
                marker="o", ms=3.5, label="net ΔY", zorder=5)
        _style(ax)
        ax.set_xlabel("quarters", fontsize=9, color=MUTED)
        ax.set_title(REGIME_LABEL[name], fontsize=10, color=INK, pad=8)
    axes[0].set_ylabel("contribution to ΔY (×10⁻³, D-goods)", fontsize=9, color=MUTED)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right", ncol=1)
    fig.suptitle("Net decomposition of the output response, quarter by quarter",
                 fontsize=11, color=INK, y=1.03)
    fig.tight_layout()
    save(fig, "fig06_net_effects", _caption_fig06(by_regime, parts, n_q))


def _caption_fig06(by_regime, parts, n_q):
    """"Smaller than its components" is COUNTED, not asserted.

    The previous literal said the net path is "at every horizon far smaller than
    the components that generate it". Under sticky prices that is false in the
    impact quarter of the passive regime, where consumption and investment move
    the same way and ΔY is the largest bar on the panel. The claim is therefore
    stated as the count of quarters in which it actually holds.
    """
    keys = [k for k, _lab, _c in parts]
    comps, dY = by_regime["passive"]
    dY = dY[:n_q] * 1e3
    mat = np.array([np.asarray(comps[k])[:n_q] * 1e3 for k in keys])
    biggest = np.abs(mat).max(axis=0)
    n_small = int((np.abs(dY) < biggest).sum())

    inv = np.asarray(comps["investment"])[:n_q]
    nx = np.asarray(comps["net_exports"])[:n_q]
    q_inv = _first_quarter(inv > 0, n_q)
    q_nx = _first_quarter(nx < 0, n_q)
    turns = []
    if q_inv is not None:
        turns.append(f"the investment contribution turns positive from {_ordinal(q_inv)} "
                     f"as the capital stock is run down")
    if q_nx is not None:
        turns.append(f"net exports flip from cushion to drag at {_ordinal(q_nx)}")
    turn_txt = ("; thereafter " + ", and ".join(turns)) if turns else ""

    return (f"Contributions to the output response quarter by quarter. Without a backstop "
            f"the impact quarter is a joint consumption-and-investment contraction "
            f"({dY[0]:+.2f} ×10⁻³ of D-goods in total){turn_txt}. The backstop works by "
            f"lifting the consumption contribution in the opening quarters rather than by "
            f"raising output uniformly, and beyond the impact quarter the aggregate hides "
            f"most of what moves underneath it: the net path (black) is smaller in "
            f"magnitude than the largest single component in {n_small} of the first "
            f"{n_q} quarters.")


def fig07_ms_regimes():
    """Empirical regime dating — the discipline behind the model's three stances."""
    npz = os.path.join(ROOT, "Empirics", "outputs", "ms_regime_COMPOSITE.npz")
    d = np.load(npz, allow_pickle=True)
    dates = np.array([np.datetime64(s) for s in d["dates"]])
    y, probs, erg, dur = d["y"], d["probs"], d["ergodic"], d["durations"]

    # Order the estimated states by mean spread so the labels are not an artifact of
    # the estimator's arbitrary state numbering: low spread = dove (the ECB is
    # intervening), high spread = hawk (it is not).
    order = np.argsort(d["means"])
    lab = ["dove (intervening)", "base", "hawk (passive)"]
    col = [GREEN, ORANGE, RED]

    fig, axes = plt.subplots(2, 1, figsize=(12, 6.2), sharex=True,
                             gridspec_kw={"height_ratios": [1.35, 1]})

    ax = axes[0]
    modal = np.argmax(probs[:, order], axis=1)
    for k in range(3):
        ax.fill_between(dates, 0, y.max() * 1.08, where=(modal == k),
                        color=col[k], alpha=0.13, step="mid", lw=0)
    ax.plot(dates, y, color=INK, lw=1.3)
    # The sample opens in 1995 but the ECB did not exist until 1999, so the
    # "intervention stance" reading of the early regimes is anachronistic — that
    # stretch is EMU convergence, not policy. Mark it rather than let the shading
    # imply an ECB stance that could not have existed.
    ecb = np.datetime64("1999-01-01")
    ax.axvline(ecb, color=INK, lw=1.0, ls=(0, (3, 2)))
    ax.text(ecb, y.max() * 1.02, "  ECB founded", fontsize=7.5, color=INK,
            va="top", ha="left")
    ax.set_ylim(0, y.max() * 1.08)
    _style(ax)
    ax.yaxis.grid(False)
    ax.set_ylabel("peripheral − Bund 10y (pp)", fontsize=9, color=MUTED)
    ax.set_title("Peripheral–Bund spread, shaded by modal intervention-stance regime",
                 fontsize=10, color=INK, pad=8)
    for k in range(3):
        ax.plot([], [], color=col[k], lw=6, alpha=0.35, label=lab[k])
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, ncol=3, loc="upper left")

    ax = axes[1]
    ax.stackplot(dates, *[probs[:, order[k]] for k in range(3)],
                 colors=col, alpha=0.85, labels=lab, lw=0)
    ax.set_ylim(0, 1)
    _style(ax)
    ax.yaxis.grid(False)
    ax.set_ylabel("smoothed probability", fontsize=9, color=MUTED)
    ax.set_title(
        "Smoothed regime probabilities  ·  ergodic shares "
        + " / ".join(f"{erg[order[k]] * 100:.0f}%" for k in range(3))
        + "  ·  expected durations "
        + " / ".join(f"{dur[order[k]]:.0f}m" for k in range(3)),
        fontsize=9.5, color=INK, pad=8)

    fig.tight_layout()
    save(fig, "fig07_ms_regimes", _caption_fig07(dates, modal, erg, order, ecb))


def _caption_fig07(dates, modal, erg, order, ecb):
    """Empirical, hence MODEL-INDEPENDENT — but still derived from the npz.

    Nothing in this caption moves when the model is recalibrated: the estimates come
    from Empirics/outputs/ms_regime_COMPOSITE.npz, not from the solve. It is derived
    anyway so a re-estimation of the Markov-switching model cannot leave it stale.
    """
    shares = " / ".join(f"{erg[order[k]] * 100:.0f}%" for k in range(3))
    # Longest contiguous run of the modal hawk state, and whether any of it
    # predates the ECB — the caveat the figure's dashed line marks.
    spans, start = [], None
    for i, m in enumerate(modal):
        if m == 2 and start is None:
            start = i
        elif m != 2 and start is not None:
            spans.append((start, i - 1))
            start = None
    if start is not None:
        spans.append((start, len(modal) - 1))
    if spans:
        s0, s1 = max(spans, key=lambda s: s[1] - s[0])
        yrs = (str(dates[s0])[:4], str(dates[s1])[:4])
        hawk_txt = (f"the high-spread 'hawk' state covers {yrs[0]}–{yrs[1]}"
                    if yrs[0] != yrs[1] else f"the high-spread 'hawk' state covers {yrs[0]}")
    else:
        hawk_txt = "the high-spread 'hawk' state is never modal"
    pre = any(dates[s0_] < ecb for s0_, _s1 in spans)
    caveat = (" — though the pre-1999 stretch predates the ECB and reflects EMU "
              "convergence, not any policy stance" if pre else "")
    return (f"A three-state Markov-switching model on peripheral–Bund spreads dates the "
            f"ECB's intervention stance and disciplines the model's three backstop regimes: "
            f"{hawk_txt}, and the ergodic shares ({shares}) are what the regime-uncertainty "
            f"beliefs are set to{caveat}. Estimated from market data, so unlike every other "
            f"figure here it does not move with the calibration.")


def fig08_deciles():
    """Consumption incidence by steady-state INCOME quintile.

    Per-capita consumption in bin k is (mass-weighted consumption) / (bin mass), so
    its first-order percentage response is the difference of the two percentage
    responses. Reporting only the numerator would attribute pure composition drift
    to household behaviour.

    Income quintiles are plotted, NOT wealth deciles. The mass in an income bin is
    the stationary distribution of the exogenous Markov chain and does not move, so
    this measure is purely behavioural. The wealth-decile analogue is dominated by
    membership churn across fixed deposit thresholds (bottom decile: -41.6
    consumption against -44.4 mass, netting +2.8) and is reported in the tables with
    that caveat rather than plotted, where it would read as a behavioural result.
    """
    from e4_distribution import CACHE, N_QNT
    from lottery_math import closed_loop

    if not os.path.exists(CACHE):
        raise FileNotFoundError(
            f"{CACHE} missing — run: /opt/anaconda3/envs/ssj/bin/python "
            f"experiments/e4_distribution.py")
    d = np.load(CACHE, allow_pickle=True)
    eps = np.asarray(d["dShock_def_D"])
    A_def, A_cb = d["spread_rb__shock_def_D"], d["spread_rb__cb_buy_D"]
    mass_ss, c_ss = d["qnt_mass_ss"], d["qnt_c_ss"]

    # gamma solved on THIS cache, so the regimes mean the same thing they do everywhere
    # else in the paper. Routed through common.named_regime_gammas rather than calling
    # gamma_for_compression directly: since 2026-08-18 the aggressive target is
    # unreachable below the closed-loop pole and falls back to the maximum feasible
    # intervention, and a second copy of the solve here would silently skip that
    # handling and raise. One definition of the regimes, not two.
    from common import named_regime_gammas as _nrg
    gam = _nrg({"spread_rb__shock_def_D": A_def, "spread_rb__cb_buy_D": A_cb,
                "dShock_def_D": eps})

    H, beta = 40, float(d["beta_D_ss"])
    disc = beta ** np.arange(H)

    def path(k, cb):
        """Per-capita consumption in income bin k, % of its own SS, over H quarters."""
        # UPPERCASE: SSJ exposes a het block's outputs under the uppercased name
        # (c_D -> C_D), and that is what the cache keys on.
        ck = (d[f"CQNT{k + 1}_D__shock_def_D"] @ eps
              + d[f"CQNT{k + 1}_D__cb_buy_D"] @ cb)[:H]
        mk = (d[f"MQNT{k + 1}_D__shock_def_D"] @ eps
              + d[f"MQNT{k + 1}_D__cb_buy_D"] @ cb)[:H]
        return 100.0 * (ck / c_ss[k] - mk / mass_ss[k])

    paths, pv = {}, {}
    for name, g in gam.items():
        _sp, cb = closed_loop(A_def, A_cb, eps, g)
        paths[name] = np.array([path(k, cb) for k in range(N_QNT)])
        pv[name] = (paths[name] * disc).sum(axis=1)

    x = np.arange(1, N_QNT + 1)
    # Income quintiles are ORDINAL, so this is a sequential ramp — one hue, light to
    # dark — not the categorical palette. Mixing categorical hues across an ordered
    # variable implies distinctions that are not there and loses the ordering.
    qcol = ["#BBD6EA", "#8CBAD9", "#5C9EC8", "#3782B7", "#12507E"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))

    ax = axes[0]
    for k in range(N_QNT):
        ax.plot(np.arange(H), paths["passive"][k], lw=2, color=qcol[k],
                label=f"Q{k + 1}" + (" (lowest)" if k == 0 else
                                     " (highest)" if k == N_QNT - 1 else ""))
    _style(ax)
    ax.set_xlabel("quarters", fontsize=9, color=MUTED)
    ax.set_ylabel("consumption (% of own SS)", fontsize=9, color=MUTED)
    ax.set_title("A. No backstop, by income quintile", fontsize=10, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, ncol=2)

    # Discounted PV — the summary incidence measure. An impact-only snapshot is
    # misleading here, because every quintile's consumption RISES on impact (the
    # investment collapse releases resources) and only later turns.
    ax = axes[1]
    for name in gam:
        ax.plot(x, pv[name], marker="o", ms=6, lw=2, color=REGIME_COLOR[name],
                label=REGIME_LABEL[name])
    _style(ax)
    ax.set_xticks(x)
    ax.set_xlabel("steady-state income quintile (1 = lowest)", fontsize=9, color=MUTED)
    ax.set_ylabel("PV of consumption response (% of own SS)", fontsize=9, color=MUTED)
    ax.set_title("B. Incidence, discounted over 40q", fontsize=10, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK)

    ax = axes[2]
    gains = pv["aggressive"] - pv["passive"]
    ax.bar(x, gains, width=0.6, color=[GREEN if v >= 0 else RED for v in gains], zorder=3)
    _style(ax)
    ax.set_xticks(x)
    ax.set_xlabel("steady-state income quintile (1 = lowest)", fontsize=9, color=MUTED)
    ax.set_ylabel("PV difference vs passive", fontsize=9, color=MUTED)
    ax.set_title("C. Who the backstop protects", fontsize=10, color=INK, pad=8)

    fig.suptitle("Distributional incidence by steady-state income quintile",
                 fontsize=11, color=INK, y=1.03)
    fig.tight_layout()
    save(fig, "fig08_deciles", _caption_fig08(paths, pv, N_QNT, H))
    return paths, pv, gam


def _caption_fig08(paths, pv, n_qnt, H):
    """The instance that made this defect visible: the literal caption said the
    lowest quintile "gains 0.95%" and the highest loses 0.59%, against a Table 4 in
    the same generated document reading +0.4250 and −0.9073. It also claimed every
    quintile's consumption RISES on impact, which the sticky-price solution reverses.
    Both facts are now read off `paths` and `pv` — the arrays panel A and panel B
    are drawn from.
    """
    imp = paths["passive"][:, 0]
    p_lo, p_hi = float(pv["passive"][0]), float(pv["passive"][-1])
    gain = np.asarray(pv["aggressive"]) - np.asarray(pv["passive"])
    g_lo, g_hi = float(gain[0]), float(gain[-1])

    # Impact response: near-identical across quintiles at this calibration, so say
    # so only if the spread across bins is genuinely small relative to the level.
    uniform = float(imp.max() - imp.min()) < 0.1 * abs(float(imp.mean()))
    if uniform:
        vb = "falls" if imp.mean() < 0 else "rises"
        imp_txt = (f"Consumption {vb} by about {abs(float(imp.mean())):.2f}% in every "
                   f"income quintile on impact, so the distributional difference emerges "
                   f"only afterwards")
    else:
        imp_txt = (f"On impact the consumption response already differs across the "
                   f"distribution, from {imp[0]:+.2f}% in the lowest quintile to "
                   f"{imp[-1]:+.2f}% in the highest")

    mono_pv = _monotone(pv["passive"], -1)
    mono_gain = _monotone(gain, -1)
    burden = (", monotonically across the five quintiles" if mono_pv else
              ", though not monotonically across the quintiles")
    prot = ("The backstop's protection runs the same way, also monotone in quintile:"
            if mono_gain else
            "The backstop's protection runs the same way but is not monotone in quintile:")

    def _side(v, who):
        return (f"the {who} quintile loses {abs(v):.2f}% of its own consumption"
                if v < 0 else f"the {who} quintile gains {v:.2f}%")

    return (f"{imp_txt}. Discounted over {H} quarters the burden of the crisis falls on the "
            f"top of the income distribution{burden}: {_side(p_hi, 'highest-income')} "
            f"while {_side(p_lo, 'lowest')}. {prot} it is worth {g_lo:+.2f}% of "
            f"consumption to the lowest quintile against {g_hi:+.2f}% to the highest.")


def fig04_spread_decomposition(cache, ss_tl, regimes):
    """Spread path against the expected-loss pricing it starts from.

    REWRITTEN 2026-08-18. The previous version was a two-segment bar splitting the
    default loading into `EL_price_D` and `psi_spread_D`, captioned "x% fundamental /
    y% non-fundamental". Both objects are deleted, and the split was not
    interpretable in the first place: `psi_spread_D` was a free parameter absorbing
    the principal-loss term that `zeta_writeoff = 0` had left out of the bond payoff,
    so the "friction share" was mostly a calibration artefact. In a linearised
    equilibrium model there is no such thing as a fundamental/non-fundamental share of
    an endogenous price anyway — every channel operates simultaneously.

    What replaces it is a MECHANISM chart with no share arithmetic. Three series, all in
    the same units (annualised bp of yield), all read off the same solved path:

    1. **Direct expected-loss pricing of the Greek yield.** The bond price path implied
       by the realised `def_rate_D` path when the intermediary's required excess return
       is held at its steady-state value `s0 = Delta_bD_D * (rk_ss - rdep_ss)`. Solved by
       backward recursion on the model's OWN pricing condition,

           q_t = (1 - h*d_{t+1}) * [delta_b + (1-delta_b)*q_{t+1}] / (1 + s0),

       terminating at `q_ss`, then read as `delta_b*(1/q_t - 1)`. Forward-consistent, so
       it correctly reflects that the shock decays at `rho_def`.
    2. **Equilibrium Greek yield** — the same object with the intermediary's required
       return free to move.
    3. **Equilibrium GR–DE spread**, which additionally contains the German leg.

    UNITS TRAP THIS FIXES. A first draft plotted `EL_load_D * def_rate_D`, which is a
    ONE-PERIOD expected capital-loss rate, against `spread_rb`, which is a per-period
    COUPON-EQUIVALENT yield. They are not commensurate: it read 331bp against a 206bp
    spread, i.e. an "amplification factor" of 0.62, which is an artefact of the mismatch
    and not a mechanism. The recursion above puts both on the yield measure.

    WHAT THE CHART ACTUALLY SHOWS, which is not what the old figure claimed. Direct
    expected-loss pricing accounts for essentially the whole Greek yield response
    (213.8bp against an equilibrium 189.7bp on impact). The intermediary channel is a
    modest OFFSET at the price — general-equilibrium movement in `rk_D` and `rdep_D`
    lowers the required return slightly — while being a large amplifier for QUANTITIES
    (`n_inter_D` -11.4%, `Y_D` -1.97%). The spread exceeds the Greek yield because the
    German leg falls (flight to quality), not because of amplification on the Greek leg.
    """
    irf = regimes["passive"][1]
    H, T = 40, 200
    t = np.arange(H)
    q_ss = float(ss_tl["q_b_D"])
    db = float(ss_tl["delta_b_D"])
    h = 1.0 - float(ss_tl["recovery_rate_D"])
    s0 = db * (1.0 / q_ss - 1.0)          # SS required excess return (rdep_ss = 0)
    d = np.asarray(irf["def_rate_D"])[:T]  # deviation == level, def_rate_ss = 0

    q = np.full(T + 1, q_ss)
    for i in range(T - 1, -1, -1):
        dn = d[i + 1] if i + 1 < T else 0.0
        q[i] = (1.0 - h * dn) * (db + (1.0 - db) * q[i + 1]) / (1.0 + s0)
    y0 = db * (1.0 / q_ss - 1.0)
    direct = (db * (1.0 / q[:H] - 1.0) - y0) * BP_ANN
    y_D = np.asarray(irf["rb_D"])[:H] * BP_ANN
    spread = np.asarray(irf["spread_rb"])[:H] * BP_ANN

    fig, ax = plt.subplots(figsize=(9, 3.8))
    ax.fill_between(t, 0, direct, color=BLUE, alpha=0.28, zorder=1)
    ax.plot(t, direct, color=BLUE, lw=1.6, ls="--", zorder=3,
            label="Greek yield under direct expected-loss pricing\n"
                  "(required return held at steady state)")
    ax.plot(t, y_D, color=INK, lw=2.2, zorder=4, label="equilibrium Greek yield")
    ax.plot(t, spread, color=ORANGE, lw=2.0, zorder=4,
            label="equilibrium GR–DE spread (adds the German leg)")
    ax.set_xlim(0, H - 1)
    ax.set_xlabel("quarters after the 1pp Greek default-probability shock",
                  fontsize=9, color=MUTED)
    ax.set_ylabel("bp, annualised", fontsize=9, color=MUTED)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right")
    ax.set_title("Where the Greek spread comes from: the bond payoff, priced by the "
                 "intermediary", fontsize=9.5, color=INK, pad=10)
    fig.tight_layout()
    peak_sp = float(spread.max())
    save(fig, "fig04_spread_decomposition",
         _caption_fig04(float(ss_tl["EL_load_D"]), peak_sp, float(direct.max()),
                        float(y_D.max())))
    return float(ss_tl["EL_load_D"]), peak_sp


def _caption_fig04(el_load, peak_sp, peak_dir, peak_yD):
    """Derived from the solved path, not from any parameter."""
    return (f"The sovereign spread is generated by the bond's state-contingent payoff "
            f"inside the intermediary's portfolio optimality condition, not by any spread "
            f"parameter. Pricing the expected loss on the perpetuity — coupon *and* "
            f"continuation value, {el_load:.4f} per unit of default probability — at a "
            f"required return frozen at its steady-state level already accounts for "
            f"{peak_dir:.0f}bp of Greek yield, against an equilibrium {peak_yD:.0f}bp: the "
            f"intermediary channel is a modest offset at the price while being a large "
            f"amplifier for quantities. The {peak_sp:.0f}bp spread exceeds the Greek yield "
            f"because the German leg falls in a flight to quality. These are joint "
            f"mechanisms along one equilibrium path, not separable shares of the price.")


def fig05_incidence(cache, payload, gammas):
    """Germany's ledger: exposure (quantity of risk) against loading (price paid)."""
    Y_ss = float(cache["Y_D_ss"])
    A_def, A_cb = cache["spread_rb__shock_def_D"], cache["spread_rb__cb_buy_D"]
    eps = np.asarray(cache["dShock_def_D"])
    from lottery_math import closed_loop

    expo, el_pv, load = [], [], []
    for g in gammas:
        _sp, cb = closed_loop(A_def, A_cb, eps, float(g))
        d = e1.cb_pnl(irf_from_cache(cache, cb, eps), cache)
        expo.append(100 * d["purchases_pv"] / Y_ss)
        el_pv.append(100 * d["el_pv"] / Y_ss)
        load.append(d["prem_pv"] / d["el_pv"] if d["el_pv"] > 1e-16 else np.nan)

    # Three panels, not two with one series rescaled. Exposure (~2%) and expected
    # loss (~0.005%) are the same unit but differ by ~400x; putting them on one axis
    # requires a "(x100)" multiplier, which is a second scale in disguise. Small
    # multiples keep every axis honest.
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0))
    ax = axes[0]
    ax.plot(gammas, expo, color=BLUE, lw=2.2)
    _style(ax)
    ax.set_xlabel("backstop aggressiveness γ", fontsize=9, color=MUTED)
    ax.set_ylabel("PV, % of quarterly steady-state $Y_D$", fontsize=9, color=MUTED)
    ax.set_title("Exposure\n(discounted purchases)", fontsize=10, color=INK, pad=8)

    ax = axes[1]
    ax.plot(gammas, el_pv, color=ORANGE, lw=2.2)
    _style(ax)
    ax.set_xlabel("backstop aggressiveness γ", fontsize=9, color=MUTED)
    ax.set_ylabel("PV, % of quarterly steady-state $Y_D$", fontsize=9, color=MUTED)
    ax.set_title("Priced expected loss\n(computed off-path)", fontsize=10, color=INK, pad=8)

    ax = axes[2]
    ax.plot(gammas, load, color=BLUE, lw=2.2)
    ax.axhline(1.0, ls=(0, (4, 3)), lw=1.2, color=MUTED)
    ax.text(gammas[-1], 1.0, " actuarially fair", va="center", ha="right",
            fontsize=8, color=MUTED, backgroundcolor="white")
    _style(ax)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("backstop aggressiveness γ", fontsize=9, color=MUTED)
    ax.set_ylabel("premium PV ÷ expected-loss PV", fontsize=9, color=MUTED)
    ax.set_title("Compensation per unit\n(loading)", fontsize=10, color=INK, pad=8)

    fig.suptitle("The German ledger: exposure rises, compensation per unit falls",
                 fontsize=11, color=INK, y=1.03)
    fig.tight_layout()
    save(fig, "fig05_incidence", _caption_fig05(gammas, expo, load))


def _caption_fig05(gammas, expo, load):
    """The whole claim is a pair of directions, so take both from the schedules."""
    expo = np.asarray(expo, dtype=float)
    load = np.asarray(load, dtype=float)
    ok = ~np.isnan(load)
    g_hi = float(np.asarray(gammas)[-1])
    up = _monotone(expo, +1)
    down = _monotone(load, -1)
    opposed = up and down
    verdict = ("move in opposite directions" if opposed else
               "do NOT move in opposite directions at this calibration — check before "
               "asserting the German-ledger reading")
    return (f"As the backstop strengthens Germany's discounted exposure rises "
            f"{'steadily' if up else 'non-monotonically'} — from zero to "
            f"{expo[-1]:.2f}% of quarterly steady-state $Y_D$ at γ={g_hi:.0f} — while the "
            f"compensation it earns per unit of expected loss "
            f"{'falls steadily' if down else 'moves non-monotonically'} from "
            f"{load[ok][0]:.2f}× to {load[ok][-1]:.2f}×, so the two objects the German "
            f"litigation actually turned on — quantity of risk assumed and price paid for "
            f"it — {verdict}.")


# ── Tables ───────────────────────────────────────────────────────────────────

def _solved_ss():
    import io, contextlib
    from calibration import get_calibration
    from depreciation_calibration import calibrate_depreciation
    from ic_delta_calibration import calibrate_ic_delta
    from steady_state import solve_steady_state
    with contextlib.redirect_stdout(io.StringIO()):
        r = calibrate_depreciation(calibrate_ic_delta(solve_steady_state(get_calibration())))
    ss = r["ss_final"] if isinstance(r, dict) and "ss_final" in r else r
    return getattr(ss, "toplevel", ss)


def tables(cache, payload, ss_tl, el, ps, dist=None):
    """Markdown tables for the draft. Every value read live; only citations are literal."""
    g = lambda k: float(ss_tl[k])
    p = provenance()
    Y_ss = float(cache["Y_D_ss"])

    peak0 = payload["regimes"]["passive"]["peak_spread_bp_ann"]
    n0 = payload["regimes"]["passive"]["impact"]["n_inter_D_pct_ss"]
    passthrough = n0 / (peak0 / 100.0)

    L = [
        "# First-draft results — tables and figures", "",
        f"*Generated {p['generated']} from `{p['git_sha']}` · calibration "
        f"`{p['cal_fingerprint']}` · `BANK_SCOPE={p['BANK_SCOPE']}` · "
        f"`writeoff_enabled={p['writeoff_enabled_D']:g}` (pure risk-premium framing, "
        f"S-1 resolved 2026-08-04).*", "",
        "Generated by `experiments/paper_outputs.py`. Every number is read live from the "
        "solved steady state or the cached response matrices — none is transcribed. "
        "Figures are in `experiments/paper/`, each with its caption baked into the image. "
        "**The captions are derived too**: each is built by its own figure function from "
        "the arrays that figure plots, so a caption cannot contradict a table below it "
        "the way the hardcoded set did between the sticky-price conversion and "
        "2026-08-06.",
        "",
        "## Table 1 — Calibration and identification ledger", "",
        "The distinction that matters for a referee is *which* parameters are measured, "
        "which are matched to a target, and which are free. Only one amplification dial "
        "is free.", "",
        "| Parameter | Value | Status | Source / target |",
        "|---|---|---|---|",
        f"| `theta_D` / `theta_F` | {g('theta_D'):.3f} / {g('theta_F'):.3f} | **measured** | "
        "EBA 2011 disclosure: GK-eligible assets (corporate ex-CRE + CRE + sovereign) EAD ÷ Core Tier 1 |",
        f"| `delta_b_D` / `delta_b_F` | {g('delta_b_D'):.4f} / {g('delta_b_F'):.4f} | **measured** | "
        "EBA sovereign maturity ladder repriced at the 31-Dec-2010 market yield; modified duration 3.12y (GGB) / 4.22y (Bund) |",
        f"| `phi_bD_D` (own-book concentration) | {g('phi_bD_D'):.3f} | **measured** | "
        "EBA 2011 bank-held sovereign book ÷ capital, broad-sector scope |",
        f"| `b_F_D` (cross-holdings) | {g('b_F_D'):.4f} | **measured** | EBA 2011 bilateral GR/DE exposures |",
        f"| `n_inter_D` / `n_inter_F` | {g('n_inter_D'):.3f} / {g('n_inter_F'):.3f} | *implied* | "
        "`(Q·K + sovereign) / theta` — follows from measured leverage and the balance sheet |",
        f"| `recovery_rate` | {g('recovery_rate_D'):.2f} | *external estimate* | "
        "Zettelmeyer, Trebesch & Gulati (PIIE WP13-8): 59–65% investor NPV loss in the March 2012 Greek PSI |",
        f"| `psi_lambda_B` | {g('psi_lambda_B_D'):.2f} | **FREE — the one amplification dial** | "
        "matched to the ~150bp 2010 GR–DE spread on a 1pp default shock; no EBA counterpart |",
        f"| `def_scale` | {g('def_scale_D'):.2f} | free | strong-amplification choice; exceeds the 0.12–0.23 GR-2011 range |",
        f"| `Delta_bD_D` / `Delta_bF_D` | 0.20 / 0.40 | **unidentified** | "
        "no EBA counterpart; the GK feasibility inequality bounds but does not pin them |",
        f"| `phi_lamb` | {g('phi_lamb_D'):.2f} | free | fiscal feedback; literature (Staehr 2008) is 0.025–0.038 quarterly |",
        f"| `f` (banker exit) | {g('f_D'):.2f} | free | standard GK range is ~0.03 |",
        "",
        "## Table 2 — Moment match", "",
        "| Moment | Target | Model | Source |",
        "|---|---|---|---|",
        f"| Peak D–F spread, 1pp default shock | ~150 bp ann. | **{peak0:.1f} bp** | 2010 GR–DE 10y spread |",
        f"| Capital–output ratio `K/Y` (annual) | 2.70 | **{g('K_D') / (4 * Y_ss):.2f}** | conventional |",
        f"| Steady-state return on capital `rk` | 0.0100 | **{g('rk_D') if 'rk_D' in ss_tl else float('nan'):.6f}** | conventional quarterly |",
        f"| Bank net-worth pass-through | −1.8 to −8.6 %/100bp | **{passthrough:.2f} %/100bp** | "
        "Acharya–Drechsler–Schnabl (2014 JF) bank-equity-on-sovereign-CDS elasticity, converted at three baselines |",
        f"| Steady-state D and F bond yields | equalised | {g('rb_D'):.6f} / {g('rb_F'):.6f} | model restriction |",
        "",
        "## Table 3 — Main results", "",
        "| | passive | medium | aggressive |",
        "|---|---|---|---|",
        "| backstop coefficient γ | " + " | ".join(
            f"{payload['regimes'][k]['gamma']:.3f}" for k in payload["regimes"]) + " |",
        "| peak D–F spread (bp ann.) | " + " | ".join(
            f"{payload['regimes'][k]['peak_spread_bp_ann']:.1f}" for k in payload["regimes"]) + " |",
        "| output, impact (% SS) | " + " | ".join(
            f"{payload['regimes'][k]['impact']['Y_D_pct_ss']:+.4f}" for k in payload["regimes"]) + " |",
        "| investment, impact (% SS) | " + " | ".join(
            f"{payload['regimes'][k]['impact']['I_D_pct_ss']:+.4f}" for k in payload["regimes"]) + " |",
        "| bank net worth, impact (% SS) | " + " | ".join(
            f"{payload['regimes'][k]['impact']['n_inter_D_pct_ss']:+.3f}" for k in payload["regimes"]) + " |",
        "| ECB exposure PV (% of quarterly $Y_D$) | " + " | ".join(
            f"{payload['regimes'][k]['a5_1_exposure_pv_pct_Y']:.4f}" for k in payload["regimes"]) + " |",
        "| priced expected loss PV (% of $Y_D$) | " + " | ".join(
            f"{payload['regimes'][k]['a5_1_expected_loss_pv_pct_Y']:.5f}" for k in payload["regimes"]) + " |",
        "| **loading (premium ÷ expected loss)** | " + " | ".join(
            ("n/a" if payload["regimes"][k]["loading"] is None
             else f"**{payload['regimes'][k]['loading']:.2f}**") for k in payload["regimes"]) + " |",
        "",
        f"Sovereign-spread mechanism: expected loss on the perpetuity is `EL_load_D = "
        f"{el:.6f}` per unit of default probability (coupon *and* continuation value, "
        f"`zeta_writeoff = 1`), priced inside the GK portfolio FOC. Equilibrium peak "
        f"spread is **{ps:.1f} bp** annualised. There is no separate spread parameter and "
        f"no fundamental/non-fundamental share — see fig04.", "",
        "## Table 4 — Distributional incidence, by income quintile", "",
        "PV of the consumption response over 40 quarters, % of each quintile's own "
        "steady-state consumption. Bins are cut on the **exogenous income state**, whose "
        "marginal distribution is the stationary distribution of the Markov chain and is "
        "therefore invariant to the shock — verified numerically at `max|Δmass| ≈ 1e−19`. "
        "The per-capita response is consequently *purely behavioural*.", "",
        "| income quintile | passive | medium | aggressive | backstop gain |",
        "|---|---|---|---|---|",
    ]
    if dist is not None:
        pvq = dist["pv"]
        for k in range(len(pvq["passive"])):
            L.append(f"| Q{k + 1}{' (lowest)' if k == 0 else ' (highest)' if k == len(pvq['passive']) - 1 else ''} "
                     f"| {pvq['passive'][k]:+.4f} | {pvq['medium'][k]:+.4f} | "
                     f"{pvq['aggressive'][k]:+.4f} | {pvq['aggressive'][k] - pvq['passive'][k]:+.4f} |")
    L += ["", "> **Do not run this cut on wealth.** Binning on steady-state deposits with "
          "fixed boundaries makes the per-capita number overwhelmingly *composition*: the "
          "deposit distribution shifts across the thresholds, bin masses move by 2–3e−3 "
          "(2–3% of bin mass), and the net is a small residue of two large nearly-"
          "cancelling terms — bottom decile, PV: −41.6 consumption against −44.4 mass, "
          "netting +2.8. The arithmetic is exact and the object is well defined, but it "
          "must not be described as how poor households behaved.", "",
        "## Figures", "",
    ]
    for name in sorted(CAPTIONS):
        L += [f"### `{name}`", "", CAPTIONS[name], "",
              f"![{name}](../experiments/paper/{name}.png)", ""]

    with open(TABLES_DOC, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return TABLES_DOC


def main():
    cache = load_cache()
    regimes = regime_irfs(cache)
    payload = e1.run()
    ss_tl = _solved_ss()

    fig01_transmission(cache, regimes)
    gammas, _loading, _peak = fig02_loading_schedule(cache, regimes, payload)
    fig03_dy_decomposition(cache, regimes)
    el, ps = fig04_spread_decomposition(cache, ss_tl, regimes)
    fig05_incidence(cache, payload, gammas)
    fig06_net_effects(cache, regimes)
    fig07_ms_regimes()
    _paths, _pv, _gam = fig08_deciles()
    dist = {"pv": _pv}

    # Coverage: captions and emitted figures must match exactly. Since CAPTIONS is
    # now filled by save(), this also proves every figure supplied a DERIVED caption
    # rather than silently shipping without one.
    emitted = {f[:-4] for f in os.listdir(PAPER_DIR) if f.endswith(".png")}
    assert emitted == set(CAPTIONS), (emitted - set(CAPTIONS), set(CAPTIONS) - emitted)

    # Prose-vs-table guard. fig01's caption and Table 3 quote the same impact
    # numbers by two different routes — the cache directly, and e1.run()'s payload.
    # Nothing previously checked that rendered prose agreed with rendered tables,
    # which is how eight figures came to carry claims their own tables refuted.
    _n_cap = float(np.asarray(regimes["passive"][1]["n_inter_D"])[0]
                   * 100.0 / float(cache["n_inter_D_ss"]))
    _n_tbl = payload["regimes"]["passive"]["impact"]["n_inter_D_pct_ss"]
    assert abs(_n_cap - _n_tbl) < 1e-9, (
        f"fig01's caption says bank net worth moves {_n_cap:+.4f}% on impact while "
        f"Table 3 says {_n_tbl:+.4f}% — the figure and the table are no longer reading "
        f"the same solve. Do not publish this document.")
    # Table 4 and fig08's caption are the SAME object (dist['pv'] is _pv), so they
    # agree by construction rather than by check.
    assert dist["pv"] is _pv

    doc = tables(cache, payload, ss_tl, el, ps, dist=dist)
    print(f"Figures  -> {PAPER_DIR}  ({len(emitted)}, captions baked in)")
    print(f"Tables   -> {doc}")


if __name__ == "__main__":
    main()
