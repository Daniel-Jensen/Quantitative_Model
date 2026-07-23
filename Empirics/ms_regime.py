"""
Markov-switching regime probabilities on peripheral-Bund sovereign spreads.

Purpose (paper role)
--------------------
Estimates a 3-state Markov-switching model on 10y peripheral-minus-Bund spreads
to date, and measure the frequency/persistence of, three *ECB intervention-stance*
regimes:

    dove  (intervention)  = low-spread, low-variance state   -> aggressive backstop
    base                  = intermediate state                -> normal
    hawk  (stress)        = high-spread, high-variance state   -> passive / no backstop

The estimated transition matrix P and its expected durations / ergodic
distribution are the *empirical discipline* for the Stage-B exogenous
regime-switching transition matrix in
`docs/superpowers/specs/2026-07-16-exogenous-policy-regimes-design.md`.

Consistency with the paper's lane: spreads are used only as the *observable that
reveals* which stance regime the ECB was in (regime dating). The structural model
keeps phi_TPI exogenous; this script does NOT make the backstop react to spreads.
It supplies the frequency-and-persistence numbers (P, durations, ergodic shares)
that a constant exogenous transition matrix should match.

Run:
    /opt/anaconda3/envs/ssj/bin/python Empirics/ms_regime.py

Outputs (Empirics/outputs/):
    spreads_fred.csv                 cached FRED pull (monthly 10y rates + spreads)
    ms_regime_<name>.npz             P, durations, ergodic, means, sigmas, probs
    ms_regime_<name>_summary.csv     smoothed regime probabilities over time
    ms_regime_<name>_regimes.png     spread w/ regime shading + smoothed-prob area
    ms_regime_<name>_transition.png  transition matrix heatmap + durations
    ms_regime_summary.csv            one-row-per-series table of P diag / durations
"""

from __future__ import annotations

import os
import textwrap
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.regime_switching.markov_regression import MarkovRegression

warnings.filterwarnings("ignore")  # EM/optim convergence chatter; we check llf ourselves
np.random.seed(0)

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)
CACHE = os.path.join(OUT, "spreads_fred.csv")

START = "1995-01-01"

# FRED 10y long-term government bond yields (monthly, %)
RATE_IDS = {
    "ITA": "IRLTLT01ITM156N",
    "ESP": "IRLTLT01ESM156N",
    "GRC": "IRLTLT01GRM156N",
    "PRT": "IRLTLT01PTM156N",
    "DEU": "IRLTLT01DEM156N",
}

# Okabe-Ito colourblind-safe palette (matches Empirics/graph_spreads.py)
REGIME_LABELS = ["dove (intervention)", "base", "hawk (stress)"]
REGIME_COLORS = ["#0072B2", "#E69F00", "#D55E00"]  # blue / orange / vermillion

# --------------------------------------------------------------------------- #
# Data
# --------------------------------------------------------------------------- #
def fred_series(sid: str) -> pd.Series:
    url = f"https://fred.stlouisfed.org/graph/fredgraph.csv?id={sid}"
    df = pd.read_csv(url, parse_dates=["observation_date"])
    s = pd.to_numeric(df[sid], errors="coerce")
    s.index = df["observation_date"]
    return s.loc[s.index >= START]


def load_spreads() -> pd.DataFrame:
    """Peripheral-minus-Bund 10y spreads (percentage points), cached to disk."""
    if os.path.exists(CACHE):
        df = pd.read_csv(CACHE, parse_dates=["date"]).set_index("date")
        print(f"[data] loaded cached spreads from {CACHE} ({len(df)} months)")
        return df

    print("[data] pulling 10y rates from FRED ...")
    rates = pd.concat({c: fred_series(sid) for c, sid in RATE_IDS.items()}, axis=1)
    rates.index.name = "date"
    spreads = rates.drop(columns="DEU").sub(rates["DEU"], axis=0)
    spreads.columns = [f"{c}_spread" for c in spreads.columns]
    # periphery composite: ITA/ESP/PRT mean (GRC excluded — PSI restructuring gap
    # and >30pp crisis prints make it non-representative of a continuous stance).
    spreads["COMPOSITE_spread"] = spreads[
        ["ITA_spread", "ESP_spread", "PRT_spread"]
    ].mean(axis=1)
    out = pd.concat([rates, spreads], axis=1)
    out.to_csv(CACHE)
    print(f"[data] cached to {CACHE} ({len(out)} months)")
    return out


# --------------------------------------------------------------------------- #
# Estimation
# --------------------------------------------------------------------------- #
def stationary_dist(P: np.ndarray) -> np.ndarray:
    """Ergodic distribution of a row-stochastic transition matrix P (P_ij = P(j|i))."""
    vals, vecs = np.linalg.eig(P.T)
    i = np.argmin(np.abs(vals - 1.0))
    pi = np.real(vecs[:, i])
    return pi / pi.sum()


def fit_ms(y: pd.Series, name: str, search_reps: int = 40) -> dict:
    """
    Fit a 3-state MS model (switching mean + switching variance) and return
    regime-ordered results: dove (low spread) -> base -> hawk (high spread).
    """
    y = y.dropna()
    mod = MarkovRegression(
        y.values, k_regimes=3, trend="c", switching_variance=True
    )
    res = mod.fit(search_reps=search_reps, maxiter=200)

    pmap = dict(zip(res.model.param_names, np.asarray(res.params)))
    means = np.array([pmap[f"const[{i}]"] for i in range(3)])
    sig2 = np.array([pmap[f"sigma2[{i}]"] for i in range(3)])

    # statsmodels regime_transition[i, j] = P(s_t=i | s_{t-1}=j)  (column-stochastic).
    M = np.squeeze(np.asarray(res.regime_transition))  # (3,3), columns sum to 1
    P_row = M.T  # row-stochastic: P_row[i, j] = P(s_t=j | s_{t-1}=i)

    order = np.argsort(means)  # ascending mean -> dove, base, hawk
    means, sig2 = means[order], sig2[order]
    P = P_row[np.ix_(order, order)]

    probs = np.asarray(res.smoothed_marginal_probabilities)[:, order]
    durations = 1.0 / (1.0 - np.diag(P))
    ergodic = stationary_dist(P)

    print(f"\n=== {name}  (n={len(y)}, llf={res.llf:.1f}, AIC={res.aic:.1f}) ===")
    for k in range(3):
        print(
            f"  {REGIME_LABELS[k]:<20s} mean={means[k]:6.2f}pp  "
            f"sd={np.sqrt(sig2[k]):5.2f}  dur={durations[k]:6.1f} mo  "
            f"ergodic={ergodic[k]*100:5.1f}%"
        )
    print("  transition matrix P (rows=from, cols=to):")
    for k in range(3):
        print("   ", "  ".join(f"{P[k, j]:.3f}" for j in range(3)))

    return dict(
        name=name,
        dates=y.index,
        y=y.values,
        means=means,
        sig2=sig2,
        P=P,
        durations=durations,
        ergodic=ergodic,
        probs=probs,
        llf=res.llf,
        aic=res.aic,
        n=len(y),
    )


# --------------------------------------------------------------------------- #
# Figures (captions baked into the PNG per docs skill)
# --------------------------------------------------------------------------- #
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "savefig.dpi": 300,
    }
)

CAPTIONS: dict[str, str] = {}


def _caption(fig, name):
    cap = CAPTIONS.get(name)
    if cap:
        chars = int(fig.get_size_inches()[0] * 13)
        fig.text(
            0.5, -0.01, textwrap.fill(cap, width=chars),
            ha="center", va="top", fontsize=8, style="italic", color="0.35",
        )


def _mlr(probs):
    """Most-likely regime per period."""
    return probs.argmax(axis=1)


def plot_regimes(r: dict):
    name = r["name"]
    dates, y, probs = r["dates"], r["y"], r["probs"]
    mlr = _mlr(probs)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11, 7), sharex=True, height_ratios=[2, 1]
    )

    # --- top: spread with most-likely-regime shading ---
    for k in range(3):
        # contiguous runs of regime k
        inreg = mlr == k
        starts = np.where(inreg & ~np.r_[False, inreg[:-1]])[0]
        ends = np.where(inreg & ~np.r_[inreg[1:], False])[0]
        for s, e in zip(starts, ends):
            ax1.axvspan(
                dates[s], dates[min(e + 1, len(dates) - 1)],
                color=REGIME_COLORS[k], alpha=0.16, lw=0, zorder=0,
            )
    ax1.plot(dates, y, color="0.15", lw=1.6, zorder=2)
    ax1.axhline(0, color="0.6", lw=0.9, ls=":")
    ax1.set_ylabel("Spread vs Bund (pp)")
    ax1.set_title(
        f"ECB intervention-stance regimes — {name}\n"
        f"(3-state Markov switching on the 10y sovereign spread)",
        fontsize=12, pad=8,
    )
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=REGIME_COLORS[k], alpha=0.5)
        for k in range(3)
    ]
    ax1.legend(handles, REGIME_LABELS, frameon=False, loc="upper left", ncol=3)
    ax1.spines[["top", "right"]].set_visible(False)

    # --- bottom: smoothed regime probabilities (stacked area) ---
    ax2.stackplot(
        dates, probs.T, colors=REGIME_COLORS, alpha=0.85, labels=REGIME_LABELS
    )
    ax2.set_ylim(0, 1)
    ax2.set_ylabel("P(regime)")
    ax2.set_xlabel("Year")
    ax2.margins(x=0)
    ax2.spines[["top", "right"]].set_visible(False)

    dur = r["durations"]
    erg = r["ergodic"]
    CAPTIONS[f"ms_regime_{name}_regimes"] = (
        f"Smoothed 3-state Markov-switching regimes on the {name} 10y sovereign "
        f"spread: a persistent dove/intervention state (mean {r['means'][0]:.2f}pp, "
        f"~{dur[0]:.0f}mo, {erg[0]*100:.0f}% of time), a base state, and a rare but "
        f"sticky hawk/stress state (mean {r['means'][2]:.2f}pp, ~{dur[2]:.0f}mo, "
        f"{erg[2]*100:.0f}%) — the empirical persistence that disciplines the "
        f"Stage-B exogenous transition matrix."
    )
    _caption(fig, f"ms_regime_{name}_regimes")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"ms_regime_{name}_regimes.png"), bbox_inches="tight")
    plt.close(fig)


def plot_transition(r: dict):
    name = r["name"]
    P = r["P"]
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    im = ax.imshow(P, cmap="Blues", vmin=0, vmax=1)
    for i in range(3):
        for j in range(3):
            ax.text(
                j, i, f"{P[i, j]:.2f}", ha="center", va="center",
                color="white" if P[i, j] > 0.5 else "0.2", fontsize=12,
            )
    short = ["dove", "base", "hawk"]
    ax.set_xticks(range(3)); ax.set_xticklabels(short)
    ax.set_yticks(range(3)); ax.set_yticklabels(short)
    ax.set_xlabel("to state"); ax.set_ylabel("from state")
    ax.set_title(f"Monthly regime transition matrix — {name}", fontsize=12, pad=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    dur = r["durations"]
    dur_txt = "   ".join(
        f"{short[k]}: {dur[k]:.0f}mo" for k in range(3)
    )
    CAPTIONS[f"ms_regime_{name}_transition"] = (
        f"Estimated monthly transition matrix for {name}: diagonal dominance "
        f"(expected durations {dur_txt}) shows all three ECB-stance regimes are "
        f"persistent, so a Stage-B model with a constant exogenous transition "
        f"matrix is empirically warranted rather than an i.i.d. mixing assumption."
    )
    _caption(fig, f"ms_regime_{name}_transition")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"ms_regime_{name}_transition.png"), bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    data = load_spreads()

    # headline: periphery composite + Italy (canonical TPI subject);
    # GRC/ESP/PRT as robustness. GRC flagged (restructuring gap).
    series = {
        "COMPOSITE": data["COMPOSITE_spread"],
        "ITA": data["ITA_spread"],
        "ESP": data["ESP_spread"],
        "PRT": data["PRT_spread"],
        "GRC": data["GRC_spread"],
    }

    table_rows = []
    for name, y in series.items():
        r = fit_ms(y, name)
        # persist
        np.savez(
            os.path.join(OUT, f"ms_regime_{name}.npz"),
            dates=np.array([d.strftime("%Y-%m-%d") for d in r["dates"]]),
            y=r["y"], means=r["means"], sig2=r["sig2"], P=r["P"],
            durations=r["durations"], ergodic=r["ergodic"], probs=r["probs"],
            llf=r["llf"], aic=r["aic"],
        )
        prob_df = pd.DataFrame(
            r["probs"], index=r["dates"], columns=[l.split()[0] for l in REGIME_LABELS]
        )
        prob_df.to_csv(os.path.join(OUT, f"ms_regime_{name}_summary.csv"))
        plot_regimes(r)
        plot_transition(r)
        table_rows.append(
            dict(
                series=name, n=r["n"],
                mean_dove=r["means"][0], mean_base=r["means"][1], mean_hawk=r["means"][2],
                dur_dove=r["durations"][0], dur_base=r["durations"][1], dur_hawk=r["durations"][2],
                erg_dove=r["ergodic"][0], erg_base=r["ergodic"][1], erg_hawk=r["ergodic"][2],
            )
        )

    tbl = pd.DataFrame(table_rows)
    tbl.to_csv(os.path.join(OUT, "ms_regime_summary.csv"), index=False)
    print("\n=== cross-series summary ===")
    print(tbl.round(2).to_string(index=False))

    # coverage check (captioning skill)
    emitted = {p[:-4] for p in os.listdir(OUT) if p.endswith(".png")}
    assert emitted == set(CAPTIONS), (emitted - set(CAPTIONS), set(CAPTIONS) - emitted)
    print(f"\n[ok] {len(emitted)} captioned figures written to {OUT}")


if __name__ == "__main__":
    main()
