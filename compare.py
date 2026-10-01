"""Linear (SSJ) vs global-projection IRFs, side by side.

Reads results/SSJ/data/comparison_irfs.json and results/GLOBAL/data/comparison_irfs.json
(written by each model's compute stage) and writes, to results/COMPARISON/:
    comparison.md        impact / peak tables for both experiments
    compare_risk.png     overlay of the sovereign-risk IRFs
    compare_tfp.png      overlay of the TFP IRFs

    python3 compare.py            (or: python3 run.py compare)

Read the magnitudes with care. The two solvers are calibrated differently and the risk
shocks differ in size and persistence (each file states its own shock), so raw levels
are not a like-for-like test. The "per 1pp" columns rescale each risk response by its
own impact move in the quarterly default probability: exact for the linear model, a
first-order reading for the nonlinear one. Signs and shapes are the robust comparison.
"""
import json
import sys
from pathlib import Path

import numpy as np

RESULTS = Path(__file__).resolve().parent / "results"
DIR = RESULTS / "COMPARISON"
SOURCES = {"ssj": RESULTS / "SSJ" / "data" / "comparison_irfs.json",
           "global": RESULTS / "GLOBAL" / "data" / "comparison_irfs.json"}
# series -> label; every series is a % deviation except pd (pp)
LABELS = {"Y_D": "D output", "C_D": "D consumption", "I_D": "D investment",
          "n_D": "D bank net worth", "Q_bD": "D bond price", "pd": "D default prob. (pp/qtr)",
          "Y_F": "F output", "C_F": "F consumption", "I_F": "F investment",
          "n_F": "F bank net worth"}
TITLES = {"risk": "Sovereign-risk shock", "tfp": "TFP shock to D"}


def load(method):
    f = SOURCES[method]
    if not f.exists():
        sys.exit(f"missing {f}\n  run first:  python3 run.py {method}")
    return json.loads(f.read_text())


def common(exp, runs):
    # series both methods report, cut to the shorter horizon
    s = [r["experiments"][exp]["series"] for r in runs]
    keys = [k for k in LABELS if all(k in x for x in s)]
    h = min(len(x[k]) for x in s for k in keys)
    return keys, [{k: np.asarray(x[k][:h]) for k in keys} for x in s], h


def table(exp, runs):
    keys, (a, b), h = common(exp, runs)
    risk = exp == "risk" and "pd" in keys
    lines = [f"### {TITLES[exp]}", ""]
    lines += [f"- **{r['method']}** ({r['label']}): {r['experiments'][exp]['shock']}" for r in runs]
    lines += ["", f"Impact = quarter 0; peak = largest |response| over quarters 0-{h - 1}.", ""]
    head = "| series | SSJ impact | global impact | SSJ peak (q) | global peak (q) |"
    rule = "|---|---:|---:|---:|---:|"
    if risk:
        head += " SSJ per 1pp | global per 1pp |"
        rule += "---:|---:|"
    lines += [head + " same sign |", rule + ":---:|"]
    for k in keys:
        ia, ib = int(np.argmax(np.abs(a[k]))), int(np.argmax(np.abs(b[k])))
        row = (f"| {LABELS[k]} | {a[k][0]:+.4f} | {b[k][0]:+.4f} | "
               f"{a[k][ia]:+.4f} ({ia}) | {b[k][ib]:+.4f} ({ib}) |")
        if risk:
            row += f" {a[k][0] / a['pd'][0]:+.4f} | {b[k][0] / b['pd'][0]:+.4f} |"
        row += " yes |" if np.sign(a[k][0]) == np.sign(b[k][0]) else " **no** |"
        lines.append(row)
    return lines + [""]


def plot(exp, runs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys, (a, b), h = common(exp, runs)
    ncol = 4
    nrow = -(-len(keys) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.6 * nrow), squeeze=False)
    q = np.arange(h)
    for ax, k in zip(axes.flat, keys):
        ax.plot(q, a[k], lw=1.8, label="SSJ (linear)")
        ax.plot(q, b[k], lw=1.8, ls="--", label="global (nonlinear)")
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_title(LABELS[k], fontsize=9)
        ax.tick_params(labelsize=7)
    for ax in axes.flat[len(keys):]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=7)
    fig.suptitle(f"{TITLES[exp]} (% deviation; quarters on x-axis)", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main():
    runs = [load(m) for m in SOURCES]
    DIR.mkdir(parents=True, exist_ok=True)
    md = ["# Linear (SSJ) vs global projection", "",
          "Calibrations differ between the two solvers, so compare signs and shapes before "
          "magnitudes (see compare.py's header).", ""]
    for exp in TITLES:
        md += table(exp, runs)
        md.append(f"![{exp}](compare_{exp}.png)\n")
        plot(exp, runs, DIR / f"compare_{exp}.png")
    text = "\n".join(md)
    (DIR / "comparison.md").write_text(text)
    print(text)
    print(f"\nwritten: {DIR / 'comparison.md'}, compare_risk.png, compare_tfp.png")


if __name__ == "__main__":
    main()
