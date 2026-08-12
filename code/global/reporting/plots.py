# IRF OVERLAY PLOT FOR THE RECURSIVE OMT/TPI ACTIVATION EXPERIMENT.
# Every figure is written to output/. The projection experiments print their
# own IRF tables; this is the one figure the comparison produces.
import os

import numpy as np
import matplotlib.pyplot as plt

# output/ lives at the package root (one level up from reporting/)
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")


def _save(fig, filename):
    # TIGHT-LAYOUT AND WRITE THE FIGURE TO output/.
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


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
