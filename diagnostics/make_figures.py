"""
Render diagnostic figures from the cached IRFs (diagnostics/*.npz).
No model solve here — reads the cache produced by solve_configs.py.

Emits (each PNG carries a baked-in one-sentence narrative caption):
  01_default_transmission_chain.png
  02_decisive_qb_networth.png
  03_tfp_control.png
"""
import os
import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams["savefig.dpi"] = 200

BASE = np.load(os.path.join(HERE, "irfs_baseline.npz"))
ZERO = np.load(os.path.join(HERE, "irfs_psilam0.npz"))

BLUE, RED = "#002147", "#8C1515"
N = 60  # quarters to plot

CAPTIONS = {
    "01_default_transmission_chain":
        "1pp sovereign-default shock, baseline psi_lambda_B=3 (blue) vs "
        "psi_lambda_B=0 (red): the identical driver def_rate_D feeds the entire "
        "transmission chain at baseline but every downstream variable collapses "
        "onto the zero axis when psi_lambda_B=0 — sovereign risk transmits ONLY "
        "through the psi_lambda_B structural terms.",
    "02_decisive_qb_networth":
        "Decisive Case-1-vs-Case-3 test: at psi_lambda_B=0 the ENDOGENOUS bond "
        "price q_b_D does not move (no forcing in its FOC), so there is no MTM "
        "loss and bank net worth n_inter_D stays flat — the price is endogenous "
        "(not stale), so this is Case 3 (missing fundamental channel), not "
        "Case 1 (stale price feeding net worth).",
    "03_tfp_control":
        "Control: a 1% TFP shock transmits normally at psi_lambda_B=0 (red "
        "tracks blue) — the model still solves and propagates a non-sovereign "
        "shock, so the null default response is specific to the sovereign-risk "
        "channel, not a globally dead solve.",
}


def g(store, shock, var):
    key = f"{shock}__{var}"
    return store[key][:N] if key in store.files else None


def save(fig, name):
    cap = CAPTIONS.get(name)
    if cap:
        chars = int(fig.get_size_inches()[0] * 14)
        fig.text(0.5, -0.02, textwrap.fill(cap, width=chars),
                 ha="center", va="top", fontsize=8, style="italic", color="0.35")
    fig.savefig(os.path.join(HERE, f"{name}.png"), bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}.png")


def chain_grid(name, shock):
    vars9 = ["def_rate_D", "q_b_D", "spread_rb", "n_inter_D", "theta_D",
             "K_D", "Y_D", "C_D", "w_D"]
    fig, axes = plt.subplots(3, 3, figsize=(12, 9))
    for ax, v in zip(axes.flat, vars9):
        b, z = g(BASE, shock, v), g(ZERO, shock, v)
        if b is not None:
            ax.plot(b, color=BLUE, lw=1.8, label="psi_lambda_B=3 (baseline)")
        if z is not None:
            ax.plot(z, color=RED, lw=1.8, ls="--", label="psi_lambda_B=0")
        ax.axhline(0, color="0.7", lw=0.6)
        ax.set_title(v, fontsize=10)
        ax.set_xlabel("quarters", fontsize=8)
    axes.flat[0].legend(fontsize=7, loc="best")
    fig.suptitle(f"Transmission chain — {shock.upper()} shock "
                 f"(baseline vs psi_lambda_B=0)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, name)


def decisive():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, v in zip(axes, ["q_b_D", "n_inter_D"]):
        b, z = g(BASE, "def", v), g(ZERO, "def", v)
        if b is not None:
            ax.plot(b, color=BLUE, lw=2, label="psi_lambda_B=3 (baseline)")
        if z is not None:
            ax.plot(z, color=RED, lw=2, ls="--", label="psi_lambda_B=0")
        ax.axhline(0, color="0.7", lw=0.6)
        ax.set_title(v, fontsize=11)
        ax.set_xlabel("quarters", fontsize=9)
    axes[0].legend(fontsize=8, loc="best")
    fig.suptitle("Decisive test — endogenous bond price & bank net worth on a "
                 "default shock", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, "02_decisive_qb_networth")


if __name__ == "__main__":
    chain_grid("01_default_transmission_chain", "def")
    decisive()
    chain_grid("03_tfp_control", "tfp")

    emitted = {f[:-4] for f in os.listdir(HERE) if f.endswith(".png")}
    missing_cap = emitted - set(CAPTIONS)
    missing_fig = set(CAPTIONS) - emitted
    assert not missing_cap and not missing_fig, (missing_cap, missing_fig)
    print("coverage check OK:", sorted(emitted))
