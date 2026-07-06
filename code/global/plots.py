"""Plotting routines for the two-country HANK-GK monetary union model.

All figures are saved to the `output/` subdirectory next to this file.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")


def plot_steady_state(ss, cal):
    """2×2 panel: savings policies and CDF of deposit holdings, D and F."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    for col, country, lbl in [(0, "D", "D (domestic)"), (1, "F", "F (foreign)")]:
        a_grid   = ss[f"a_grid_{country}"]
        a_pol    = ss[f"a_pol_{country}_ss"]
        D        = ss[f"D_{country}_ss"]
        A_ss     = ss[f"A_{country}_ss"]
        a_lim    = min(cal[f"a_max_{country}"], 4 * A_ss + 1)

        # Savings policy
        ax = axes[0, col]
        ax.plot(a_grid, a_pol[:, 0], label="low income")
        ax.plot(a_grid, a_pol[:, 1], label="high income")
        ax.plot(a_grid, a_grid, "k--", lw=0.7)
        ax.set_xlim(0, a_lim); ax.set_ylim(0, a_lim)
        ax.set_xlabel("a (deposits)"); ax.set_ylabel("a'(a, e)")
        ax.set_title(f"Savings policy — country {lbl}")
        ax.legend(fontsize=8)

        # CDF
        ax = axes[1, col]
        marg = D.sum(axis=1)
        ax.plot(a_grid, np.cumsum(marg))
        ax.set_xlim(0, a_lim); ax.set_xlabel("a (deposits)")
        ax.set_title(f"CDF of deposits — country {lbl}")

    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, "steady_state.png"), dpi=150)
    plt.close(fig)


def _panel(ax, t, data_D, data_F, title, ylabel,
           label_D="D (domestic)", label_F="F (foreign)"):
    """Two-line panel: country D solid blue, country F dashed red."""
    c_D, c_F = "#1f77b4", "#d62728"
    ax.plot(t, data_D, color=c_D, lw=1.5, label=label_D)
    if data_F is not None:
        ax.plot(t, data_F, color=c_F, lw=1.5, ls="--", label=label_F)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title(title, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)
    ax.legend(fontsize=7)


def plot_default_irf(out, ss, cal, def_D_path, def_F_path,
                     T_plot=40, filename="default_irf.png"):
    """3×4 panel IRF: response to a default-risk shock in country D.

    Each panel overlays country D (solid blue) and F (dashed red) so that
    cross-border spillovers through the GK financial channel are visible.
    Panels: real economy (Y, C, I, p), financial rates (rdep, rk, Q_b,
    excess return), balance sheets (θ, n, bond spread, default rate).
    """
    t = np.arange(T_plot)

    Y_D_ss     = ss["ss_firm_D"]["Y_ss"]
    Y_F_ss     = ss["ss_firm_F"]["Y_ss"]
    C_D_ss     = ss["C_D_ss"]
    C_F_ss     = ss["C_F_ss"]
    I_D_ss     = ss["ss_firm_D"]["I_ss"]
    I_F_ss     = ss["ss_firm_F"]["I_ss"]
    rdep_D_ss  = cal["r_dep_D_target"]
    rdep_F_ss  = cal["r_dep_F_target"]
    rk_D_ss    = ss["rk_D_ss"]
    rk_F_ss    = ss["rk_F_ss"]
    theta_D_ss = ss["ss_bank_D"]["theta_ss"]
    theta_F_ss = ss["ss_bank_F"]["theta_ss"]
    n_D_ss     = ss["ss_bank_D"]["n_ss"]
    n_F_ss     = ss["ss_bank_F"]["n_ss"]
    Q_bD_ss    = ss["Q_bD_ss"]
    Q_bF_ss    = ss["Q_bF_ss"]
    p_ss       = ss["p_ss"]

    def pct(series, ref):
        return 100.0 * (np.asarray(series)[:T_plot] / ref - 1.0)

    def bps(series, ref=0.0):
        return 10000.0 * (np.asarray(series)[:T_plot] - ref)

    fig, axes = plt.subplots(3, 4, figsize=(16, 10))
    fig.suptitle("IRF: default-risk shock in country D (Greece analogue)",
                 fontsize=11, y=1.01)

    # ── Row 0: real economy ───────────────────────────────────────────────────
    _panel(axes[0, 0], t, pct(out["Y_D"], Y_D_ss), pct(out["Y_F"], Y_F_ss),
           "Output", "% dev. from SS")
    _panel(axes[0, 1], t, pct(out["C_D"], C_D_ss), pct(out["C_F"], C_F_ss),
           "Consumption", "% dev.")
    _panel(axes[0, 2], t, pct(out["I_D"], I_D_ss), pct(out["I_F"], I_F_ss),
           "Investment", "% dev.")

    ax = axes[0, 3]
    ax.plot(t, pct(out["p"], p_ss), color="#7f2be8", lw=1.5)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title("Real exch. rate p", fontsize=9)
    ax.set_ylabel("% dev.", fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)

    # ── Row 1: financial rates ────────────────────────────────────────────────
    _panel(axes[1, 0], t, bps(out["rdep_D"], rdep_D_ss), bps(out["rdep_F"], rdep_F_ss),
           "Deposit rate rdep", "bps")
    _panel(axes[1, 1], t, bps(out["rk_D"], rk_D_ss), bps(out["rk_F"], rk_F_ss),
           "Capital return rk", "bps")
    _panel(axes[1, 2], t,
           pct(out["Q_bD"], Q_bD_ss), pct(out["Q_bF"], Q_bF_ss),
           "Bond price Q_b", "% dev.", label_D="Q_bD", label_F="Q_bF")
    _panel(axes[1, 3], t,
           bps(np.asarray(out["rk_D"]) - np.asarray(out["rdep_D"])),
           bps(np.asarray(out["rk_F"]) - np.asarray(out["rdep_F"])),
           "Excess return rk − rdep", "bps")

    # ── Row 2: balance sheets and spreads ─────────────────────────────────────
    _panel(axes[2, 0], t, pct(out["theta_D"], theta_D_ss), pct(out["theta_F"], theta_F_ss),
           "Bank leverage θ", "% dev.")
    _panel(axes[2, 1], t, pct(out["n_D"], n_D_ss), pct(out["n_F"], n_F_ss),
           "Bank net worth n", "% dev.")

    ax = axes[2, 2]
    ax.plot(t, bps(np.asarray(out["Q_bD"]) - np.asarray(out["Q_bF"])),
            color="#2ca02c", lw=1.5)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title("Bond spread Q_bD − Q_bF", fontsize=9)
    ax.set_ylabel("bps", fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)

    c_D, c_F = "#1f77b4", "#d62728"
    ax = axes[2, 3]
    ax.plot(t, 100.0 * np.asarray(def_D_path)[:T_plot],
            color=c_D, lw=1.5, label="def D (shock)")
    ax.plot(t, 100.0 * np.asarray(def_F_path)[:T_plot],
            color=c_F, lw=1.5, ls="--", label="def F (none)")
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title("Default risk", fontsize=9)
    ax.set_ylabel("%", fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)
    ax.legend(fontsize=7)

    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_irf(out, ss, cal, T_plot=60):
    """3×4 panel IRF: D-country, F-country, financial, and global variables."""
    t = np.arange(T_plot)

    Y_D_ss  = ss["ss_firm_D"]["Y_ss"]
    Y_F_ss  = ss["ss_firm_F"]["Y_ss"]
    C_D_ss  = ss["C_D_ss"]
    C_F_ss  = ss["C_F_ss"]
    rk_D_ss = ss["rk_D_ss"]
    rk_F_ss = ss["rk_F_ss"]
    rdep_D_ss = cal["r_dep_D_target"]
    rdep_F_ss = cal["r_dep_F_target"]
    thetaD_ss = ss["ss_bank_D"]["theta_ss"]
    thetaF_ss = ss["ss_bank_F"]["theta_ss"]
    Q_bD_ss = ss["Q_bD_ss"]
    Q_bF_ss = ss["Q_bF_ss"]
    p_ss    = ss["p_ss"]

    def pct(series, ss_val):
        return 100.0 * (series[:T_plot] / ss_val - 1.0)

    def bps(series, ss_val):
        return 10000.0 * (series[:T_plot] - ss_val)

    fig, axes = plt.subplots(3, 4, figsize=(16, 10))

    panels = [
        (axes[0, 0], pct(out["Y_D"], Y_D_ss),    "Output D",       "% dev."),
        (axes[0, 1], pct(out["C_D"], C_D_ss),    "Consumption D",   "% dev."),
        (axes[0, 2], pct(out["Y_F"], Y_F_ss),    "Output F",       "% dev."),
        (axes[0, 3], pct(out["C_F"], C_F_ss),    "Consumption F",   "% dev."),

        (axes[1, 0], bps(out["rk_D"], rk_D_ss),   "rk D",          "bps"),
        (axes[1, 1], bps(out["rdep_D"], rdep_D_ss), "rdep D",       "bps"),
        (axes[1, 2], bps(out["rk_F"], rk_F_ss),   "rk F",          "bps"),
        (axes[1, 3], bps(out["rdep_F"], rdep_F_ss), "rdep F",       "bps"),

        (axes[2, 0], pct(out["theta_D"], thetaD_ss),  "Bank leverage D", "% dev."),
        (axes[2, 1], pct(out["theta_F"], thetaF_ss),  "Bank leverage F", "% dev."),
        (axes[2, 2], bps(out["Q_bD"] - out["Q_bF"], 0.0), "Spread Q_bD−Q_bF", "bps"),
        (axes[2, 3], pct(out["p"], p_ss),          "Real exch. rate p", "% dev."),
    ]

    for ax, data, title, ylabel in panels:
        ax.plot(t, data)
        ax.axhline(0.0, color="k", lw=0.7, ls="--")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("quarter", fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)

    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, "tfp_irf.png"), dpi=150)
    plt.close(fig)
