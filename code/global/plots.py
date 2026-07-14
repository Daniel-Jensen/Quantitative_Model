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


def plot_default_irf(out, ss, cal, T_plot=100, filename="default_irf.png"):
    """3×4 panel IRF: Cole-Kehoe risk-only sunspot shock in country D
    (Bocola 2016 pass-through experiment).

    Each panel overlays country D (solid blue) and F (dashed red) so that
    cross-border spillovers through the GK financial channel are visible.
    Panels: real economy (Y, C, I, p), financial rates (rdep, rk, Q_b,
    excess return), balance sheets (θ, n), sovereign-spread decomposition
    (default compensation vs liquidity premium), and priced risk + debt.
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
    fig.suptitle("IRF: Cole-Kehoe sunspot (risk-only) in country D — "
                 "Bocola pass-through", fontsize=11, y=1.01)

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

    # Sovereign yield-spread decomposition (annualized bps), Bocola style:
    # total promised-yield spread; expected-return components: default
    # compensation, RISK premium (Bocola channel; 0 in risk-neutral mode)
    # and liquidity premium (IC component λμ/Ω̃).
    from risk_branch import bond_decomposition
    dec = bond_decomposition(out, ss, cal)
    ax = axes[2, 2]
    ax.plot(t, dec["total_yield"][:T_plot], color="#2ca02c", lw=1.8,
            label="sov yield spread")
    ax.plot(t, dec["defcomp"][:T_plot], color="#1f77b4", lw=1.2, ls="--",
            label="default compensation")
    ax.plot(t, dec["risk"][:T_plot], color="#9467bd", lw=1.4, ls="-.",
            label="risk premium (Bocola)")
    ax.plot(t, dec["liquidity"][:T_plot], color="#d62728", lw=1.2, ls=":",
            label="liquidity premium λμ/Ω̃")
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title("Sovereign spread decomposition (D)", fontsize=9)
    ax.set_ylabel("bps ann. dev.", fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)
    ax.legend(fontsize=7)

    c_D = "#1f77b4"
    ax = axes[2, 3]
    ax.plot(t, 100.0 * np.asarray(out["def_price_D"])[:T_plot],
            color=c_D, lw=1.5, label="priced def. prob (sunspot)")
    ax.plot(t, 100.0 * (np.asarray(out["b_gov_D"])[:T_plot] / cal["B_gov_D_ss"] - 1.0),
            color="#ff7f0e", lw=1.5, ls="--", label="b_gov_D (% dev)")
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title("Priced default risk & public debt", fontsize=9)
    ax.set_ylabel("%", fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)
    ax.legend(fontsize=7)

    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_risk_comparison(out_on, out_off, ss, cal, T_plot=40,
                         filename="risk_channel_irf.png"):
    """Risk channel on vs off (Bocola liquidity-only counterfactual):
    2×3 panels for country D — Y, I, C, bank net worth, lending spread,
    sovereign yield spread."""
    t = np.arange(T_plot)
    Y_ss, I_ss, C_ss = (ss["ss_firm_D"]["Y_ss"], ss["ss_firm_D"]["I_ss"],
                        ss["C_D_ss"])
    n_ss = ss["ss_bank_D"]["n_ss"]
    db = cal["delta_b_D"]

    def pct(x, ref):
        return 100.0 * (np.asarray(x)[:T_plot] / ref - 1.0)

    def lend(out):
        rdep_lag = np.concatenate([[cal["r_dep_D_target"]], out["rdep_D"][:-1]])
        return 1e4 * 4 * ((np.asarray(out["rk_D"]) - rdep_lag) - ss["rk_D_ss"])[:T_plot]

    def sov(out):
        y = db / np.asarray(out["Q_bD"]) - db
        return 4e4 * (y - (db / ss["Q_bD_ss"] - db))[:T_plot]

    panels = [
        ("Output D", "% dev.", pct(out_on["Y_D"], Y_ss), pct(out_off["Y_D"], Y_ss)),
        ("Investment D", "% dev.", pct(out_on["I_D"], I_ss), pct(out_off["I_D"], I_ss)),
        ("Consumption D", "% dev.", pct(out_on["C_D"], C_ss), pct(out_off["C_D"], C_ss)),
        ("Bank net worth D", "% dev.", pct(out_on["n_D"], n_ss), pct(out_off["n_D"], n_ss)),
        ("Lending spread rk − rdep", "bps ann. dev.", lend(out_on), lend(out_off)),
        ("Sovereign yield spread", "bps ann. dev.", sov(out_on), sov(out_off)),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    fig.suptitle("Bocola risk channel: two-branch pricing (on) vs "
                 "liquidity channel only (off)", fontsize=11, y=1.01)
    for ax, (title, ylab, on, off) in zip(axes.flat, panels):
        ax.plot(t, on, color="#9467bd", lw=1.8, label="risk channel ON")
        ax.plot(t, off, color="#7f7f7f", lw=1.4, ls="--", label="risk channel OFF")
        ax.axhline(0, color="k", lw=0.7, ls=":")
        ax.set_title(title, fontsize=9)
        ax.set_ylabel(ylab, fontsize=8)
        ax.set_xlabel("quarter", fontsize=8)
        ax.legend(fontsize=7)

    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_household_policies(ss, cal, filename="household_policies.png"):
    """2×2 panel: consumption policy, savings policy, MPC with quartile
    averages, and a bar-chart summary of average MPC by wealth quartile.

    Country D only (symmetric SS ⟹ F is identical).
    """
    a_grid = ss["a_grid_D"]
    c      = ss["c_D_ss"]
    a_pol  = ss["a_pol_D_ss"]
    D      = ss["D_D_ss"]
    e      = ss["e_D"]
    A_ss   = ss["A_D_ss"]
    a_lim  = min(cal["a_max_D"], 4 * A_ss + 1)

    # ── MPC: forward finite differences on consumption policy ─────────────────
    mpc   = np.clip(np.diff(c, axis=0) / np.diff(a_grid)[:, None], 0.0, 1.0)
    a_mid = 0.5 * (a_grid[:-1] + a_grid[1:])

    # ── Wealth quartile boundaries ────────────────────────────────────────────
    marg = D.sum(axis=1)
    cdf  = np.cumsum(marg)
    # q_bounds_idx: 5 index values [0, Q1, Q2, Q3, n_a]
    q_bounds_idx = ([0]
                    + [int(np.searchsorted(cdf, q)) for q in (0.25, 0.50, 0.75)]
                    + [len(a_grid) - 1])
    q_bounds_a = a_grid[q_bounds_idx]          # asset levels at each boundary

    # ── Quartile-average MPC (distribution-weighted within each quartile) ─────
    # weights[i, j] = D[i, j] aligned to mpc midpoint i  (use D[:-1])
    weights = D[:-1, :]
    mpc_q   = []
    for lo, hi in zip(q_bounds_idx[:-1], q_bounds_idx[1:]):
        hi_m = min(hi, len(mpc))
        w = weights[lo:hi_m, :]
        m = mpc[lo:hi_m, :]
        wsum = w.sum()
        mpc_q.append(float(np.sum(w * m) / wsum) if wsum > 0 else 0.0)

    # ── Overall distribution-weighted average MPC ─────────────────────────────
    mpc_avg = float(np.sum(weights * mpc) / weights.sum())

    # ── Style ─────────────────────────────────────────────────────────────────
    c_lo = "#1f77b4"
    c_hi = "#d62728"
    q_colors = ["#2ca02c", "#bcbd22", "#ff7f0e", "#9467bd"]
    q_labels = ["Q1\n(0–25%)", "Q2\n(25–50%)", "Q3\n(50–75%)", "Q4\n(75–100%)"]
    lbl_lo = f"e = {e[0]:.2f}  (low income)"
    lbl_hi = f"e = {e[1]:.2f}  (high income)"

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("Household policy functions — country D (Greece) steady state",
                 fontsize=10, y=1.01)

    # ── [0,0] Consumption policy ──────────────────────────────────────────────
    ax = axes[0, 0]
    ax.plot(a_grid, c[:, 0], color=c_lo, lw=1.5, label=lbl_lo)
    ax.plot(a_grid, c[:, 1], color=c_hi, lw=1.5, ls="--", label=lbl_hi)
    ax.set_xlim(0, a_lim)
    ax.set_title("Consumption policy  c(a, e)", fontsize=9)
    ax.set_xlabel("a  (wealth)", fontsize=8)
    ax.set_ylabel("c", fontsize=8)
    ax.legend(fontsize=7)

    # ── [0,1] Savings policy + 45° line ───────────────────────────────────────
    ax = axes[0, 1]
    ax.plot(a_grid, a_pol[:, 0], color=c_lo, lw=1.5, label=lbl_lo)
    ax.plot(a_grid, a_pol[:, 1], color=c_hi, lw=1.5, ls="--", label=lbl_hi)
    ax.plot(a_grid, a_grid, "k--", lw=0.8, label="45°  (a′ = a)")
    ax.set_xlim(0, a_lim); ax.set_ylim(0, a_lim)
    ax.set_title("Savings policy  a′(a, e)", fontsize=9)
    ax.set_xlabel("a  (wealth)", fontsize=8)
    ax.set_ylabel("a′", fontsize=8)
    ax.legend(fontsize=7)

    # ── [1,0] MPC curves + quartile-average steps + overall average ───────────
    ax = axes[1, 0]
    # Shaded quartile bands
    for i, (a_lo, a_hi, qc) in enumerate(
            zip(q_bounds_a[:-1], q_bounds_a[1:], q_colors)):
        ax.axvspan(a_lo, a_hi, alpha=0.07, color=qc)
    # MPC curves (thin, for reference)
    mask = a_mid <= a_lim
    ax.plot(a_mid[mask], mpc[mask, 0], color=c_lo, lw=1.0, alpha=0.6, label=lbl_lo)
    ax.plot(a_mid[mask], mpc[mask, 1], color=c_hi, lw=1.0, alpha=0.6,
            ls="--", label=lbl_hi)
    # Quartile-average MPC: bold horizontal step across each band
    for i, (a_lo, a_hi, qc, mq) in enumerate(
            zip(q_bounds_a[:-1], q_bounds_a[1:], q_colors, mpc_q)):
        ax.hlines(mq, a_lo, a_hi, colors=qc, lw=3.0, zorder=4)
        mid = (a_lo + a_hi) / 2
        va  = "bottom" if mq < 0.85 else "top"
        off = 0.025 if va == "bottom" else -0.025
        ax.text(mid, mq + off, f"Q{i+1}: {mq:.3f}",
                ha="center", va=va, fontsize=7.5,
                color=qc, fontweight="bold",
                bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))
    # Overall average MPC
    ax.axhline(mpc_avg, color="k", lw=1.2, ls=":",
               label=f"avg MPC = {mpc_avg:.3f}")
    ax.set_xlim(0, a_lim); ax.set_ylim(0, 1.08)
    ax.set_title("MPC(a, e)  —  distribution-weighted quartile averages", fontsize=9)
    ax.set_xlabel("a  (wealth)", fontsize=8)
    ax.set_ylabel("MPC", fontsize=8)
    ax.legend(fontsize=7, loc="upper right")

    # ── [1,1] Bar chart: average MPC by wealth quartile ───────────────────────
    ax = axes[1, 1]
    x = np.arange(4)
    bars = ax.bar(x, mpc_q, color=q_colors, width=0.6, zorder=3,
                  edgecolor="white", linewidth=0.5)
    # Value labels on bars
    for bar, mq in zip(bars, mpc_q):
        ax.text(bar.get_x() + bar.get_width() / 2,
                mq + 0.008, f"{mq:.3f}",
                ha="center", va="bottom", fontsize=8, fontweight="bold")
    # Overall average line
    ax.axhline(mpc_avg, color="k", lw=1.2, ls="--",
               label=f"avg MPC = {mpc_avg:.3f}", zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(q_labels, fontsize=8)
    ax.set_ylim(0, max(mpc_q) * 1.25)
    ax.set_title("Average MPC by wealth quartile", fontsize=9)
    ax.set_ylabel("MPC  (distribution-weighted mean)", fontsize=8)
    ax.legend(fontsize=7)
    ax.yaxis.grid(True, lw=0.5, alpha=0.5)
    ax.set_axisbelow(True)

    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_irf(out, ss, cal, T_plot=200):
    """4×4 panel IRF to a TFP shock, country D (solid blue) vs F (dashed red).

    Rows: (0) real economy — Y, C, I, wage; (1) financial rates — rk, rdep,
    capital price Q, excess return rk−rdep; (2) banking — net worth n,
    leverage θ, sovereign bond price Q_b, bank dividends; (3) external &
    fiscal — real exchange rate p, net exports, public debt, household
    deposits.
    """
    t = np.arange(T_plot)

    Y_D_ss  = ss["ss_firm_D"]["Y_ss"]
    Y_F_ss  = ss["ss_firm_F"]["Y_ss"]
    C_D_ss  = ss["C_D_ss"]
    C_F_ss  = ss["C_F_ss"]
    I_D_ss  = ss["ss_firm_D"]["I_ss"]
    I_F_ss  = ss["ss_firm_F"]["I_ss"]
    w_D_ss  = ss["ss_firm_D"]["w_ss"]
    w_F_ss  = ss["ss_firm_F"]["w_ss"]
    rk_D_ss = ss["rk_D_ss"]
    rk_F_ss = ss["rk_F_ss"]
    rdep_D_ss = cal["r_dep_D_target"]
    rdep_F_ss = cal["r_dep_F_target"]
    Q_D_ss  = ss["ss_cap_D"]["Q_ss"] if "ss_cap_D" in ss else 1.0
    Q_F_ss  = ss["ss_cap_F"]["Q_ss"] if "ss_cap_F" in ss else 1.0
    n_D_ss  = ss["ss_bank_D"]["n_ss"]
    n_F_ss  = ss["ss_bank_F"]["n_ss"]
    thetaD_ss = ss["ss_bank_D"]["theta_ss"]
    thetaF_ss = ss["ss_bank_F"]["theta_ss"]
    div_D_ss  = ss["ss_bank_D"]["div_ss"] if "div_ss" in ss["ss_bank_D"] else None
    div_F_ss  = ss["ss_bank_F"]["div_ss"] if "div_ss" in ss["ss_bank_F"] else None
    Q_bD_ss = ss["Q_bD_ss"]
    Q_bF_ss = ss["Q_bF_ss"]
    A_D_ss  = ss["A_D_ss"]
    A_F_ss  = ss["A_F_ss"]
    p_ss    = ss["p_ss"]

    # Steady-state sovereign-bond positions (naming: b_{issuer}_{holder}).
    # Home leg = total issued stock − cross-border leg (see CLAUDE.md).
    b_D_F_ss = cal["b_D_F_ss"]                       # D-bonds held by F banks
    b_F_D_ss = cal["b_F_D_ss"]                       # F-bonds held by D banks
    b_D_D_ss = cal["B_gov_D_ss"] - b_D_F_ss          # D-bonds held by D banks
    b_F_F_ss = cal["B_gov_F_ss"] - b_F_D_ss          # F-bonds held by F banks

    def pct(series, ref):
        return 100.0 * (np.asarray(series)[:T_plot] / ref - 1.0)

    def bps(series, ref=0.0):
        return 10000.0 * (np.asarray(series)[:T_plot] - ref)

    fig, axes = plt.subplots(5, 4, figsize=(16, 16))
    fig.suptitle("IRF: TFP shock — country D (solid blue) vs F (dashed red)",
                 fontsize=11, y=1.005)

    # ── Row 0: real economy ───────────────────────────────────────────────────
    _panel(axes[0, 0], t, pct(out["Y_D"], Y_D_ss), pct(out["Y_F"], Y_F_ss),
           "Output", "% dev. from SS")
    _panel(axes[0, 1], t, pct(out["C_D"], C_D_ss), pct(out["C_F"], C_F_ss),
           "Consumption", "% dev.")
    _panel(axes[0, 2], t, pct(out["I_D"], I_D_ss), pct(out["I_F"], I_F_ss),
           "Investment", "% dev.")
    _panel(axes[0, 3], t, pct(out["w_D"], w_D_ss), pct(out["w_F"], w_F_ss),
           "Real wage", "% dev.")

    # ── Row 1: financial rates ────────────────────────────────────────────────
    _panel(axes[1, 0], t, bps(out["rk_D"], rk_D_ss), bps(out["rk_F"], rk_F_ss),
           "Capital return rk", "bps")
    _panel(axes[1, 1], t, bps(out["rdep_D"], rdep_D_ss), bps(out["rdep_F"], rdep_F_ss),
           "Deposit rate rdep", "bps")
    _panel(axes[1, 2], t, pct(out["Q_D"], Q_D_ss), pct(out["Q_F"], Q_F_ss),
           "Capital price Q", "% dev.")
    _panel(axes[1, 3], t,
           bps(np.asarray(out["rk_D"]) - np.asarray(out["rdep_D"])),
           bps(np.asarray(out["rk_F"]) - np.asarray(out["rdep_F"])),
           "Excess return rk − rdep", "bps")

    # ── Row 2: banking block ──────────────────────────────────────────────────
    _panel(axes[2, 0], t, pct(out["n_D"], n_D_ss), pct(out["n_F"], n_F_ss),
           "Bank net worth n", "% dev.")
    _panel(axes[2, 1], t, pct(out["theta_D"], thetaD_ss), pct(out["theta_F"], thetaF_ss),
           "Bank leverage θ", "% dev.")
    _panel(axes[2, 2], t, pct(out["Q_bD"], Q_bD_ss), pct(out["Q_bF"], Q_bF_ss),
           "Sovereign bond price Q_b", "% dev.", label_D="Q_bD", label_F="Q_bF")
    if div_D_ss:
        _panel(axes[2, 3], t, pct(out["div_D"], div_D_ss), pct(out["div_F"], div_F_ss),
               "Bank dividends", "% dev.")
    else:
        _panel(axes[2, 3], t, bps(out["div_D"] - ss.get("div_D_ss", out["div_D"][0])),
               bps(out["div_F"] - ss.get("div_F_ss", out["div_F"][0])),
               "Bank dividends", "dev. (level)")

    # ── Row 3: external & fiscal ──────────────────────────────────────────────
    ax = axes[3, 0]
    ax.plot(t, pct(out["p"], p_ss), color="#7f2be8", lw=1.5)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title("Real exch. rate p", fontsize=9)
    ax.set_ylabel("% dev.", fontsize=8); ax.set_xlabel("quarter", fontsize=8)

    _panel(axes[3, 1], t,
           100.0 * np.asarray(out["NX_D"])[:T_plot] / Y_D_ss,
           100.0 * np.asarray(out["NX_F"])[:T_plot] / Y_F_ss,
           "Net exports", "% of SS GDP")
    _panel(axes[3, 2], t,
           pct(out["b_gov_D"], cal["B_gov_D_ss"]), pct(out["b_gov_F"], cal["B_gov_F_ss"]),
           "Public debt b_gov", "% dev.")
    _panel(axes[3, 3], t, pct(out["A_D"], A_D_ss), pct(out["A_F"], A_F_ss),
           "Household deposits A", "% dev.")

    # ── Row 4: sovereign-bond positions (cross-border) ────────────────────────
    _panel(axes[4, 0], t, pct(out["b_D_D"], b_D_D_ss), pct(out["b_D_F"], b_D_F_ss),
           "D sovereign bonds by holder", "% dev.",
           label_D="held by D banks (home)", label_F="held by F banks (cross-border)")
    _panel(axes[4, 1], t, pct(out["b_F_F"], b_F_F_ss), pct(out["b_F_D"], b_F_D_ss),
           "F sovereign bonds by holder", "% dev.",
           label_D="held by F banks (home)", label_F="held by D banks (cross-border)")
    _panel(axes[4, 2], t, pct(out["b_F_D"], b_F_D_ss), pct(out["b_D_F"], b_D_F_ss),
           "Cross-border bond holdings", "% dev.",
           label_D="D banks' F-bonds", label_F="F banks' D-bonds")
    _panel(axes[4, 3], t,
           pct(np.asarray(out["b_D_D"]) + np.asarray(out["b_F_D"]), b_D_D_ss + b_F_D_ss),
           pct(np.asarray(out["b_F_F"]) + np.asarray(out["b_D_F"]), b_F_F_ss + b_D_F_ss),
           "Total bank bond book", "% dev.",
           label_D="D banks (home+foreign)", label_F="F banks (home+foreign)")

    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, "tfp_irf.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
