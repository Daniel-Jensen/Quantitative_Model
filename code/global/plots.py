# **IRF panels for the TFP and sunspot experiments + the spread decomposition (saved to output/).**
import os
import numpy as np
import matplotlib.pyplot as plt

from risk_branch import bond_decomposition

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")


def _panel(ax, t, data_D, data_F, title, ylabel,
           label_D="D (domestic)", label_F="F (foreign)"):
    # **Two-line panel: country D solid blue, F dashed red.**
    c_D, c_F = "#1f77b4", "#d62728"
    ax.plot(t, data_D, color=c_D, lw=1.5, label=label_D)
    if data_F is not None:
        ax.plot(t, data_F, color=c_F, lw=1.5, ls="--", label=label_F)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title(title, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)
    ax.legend(fontsize=7)


def _save(fig, filename):
    # **tight-layout and write the figure to output/.**
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_irf(out, ss, cal, T_plot=20, filename="tfp_irf.png"):
    # **5×4 IRF panel for the TFP shock (real economy, rates, banking, external, bond positions).**
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
    Q_D_ss  = 1.0   # SS capital price
    Q_F_ss  = 1.0
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

    _panel(axes[0, 0], t, pct(out["Y_D"], Y_D_ss), pct(out["Y_F"], Y_F_ss),
           "Output", "% dev. from SS")
    _panel(axes[0, 1], t, pct(out["C_D"], C_D_ss), pct(out["C_F"], C_F_ss),
           "Consumption", "% dev.")
    _panel(axes[0, 2], t, pct(out["I_D"], I_D_ss), pct(out["I_F"], I_F_ss),
           "Investment", "% dev.")
    _panel(axes[0, 3], t, pct(out["w_D"], w_D_ss), pct(out["w_F"], w_F_ss),
           "Real wage", "% dev.")

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

    _save(fig, filename)


def plot_default_irf(out, ss, cal, T_plot=100, filename="default_irf.png"):
    # **3×4 IRF panel for the Cole-Kehoe risk-only sunspot (Bocola pass-through).**
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

    _panel(axes[2, 0], t, pct(out["theta_D"], theta_D_ss), pct(out["theta_F"], theta_F_ss),
           "Bank leverage θ", "% dev.")
    _panel(axes[2, 1], t, pct(out["n_D"], n_D_ss), pct(out["n_F"], n_F_ss),
           "Bank net worth n", "% dev.")

    # Lending spread rk_t − rdep_{t−1} (the deposit rate paid at t was locked
    # at t−1), annualized bps deviation from SS — the pass-through headline.
    rdep_lag_D = np.concatenate([[rdep_D_ss], np.asarray(out["rdep_D"])[:-1]])
    rdep_lag_F = np.concatenate([[rdep_F_ss], np.asarray(out["rdep_F"])[:-1]])
    lend_D = 4e4 * ((np.asarray(out["rk_D"]) - rdep_lag_D) - (rk_D_ss - rdep_D_ss))
    lend_F = 4e4 * ((np.asarray(out["rk_F"]) - rdep_lag_F) - (rk_F_ss - rdep_F_ss))
    _panel(axes[2, 2], t, lend_D[:T_plot], lend_F[:T_plot],
           "Lending spread rk − rdep(−1)", "bps ann. dev.")

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

    _save(fig, filename)


def plot_bond_decomposition(out, ss, cal, T_plot=100,
                            filename="bond_decomposition.png"):
    # **Single-panel D-bond spread decomposition: default comp + risk premium + liquidity.**
    dec = bond_decomposition(out, ss, cal)
    t = np.arange(T_plot)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, dec["total_yield"][:T_plot], color="#2ca02c", lw=2.0,
            label="sovereign yield spread")
    ax.plot(t, dec["defcomp"][:T_plot], color="#1f77b4", lw=1.3, ls="--",
            label="default compensation")
    ax.plot(t, dec["risk"][:T_plot], color="#9467bd", lw=1.6, ls="-.",
            label="risk premium (Bocola channel)")
    ax.plot(t, dec["liquidity"][:T_plot], color="#d62728", lw=1.3, ls=":",
            label="liquidity premium λμ/Ω̃")
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title("Sovereign spread decomposition — D bonds", fontsize=10)
    ax.set_ylabel("bps ann., dev. from SS", fontsize=9)
    ax.set_xlabel("quarter", fontsize=9)
    ax.legend(fontsize=8)

    _save(fig, filename)
