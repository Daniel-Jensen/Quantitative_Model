# IRF PANELS FOR THE TFP AND SOVEREIGN-RISK EXPERIMENTS, PLUS THE SPREAD DECOMPOSITION.
# Every figure is written to output/.
import os

import numpy as np
import matplotlib.pyplot as plt

from risk_branch import bond_decomposition
from prints import lending_spread_bps

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

_C_D, _C_F = "#1f77b4", "#d62728"
_C_P = "#7f2be8"   # real exchange rate


def _panel(ax, t, data_D, data_F, title, ylabel,
           label_D="D (domestic)", label_F="F (foreign)"):
    # TWO-LINE PANEL: COUNTRY D SOLID BLUE, F DASHED RED.
    ax.plot(t, data_D, color=_C_D, lw=1.5, label=label_D)
    if data_F is not None:
        ax.plot(t, data_F, color=_C_F, lw=1.5, ls="--", label=label_F)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title(title, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)
    ax.legend(fontsize=7)


def _single_panel(ax, t, data, title, ylabel, color=_C_P):
    # ONE-LINE PANEL FOR A VARIABLE WITH NO COUNTRY COUNTERPART.
    ax.plot(t, data, color=color, lw=1.5)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title(title, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)


def _save(fig, filename):
    # TIGHT-LAYOUT AND WRITE THE FIGURE TO output/.
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(os.path.join(OUTDIR, filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _scalers(T_plot):
    # PERCENT- AND BASIS-POINT-DEVIATION HELPERS TRUNCATED TO THE PLOT WINDOW.
    def pct(series, ref):
        # PERCENT DEVIATION FROM A REFERENCE LEVEL.
        return 100.0 * (np.asarray(series)[:T_plot] / ref - 1.0)

    def bps(series, ref=0.0):
        # BASIS-POINT DEVIATION FROM A REFERENCE RATE.
        return 10000.0 * (np.asarray(series)[:T_plot] - ref)

    return pct, bps


def plot_irf(out, ss, cal, T_plot=20, filename="tfp_irf.png"):
    # 5x4 IRF PANEL FOR THE TFP SHOCK (REAL, RATES, BANKING, EXTERNAL, BOND POSITIONS).
    t = np.arange(T_plot)
    pct, bps = _scalers(T_plot)
    fm_D = ss["ss_firm_D"];  fm_F = ss["ss_firm_F"]
    bk_D = ss["ss_bank_D"];  bk_F = ss["ss_bank_F"]

    # SS bond positions, naming b_{issuer}_{holder}: the home leg is the total
    # issued stock minus the cross-border leg
    b_D_F_ss = cal["b_D_F_ss"]
    b_F_D_ss = cal["b_F_D_ss"]
    b_D_D_ss = cal["B_gov_D_ss"] - b_D_F_ss
    b_F_F_ss = cal["B_gov_F_ss"] - b_F_D_ss

    fig, axes = plt.subplots(5, 4, figsize=(16, 16))
    fig.suptitle("IRF: TFP shock — country D (solid blue) vs F (dashed red)",
                 fontsize=11, y=1.005)

    _panel(axes[0, 0], t, pct(out["Y_D"], fm_D["Y_ss"]), pct(out["Y_F"], fm_F["Y_ss"]),
           "Output", "% dev. from SS")
    _panel(axes[0, 1], t, pct(out["C_D"], ss["C_D_ss"]), pct(out["C_F"], ss["C_F_ss"]),
           "Consumption", "% dev.")
    _panel(axes[0, 2], t, pct(out["I_D"], fm_D["I_ss"]), pct(out["I_F"], fm_F["I_ss"]),
           "Investment", "% dev.")
    _panel(axes[0, 3], t, pct(out["w_D"], fm_D["w_ss"]), pct(out["w_F"], fm_F["w_ss"]),
           "Real wage", "% dev.")

    _panel(axes[1, 0], t, bps(out["rk_D"], ss["rk_D_ss"]), bps(out["rk_F"], ss["rk_F_ss"]),
           "Capital return rk", "bps")
    _panel(axes[1, 1], t, bps(out["rdep_D"], cal["r_dep_D_target"]),
           bps(out["rdep_F"], cal["r_dep_F_target"]), "Deposit rate rdep", "bps")
    _panel(axes[1, 2], t, pct(out["Q_D"], 1.0), pct(out["Q_F"], 1.0),
           "Capital price Q", "% dev.")
    _panel(axes[1, 3], t,
           bps(np.asarray(out["rk_D"]) - np.asarray(out["rdep_D"])),
           bps(np.asarray(out["rk_F"]) - np.asarray(out["rdep_F"])),
           "Excess return rk - rdep", "bps")

    _panel(axes[2, 0], t, pct(out["n_D"], bk_D["n_ss"]), pct(out["n_F"], bk_F["n_ss"]),
           "Bank net worth n", "% dev.")
    _panel(axes[2, 1], t, pct(out["theta_D"], bk_D["theta_ss"]),
           pct(out["theta_F"], bk_F["theta_ss"]), "Bank leverage theta", "% dev.")
    _panel(axes[2, 2], t, pct(out["Q_bD"], ss["Q_bD_ss"]), pct(out["Q_bF"], ss["Q_bF_ss"]),
           "Sovereign bond price Q_b", "% dev.", label_D="Q_bD", label_F="Q_bF")
    _panel(axes[2, 3], t, pct(out["div_D"], bk_D["div_ss"]), pct(out["div_F"], bk_F["div_ss"]),
           "Bank dividends", "% dev.")

    _single_panel(axes[3, 0], t, pct(out["p"], ss["p_ss"]), "Real exch. rate p", "% dev.")
    _panel(axes[3, 1], t,
           100.0 * np.asarray(out["NX_D"])[:T_plot] / fm_D["Y_ss"],
           100.0 * np.asarray(out["NX_F"])[:T_plot] / fm_F["Y_ss"],
           "Net exports", "% of SS GDP")
    _panel(axes[3, 2], t, pct(out["b_gov_D"], cal["B_gov_D_ss"]),
           pct(out["b_gov_F"], cal["B_gov_F_ss"]), "Public debt b_gov", "% dev.")
    _panel(axes[3, 3], t, pct(out["A_D"], ss["A_D_ss"]), pct(out["A_F"], ss["A_F_ss"]),
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
    # 3x4 IRF PANEL FOR THE SOVEREIGN-RISK SHOCK (BOCOLA PASS-THROUGH).
    t = np.arange(T_plot)
    pct, bps = _scalers(T_plot)
    fm_D = ss["ss_firm_D"];  fm_F = ss["ss_firm_F"]
    bk_D = ss["ss_bank_D"];  bk_F = ss["ss_bank_F"]

    fig, axes = plt.subplots(3, 4, figsize=(16, 10))
    fig.suptitle("IRF: sovereign-risk shock (risk-only, priced never realized) "
                 "in country D — Bocola pass-through", fontsize=11, y=1.01)

    _panel(axes[0, 0], t, pct(out["Y_D"], fm_D["Y_ss"]), pct(out["Y_F"], fm_F["Y_ss"]),
           "Output", "% dev. from SS")
    _panel(axes[0, 1], t, pct(out["C_D"], ss["C_D_ss"]), pct(out["C_F"], ss["C_F_ss"]),
           "Consumption", "% dev.")
    _panel(axes[0, 2], t, pct(out["I_D"], fm_D["I_ss"]), pct(out["I_F"], fm_F["I_ss"]),
           "Investment", "% dev.")
    _single_panel(axes[0, 3], t, pct(out["p"], ss["p_ss"]), "Real exch. rate p", "% dev.")

    _panel(axes[1, 0], t, bps(out["rdep_D"], cal["r_dep_D_target"]),
           bps(out["rdep_F"], cal["r_dep_F_target"]), "Deposit rate rdep", "bps")
    _panel(axes[1, 1], t, bps(out["rk_D"], ss["rk_D_ss"]), bps(out["rk_F"], ss["rk_F_ss"]),
           "Capital return rk", "bps")
    _panel(axes[1, 2], t, pct(out["Q_bD"], ss["Q_bD_ss"]), pct(out["Q_bF"], ss["Q_bF_ss"]),
           "Bond price Q_b", "% dev.", label_D="Q_bD", label_F="Q_bF")
    _panel(axes[1, 3], t,
           bps(np.asarray(out["rk_D"]) - np.asarray(out["rdep_D"])),
           bps(np.asarray(out["rk_F"]) - np.asarray(out["rdep_F"])),
           "Excess return rk - rdep", "bps")

    _panel(axes[2, 0], t, pct(out["theta_D"], bk_D["theta_ss"]),
           pct(out["theta_F"], bk_F["theta_ss"]), "Bank leverage theta", "% dev.")
    _panel(axes[2, 1], t, pct(out["n_D"], bk_D["n_ss"]), pct(out["n_F"], bk_F["n_ss"]),
           "Bank net worth n", "% dev.")
    _panel(axes[2, 2], t,
           lending_spread_bps(out, ss, cal, "D")[:T_plot],
           lending_spread_bps(out, ss, cal, "F")[:T_plot],
           "Lending spread rk - rdep(-1)", "bps ann. dev.")

    ax = axes[2, 3]
    ax.plot(t, 100.0 * np.asarray(out["def_price_D"])[:T_plot],
            color=_C_D, lw=1.5, label="priced def. prob pi")
    ax.plot(t, pct(out["b_gov_D"], cal["B_gov_D_ss"]),
            color="#ff7f0e", lw=1.5, ls="--", label="b_gov_D (% dev)")
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title("Priced default risk & public debt", fontsize=9)
    ax.set_ylabel("%", fontsize=8)
    ax.set_xlabel("quarter", fontsize=8)
    ax.legend(fontsize=7)

    _save(fig, filename)


def plot_bond_decomposition(out, ss, cal, T_plot=100,
                            filename="bond_decomposition.png"):
    # SINGLE-PANEL D-BOND SPREAD DECOMPOSITION: DEFAULT COMP + RISK PREMIUM + LIQUIDITY.
    dec = bond_decomposition(out, ss, cal)
    t = np.arange(T_plot)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, dec["total_yield"][:T_plot], color="#2ca02c", lw=2.0,
            label="sovereign yield spread")
    ax.plot(t, dec["defcomp"][:T_plot], color=_C_D, lw=1.3, ls="--",
            label="default compensation")
    ax.plot(t, dec["risk"][:T_plot], color="#9467bd", lw=1.6, ls="-.",
            label="risk premium (Bocola channel)")
    ax.plot(t, dec["liquidity"][:T_plot], color=_C_F, lw=1.3, ls=":",
            label="liquidity premium lambda*mu/Omega")
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_title("Sovereign spread decomposition — D bonds", fontsize=10)
    ax.set_ylabel("bps ann., dev. from SS", fontsize=9)
    ax.set_xlabel("quarter", fontsize=9)
    ax.legend(fontsize=8)

    _save(fig, filename)
