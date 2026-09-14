# D-vs-F SOVEREIGN BOND SPREAD UNDER THE RISK SHOCK, ACROSS OMT/TPI ACTIVATION.
# Answers "why does the no-TPI bond fall ~1% while the standard exercise falls ~8%?":
# the no-TPI (a=0) case IS the standard risk pass-through (a=0 nests the two-branch),
# so they are identical -- the ~8% is the OLD perfect-foresight liquidity channel,
# which the mu=1 recursive solver understates (documented indicative magnitudes).
# Extracts the D-bond price Q_bD, the (safe) F-bond price Q_bF, and the D-F YIELD
# SPREAD y_D - y_F (the sovereign risk premium of D over F, annualised bps), and
# plots: (1) the baseline (a=0) bond prices + spread; (2) Q_bD across 10 activation
# probabilities; (3) the D-F spread across those 10 probabilities.
import numpy as np

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import s_process_params, default_prob, IS
from solver_recursive.recursive_experiment import s_from_pd
from solver_recursive.recursive_main import ss_state, calibrate_household_anchors
from solver_recursive.recursive_experiment import solve_recursive, read_at

ACTIVATIONS = np.round(np.arange(0.0, 0.95, 0.1), 2)   # phi = 0,10,...,90 %
# The shock is a TARGET one-quarter-ahead default probability (main.py's
# RISK_SHOCK_PD), not a hard-wired s. The old constant -3.9 was labelled
# "+2 sigma": at the calibrated sigma_s = 0.63 it is +4.77 sigma, and a genuine
# +2 sigma shock is p^d = 0.35%, not 2%.
PD_SHOCK, T_IRF = 0.0198, 21


def bond_irf(rules, cal, ss, sproc):
    # Q_bD, Q_bF (% dev) AND the D-F YIELD SPREAD (abs, bps ann) ALONG THE s-DECAY.
    dbD, dbF = cal["delta_b_D"], cal["delta_b_F"]
    S0 = ss_state(ss, cal, sproc)
    base = read_at(rules, cal, ss, sproc, S0.copy())
    QbD0, QbF0 = base["Q_bD"], base["Q_bF"]
    P = {k: np.empty(T_IRF) for k in ("pd", "QbD_pct", "QbF_pct", "spread_bp")}
    for t in range(T_IRF):
        s_t = (sproc["s_star"]
               + sproc["rho_s"] ** t * (s_from_pd(PD_SHOCK) - sproc["s_star"]))
        S = S0.copy(); S[IS] = s_t
        o = read_at(rules, cal, ss, sproc, S)
        QbD, QbF = o["Q_bD"], o["Q_bF"]
        yD = dbD * (1.0 - QbD) / QbD                    # flow yield = coupon(1-Q)/Q
        yF = dbF * (1.0 - QbF) / QbF
        P["pd"][t] = 100 * default_prob(s_t)
        P["QbD_pct"][t] = 100 * (QbD / QbD0 - 1)
        P["QbF_pct"][t] = 100 * (QbF / QbF0 - 1)
        P["spread_bp"][t] = 4e4 * (yD - yF)             # D-F sovereign spread, bps ann
    return P


def main():
    import os, time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from reporting.plots import OUTDIR
    os.makedirs(OUTDIR, exist_ok=True)

    cal = get_calibration()
    cal["nw_floor_frac"] = 0.15      # match main.py: without it this is a different model
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print(f"=== D-F bond spread across OMT/TPI activation (SS Q_bD=Q_bF={ss['Q_bD_ss']:.3f}) ===",
          flush=True)

    t0 = time.perf_counter()
    irfs = {}
    for a in ACTIVATIONS:
        cal["phi_ltro"] = float(a)
        # coarse grid: this is a parameter sweep, not a headline result (see
        # calibrate_stochastic.py for the same reasoning)
        rules = solve_recursive(cal, ss, sproc, mu=1, verbose=False, s_refine=0,
                                with_cb=(a > 0.0))
        irfs[a] = bond_irf(rules, cal, ss, sproc)
        P = irfs[a]
        print(f"  a={a:.1f}: impact Q_bD={P['QbD_pct'][0]:+.2f}%  Q_bF={P['QbF_pct'][0]:+.2f}%"
              f"  D-F spread={P['spread_bp'][0]:+.0f}bp  ({time.perf_counter()-t0:.0f}s)",
              flush=True)

    t = np.arange(T_IRF)
    base = irfs[ACTIVATIONS[0]]                          # a=0 == standard risk exercise

    # FIG 1: baseline (no TPI = standard exercise) -- Q_bD, Q_bF, and the D-F spread
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Baseline sovereign-risk shock (no TPI = standard exercise): "
                 "D-bond, F-bond, and their spread", fontsize=12)
    ax[0].plot(t, base["QbD_pct"], color="#d62728", lw=2, label="D bond (risky)")
    ax[0].plot(t, base["QbF_pct"], color="#1f77b4", lw=2, ls="--", label="F bond (safe)")
    ax[0].axhline(0, color="k", lw=0.7, ls=":"); ax[0].legend(fontsize=9)
    ax[0].set_title("Bond price (% dev from low-risk)"); ax[0].set_xlabel("quarter")
    ax[0].set_ylabel("% dev")
    ax[1].plot(t, base["spread_bp"], color="#2ca02c", lw=2)
    ax[1].axhline(0, color="k", lw=0.7, ls=":")
    ax[1].set_title("D-F yield spread (sovereign risk premium)")
    ax[1].set_xlabel("quarter"); ax[1].set_ylabel("bps ann.")
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "bond_spread_baseline.png"),
                                    dpi=150, bbox_inches="tight"); plt.close(fig)

    # FIG 2 & 3: Q_bD and the D-F spread across the 10 activation probabilities
    cmap = plt.cm.viridis(np.linspace(0, 0.92, len(ACTIVATIONS)))
    for key, ylab, title, fname in (
            ("QbD_pct", "% dev", "D-bond price Q_bD under the risk shock, by OMT/TPI activation",
             "bond_spread_tpi_qbd.png"),
            ("spread_bp", "bps ann.", "D-F sovereign spread under the risk shock, by OMT/TPI activation",
             "bond_spread_tpi_spread.png")):
        fig, ax = plt.subplots(figsize=(9, 6))
        for a, c in zip(ACTIVATIONS, cmap):
            ax.plot(t, irfs[a][key], color=c, lw=1.8, label=f"{int(a*100)}%")
        ax.axhline(0, color="k", lw=0.7, ls=":")
        ax.set_title(title, fontsize=11); ax.set_xlabel("quarter"); ax.set_ylabel(ylab)
        ax.legend(title="TPI activation", fontsize=8, ncol=2)
        fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, fname), dpi=150,
                                        bbox_inches="tight"); plt.close(fig)

    print(f"\n  a=0 IS the standard risk exercise (no TPI). Impact Q_bD fall = "
          f"{base['QbD_pct'][0]:+.2f}% at mu=1 (indicative; the ~8% you recall is the old "
          f"perfect-foresight liquidity channel).", flush=True)
    print(f"  figures -> {OUTDIR}/bond_spread_baseline.png, bond_spread_tpi_qbd.png, "
          f"bond_spread_tpi_spread.png", flush=True)


if __name__ == "__main__":
    main()
