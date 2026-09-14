"""E1 — the core backstop schedule at the three named regimes.

Canonical parameterisation is the NAMED REGIMES (passive / medium / aggressive),
with gamma SOLVED for 0/25/50% peak-spread compression, not code/tpi.py's round
gamma in {0,2,5,10}. A solved gamma keeps its meaning across recalibrations; a
round number silently drifts into a different policy stance every time the model
is re-tuned.

A5-1: the German fiscal object is reported as THREE SEPARATE quantities and is
never summed into an "implicit transfer":

  1. exposure            — discounted purchases; what capital-key sharing acts on
  2. expected loss       — priced, computed OFF-PATH (see below)
  3. Greek fiscal saving — pd_D differential vs passive; no pricing assumption,
                           which makes it the cleanest headline of the three

The off-path requirement is not a stylistic choice. The excess-return flow
EL_price * def_rate_t * b_ss is a first-order deviation times a STEADY-STATE
level, so it does NOT vanish along the computed path: bondholders earn the premium
with no offsetting loss, because writeoff_enabled=0 and the IRF traces the
no-default branch. Reading expected loss off the realised path therefore shows the
CB mechanically profiting, which is an artifact of inserting a premium without the
compensating branch. It must be summed by hand as

    Sum_t beta^t * EL_price * def_rate_t * q_b * cb_buy_t

Welfare is computed but SECONDARY, and labelled so in the payload: docs/SPEC.md
says not to lead with it. It is decomposition-sensitive and comes out near-exactly
zero-sum, so it is the wrong thing for the paper to headline.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (BP_ANN, FIGURES_DIR, irf_from_cache, load_cache, pct_of_ss,
                    provenance, regime_irfs, write_results)

T_PNL = 100          # PV horizon, matches code/tpi.py's cb_pnl
T_WELFARE = 100      # discounted welfare horizon, matches run_tpi
PLOT_N = 60

# SS yields are equalised in this model, but only to solver tolerance (~1e-17 on
# the current solve). This guard must sit well above that and well below the
# 9.2e-04 that a wrong-duration bug produces — see cb_pnl.
SS_SPREAD_TOL = 1e-6

BLUE, RED, GREEN, ORANGE = "#002147", "#8C1515", "#1a6e3a", "#c87941"
REGIME_COLORS = {"passive": RED, "medium": ORANGE, "aggressive": GREEN}


def cb_pnl(irf, cache, T_pnl=T_PNL):
    """PV decomposition of the CB's D-bond position, in D-goods units.

    Ported from code/tpi.py's cb_pnl, reading steady-state levels from the cache
    rather than re-deriving them, so the two cannot drift apart.
    """
    beta_F = float(cache["beta_F"])
    q_b_D_ss = float(cache["q_b_D_ss"])
    q_b_F_ss = float(cache["q_b_F_ss"])
    delta_b_D = float(cache["delta_b_D_ss"])
    delta_b_F = float(cache["delta_b_F_ss"])
    EL_load_D = float(cache["EL_load_D"])

    disc = beta_F ** np.arange(T_pnl)
    cb = np.asarray(irf["cb_buy_D"])[:T_pnl]
    cb_l = np.concatenate([[0.0], cb[:-1]])
    dq = np.asarray(irf["q_b_D"])[:T_pnl]
    dq_l = np.concatenate([[0.0], dq[:-1]])
    defr = np.asarray(irf["def_rate_D"])[:T_pnl]
    dspr = np.asarray(irf["spread_rb"])[:T_pnl]

    # Each leg uses ITS OWN duration. delta_b_F (0.0568) != delta_b_D (0.0777) —
    # the two countries' bank books have different measured maturity ladders. An
    # earlier draft of this function used D's duration on both legs, which puts
    # the SS spread at -9.2e-04 instead of its true ~1e-17 and silently
    # contaminates carry_ss_pv. delta_b_F_ss was added to the cache (schema 3)
    # specifically so this could be done right.
    rb_D_ss = delta_b_D * (1.0 / q_b_D_ss - 1.0)
    rb_F_ss = delta_b_F * (1.0 / q_b_F_ss - 1.0)
    spread_ss = rb_D_ss - rb_F_ss
    assert abs(spread_ss) < SS_SPREAD_TOL, (
        f"SS spread {spread_ss:.3e} exceeds {SS_SPREAD_TOL:.0e}. SS yields are "
        f"equalised in this model, so this means a duration/price mismatch in the "
        f"carry legs — check delta_b_F_ss and q_b_F_ss came from the cache.")

    purchases = cb - (1.0 - delta_b_D) * cb_l
    return {
        "peak_exposure": float(np.max(q_b_D_ss * cb)),
        "purchases_pv": float((disc * q_b_D_ss * purchases).sum()),
        "el_pv": float((disc * EL_load_D * defr * q_b_D_ss * cb).sum()),
        "prem_pv": float((disc * dspr * q_b_D_ss * cb_l).sum()),
        "carry_ss_pv": float((disc * spread_ss * q_b_D_ss * cb_l).sum()),
        "mtm_pv": float((disc * (1.0 - delta_b_D) * cb_l * (dq - dq_l)).sum()),
    }


def primary_deficit(irf, cache):
    """pd_D = dG - P_CES_ss*dTAX - TAX_ss*dP_CES  (the austerity channel).

    Same construction as run_regimes.irf_all, kept identical on purpose. G_D is
    absent from the Jacobian (government spending is constant), so dG = 0.
    """
    T = len(np.asarray(irf["Y_D"]))
    dG = np.asarray(irf["G_D"]) if "G_D" in irf else np.zeros(T)
    return (dG - float(cache["P_CES_D_ss"]) * np.asarray(irf["TAX_D"])
            - float(cache["TAX_D_ss"]) * np.asarray(irf["P_CES_D"]))


def welfare(irf, cache, T_w=T_WELFARE):
    """Discounted utility deviation, % of quarterly SS consumption. SECONDARY."""
    beta_D, beta_F = float(cache["beta_D"]), float(cache["beta_F"])
    W_D = float((np.asarray(irf["U_D"])[:T_w] * beta_D ** np.arange(T_w) * 100).sum())
    W_F = float((np.asarray(irf["U_F"])[:T_w] * beta_F ** np.arange(T_w) * 100).sum())
    return W_D, W_F


def loading_schedule(cache, gamma_max=30.0, n=60):
    """Loading = premium PV / expected-loss PV over a fine gamma grid.

    THE key figure. The paper's self-extinguishing-premium claim is the DECLINE —
    the wedge exists because the marginal holder is balance-sheet constrained, and
    the backstop relieves that constraint, so intervention erodes its own profit
    source. The schedule, not any single point, is therefore the object.

    THE GRID STOPS BELOW THE CLOSED-LOOP POLE (2026-08-18). `(I - gamma*A_cb)` goes
    singular at `gamma ~ 27.3` on the post-GK-refactor cache, and the default
    `gamma_max = 30` ran straight through it: the loading spiked to 1.17 and collapsed
    to 0.38 across two grid points, the peak-spread panel showed a spurious dip to 82bp,
    and the caption's monotonicity test read the artefact as a real non-monotonicity.
    Everything past the pole is a DIFFERENT BRANCH of the closed loop, not a stronger
    version of the same policy, so it must not be plotted on the same axis.
    """
    from lottery_math import closed_loop, closed_loop_pole, POLE_SAFETY_FRACTION
    A_def, A_cb = cache["spread_rb__shock_def_D"], cache["spread_rb__cb_buy_D"]
    eps = np.asarray(cache["dShock_def_D"])
    pole = closed_loop_pole(A_cb, hi=max(gamma_max, 60.0))
    if pole is not None and POLE_SAFETY_FRACTION * pole <= gamma_max:
        gamma_max = POLE_SAFETY_FRACTION * pole
        print(f"  [loading_schedule] closed-loop pole at gamma = {pole:.2f}; "
              f"capping the grid at {gamma_max:.2f} "
              f"({POLE_SAFETY_FRACTION:g} x pole, lottery_math.POLE_SAFETY_FRACTION)")
    gammas = np.linspace(0.0, gamma_max, n)
    loading, peak_bp = np.full(n, np.nan), np.empty(n)
    for i, g in enumerate(gammas):
        spread, cb = closed_loop(A_def, A_cb, eps, float(g))
        irf = irf_from_cache(cache, cb, eps)
        d = cb_pnl(irf, cache)
        peak_bp[i] = float(np.max(spread[:T_PNL])) * BP_ANN
        if d["el_pv"] > 1e-16:
            loading[i] = d["prem_pv"] / d["el_pv"]
    return gammas, loading, peak_bp


def run():
    cache = load_cache()
    Y_ss = float(cache["Y_D_ss"])
    n_ss = float(cache["n_inter_D_ss"])
    K_ss = float(cache["K_D_ss"])
    C_ss = float(cache["C_D_ss"])
    I_ss = float(cache["I_D_ss"])
    beta_D = float(cache["beta_D"])

    regimes = regime_irfs(cache)
    payload = {
        "provenance": provenance(),
        "gamma_selection_rule":
            "gamma SOLVED for peak-spread compression, not chosen. medium = 25% (spec "
            "section 7). aggressive was 50%, but since the 2026-08-18 GK structural "
            "refactor that target lies beyond a closed-loop pole at gamma ~ 27.3 and is "
            "unreachable; it falls back to the strongest intervention the model can "
            "represent, gamma just below the pole, achieving ~46.6%. DO NOT describe the "
            "aggressive regime as 50% compression -- see common.named_regime_gammas and "
            "lottery_math.closed_loop_pole.",
        "welfare_caveat": "SECONDARY. SPEC: do not lead with welfare — it is a "
                          "delicate decomposition-dependent object and comes out "
                          "near-exactly zero-sum.",
        "regimes": {},
    }

    pd_passive = primary_deficit(regimes["passive"][1], cache)
    disc = beta_D ** np.arange(T_PNL)

    for name, (gamma, irf) in regimes.items():
        pnl = cb_pnl(irf, cache)
        W_D, W_F = welfare(irf, cache)
        spread = np.asarray(irf["spread_rb"])
        pd_here = primary_deficit(irf, cache)
        # A5-1 object 3: Greek fiscal saving vs the passive counterfactual.
        fiscal_saving_pv = float((disc * (pd_passive - pd_here)[:T_PNL]).sum())

        payload["regimes"][name] = {
            "gamma": gamma,
            "peak_spread_bp_ann": float(np.max(spread[:T_PNL]) * BP_ANN),
            # Every percentage divides by its OWN steady-state level. SSJ IRFs are
            # LEVEL deviations, so x100 is a percentage only where the SS level is
            # ~1. Y_D_ss~1 passes; n_inter_D_ss=2.138 and K_D_ss=10.8 do NOT, and
            # a past bug mislabelled exactly these two by 2.1x and 10x.
            "impact": {
                "Y_D_pct_ss": float(pct_of_ss(np.asarray(irf["Y_D"])[:1], Y_ss)[0]),
                "C_D_pct_ss": float(pct_of_ss(np.asarray(irf["C_D"])[:1], C_ss)[0]),
                "I_D_pct_ss": float(pct_of_ss(np.asarray(irf["I_D"])[:1], I_ss)[0]),
                "n_inter_D_pct_ss": float(pct_of_ss(np.asarray(irf["n_inter_D"])[:1], n_ss)[0]),
                "K_D_pct_ss": float(pct_of_ss(np.asarray(irf["K_D"])[:1], K_ss)[0]),
            },
            "trough": {
                "Y_D_pct_ss": float(np.asarray(irf["Y_D"])[:T_PNL].min() * 100.0 / Y_ss),
                "n_inter_D_pct_ss": float(np.asarray(irf["n_inter_D"])[:T_PNL].min() * 100.0 / n_ss),
            },
            # ---- A5-1: three separate objects, never summed ----
            "a5_1_exposure_pv_pct_Y": 100.0 * pnl["purchases_pv"] / Y_ss,
            "a5_1_expected_loss_pv_pct_Y": 100.0 * pnl["el_pv"] / Y_ss,
            "a5_1_greek_fiscal_saving_pv": fiscal_saving_pv,
            "peak_exposure_pct_Y": 100.0 * pnl["peak_exposure"] / Y_ss,
            "premium_pv_pct_Y": 100.0 * pnl["prem_pv"] / Y_ss,
            "mtm_pv_pct_Y": 100.0 * pnl["mtm_pv"] / Y_ss,
            "carry_ss_pv_pct_Y": 100.0 * pnl["carry_ss_pv"] / Y_ss,
            # None, not NaN: write_results uses allow_nan=False, and None is the
            # deliberate "not applicable" encoding — the passive regime buys
            # nothing, so the expected-loss denominator is zero by construction.
            "loading": (pnl["prem_pv"] / pnl["el_pv"]) if pnl["el_pv"] > 1e-16 else None,
            "welfare_W_D_secondary": W_D,
            "welfare_W_F_secondary": W_F,
        }

    assert payload["regimes"]["passive"]["a5_1_exposure_pv_pct_Y"] == 0.0, \
        "passive regime must have zero purchases by construction"

    gammas, loading, peak_bp = loading_schedule(cache)
    payload["loading_schedule"] = {
        "gamma": gammas.tolist(),
        "loading": [None if np.isnan(x) else float(x) for x in loading],
        "peak_spread_bp_ann": peak_bp.tolist(),
    }

    _plot(payload, gammas, loading, peak_bp)
    write_results("e1_backstop_schedule", payload)
    return payload


def _plot(payload, gammas, loading, peak_bp):
    os.makedirs(FIGURES_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    axes[0].plot(gammas, loading, color=BLUE, lw=2)
    for name, r in payload["regimes"].items():
        if r["loading"] is not None:
            axes[0].plot(r["gamma"], r["loading"], "o", color=REGIME_COLORS[name],
                         ms=8, label=f"{name} (γ={r['gamma']:.2f})")
    axes[0].axhline(1.0, ls="--", lw=0.8, color="gray")
    axes[0].set_xlabel("backstop aggressiveness γ")
    axes[0].set_ylabel("loading = premium PV / expected-loss PV")
    axes[0].set_title("The premium self-extinguishes as the backstop strengthens")
    axes[0].legend(fontsize=8)

    axes[1].plot(gammas, peak_bp, color=BLUE, lw=2)
    for name, r in payload["regimes"].items():
        axes[1].plot(r["gamma"], r["peak_spread_bp_ann"], "o",
                     color=REGIME_COLORS[name], ms=8, label=name)
    axes[1].set_xlabel("backstop aggressiveness γ")
    axes[1].set_ylabel("peak D–F spread (bp, annualised)")
    axes[1].set_title("Spread compression")
    axes[1].legend(fontsize=8)

    # Provenance read LIVE, never hardcoded: run_regimes.py once shipped a
    # hardcoded "market-value rule" suptitle while running at mv_rule=0.
    p = payload["provenance"]
    fig.suptitle(f"E1 — backstop schedule (psi_lambda_B={p['psi_lambda_B_D']}, "
                 f"{'market-value' if p['mv_rule_D'] else 'par-value'} rule, "
                 f"writeoff={'on' if p['writeoff_enabled_D'] else 'off'}, "
                 f"scope={p['BANK_SCOPE']}, {p['git_sha']})", fontsize=10)
    fig.text(0.5, 0.01,
             "Loading is premium PV over expected-loss PV, with expected loss computed "
             "OFF-PATH. The decline in γ is the self-extinguishing-premium result: the "
             "wedge exists because the marginal holder is balance-sheet constrained, and "
             "the backstop relieves that constraint, so intervention erodes its own "
             "profit source.", ha="center", fontsize=7.5, style="italic", wrap=True)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(os.path.join(FIGURES_DIR, "fig_e1_loading_schedule.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    res = run()
    print(f"{'regime':>12} {'gamma':>9} {'peak bp':>9} {'Y[0] %SS':>10} "
          f"{'n_int[0] %SS':>13} {'exposure %Y':>12} {'EL PV %Y':>10} {'loading':>8}")
    print("-" * 92)
    for name, r in res["regimes"].items():
        ld = "n/a" if r["loading"] is None else f"{r['loading']:.2f}"
        print(f"{name:>12} {r['gamma']:>9.4f} {r['peak_spread_bp_ann']:>9.1f} "
              f"{r['impact']['Y_D_pct_ss']:>+10.4f} {r['impact']['n_inter_D_pct_ss']:>+13.3f} "
              f"{r['a5_1_exposure_pv_pct_Y']:>12.3f} "
              f"{r['a5_1_expected_loss_pv_pct_Y']:>10.4f} {ld:>8}")
    print("-" * 92)
    print("A5-1: exposure / expected loss / Greek fiscal saving are SEPARATE objects — "
          "never sum them into an 'implicit transfer'.")
    print("Welfare is reported in the JSON as SECONDARY only (SPEC: do not lead with it).")
