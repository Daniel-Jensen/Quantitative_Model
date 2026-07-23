"""
run_regimes.py — Stage A: three deterministic backstop regimes (MAIN model).

Regime gammas: passive = 0; aggressive = ~50% peak-SPREAD compression;
medium = ~25% (spec §7, computed not chosen). On main's model the ECB capital-key
conduit (kappa_cb_F=0.929, funding socialised to F) plus the market-value fiscal
rule make CB purchases COMPRESS the D-F spread (A_cb < 0) — the ms-regime SA-1
pathology (purchases widen the spread) is ABSENT here, so spread-compression
targeting is feasible and is used, as the spec originally intended.

Outputs: three-regime IRF set (§8.1), A5 implicit-transfer schedule, A6
psi_lambda_B=0 amplifier-invariance table, A7 blow-up guard. Pure numpy on
regime_model caches; writes regimes_calibration.json.
"""
import os, sys, json, datetime
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from regime_model import build_caches, load_cache, PSILAM_MAIN
from lottery_math import closed_loop, gamma_for_compression, peak

LOG = os.path.join(HERE, "regimes_log.md")
BLUE, RED, GREEN, ORANGE = "#002147", "#8C1515", "#1a6e3a", "#c87941"
BP_ANN = 4.0 * 1e4   # quarterly rate deviation -> annualised basis points
N = 60               # plot horizon


def log(m=""):
    print(m, flush=True)
    with open(LOG, "a") as f:
        f.write(m + "\n")


def irf_all(cache, cb_path, eps, T):
    """Apply cached matrices: output deviation = M_def@eps + M_cb@cb."""
    outs = sorted({k.split("__")[0] for k in cache if k.endswith("__cb_buy_D")})
    d = {}
    for o in outs:
        d[o] = cache[f"{o}__shock_def_D"] @ eps + cache[f"{o}__cb_buy_D"] @ cb_path
    d["cb_buy_D"] = cb_path
    dG = d["G_D"] if "G_D" in d else np.zeros(T)
    d["pd_D"] = dG - float(cache["P_CES_D_ss"]) * d["TAX_D"] - float(cache["TAX_D_ss"]) * d["P_CES_D"]
    d["mv_b_D_D"] = float(cache["q_b_D_ss"]) * d["b_D_D"] + float(cache["b_D_D_ss"]) * d["q_b_D"]
    return d


def main():
    open(LOG, "a").write(f"\n\n## Stage A run (MAIN model) — {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n")
    build_caches()
    cM, c0 = load_cache(PSILAM_MAIN), load_cache(0.0)
    T = int(cM["T"]); eps = np.asarray(cM["dShock_def_D"])
    A_def, A_cb = cM["spread_rb__shock_def_D"], cM["spread_rb__cb_buy_D"]

    # ── Regime calibration (spec §7): 50%/25% peak-spread compression ─────────
    sp_passive, _ = closed_loop(A_def, A_cb, eps, 0.0)
    pk_passive_bp = peak(sp_passive) * BP_ANN
    log(f"- passive peak spread: {pk_passive_bp:.1f} bp ann "
        f"(main is calibrated to ~151bp at psi_lambda_B={PSILAM_MAIN}; investigate if outside 120-180)")
    if pk_passive_bp > 1000:
        log("  **A7 FLAG: passive peak spread > 1000bp — LINEAR-APPROXIMATION BREAKDOWN, "
            "not an economic finding.**")
    log(f"- impact A_cb = d(spread)/d(cb_buy)[0,0] = {A_cb[0,0]:+.3e} "
        f"({'COMPRESSES (backstop works — SA-1 absent on main)' if A_cb[0,0] < 0 else 'WIDENS — unexpected on main, investigate'})")
    g_A = gamma_for_compression(A_def, A_cb, eps, target=0.50)
    g_M = gamma_for_compression(A_def, A_cb, eps, target=0.25)
    gammas = {"aggressive": g_A, "medium": g_M, "passive": 0.0}
    log(f"- gamma_aggressive = {g_A:.4f} (50% peak-spread compression), "
        f"gamma_medium = {g_M:.4f} (25%), gamma_passive = 0 (fixed anchor)")
    json.dump({"gammas": gammas, "pk_passive_bp": pk_passive_bp, "A_cb_impact": float(A_cb[0, 0]),
               "selection_rule": "peak-spread compression 50%/25% (spec §7; feasible on main, "
                                 "A_cb<0 — capital-key conduit; ms-regime SA-1 pathology absent)",
               "model": "main (psi_lambda_B=%.4f, mv_rule=1, capital-key conduit)" % PSILAM_MAIN,
               "generated": str(datetime.datetime.now())},
              open(os.path.join(HERE, "regimes_calibration.json"), "w"), indent=2)

    # ── Closed-loop IRFs per regime ──────────────────────────────────────────
    irfs = {}
    for name, g in gammas.items():
        sp, cb = closed_loop(A_def, A_cb, eps, g)
        irfs[name] = irf_all(cM, cb, eps, T)
        log(f"  {name:>10}: peak spread {peak(sp)*BP_ANN:7.1f} bp, "
            f"Y_D[0] {irfs[name]['Y_D'][0]*100:+.4f}%, "
            f"n_inter_D[0] {irfs[name]['n_inter_D'][0]*100:+.3f}%")

    # ── §8.1 IRF figure ──────────────────────────────────────────────────────
    panels = [("spread_rb", "D-F yield spread (bp ann)", BP_ANN),
              ("q_b_D", "q_b_D", 1.0), ("Y_D", "Y_D", 100.0), ("I_D", "I_D", 100.0),
              ("C_D", "C_D", 100.0), ("NX_D", "NX_D", 100.0),
              ("n_inter_D", "bank net worth n_inter_D (%)", 100.0),
              ("K_D", "capital level K_D (%)", 100.0),
              ("b_D_D", "b_D_D quantity", 1.0), ("mv_b_D_D", "b_D_D market value", 1.0),
              ("pd_D", "primary deficit (austerity channel)", 1.0),
              ("cb_buy_D", "CB purchases", 1.0)]
    colors = {"aggressive": GREEN, "medium": ORANGE, "passive": RED}
    fig, axes = plt.subplots(3, 4, figsize=(16, 10))
    for ax, (v, title, scale) in zip(axes.flat, panels):
        for name in gammas:
            ax.plot(np.arange(N), irfs[name][v][:N] * scale, color=colors[name], label=name)
        ax.set_title(title, fontsize=9); ax.axhline(0, lw=0.5, color="gray")
    axes.flat[0].legend(fontsize=8)
    fig.suptitle(f"Stage A (main): 1pp default shock under three exogenous backstop regimes "
                 f"(psi_lambda_B={PSILAM_MAIN}, market-value rule, capital-key conduit)")
    fig.text(0.5, 0.02,
             "Three exogenous ECB backstop regimes (deviation-form rule TPI=gamma*(spread-ss)): "
             "more aggressive regimes COMPRESS the D-F spread and cushion investment, output and "
             "bank net worth. On main the capital-key conduit socialises the purchase funding to "
             "the core (F), so a periphery-bond purchase relieves periphery banks and narrows the "
             "spread as intended — the ms-regime SA-1 spread-widening pathology is absent. "
             "G_D absent from the Jacobian (constant), so the primary-deficit panel uses dG=0.",
             ha="center", fontsize=7.5, style="italic", wrap=True)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(HERE, "fig_stageA_irfs.png"), dpi=200)

    # ── A5: implicit-transfer schedule ───────────────────────────────────────
    beta = float(cM["beta_D"]); disc = beta ** np.arange(100)
    qb = float(cM["q_b_D_ss"])
    log("\n| regime | discounted CB purchases (Sum beta^t q_b cb_t) | dY_D peak (%) | dC_D peak (%) | pd_D peak |")
    log("|---|---|---|---|---|")
    rows = {}
    for name in gammas:
        buy = (irfs[name]["cb_buy_D"][:100] * qb * disc).sum()
        rows[name] = buy
        log(f"| {name} | {buy:.5f} | {irfs[name]['Y_D'][:100].min()*100:+.4f} | "
            f"{irfs[name]['C_D'][:100].min()*100:+.4f} | {irfs[name]['pd_D'][:100].max():+.5f} |")
    assert rows["passive"] == 0.0, "A5: passive regime must have zero purchases by construction"

    # ── A6: psi_lambda_B=0 amplifier-invariance (crisis-severity ranking) ─────
    A_def0, A_cb0 = c0["spread_rb__shock_def_D"], c0["spread_rb__cb_buy_D"]
    log("\n### A6 — ranking at psi_lambda_B = 0 (fundamental floor)")
    sev = {}
    for name, g in gammas.items():
        sp0, cb0 = closed_loop(A_def0, A_cb0, eps, g)
        y0 = c0["Y_D__shock_def_D"] @ eps + c0["Y_D__cb_buy_D"] @ cb0
        sev[name] = (peak(sp0) * BP_ANN, -y0[:100].min() * 100)
        log(f"  {name:>10}: peak spread {sev[name][0]:6.2f} bp, output loss {sev[name][1]:.5f}%")
    rank_ok = sev["aggressive"][0] < sev["medium"][0] < sev["passive"][0]
    log(f"- **A6 spread ranking survives at psi_lambda_B=0: {'YES' if rank_ok else 'NO'}** "
        f"(aggressive < medium < passive in peak spread = crisis severity)")

    log("\nStage A (main) complete.")


if __name__ == "__main__":
    main()
