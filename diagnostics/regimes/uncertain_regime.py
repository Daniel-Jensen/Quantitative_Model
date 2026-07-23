"""
uncertain_regime.py — Stage B-lite: the CB's type {aggressive, medium, passive}
is unknown when the default shock hits; banks price the lottery under beliefs pi;
the CB is silent until revelation date k (spec §10). Pure numpy on the
regime_model cache (MAIN model) + Stage A's regimes_calibration.json + beliefs.json.

Assertion boundary:
  * HARD (structural, blocking): the §10.3 checks — pre-k branch identity, the
    pi-weighted revelation jump = 0 (rational expectations), and k=0 nesting of
    Stage A. Mathematical identities of the construction, independent of A_cb's
    sign; they MUST hold.
  * REPORTED (computed findings, NOT pass/fail): every economic SIGN — whether the
    impact spread rises or falls in pi_passive, the sign of the uncertainty premium,
    the A6 lottery ranking. On MAIN A_cb<0 (purchases COMPRESS the spread), so the
    intuitive spec-§10 prediction is expected — fear of a passive CB raises today's
    spread — but the sign is still computed, and captions are generated from it so
    the narrative is honest whichever way it falls.
"""
import os, sys, json, datetime
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from regime_model import build_caches, load_cache, PSILAM_MAIN
from lottery_math import shift_k, closed_loop, solve_lottery, branch_output, peak

LOG = os.path.join(HERE, "regimes_log.md")
BLUE, RED, GREEN, ORANGE = "#002147", "#8C1515", "#1a6e3a", "#c87941"
BP_ANN = 4.0 * 1e4
NAMES = ("aggressive", "medium", "passive")
COLORS = dict(zip(NAMES, (GREEN, ORANGE, RED)))
K_DEFAULT, K_SWEEP = 2, (1, 2, 4)


def log(m=""):
    print(m, flush=True)
    with open(LOG, "a") as f:
        f.write(m + "\n")


def assert_lottery(spreads, cbs, cb_e, pi, k, A_cb):
    """§10.3 blocking assertions — structural identities, independent of A_cb sign."""
    assert np.max(np.abs(spreads[:, :k] - spreads[0, :k])) < 1e-9, \
        "BLOCKING: branch paths differ pre-revelation (shared-information violation)"
    jump = np.stack([shift_k(A_cb, k) @ (cbs[s] - cb_e) for s in range(len(pi))])
    assert np.max(np.abs((pi[:, None] * jump).sum(0))) < 1e-9, \
        "BLOCKING: pi-weighted revelation jump != 0 (rational-expectations violation)"


def main():
    open(LOG, "a").write(f"\n\n## Stage B-lite run (MAIN model) — {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n")
    build_caches()
    c = load_cache(PSILAM_MAIN)
    T = int(c["T"]); eps = np.asarray(c["dShock_def_D"])
    A_def, A_cb = c["spread_rb__shock_def_D"], c["spread_rb__cb_buy_D"]
    cal = json.load(open(os.path.join(HERE, "regimes_calibration.json")))
    gammas = np.array([cal["gammas"][n] for n in NAMES])
    B = json.load(open(os.path.join(HERE, "beliefs.json")))
    pi_onset = np.array([B["pi_onset_1q"][n] for n in NAMES])
    pi_erg   = np.array([B["pi_ergodic"][n]  for n in NAMES])
    log(f"- gammas (spread-compression, spec §7): {dict(zip(NAMES, np.round(gammas, 3)))}")
    log(f"- impact A_cb = {A_cb[0,0]:+.3e} ({'compresses' if A_cb[0,0] < 0 else 'widens'} on main)")
    log(f"- pi_onset = {dict(zip(NAMES, np.round(pi_onset, 4)))}, "
        f"pi_ergodic = {dict(zip(NAMES, np.round(pi_erg, 4)))}")

    # ── Baseline lottery: onset beliefs, k=2 ─────────────────────────────────
    k = K_DEFAULT
    spreads, cbs, cb_e = solve_lottery(A_def, A_cb, eps, gammas, pi_onset, k)
    assert_lottery(spreads, cbs, cb_e, pi_onset, k, A_cb)
    sp_k0, _, _ = solve_lottery(A_def, A_cb, eps, gammas, pi_onset, 0)
    for s, g in enumerate(gammas):
        sp_known, _ = closed_loop(A_def, A_cb, eps, g)
        assert np.allclose(sp_k0[s], sp_known, atol=1e-9), "BLOCKING: k=0 does not nest Stage A"
    log(f"- §10.3 assertions PASS (pre-k identity, RE jump=0, Stage A nesting) at pi_onset, k={k}")
    if peak(spreads[NAMES.index("passive")]) * BP_ANN > 1000:
        log("  **A7 FLAG (passive branch): >1000bp — linear-approximation breakdown, not a finding.**")

    # ── Decomposition: delay cost vs uncertainty premium (peak spread, bp) ────
    log("\n| branch | known-immediate | known-delayed(k) | lottery | delay cost | uncertainty premium |")
    log("|---|---|---|---|---|---|")
    Pi_k = np.diag((np.arange(T) >= k).astype(float))
    for s, n in enumerate(NAMES):
        sp_imm, _ = closed_loop(A_def, A_cb, eps, gammas[s])
        sp_del = np.linalg.solve(np.eye(T) - gammas[s] * A_cb @ Pi_k, A_def @ eps)
        imm, dl, lot = peak(sp_imm), peak(sp_del), peak(spreads[s])
        log(f"| {n} | {imm*BP_ANN:.1f} | {dl*BP_ANN:.1f} | {lot*BP_ANN:.1f} | "
            f"{(dl-imm)*BP_ANN:+.1f} | {(lot-dl)*BP_ANN:+.1f} |")
    log("(all peak spread, bp ann; delay cost = known-delayed - known-immediate; "
        "uncertainty premium = lottery - known-delayed)")

    # ── Welfare per branch + expected (tpi.py convention: 100q, beta^t) ──────
    bD, bF = float(c["beta_D"]), float(c["beta_F"])
    dD, dF = bD ** np.arange(100), bF ** np.arange(100)
    qb = float(c["q_b_D_ss"])
    W = {}
    log("\n| branch | W_D | W_F | discounted CB purchases (A5, per branch) |")
    log("|---|---|---|---|")
    for s, n in enumerate(NAMES):
        UD = branch_output(c["U_D__shock_def_D"], c["U_D__cb_buy_D"], eps, cb_e, cbs[s], k)
        UF = branch_output(c["U_F__shock_def_D"], c["U_F__cb_buy_D"], eps, cb_e, cbs[s], k)
        buy = (cbs[s][:100] * qb * dD).sum()
        W[n] = ((UD[:100] * dD * 100).sum(), (UF[:100] * dF * 100).sum(), buy)
        log(f"| {n} | {W[n][0]:+.4f} | {W[n][1]:+.4f} | {W[n][2]:.5f} |")
    EW = tuple((pi_onset * np.array([W[n][i] for n in NAMES])).sum() for i in (0, 1, 2))
    log(f"| **E_pi** | {EW[0]:+.4f} | {EW[1]:+.4f} | {EW[2]:.5f} |")
    log("(welfare: % SS cons., 100q, tpi.py convention; purchases: Sigma beta^t q_b cb_t)")

    # ── Figure: branch paths + revelation jumps ──────────────────────────────
    Nplot = 24
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for s, n in enumerate(NAMES):
        axes[0].plot(np.arange(Nplot), spreads[s][:Nplot] * BP_ANN, color=COLORS[n], label=n)
        axes[1].plot(np.arange(Nplot), cbs[s][:Nplot], color=COLORS[n])
        yD = branch_output(c["Y_D__shock_def_D"], c["Y_D__cb_buy_D"], eps, cb_e, cbs[s], k)
        axes[2].plot(np.arange(Nplot), yD[:Nplot] * 100, color=COLORS[n])
    for ax, t in zip(axes, ("spread_rb (bp ann)", "CB purchases", "Y_D (%)")):
        ax.axvline(k - 0.5, color="gray", lw=0.8, ls="--"); ax.set_title(t, fontsize=10)
        ax.axhline(0, lw=0.5, color="gray")
    axes[0].legend(fontsize=8)
    fig.suptitle(f"Stage B-lite (main): CB-type lottery, onset beliefs, revelation at k={k}")
    fig.text(0.5, 0.02,
             "Before the CB's type is revealed (dashed line) all branches share one path in which "
             "the spread already prices the policy lottery; at revelation the branches split into "
             "announcement-effect jumps whose belief-weighted average is exactly zero (rational "
             "expectations). The pre-revelation level is the regime-uncertainty object (Stage B-only).",
             ha="center", fontsize=8, style="italic", wrap=True)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.savefig(os.path.join(HERE, "fig_stageB_lottery.png"), dpi=200)

    # ── Sweep: impact/peak spread vs pi_passive (ergodic -> pure passive) ────
    # SIGN IS A FINDING. Caption generated from the computed direction.
    ts = np.linspace(0, 1, 11)
    pure_passive = np.array([0.0, 0.0, 1.0])
    imp, pk, pi_pass = [], [], []
    for t in ts:
        pi_t = (1 - t) * pi_erg + t * pure_passive
        pi_t = pi_t / pi_t.sum()
        sp_t, cb_t, ce_t = solve_lottery(A_def, A_cb, eps, gammas, pi_t, k)
        imp.append(sp_t[0][0] * BP_ANN)                          # t=0, common pre-k across branches
        pk.append(peak(sp_t[NAMES.index("passive")]) * BP_ANN)
        pi_pass.append(pi_t[2])
    imp, pk, pi_pass = np.array(imp), np.array(pk), np.array(pi_pass)
    d_imp = np.diff(imp)
    if np.all(d_imp > 1e-12):
        direction, story = "rises", ("fear of a passive central bank is priced before any policy "
                                     "acts — the impact spread rises with belief weight on the passive "
                                     "type (the intuitive spec-§10 sign: on main purchases compress the "
                                     "spread, so expecting no backstop raises today's spread)")
    elif np.all(d_imp < -1e-12):
        direction, story = "falls", ("the impact spread falls as belief shifts to the passive type — "
                                     "the counter-intuitive sign; investigate before reporting (would "
                                     "require the anticipation channel to dominate with A_cb<0)")
    else:
        direction, story = "non-monotone", ("the impact spread is non-monotone in the belief weight "
                                            "on the passive type (see the run log for the profile)")
    log(f"\n- impact spread vs pi_passive: **{direction}** "
        f"(from {imp[0]:.2f} bp at pi_passive={pi_pass[0]:.3f} to {imp[-1]:.2f} bp at "
        f"pi_passive={pi_pass[-1]:.3f}); this is the regime-uncertainty price, sign computed not targeted.")

    fig2, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.plot(pi_pass, imp, color=BLUE, marker="o", lw=3, label="impact spread (t=0, pre-revelation)")
    ax.plot(pi_pass, pk, color=RED, marker="s", ls="--", label="peak spread, passive branch")
    ax.set_xlabel("belief weight on passive CB (pi_passive)"); ax.set_ylabel("bp ann")
    ax.legend(fontsize=8)
    ax.annotate("(impact spread and passive-branch peak coincide:\nthe peak is pre-revelation)",
                xy=(0.5, 0.12), xycoords="axes fraction", ha="center", fontsize=7, color="gray")
    fig2.suptitle("Stage B-lite (main): regime-uncertainty price vs belief on the passive type")
    fig2.text(0.5, 0.02,
              f"As beliefs shift from the ergodic mix toward certainty on the passive type, the "
              f"pre-revelation impact spread {direction}: {story}. "
              f"No steady-state spread differs across types — this is a Stage-B-only object.",
              ha="center", fontsize=8, style="italic", wrap=True)
    fig2.tight_layout(rect=(0, 0.20, 1, 1))
    fig2.savefig(os.path.join(HERE, "fig_stageB_premium_vs_pi.png"), dpi=200)

    # ── Sweep: k (announcement-delay cost), onset-k-consistent beliefs ───────
    log("\n| k | impact spread (bp) | passive-branch peak (bp) | E_pi[W_D] |")
    log("|---|---|---|---|")
    for kk in K_SWEEP:
        pi_k = np.array([B["pi_onset_k"][str(kk)][n] for n in NAMES])
        sp_k, cb_k, ce_k = solve_lottery(A_def, A_cb, eps, gammas, pi_k, kk)
        assert_lottery(sp_k, cb_k, ce_k, pi_k, kk, A_cb)
        WDk = []
        for s in range(3):
            UD = branch_output(c["U_D__shock_def_D"], c["U_D__cb_buy_D"], eps, ce_k, cb_k[s], kk)
            WDk.append((UD[:100] * dD * 100).sum())
        log(f"| {kk} | {sp_k[0][0]*BP_ANN:.1f} | {peak(sp_k[2])*BP_ANN:.1f} | "
            f"{(pi_k*np.array(WDk)).sum():+.4f} |")

    # ── A6 extension: lottery ranking at psi_lambda_B = 0 (REPORTED) ─────────
    c0 = load_cache(0.0)
    sp0, cb0, ce0 = solve_lottery(c0["spread_rb__shock_def_D"], c0["spread_rb__cb_buy_D"],
                                  eps, gammas, pi_onset, k)
    pk0 = [peak(sp0[s]) * BP_ANN for s in range(3)]
    surv = pk0[0] < pk0[1] < pk0[2]
    log(f"\n- A6 (lottery, psi_lambda_B=0): branch peaks {dict(zip(NAMES, [round(p, 2) for p in pk0]))} bp — "
        f"aggressive<medium<passive ordering survives: {'YES' if surv else 'NO'} (reported)")

    log("\nStage B-lite (main) complete.")


if __name__ == "__main__":
    main()
