"""
psi_lambda_B sweep with EL_price ANCHORED (macro-pru-fix), to discipline the
amplification dial by external moments rather than by continuity with the old
calibration.

For each psi_lambda_B on a grid: hold EL_price fixed at its SS-anchored value,
set psi_spread = psi_spread_ss * (psi_lambda_B / psi_lambda_B_base) (exact —
psi_spread is linear in psi_lambda_B), re-solve the Jacobian, apply the 1pp
default shock, and report the model-implied moments:

  - peak sovereign spread response (bp, annualised ×4)      <- match to observed GR spread
  - peak bank net worth response (% of SS)
  - peak output response (% of SS)
  - bank-networth-to-spread pass-through  Δn(%SS) / Δspread(100bp)   <- doom-loop moment
  - spread amplification vs the psi_lambda_B=0 fundamental floor

The psi_lambda_B=0 row is ALSO the Case-3-resolved reproduction (fundamental
floor from EL_price alone).

Outputs: diagnostics/psilam_moment_sweep.md, .npz, 05_psilam_moment_sweep.png
"""
import os, sys, json, copy, datetime, textwrap
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "code"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["savefig.dpi"] = 200


def ts():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


GRID = [0.0, 0.5, 1.0, 1.5, 2.0, 2.6, 3.0, 4.0, 5.0]
ANNUAL = 4.0          # quarterly holding-yield -> annualised
BP = 1e4


def main():
    from calibration import get_calibration
    from steady_state import solve_steady_state
    from ic_delta_calibration import calibrate_ic_delta
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve

    print(f"[{ts()}] pipeline: calibration -> SS -> ic_delta -> depreciation")
    cal = get_calibration()
    ssr = solve_steady_state(cal)
    ssr = calibrate_ic_delta(ssr)
    ssr = calibrate_depreciation(ssr)

    print(f"[{ts()}] build_and_solve (baseline psi_lambda_B=3)")
    res = build_and_solve(ssr)
    ha, ss = res["ha_full"], res["ss_final"]
    unk, tgt, T, dshock = res["unknowns_tp"], res["targets_tp"], res["T"], res["dShock_def_D"]

    psi_spread_base = float(ss["psi_spread_D"])
    psilam_base     = float(ss["psi_lambda_B_D"])
    EL_D            = float(ss["EL_price_D"])
    n_ss            = float(ss["n_inter_D"])
    Y_ss            = float(ss["Y_D"])
    print(f"[{ts()}] anchors: EL_price_D={EL_D:.6f} (FIXED), "
          f"psi_spread_D(base)={psi_spread_base:.6f} at psi_lambda_B={psilam_base}")

    zero = np.zeros(T)
    rows = []
    for g in GRID:
        ssg = copy.deepcopy(ss)
        ssg.toplevel["psi_lambda_B_D"] = g
        ssg.toplevel["psi_lambda_B_F"] = g
        ssg.toplevel["psi_spread_D"]   = psi_spread_base * g / psilam_base
        ssg.toplevel["psi_spread_F"]   = psi_spread_base * g / psilam_base
        # EL_price_D/F untouched (anchored)
        Gg = ha.solve_jacobian(ssg, unknowns=unk, targets=tgt,
                               inputs=["Z_D", "shock_def_D", "Z_F", "shock_def_F"], T=T)
        irf = Gg @ {"Z_D": zero, "Z_F": zero, "shock_def_D": dshock, "shock_def_F": zero}
        sp   = np.asarray(irf["spread_rb"])[:100]
        nD   = np.asarray(irf["n_inter_D"])[:100]
        YD   = np.asarray(irf["Y_D"])[:100]
        spread_bp = float(sp.max()) * ANNUAL * BP
        n_pct     = float(nD.min()) / n_ss * 100.0
        Y_pct     = float(YD.min()) / Y_ss * 100.0
        passthru  = n_pct / (spread_bp / 100.0) if abs(spread_bp) > 1e-9 else float("nan")
        rows.append(dict(psilam=g, psi_spread=psi_spread_base * g / psilam_base,
                         foc_load=EL_D + psi_spread_base * g / psilam_base,
                         spread_bp=spread_bp, n_pct=n_pct, Y_pct=Y_pct, passthru=passthru))
        print(f"[{ts()}]  psi_lambda_B={g:4.2f}  spread={spread_bp:7.2f}bp  "
              f"dn={n_pct:+.3f}%SS  dY={Y_pct:+.4f}%SS  passthru={passthru:+.3f}%/100bp")

    floor_bp = rows[0]["spread_bp"]  # psi_lambda_B=0
    for r in rows:
        r["amp_vs_floor"] = r["spread_bp"] / floor_bp if abs(floor_bp) > 1e-12 else float("nan")

    # ---- markdown table ----
    lines = [
        "# psi_lambda_B moment sweep (EL_price anchored) — macro-pru-fix",
        "",
        f"Generated {ts()}. `EL_price_D = {EL_D:.6f}` held FIXED (the empirically-anchored",
        "fundamental expected-loss loading, from recovery/duration). `psi_spread` scales",
        "linearly with `psi_lambda_B`. Moments are the response to a **1pp sovereign default",
        "shock**. Choose `psi_lambda_B` to match an external moment, NOT the old IRFs.",
        "",
        "| psi_lambda_B | psi_spread | FOC load (EL+ψs) | peak spread (bp, ann) | peak Δn_inter (%SS) | peak ΔY (%SS) | Δn per 100bp | spread amp vs ψλ=0 |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        tag = ""
        if r["psilam"] == 0.0:   tag = " ← fundamental floor / Case-3 check"
        if r["psilam"] == 3.0:   tag = " ← current baseline"
        if r["psilam"] == 2.6:   tag = " ← loading-matched to old"
        lines.append(
            f"| {r['psilam']:.2f}{tag} | {r['psi_spread']:.3f} | {r['foc_load']:.3f} | "
            f"{r['spread_bp']:.1f} | {r['n_pct']:+.3f} | {r['Y_pct']:+.4f} | "
            f"{r['passthru']:+.3f} | {r['amp_vs_floor']:.2f}× |")
    lines += [
        "",
        "## How to read this",
        "- **Disciplining moment = the spread level / amplification** (highly `psi_lambda_B`-sensitive):",
        "  pick `psi_lambda_B` so the peak spread matches the observed Greek response to a",
        "  comparable default-probability move (2010 GR–DE spread ≈ 150bp is the paper's target).",
        "- **Consistency moment = the bank-networth-to-spread pass-through** `Δn per 100bp`",
        "  (roughly `psi_lambda_B`-robust): compare to bank-equity/sovereign-spread event studies",
        "  (e.g. Altavilla–Pagano–Simonelli; Acharya–Drechsler–Schnabl). If it sits in the",
        "  empirical range across the sweep, the doom-loop transmission is the right size.",
        "- **`psi_lambda_B = 0` row is the fundamental floor**: nonzero, correctly-signed response",
        "  from `EL_price` alone — the Case-3 null is resolved independent of the dial.",
        "",
        "NOTE: exact empirical target values are a literature-retrieval task (flagged in the",
        "handoff); this table supplies the model side of the mapping so the dial is set by data.",
    ]
    with open(os.path.join(HERE, "psilam_moment_sweep.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    np.savez(os.path.join(HERE, "psilam_moment_sweep.npz"),
             psilam=np.array([r["psilam"] for r in rows]),
             spread_bp=np.array([r["spread_bp"] for r in rows]),
             n_pct=np.array([r["n_pct"] for r in rows]),
             Y_pct=np.array([r["Y_pct"] for r in rows]),
             passthru=np.array([r["passthru"] for r in rows]),
             EL_price_D=EL_D)

    # ---- figure ----
    g   = np.array([r["psilam"] for r in rows])
    spb = np.array([r["spread_bp"] for r in rows])
    npc = np.array([r["n_pct"] for r in rows])
    pt  = np.array([r["passthru"] for r in rows])
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].plot(g, spb, "o-", color="#8C1515"); ax[0].set_title("peak spread (bp, ann)")
    ax[1].plot(g, npc, "o-", color="#002147"); ax[1].set_title("peak Δ net worth (% SS)")
    ax[2].plot(g, pt,  "o-", color="#1a6e3a"); ax[2].set_title("Δn per 100bp (pass-through)")
    for a in ax:
        a.axvline(3.0, color="0.7", ls="--", lw=1, label="current (3.0)")
        a.axvline(0.0, color="0.85", ls=":", lw=1)
        a.set_xlabel("psi_lambda_B (EL_price fixed)")
    ax[0].legend(fontsize=8)
    fig.suptitle("psi_lambda_B disciplined by moments (EL_price anchored) — macro-pru-fix", fontsize=12)
    cap = ("With the fundamental expected-loss loading EL_price held fixed, the peak sovereign "
           "spread is strongly increasing in the amplification dial psi_lambda_B (left) while the "
           "bank-networth-to-spread pass-through (right) is comparatively flat — so the spread "
           "level disciplines psi_lambda_B against data and the pass-through is a consistency check; "
           "psi_lambda_B=0 is the nonzero fundamental floor (Case-3 resolved).")
    chars = int(fig.get_size_inches()[0] * 14)
    fig.text(0.5, -0.04, textwrap.fill(cap, width=chars), ha="center", va="top",
             fontsize=8, style="italic", color="0.35")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(HERE, "05_psilam_moment_sweep.png"), bbox_inches="tight")
    print(f"[{ts()}] wrote psilam_moment_sweep.md/.npz and 05_psilam_moment_sweep.png")


if __name__ == "__main__":
    main()
