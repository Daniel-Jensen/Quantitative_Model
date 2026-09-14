"""
psilam_breakdown_sweep.py — re-derive PSILAM_BREAKDOWN for the broad-scope EBA
calibration (2026-07-31).

`regime_model.PSILAM_BREAKDOWN` is a hard guard: the regimes diagnostics refuse
to run at a psi_lambda_B where the linear approximation has broken down. The
threshold is CALIBRATION-DEPENDENT (it scales with bank net worth), and the
committed value 2.5 was derived for the CT1-scope net worth n_inter_D=0.408.
Under BANK_SCOPE="broad" the live calibration is psi_lambda_B=8.5 with
n_inter_D=2.138, so the old threshold both blocks the run and carries no
information about where breakdown actually is.

Method. psi_spread is EXACTLY linear in psi_lambda_B —

    psi_spread_D = lambda_gk_D * psi_lambda_B_D / (beta_inter_D * Omega_D)

(steady_state._apply_ss_anchors) — and lambda_gk / Omega / beta_inter are solved
without reference to psi_lambda_B. So the steady state is invariant and only the
Jacobian has to be re-solved per grid point. Both dials MUST move together:
patching psi_lambda_B alone leaves psi_spread stale and inverts the apparent sign
of the spread response (the void-sweep caveat in docs/eba_calibration.md).

Breakdown criteria, in order of severity:
  1. peak spread > 1000bp ann      — run_regimes.py's A7 flag; linearisation gone
  2. b_gov_D[499] not ~ 0          — debt mode non-stationary
  3. sign flip: n_inter_D[0] or Y_D[0] positive on a default shock
  4. non-monotone peak spread in psi_lambda_B

Outputs: diagnostics/psilam_breakdown_sweep.md + .npz
"""
import os, sys, copy, datetime
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "code"))

# Above 9 is the unexplored region: docs/eba_calibration.md already records
# 5.4 / 23.1 / 58.6 / 111.0 / 142.5 / 157.8 bp at 0 / 1 / 3 / 6 / 8 / 9.
# Refined 26-29 after the coarse pass put a pole between 25 and 30, and 15-18
# after n_inter_D[0] was seen to reverse there while the spread kept rising.
GRID = [8.5, 10.0, 12.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 22.0,
        25.0, 26.0, 27.0, 28.0, 29.0, 30.0]
BP_ANN = 4.0 * 1e4


def ts():
    return datetime.datetime.now().strftime("%H:%M:%S")


def main():
    from calibration import get_calibration
    from steady_state import solve_steady_state
    from ic_delta_calibration import calibrate_ic_delta
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve, solve_jacobian_padded

    cal = get_calibration()
    print(f"[{ts()}] live psi_lambda_B_D = {cal['psi_lambda_B_D']}, "
          f"BANK_SCOPE via n_inter_D = {cal['n_inter_D']}")
    ssr = calibrate_depreciation(calibrate_ic_delta(solve_steady_state(cal)))
    res = build_and_solve(ssr)
    ha, ss = res["ha_full"], res["ss_final"]
    unk, tgt, T, dshock = res["unknowns_tp"], res["targets_tp"], res["T"], res["dShock_def_D"]

    # The exact anchor formula, re-evaluated per grid point (never rescaled by hand).
    def psi_spread(country, g):
        return (float(ss[f"lambda_gk_{country}"]) * g
                / (float(ss[f"beta_inter_{country}"]) * float(ss[f"Omega_{country}"])))

    base = float(ss["psi_spread_D"])
    assert abs(psi_spread("D", cal["psi_lambda_B_D"]) - base) < 1e-10 * max(1.0, abs(base)), \
        "psi_spread reconstruction disagrees with the solved SS — anchor formula drifted"
    print(f"[{ts()}] psi_spread_D reconstruction OK ({base:.6f}); "
          f"n_inter_D_ss={float(ss['n_inter_D']):.4f}, lambda_gk_D={float(ss['lambda_gk_D']):.4f}")

    zero = np.zeros(T)
    rows = []
    for g in GRID:
        ssg = copy.deepcopy(ss)
        for c in ("D", "F"):
            ssg.toplevel[f"psi_lambda_B_{c}"] = g
            ssg.toplevel[f"psi_spread_{c}"] = psi_spread(c, g)
        Gg = solve_jacobian_padded(ha, ssg, unk, tgt,
                                   ["Z_D", "shock_def_D", "Z_F", "shock_def_F"], T)
        irf = Gg @ {"Z_D": zero, "Z_F": zero, "shock_def_D": dshock, "shock_def_F": zero}
        sp = np.asarray(irf["spread_rb"])
        r = dict(psilam=g,
                 psi_spread=psi_spread("D", g),
                 spread_bp=float(sp[:100].max()) * BP_ANN,
                 b_gov_end=float(np.asarray(irf["b_gov_D"])[T - 1]),
                 n0=float(np.asarray(irf["n_inter_D"])[0]) / float(ss["n_inter_D"]) * 100.0,
                 Y0=float(np.asarray(irf["Y_D"])[0]) / float(ss["Y_D"]) * 100.0)
        r["ok"] = (r["spread_bp"] < 1000.0 and abs(r["b_gov_end"]) < 1e-2
                   and r["n0"] < 0 and r["Y0"] < 0)
        rows.append(r)
        print(f"[{ts()}]  psi_lambda_B={g:5.1f}  spread={r['spread_bp']:9.1f}bp  "
              f"b_gov[T-1]={r['b_gov_end']:+.3e}  n_inter[0]={r['n0']:+.3f}%  "
              f"Y[0]={r['Y0']:+.4f}%  {'ok' if r['ok'] else 'BREAKDOWN'}", flush=True)

    sp = np.array([r["spread_bp"] for r in rows])
    monotone_upto = GRID[-1]
    for i in range(1, len(sp)):
        if sp[i] <= sp[i - 1]:
            monotone_upto = GRID[i - 1]
            break
    first_bad = next((r["psilam"] for r in rows if not r["ok"]), None)
    # The net-worth response can reverse while the spread still rises — a milder
    # pathology than the pole, and the one that actually sets the guard.
    n0 = np.array([r["n0"] for r in rows])
    rev = [GRID[i] for i in range(1, len(n0)) if n0[i] > n0[i - 1]]
    first_rev = rev[0] if rev else None

    lines = [
        "# psi_lambda_B breakdown sweep — broad-scope EBA calibration",
        "",
        f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M:%S}. Steady state solved once "
        f"(psi_lambda_B enters only via the psi_spread anchor, which is exactly linear in it); "
        f"Jacobian re-solved per point with BOTH dials moved together.",
        "",
        f"`n_inter_D_ss = {float(ss['n_inter_D']):.4f}`, `lambda_gk_D = {float(ss['lambda_gk_D']):.4f}`, "
        f"`Omega_D = {float(ss['Omega_D']):.4f}`.",
        "",
        "| psi_lambda_B | psi_spread_D | peak spread (bp ann) | b_gov_D[T-1] | n_inter_D[0] (%) | Y_D[0] (%) | verdict |",
        "|---:|---:|---:|---:|---:|---:|:--|",
    ]
    for r in rows:
        lines.append(f"| {r['psilam']:.1f} | {r['psi_spread']:.3f} | {r['spread_bp']:.1f} | "
                     f"{r['b_gov_end']:+.2e} | {r['n0']:+.3f} | {r['Y0']:+.4f} | "
                     f"{'ok' if r['ok'] else '**BREAKDOWN**'} |")
    lines += [
        "",
        f"- peak spread monotone increasing up to psi_lambda_B = **{monotone_upto}**",
        f"- first psi_lambda_B at which n_inter_D[0] REVERSES (shrinks while the spread "
        f"still rises): **{first_rev if first_rev is not None else 'none on this grid'}**",
        f"- first breakdown row (A7 / stationarity / sign): "
        f"**{first_bad if first_bad is not None else 'none on this grid'}**",
        "",
        "Breakdown criteria: peak spread > 1000bp (A7), |b_gov_D[T-1]| > 1e-2, or a sign flip in "
        "n_inter_D[0] / Y_D[0] on a default shock.",
        "",
        "`regime_model.PSILAM_BREAKDOWN` is set from the FIRST pathology (the n_inter_D[0] "
        "reversal), not from the pole, so the guard has real margin rather than sitting on the "
        "edge of the singularity.",
    ]
    with open(os.path.join(HERE, "psilam_breakdown_sweep.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    np.savez(os.path.join(HERE, "psilam_breakdown_sweep.npz"),
             **{k: np.array([r[k] for r in rows]) for k in
                ("psilam", "psi_spread", "spread_bp", "b_gov_end", "n0", "Y0")})
    print(f"[{ts()}] monotone up to {monotone_upto}; first breakdown {first_bad}")


if __name__ == "__main__":
    main()
