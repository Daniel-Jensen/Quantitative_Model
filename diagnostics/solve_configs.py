"""
Diagnostic solve: baseline (psi_lambda_B=3.0) vs psi_lambda_B=0.0.

Reproduces the reported symptom (null default-shock IRFs at psi_lambda=0) and
caches all IRF series + summary statistics so downstream analysis needs no
re-solve. DIAGNOSTIC ONLY — imports the unmodified model from code/, writes
only to diagnostics/.

Runs, for each of {baseline, psi_lambda0}:
  - a 1pp sovereign default shock to D  (the shock under investigation)
  - a 1% TFP shock to D                 (control: proves the model still
                                          transmits *something* at psi_lambda=0)

Key comparison variables (the transmission chain):
  def_rate_D  -> q_b_D -> rb_actual_D/spread_rb -> n_inter_D -> theta_D
              -> K_D/I_D -> Y_D/C_D/w_D ; plus fiscal b_gov_D/TAX_D.

Outputs:
  diagnostics/irfs_baseline.npz, diagnostics/irfs_psilam0.npz
  diagnostics/summary.json
  appends a timestamped section to diagnostics/run_log.md
"""
import os, sys, json, copy, datetime, traceback
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CODE = os.path.join(ROOT, "code")
sys.path.insert(0, CODE)

RUN_LOG = os.path.join(HERE, "run_log.md")

# Optional output tag: `python solve_configs.py postfix` writes *_postfix artifacts
# so pre-fix (committed) and post-fix evidence coexist without clobbering.
TAG = sys.argv[1] if len(sys.argv) > 1 else ""
SUF = f"_{TAG}" if TAG else ""


def ts():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg=""):
    line = f"{msg}"
    print(line, flush=True)
    with open(RUN_LOG, "a") as f:
        f.write(line + "\n")


def logts(msg):
    log(f"- `{ts()}` {msg}")


# core series to persist in full (length T) for both configs / both shocks
CORE_VARS = [
    "def_rate_D", "q_b_D", "q_b_F", "rb_actual_D", "rb_actual_F", "spread_rb",
    "n_inter_D", "n_inter_F", "theta_D", "theta_F", "K_D", "I_D",
    "Y_D", "Y_F", "C_D", "C_F", "w_D", "b_gov_D", "TAX_D", "rdep_D", "div_D",
]


def peak(x, n=100):
    x = np.asarray(x)
    return float(np.max(np.abs(x[:n])))


def impact(x):
    return float(np.asarray(x)[0])


def irf_to_dict(irf):
    """Convert an SSJ ImpulseDict / dict result to a plain {name: 1d np.array}."""
    out = {}
    for k in irf.keys():
        try:
            v = np.asarray(irf[k])
        except Exception:
            continue
        if v.ndim == 1:
            out[k] = v
    return out


def summarize(irf, shock_name):
    d = irf_to_dict(irf)
    rows = {}
    for v in CORE_VARS:
        if v in d:
            rows[v] = {"impact": impact(d[v]), "peak_abs_100": peak(d[v])}
        else:
            rows[v] = {"impact": None, "peak_abs_100": None, "missing": True}
    return rows


def main():
    log("\n---\n")
    log(f"## Solve run{(' ['+TAG+']') if TAG else ''} — {ts()}")
    logts("Importing pipeline modules from code/ (unmodified).")

    from calibration import get_calibration
    from steady_state import solve_steady_state
    from ic_delta_calibration import calibrate_ic_delta
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve

    # ---- Steps 1-4: calibration -> SS -> ic_delta -> depreciation (once) ----
    logts("Step 1: get_calibration()")
    cal = get_calibration()
    logts(f"psi_lambda_B_D (baseline) = {cal['psi_lambda_B_D']}, "
          f"writeoff_enabled_D = {cal['writeoff_enabled_D']}, T1_D = {cal['T1_D']}, "
          f"mv_rule_D = {cal['mv_rule_D']}")

    logts("Step 2: solve_steady_state()")
    ss_results = solve_steady_state(cal)
    logts("Step 3: calibrate_ic_delta()")
    ss_results = calibrate_ic_delta(ss_results)
    logts("Step 4: calibrate_depreciation()  (final SS re-solve)")
    ss_results = calibrate_depreciation(ss_results)

    ss_final = ss_results["ss_final"]

    # ---- SS sanity: market clearing + psi params present (Step 5 sanity) ----
    ss_checks = {}
    for k in ["goods_mkt_D", "goods_mkt_F", "ca_res_D", "deposit_mkt_D",
              "deposit_mkt_F"]:
        try:
            ss_checks[k] = float(ss_final[k])
        except Exception as e:
            ss_checks[k] = f"MISSING ({e})"
    for k in ["psi_lambda_B_D", "psi_lambda_B_F", "psi_spread_D", "psi_spread_F",
              "q_b_D", "q_b_F", "rb_D", "rb_F", "n_inter_D", "theta_D",
              "beta_D", "beta_F", "p"]:
        try:
            ss_checks[k] = float(ss_final[k])
        except Exception as e:
            ss_checks[k] = f"MISSING ({e})"
    log("")
    log("### SS sanity checks (market clearing + key SS values)")
    for k, v in ss_checks.items():
        log(f"    {k:16s} = {v}")
    # spread at SS should be ~0 (equal duration/price) -> dormant TPI premise
    try:
        ss_spread = float(ss_final["spread_rb"])
        log(f"    spread_rb (SS)   = {ss_spread}   (expect ~0)")
        ss_checks["spread_rb"] = ss_spread
    except Exception as e:
        log(f"    spread_rb (SS)   = MISSING ({e})")

    # ---- Baseline Jacobian (psi_lambda_B = 3.0) via build_and_solve ----
    log("")
    logts("Building baseline Jacobian G (psi_lambda_B = 3.0) via build_and_solve()...")
    res = build_and_solve(ss_results)
    ha_full = res["ha_full"]
    unknowns_tp = res["unknowns_tp"]
    targets_tp = res["targets_tp"]
    T = res["T"]
    dShock_def_D = res["dShock_def_D"]
    dZ_D = res["dZ_D"]
    logts(f"Baseline G computed. T={T}. shock peaks: def={dShock_def_D[0]:.4g}, "
          f"TFP={dZ_D[0]:.4g}")

    irf_def_base = res["irfs_def_D"]           # default shock, baseline
    irf_tfp_base = res["irfs_Z_D"]             # TFP shock, baseline

    # ---- psi_lambda_B = 0 Jacobian (same SS, same model) ----
    log("")
    logts("Building psi_lambda_B = 0 Jacobian G0 (zero psi_lambda_B_* and psi_spread_*)...")
    ss0 = copy.deepcopy(ss_final)
    for k in ["psi_lambda_B_D", "psi_lambda_B_F", "psi_spread_D", "psi_spread_F"]:
        ss0.toplevel[k] = 0.0
    logts("Confirming ss0 overrides: " + ", ".join(
        f"{k}={float(ss0[k])}" for k in
        ["psi_lambda_B_D", "psi_lambda_B_F", "psi_spread_D", "psi_spread_F"]))

    exogenous = ["Z_D", "shock_def_D", "Z_F", "shock_def_F"]
    G0 = ha_full.solve_jacobian(
        ss0, unknowns=unknowns_tp, targets=targets_tp, inputs=exogenous, T=T)
    logts("G0 computed.")

    zero = np.zeros(T)
    irf_def_0 = G0 @ {"Z_D": zero, "Z_F": zero,
                      "shock_def_D": dShock_def_D, "shock_def_F": zero}
    irf_tfp_0 = G0 @ {"Z_D": dZ_D, "Z_F": zero,
                      "shock_def_D": zero, "shock_def_F": zero}

    # ---- Persist full core series ----
    def dump(path, irf_def, irf_tfp):
        dd = irf_to_dict(irf_def)
        dt = irf_to_dict(irf_tfp)
        payload = {}
        for v in CORE_VARS:
            if v in dd:
                payload[f"def__{v}"] = dd[v]
            if v in dt:
                payload[f"tfp__{v}"] = dt[v]
        np.savez(path, **payload)
        return sorted(payload.keys())

    kb = dump(os.path.join(HERE, f"irfs_baseline{SUF}.npz"), irf_def_base, irf_tfp_base)
    k0 = dump(os.path.join(HERE, f"irfs_psilam0{SUF}.npz"), irf_def_0, irf_tfp_0)
    logts(f"Saved irfs_baseline{SUF}.npz ({len(kb)} series), irfs_psilam0{SUF}.npz ({len(k0)} series).")

    # ---- Summary statistics ----
    summary = {
        "timestamp": ts(),
        "git_note": "see diagnostics/env.txt",
        "T": T,
        "shock_sizes": {"def_impact": float(dShock_def_D[0]),
                        "tfp_impact": float(dZ_D[0])},
        "ss_checks": ss_checks,
        "baseline_psi_lambda_B": float(cal["psi_lambda_B_D"]),
        "default_shock": {
            "baseline": summarize(irf_def_base, "def"),
            "psi_lambda0": summarize(irf_def_0, "def"),
        },
        "tfp_shock_control": {
            "baseline": summarize(irf_tfp_base, "tfp"),
            "psi_lambda0": summarize(irf_tfp_0, "tfp"),
        },
    }
    with open(os.path.join(HERE, f"summary{SUF}.json"), "w") as f:
        json.dump(summary, f, indent=2)
    logts(f"Wrote summary{SUF}.json")

    # ---- Human-readable comparison tables in the log ----
    def table(title, base_rows, zero_rows):
        log("")
        log(f"### {title}")
        log("")
        log(f"| var | impact (base) | peak|·| (base) | impact (ψλ=0) | peak|·| (ψλ=0) | peak ratio ψλ0/base |")
        log(f"|-----|--------------:|---------------:|--------------:|---------------:|--------------------:|")
        for v in CORE_VARS:
            b = base_rows.get(v, {})
            z = zero_rows.get(v, {})
            bi, bp = b.get("impact"), b.get("peak_abs_100")
            zi, zp = z.get("impact"), z.get("peak_abs_100")
            if bp not in (None, 0) and zp is not None:
                ratio = f"{zp / bp:.3e}"
            else:
                ratio = "n/a"
            fmt = lambda x: "     ---" if x is None else f"{x:.4e}"
            log(f"| {v} | {fmt(bi)} | {fmt(bp)} | {fmt(zi)} | {fmt(zp)} | {ratio} |")

    table("Default shock (1pp): baseline vs psi_lambda_B=0",
          summary["default_shock"]["baseline"], summary["default_shock"]["psi_lambda0"])
    table("TFP shock (1%, CONTROL): baseline vs psi_lambda_B=0",
          summary["tfp_shock_control"]["baseline"], summary["tfp_shock_control"]["psi_lambda0"])

    log("")
    logts("SOLVE RUN COMPLETE.")
    log("\n---\n")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        log("")
        log("### !!! SOLVE RUN FAILED !!!")
        log("```")
        log(traceback.format_exc())
        log("```")
        raise
