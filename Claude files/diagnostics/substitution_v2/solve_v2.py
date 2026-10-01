"""
Diagnostic v2 solve: psi_lambda_B = 2.8, EL_price on vs off, 1pp default shock.

SS is invariant to psi_lambda_B and EL_price (both ∝ def_rate(+1)=0 at SS), so we
solve the SS once and compute two Jacobians. Saves the full IRF set for each config
plus SS levels needed to convert deviations to %-of-SS and to build balance-sheet
identities. Diagnose only — imports unmodified code/, writes only to
diagnostics/substitution_v2/.
"""
import os, sys, json, copy, datetime
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "code"))
RUN_LOG = os.path.join(HERE, "run_log.md")

PSILAM = 2.8
EL_ANCHOR = None  # read from SS


def ts():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(m=""):
    print(m, flush=True)
    with open(RUN_LOG, "a") as f:
        f.write(m + "\n")


def irf_to_dict(irf):
    out = {}
    for k in irf.keys():
        try:
            v = np.asarray(irf[k])
        except Exception:
            continue
        if v.ndim == 1:
            out[k] = v
    return out


def main():
    from calibration import get_calibration
    from steady_state import solve_steady_state
    from ic_delta_calibration import calibrate_ic_delta
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve, solve_jacobian_padded

    log(f"\n## Step 0 — solve — {ts()}")
    cal = get_calibration()
    ssr = solve_steady_state(cal)
    ssr = calibrate_ic_delta(ssr)
    ssr = calibrate_depreciation(ssr)
    res = build_and_solve(ssr)
    ha, ss = res["ha_full"], res["ss_final"]
    unk, tgt, T, dshock = res["unknowns_tp"], res["targets_tp"], res["T"], res["dShock_def_D"]

    psi_spread_base = float(ss["psi_spread_D"])
    psilam_base     = float(ss["psi_lambda_B_D"])
    EL_anchor       = float(ss["EL_price_D"])
    psi_spread_28   = psi_spread_base * PSILAM / psilam_base

    # SS levels for normalisation / balance-sheet identities
    SS = {k: float(ss[k]) for k in [
        "n_inter_D", "K_D", "b_D_D", "b_F_D", "q_b_D", "q_b_F", "Q_D",
        "C_D", "I_D", "Y_D", "theta_D", "P_CES_D", "w_D", "N_D",
        "EL_price_D", "psi_spread_D", "rk_D", "mpk_D", "delta_D"]}
    SS["psi_spread_28"] = psi_spread_28
    try:
        SS["G_D"] = float(ss["G_D"])
    except Exception:
        SS["G_D"] = None
    try:
        SS["NX_D"] = float(ss["NX_D"])
    except Exception:
        SS["NX_D"] = None
    with open(os.path.join(HERE, "ss_values.json"), "w") as f:
        json.dump(SS, f, indent=2)

    log(f"- {ts()} anchors: EL_price_D(anchor)={EL_anchor:.6f}; "
        f"psi_spread_D(base@3)={psi_spread_base:.6f}; psi_spread@2.8={psi_spread_28:.6f}")
    log(f"- SS levels: n_inter={SS['n_inter_D']:.4f} K_D={SS['K_D']:.4f} "
        f"b_D_D={SS['b_D_D']:.4f} q_b_D={SS['q_b_D']:.4f} theta={SS['theta_D']:.3f} "
        f"Y_D={SS['Y_D']:.4f} I_D={SS['I_D']:.4f} C_D={SS['C_D']:.4f}")

    # SS market clearing (confirm we diagnose the intended SS)
    log("- SS market clearing: " + "  ".join(
        f"{k}={float(ss[k]):.2e}" for k in
        ["goods_mkt_D", "goods_mkt_F", "ca_res_D", "deposit_mkt_D"]))

    inputs = ["Z_D", "shock_def_D", "Z_F", "shock_def_F"]
    zero = np.zeros(T)

    def solve(el_price, tag):
        ssg = copy.deepcopy(ss)
        ssg.toplevel["psi_lambda_B_D"] = PSILAM
        ssg.toplevel["psi_lambda_B_F"] = PSILAM
        ssg.toplevel["psi_spread_D"]   = psi_spread_28
        ssg.toplevel["psi_spread_F"]   = psi_spread_28
        ssg.toplevel["EL_price_D"]     = el_price
        ssg.toplevel["EL_price_F"]     = el_price
        log(f"- {ts()} solving [{tag}] psi_lambda_B={PSILAM}, psi_spread={psi_spread_28:.6f}, "
            f"EL_price={el_price:.6f}")
        G = solve_jacobian_padded(ha, ssg, unk, tgt, inputs, T)
        irf = G @ {"Z_D": zero, "Z_F": zero, "shock_def_D": dshock, "shock_def_F": zero}
        d = irf_to_dict(irf)
        np.savez(os.path.join(HERE, f"irfs_2p8_{tag}.npz"), **d)
        # Step 0/1 quick checks
        sp = d.get("spread_rb"); qb = d.get("q_b_D"); Y = d.get("Y_D")
        K = d.get("K_D"); n = d.get("n_inter_D"); I = d.get("I_D"); b = d.get("b_D_D")
        def pk(x, mode="absmax"):
            if x is None:
                return float("nan")
            x = x[:100]
            return float(x[np.argmax(np.abs(x))])
        log(f"    [{tag}] shock reaches q_b: spread_rb impact={sp[0]:.4e} "
            f"peak(bp,ann)={sp[:100].max()*4*1e4:.1f}; q_b_D impact={qb[0]:.4e}")
        log(f"    [{tag}] Y_D impact={Y[0]:.4e} ({Y[0]/SS['Y_D']*100:+.4f}%SS)  "
            f"peakabs={pk(Y):.4e}")
        log(f"    [{tag}] K_D impact={K[0]:.4e} ({K[0]/SS['K_D']*100:+.4f}%SS)  "
            f"peakabs={pk(K):.4e} ({pk(K)/SS['K_D']*100:+.4f}%SS)")
        log(f"    [{tag}] n_inter_D impact={n[0]:.4e} ({n[0]/SS['n_inter_D']*100:+.4f}%SS); "
            f"b_D_D impact={b[0]:.4e} ({b[0]/SS['b_D_D']*100:+.4f}%SS); I_D impact={I[0]:.4e}")
        return d

    log("\n### Step 0/1 checks")
    solve(EL_anchor, "ELon")
    solve(0.0, "ELoff")
    log(f"\n- {ts()} SOLVE v2 COMPLETE. saved irfs_2p8_ELon.npz, irfs_2p8_ELoff.npz, ss_values.json")


if __name__ == "__main__":
    main()
