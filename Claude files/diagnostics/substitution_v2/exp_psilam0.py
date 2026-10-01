"""
Experiment: output response at psi_lambda_B = 0 (fundamental floor, EL_price ON).

Solves one Jacobian at psi_lambda_B=0 (psi_spread=0, EL_price anchored), applies the
1pp default shock, and characterises the OUTPUT response — magnitude, I/C/NX
decomposition, and the deleveraging-vs-substitution regime (sign of capital LEVEL) —
against the calibrated psi_lambda_B=2.8 case (cached irfs_2p8_ELon.npz).

Diagnose only; imports unmodified code/, writes only to diagnostics/substitution_v2/.
"""
import os, sys, json, copy, datetime, textwrap
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "code"))
RUN_LOG = os.path.join(HERE, "run_log.md")
plt.rcParams["savefig.dpi"] = 200
BLUE, RED, GREEN = "#002147", "#8C1515", "#1a6e3a"
N = 60


def ts(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
def log(m=""):
    print(m, flush=True)
    with open(RUN_LOG, "a") as f: f.write(m + "\n")


def irf_to_dict(irf):
    out = {}
    for k in irf.keys():
        try:
            v = np.asarray(irf[k])
        except Exception:
            continue
        if v.ndim == 1: out[k] = v
    return out


def extremum(x, n=100):
    x = np.asarray(x)[:n]; i = int(np.argmax(np.abs(x))); return float(x[i]), i


def main():
    from calibration import get_calibration
    from steady_state import solve_steady_state
    from ic_delta_calibration import calibrate_ic_delta
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve, solve_jacobian_padded

    log(f"\n## Experiment — output response at psi_lambda_B = 0 (EL_price ON) — {ts()}")
    cal = get_calibration()
    ssr = calibrate_depreciation(calibrate_ic_delta(solve_steady_state(cal)))
    res = build_and_solve(ssr)
    ha, ss = res["ha_full"], res["ss_final"]
    unk, tgt, T, dshock = res["unknowns_tp"], res["targets_tp"], res["T"], res["dShock_def_D"]

    ss0 = copy.deepcopy(ss)
    ss0.toplevel["psi_lambda_B_D"] = 0.0; ss0.toplevel["psi_lambda_B_F"] = 0.0
    ss0.toplevel["psi_spread_D"]   = 0.0; ss0.toplevel["psi_spread_F"]   = 0.0
    # EL_price stays anchored (~0.1025)
    log(f"- {ts()} solving psi_lambda_B=0, psi_spread=0, EL_price={float(ss0['EL_price_D']):.6f}")
    G = solve_jacobian_padded(ha, ss0, unk, tgt,
                              ["Z_D", "shock_def_D", "Z_F", "shock_def_F"], T)
    d = irf_to_dict(G @ {"Z_D": np.zeros(T), "Z_F": np.zeros(T),
                         "shock_def_D": dshock, "shock_def_F": np.zeros(T)})
    np.savez(os.path.join(HERE, "irfs_psilam0_full.npz"), **d)

    SS = json.load(open(os.path.join(HERE, "ss_values.json")))
    C28 = np.load(os.path.join(HERE, "irfs_2p8_ELon.npz"))
    Yss, Kss, bss, nss = SS["Y_D"], SS["K_D"], SS["b_D_D"], SS["n_inter_D"]
    Pss, Css = SS["P_CES_D"], SS["C_D"]

    def gg(store, k):
        return store[k] if k in store.files else None

    Y, K, b, n, I, C = (d["Y_D"], d["K_D"], d["b_D_D"], d["n_inter_D"], d["I_D"], d["C_D"])
    NX = gg(d, "NX_D"); P = gg(d, "P_CES_D"); sp = d["spread_rb"]

    # ---- output response ----
    ye, yi = extremum(Y); ye28, _ = extremum(C28["Y_D"])
    log("\n### Output response")
    log(f"- spread peak = {sp[:100].max()*4*1e4:.1f}bp (fundamental floor)")
    log(f"- Y_D: impact {Y[0]/Yss*100:+.4f}%SS; trough {ye/Yss*100:+.4f}%SS at t={yi}")
    log(f"- vs psi_lambda_B=2.8: Y trough {ye28/Yss*100:+.4f}%SS  "
        f"(amplification 2.8/0 = {ye28/ye:.1f}x)")
    log(f"- SIGN: Y_D {'FALLS (contraction, no perverse rise)' if Y[:100].min()<0 else 'RISES — PERVERSE'}")

    # ---- decomposition ----
    contrib_I = I
    contrib_C = Pss * C + (Css * P if P is not None else 0.0)
    contrib_NX = NX if NX is not None else np.zeros_like(Y)
    log("\n### ΔY decomposition at psi_lambda_B=0 (impact)")
    for nm, x in [("I_D", contrib_I), ("P·C_D", contrib_C), ("NX_D", contrib_NX)]:
        e, ei = extremum(x)
        log(f"- {nm:7s}: impact {x[0]:+.4e}  extremum {e:+.4e} at t={ei}")

    # ---- regime: sign of capital LEVEL ----
    Ke, Ki = extremum(K); be, bi = extremum(b); ne, ni = extremum(n)
    log("\n### Regime at psi_lambda_B=0 (sign of capital LEVEL)")
    log(f"- K_D (capital LEVEL): impact {K[0]/Kss*100:+.4f}%SS, extremum {Ke/Kss*100:+.4f}%SS at t={Ki}")
    log(f"- b_D_D (sov quantity): impact {b[0]/bss*100:+.4f}%SS, extremum {be/bss*100:+.4f}%SS at t={bi}")
    log(f"- n_inter_D (net worth): impact {n[0]/nss*100:+.4f}%SS")
    Kmin, Kmax = K[:100].min(), K[:100].max()
    if Kmax <= 1e-12:
        reg = "DELEVERAGING-DOMINANT (capital level falls, never rises)"
    elif Kmin >= -1e-12:
        reg = "SUBSTITUTION-DOMINANT (capital level rises)"
    else:
        reg = (f"MIXED (K_D min {Kmin/Kss*100:+.3f}%SS, max {Kmax/Kss*100:+.3f}%SS)")
    log(f"- REGIME at psi_lambda_B=0: {reg}")
    log(f"- vs 2.8: K_D there fell to {extremum(C28['K_D'])[0]/Kss*100:+.4f}%SS (deleveraging-dominant)")

    # ---- figure ----
    t = np.arange(N)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].plot(t, C28["Y_D"][:N]/Yss*100, color=BLUE, lw=2, label="ψλ=2.8 (calibrated)")
    ax[0].plot(t, Y[:N]/Yss*100, color=RED, lw=2, label="ψλ=0 (fundamental floor)")
    ax[0].axhline(0, color="0.7", lw=0.6); ax[0].set_title("output Y_D"); ax[0].legend(fontsize=8)
    ax[1].plot(t, contrib_I[:N]/Yss*100, color=RED, label="I_D")
    ax[1].plot(t, contrib_C[:N]/Yss*100, color=BLUE, label="P·C_D")
    ax[1].plot(t, contrib_NX[:N]/Yss*100, color=GREEN, label="NX_D")
    ax[1].plot(t, Y[:N]/Yss*100, color="k", lw=1.8, label="ΔY_D")
    ax[1].axhline(0, color="0.7", lw=0.6); ax[1].set_title("ΔY decomposition @ ψλ=0"); ax[1].legend(fontsize=8)
    ax[2].plot(t, K[:N]/Kss*100, color=RED, lw=2, label="capital K_D (level)")
    ax[2].plot(t, b[:N]/bss*100, color=BLUE, lw=2, label="sov quantity b_D_D")
    ax[2].plot(t, n[:N]/nss*100, color=GREEN, lw=1.6, label="net worth n_inter_D")
    ax[2].axhline(0, color="0.7", lw=0.6); ax[2].set_title("balance sheet @ ψλ=0"); ax[2].legend(fontsize=8)
    for a in ax: a.set_xlabel("quarters")
    fig.suptitle("Output response at psi_lambda_B = 0 (EL_price fundamental floor)", fontsize=12)
    cap = ("At psi_lambda_B=0 the collateral-amplification is off, so only the fundamental "
           "expected-loss channel (EL_price) operates: output still falls (no perverse rise) but "
           "is ~an order of magnitude smaller than at the calibrated 2.8; the regime (sign of the "
           "capital LEVEL) shows whether deleveraging still dominates once the amplifier is removed.")
    chars = int(fig.get_size_inches()[0] * 14)
    fig.text(0.5, -0.04, textwrap.fill(cap, width=chars), ha="center", va="top",
             fontsize=8, style="italic", color="0.35")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(HERE, "exp_psilam0_output.png"), bbox_inches="tight")
    log(f"\n- {ts()} wrote irfs_psilam0_full.npz and exp_psilam0_output.png")


if __name__ == "__main__":
    main()
