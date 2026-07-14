"""
Analysis of the psi_lambda_B=0 output experiment (reads cached npz; no solve).
Fixes the .files-on-dict bug in exp_psilam0.py by loading the saved npz.
"""
import os, json, textwrap
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN_LOG = os.path.join(HERE, "run_log.md")
plt.rcParams["savefig.dpi"] = 200
BLUE, RED, GREEN = "#002147", "#8C1515", "#1a6e3a"
N = 60


def log(m=""):
    print(m, flush=True)
    with open(RUN_LOG, "a") as f:
        f.write(m + "\n")


def gg(store, k):
    return store[k] if k in store.files else None


def extremum(x, n=100):
    x = np.asarray(x)[:n]; i = int(np.argmax(np.abs(x))); return float(x[i]), i


F0  = np.load(os.path.join(HERE, "irfs_psilam0_full.npz"))
F28 = np.load(os.path.join(HERE, "irfs_2p8_ELon.npz"))
SS  = json.load(open(os.path.join(HERE, "ss_values.json")))
Yss, Kss, bss, nss = SS["Y_D"], SS["K_D"], SS["b_D_D"], SS["n_inter_D"]
Pss, Css = SS["P_CES_D"], SS["C_D"]

Y, K, b, n = F0["Y_D"], F0["K_D"], F0["b_D_D"], F0["n_inter_D"]
I, C, sp = F0["I_D"], F0["C_D"], F0["spread_rb"]
NX, P = gg(F0, "NX_D"), gg(F0, "P_CES_D")

log("\n### Output response (fixed analysis)")
ye, yi = extremum(Y); ye28, _ = extremum(F28["Y_D"])
log(f"- spread peak = {sp[:100].max()*4*1e4:.1f}bp (fundamental floor, EL_price only)")
log(f"- Y_D: impact {Y[0]/Yss*100:+.4f}%SS; trough {ye/Yss*100:+.4f}%SS at t={yi}")
log(f"- vs psi_lambda_B=2.8: Y trough {ye28/Yss*100:+.4f}%SS (amplification 2.8/0 = {ye28/ye:.1f}x)")
log(f"- SIGN: Y_D {'FALLS (contraction, no perverse rise)' if Y[:100].min()<0 else 'RISES — PERVERSE'}")

log("\n### ΔY decomposition at psi_lambda_B=0 (impact)")
cI = I
cC = Pss * C + (Css * P if P is not None else 0.0)
cNX = NX if NX is not None else np.zeros_like(Y)
for nm, x in [("I_D", cI), ("P·C_D", cC), ("NX_D", cNX)]:
    e, ei = extremum(x)
    log(f"- {nm:7s}: impact {x[0]:+.4e}  extremum {e:+.4e} at t={ei}")
resid = Y - (cI + cC + cNX)
log(f"- identity residual max|·| = {np.max(np.abs(resid[:100])):.2e}")

log("\n### Regime at psi_lambda_B=0 (sign of capital LEVEL)")
Ke, Ki = extremum(K); be, bi = extremum(b)
log(f"- K_D (capital LEVEL): impact {K[0]/Kss*100:+.4f}%SS, extremum {Ke/Kss*100:+.4f}%SS at t={Ki}")
log(f"- b_D_D (sov quantity): impact {b[0]/bss*100:+.4f}%SS, extremum {be/bss*100:+.4f}%SS at t={bi}")
log(f"- n_inter_D (net worth): impact {n[0]/nss*100:+.4f}%SS")
Kmin, Kmax = K[:100].min(), K[:100].max()
if Kmax <= 1e-12:
    reg = "DELEVERAGING-DOMINANT (capital level falls, never rises)"
elif Kmin >= -1e-12:
    reg = "SUBSTITUTION-DOMINANT (capital level rises)"
else:
    reg = f"MIXED (K_D min {Kmin/Kss*100:+.4f}%SS, max {Kmax/Kss*100:+.4f}%SS)"
log(f"- REGIME at psi_lambda_B=0: {reg}")
log(f"- vs 2.8: K_D fell to {extremum(F28['K_D'])[0]/Kss*100:+.4f}%SS (deleveraging-dominant)")

# figure
t = np.arange(N)
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
ax[0].plot(t, F28["Y_D"][:N]/Yss*100, color=BLUE, lw=2, label="ψλ=2.8 (calibrated, 147bp)")
ax[0].plot(t, Y[:N]/Yss*100, color=RED, lw=2, label="ψλ=0 (floor, 13bp)")
ax[0].axhline(0, color="0.7", lw=0.6); ax[0].set_title("output Y_D (%SS)"); ax[0].legend(fontsize=8)
ax[1].plot(t, cI[:N]/Yss*100, color=RED, label="I_D")
ax[1].plot(t, cC[:N]/Yss*100, color=BLUE, label="P·C_D")
ax[1].plot(t, cNX[:N]/Yss*100, color=GREEN, label="NX_D")
ax[1].plot(t, Y[:N]/Yss*100, color="k", lw=1.8, label="ΔY_D")
ax[1].axhline(0, color="0.7", lw=0.6); ax[1].set_title("ΔY decomposition @ ψλ=0"); ax[1].legend(fontsize=8)
ax[2].plot(t, K[:N]/Kss*100, color=RED, lw=2, label="capital K_D (level)")
ax[2].plot(t, b[:N]/bss*100, color=BLUE, lw=2, label="sov quantity b_D_D")
ax[2].plot(t, n[:N]/nss*100, color=GREEN, lw=1.6, label="net worth n_inter_D")
ax[2].axhline(0, color="0.7", lw=0.6); ax[2].set_title("balance sheet @ ψλ=0 (%SS)"); ax[2].legend(fontsize=8)
for a in ax:
    a.set_xlabel("quarters")
fig.suptitle("Output response at psi_lambda_B = 0 (EL_price fundamental floor)", fontsize=12)
cap = ("At psi_lambda_B=0 the collateral amplification is off, so only the fundamental "
       "expected-loss channel (EL_price) operates: output still falls (no perverse rise) but is "
       f"~{extremum(F28['Y_D'])[0]/extremum(Y)[0]:.0f}x smaller than at the calibrated 2.8; the sign "
       "of the capital LEVEL shows whether deleveraging still dominates once the amplifier is removed.")
chars = int(fig.get_size_inches()[0] * 14)
fig.text(0.5, -0.04, textwrap.fill(cap, width=chars), ha="center", va="top",
         fontsize=8, style="italic", color="0.35")
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(HERE, "exp_psilam0_output.png"), bbox_inches="tight")
log("\n- wrote exp_psilam0_output.png (analysis from cached irfs_psilam0_full.npz)")
