"""
Diagnostic v2 analysis (reads cached IRFs; no solve).
Steps 1-4 + Step 6 note. Classifies regime by the SIGN OF THE CAPITAL LEVEL K_D.
Figures carry baked-in captions.
"""
import os, json, textwrap
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN_LOG = os.path.join(HERE, "run_log.md")
plt.rcParams["savefig.dpi"] = 200
BLUE, RED, GREEN, ORANGE = "#002147", "#8C1515", "#1a6e3a", "#c87941"
N = 60


def log(m=""):
    print(m, flush=True)
    with open(RUN_LOG, "a") as f:
        f.write(m + "\n")


ON  = np.load(os.path.join(HERE, "irfs_2p8_ELon.npz"))
OFF = np.load(os.path.join(HERE, "irfs_2p8_ELoff.npz"))
SS  = json.load(open(os.path.join(HERE, "ss_values.json")))


def g(store, k):
    return store[k] if k in store.files else None


def extremum(x, n=100):
    x = np.asarray(x)[:n]
    i = int(np.argmax(np.abs(x)))
    return float(x[i]), i


def cap_fig(fig, name, cap):
    chars = int(fig.get_size_inches()[0] * 14)
    fig.text(0.5, -0.03, textwrap.fill(cap, width=chars), ha="center", va="top",
             fontsize=8, style="italic", color="0.35")
    fig.savefig(os.path.join(HERE, name), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------- Step 1: aggregate sign ----------------
log("\n## Step 1 — aggregate output sign at psi_lambda_B=2.8 (ELon)")
Y = g(ON, "Y_D"); Yss = SS["Y_D"]
yext, yi = extremum(Y)
log(f"- Y_D impact = {Y[0]:.4e} ({Y[0]/Yss*100:+.4f}%SS); "
    f"extremum(100q) = {yext:.4e} ({yext/Yss*100:+.4f}%SS) at t={yi}")
log(f"- Sweep implied ΔY at 2.8 ≈ interpolate(2.6:-0.0271%, 3.0:-0.0389%) ≈ -0.033%SS. "
    f"Consistent: {'YES' if Y[:100].min() < 0 else 'NO — FLAG'}")

# ---------------- Step 2: bank balance sheet, classify by K_D level ----------------
log("\n## Step 2 — bank balance sheet (LEVELS), classify by sign of capital LEVEL K_D")
K = g(ON, "K_D"); bDD = g(ON, "b_D_D"); n = g(ON, "n_inter_D"); th = g(ON, "theta_D")
Q = g(ON, "Q_D"); qbD = g(ON, "q_b_D"); qbF = g(ON, "q_b_F"); bFD = g(ON, "b_F_D")
Kss, bss, nss, thss = SS["K_D"], SS["b_D_D"], SS["n_inter_D"], SS["theta_D"]
# total assets (market value) = theta*n  (k_balance_sheet identity); dTA = th_ss*dn + n_ss*dth
dTA = thss * n + nss * th
TAss = thss * nss
# bond market-value share of the book: s = q_b_D*b_D_D / (theta*n); first-order ds
mvbond_ss = SS["q_b_D"] * bss
dmvbond = SS["q_b_D"] * bDD + bss * qbD
d_share = (dmvbond - (mvbond_ss / TAss) * dTA) / TAss

Kext, Ki = extremum(K); bext, bi = extremum(bDD); next_, ni = extremum(n)
log(f"- K_D (capital LEVEL): impact {K[0]/Kss*100:+.4f}%SS, extremum {Kext/Kss*100:+.4f}%SS at t={Ki}")
log(f"- b_D_D (Greek sov QUANTITY): impact {bDD[0]/bss*100:+.4f}%SS, extremum {bext/bss*100:+.4f}%SS at t={bi}")
log(f"- n_inter_D (net worth): impact {n[0]/nss*100:+.4f}%SS, extremum {next_/nss*100:+.4f}%SS at t={ni}")
log(f"- total assets (θ·n): impact {dTA[0]/TAss*100:+.4f}%SS, min {dTA[:100].min()/TAss*100:+.4f}%SS")
log(f"- bond MV share of book: impact Δ {d_share[0]*100:+.4f}pp, extremum {extremum(d_share)[0]*100:+.4f}pp "
    f"(>0 = tilt TOWARD bonds)")
K_min = K[:100].min(); K_max = K[:100].max()
if K_max <= 1e-12:
    regime = "DELEVERAGING-DOMINANT: capital LEVEL falls (K_D ≤ 0 throughout)"
elif K_min >= -1e-12:
    regime = "SUBSTITUTION-DOMINANT (RED FLAG vs ΔY<0): capital LEVEL rises (K_D ≥ 0)"
else:
    regime = (f"MIXED: K_D crosses zero (min {K_min/Kss*100:+.3f}%SS, max {K_max/Kss*100:+.3f}%SS) "
              "— report both; classify by dominant/persistent sign")
log(f"- REGIME: {regime}")

# ---------------- Step 3a: Y decomposition ----------------
log("\n## Step 3a — ΔY_D decomposition (goods-market identity)")
I = g(ON, "I_D"); C = g(ON, "C_D"); NX = g(ON, "NX_D")
P = g(ON, "P_CES_D"); G = g(ON, "G_D"); Phi = g(ON, "Phi_D"); Tt = g(ON, "T_D")
Css, Pss = SS["C_D"], SS["P_CES_D"]
cons_bundle = (P if P is not None else 0*C) * Css + Pss * C if P is not None else Pss * C
# contribution of bundle consumption: d(P_CES*C) = P_ss*dC + C_ss*dP
contrib_C = Pss * C + (Css * P if P is not None else 0.0)
contrib_I = I
contrib_NX = NX if NX is not None else np.zeros_like(Y)
contrib_G = G if G is not None else np.zeros_like(Y)
for nm, x in [("I_D", contrib_I), ("P_CES·C_D", contrib_C), ("NX_D", contrib_NX), ("G_D", contrib_G)]:
    e, ei = extremum(x)
    log(f"- contribution {nm:10s}: impact {x[0]:+.4e}  extremum {e:+.4e} at t={ei}")
resid = Y - (contrib_I + contrib_C + contrib_NX + contrib_G
             + (Phi if Phi is not None else 0) + (Tt if Tt is not None else 0))
log(f"- identity residual max|·| = {np.max(np.abs(resid[:100])):.2e} (should be ~0)")

# ---------------- Step 3b: EL on/off capital decomposition ----------------
log("\n## Step 3b — capital: substitution push vs deleveraging pull (EL_price on/off)")
Kon = g(ON, "K_D"); Koff = g(OFF, "K_D")
deleverage = Kon - Koff            # what EL_price adds (expect < 0)
substitution = Koff                # residual, substitution-leaning (impure)
for label, x in [("net K_D (ELon)", Kon), ("substitution-leaning (ELoff)", substitution),
                 ("deleveraging pull (ELon-ELoff)", deleverage)]:
    e, ei = extremum(x)
    log(f"- {label:34s}: impact {x[0]/Kss*100:+.4f}%SS  extremum {e/Kss*100:+.4f}%SS at t={ei}")

# ---------------- Step 4: pass-through ----------------
log("\n## Step 4 — net-worth-to-spread pass-through at 2.8 (validation moment)")
sp = g(ON, "spread_rb")
spread_bp = float(sp[:100].max()) * 4 * 1e4
n_pct = float(n[:100].min()) / nss * 100
passthru = n_pct / (spread_bp / 100)
log(f"- peak spread = {spread_bp:.1f}bp; peak Δn_inter = {n_pct:+.3f}%SS; "
    f"pass-through = {passthru:+.3f}%/100bp  (sweep flag ≈ -0.85)")

# ---------------- Step 6: consistency note ----------------
log("\n## Step 6 — empirical-prior consistency")
credit_sign = "contracting" if (Kon[:100].min() < 0 and I[:100].min() < 0) else "NOT contracting"
bond_sign = "down" if bDD[:100].min() < 0 and bDD[:100].max() <= abs(bDD[:100].min()) else "mixed/up"
log(f"- capital/credit: {credit_sign}; Greek-sov quantity b_D_D: {bond_sign}; "
    f"bond MV share tilt: {'toward bonds' if extremum(d_share)[0] > 0 else 'toward capital'}")
log("- Model omits renationalisation motive (always-binding IC). Check it does NOT produce the "
    "doubly-counterfactual bonds↓/capital↑ WITH credit expanding.")
doubly_cf = (bDD[:100].min() < 0) and (Kon[:100].max() > 1e-9 and Kon[:100].min() >= -1e-12) and (I[:100].min() >= 0)
log(f"- doubly-counterfactual (bonds↓ & capital↑ & credit↑)? {'YES — PROBLEM' if doubly_cf else 'NO'}")

# ---------------- figures ----------------
t = np.arange(N)
# Fig A: Y and decomposition
figA, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(t, Y[:N]/Yss*100, color="k", lw=2.2, label="ΔY_D")
ax.plot(t, contrib_I[:N]/Yss*100, color=RED, lw=1.6, label="investment I_D")
ax.plot(t, contrib_C[:N]/Yss*100, color=BLUE, lw=1.6, label="bundle cons P·C_D")
ax.plot(t, contrib_NX[:N]/Yss*100, color=GREEN, lw=1.6, label="net exports NX_D")
if G is not None and np.any(np.abs(G[:N]) > 1e-12):
    ax.plot(t, contrib_G[:N]/Yss*100, color=ORANGE, lw=1.2, label="gov G_D")
ax.axhline(0, color="0.7", lw=0.6); ax.set_xlabel("quarters"); ax.set_ylabel("% of SS Y_D")
ax.legend(fontsize=8); ax.set_title("ΔY_D decomposition (psi_lambda_B=2.8, ELon)")
cap_fig(figA, "v2_01_Y_decomposition.png",
        "Output falls on the Greek default shock; the decomposition shows which of investment, "
        "consumption and net exports drives the negative sign — confirming (or not) that it is the "
        "investment/deleveraging channel rather than a terms-of-trade NX effect.")

# Fig B: balance sheet levels
figB, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(t, K[:N]/Kss*100, color=RED, lw=2, label="capital K_D (level)")
ax.plot(t, bDD[:N]/bss*100, color=BLUE, lw=2, label="Greek-sov quantity b_D_D")
ax.plot(t, dTA[:N]/TAss*100, color="k", lw=1.6, ls="--", label="total assets θ·n")
ax.plot(t, n[:N]/nss*100, color=GREEN, lw=1.6, label="net worth n_inter_D")
ax.axhline(0, color="0.7", lw=0.6); ax.set_xlabel("quarters"); ax.set_ylabel("% of SS")
ax.legend(fontsize=8); ax.set_title("Bank balance sheet — LEVELS (psi_lambda_B=2.8, ELon)")
cap_fig(figB, "v2_02_balance_sheet.png",
        "Decisive probe: the SIGN of the capital LEVEL K_D classifies the regime. Deleveraging "
        "dominates iff total assets, net worth, bonds AND capital all fall (bonds by more); a rising "
        "bond share with a falling capital level is still deleveraging-dominant.")

# Fig C: EL on/off capital overlay
figC, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(t, Kon[:N]/Kss*100, color=RED, lw=2, label="net K_D (EL_price ON)")
ax.plot(t, Koff[:N]/Kss*100, color=BLUE, lw=2, ls="--", label="EL_price OFF (substitution-leaning)")
ax.plot(t, deleverage[:N]/Kss*100, color=GREEN, lw=1.6, label="deleveraging pull (ON−OFF)")
ax.axhline(0, color="0.7", lw=0.6); ax.set_xlabel("quarters"); ax.set_ylabel("% of SS K_D")
ax.legend(fontsize=8); ax.set_title("Capital: substitution push vs deleveraging pull")
cap_fig(figC, "v2_03_EL_on_off_capital.png",
        "EL_price on/off decomposition of the capital response: EL_price OFF isolates the "
        "substitution-leaning residual (agency+IC channels only), the ON−OFF gap is the "
        "deleveraging pull the fundamental expected-loss/net-worth channel supplies.")

log(f"\n- analysis complete; figures v2_01/02/03 written.")
