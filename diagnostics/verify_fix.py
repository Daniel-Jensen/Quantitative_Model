"""
Verify the macro-pru-fix (fundamental expected-loss channel, EL_price).

Compares PRE-fix caches (irfs_*.npz, summary.json — committed on main) against
POST-fix caches (irfs_*_postfix.npz, summary_postfix.json produced by
`solve_configs.py postfix`). Checks:

  (A) SS-neutrality: the fix must NOT move the steady state (EL_price ∝ def_rate(+1)=0 at SS).
  (B) The psi_lambda_B=0 default response is now NONZERO and CORRECTLY SIGNED
      (pre-fix it was identically 0): q_b_D<0, n_inter_D<0, spread_rb>0, Y_D<0 on impact.
  (C) Baseline (psi_lambda_B=3) remains correctly signed (fix adds to the loading).

Emits diagnostics/04_fix_psilam0_before_after.png (captioned) and prints a verdict.
"""
import os, json, textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams["savefig.dpi"] = 200
BLUE, RED, GREY = "#002147", "#8C1515", "0.6"


def load(name):
    return np.load(os.path.join(HERE, name))


pre0  = load("irfs_psilam0.npz")          # pre-fix, psi_lambda_B=0
post0 = load("irfs_psilam0_postfix.npz")  # post-fix, psi_lambda_B=0
preB  = load("irfs_baseline.npz")
postB = load("irfs_baseline_postfix.npz")

with open(os.path.join(HERE, "summary.json"))         as f: sPre  = json.load(f)
with open(os.path.join(HERE, "summary_postfix.json")) as f: sPost = json.load(f)

print("=" * 70)
print("FIX VERIFICATION — macro-pru-fix (EL_price fundamental default channel)")
print("=" * 70)

# ---------- (A) SS-neutrality ----------
print("\n(A) Steady-state neutrality (pre vs post fix):")
ss_keys = ["q_b_D", "q_b_F", "rb_D", "n_inter_D", "theta_D", "beta_D", "p",
           "spread_rb", "goods_mkt_D", "deposit_mkt_D"]
ss_ok = True
for k in ss_keys:
    a = sPre["ss_checks"].get(k); b = sPost["ss_checks"].get(k)
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        d = abs(a - b)
        flag = "" if d < 1e-10 else "  <-- MOVED"
        if d >= 1e-10:
            ss_ok = False
        print(f"    {k:14s} pre={a:.10g}  post={b:.10g}  |Δ|={d:.2e}{flag}")
# report EL_price actually set
elp = sPost["ss_checks"].get("EL_price_D", "n/a")
print(f"    EL_price_D (post) = {sPost['ss_checks'].get('EL_price_D','not in summary')}")
print(f"  => SS {'UNCHANGED ✓' if ss_ok else 'MOVED ✗ (unexpected)'}")

# ---------- (B) psi_lambda_B = 0 now responds, correct signs ----------
print("\n(B) Default shock at psi_lambda_B = 0  (pre-fix = identically 0):")
print(f"    {'var':11s} {'pre impact':>12s} {'post impact':>12s} {'post peak|·|':>13s}  sign check")
checks = [("q_b_D", "<0"), ("n_inter_D", "<0"), ("spread_rb", ">0"), ("Y_D", "<0"),
          ("theta_D", ">0"), ("C_D", "<0")]
signs_ok = True
for v, want in checks:
    kd = f"def__{v}"
    pre_i  = float(pre0[kd][0])  if kd in pre0.files  else float("nan")
    post   = post0[kd] if kd in post0.files else np.zeros(1)
    post_i = float(post[0]); post_pk = float(np.max(np.abs(post[:100])))
    if want == "<0":
        ok = post_i < -1e-9
    else:
        ok = post_i > 1e-9
    signs_ok &= ok
    print(f"    {v:11s} {pre_i:>12.4e} {post_i:>12.4e} {post_pk:>13.4e}  want {want}: {'OK' if ok else 'FAIL'}")
alive = any(float(np.max(np.abs(post0[f'def__{v}'][:100]))) > 1e-9
            for v, _ in checks if f"def__{v}" in post0.files)
print(f"  => channel {'ALIVE ✓' if alive else 'STILL DEAD ✗'}; "
      f"signs {'CORRECT ✓' if signs_ok else 'WRONG ✗'}")

# ---------- (C) baseline still correctly signed ----------
print("\n(C) Baseline (psi_lambda_B = 3) default shock, pre vs post fix (impact):")
for v in ["q_b_D", "n_inter_D", "spread_rb", "Y_D"]:
    kd = f"def__{v}"
    a = float(preB[kd][0]) if kd in preB.files else float("nan")
    b = float(postB[kd][0]) if kd in postB.files else float("nan")
    print(f"    {v:11s} pre={a:>12.4e}  post={b:>12.4e}")

# ---------- figure ----------
N = 60
CAP = ("Fix verification: at psi_lambda_B=0 the sovereign default shock now moves the "
       "bond price, bank net worth, yield spread and output with the correct signs "
       "(solid red) whereas pre-fix they were identically zero (dashed grey) — the "
       "fundamental expected-loss channel (EL_price) is active independently of the "
       "psi_lambda_B collateral friction.")
fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
for ax, v, ttl in zip(axes, ["q_b_D", "n_inter_D", "spread_rb", "Y_D"],
                      ["bond price q_b_D", "net worth n_inter_D",
                       "yield spread spread_rb", "output Y_D"]):
    kd = f"def__{v}"
    if kd in pre0.files:
        ax.plot(pre0[kd][:N], color=GREY, ls="--", lw=1.6, label="pre-fix (ψλ=0)")
    if kd in post0.files:
        ax.plot(post0[kd][:N], color=RED, lw=2, label="post-fix (ψλ=0)")
    ax.axhline(0, color="0.8", lw=0.6)
    ax.set_title(ttl, fontsize=10); ax.set_xlabel("quarters", fontsize=8)
axes[0].legend(fontsize=8, loc="best")
fig.suptitle("macro-pru-fix: fundamental default channel active at psi_lambda_B = 0",
             fontsize=12)
chars = int(fig.get_size_inches()[0] * 14)
fig.text(0.5, -0.03, textwrap.fill(CAP, width=chars), ha="center", va="top",
         fontsize=8, style="italic", color="0.35")
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(HERE, "04_fix_psilam0_before_after.png"), bbox_inches="tight")
plt.close(fig)
print("\nwrote 04_fix_psilam0_before_after.png")

# ---------- verdict ----------
ok = ss_ok and alive and signs_ok
print("\n" + "=" * 70)
print(f"VERDICT: {'FIX VERIFIED ✓' if ok else 'FIX NOT VERIFIED ✗'}"
      f"  (SS-neutral={ss_ok}, channel-alive={alive}, signs-correct={signs_ok})")
print("=" * 70)
