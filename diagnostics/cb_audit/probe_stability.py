"""CB-block audit — Step 4 analysis (feedback sign + closed-loop stability).

Reads diagnostics/cb_audit/probe_pipeline.npz (written by probe_pipeline.py) and
reports:
  * the sign and shape of A_cb = d(spread_rb)/d(cb_buy_D);
  * the closed-loop pole and the safety margin against the intended gamma range;
  * per-gamma stability: b_gov_D[499], and the dominant modulus of the spread,
    b_gov_D and n_inter_D tails via the validated Prony estimator;
  * the CB conduit's on-path fiscal incidence (rem_cb_D / rem_cb_F signs).

Audit-only; writes probe_stability.md.
"""
import os, json
import numpy as np
from prony import prony_modulus

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, "probe_pipeline.npz"))
J = json.load(open(os.path.join(HERE, "probe_pipeline.json")))

L = []
def p(s=""):
    print(s)
    L.append(s)

A_cb = d["A_cb"]
A_def = d["A_def"]
T = A_cb.shape[0]

p("# Step 4 — CB feedback sign and closed-loop stability\n")
p(f"A_cb = d(spread_rb)/d(cb_buy_D), T = {T}\n")
p("## 4a. Sign of the feedback\n")
p("| lag h | A_cb[h,0]  (impulse at t=0) | A_cb[h,h] (diagonal) |")
p("|---|---|---|")
for h in range(8):
    p(f"| {h} | {A_cb[h,0]:+.6e} | {A_cb[h,h]:+.6e} |")
neg0 = A_cb[0, 0] < 0
colsum = A_cb[:, 0].sum()
p("")
p(f"- `A_cb[0,0] = {A_cb[0,0]:+.6e}`  -> impact effect of a unit purchase on the "
  f"impact spread is **{'NEGATIVE (compresses)' if neg0 else 'POSITIVE (WIDENS — FAIL)'}**")
p(f"- column sum `sum_h A_cb[h,0] = {colsum:+.6e}`  (cumulative spread response to a "
  f"one-period purchase at t=0)")
p(f"- fraction of the t=0 column that is negative: "
  f"{100.0*np.mean(A_cb[:,0] < 0):.1f}%")
p(f"- most positive entry in the t=0 column: {A_cb[:,0].max():+.3e} at lag "
  f"{int(np.argmax(A_cb[:,0]))}")
ev = np.linalg.eigvals(A_cb)
p(f"- spectral radius of A_cb = {np.max(np.abs(ev)):.6f}; "
  f"max real eigenvalue = {np.max(ev.real):+.6f}; "
  f"min real eigenvalue = {np.min(ev.real):+.6f}")
p(f"- 1 / max real eigenvalue = {1.0/np.max(ev.real):+.4f} "
  f"(the gamma at which I - gamma*A_cb becomes singular, if the max real "
  f"eigenvalue is the binding one)")

p("\n## 4b. Closed-loop pole and margin\n")
pole = J.get("closed_loop_pole")
p(f"- condition-number scan: {pole}")
p(f"- intended gamma range: [0, 10] (code/tpi.py gamma_values = [0, 2, 5, 10])")
if pole:
    p(f"- margin: pole / gamma_max = {pole['gamma_pole']/10.0:.2f}x; "
      f"the 0.75-safety cap sits at gamma = {pole['gamma_safe_max']:.2f}")

p("\n## 4c. Per-gamma stability\n")
p("| gamma | peak spread (bp ann) | b_gov_D[499] | n_inter_D[0] | Y_D[0] | "
  "|lam| spread | |lam| b_gov_D | |lam| n_inter_D |")
p("|---|---|---|---|---|---|---|---|")
for g in (0, 2, 5, 10):
    sp = d[f"irf_g{g}_spread_rb"]
    bg = d[f"irf_g{g}_b_gov_D"]
    ni = d[f"irf_g{g}_n_inter_D"]
    yd = d[f"irf_g{g}_Y_D"]
    m_sp = prony_modulus(sp[:400])[0]
    m_bg = prony_modulus(bg[:400])[0]
    m_ni = prony_modulus(ni[:400])[0]
    p(f"| {g} | {J[f'gamma{g}.peak_spread_bp']:.2f} | "
      f"{J[f'gamma{g}.b_gov_D_499']:+.3e} | {ni[0]:+.4e} | {yd[0]:+.4e} | "
      f"{m_sp:.6f} | {m_bg:.6f} | {m_ni:.6f} |")
p("")
p("All moduli must be < 1 for a stationary closed loop.")

p("\n## 4d. Breakdown scan over gamma (open-loop grid, up to the pole)\n")
eps = d["irf_g10_cb_buy_D"] * 0  # placeholder shape
# reconstruct the shock from the gamma=0 spread: spread_0 = A_def @ eps
eps = np.linalg.lstsq(A_def, d["irf_g0_spread_rb"], rcond=None)[0]
p("| gamma | peak spread (bp ann) | compression vs g=0 | cond(I - g A_cb) | |lam| spread |")
p("|---|---|---|---|---|")
base = None
I = np.eye(T)
for g in [0, 1, 2, 5, 10, 15, 19.88, 22, 25, 26.5]:
    M = I - g * A_cb
    c = np.linalg.cond(M)
    sp = np.linalg.solve(M, A_def @ eps)
    pk = sp[:100].max() * 4e4
    if base is None:
        base = pk
    m = prony_modulus(sp[:400])[0]
    p(f"| {g:g} | {pk:.2f} | {100*(1-pk/base):+.1f}% | {c:.3e} | {m:.6f} |")

p("\n## 4e. On-path fiscal incidence of the conduit\n")
have = [k for k in d.files if "rem_cb" in k or "cb_flow" in k]
p(f"conduit series present in the dump: {sorted(have)}")
p(f"G columns that do not exist (pure-CB objects have no shock_def_D loading): "
  f"{J.get('G_missing_columns')}")
for g in (2, 5, 10):
    for nm in ("rem_cb_D", "rem_cb_F", "cb_flow_D"):
        k = f"irf_g{g}_{nm}"
        if k in d.files:
            x = d[k]
            p(f"- g={g:<2} {nm:10s}: t0 = {x[0]:+.6e}, t1 = {x[1]:+.6e}, "
              f"min = {x.min():+.6e}, max = {x.max():+.6e}, "
              f"undiscounted sum(0:100) = {x[:100].sum():+.6e}")
for g in (2, 5, 10):
    for nm in ("TAX_D", "TAX_F"):
        k = f"irf_g{g}_{nm}"
        if k in d.files:
            x = d[k]
            p(f"- g={g:<2} {nm:10s}: t0 = {x[0]:+.6e}, peak|.| = {np.abs(x).max():+.6e}")

p("\n## 4f. Walras residuals along the closed loop\n")
p("| gamma | max|ca_res_D| | max|goods_mkt_D| | max|goods_mkt_F| |")
p("|---|---|---|---|")
for g in (0, 2, 5, 10):
    p(f"| {g} | {J[f'gamma{g}.max_ca_res_D']:.2e} | "
      f"{J[f'gamma{g}.max_goods_mkt_D']:.2e} | {J[f'gamma{g}.max_goods_mkt_F']:.2e} |")

open(os.path.join(HERE, "probe_stability.md"), "w").write("\n".join(L) + "\n")
print("\nWROTE probe_stability.md")
