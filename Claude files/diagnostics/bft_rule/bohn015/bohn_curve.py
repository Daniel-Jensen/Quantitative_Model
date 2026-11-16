# THE BINDING THRESHOLD AND THE 1.5% IMPACT ALONG THE SAVED, CONVERGED CONTINUATION IN THE BOHN COEFFICIENT
import json
from bohn_compare import Econ, RUNS, OUT, clear, fitted, shocked, summ, threshold, RE
ROOT = RUNS.parents[3]
files = {0.06: ROOT / "data" / "global_warm_start_bohn006.pkl", 0.10: RUNS / "rules_12pct_bohn0.100.pkl",
         0.15: RUNS / "rules_12pct_bohn0.150.pkl", 0.20: RUNS / "rules_12pct_bohn0.200.pkl",
         0.25: RUNS / "rules_12pct_bohn0.250.pkl"}
out = {}
for b, fpath in files.items():
    e = Econ(b, fpath)
    S = RE.stochastic_rest_point(e.rules, e.cal, e.ss, e.sp, verbose=False)
    rest = summ(e, clear(e, S))
    thr = threshold(e, S, "exact", hi=0.12)
    imp = {pd: summ(e, clear(e, shocked(S, pd))) for pd in (0.015, 0.03, 0.06)}
    out[str(b)] = dict(root=1 - e.cal["delta_b_D"] + (e.cal["delta_b_D"] - b / 4) / e.ss["gs_D"]["Q_B_ss"],
                       thr=thr["pd"], rest=rest, imp=imp)
    print(b, f"root {out[str(b)]['root']:.4f} thr {thr['pd']}", f"rest n {rest['n']:.4f} b_DD {rest['b_DD']:.4f} debt {rest['Bp']:.4f} psi {rest['psi']:.4f} EOm {rest['EOm']:.4f}",
          {pd: (round(v['ls_bp'], 1), round(100 * (v['Y'] / rest['Y'] - 1), 4)) for pd, v in imp.items()}, flush=True)
json.dump(out, open(OUT / "bohn_curve.json", "w"), indent=1)
