# finer continuation groups for the source of the Q_K rise (one-quarter exact clearing, Shapley)
import json, time
from math import factorial
import forensic3 as F3
from forensic_core import *
F3.GROUPS.clear()
F3.GROUPS.update(dict(Kpp=("Kp_D",), Hours=("N_D",), Bond=("Q_bD",), Fran=("alpha_D",), Cons=("C_D",), ToT=("p",),
                      Fside=("C_F", "N_F", "Q_bF", "Kp_F", "alpha_F")))
G = list(F3.GROUPS) + ["Today", "Buy"]
F3.GNAMES[:] = G
out = {}
for pd in (0.015, 0.06):
    So = shocked(F3.S_off, pd)
    vA = clear(rules_off, cal0, So); vDc = clear(rules_pol, cal_pol, So); share = vDc.m_cb / B
    eq = {}
    for mask in range(1 << len(G)):
        on = [g for i, g in enumerate(G) if mask >> i & 1]
        w = clear(rules_off, cal0, So, share=share if "Buy" in on else 0.0, cont=F3.Hybrid(on),
                  x0=(vDc.x if len(on) > 4 else vA.x).copy())
        eq[mask] = F3.summ(w)
    n = len(G)
    def shap(q):
        phi = [0.0] * n
        for i in range(n):
            for m in range(1 << n):
                if m >> i & 1: continue
                c = bin(m).count("1")
                phi[i] += factorial(c) * factorial(n - c - 1) / factorial(n) * (q(eq[m | 1 << i]) - q(eq[m]))
        return dict(zip(G, phi))
    out[str(pd)] = dict(QK=shap(lambda e: 100 * np.log(e["QK"])), ls=shap(lambda e: e["ls_bp"]),
                        rdep=shap(lambda e: e["rdep_bp"]), n=shap(lambda e: 100 * np.log(e["n"])),
                        solo={g: {k: eq[1 << i][k] - eq[0][k] for k in ("QK", "ls_bp", "rdep_bp")} for i, g in enumerate(G)},
                        maxerr=max(e["err"] for e in eq.values()), full=eq[(1 << n) - 1], base=eq[0])
    print(pd, "done", flush=True)
json.dump(out, open(OUT / "forensic" / "forensic3b_out.json", "w"), default=float, indent=1)
