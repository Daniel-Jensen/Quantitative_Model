# WHY THE THRESHOLD MOVES: AT A COMMON STATE (THE 0.25 REST POINT, SHOCKED), EACH ECONOMY'S EXACT CLEARING, WITH THE
# BANKER'S VALUATION OF NET WORTH E[Omega] = E[Lam (f + (1-f) alpha')] SPLIT INTO THE HOUSEHOLD DISCOUNT FACTOR AND THE
# FRANCHISE VALUE, AND THE NEXT-QUARTER FISCAL AND BOND-MARKET OBJECTS BEHIND alpha'.
import json, sys
import numpy as np
from bohn_compare import Econ, RUNS, OUT, clear, shocked, summ, PM, RE
E = {0.25: Econ(0.25, RUNS / "rules_12pct_bohn0.250.pkl"), 0.15: Econ(0.15, RUNS / "rules_12pct_bohn0.150.pkl")}
S25 = RE.stochastic_rest_point(E[0.25].rules, E[0.25].cal, E[0.25].ss, E[0.25].sp, verbose=False)
out = {}
for pd in (0.012, 0.015, 0.03):
    S = shocked(S25, pd)
    for b, e in E.items():
        v = clear(e, S)
        cal, f = e.cal, e.cal["f_D"]
        xc = PM._ghh(v.C_D, v.N_D, cal["chi_D"], cal["frisch_D"], v.X_FLOOR_D, v.X_EPS)
        lam = [cal["beta_inter_D"] * (xc / v.XN_D[j]) ** cal["sigma_D"] for j in range(v.nreg)]
        alp = [np.real(v.R[j]["alpha_D"]) for j in range(v.nreg)]
        Elam = float(PM._expect(v.wgt, lam)); Ealp = float(PM._expect(v.wgt, alp))
        Ealp_nd = float(np.mean(alp[0])) if np.ndim(alp[0]) else float(alp[0])
        Qn = [np.real(v.R[j]["Q_bD"]) for j in range(v.nreg)]
        # the quarter after: debt and taxes at the next no-default quadrature states, from the exact next-state clearing
        eps, _ = PM.gh_nodes(e.ngh); Sn, _ = PM._next_states(v, e.sp, eps)
        nxt = [summ(e, clear(e, Sn[i], d=0)) for i in range(Sn.shape[0])]
        s = summ(e, v)
        out[f"{pd}_{b}"] = dict(psi=s["psi"], EOm=s["EOm"], Elam=Elam, Ealpha=Ealp, Ealpha_nodef=Ealp_nd,
                                EOm_rebuilt=float(PM._expect(v.wgt, [lam[j] * (f + (1 - f) * alp[j]) for j in range(v.nreg)])),
                                R=s["R"], n=s["n"], lev=s["lev"], Tax=s["Tax"], new=s["new"], Bp=s["Bp"], C=s["C"],
                                er_bp=float(4e4 * v.er_D), EQ_next=float(PM._expect(v.wgt, Qn)),
                                next_Bp=float(np.mean([x["Bp"] for x in nxt])), next_Tax=float(np.mean([x["Tax"] for x in nxt])),
                                next_mu=float(np.mean([x["mu"] for x in nxt])), next_er=None,
                                next_alpha=float(np.mean([x["alpha"] for x in nxt])), next_n=float(np.mean([x["n"] for x in nxt])),
                                next_ls=float(np.mean([x["ls_bp"] for x in nxt])))
json.dump(out, open(OUT / "bohn_why.json", "w"), indent=1)
for k, d in out.items():
    print(k, {q: (round(x, 5) if isinstance(x, float) else x) for q, x in d.items()})
