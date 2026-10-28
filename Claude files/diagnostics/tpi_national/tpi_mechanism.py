# WHY A PURCHASE MOVES (OR DOES NOT MOVE) THE D BOND PRICE: THE RISK-TRANSFER FORMULA, MEASURED.
#
# At fixed continuation rules a purchase dm (paid with Z = q dm at the bond's lambda) leaves today's
# balance sheet and mu unchanged and moves the price only through next period's kernel Omega_j:
#   d ln q = sum_j w_j dOmega_j [ e_j + lam (1-mu)/E[Om] ] / (E[Om] R + lam mu)        (risk, franchise)
#          + sum_j w_j Omega_j dXi_j / E[Om Xi]                                       (continuation)
#   with e_j = Xi_j/q - R
# This script measures dOmega_j per regime from the solved rules and checks the formula against the
# exactly re-cleared price, for each economy, at rest and at the headline shock.
#   python3 tpi_mechanism.py nat uni nat_lz0
import json
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / __import__("os").environ.get("WT", "dff3")))
sys.path.insert(0, str(HERE))

import numpy as np

DM = 0.05 * 0.98  # 5% of the SS stock


def legs(v, cal):
    # PER-REGIME PIECES OF THE D-BANK D-BOND FOC.
    db, rec = cal["delta_b_D"], cal["recovery_rate_D"]
    W = [np.asarray(v.wgt[j]) for j in range(v.nreg)]
    Om = [np.asarray(v.Om_D[j]) for j in range(v.nreg)]
    Xi = [(rec if v.reg[j] else 1.0) * (db + (1 - db) * np.asarray(v.R[j]["Q_bD"]))
          for j in range(v.nreg)]
    return W, Om, Xi


def decompose(tag, rules, cal, ss, sproc, S):
    # FIXED RULES, FIXED OTHER UNKNOWNS: dOmega_j FROM A FORCED PURCHASE, AND THE FORMULA.
    import tpi_national as tn
    from global_projection.solver_recursive.recursive_experiment import read_exact
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive import tpi_experiment as te
    c = dict(cal, tpi_on=False)  # with the TPI off the last unknown is the purchase itself
    o = read_exact(rules, cal, ss, sproc, S)
    x0 = np.array(o["_x"], dtype=float)
    x0[-1] = 0.0
    v0 = tn.internals(S, x0, rules, c, ss, sproc)
    x1 = x0.copy()
    x1[-1] = DM / cal["B_gov_D_ss"]  # x_cb is the purchase share of the SS stock
    v1 = tn.internals(S, x1, rules, c, ss, sproc)
    W, Om0, Xi0 = legs(v0, cal)
    _, Om1, Xi1 = legs(v1, cal)
    q, R = float(v0.Q_bD), 1.0 + float(v0.rdep_D)
    lam = cal["lambda_bD_D"]
    EOm = float(sum(np.dot(W[j], Om0[j]) for j in range(v0.nreg)))
    den = EOm * R + lam * float(v0.mu_D)
    fr = lam * (1 - float(v0.mu_D)) / EOm
    # per-regime pieces: probability, mean dOmega, mean excess return
    EOmXi0 = float(sum(np.dot(W[j], Om0[j] * Xi0[j]) for j in range(v0.nreg)))
    rows, risk, franchise, cont = [], 0.0, 0.0, 0.0
    for j in range(v0.nreg):
        w, dOm, e = W[j], Om1[j] - Om0[j], Xi0[j] / q - R
        if w.sum() <= 0:
            continue
        risk += float(np.dot(w, dOm * e)) / den
        franchise += float(np.dot(w, dOm)) * fr / den
        # the third leg: tomorrow's price itself moves with the smaller private float
        cont += float(np.dot(w, Om0[j] * (Xi1[j] - Xi0[j]))) / EOmXi0
        f_ = cal["f_D"]
        al0, al1 = np.asarray(v0.R[j]["alpha_D"]), np.asarray(v1.R[j]["alpha_D"])
        La0, La1 = Om0[j] / (f_ + (1 - f_) * al0), Om1[j] / (f_ + (1 - f_) * al1)
        rd0, rd1 = np.asarray(v0.R[j]["rdep_D"]), np.asarray(v1.R[j]["rdep_D"])
        rows.append(dict(regime=int(v0.reg[j]), prob=float(w.sum()),
                         dalpha_rel=float(np.dot(w, al1 / al0 - 1) / w.sum()),
                         dLam_rel=float(np.dot(w, La1 / La0 - 1) / w.sum()),
                         drdep_bp=4e4 * float(np.dot(w, rd1 - rd0) / w.sum()),
                         Om=float(np.dot(w, Om0[j]) / w.sum()),
                         dOm=float(np.dot(w, dOm) / w.sum()),
                         dOm_rel=float(np.dot(w, dOm / Om0[j]) / w.sum()),
                         e_bp=4e4 * float(np.dot(w, e) / w.sum())))
    # the price the FOC implies at the perturbed continuation (other unknowns held)
    EOmXi1 = float(sum(np.dot(W[j], Om1[j] * Xi1[j]) for j in range(v1.nreg)))
    EOm1 = float(sum(np.dot(W[j], Om1[j]) for j in range(v1.nreg)))
    mu1 = 1.0 - EOm1 * R * float(v0.n_D) / float(v0.lev_D)
    q_implied = EOmXi1 / (EOm1 * R + lam * max(mu1, 0.0))
    # the fully re-cleared price at fixed rules (the other 13 unknowns re-solved)
    full = te.price_impact(rules, cal, ss, sproc, S, ms=(DM / cal["B_gov_D_ss"],))
    out = dict(tag=tag, q=q, mu=float(v0.mu_D), dlnq_formula_pct=100 * (risk + franchise + cont),
               risk_pct=100 * risk, franchise_pct=100 * franchise, cont_pct=100 * cont,
               dlnq_implied_pct=100 * (q_implied / q - 1),
               dlnq_full_pct=100 * (full[1]["Q_bD"] / full[0]["Q_bD"] - 1),
               dsov_full_bp=full[1]["sov_bp"] - full[0]["sov_bp"], regimes=rows,
               dn_next_nd=None)
    return out


def main():
    # EVERY RUN NAMED ON THE COMMAND LINE, AT REST AND AT THE SHOCK, IN ITS NO-TPI AND TPI RULES.
    import argparse
    import tpi_national as tn
    from global_projection.solver_recursive.recursive_experiment import (
        stochastic_rest_point, s_from_pd)
    from global_projection.solver_recursive.state_grid import IS
    res = {}
    for tag in sys.argv[1:]:
        meta = json.loads((HERE / "tpi_runs" / f"{tag}.json").read_text())["args"]
        a = argparse.Namespace(**{k: meta.get(k) for k in ("market", "n_jobs", "cap", "lamZ", "phi", "keyD", "theta", "warranted", "keyLoss")})
        cal, ss, sproc = tn.setup(a)
        rr = pickle.load(open(HERE / "tpi_runs" / f"{tag}_rules.pkl", "rb"))
        for name in (("off",) if len(sys.argv) > 1 else ("off", "on")):
            c = dict(cal, tpi_on=(name == "on"))
            S0 = stochastic_rest_point(rr[name], c, ss, sproc, verbose=False)
            S1 = S0.copy()
            S1[IS] = s_from_pd(meta["pd"])
            for where, S in (("rest", S0), ("shock", S1)):
                d = decompose(tag, rr[name], c, ss, sproc, S)
                res[f"{tag}/{name}/{where}"] = d
                print(f"  {tag:8s} {name:3s} {where:5s}  q={d['q']:.4f} mu={d['mu']:.4f}  "
                      f"dlnq per 5% of stock: formula {d['dlnq_formula_pct']:+.4f}% "
                      f"(risk {d['risk_pct']:+.4f}, franchise {d['franchise_pct']:+.4f}, "
                      f"continuation {d['cont_pct']:+.4f}), "
                      f"implied {d['dlnq_implied_pct']:+.4f}%, full {d['dlnq_full_pct']:+.4f}% "
                      f"= {d['dsov_full_bp']:+.2f} bp", flush=True)
                for r in d["regimes"]:
                    print(f"        d'={r['regime']} prob {r['prob']:.4f}  Om {r['Om']:.4f}  "
                          f"dOm/Om {100*r['dOm_rel']:+.4f}%  (dLam {100*r['dLam_rel']:+.4f}%, "
                          f"dalpha {100*r['dalpha_rel']:+.4f}%, drdep' {r['drdep_bp']:+.2f} bp)  "
                          f"excess return {r['e_bp']:+.0f} bp/yr")
    (HERE / "tpi_runs" / "mechanism.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
