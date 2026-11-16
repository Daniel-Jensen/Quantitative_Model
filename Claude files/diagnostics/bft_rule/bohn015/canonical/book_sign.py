# WHERE THE EUROSYSTEM BOOK CAN TURN NEGATIVE (2026-11-13): AT THE GRID NODES (BY THE BOOK CARRIED IN), ALONG THE EXACT
# IMPULSE RESPONSES, AND ON THE SIMULATED ERGODIC SET OF THE FITTED RULES; AND THE OFF-GRID ERRORS ON THE ECONOMICALLY
# ADMISSIBLE PART OF THAT SET (BOOK >= 0). NO SOLVE.
import json, pickle, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE.parents[1])); sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "k2"))
from results_io import load
from global_projection.solver_recursive.state_grid import IM, IS, default_prob
from global_projection.solver_recursive import accuracy as AC
from run_bft import node_checks
import headline as H

d = json.load(open(H.OUT / "headline.json"))
cal = load(H.DATA / "calibration"); ss = load(H.DATA / "steady_state"); cal["tpi_on"] = False
info = json.load(open(H.OUT / "solve_k2_info.json"))
cal_k = dict(cal, tpi_on=True, tpi_rule="bft", tpi_gz=False, tpi_bft_er_ref=info["er_ref_bp"] / 4e4, tpi_bft_k=2.0)
B = cal["B_gov_D_ss"]
out = {}
for tag, c, rf in (("canonical_0.15", cal_k, H.OUT / "rules_bft_k2.pkl"),
                   ("pr38_0.25", None, H.ROOT.parent / "QUANTITATIVE_MODEL" / "results" / "GLOBAL_2026-11-12_pr38_bft" / "rules_bft_k2.pkl")):
    if c is None:
        import compare_k2 as C
        off = C.Econ(0.25, C.RUNS / "rules_12pct_bohn0.250.pkl")
        E = C.policy(off, rf, 15.376861746491961)
        S0 = np.array(json.load(open(H.DIAG / "compare_k2.json"))["0.25_pol"]["rest"]["state"])
    else:
        E = H.Eco(c, ss, pickle.load(open(rf, "rb")))
        S0 = np.array(d["pol"]["rest"]["state"])
    rows = node_checks(E.rules, E.cal, E.ss, E.sp)["rows"]
    nonneg_in = rows["M_lag"] >= -1e-12
    o = dict(nodes=int(len(rows["M_lag"])), nodes_negative_book_in=int((~nonneg_in).sum()),
             min_book_out_given_nonneg_in_pct=float(100 * rows["M_new"][nonneg_in].min()),
             min_purchase_given_nonneg_in_pct=float(100 * rows["m"][nonneg_in].min()))
    states, _ = AC.simulate(E.rules, E.cal, E.ss, E.sp, S0, T=1200, burn=200, seed=0)
    neg = states[:, IM] < -1e-10
    first = int(np.argmax(neg)) if neg.any() else None
    err = AC.euler_errors(E.rules, E.cal, E.ss, E.sp, states)
    ok = np.isfinite(err).all(axis=1)
    l10 = np.log10(np.maximum(err, 1e-16))
    adm = ok & ~neg
    o.update(sim_share_negative=float(neg.mean()), sim_first_negative=first,
             sim_book_min_pct=float(100 * states[:, IM].min() / B),
             pd_at_first_negative=float(default_prob(states[first, IS])) if first is not None else None,
             l10_all=(float(l10[ok].mean()), float(l10[ok].max())),
             l10_admissible=(float(l10[adm].mean()), float(l10[adm].max())) if adm.any() else None,
             n_admissible=int(adm.sum()))
    out[tag] = o
    print(tag, o, flush=True)
for arm in ("off", "pol"):
    for pd in ("0.02", "0.015"):
        print(arm, pd, "IRF book min (% of debt)", min(d[arm]["irf"][pd]["M_cb"]), "purchase min", min(d[arm]["irf"][pd]["m_cb"]))
json.dump(out, open(H.OUT / "book_sign.json", "w"), indent=1)
