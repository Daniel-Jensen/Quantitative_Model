# NO POLICY AGAINST THE ANNOUNCED k = 2 PROPORTIONAL RULE, AT BOHN 0.25 AND 0.15, AT 1.5% AND 2% (AND 6%, EXACT ONLY).
# NO GLOBAL SOLVE: THE FOUR SAVED SOLUTIONS ARE CLEARED EXACTLY (ONE QUARTER), JOINTLY OVER TWO QUARTERS, AND READ OFF
# THE FITTED RULES. THE CRISIS RESPONSE OF POLICY OVER NO POLICY SPLITS EXACTLY INTO THE REST-POINT ANNOUNCEMENT
# EFFECT AND THE INCREMENTAL RESPONSE FROM EACH ECONOMY'S OWN REST POINT; A MATCHED-STATE COMPARISON IS ALSO KEPT.
#   python3 -u "Claude files/diagnostics/bft_rule/bohn015/k2/compare_k2.py" > <out>/compare_k2.log
import copy, json, pickle, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE.parent))
from global_projection.solver_recursive.state_grid import default_prob, IS, IM, IO
from global_projection.solver_recursive import recursive_experiment as RE
from global_projection.solver_recursive import accuracy as AC
import bohn_compare as BC
from bohn_compare import Econ, RUNS, fitted, clear, two_quarter, shocked, cont_accuracy, quiet

BFT = ROOT.parent / "QUANTITATIVE_MODEL" / "results" / "GLOBAL_2026-11-12_pr38_bft"
OUT = BFT / "bohn015" / "k2"
SHOCKS = (0.015, 0.02)
IRF_KEYS = ("Y", "C", "I", "C_F", "Y_F", "n", "n_F", "mu", "spread", "sov_bp", "Q_bD", "m_cb", "M_cb", "rdep",
            "at_wall")
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:6.0f}s]", *a, flush=True)


def summ(E, v):
    s = BC.summ(E, v)
    B = E.cal["B_gov_D_ss"]
    s.update(m_pct=float(100 * v.m_cb / B), M_pct=float(100 * v.M_cb_new / B), n_F=float(v.n_F), C_F=float(v.C_F),
             Y_F=float(v.Y_F), er_bp=float(4e4 * v.er_D))
    return s


def policy(off, rules_file, er_ref_bp):
    # THE SAME ECONOMY WITH THE ANNOUNCED RULE: ITS CALIBRATION PLUS THE RULE, AND THE SOLVED RULES
    E = copy.copy(off)
    E.cal = dict(off.cal, tpi_on=True, tpi_rule="bft", tpi_gz=False, tpi_bft_er_ref=er_ref_bp / 4e4, tpi_bft_k=2.0)
    E.rules = pickle.load(open(rules_file, "rb"))
    for tag in ("_rest_point_exact", "_rest_point"):
        if hasattr(E.rules, tag):
            delattr(E.rules, tag)
    E.ngh = E.rules.n_gh or RE.N_GH
    return E


def psi_at(E, S_rest, pd, mode):
    v = clear(E, shocked(S_rest, pd))
    if mode == "2q":
        v = two_quarter(E, shocked(S_rest, pd), v.x)
    return summ(E, v)["psi"]


def threshold(E, S_rest, mode, lo=0.004, hi=0.06):
    g = lambda pd: psi_at(E, S_rest, pd, mode)
    glo, ghi = g(lo), g(hi)
    if glo > 0 or ghi < 0:
        return dict(pd=None, psi_lo=glo, psi_hi=ghi)
    return dict(pd=float(brentq(g, lo, hi, xtol=5e-5)), psi_lo=glo, psi_hi=ghi)


def pct(a, b):
    return 100 * (a / b - 1)


KEYS_PCT = ("Q_bD", "n", "n_F", "I", "Y", "C", "C_F", "QK")
KEYS_LVL = ("sov_bp", "ls_bp", "mu", "m_pct", "rdep_bp")


def diff(a, b):
    # a relative to b: percent for quantities and prices, level difference for rates, mu and purchases
    d = {k: pct(a[k], b[k]) for k in KEYS_PCT}
    d.update({k: a[k] - b[k] for k in KEYS_LVL})
    d["CD_over_CF"] = pct(a["C"] / a["C_F"], b["C"] / b["C_F"])
    return d


def main():
    er25 = 15.376861746491961  # the 0.25 run's reference (results/GLOBAL_2026-11-12_pr38_bft/summary)
    er15 = json.load(open(OUT / "solve_info.json"))["er_ref_bp"]
    eco = {}
    for b, rf, pf, er in ((0.25, "rules_12pct_bohn0.250.pkl", BFT / "rules_bft_k2.pkl", er25),
                          (0.15, "rules_12pct_bohn0.150.pkl", OUT / "rules_bft_k2_bohn015.pkl", er15)):
        off = Econ(b, RUNS / rf)
        eco[(b, "off")] = off
        eco[(b, "pol")] = policy(off, pf, er)
    res = {}
    for (b, arm), E in eco.items():
        key = f"{b}_{arm}"
        r = res[key] = {}
        S_rest = RE.stochastic_rest_point(E.rules, E.cal, E.ss, E.sp, verbose=False, exact=True)
        if arm == "pol":
            # the fitted iteration can carry the book below zero, where the rule switches itself off (at Bohn 0.15
            # it settled at -3.4% of the debt): the rest point is the exact fixed point from a zero book at the
            # no-policy rest point, which the exact no-shock path also converges to
            r["rest_fitted_iteration_M"] = float(S_rest[IM])
            S0 = eco[(b, "off")].S_rest.copy()
            S0[IM] = S0[IO] = 0.0
            S_rest, err = RE._exact_rest(E.rules, E.cal, E.ss, E.sp, S0)
            assert err < 1e-9 and S_rest[IM] >= 0.0, (err, S_rest[IM])
            E.rules._rest_point_exact = S_rest.copy()
        E.S_rest = S_rest
        vx = clear(E, S_rest)
        r["rest"] = dict(fitted=summ(E, fitted(E, S_rest)), exact=summ(E, vx),
                         two_quarter=summ(E, two_quarter(E, S_rest, vx.x)), pd=float(default_prob(S_rest[IS])),
                         state=S_rest.tolist())
        log(key, "rest", {m: (round(x["sov_bp"], 2), round(x["mu"], 5), round(x["Y"], 6))
                          for m, x in r["rest"].items() if isinstance(x, dict)})
        r["impact"] = {}
        for pd in SHOCKS + (0.06,):
            S = shocked(S_rest, pd)
            vx = clear(E, S)
            lv = dict(fitted=summ(E, fitted(E, S)), exact=summ(E, vx))
            if pd in SHOCKS:
                lv["two_quarter"] = summ(E, two_quarter(E, S, vx.x))
                r.setdefault("cont_accuracy", {})[str(pd)] = cont_accuracy(E, vx)
            r["impact"][str(pd)] = lv
            log(key, pd, {m: (round(x["ls_bp"], 1), round(x["mu"], 4), round(pct(x["Y"], r["rest"][m]["Y"]), 4))
                          for m, x in lv.items()})
        # the binding threshold (exact one-quarter and joint two-quarter clearing) and a scan around it
        r["threshold"] = {m: threshold(E, S_rest, m) for m in ("exact", "2q")}
        log(key, "threshold", {m: t["pd"] for m, t in r["threshold"].items()})
        r["scan"] = {str(pd): summ(E, clear(E, shocked(S_rest, pd)))
                     for pd in (0.008, 0.010, 0.012, 0.013, 0.014, 0.015, 0.0175, 0.02, 0.025, 0.03, 0.04)}
        # exact dynamic responses
        r["irf"] = {}
        for pd in SHOCKS:
            p, txt = quiet(RE.dynamic_irf, E.rules, E.cal, E.ss, E.sp, pd_shock=pd, T=25, rest_verbose=False,
                           exact=True)
            r["irf"][str(pd)] = {k: np.asarray(p[k]).tolist() for k in IRF_KEYS if k in p}
            r["irf"][str(pd)]["_text"] = txt
        log(key, "irfs done")
    # matched states: the policy rules at the no-policy rest point with the risk shock
    for b in (0.25, 0.15):
        off, pol = eco[(b, "off")], eco[(b, "pol")]
        res[f"{b}_matched"] = {str(pd): summ(pol, clear(pol, shocked(off.S_rest, pd))) for pd in SHOCKS + (0.06,)}
    # off-grid accuracy of the new solution on a simulated ergodic set (fitted rules, as the solver's own report)
    E = eco[(0.15, "pol")]
    states, off_box = AC.simulate(E.rules, E.cal, E.ss, E.sp, E.S_rest, T=1200, burn=200, seed=0)
    err = AC.euler_errors(E.rules, E.cal, E.ss, E.sp, states)
    ok = np.isfinite(err).all(axis=1)
    l10 = np.log10(np.maximum(err[ok], 1e-16))
    pdv = default_prob(states[:, IS])
    res["accuracy_0.15_pol"] = dict(n=int(ok.sum()), l10_mean=float(l10.mean()), l10_max=float(l10.max()),
                                    off_frac=float((off_box > 0).any(axis=1).mean()),
                                    pd_p95=float(np.percentile(pdv, 95)), pd_max=float(pdv.max()))
    log("accuracy", res["accuracy_0.15_pol"])
    # decomposition: policy minus no policy at the shock = announcement (rest points) + incremental (from own rest)
    dec = {}
    for b in (0.25, 0.15):
        for m in ("fitted", "exact", "two_quarter"):
            ro, rp = res[f"{b}_off"]["rest"][m], res[f"{b}_pol"]["rest"][m]
            for pd in SHOCKS:
                so, sp_ = res[f"{b}_off"]["impact"][str(pd)][m], res[f"{b}_pol"]["impact"][str(pd)][m]
                dec[f"{b}_{m}_{pd}"] = dict(no_policy_response=diff(so, ro), policy_response=diff(sp_, rp),
                                            announcement=diff(rp, ro), total=diff(sp_, so),
                                            matched=diff(res[f"{b}_matched"][str(pd)], so) if m == "exact" else None)
    res["decomposition"] = dec
    json.dump(res, open(OUT / "compare_k2.json", "w"), indent=1,
              default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    log("saved", OUT)


if __name__ == "__main__":
    main()
