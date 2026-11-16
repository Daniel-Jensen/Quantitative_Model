# THE CANONICAL HEADLINE RESULTS (BOHN 0.15; SHOCKS 2% AND 1.5%), NO POLICY AGAINST THE ANNOUNCED k = 2 RULE, 2026-11-13.
# NO SOLVE: THE TWO SAVED CANONICAL SOLUTIONS ARE CLEARED EXACTLY (ONE QUARTER), JOINTLY OVER TWO QUARTERS, AND READ
# OFF THE FITTED RULES. ALSO: REPRODUCTION OF THE DIAGNOSTIC SOLUTIONS, THE ACCOUNTING IDENTITIES AT EVERY REPORTED
# STATE, THE SIGN AND NORMALISATION OF THE EUROSYSTEM BOOK, THE BINDING THRESHOLDS AND OFF-GRID ACCURACY.
#   python3 -u "Claude files/diagnostics/bft_rule/bohn015/canonical/headline.py"   (from the worktree root)
import json, pickle, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE.parent / "k2"))
from results_io import load
from global_projection.solver_recursive.state_grid import s_process_params, default_prob, IS, IM, IO
from global_projection.solver_recursive.decision_rules import SOLVE
from global_projection.solver_recursive import recursive_experiment as RE
from global_projection.solver_recursive import accuracy as AC
from bohn_compare import fitted, clear, two_quarter, shocked, cont_accuracy, quiet, RUNS
from compare_k2 import summ, threshold, diff, pct, IRF_KEYS

DATA = ROOT / "results" / "GLOBAL" / "data"
OUT = ROOT.parent / "QUANTITATIVE_MODEL" / "results" / "GLOBAL_2026-11-13_bohn015"
DIAG = ROOT.parent / "QUANTITATIVE_MODEL" / "results" / "GLOBAL_2026-11-12_pr38_bft" / "bohn015" / "k2"
SHOCKS = (0.02, 0.015)
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:6.0f}s]", *a, flush=True)


class Eco:
    def __init__(self, cal, ss, rules):
        self.cal, self.ss, self.rules = cal, ss, rules
        self.sp = s_process_params(cal)
        self.ngh = rules.n_gh or RE.N_GH
        for tag in ("_rest_point_exact", "_rest_point"):
            if hasattr(rules, tag):
                delattr(rules, tag)


def identities(E, v):
    # THE ACCOUNTING IDENTITIES AT ONE CLEARED STATE (EACH SHOULD BE ZERO TO ROUNDING)
    c = E.cal
    f, om, dK = c["f_D"], c["omega_ent_D"], c["delta_D"]
    payF = c["delta_b_F"] + (1 - c["delta_b_F"]) * v.Q_bF
    n_sum = ((1 - f) * (v.mpk_D + (1 - dK) * v.Q_D) * v.K_D + (1 - f) * v.Xi_D * v.b_D_D_lag
             + (1 - f) * v.p * payF * v.b_F_D_lag + (1 - f) * v.O_lag - (1 - f) * v.P_D + om * v.assets_D + v.recap_D)
    psi = 1 - v.E_Om_D * (1 + v.rdep_D) * v.n_D / v.lev_D
    return dict(
        net_worth=float(n_sum - v.n_D),  # n = (1-f)(returns on capital, bonds, claim - deposits) + start-up + recap
        bond_market=float(v.Bp_D - v.b_D_D_new - v.b_D_F_new - v.M_cb_new),  # D bonds: issued = held
        multiplier=float(v.mu_D - max(psi, 0.0)),  # mu = max(psi, 0) (the fitted kink is not used: exact)
        lending_spread=float((v.r_wc_D - v.rdep_D) - c["lambda_K_D"] * v.mu_D / v.E_Om_D),
        cb_claim=float(v.Z_cb - v.Q_bD * v.M_cb_new),  # the Eurosystem's claim = market value of its book
        goods_F=float(v.goods_F),  # Walras: the F goods market clears when every other market does
        book_runoff=float(v.M_cb_new - (1 - c["delta_b_D"]) * v.surv * v.M_lag - v.m_cb),  # book = runoff + purchase
    )


def purchases(E, v):
    B, Qss = E.cal["B_gov_D_ss"], E.ss["Q_bD_ss"]
    return dict(face_pct_ss_debt=float(100 * v.m_cb / B), mkt_pct_ss_debt=float(100 * v.Q_bD * v.m_cb / (Qss * B)),
                pct_current_debt=float(100 * v.m_cb / v.Bp_D), pct_private_stock=float(
                    100 * v.m_cb / ((1 - E.cal["delta_b_D"]) * v.surv * (v.b_D_D_lag + v.b_D_F_lag))),
                book_pct_ss_debt=float(100 * v.M_cb_new / B), M_lag_pct=float(100 * v.M_lag / B))


def main():
    cal = load(DATA / "calibration"); ss = load(DATA / "steady_state")
    cal["tpi_on"] = False
    info = json.load(open(OUT / "solve_k2_info.json"))
    off = Eco(cal, ss, pickle.load(open(DATA / "rules_risk.pkl", "rb")))
    cal_k = dict(cal, tpi_on=True, tpi_rule="bft", tpi_gz=False, tpi_bft_er_ref=info["er_ref_bp"] / 4e4, tpi_bft_k=2.0)
    pol = Eco(cal_k, ss, pickle.load(open(OUT / "rules_bft_k2.pkl", "rb")))
    eco = dict(off=off, pol=pol)
    res = dict(settings=dict(bohn=cal["bohn_annual"], B_gov_D_ss=cal["B_gov_D_ss"], Q_bD_ss=ss["Q_bD_ss"],
                             er_ref_bp=info["er_ref_bp"], solve_k2=info))
    # 1. reproduction: the canonical solutions against the diagnostic ones of the same day
    diag = json.load(open(DIAG / "compare_k2.json"))
    rep = {}
    for arm, rf in (("off", RUNS / "rules_12pct_bohn0.150.pkl"), ("pol", DIAG / "rules_bft_k2_bohn015.pkl")):
        old = pickle.load(open(rf, "rb"))
        new = eco[arm].rules
        same_grid = bool(np.allclose(old.grid.points, new.grid.points))
        gap = {k: float(np.max(np.abs(np.asarray(new.vals[k][0]) - np.asarray(old.vals[k][0]))
                            / np.maximum(np.abs(np.asarray(old.vals[k][0])), 1e-3))) for k in SOLVE}
        rep[arm] = dict(same_grid=same_grid, max_rel_gap=max(gap.values()), worst=max(gap, key=gap.get))
    res["reproduction"] = rep
    log("reproduction", rep)
    # 2. rest points (the policy economy's from a zero book: the fitted iteration can push the book below zero,
    # where the rule switches itself off)
    for arm, E in eco.items():
        r = res[arm] = {}
        S = RE.stochastic_rest_point(E.rules, E.cal, E.ss, E.sp, verbose=False, exact=True)
        r["rest_book_from_fitted_iteration"] = float(S[IM])
        if arm == "pol":
            S0 = eco["off"].S_rest.copy(); S0[IM] = S0[IO] = 0.0
            S, err = RE._exact_rest(E.rules, E.cal, E.ss, E.sp, S0)
            assert err < 1e-9 and S[IM] >= 0.0, (err, S[IM])
            E.rules._rest_point_exact = S.copy()
        E.S_rest = S
        vx = clear(E, S)
        r["rest"] = dict(fitted=summ(E, fitted(E, S)), exact=summ(E, vx), two_quarter=summ(E, two_quarter(E, S, vx.x)),
                         pd=float(default_prob(S[IS])), state=S.tolist(), identities=identities(E, vx),
                         purchases=purchases(E, vx))
        log(arm, "rest", {m: (round(x["sov_bp"], 2), round(x["mu"], 5)) for m, x in r["rest"].items()
                          if m in ("fitted", "exact", "two_quarter")}, "book", r["rest"]["purchases"]["book_pct_ss_debt"])
        r["impact"], r["cont_accuracy"] = {}, {}
        for pd in SHOCKS + (0.06,):
            Sx = shocked(S, pd)
            vx = clear(E, Sx)
            lv = dict(fitted=summ(E, fitted(E, Sx)), exact=summ(E, vx), identities=identities(E, vx),
                      purchases=purchases(E, vx))
            if pd in SHOCKS:
                v2 = two_quarter(E, Sx, vx.x)
                lv["two_quarter"] = summ(E, v2)
                lv["identities_2q"] = identities(E, v2)
                lv["purchases_2q"] = purchases(E, v2)
                r["cont_accuracy"][str(pd)] = cont_accuracy(E, vx)
            r["impact"][str(pd)] = lv
            log(arm, pd, {m: (round(lv[m]["ls_bp"], 1), round(lv[m]["mu"], 4))
                          for m in ("fitted", "exact", "two_quarter") if m in lv},
                "max identity gap", max(abs(x) for x in lv["identities"].values()))
        r["threshold"] = {m: threshold(E, S, m) for m in ("exact", "2q")}
        log(arm, "threshold", {m: t["pd"] for m, t in r["threshold"].items()})
        r["scan"] = {str(pd): summ(E, clear(E, shocked(S, pd)))
                     for pd in (0.010, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025, 0.03, 0.04, 0.06)}
        r["irf"] = {}
        for pd in SHOCKS:
            p, txt = quiet(RE.dynamic_irf, E.rules, E.cal, E.ss, E.sp, pd_shock=pd, T=25, rest_verbose=False,
                           exact=True)
            r["irf"][str(pd)] = {k: np.asarray(p[k]).tolist() for k in IRF_KEYS if k in p}
            r["irf"][str(pd)]["_text"] = txt
        # off-grid accuracy on a simulated ergodic set, and the Eurosystem book along it
        states, off_box = AC.simulate(E.rules, E.cal, E.ss, E.sp, S, T=1200, burn=200, seed=0)
        err = AC.euler_errors(E.rules, E.cal, E.ss, E.sp, states)
        ok = np.isfinite(err).all(axis=1)
        l10 = np.log10(np.maximum(err[ok], 1e-16))
        pdv = default_prob(states[:, IS])
        r["accuracy"] = dict(n=int(ok.sum()), l10_mean=float(l10.mean()), l10_max=float(l10.max()),
                             l10_mean_by_eq=l10.mean(axis=0).tolist(), l10_max_by_eq=l10.max(axis=0).tolist(),
                             off_frac=float((np.asarray(off_box) > 0).any(axis=-1).mean()) if np.ndim(off_box) else
                             float(np.mean(off_box)),
                             pd_p95=float(np.percentile(pdv, 95)), pd_max=float(pdv.max()),
                             book_min_pct=float(100 * states[:, IM].min() / cal["B_gov_D_ss"]),
                             book_max_pct=float(100 * states[:, IM].max() / cal["B_gov_D_ss"]))
        log(arm, "accuracy", r["accuracy"]["l10_mean"], r["accuracy"]["l10_max"], "book",
            r["accuracy"]["book_min_pct"], r["accuracy"]["book_max_pct"])
    # 3. reproduction of the diagnostic numbers (exact one-quarter clearing)
    res["reproduction"]["numbers"] = {
        arm: {pd: dict(ls_new=res[arm]["impact"][pd]["exact"]["ls_bp"], ls_diag=diag[f"0.15_{arm}"]["impact"][pd]["exact"]["ls_bp"],
                       Y_new=res[arm]["impact"][pd]["exact"]["Y"], Y_diag=diag[f"0.15_{arm}"]["impact"][pd]["exact"]["Y"])
              for pd in ("0.015", "0.02")} for arm in ("off", "pol")}
    # 4. policy over no policy at the shock = announcement (rest points) + incremental crisis response
    dec = {}
    for m in ("fitted", "exact", "two_quarter"):
        ro, rp = res["off"]["rest"][m], res["pol"]["rest"][m]
        for pd in SHOCKS:
            so, sp_ = res["off"]["impact"][str(pd)][m], res["pol"]["impact"][str(pd)][m]
            a, nr, prr, tot = diff(rp, ro), diff(so, ro), diff(sp_, rp), diff(sp_, so)
            dec[f"{m}_{pd}"] = dict(no_policy_response=nr, policy_response=prr, announcement=a, total=tot,
                                    incremental={k: prr[k] - nr[k] for k in prr},
                                    # for rates, mu and purchases the split is exact; for percent changes to first order
                                    split_gap={k: tot[k] - a[k] - (prr[k] - nr[k]) for k in tot})
    res["decomposition"] = dec
    json.dump(res, open(OUT / "headline.json", "w"), indent=1,
              default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    log("saved", OUT / "headline.json")


if __name__ == "__main__":
    main()
