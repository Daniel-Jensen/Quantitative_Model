# THE BI-FOERSTER-TRAUM PURCHASE RULE AS AN ANNOUNCED COMMITMENT ON THE PR #38 SOLUTION (2026-11-12).
# Starts from the saved no-TPI rules of results/GLOBAL_2026-11-12_pr38 (main repo), walks the rule's coefficient k
# through the pre-specified set K_SET (stopping at the first that fails), and for each solved k reports accuracy, the
# Eurosystem book against the box, the rest point (announcement) and the crisis responses (intervention), all cleared
# exactly. Nothing in the private economy, the calibration, the reserve treatment or the solver is changed.
#   python3 "Claude files/diagnostics/bft_rule/run_bft.py"    (from the PR #38 worktree root)
import copy
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from results_io import load, save
from global_projection.solver_recursive.state_grid import s_process_params, default_prob, IS, IM, STATE_NAMES
from global_projection.solver_recursive.decision_rules import SOLVE
from global_projection.solver_recursive.point_map import point_residuals
from global_projection.solver_recursive.accuracy import accuracy_report
from global_projection.solver_recursive import recursive_experiment as RE

MAIN = ROOT.parent / "QUANTITATIVE_MODEL"
DATA = MAIN / "results" / "GLOBAL_2026-11-12_pr38" / "data"
OUT = MAIN / "results" / "GLOBAL_2026-11-12_pr38_bft"
# the pre-specified coefficients, in BFT units (phi_cb = k * Q_ss * B_ss; BFT's own value is 16)
K_SET = (2.0, 4.0, 8.0, 16.0)
# the PR #38 headline shock and the deeper crisis shock of 2026-11-06, p^d per quarter
SHOCKS = (0.015, 0.06)
T_IRF, T_ACC = 25, 1200
IRF_KEYS = ("Y", "C", "I", "C_F", "Y_F", "n", "n_F", "mu", "loan_D", "spread", "sov_bp", "Q_bD", "dQ_bD",
            "m_cb", "M_cb", "Pi_cb", "rdep", "at_wall")


def bp(x):
    return 4e4 * x


def rest_row(rules, cal, ss, sproc):
    # THE STOCHASTIC REST POINT, CLEARED EXACTLY.
    S = RE.stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    o = RE.read_exact(rules, cal, ss, sproc, S.copy())
    B = cal["B_gov_D_ss"]
    row = dict(Q_bD=o["Q_bD"], sov_bp=RE._sov_bp(o, cal), er_bp=bp(o["er_D"]), mu_D=o["mu_D"], n_D=o["n_D"],
               n_F=o["n_F"], loan_bp=bp(o["r_wc_D"] - o["rdep_D"]), rdep_bp=bp(o["rdep_D"]), Y_D=o["Y_D"],
               C_D=o["C_D"], I_D=o["I_D"], Y_F=o["Y_F"], C_F=o["C_F"], alpha_D=o["alpha_D"],
               M_pct=100 * o["M_cb_new"] / B, m_pct=100 * o["m_cb"] / B, resid=o["_resid"])
    return row, S


def node_checks(rules, cal, ss, sproc):
    # AT EVERY NO-DEFAULT NODE: THE RULE, THE BOOK, AND WHETHER NEXT QUARTER'S STATE LEAVES THE BOX.
    B, rows = cal["B_gov_D_ss"], []
    for i in range(rules.grid.n):
        S = rules.grid.points[i]
        x = np.array([rules.vals[k][0][i] for k in SOLVE])
        res, o = point_residuals(S, 0, x, rules, cal, ss, sproc, n_gh=rules.n_gh or 5)
        o["_x"] = x
        out = rules.grid.outside(np.atleast_2d(RE.advance(o, S, sproc)))[0]
        rows.append(dict(pd=float(default_prob(S[IS])), M_lag=S[IM] / B, M_new=o["M_cb_new"] / B,
                         m=o["m_cb"] / B, er_bp=bp(o["er_D"]), mu=o["mu_D"], b_DD=o["b_D_D_new"] / B,
                         resid=float(np.max(np.abs(res))), out_M=out[IM],
                         out_other=float(np.max(np.delete(out, [IS, IM])))))
    r = {k: np.array([row[k] for row in rows]) for k in rows[0]}
    return dict(n=len(rows), max_resid=float(r["resid"].max()), n_buy=int((r["m"] > 1e-4).sum()),
                max_M_new=float(r["M_new"].max()), band_M=float(rules.grid.hi[IM] / B),
                n_out_M=int((r["out_M"] > 0).sum()), max_out_M=float(r["out_M"].max()),
                n_out_other=int((r["out_other"] > 0).sum()), min_b_DD=float(r["b_DD"].min()), rows=r)


def path_checks(rules, cal, ss, sproc, irf, S_rest, pd):
    # ALONG THE EXACT CRISIS PATH: FITTED-vs-EXACT GAP, THE F GOODS MARKET (WALRAS) AND THE BOOK'S IDENTITIES.
    S = S_rest.copy()
    S[IS] = RE.s_from_pd(pd)
    o_x = RE.read_exact(rules, cal, ss, sproc, S.copy())
    o_f = RE.read_at(rules, cal, ss, sproc, S.copy())
    B = cal["B_gov_D_ss"]
    clear = o_x["Bp_D"] - o_x["b_D_D_new"] - o_x["b_D_F_new"] - o_x["M_cb_new"]
    return dict(Y_exact=o_x["Y_D"], Y_fit=o_f["Y_D"], mu_exact=o_x["mu_D"], mu_fit=o_f["mu_D"],
                m_exact=100 * o_x["m_cb"] / B, m_fit=100 * o_f["m_cb"] / B, fit_resid=o_f["_resid"],
                goods_F=o_x["goods_F"], clearing=clear, claim=o_x["Z_cb"] - o_x["Q_bD"] * o_x["M_cb_new"],
                walled=int(np.sum(irf["at_wall"])))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cal0 = load(DATA / "calibration")
    ss0 = load(DATA / "steady_state")
    cal0["tpi_on"] = False
    sproc = s_process_params(cal0)
    with open(DATA / "rules_risk.pkl", "rb") as fh:
        rules_off = pickle.load(fh)
    print(f"no-TPI rules: {rules_off.grid.n} points, solve_ok {rules_off.solve_ok}, smoothing "
          f"{getattr(rules_off, 'smoothing', None)}")
    chk = node_checks(rules_off, cal0, ss0, sproc)
    print(f"  saved solution reproduces at the nodes: max|F| {chk['max_resid']:.1e}")

    rest_off, S_off = rest_row(rules_off, cal0, ss0, sproc)
    er_ref = rest_off["er_bp"] / 4e4
    print(f"  no-TPI rest point: ER {rest_off['er_bp']:.2f} bp/yr (the rule's reference), sovereign spread "
          f"{rest_off['sov_bp']:.1f} bp, mu_D {rest_off['mu_D']:.2e}")
    irf_off = {pd: RE.dynamic_irf(rules_off, cal0, ss0, sproc, pd_shock=pd, T=T_IRF, rest_verbose=False, exact=True)
               for pd in SHOCKS}
    summary = dict(er_ref_bp=rest_off["er_bp"], K_SET=K_SET, SHOCKS=SHOCKS, rest_off=rest_off,
                   irf_off={str(pd): {k: irf_off[pd][k] for k in IRF_KEYS} for pd in SHOCKS}, runs={})

    cal_b = dict(cal0, tpi_on=True, tpi_rule="bft", tpi_gz=False, tpi_bft_er_ref=er_ref)
    prev, k_prev = rules_off, 0.0
    for k in K_SET:
        t0 = time.perf_counter()
        print(f"\n===== BFT rule, k = {k:g} (from k = {k_prev:g}) =====", flush=True)
        cal_k = dict(cal_b, tpi_bft_k=k, tpi_bft_from=k_prev)
        ss_k = copy.deepcopy(ss0)
        r = RE.solve_recursive(cal_k, ss_k, sproc, s_refine=0, base=prev)
        ok = bool(r.solve_ok) and r.tpi_bft_solved == k
        run = dict(k=k, ok=ok, solved=r.tpi_bft_solved, worst=float(r.solve_worst), minutes=0.0,
                   smoothing=getattr(r, "smoothing", None))
        with open(OUT / f"rules_bft_k{k:g}.pkl", "wb") as fh:
            pickle.dump(r, fh)
        run["minutes"] = (time.perf_counter() - t0) / 60
        print(f"  k = {k:g}: ok {ok}, last k solved {r.tpi_bft_solved}, max|F| {r.solve_worst:.1e}, "
              f"{run['minutes']:.0f} min", flush=True)
        summary["runs"][str(k)] = run
        if not ok:
            save(OUT / "summary", summary)
            print(f"  STOP: k = {k:g} did not solve; no larger coefficient is attempted")
            break
        nodes = node_checks(r, cal_k, ss_k, sproc)
        run["nodes"] = {kk: v for kk, v in nodes.items() if kk != "rows"}
        run["node_rows"] = nodes["rows"]
        print(f"  nodes: max|F| {nodes['max_resid']:.1e}; buying at {nodes['n_buy']}/{nodes['n']}; largest next "
              f"book {100 * nodes['max_M_new']:.1f}% of the debt against a box band of {100 * nodes['band_M']:.1f}%; "
              f"next book outside the box at {nodes['n_out_M']} nodes (worst {100 * nodes['max_out_M']:.1f}% of "
              f"width); other states outside at {nodes['n_out_other']}; min D-bank book {nodes['min_b_DD']:.3f}")
        run["rest_on"], S_on = rest_row(r, cal_k, ss_k, sproc)
        run["irf_on"], run["paths"] = {}, {}
        for pd in SHOCKS:
            irf = RE.dynamic_irf(r, cal_k, ss_k, sproc, pd_shock=pd, T=T_IRF, rest_verbose=False, exact=True)
            run["irf_on"][str(pd)] = {kk: irf[kk] for kk in IRF_KEYS}
            run["paths"][str(pd)] = path_checks(r, cal_k, ss_k, sproc, irf, S_on, pd)
            print(f"  shock {pd}: {run['paths'][str(pd)]}", flush=True)
        acc = accuracy_report(r, cal_k, ss_k, sproc, S_on, T=T_ACC, label=f"BFT k={k:g}")
        if acc is not None:
            run["accuracy"] = dict(l10_mean=acc["l10_mean"], l10_max=acc["l10_max"], off_frac=acc["off_frac"])
        save(OUT / "summary", summary)
        prev, k_prev = r, k
    save(OUT / "summary", summary)
    print("\ndone")


if __name__ == "__main__":
    main()
