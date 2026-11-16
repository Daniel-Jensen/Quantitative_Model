# OFF-GRID ACCURACY ON AN EXACTLY CLEARED SIMULATION (2026-11-13). THE SOLVER'S OWN REPORT SIMULATES THE FITTED RULES;
# UNDER THE ANNOUNCED RULE THE FITTED PURCHASE DIPS BELOW ZERO AT ONCE, THE BOOK TURNS NEGATIVE AND THE RULE SWITCHES
# ITSELF OFF, SO 89% OF THAT SAMPLE IS AN ECONOMY WITHOUT THE BACKSTOP. HERE EVERY QUARTER IS CLEARED EXACTLY (THE
# BOOK STAYS >= 0), s DRAWS AS IN accuracy.simulate (SAME SEED, NO-DEFAULT HISTORY), AND THE FITTED RULES' EULER
# ERRORS ARE EVALUATED ON THE VISITED STATES. NO SOLVE.
import json, pickle, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE.parent / "k2"))
from results_io import load
from global_projection.solver_recursive.state_grid import IM, IS, default_prob
from global_projection.solver_recursive import accuracy as AC
from global_projection.solver_recursive import recursive_experiment as RE
import headline as H

T, BURN = 1200, 200


def main():
    d = json.load(open(H.OUT / "headline.json"))
    cal = load(H.DATA / "calibration"); ss = load(H.DATA / "steady_state"); cal["tpi_on"] = False
    info = json.load(open(H.OUT / "solve_k2_info.json"))
    cal_k = dict(cal, tpi_on=True, tpi_rule="bft", tpi_gz=False, tpi_bft_er_ref=info["er_ref_bp"] / 4e4,
                 tpi_bft_k=2.0)
    B = cal["B_gov_D_ss"]
    out = {}
    for arm, c, rf in (("off", cal, H.DATA / "rules_risk.pkl"), ("pol", cal_k, H.OUT / "rules_bft_k2.pkl")):
        E = H.Eco(c, ss, pickle.load(open(rf, "rb")))
        rng = np.random.default_rng(0)
        S = np.array(d[arm]["rest"]["state"])
        states, mu, ls, book, buy = [], [], [], [], []
        for t in range(T + BURN):
            o = RE.read_exact(E.rules, E.cal, E.ss, E.sp, S.copy())
            if t >= BURN:
                states.append(S.copy()); mu.append(o["mu_D"]); ls.append(4e4 * (o["r_wc_D"] - o["rdep_D"]))
                book.append(o["M_cb_new"] / B); buy.append(o["m_cb"] / B)
            Sn = RE.advance(o, S, E.sp)
            Sn[IS] = (1 - E.sp["rho_s"]) * E.sp["s_star"] + E.sp["rho_s"] * S[IS] + E.sp["sigma_s"] * rng.standard_normal()
            S = E.rules.grid.clip(Sn)[0]
        states = np.array(states)
        err = AC.euler_errors(E.rules, E.cal, E.ss, E.sp, states)
        ok = np.isfinite(err).all(axis=1)
        l10 = np.log10(np.maximum(err[ok], 1e-16))
        pdv = default_prob(states[:, IS])
        hi = pdv > 0.01
        out[arm] = dict(n=int(ok.sum()), l10_mean=float(l10.mean()), l10_max=float(l10.max()),
                        l10_max_by_eq=dict(zip(AC._EQN, l10.max(axis=0).tolist())),
                        l10_mean_pd_above_1pct=float(l10[hi[ok]].mean()) if hi.any() else None,
                        pd_p95=float(np.percentile(pdv, 95)), pd_max=float(pdv.max()),
                        share_binding=float(np.mean(np.array(mu) > 1e-8)), ls_max_bp=float(max(ls)),
                        book_min_pct=float(100 * min(book)), book_max_pct=float(100 * max(book)),
                        buy_min_pct=float(100 * min(buy)), outside_box=float(np.mean(
                            E.rules.grid.outside(states).max(axis=1) > 0)))
        print(arm, out[arm], flush=True)
    json.dump(out, open(H.OUT / "accuracy_exact.json", "w"), indent=1)


if __name__ == "__main__":
    main()
