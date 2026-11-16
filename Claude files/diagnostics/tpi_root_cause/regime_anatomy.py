# THE BANKER'S KERNEL AT IMPACT, REGIME BY REGIME, AND THE DEFAULT-REGIME RULES WHERE THE EUROSYSTEM HOLDS BONDS
# (2026-11-06, READ-ONLY).
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import root

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

PEEK = ("alpha_D", "C_D", "N_D", "Q_bD", "Kp_D", "p", "rdep_D", "A_D")


def internals(rules, cal, ss, sproc, S, x):
    # THE PERIOD MAP'S INTERNALS AT (S, x) IN THE NO-DEFAULT REGIME.
    from global_projection.solver_recursive import point_map as PM
    v = PM._read_state_and_regime(S, 0, x, rules, cal, ss)
    PM._production(v, cal)
    PM._government(v, cal, ss)
    PM._balance_sheets(v, rules, cal)
    PM._continuation(v, rules, cal, sproc, rules.n_gh or 5, False)
    PM._discount_kernels(v, rules, cal)
    PM._incentive_constraint(v, cal)
    return v


def solve_regime(rules, cal, ss, sproc, S, d):
    # CLEAR THE PERIOD MAP EXACTLY AT S IN REGIME d, CONTINUATION FROM THE FITTED RULES.
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    x0 = np.array([float(rules.eval(k, d, np.atleast_2d(S))[0]) for k in SOLVE])
    f = lambda x: point_residuals(S, d, x, rules, cal, ss, sproc, n_gh=rules.n_gh or 5)[0]
    sol = root(f, x0, method="hybr", tol=1e-13)
    _, o = point_residuals(S, d, sol.x, rules, cal, ss, sproc, n_gh=rules.n_gh or 5)
    o["_x"] = sol.x
    o["_resid"] = float(np.max(np.abs(sol.fun)))
    return o


def exact_val(o, k):
    # THE EXACT COUNTERPART OF RULE k FROM A SOLVED POINT.
    from global_projection.solver_recursive.decision_rules import SOLVE
    return float(o[k]) if k in ("alpha_D", "C_D") else float(o["_x"][SOLVE.index(k)])


if __name__ == "__main__":
    from post_impact import setup
    from one_off_purchase import read_forced
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS, IM, IO, IBDD
    cal, ss, sproc, r_off, r_on, c_off, c_on, run = setup()
    rests = np.load(HERE / "runs" / "exact_rest_points.npy", allow_pickle=True).item()
    s_sh = RE.s_from_pd(run.RISK_SHOCK_PD)
    B = cal["B_gov_D_ss"]
    cases = {}
    for lab, rules, c, m in (("A no TPI", r_off, c_off, None), ("B TPI", r_on, c_on, None),
                             ("C no TPI + 41% buy", r_off, c_off, 0.4071)):
        S = rests["off" if rules is r_off else "on"].copy(); S[IS] = s_sh
        o = RE.read_exact(rules, c, ss, sproc, S) if m is None else read_forced(rules, c, ss, sproc, S, m)
        cases[lab] = (rules, c, S, o)

    print("THE BANKER'S KERNEL AT IMPACT: E_Om = sum_j w_j E[Lam_j (f + (1-f) alpha'_j)]")
    for lab, (rules, c, S, o) in cases.items():
        v = internals(rules, c, ss, sproc, S, o["_x"])
        print(f"\n [{lab}]  E_Om {v.E_Om_D:.4f}  pd {100*sum(w.sum() for w in v.wgt[1:]):.2f}%  "
              f"b_DD' {v.b_D_D_new:.4f}  M' {v.M_cb_new/B:.4f}B  O' {v.Op_cb:.4f}  Q_bD {v.Q_bD:.4f}  Q_K {v.Q_D:.5f}")
        for j in range(v.nreg):
            w = v.wgt[j]; W = w.sum()
            lam = v.Om_D[j] / (c["f_D"] + (1 - c["f_D"]) * v.R[j]["alpha_D"])
            print(f"   regime {j}: weight {W:.4f}  E[Lam] {np.dot(w, lam)/W:.4f}  E[alpha'] {np.dot(w, v.R[j]['alpha_D'])/W:.4f}"
                  f"  E[Om] {np.dot(w, v.Om_D[j])/W:.4f}  contribution {np.dot(w, v.Om_D[j]):.4f}"
                  f"  E[C'] {np.dot(w, v.R[j]['C_D'])/W:.4f}  E[Q_bD'] {np.dot(w, v.R[j]['Q_bD'])/W:.4f}"
                  f"  E[R_K'] {4e4*np.dot(w, v.RKD[j])/W:8.1f}bp")
            print(f"      alpha' by GH node: {np.round(v.R[j]['alpha_D'], 4)}   Lam: {np.round(lam, 4)}")

    # the default regime's rules against an exact default-regime solve at the next states (the GH nodes' centre)
    print("\nRULE vs EXACT AT THE NEXT STATE (centre GH node), BOTH REGIMES: % gaps (rdep in bp/yr)")
    for lab, (rules, c, S, o) in cases.items():
        Sn = RE.advance(o, S, sproc)
        Sn = rules.grid.clip(Sn)[0]
        u = rules.grid.to_unit(np.atleast_2d(Sn))[0]
        print(f" [{lab}]  next state in box units: s {u[IS]:+.3f}  b_DD+M {u[IBDD]:+.3f}  M {u[IM]:+.3f}  O {u[IO]:+.3f}")
        for d in (0, 1):
            ox = solve_regime(rules, c, ss, sproc, Sn, d)
            gaps = []
            for k in PEEK:
                fit = float(rules.eval(k, d, np.atleast_2d(Sn))[0])
                ex = exact_val(ox, k)
                gaps.append(4e4 * (fit - ex) if k == "rdep_D" else 100 * (fit / ex - 1))
            print(f"   regime {d} (resid {ox['_resid']:.0e}): " + "  ".join(f"{k} {g:+.3f}" for k, g in zip(PEEK, gaps))
                  + f"   exact alpha {ox['alpha_D']:.4f} mu {ox['mu_D']:.4f} n_D {ox['n_D']:.4f}")

    # the M axis of the rules: every other state at the box centre, s at the centre and at the top
    print("\nTHE RULES ALONG THE M AXIS (other states at the box centre)")
    g = r_on.grid
    for rules, nm in ((r_off, "no TPI"), (r_on, "TPI")):
        for s_u in (0.0, 1.0):
            for d in (0, 1):
                line = []
                for uM in (-1.0, -0.5, 0.0, 0.47, 0.87, 1.0):
                    uu = np.zeros(g.d); uu[IS] = s_u; uu[IM] = uM
                    X = g.from_unit(uu)
                    line.append(" ".join(f"{float(rules.eval(k, d, X)[0]):7.4f}" for k in ("alpha_D", "C_D", "Q_bD")))
                print(f"  {nm:6s} s_u {s_u:+.0f} regime {d}: (alpha, C_D, Q_bD) at u_M -1,-.5,0,.47,.87,1:  " + " | ".join(line))
