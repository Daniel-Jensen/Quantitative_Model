# THE D BOND'S FUNDAMENTAL VALUE AT THE MODEL'S OWN DISCOUNT RATES, AGAINST THE FITTED PRICE RULE (2026-10-23).
# On a fine grid in s (the other states at the rest point) the script reads, from the solved rules, the bank's required
# return req(s) = 1 + rdep + lambda_bD mu / E[Om], the Om-weighted no-default weights on each s' quadrature node and the
# default leg. Holding those fixed it solves the perpetuity recursion Q(s) = (sum_k w_k(s) (delta + (1 - delta) Q(s'_k))
# + default leg(s)) / req(s) on the fine grid (a contraction, linear interpolation in s, s' clamped to the box). If the
# rule's price is far above this fixed point, the rule is not what the model's own discount rates price.
# It then solves the same recursion by Chebyshev collocation in s alone at m = 3, 5, 9, 17 nodes (the inputs read off the
# fine grid), to show what a degree m-1 rule in s does to the price with nothing else approximated.
# usage: python3 diag_fundamental.py runs/pd002.pkl
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def main():
    # READ THE DISCOUNT RATES, SOLVE THE RECURSION, COMPARE.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS, default_prob
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.decision_rules import SOLVE
    import global_projection.solver_recursive.point_map as PM
    import global_projection.solver_recursive.recursive_experiment as RE
    from diag_slack import pieces
    d_ = pickle.load(open(os.path.join(HERE, sys.argv[1]), "rb"))
    rules = d_["rules"]
    cal = get_calibration()
    cal["nw_floor_frac"] = d_["nw_floor_frac"]
    cal.update(d_.get("cal_over", {}))
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    sp.update(d_.get("sp_over", {}))
    calibrate_household_anchors(cal, ss, sp)
    g, db = rules.grid, cal["delta_b_D"]
    S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)
    sg = np.linspace(g.lo[IS], g.hi[IS], 81)
    W, SN, DL, REQ, QR, MU, RD = [], [], [], [], [], [], []
    for s in sg:
        S = S0.copy()
        S[IS] = s
        x = np.array([float(rules.eval_read(k, 0, S[None])[0]) for k in SOLVE])
        v = pieces(PM, S, 0, x, rules, cal, ss, sp, False)
        eps, _ = PM.gh_nodes(len(v.wgt[0]))
        Sn, _ = PM._next_states(v, sp, eps)
        nd = [j for j in range(v.nreg) if not v.reg[j]][0]
        W.append(v.wgt[nd] * np.asarray(v.Om_D[nd]) / v.E_Om_D)
        SN.append(np.clip(Sn[:, IS], g.lo[IS], g.hi[IS]))
        DL.append(sum(float(np.dot(v.wgt[j], np.asarray(v.Om_D[j]) * cal["recovery_rate_D"]
                                   * (db + (1 - db) * np.asarray(v.R[j]["Q_bD"]))))
                      for j in range(v.nreg) if v.reg[j]) / v.E_Om_D)
        REQ.append(1 + v.rdep_D + cal["lambda_bD_D"] * v.mu_D / v.E_Om_D)
        QR.append(v.Q_bD); MU.append(v.mu_D); RD.append(v.rdep_D)
    W, SN, DL, REQ, QR = map(np.array, (W, SN, DL, REQ, QR))
    Q = np.full(sg.size, 0.9)
    for it in range(5000):
        Qn = (np.sum(W * (db + (1 - db) * np.interp(SN, sg, Q)), axis=1) + DL) / REQ
        if np.max(np.abs(Qn - Q)) < 1e-13:
            break
        Q = Qn
    # the same recursion with the rule's own Q' (one step): equals the rule where the FOC holds
    Q1 = (np.sum(W * (db + (1 - db) * np.interp(SN, sg, QR)), axis=1) + DL) / REQ
    u_rest = float(g.to_unit(S0)[0][IS])
    print(f"fixed point in {it} iterations; rest point s {S0[IS]:.3f} (u {u_rest:+.2f})")
    print(f"{'u_s':>6s} {'p^d %':>8s} {'req bp':>7s} {'rdep bp':>7s} {'mu':>7s} {'Q rule':>7s} {'Q 1-step':>8s} {'Q fundam':>8s} {'gap':>7s}")
    for i in range(0, sg.size, 5):
        u = float(g.to_unit(np.r_[S0[:IS], sg[i], S0[IS + 1:]][None])[0][IS])
        print(f"{u:+6.2f} {100 * default_prob(sg[i]):8.4f} {4e4 * (REQ[i] - 1):7.1f} {4e4 * RD[i]:7.1f} {MU[i]:7.4f} "
              f"{QR[i]:7.4f} {Q1[i]:8.4f} {Q[i]:8.4f} {QR[i] - Q[i]:+7.4f}")
    # the recursion collocated in s alone, the inputs interpolated from the fine grid
    to_u = lambda s_: 2 * (s_ - g.lo[IS]) / (g.hi[IS] - g.lo[IS]) - 1
    T = lambda u, m: np.cos(np.arange(m) * np.arccos(np.clip(u, -1, 1))[..., None])
    ug, ur = to_u(sg), to_u(S0[IS])
    print(f"{'nodes':>5s} {'Q at rest':>9s} {'max gap to fundamental':>23s} {'Q at u=-1':>9s} {'monotone':>8s}")
    for m in (3, 5, 9, 17):
        un = -np.cos(np.pi * np.arange(m) / (m - 1))
        Wn = np.array([[np.interp(u, ug, W[:, k]) for k in range(W.shape[1])] for u in un])
        SNn = np.array([[np.interp(u, ug, SN[:, k]) for k in range(W.shape[1])] for u in un])
        DLn, REQn = np.interp(un, ug, DL), np.interp(un, ug, REQ)
        # Q(u_i) = sum_k W_ik (db + (1-db) T(u'_ik) c) / REQ_i + DL_i / REQ_i, linear in c
        A = T(un, m) - (1 - db) * np.einsum("ik,ikn->in", Wn, T(to_u(SNn), m)) / REQn[:, None]
        c = np.linalg.solve(A, (db * Wn.sum(1) + DLn) / REQn)
        Qm = T(ug, m) @ c
        print(f"{m:5d} {float(T(np.array(ur), m) @ c):9.4f} {np.max(np.abs(Qm - Q)):23.4f} {Qm[0]:9.4f} "
              f"{str(bool(np.all(np.diff(Qm) <= 1e-12))):>8s}")
    print(f"at rest: rule {np.interp(S0[IS], sg, QR):.4f}  fundamental {np.interp(S0[IS], sg, Q):.4f}  "
          f"(yield {4e4 * db * (1 - np.interp(S0[IS], sg, Q)) / np.interp(S0[IS], sg, Q):.1f} bp vs rule "
          f"{4e4 * db * (1 - np.interp(S0[IS], sg, QR)) / np.interp(S0[IS], sg, QR):.1f})")


if __name__ == "__main__":
    main()
