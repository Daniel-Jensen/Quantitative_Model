# THE SOVEREIGN-SPREAD RESPONSE TO THE HEADLINE RISK SHOCK ON A SAVED SOLVE, FITTED AND EXACTLY CLEARED, WITH THE D-BOND
# FOC ERROR OF THE FITTED RULES ALONG THE PATH (2026-10-23). The response is the shocked path minus the unshocked one,
# both from the rest point (the pipeline's dynamic_irf convention). The FOC error is in price units: the rule's Q_bD
# minus the price its own FOC gives at that state.
# usage: python3 -u irf_spread.py runs/x1_it2.pkl [PD_SHOCK=0.02] [T=25]
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def main():
    # SIMULATE BOTH PATHS, BOTH READS.
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
    pd_shock = float(sys.argv[2]) if len(sys.argv) > 2 else 0.02
    T = int(sys.argv[3]) if len(sys.argv) > 3 else 25
    rules = d_["rules"]
    cal = get_calibration()
    cal["nw_floor_frac"] = d_["nw_floor_frac"]
    cal.update(d_.get("cal_over", {}))
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    sp.update(d_.get("sp_over", {}))
    calibrate_household_anchors(cal, ss, sp)
    g = rules.grid
    yld = lambda Q, db: 4e4 * db * (1.0 - Q) / Q
    sov = lambda o: yld(o["Q_bD"], cal["delta_b_D"]) - yld(o["Q_bF"], cal["delta_b_F"])
    S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)

    def path(shock, read):
        # ONE PATH FROM THE REST POINT: READS, STATES, WALL FLAGS.
        S, out = S0.copy(), []
        if shock:
            S[IS] = RE.s_from_pd(pd_shock)
        for t in range(T):
            o = read(rules, cal, ss, sp, S.copy())
            Sn = RE.advance(o, S, sp)
            Sc = g.clip(Sn)[0]
            out.append((S.copy(), o, bool(np.any(np.abs(Sn - Sc) > 1e-12))))
            S = Sc
        return out

    fit, fit0 = path(True, RE.read_at), path(False, RE.read_at)
    ex, ex0 = path(True, RE.read_exact), path(False, RE.read_exact)
    print(f"rest point: p^d {100 * default_prob(S0[IS]):.4f}%, sov {sov(fit0[0][1]):.1f} bp (fitted) {sov(ex0[0][1]):.1f} "
          f"(exact); max|u| bank stocks {np.max(np.abs(g.to_unit(S0)[0][:8])):.3f}; shock to p^d {100 * pd_shock:.2f}%")
    print(f"{'q':>3s} {'p^d%':>6s} | {'sov fit':>7s} {'sov ex':>7s} | {'dsov fit':>8s} {'dsov ex':>8s} | {'dQ% fit':>7s} "
          f"{'dQ% ex':>7s} | {'FOC err':>8s} | {'dn% fit':>7s} {'mu fit':>7s} {'mu ex':>7s} | {'dloan fit':>9s} "
          f"{'dloan ex':>8s} | wall")
    errs = []
    for t in range(T):
        S, o, wall = fit[t]
        o0, oe, oe0 = fit0[t][1], ex[t][1], ex0[t][1]
        x = np.array([float(rules.eval_read(k, 0, S[None])[0]) for k in SOLVE])
        v = pieces(PM, S, 0, x, rules, cal, ss, sp, False)
        err = v.Q_bD - v.E_Om_payD / v.dmd_D
        errs.append(err)
        loan = lambda r: 4e4 * (r["r_wc_D"] - r["rdep_D"])
        if t < 13 or t % 4 == 0:
            print(f"{t:3d} {100 * default_prob(S[IS]):6.3f} | {sov(o):7.1f} {sov(oe):7.1f} | {sov(o) - sov(o0):+8.1f} "
                  f"{sov(oe) - sov(oe0):+8.1f} | {100 * (o['Q_bD'] / o0['Q_bD'] - 1):+7.2f} "
                  f"{100 * (oe['Q_bD'] / oe0['Q_bD'] - 1):+7.2f} | {err:+8.1e} | {100 * (o['n_D'] / o0['n_D'] - 1):+7.2f} "
                  f"{o['mu_D']:7.4f} {oe['mu_D']:7.4f} | {loan(o) - loan(o0):+9.1f} {loan(oe) - loan(oe0):+8.1f} | "
                  f"{int(wall)}")
    print(f"max |bond FOC error| along the path {np.max(np.abs(errs)):.1e}; quarters at the box wall "
          f"{sum(w for _, _, w in fit)} (shocked), {sum(w for _, _, w in fit0)} (unshocked)")


if __name__ == "__main__":
    main()
