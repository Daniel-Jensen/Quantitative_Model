# THE 0.02% ANOMALY: THE D BOND PRICE RULE ALONG s AT THE REST POINT'S BANK STATE, THE NODE VALUES IT INTERPOLATES, AND
# THE BOND'S RETURN AT REST FROM THE RULE ITSELF (A STATIONARY POINT: E[Q'] = Q) VS THE DEPOSIT RATE.
# usage: python3 anomaly_q.py runs/pd002.pkl
import os, pickle, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def main():
    # PROBE.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.point_map import default_prob
    import global_projection.solver_recursive.recursive_experiment as RE
    d_ = pickle.load(open(os.path.join(HERE, sys.argv[1]), "rb"))
    rules, g = d_["rules"], d_["rules"].grid
    cal = get_calibration(); cal["nw_floor_frac"] = d_["nw_floor_frac"]; cal.update(d_.get("cal_over", {}))
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); sp.update(d_.get("sp_over", {}))
    calibrate_household_anchors(cal, ss, sp)
    S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)
    db = cal["delta_b_D"]
    print(f"Q_bD_ss {ss['Q_bD_ss']:.4f}; rest s {S0[IS]:+.3f} (p^d {100 * default_prob(S0[IS]):.4f}%/qtr), u_s {g.to_unit(S0)[0][IS]:+.3f}")
    for u in (-1.0, -0.6, -0.41, -0.2, 0.0, 0.5, 1.0):
        S = S0.copy(); S[IS] = g.from_unit(np.where(np.arange(g.d) == IS, u, g.to_unit(S0)[0])[None])[0][IS]
        q0 = float(rules.eval_read("Q_bD", 0, S[None])[0]); r0 = float(rules.eval_read("rdep_D", 0, S[None])[0])
        y = db * (1 - q0) / q0
        print(f"  u_s {u:+.2f} s {S[IS]:+.3f} p^d {100 * default_prob(S[IS]):7.4f}%/qtr: Q_bD rule {q0:.4f}, yield {400 * y:5.2f}%/yr, "
              f"rdep {400 * r0:5.2f}%/yr")
    nodes = np.where(np.all(np.abs(np.delete(g.to_unit(g.points), IS, axis=1)) < 1e-12, axis=1))[0]
    for i in nodes:
        print(f"  centre node s {g.points[i, IS]:+.3f}: Q_bD node value d0 {rules.vals['Q_bD'][0][i]:.4f}")


if __name__ == "__main__":
    main()
