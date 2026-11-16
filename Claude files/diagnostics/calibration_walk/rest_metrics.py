# THE FLOOR SWEEP'S REST-POINT NUMBERS FOR ANY SAVED SOLVE (probe_solve / floor_sweep pickles), 2026-10-23: THE TARGETED
# MOMENTS (FITTED AND EXACT), mu AT REST, THE SOVEREIGN SPREAD, A DEFAULT AT REST, AND THE REST POINT FROM OTHER STARTS.
# usage: python3 -u rest_metrics.py runs/pd002.pkl
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def main():
    # READ THE REST POINT.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES, IS, IM, IO
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors, ss_state
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    import global_projection.solver_recursive.recursive_experiment as RE
    import global_projection.solver_recursive.calibrate_stochastic as CS
    t0 = time.perf_counter()
    d_ = pickle.load(open(os.path.join(HERE, sys.argv[1]), "rb"))
    rules = d_["rules"]
    rules.solve_ok, rules.solve_worst = d_.get("ok"), d_.get("worst", float("nan"))
    cal = get_calibration()
    cal["nw_floor_frac"] = d_["nw_floor_frac"]
    cal.update(d_.get("cal_over", {}))
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    sp.update(d_.get("sp_over", {}))
    calibrate_household_anchors(cal, ss, sp)
    g, n_ss = rules.grid, ss["ss_bank_D"]["n_ss"]
    S0, o = CS.rest_point_read(rules, cal, ss, sp)
    oe = RE.read_exact(rules, cal, ss, sp, S0.copy())
    mf, me = CS.rest_point_moments(o), CS.rest_point_moments(oe)
    yield_ = lambda Q, db: 4e4 * db * (1.0 - Q) / Q
    sov = lambda r: yield_(r["Q_bD"], cal["delta_b_D"]) - yield_(r["Q_bF"], cal["delta_b_F"])
    # a default at the rest point: the default-regime rules there, and the state it leads to
    x1 = np.array([float(rules.eval_read(k, 1, S0[None])[0]) for k in SOLVE])
    _, o1 = point_residuals(S0, 1, x1, rules, cal, ss, sp, n_gh=rules.n_gh or RE.N_GH, no_default=False)
    o1["_x"] = x1
    S1 = RE.advance(o1, S0.copy(), sp)
    u0, u1 = g.to_unit(S0)[0], g.to_unit(S1)[0]
    far = int(np.argmax(np.abs(u1[:8])))
    # the rest point from other starts (no shocks): its spread across starts in unit coordinates
    rng, ends = np.random.default_rng(1), []
    for _ in range(4):
        u = rng.uniform(-0.8, 0.8, g.d)
        u[[IM, IO]] = 0.0
        S = g.from_unit(u[None])[0]
        S[IS] = ss_state(ss, cal, sp)[IS]
        for _ in range(4000):
            Sn = RE.advance(RE.read_at(rules, cal, ss, sp, S), S, sp, g)
            if np.max(np.abs(Sn - S)) < 1e-11:
                break
            S = Sn
        ends.append(g.to_unit(S)[0])
    spread_starts = float(np.max(np.ptp(np.array(ends + [u0]), axis=0)))
    row = dict(ok=bool(rules.solve_ok), worst=float(rules.solve_worst),
               lev_D=me["lev_D"], lev_F=me["lev_F"],
               spread_D_fit=4e4 * mf["spread_D"], spread_D_exact=4e4 * me["spread_D"],
               spread_F_exact=4e4 * me["spread_F"], exp_D=100 * me["exp_D"], exp_F=100 * me["exp_F"],
               mu_D_fit=float(o["mu_D"]), mu_D_exact=float(oe["mu_D"]),
               sov_spread_fit=sov(o), sov_spread_exact=sov(oe),
               n_D_rest=float(o["n_D"] / n_ss), alpha_D_rest=float(o["alpha_D"]),
               n_D_default=float(o1["n_D"] / n_ss), alpha_D_default=float(o1["alpha_D"]),
               rest_u_max=float(np.max(np.abs(u0[:8]))),
               post_default_u_max=float(np.abs(u1[far])), post_default_far=STATE_NAMES[far],
               rest_spread_across_starts=spread_starts, read_resid=float(o["_resid"]),
               seconds=time.perf_counter() - t0)
    print("  ".join(f"{k} {v:.4g}" if isinstance(v, float) else f"{k} {v}" for k, v in row.items()), flush=True)


if __name__ == "__main__":
    main()
