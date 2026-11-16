# WHERE THE TPI IRF HITS THE BOX (2026-10-25): AT THE REST POINT AND ON IMPACT OF THE 2% SHOCK, THE PURCHASE, THE
# EUROSYSTEM BOOK AND THE D BANKS' D BONDS THE LAW OF MOTION WANTS, AGAINST THE BOX BOUNDS, IN % OF THE D DEBT STOCK.
# usage: python3 tpi_wall.py   (reads results/GLOBAL/data)
import json
import os
import pickle
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

if __name__ == "__main__":
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS, IM, IO, IBDD
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    import global_projection.solver_recursive.recursive_experiment as RE
    data = os.path.join(ROOT, "results", "GLOBAL", "data")
    rules = pickle.load(open(os.path.join(data, "rules_tpi.pkl"), "rb"))
    cal = get_calibration(); cal["nw_floor_frac"] = 0.15; cal["tpi_on"] = True; cal["tpi_cap_bp"] = 200.0
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    B = cal["B_gov_D_ss"]; g = rules.grid
    S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)
    corners = g.from_unit(np.array([[s * (k == j) for k in range(g.d)] for j in (IM, IBDD) for s in (-1, 1)]))
    print(f"solve_ok {rules.solve_ok} worst {rules.solve_worst:.1e} cap solved {rules.tpi_cap_solved}")
    print(f"box: M_cb [{100 * corners[0, IM] / B:+.2f}, {100 * corners[1, IM] / B:+.2f}]% of B (axis through the centre); "
          f"b_DD [{100 * corners[2, IBDD] / B:.1f}, {100 * corners[3, IBDD] / B:.1f}]% of B")
    for lab, s in (("rest", None), ("impact, p^d 2%", RE.s_from_pd(0.02))):
        S = S0.copy()
        if s is not None:
            S[IS] = s
        o = RE.read_exact(rules, cal, ss, sp, S)
        Sn = RE.advance(o, S, sp); Sc = g.clip(Sn)[0]
        print(f"{lab:15s} Q_bD {o['Q_bD']:.4f}  purchase m {100 * o['m_cb'] / B:6.2f}% of B  -> book M' {100 * Sn[IM] / B:6.2f}% "
              f"(clipped to {100 * Sc[IM] / B:6.2f}%)  D banks' b_DD' {100 * Sn[IBDD] / B:6.2f}% (clipped {100 * Sc[IBDD] / B:6.2f}%)  "
              f"mu_D {o['mu_D']:.4f}")
