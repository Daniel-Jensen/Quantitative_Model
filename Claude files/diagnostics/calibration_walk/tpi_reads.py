# READS ON THE SAVED RUN (2026-10-25): FOR SHOCKS OF SEVERAL SIZES, THE NO-TPI MARKET SOVEREIGN SPREAD AND, IN THE TPI
# ECONOMY (CAP 200 BP), THE IMPACT PURCHASE AND WHETHER THE NEXT-QUARTER BOOK FITS THE BOX. Reads only: the TPI rules are
# the saved, flagged solve; impact states lie inside the box, so the purchase is a rule read, the path after it is not.
# usage: python3 tpi_reads.py
import os
import pickle
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

if __name__ == "__main__":
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS, IM
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    import global_projection.solver_recursive.recursive_experiment as RE
    data = os.path.join(ROOT, "results", "GLOBAL", "data")
    on = pickle.load(open(os.path.join(data, "rules_tpi.pkl"), "rb"))
    off = pickle.load(open(os.path.join(data, "rules_risk.pkl"), "rb"))
    cal = get_calibration(); cal["nw_floor_frac"] = 0.15
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    c_on, c_off = dict(cal, tpi_on=True, tpi_cap_bp=200.0), dict(cal, tpi_on=False)
    B = cal["B_gov_D_ss"]
    M_top = on.grid.from_unit(np.eye(on.grid.d)[IM][None])[0][IM]
    yld = lambda Q, db: 4e4 * db * (1 - Q) / Q
    sov = lambda o: yld(o["Q_bD"], cal["delta_b_D"]) - yld(o["Q_bF"], cal["delta_b_F"])
    S_on = RE.stochastic_rest_point(on, c_on, ss, sp, verbose=False)
    S_off = RE.stochastic_rest_point(off, c_off, ss, sp, verbose=False)
    print(f"box top of the Eurosystem book: {100 * M_top / B:.1f}% of B (+5% of a half-width by the C1 saturation)")
    print(f"{'shock p^d':>9s} {'no-TPI sov bp':>13s} {'TPI sov bp':>10s} {'purchase %B':>11s} {'fits box':>8s}")
    for pd in (0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02):
        So, Sf = S_on.copy(), S_off.copy()
        So[IS] = Sf[IS] = RE.s_from_pd(pd)
        of, oo = RE.read_exact(off, c_off, ss, sp, Sf), RE.read_exact(on, c_on, ss, sp, So)
        m = oo["m_cb"]
        print(f"{100 * pd:8.2f}% {sov(of):13.1f} {sov(oo):10.1f} {100 * m / B:11.1f} {str(bool(m <= M_top)):>8s}")
