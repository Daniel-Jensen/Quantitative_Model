# WHY THE TPI ECONOMY DOES WORSE AFTER IMPACT (2026-11-06): REST POINTS, SAME-STATE READS AND RULE ACCURACY ALONG BOTH
# IRF PATHS, ON THE SAVED RULES (READ-ONLY).
import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

KEYS = ("Q_bD", "n_D", "E_Om_D", "alpha_D", "rdep_D", "C_D", "N_D", "Y_D", "I_D", "p", "A_D", "Tax_D",
        "inc_D", "Div_D", "W_D", "dep_D", "Pi_cb", "m_cb", "M_cb_new", "b_D_D_new", "Bp_D", "C_F", "n_F")


def setup():
    # THE CALIBRATION AND BOTH SOLVED ECONOMIES, AS run.py's COMPARE SEES THEM.
    import run
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    data = ROOT / "results" / "GLOBAL" / "data"
    cal = get_calibration()
    cal.update(nw_floor_frac=run.NW_FLOOR, jacobian=run.JACOBIAN, pd_box_top=run.PD_BOX_TOP)
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    cal["tpi_cap_bp"] = float(run.TPI_CAP)
    r_off = pickle.load(open(data / "rules_risk.pkl", "rb"))
    r_on = pickle.load(open(data / "rules_tpi.pkl", "rb"))
    c_on = dict(cal, tpi_on=True, tpi_eps=getattr(r_on, "tpi_eps_solved", cal["tpi_eps"]))
    c_off = dict(cal, tpi_on=False)
    return cal, ss, sproc, r_off, r_on, c_off, c_on, run


def irf_states(rules, cal, ss, sproc, pd_shock, T):
    # THE EXACT IRF PATH'S STATES AND READS (AS dynamic_irf(exact=True) WALKS IT).
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    S = RE.stochastic_rest_point(rules, cal, ss, sproc, verbose=False).copy()
    S[IS] = RE.s_from_pd(pd_shock)
    Ss, os_ = [], []
    for _ in range(T):
        o = RE.read_exact(rules, cal, ss, sproc, S)
        Ss.append(S.copy()); os_.append(o)
        S = rules.grid.clip(RE.advance(o, S, sproc))[0]
    return Ss, os_


def rule_gap(rules, cal, ss, sproc, S, o_exact):
    # FITTED RULE MINUS EXACT SOLVE AT S, FOR THE RULES THE CONTINUATION READS.
    from global_projection.solver_recursive.decision_rules import SOLVE
    Sm = np.atleast_2d(S)
    out = {}
    for k, ex in (("Q_bD", o_exact["Q_bD"]), ("C_D", o_exact["C_D"]), ("alpha_D", o_exact["alpha_D"]),
                  ("N_D", o_exact["_x"][0]), ("p", o_exact["_x"][6]), ("Kp_D", o_exact["_x"][2]),
                  ("rdep_D", o_exact["_x"][4]), ("C_F", o_exact["C_F"])):
        fit = float(rules.eval(k, 0, Sm)[0])
        out[k] = 100.0 * (fit / ex - 1.0) if k != "rdep_D" else 4e4 * (fit - ex)
    return out


if __name__ == "__main__":
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import STATE_NAMES, IS, IM
    cal, ss, sproc, r_off, r_on, c_off, c_on, run = setup()
    np.set_printoptions(linewidth=200, precision=5, suppress=True)

    # 1. the two rest points, state by state
    S_off = RE.stochastic_rest_point(r_off, c_off, ss, sproc, verbose=False)
    S_on = RE.stochastic_rest_point(r_on, c_on, ss, sproc, verbose=False)
    print("REST STATES (on vs off)")
    for i, nm in enumerate(STATE_NAMES):
        d = (100 * (S_on[i] / S_off[i] - 1)) if abs(S_off[i]) > 1e-9 else float("nan")
        print(f"  {nm:6s} off {S_off[i]: .6e}  on {S_on[i]: .6e}  diff% {d:+.4f}")

    # 2. both economies read exactly at the SAME state (the no-TPI rest point): only the continuation rules differ
    print("\nSAME STATE (the no-TPI rest point): TPI rules vs no-TPI rules, exact reads")
    for lab, S in (("S_off", S_off), ("S_on", S_on)):
        a = RE.read_exact(r_off, c_off, ss, sproc, S.copy())
        b = RE.read_exact(r_on, c_on, ss, sproc, S.copy())
        print(f"  at {lab}: resid off {a['_resid']:.1e} on {b['_resid']:.1e}")
        for k in KEYS:
            va, vb = float(a[k]), float(b[k])
            d = 100 * (vb / va - 1) if abs(va) > 1e-9 else vb - va
            print(f"    {k:10s} off {va: .6f}  on {vb: .6f}  diff {d:+.4f}{'%' if abs(va) > 1e-9 else ''}")

    # 3. rule accuracy along both exact IRF paths: what quarter t's expectations read at S_{t+1} vs the exact solve
    T = 13
    for lab, rules, c in (("NO TPI", r_off, c_off), ("TPI", r_on, c_on)):
        Ss, os_ = irf_states(rules, c, ss, sproc, run.RISK_SHOCK_PD, T)
        print(f"\n{lab} PATH: fitted-rule minus exact (%; rdep in bp/yr), and the state's TPI book")
        print("   t   M%B    b_DD     s      " + "  ".join(f"{k:>8s}" for k in
                                                       ("Q_bD", "C_D", "alpha_D", "N_D", "p", "Kp_D", "rdep_D", "C_F")))
        for t in range(T):
            g = rule_gap(rules, c, ss, sproc, Ss[t], os_[t])
            print(f"  {t:2d} {100*Ss[t][IM]/cal['B_gov_D_ss']:6.2f} {Ss[t][4]:7.4f} {Ss[t][IS]:6.3f} "
                  + "  ".join(f"{g[k]:+8.4f}" for k in ("Q_bD", "C_D", "alpha_D", "N_D", "p", "Kp_D", "rdep_D", "C_F"))
                  + f"   resid {os_[t]['_resid']:.0e}")
