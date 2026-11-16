# A ONE-OFF PURCHASE IN THE NO-TPI ECONOMY (2026-11-06, READ-ONLY): THE SAME 41% BOUGHT AT IMPACT, THEN THE NO-TPI RULES
# RUN (NO CAP, SO NO CAP KINK IN THE RULES). SEPARATES "THE EUROSYSTEM HOLDS THE BONDS" FROM THE TPI RULES' OWN FIT.
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import root

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

T = 25
SHOW = (0, 1, 2, 4, 8, 12, 16, 24)


def read_forced(rules, cal, ss, sproc, S, m):
    # THE PERIOD MAP CLEARED WITH THE PURCHASE FIXED AT m (TPI OFF: THE LAST UNKNOWN IS THE PURCHASE ITSELF).
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    ngh = rules.n_gh or RE.N_GH
    x0 = np.array([float(rules.eval(k, 0, np.atleast_2d(S))[0]) for k in SOLVE])[:-1]
    f = lambda z: point_residuals(S, 0, np.append(z, m), rules, cal, ss, sproc, n_gh=ngh)[0][:-1]
    sol = root(f, x0, method="hybr", tol=1e-13)
    _, o = point_residuals(S, 0, np.append(sol.x, m), rules, cal, ss, sproc, n_gh=ngh)
    o["_x"] = np.append(sol.x, m)
    o["_resid"] = float(np.max(np.abs(sol.fun)))
    return o


def capital_price(rules, cal, ss, S, o):
    # THE D CAPITAL PRICE Q_K AT THE SOLVED POINT.
    from global_projection.solver_recursive import point_map as PM
    v = PM._read_state_and_regime(S, 0, o["_x"], rules, cal, ss)
    PM._production(v, cal)
    return float(v.Q_D)


def record(rules, cal, ss, S, o):
    # ONE QUARTER'S OBJECTS.
    from global_projection.solver_recursive.state_grid import IM
    yD = cal["delta_b_D"] * (1 - o["Q_bD"]) / o["Q_bD"]
    yF = cal["delta_b_F"] * (1 - o["Q_bF"]) / o["Q_bF"]
    Y = ss["ss_firm_D"]["Y_ss"]
    return dict(sov_bp=4e4 * (yD - yF), Q_bD=o["Q_bD"], Q_K=capital_price(rules, cal, ss, S, o), n_D=o["n_D"],
                I_D=o["I_D"], Y_D=o["Y_D"], C_D=o["C_D"], N_D=o["_x"][0], rdep_D=4e4 * o["rdep_D"],
                rdep_F=4e4 * o["rdep_F"], kap_nfa=4e4 * cal.get("kappa_nfa", 0.0) * o["nfa_dep_D"] / Y,
                nfa=o["nfa_dep_D"] / Y, B=o["Bp_D"], Tax_D=o["Tax_D"], Div_D=o["Div_D"], inc_D=o["inc_D"],
                Pi_cb=o["Pi_cb"] / Y, M=S[IM] / cal["B_gov_D_ss"], E_Om=o["E_Om_D"], C_F=o["C_F"], p=o["p"],
                A_D=o["A_D"], resid=o["_resid"])


def path(rules, cal, ss, sproc, S0, s_shock=None, m0=None):
    # T QUARTERS FROM S0: OPTIONAL RISK SHOCK AND OPTIONAL FORCED PURCHASE AT t = 0.
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    S = S0.copy()
    if s_shock is not None:
        S[IS] = s_shock
    rows = []
    for t in range(T):
        o = (read_forced(rules, cal, ss, sproc, S, m0) if (t == 0 and m0 is not None)
             else RE.read_exact(rules, cal, ss, sproc, S))
        rows.append(record(rules, cal, ss, S, o))
        S = rules.grid.clip(RE.advance(o, S, sproc))[0]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


LEVEL = ("sov_bp", "rdep_D", "rdep_F", "kap_nfa", "nfa", "M", "Pi_cb", "resid")


def dev(a, ref):
    # DEVIATIONS FROM THE REFERENCE PATH: % FOR QUANTITIES, LEVEL DIFFERENCES FOR RATES AND SHARES.
    return {k: (a[k] - ref[k]) if k in LEVEL else 100 * (a[k] / ref[k] - 1) for k in a}


if __name__ == "__main__":
    from exact_rest import exact_rest
    from post_impact import setup
    from global_projection.solver_recursive import recursive_experiment as RE
    cal, ss, sproc, r_off, r_on, c_off, c_on, run = setup()
    cache = HERE / "runs" / "exact_rest_points.npy"
    if cache.exists():
        rests = np.load(cache, allow_pickle=True).item()
    else:
        rests = {lab: exact_rest(r, c, ss, sproc, RE.stochastic_rest_point(r, c, ss, sproc, verbose=False))[0]
                 for lab, r, c in (("off", r_off, c_off), ("on", r_on, c_on))}
        np.save(cache, rests, allow_pickle=True)
    s_sh = RE.s_from_pd(run.RISK_SHOCK_PD)
    B = cal["B_gov_D_ss"]
    ref_off = path(r_off, c_off, ss, sproc, rests["off"])
    ref_on = path(r_on, c_on, ss, sproc, rests["on"])
    runs = {
        "A no TPI, risk shock": dev(path(r_off, c_off, ss, sproc, rests["off"], s_sh), ref_off),
        "B TPI, risk shock": dev(path(r_on, c_on, ss, sproc, rests["on"], s_sh), ref_on),
        # x_cb is the purchase as a share of the SS debt (point_map._tpi_purchase_and_price)
        "C no TPI, risk shock + one-off 41% buy": dev(path(r_off, c_off, ss, sproc, rests["off"], s_sh, 0.4071), ref_off),
        "D no TPI, one-off 41% buy, no shock": dev(path(r_off, c_off, ss, sproc, rests["off"], None, 0.4071), ref_off),
    }
    np.save(HERE / "runs" / "one_off_purchase.npy", dict(runs=runs, ref_off=ref_off, ref_on=ref_on), allow_pickle=True)
    print("quarters:", SHOW, "  (% dev from own unshocked path; rates bp/yr and shares as level differences)")
    for k in ("sov_bp", "Q_bD", "Q_K", "n_D", "I_D", "Y_D", "C_D", "N_D", "rdep_D", "rdep_F", "kap_nfa", "nfa", "B",
              "Div_D", "inc_D", "Pi_cb", "M", "E_Om", "C_F", "p", "resid"):
        print(f"\n {k}")
        for lab, d in runs.items():
            print(f"   {lab:42s} " + " ".join(f"{d[k][t]:+9.4f}" for t in SHOW))
    print("\nREFERENCE LEVELS at t = 0 / 24:")
    for k in ("sov_bp", "rdep_D", "nfa", "B", "n_D"):
        print(f"   {k:8s} off {ref_off[k][0]:.5f} / {ref_off[k][-1]:.5f}   on {ref_on[k][0]:.5f} / {ref_on[k][-1]:.5f}")
