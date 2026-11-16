# CRISIS-SIZED SHOCKS ON THE APPENDIX-B 12% BOX (2026-11-06, READ-ONLY; data/global_warm_start.pkl, NO TPI):
# which p^d shock gives a recession with the acceptance signs, how much a one-off purchase must buy to switch the D
# constraint off, and what that purchase does after impact with the no-TPI rules continuing.
import pickle
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

SHOCKS = (0.04, 0.05, 0.06, 0.08, 0.09)
SHARES = (0.05, 0.10, 0.15, 0.20, 0.25)
T = 13
SHOW = (0, 1, 2, 4, 8, 12)


def setup12():
    # THE APPENDIX-B ECONOMY ON ITS 12% BOX, AS SOLVED (FIX A).
    import run
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    cal = get_calibration()
    cal.update(nw_floor_frac=run.NW_FLOOR, jacobian=run.JACOBIAN, pd_box_top=0.12, tpi_on=False)
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    rules = pickle.load(open(ROOT / "data" / "global_warm_start.pkl", "rb"))
    return cal, ss, sproc, rules


def row(rules, cal, ss, S, o):
    # ONE QUARTER'S OBJECTS (LEVELS).
    from one_off_purchase import capital_price
    yD = cal["delta_b_D"] * (1 - o["Q_bD"]) / o["Q_bD"]
    yF = cal["delta_b_F"] * (1 - o["Q_bF"]) / o["Q_bF"]
    return dict(sov=4e4 * (yD - yF), Q=o["Q_bD"], n=o["n_D"], nF=o["n_F"], Y=o["Y_D"], C=o["C_D"], I=o["I_D"],
                N=o["_x"][0], mu=o["mu_D"], lend=4e4 * (o["r_wc_D"] - o["rdep_D"]), rdep=4e4 * o["rdep_D"],
                QK=capital_price(rules, cal, ss, S, o), Pi=100 * o["Pi_cb"] / ss["ss_firm_D"]["Y_ss"],
                M=100 * o["M_cb_new"] / cal["B_gov_D_ss"])


def walk(rules, cal, ss, sproc, S0, s_shock=None, m0=None):
    # T QUARTERS FROM S0, EXACT EACH QUARTER; OPTIONAL SHOCK AND A ONE-OFF PURCHASE (SHARE OF THE DEBT) AT t = 0.
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    from one_off_purchase import read_forced
    S = S0.copy()
    if s_shock is not None:
        S[IS] = s_shock
    out = []
    for t in range(T):
        o = (read_forced(rules, cal, ss, sproc, S, m0) if (t == 0 and m0) else RE.read_exact(rules, cal, ss, sproc, S))
        out.append(row(rules, cal, ss, S, o))
        S = rules.grid.clip(RE.advance(o, S, sproc))[0]
    return {k: np.array([r[k] for r in out]) for k in out[0]}


LEVELS = ("sov", "mu", "lend", "rdep", "Pi", "M")


def dev(a, ref):
    # DEVIATION FROM THE UNSHOCKED PATH: % FOR QUANTITIES, DIFFERENCES FOR RATES, SPREADS, mu AND SHARES.
    return {k: (a[k] - ref[k]) if k in LEVELS else 100 * (a[k] / ref[k] - 1) for k in a}


if __name__ == "__main__":
    from exact_rest import exact_rest
    from one_off_purchase import read_forced
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    cal, ss, sproc, rules = setup12()
    S_fit = RE.stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    S0, err = exact_rest(rules, cal, ss, sproc, S_fit)
    np.save(HERE / "runs" / "crisis_rest_12pct.npy", S0)
    o0 = RE.read_exact(rules, cal, ss, sproc, S0.copy())
    r0 = row(rules, cal, ss, S0, o0)
    print(f"REST (exact fixed point, {err:.0e}): sov {r0['sov']:.1f} bp, mu_D {r0['mu']:.4f}, lending spread "
          f"{r0['lend']:.1f} bp, n_D {r0['n']:.4f}")
    ref = walk(rules, cal, ss, sproc, S0)

    print("\n1. NO-TPI IRFs (exact): % vs the unshocked path; spreads/rates bp/yr and mu as differences; quarters", SHOW)
    for pd in SHOCKS:
        d = dev(walk(rules, cal, ss, sproc, S0, RE.s_from_pd(pd)), ref)
        signs = dict(Q_bD=d["Q"][0] < 0, n_D=d["n"][0] < 0, n_F=d["nF"][0] < 0, Y_D=d["Y"][0] < 0,
                     C_D=d["C"][0] < 0, lending=d["lend"][0] > 0)
        print(f"\n p^d {100*pd:.0f}%   acceptance signs: " + ", ".join(f"{k} {'ok' if v else 'WRONG'}" for k, v in signs.items()))
        for k in ("sov", "Q", "n", "nF", "Y", "C", "I", "N", "mu", "lend", "rdep", "QK"):
            print(f"   {k:5s} " + " ".join(f"{d[k][t]:+9.3f}" for t in SHOW))

    print("\n2. ONE-OFF PURCHASE AT IMPACT, RULES HELD FIXED (% of the D debt): mu_D, lending spread, sov spread, D output")
    m_star = {}
    for pd in SHOCKS:
        S = S0.copy(); S[IS] = RE.s_from_pd(pd)
        o_a = RE.read_exact(rules, cal, ss, sproc, S)
        cells = []
        for m in SHARES:
            o_b = read_forced(rules, cal, ss, sproc, S, m)
            cells.append(f"{100*m:.0f}%: mu {o_b['mu_D']:.3f}, lend {4e4*(o_b['r_wc_D']-o_b['rdep_D']):5.0f}, sov "
                         f"{4e4*(cal['delta_b_D']*(1-o_b['Q_bD'])/o_b['Q_bD'] - cal['delta_b_F']*(1-o_b['Q_bF'])/o_b['Q_bF']):5.0f}, "
                         f"Y {100*(o_b['Y_D']/o_a['Y_D']-1):+.3f}%")
        print(f" p^d {100*pd:.0f}% (no buy: mu {o_a['mu_D']:.3f}, lend {4e4*(o_a['r_wc_D']-o_a['rdep_D']):.0f} bp)")
        for c in cells:
            print("    " + c)
        # the purchase that just switches the constraint off (mu crosses 0): rules held fixed
        if o_a["mu_D"] > 1e-6:
            g = lambda m: read_forced(rules, cal, ss, sproc, S, m)["mu_D"] - 1e-6
            hi = next((m for m in SHARES if g(m) < 0), None)
            m_star[pd] = brentq(g, 0.0, hi, xtol=1e-5) if hi else None
            print(f"    purchase that switches the constraint off: "
                  + (f"{100*m_star[pd]:.1f}% of the debt" if m_star[pd] else "more than 25%"))

    print("\n3. THE PATH AFTER A ONE-OFF PURCHASE THAT SWITCHES THE CONSTRAINT OFF (no-TPI rules continue, key_D 0.071)")
    for pd, m in m_star.items():
        if not m:
            continue
        a = dev(walk(rules, cal, ss, sproc, S0, RE.s_from_pd(pd)), ref)
        b = dev(walk(rules, cal, ss, sproc, S0, RE.s_from_pd(pd), m), ref)
        print(f"\n p^d {100*pd:.0f}%, buy {100*m:.1f}%:   no purchase | one-off purchase, quarters {SHOW}")
        for k in ("sov", "n", "Y", "C", "I", "mu", "lend", "Pi", "M"):
            print(f"   {k:5s} " + " ".join(f"{a[k][t]:+8.3f}" for t in SHOW) + "  |  "
                  + " ".join(f"{b[k][t]:+8.3f}" for t in SHOW))
