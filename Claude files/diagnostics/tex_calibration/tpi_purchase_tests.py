# THE PURCHASE TESTS OF 2026-10-30, IN ONE PLACE (READ-ONLY; THE SOLVED RULES FROM results/GLOBAL/data/rules_risk.pkl):
#   1. a one-off purchase of 5% / 15% of the D debt stock (rules held fixed) at p^d 2/5/12%: bond price, lending spread,
#      output, investment, for the current design (the TPI claim counts against the D bank's IC at the bond's lambda) and
#      a Gertler-Karadi variant (central-bank reserves not divertable; a diagnostic patch, no model edit);
#   2. Bocola's capital (lending) premia, E[r^k'] - r = lambda_K*mu/E[Om] - Cov(Om, r^k')/E[Om], at rest and p^d 2/5/12%,
#      with D output against the rest point (the damage the purchases of test 1 repair).
#   python3 -u tpi_purchase_tests.py
import inspect, os, pickle, sys, textwrap
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
import numpy as np
from scipy.optimize import root


def setup():
    # THE PIPELINE'S ECONOMY AND THE SOLVED RULES.
    import run
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    cal = get_calibration(); cal.update(nw_floor_frac=run.NW_FLOOR, jacobian=run.JACOBIAN)
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    rules = pickle.load(open(ROOT / "results/GLOBAL/data/rules_risk.pkl", "rb"))
    return cal, ss, sp, rules


def purchases(cal, ss, sp, r):
    # TEST 1: THE ONE-OFF PURCHASE, BOTH DESIGNS.
    import global_projection.solver_recursive.point_map as pm
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.state_grid import IS
    c, B = dict(cal, tpi_on=False), cal["B_gov_D_ss"]
    S0 = RE.stochastic_rest_point(r, c, ss, sp, verbose=False)
    orig = pm._incentive_constraint
    ns = {}
    exec(textwrap.dedent(inspect.getsource(orig)).replace("lbDD * v.Z_cb", "0.0 * v.Z_cb"), pm.__dict__, ns)

    def at(S, m, x0):
        # THE PERIOD MAP CLEARED AT S WITH A PURCHASE m (TPI OFF: THE LAST UNKNOWN IS THE PURCHASE).
        f = lambda z: pm.point_residuals(S, 0, np.append(z, m), r, c, ss, sp, n_gh=5)[0][:-1]
        sol = root(f, x0, method="hybr", tol=1e-12)
        return sol.x, pm.point_residuals(S, 0, np.append(sol.x, m), r, c, ss, sp, n_gh=5)[1]
    spread = lambda o: 4e4 * (o["r_wc_D"] - o["rdep_D"])
    for name, fn in (("reserves count like bonds (current)", orig), ("reserves not divertable (Gertler-Karadi)", ns["_incentive_constraint"])):
        pm._incentive_constraint = fn
        print(name)
        for pd in (0.02, 0.05, 0.12):
            S = S0.copy(); S[IS] = RE.s_from_pd(pd)
            x, o = at(S, 0.0, np.array([float(r.eval_read(k, 0, S[None])[0]) for k in SOLVE])[:-1])
            rows = []
            for share in (0.05, 0.15):
                _, o1 = at(S, share * B, x)
                rows.append(f"buy {100 * share:.0f}%: price {100 * (o1['Q_bD'] / o['Q_bD'] - 1):+.2f}%, lending spread "
                            f"{spread(o1) - spread(o):+.0f} bp, output {100 * (o1['Y_D'] / o['Y_D'] - 1):+.3f}%, investment "
                            f"{100 * (o1['I_D'] / o['I_D'] - 1):+.2f}%, mu {o['mu_D']:.3f}->{o1['mu_D']:.3f}")
            print(f"   p^d {100 * pd:.0f}%: " + " | ".join(rows))
    pm._incentive_constraint = orig


def capital_premia(cal, ss, sp, r):
    # TEST 2: BOCOLA'S CAPITAL PREMIA AT EXACT READS.
    import global_projection.solver_recursive.point_map as pm
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    S0 = RE.stochastic_rest_point(r, cal, ss, sp, verbose=False)
    print("Bocola's capital premia, E[r^k'] - r (bp/yr); D output against the rest point")
    y_rest = None
    for pd in (None, 0.02, 0.05, 0.12):
        S = S0.copy()
        if pd:
            S[IS] = RE.s_from_pd(pd)
        o = RE.read_exact(r, cal, ss, sp, S)
        y_rest = y_rest or o["Y_D"]
        v = pm._read_state_and_regime(S, 0, o["_x"], r, cal, ss)
        pm._production(v, cal); pm._government(v, cal, ss); pm._balance_sheets(v, r, cal)
        pm._continuation(v, r, cal, sp, r.n_gh or 5, False); pm._discount_kernels(v, r, cal); pm._incentive_constraint(v, cal)
        E = lambda arrs: sum(float(np.dot(v.wgt[j], arrs[j])) for j in range(v.nreg))
        erk, eom = E(v.RKD), v.E_Om_D
        risk = 4e4 * (erk - E([v.Om_D[j] * v.RKD[j] for j in range(v.nreg)]) / eom)
        liq = 4e4 * cal["lambda_K_D"] * v.mu_D / eom
        print(f"   p^d {'rest' if pd is None else f'{100 * pd:.0f}%':>5s}: {4e4 * (erk - v.rdep_D):7.1f} = risk {risk:6.1f} + "
              f"liquidity {liq:6.1f}  (mu_D {v.mu_D:.4f}), Y_D {100 * (o['Y_D'] / y_rest - 1):+.3f}%")


if __name__ == "__main__":
    cal, ss, sp, r = setup()
    purchases(cal, ss, sp, r)
    capital_premia(cal, ss, sp, r)
