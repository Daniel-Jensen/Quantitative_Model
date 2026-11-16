# WHY IS THE D BANK SLACK AT THE REST POINT? ACCOUNTING ON A SAVED SOLVE, NO NEW SOLVES (2026-10-22).
# The multiplier is mu = kink(z), z = 1 - x, x = E[Omega] (1 + rdep) n / lev (lev = divertable assets). The script
# reads every piece at the stochastic rest point (fitted and exactly cleared) and at the deterministic SS (pi = 0,
# SS-constant rules: the model nests the SS there), and splits log x - log x_ss into the four pieces. It also prices
# the D bond's no-default excess return, the income the bank keeps on a path where default is priced but never comes.
# usage: python3 diag_slack.py runs/floor_0.15.pkl
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def pieces(PM, S, d, x, rules, cal, ss, sp, no_default):
    # THE PERIOD MAP UP TO THE INCENTIVE CONSTRAINT AND THE BOND FOCS, RETURNING ITS STATE OBJECT.
    v = PM._read_state_and_regime(S, d, x, rules, cal, ss)
    PM._production(v, cal)
    PM._government(v, cal, ss)
    PM._balance_sheets(v, rules, cal)
    PM._continuation(v, rules, cal, sp, rules.n_gh or 5, no_default)
    PM._discount_kernels(v, rules, cal)
    PM._incentive_constraint(v, cal)
    PM._bond_demands(v, cal)
    return v


def main():
    # DECOMPOSE.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors, ss_state, ss_x
    from global_projection.solver_recursive.decision_rules import RuleSet, SOLVE
    import global_projection.solver_recursive.point_map as PM
    import global_projection.solver_recursive.recursive_experiment as RE
    d_ = pickle.load(open(os.path.join(HERE, sys.argv[1]), "rb"))
    rules = d_["rules"]
    cal = get_calibration()
    cal["nw_floor_frac"] = d_["nw_floor_frac"]
    for k, v in d_.get("cal_over", {}).items():
        cal[k] = v
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    sp.update(d_.get("sp_over", {}))
    calibrate_household_anchors(cal, ss, sp)
    S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)
    xf = np.array([float(rules.eval_read(k, 0, S0[None])[0]) for k in SOLVE])
    xe = RE.read_exact(rules, cal, ss, sp, S0.copy(), x0=xf)["_x"]
    Sss = ss_state(ss, cal, sp)
    rss = RuleSet.from_ss(rules.grid, ss, cal)
    rss.n_gh = rules.n_gh
    vs = pieces(PM, Sss, 0, ss_x(ss, cal), rss, cal, ss, sp, True)
    rows = []
    for lab, x in (("SS (pi=0)", None), ("rest, fitted", xf), ("rest, exact", xe)):
        v = vs if x is None else pieces(PM, S0, 0, x, rules, cal, ss, sp, False)
        xx = v.E_Om_D * (1 + v.rdep_D) * v.n_D / v.lev_D
        # the D bond's quarterly return on the no-default branch, over the deposit rate
        exc = v.E_payD_nodef / v.Q_bD - (1 + v.rdep_D) if hasattr(v, "E_payD_nodef") else float("nan")
        rows.append((lab, v, xx, exc))
    print("state at rest vs SS (%): " + "  ".join(f"{n} {100 * (S0[i] / Sss[i] - 1):+.2f}" for i, n in enumerate(STATE_NAMES[:5])))
    print(f"{'':14s} {'mu_D':>8s} {'x':>8s} {'E[Om]':>8s} {'1+rdep':>8s} {'n_D':>8s} {'lev_D':>8s} {'Q_bD':>7s} {'bond exc/yr':>11s} {'alpha_D':>8s}")
    for lab, v, xx, exc in rows:
        print(f"{lab:14s} {v.mu_D:8.5f} {xx:8.5f} {v.E_Om_D:8.4f} {1 + v.rdep_D:8.5f} {v.n_D:8.4f} {v.lev_D:8.4f} "
              f"{v.Q_bD:7.4f} {4e4 * exc:9.1f}bp {v.alpha_D_new:8.4f}")
    _, b, xb, _ = rows[0]
    for lab, v, xx, _ in rows[1:]:
        parts = dict(E_Om=np.log(v.E_Om_D / b.E_Om_D), R=np.log((1 + v.rdep_D) / (1 + b.rdep_D)),
                     n=np.log(v.n_D / b.n_D), lev=-np.log(v.lev_D / b.lev_D))
        print(f"{lab}: log x - log x_ss = {np.log(xx / xb):+.4f} = " + "  ".join(f"{k} {p:+.4f}" for k, p in parts.items()))
    # the bank's book at rest: shares of assets
    v = rows[2][1]
    A = v.Q_D * v.Kp_D + v.Q_bD * v.b_D_D_new + v.p * v.Q_bF * v.b_F_D_new + v.L_wc_D
    print(f"assets at rest: capital {v.Q_D * v.Kp_D / A:.3f}, D bonds {v.Q_bD * v.b_D_D_new / A:.3f}, F bonds "
          f"{v.p * v.Q_bF * v.b_F_D_new / A:.3f}, wc loans {v.L_wc_D / A:.3f}; assets/n {A / v.n_D:.3f}")


if __name__ == "__main__":
    main()
