# APPENDIX-B MEASUREMENT RUNS (2026-10-26): SOLVE THE GLOBAL MODEL AT THE PORTED CALIBRATION, THEN READ
# SOLVABILITY, THE REST POINT AND THE RISK-SHOCK SIGNS (FITTED AND EXACT) WITH THE E[OMEGA] DECOMPOSITION.
#   python3 "Claude files/diagnostics/tex_calibration/measure.py" --tag m1
#   ... --tag m4 --floor 0.5 --instruments runs/m2_instruments.json
#   ... --tag m3 --s-refine 9 --instruments runs/m2_instruments.json
# Writes runs/<tag>_rules.pkl and runs/<tag>_summary.json; the log is the caller's stdout.
import argparse
import json
import os
import pickle
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

import numpy as np

PDS = (0.0025, 0.005, 0.01, 0.02)


def omega_legs(rules, cal, ss, sproc, S, read):
    # THE D BANK'S KERNEL BY REGIME AT STATE S, FROM THE PERIOD MAP CLEARED BY read.
    import global_projection.solver_recursive.point_map as pm
    o = read(rules, cal, ss, sproc, S.copy())
    v = pm._read_state_and_regime(S, 0, o["_x"], rules, cal, ss)
    pm._production(v, cal); pm._government(v, cal, ss); pm._balance_sheets(v, rules, cal)
    pm._continuation(v, rules, cal, sproc, rules.n_gh or 5, False)
    pm._discount_kernels(v, rules, cal); pm._incentive_constraint(v, cal)
    w = [float(np.sum(v.wgt[j])) for j in range(v.nreg)]
    leg = lambda arr, j: float(np.dot(v.wgt[j], arr)) / max(w[j], 1e-300)
    return dict(pd=w[1], E_Om=v.E_Om_D, Om_nd=leg(v.Om_D[0], 0), Om_d=leg(v.Om_D[1], 1),
                alpha_nd=leg(v.R[0]["alpha_D"], 0), alpha_d=leg(v.R[1]["alpha_D"], 1),
                n_D=v.n_D, lev_D=v.lev_D, mu_D=v.mu_D, resid=o["_resid"], o=o)


def row(o, base):
    # THE ACCEPTANCE OBJECTS AS DEVIATIONS FROM THE REST-POINT READ.
    pct = lambda k: 100.0 * (o[k] / base[k] - 1.0)
    return dict(Q_bD=pct("Q_bD"), n_D=pct("n_D"), n_F=pct("n_F"), Y_D=pct("Y_D"), C_D=pct("C_D"),
                spread_bp=4e4 * ((o["r_wc_D"] - o["rdep_D"]) - (base["r_wc_D"] - base["rdep_D"])),
                mu_D=o["mu_D"], goods_F=o["goods_F"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--floor", type=float, default=0.15)
    ap.add_argument("--s-refine", type=int, default=0)
    ap.add_argument("--instruments", default=None)
    a = ap.parse_args()

    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES, IS
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.calibrate_stochastic import rest_point_moments

    t0 = time.perf_counter()
    cal = get_calibration()
    cal["nw_floor_frac"] = a.floor
    if a.instruments:
        cal.update(json.loads((HERE / a.instruments).read_text()))
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print(f"[{a.tag}] SS: lambda {cal['lambda_K_D']:.5f}  mu_ss {ss['ss_bank_D']['mu_ss']:.5f}  "
          f"q_D {ss['Q_bD_ss']:.4f}  p_ss {ss['p_ss']:.5f}  floor {a.floor}", flush=True)

    rules = RE.solve_recursive(cal, ss, sproc, s_refine=a.s_refine)
    with open(HERE / "runs" / f"{a.tag}_rules.pkl", "wb") as fh:
        pickle.dump(dict(rules=rules, cal=cal, ss=ss, sproc=sproc), fh)
    print(f"[{a.tag}] solve: ok {rules.solve_ok}  worst {rules.solve_worst:.2e}  smoothing "
          f"{rules.smoothing}  ({time.perf_counter() - t0:.0f}s)", flush=True)

    S0 = RE.stochastic_rest_point(rules, cal, ss, sproc, verbose=True)
    u = rules.grid.to_unit(np.atleast_2d(S0))[0]
    out = dict(tag=a.tag, floor=a.floor, solve_ok=bool(rules.solve_ok),
               solve_worst=float(rules.solve_worst), smoothing=rules.smoothing,
               rest_state=dict(zip(STATE_NAMES, map(float, S0))),
               rest_u=dict(zip(STATE_NAMES, map(float, u))), reads={})
    for name, read in (("fitted", RE.read_at), ("exact", RE.read_exact)):
        base = omega_legs(rules, cal, ss, sproc, S0, read)
        m = rest_point_moments(base["o"])
        rec = dict(rest={k: v for k, v in base.items() if k != "o"},
                   rest_moments={k: float(v) for k, v in m.items()}, shocks={})
        for pd in PDS:
            S = S0.copy(); S[IS] = RE.s_from_pd(pd)
            leg = omega_legs(rules, cal, ss, sproc, S, read)
            rec["shocks"][str(pd)] = dict(row(leg["o"], base["o"]),
                                          **{k: v for k, v in leg.items() if k != "o"})
        out["reads"][name] = rec
    RE.impact_table(rules, cal, ss, sproc)
    irf = RE.dynamic_irf(rules, cal, ss, sproc, pd_shock=0.02, T=25)
    out["irf_at_wall"] = [bool(x) for x in np.atleast_1d(irf.get("at_wall", []))]
    (HERE / "runs" / f"{a.tag}_summary.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: out[k] for k in ("solve_ok", "solve_worst", "smoothing", "rest_u")},
                     indent=1, default=float))
    for name, rec in out["reads"].items():
        print(f"\n[{a.tag}] {name} read: rest moments {rec['rest_moments']}")
        print(f"   rest: mu_D {rec['rest']['mu_D']:.5f}  E[Om] {rec['rest']['E_Om']:.4f}")
        for pd, r in rec["shocks"].items():
            print(f"   p^d {pd}: dQ_bD {r['Q_bD']:+.2f}%  dn_D {r['n_D']:+.2f}%  dn_F {r['n_F']:+.2f}%  "
                  f"dY_D {r['Y_D']:+.4f}%  dC_D {r['C_D']:+.4f}%  dspread {r['spread_bp']:+.1f}bp  "
                  f"mu_D {r['mu_D']:.5f}  E[Om] {r['E_Om']:.4f} = nd {r['Om_nd']:.4f} / d {r['Om_d']:.3f} "
                  f"(alpha' nd {r['alpha_nd']:.3f}, d {r['alpha_d']:.3f})  goods_F {r['goods_F']:+.1e}")
    print(f"\n[{a.tag}] done in {time.perf_counter() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
