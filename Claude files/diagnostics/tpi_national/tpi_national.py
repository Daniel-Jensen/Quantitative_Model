# NATIONAL VS UNION DEPOSITS WITH THE TPI: SOLVE, COMPARE, AND DECOMPOSE THE D-BOND PRICE.
#   python3 tpi_national.py --market national --tag nat
#   python3 tpi_national.py --market union --tag uni
#   python3 tpi_national.py --market national --lamZ 0 --tag nat_lz0
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
import os
sys.path.insert(0, str(HERE / os.environ.get("WT", "dff")))  # the patched domestic-funding-fix worktree
OUT = HERE / "tpi_runs"

import numpy as np


def fundamental_spread(cal, sproc, n=600, n_gh=9):
    # THE RISK-NEUTRAL (EXPECTED-LOSS) D-F SPREAD AS A FUNCTION OF s, FROM THE EXOGENOUS p^d(s) ONLY.
    from global_projection.solver_recursive.state_grid import default_prob
    sg = np.linspace(-14.0, 3.0, n)
    x, w = np.polynomial.hermite_e.hermegauss(n_gh)
    w = w / w.sum()
    db, h = cal["delta_b_D"], cal["recovery_rate_D"]
    beta = 1.0 / (1.0 + cal["r_dep_D_target"])
    pd = default_prob(sg)
    q = np.full(n, beta * db / (1 - beta * (1 - db)))
    for _ in range(4000):
        sn = (1 - sproc["rho_s"]) * sproc["s_star"] + sproc["rho_s"] * sg[:, None] + sproc["sigma_s"] * x
        Eq = (np.interp(sn, sg, q) * w).sum(axis=1)
        qn = beta * (1 - pd * (1 - h)) * (db + (1 - db) * Eq)
        if np.max(np.abs(qn - q)) < 1e-14:
            break
        q = qn
    q_rf = beta * db / (1 - beta * (1 - db))
    spread = db * (1 - q) / q - db * (1 - q_rf) / q_rf  # quarterly yield units
    return sg, spread


def setup(args):
    # CALIBRATION, SS AND ANCHORS FOR ONE ECONOMY.
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    cal = get_calibration()
    cal.update(nw_floor_frac=0.15, deposit_market=args.market, national_m_band=False,
               n_jobs=args.n_jobs, tpi_cap_bp=args.cap)
    if args.lamZ is not None:
        cal["lambda_Z_D"] = args.lamZ
    if getattr(args, "phi", 0.0):
        cal["phi_float"] = args.phi
    if getattr(args, "keyD", None) is not None:
        cal["tpi_key_D"] = args.keyD
    if getattr(args, "theta", 0.0):
        cal["theta_cnt"] = args.theta
    if getattr(args, "keyLoss", None) is not None:
        cal["tpi_key_D_loss"] = args.keyLoss
    ss = solve_steady_state(cal)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    if getattr(args, "warranted", False):
        cal["tpi_fund_s"], cal["tpi_fund_spread"] = fundamental_spread(cal, sproc)
    return cal, ss, sproc


def internals(S, x, rules, cal, ss, sproc):
    # THE PERIOD MAP'S INTERNAL OBJECTS AT ONE STATE (SAME SEQUENCE AS point_residuals).
    from global_projection.solver_recursive import point_map as pm
    v = pm._read_state_and_regime(S, 0, x, rules, cal, ss)
    pm._production(v, cal)
    pm._government(v, cal, ss)
    pm._balance_sheets(v, rules, cal)
    pm._continuation(v, rules, cal, sproc, rules.n_gh or 5, False)
    pm._discount_kernels(v, rules, cal)
    pm._incentive_constraint(v, cal)
    pm._bond_demands(v, cal)
    return v


def price_legs(rules, cal, ss, sproc, S):
    # Q = E[Om pay] / (E[Om] R + lam mu): EXPECTED LOSS, RISK PREMIUM, LIQUIDITY PREMIUM (bp/yr).
    from global_projection.solver_recursive.recursive_experiment import read_exact
    o = read_exact(rules, cal, ss, sproc, S)
    v = internals(S, o["_x"], rules, cal, ss, sproc)
    db, rec = cal["delta_b_D"], cal["recovery_rate_D"]
    w, Om, reg = v.wgt, v.Om_D, v.reg
    pay = [(rec if reg[j] else 1.0) * (db + (1 - db) * v.R[j]["Q_bD"]) for j in range(v.nreg)]
    pay_nd = [db + (1 - db) * v.R[j]["Q_bD"] for j in range(v.nreg)]
    E = lambda z: float(sum(np.dot(w[j], z[j]) for j in range(v.nreg)))
    EOm, Epay, Epay_nd = E(Om), E(pay), E(pay_nd)
    EOmpay = E([Om[j] * pay[j] for j in range(v.nreg)])
    cov = EOmpay - EOm * Epay
    R = 1.0 + v.rdep_D
    lam_mu = cal["lambda_bD_D"] * v.mu_D
    wsum = lambda j: float(np.sum(w[j]))
    Om_reg = [float(np.dot(w[j], Om[j]) / wsum(j)) if wsum(j) > 0 else np.nan
              for j in range(v.nreg)]
    rd_reg = [float(np.dot(w[j], v.R[j]["rdep_D"]) / wsum(j)) if wsum(j) > 0 else np.nan
              for j in range(v.nreg)]
    return dict(
        Q_bD=float(o["Q_bD"]), pd=float(sum(wsum(j) for j in range(v.nreg) if reg[j])),
        exp_loss_bp=4e4 * (Epay_nd - Epay) / Epay_nd,
        risk_prem_bp=-4e4 * cov / (EOm * Epay),
        liq_prem_bp=4e4 * lam_mu / (EOm * R),
        Om_nd=Om_reg[0], Om_d=Om_reg[1], Om_ratio=Om_reg[1] / Om_reg[0],
        rdep_now_bp=4e4 * v.rdep_D, rdep_next_nd_bp=4e4 * rd_reg[0], rdep_next_d_bp=4e4 * rd_reg[1],
        mu=float(v.mu_D), b_DD=float(o["b_D_D_new"]), b_DF=float(o["b_D_F_new"]),
        m_cb=float(o["m_cb"]), M_cb=float(o["M_cb_new"]), p=float(o["p"]))


def kernel_report(tag, rules_off, rules_on, c_off, c_on, ss, sproc, pd):
    # THE PRICE LEGS AT REST AND ON IMPACT, IN BOTH ECONOMIES.
    from global_projection.solver_recursive.recursive_experiment import (
        stochastic_rest_point, s_from_pd)
    from global_projection.solver_recursive.state_grid import IS
    rep = {}
    for name, r, c in (("off", rules_off, c_off), ("on", rules_on, c_on)):
        if r is None:
            continue
        S0 = stochastic_rest_point(r, c, ss, sproc, verbose=False)
        S1 = S0.copy()
        S1[IS] = s_from_pd(pd)
        rep[f"rest_{name}"] = price_legs(r, c, ss, sproc, S0)
        rep[f"shock_{name}"] = price_legs(r, c, ss, sproc, S1)
    print(f"\n  [{tag}] D-BOND PRICE LEGS (one-period premia, bp/yr)")
    print("   state        Q_bD   exp.loss  risk.prem  liq.prem  Om_d/Om_nd  rdep  rdep'_nd  rdep'_d    mu    m_cb")
    for k, d in rep.items():
        print(f"   {k:10s} {d['Q_bD']:.4f} {d['exp_loss_bp']:9.1f} {d['risk_prem_bp']:10.1f} "
              f"{d['liq_prem_bp']:9.1f} {d['Om_ratio']:11.4f} {d['rdep_now_bp']:6.0f} "
              f"{d['rdep_next_nd_bp']:9.0f} {d['rdep_next_d_bp']:8.0f} {d['mu']:7.4f} {d['m_cb']:7.4f}")
    return rep


def _clean(x):
    # JSON-SAFE COPY.
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    return x


def main():
    # SOLVE THE NO-TPI BASELINE AND THE TPI ECONOMY, COMPARE, SAVE.
    ap = argparse.ArgumentParser()
    ap.add_argument("--market", default="national")
    ap.add_argument("--lamZ", type=float, default=None)
    ap.add_argument("--phi", type=float, default=0.0)
    ap.add_argument("--keyD", type=float, default=None)
    ap.add_argument("--theta", type=float, default=0.0)
    ap.add_argument("--warranted", action="store_true")
    ap.add_argument("--keyLoss", type=float, default=None)
    ap.add_argument("--cap", type=float, default=200.0)
    ap.add_argument("--s-refine", type=int, default=0)
    ap.add_argument("--n_jobs", type=int, default=0)
    ap.add_argument("--pd", type=float, default=0.02)
    ap.add_argument("--tag", required=True)
    args = ap.parse_args()
    from global_projection.solver_recursive.recursive_experiment import solve_recursive
    from global_projection.solver_recursive import tpi_experiment
    OUT.mkdir(exist_ok=True)
    t0 = time.perf_counter()
    cal, ss, sproc = setup(args)
    print(f"[{args.tag}] market={args.market} lambda_Z={cal.get('lambda_Z_D', 'lambda_bD')} "
          f"cap={args.cap} s_refine={args.s_refine}", flush=True)
    base = []
    rules_off = solve_recursive(dict(cal, tpi_on=False), ss, sproc, s_refine=args.s_refine,
                                base_out=base)
    print(f"[{args.tag}] no-TPI solved ok={rules_off.solve_ok} ({time.perf_counter()-t0:.0f}s)",
          flush=True)
    c_on, c_off = dict(cal, tpi_on=True), dict(cal, tpi_on=False)
    rules_on = solve_recursive(c_on, ss, sproc, s_refine=args.s_refine, base=base[0])
    cap_solved = getattr(rules_on, "tpi_cap_solved", None)
    print(f"[{args.tag}] TPI solved ok={rules_on.solve_ok} cap reached={cap_solved} "
          f"({time.perf_counter()-t0:.0f}s)", flush=True)
    # the comparison at the cap actually reached
    if cap_solved is not None:
        c_on["tpi_cap_bp"] = cap_solved
    out = tpi_experiment.compare(dict(cal, tpi_cap_bp=c_on["tpi_cap_bp"]), ss, sproc,
                                 rules_on, rules_off, args.pd)
    out["kernel"] = kernel_report(args.tag, rules_off, rules_on, c_off, c_on, ss, sproc, args.pd)
    out["args"] = vars(args)
    out["seconds"] = time.perf_counter() - t0
    (OUT / f"{args.tag}.json").write_text(json.dumps(_clean(out), indent=1))
    with open(OUT / f"{args.tag}_rules.pkl", "wb") as fh:
        pickle.dump(dict(on=rules_on, off=rules_off), fh)
    r0, r1 = out["rest_off"], out["rest_on"]
    print(f"\n[{args.tag}] ANNOUNCEMENT (rest point, nothing bought): sov spread "
          f"{r0['sov_bp']:.1f} -> {r1['sov_bp']:.1f} bp/yr, lending spread {r0['spread_bp']:.1f} -> "
          f"{r1['spread_bp']:.1f} bp/yr, M at rest {r1['M']:.2f}% of debt")
    for k in ("impact_rest", "impact_shock"):
        print(f"[{args.tag}] {k}: " + ", ".join(f"m={r['m']:.2f}: Q {r['Q_bD']:.5f} sov {r['sov_bp']:.1f}bp"
                                             for r in out[k]))
    print(f"[{args.tag}] KT: {out['kt']}")
    print(f"[{args.tag}] done in {out['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    # the collocation Jacobian pool re-imports this file in every worker
    main()
