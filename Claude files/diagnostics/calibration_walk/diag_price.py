# WHY IS THE D BOND DEAR AT A LOW PRICED DEFAULT PROBABILITY? THE BOND FOC TERM BY TERM ON A SAVED SOLVE (2026-10-23).
# The D bank's FOC is Q = E[Om payD] / (E[Om](1 + rdep) + lambda_bD mu), payD = hc (delta + (1 - delta) Q'). Written as a
# required return, Q = Ehat[payD] / (1 + rdep + lambda_bD mu / E[Om]) with Ehat the Om-weighted mean. The script prints each
# term at the s nodes (other states at the box centre and at the rest point's bank state), at the rest point (fitted and
# exact) and at the SS; the price a stationary perpetuity would have at that required return; and the pieces of Ehat[payD]:
# the plain mean of Q', its Om tilt, and the default leg.
# usage: python3 diag_price.py runs/pd002.pkl
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def terms(PM, v, cal, sp):
    # THE FOC'S PIECES AT ONE EVALUATED POINT.
    import numpy as np
    from global_projection.solver_recursive.state_grid import default_prob
    db = cal["delta_b_D"]
    eps, wq = PM.gh_nodes(len(v.wgt[0]))
    Sn, _ = PM._next_states(v, sp, eps)
    nd = [j for j in range(v.nreg) if not v.reg[j]][0]
    Qn = np.asarray(v.R[nd]["Q_bD"], dtype=float)
    Om = np.asarray(v.Om_D[nd], dtype=float)
    pd = float(default_prob(v.s))
    EQ = float(wq @ Qn)
    EQ_om = float((wq * Om) @ Qn / (wq @ Om))
    wedge = cal["lambda_bD_D"] * v.mu_D / v.E_Om_D
    req = (1 + v.rdep_D + wedge)
    pay_hat = v.E_Om_payD / v.E_Om_D
    # stationary perpetuity: Q = hc_bar (db + (1 - db) Q) / req, hc_bar the Om-weighted repayment share
    gross = sum(float(np.dot(v.wgt[j], v.Om_D[j] * (db + (1 - db) * np.asarray(v.R[j]["Q_bD"])))) for j in range(v.nreg))
    hc_bar = v.E_Om_payD / gross
    Q_stat = hc_bar * db / (req - hc_bar * (1 - db))
    return dict(s=v.s, pd=pd, Q=v.Q_bD, foc=v.foc_D if hasattr(v, "foc_D") else (v.E_Om_payD - v.dmd_D * v.Q_bD) / v.E_Om_D,
                rdep=4e4 * v.rdep_D, mu=v.mu_D, wedge=4e4 * wedge, EQ=EQ, EQ_om=EQ_om,
                sn_lo=float(Sn[0, PM.IS]), sn_hi=float(Sn[-1, PM.IS]), Qn_lo=float(Qn[0]), Qn_hi=float(Qn[-1]),
                hc_bar=hc_bar, Q_implied=pay_hat / req, Q_stat=Q_stat,
                yld=4e4 * db * (1 - v.Q_bD) / v.Q_bD, req_ann=4e4 * (req - 1))


def show(lab, t):
    print(f"{lab:22s} s {t['s']:7.3f} p^d {100 * t['pd']:7.4f}%  Q {t['Q']:.4f} (FOC {t['Q_implied']:.4f}, res {t['foc']:+.1e}) "
          f"yield {t['yld']:6.1f}  rdep {t['rdep']:6.1f}  lam.mu/EOm {t['wedge']:6.1f}  req {t['req_ann']:6.1f}bp | "
          f"E[Q'] {t['EQ']:.4f} Om-tilt {t['EQ_om'] - t['EQ']:+.4f}  Q'(s'={t['sn_lo']:.2f}) {t['Qn_lo']:.4f} "
          f"Q'(s'={t['sn_hi']:.2f}) {t['Qn_hi']:.4f}  hc_bar {t['hc_bar']:.5f}  stationary Q {t['Q_stat']:.4f}", flush=True)


def main():
    # EVALUATE THE FOC ACROSS s AND AT THE REST POINT.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors, ss_state, ss_x
    from global_projection.solver_recursive.decision_rules import RuleSet, SOLVE
    import global_projection.solver_recursive.point_map as PM
    import global_projection.solver_recursive.recursive_experiment as RE
    from diag_slack import pieces
    d_ = pickle.load(open(os.path.join(HERE, sys.argv[1]), "rb"))
    rules = d_["rules"]
    cal = get_calibration()
    cal["nw_floor_frac"] = d_["nw_floor_frac"]
    cal.update(d_.get("cal_over", {}))
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    sp.update(d_.get("sp_over", {}))
    calibrate_household_anchors(cal, ss, sp)
    g = rules.grid
    print(f"s_star {sp['s_star']:.3f} (p^d {100 / (1 + np.exp(-sp['s_star'])):.4f}%), rho_s {sp['rho_s']}, sigma_s {sp['sigma_s']}; "
          f"box s [{g.lo[IS]:.3f}, {g.hi[IS]:.3f}] = p^d [{100 / (1 + np.exp(-g.lo[IS])):.4f}, {100 / (1 + np.exp(-g.hi[IS])):.3f}]%; "
          f"Q_ss {ss['ss_bank_D']['Q_bD'] if 'Q_bD' in ss['ss_bank_D'] else float('nan')}; points {g.n}", flush=True)
    S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)
    xf = np.array([float(rules.eval_read(k, 0, S0[None])[0]) for k in SOLVE])
    xe = RE.read_exact(rules, cal, ss, sp, S0.copy(), x0=xf)["_x"]
    Sss = ss_state(ss, cal, sp)
    rss = RuleSet.from_ss(rules.grid, ss, cal)
    rss.n_gh = rules.n_gh
    show("SS (pi=0, SS rules)", terms(PM, pieces(PM, Sss, 0, ss_x(ss, cal), rss, cal, ss, sp, True), cal, sp))
    show("rest, fitted", terms(PM, pieces(PM, S0, 0, xf, rules, cal, ss, sp, False), cal, sp))
    show("rest, exact", terms(PM, pieces(PM, S0, 0, xe, rules, cal, ss, sp, False), cal, sp))
    print(f"rest point u: {np.array2string(g.to_unit(S0)[0], precision=2)}", flush=True)
    for base, Sb in (("centre", g.from_unit(np.zeros((1, g.d)))[0]), ("rest bank", S0)):
        for u in sorted(set(np.round(g.points_unit[:, IS], 6)) | {-0.5, 0.5}):
            S = Sb.copy()
            S[IS] = g.from_unit(np.full((1, g.d), u))[0][IS]
            x = np.array([float(rules.eval_read(k, 0, S[None])[0]) for k in SOLVE])
            show(f"{base} u_s={u:+.2f}", terms(PM, pieces(PM, S, 0, x, rules, cal, ss, sp, False), cal, sp))


if __name__ == "__main__":
    main()
