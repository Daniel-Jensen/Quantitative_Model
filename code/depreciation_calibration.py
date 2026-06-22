"""
Capital depreciation rate calibration and final steady-state re-solve.

Targets rk = 0.01 per quarter for both countries, then does one final
SS solve with the calibrated delta values. Also runs the post-SS
residual diagnostic.
"""
import copy

from steady_state import _apply_ss_anchors


def calibrate_depreciation(ss_results):
    ss                = ss_results['ss']
    ha                = ss_results['ha']
    calibration_start = ss_results['calibration_start']
    unknowns_ss       = ss_results['unknowns_ss']
    targets_ss        = ss_results['targets_ss']

    rk_D_target = 0.01
    rk_F_target = 0.01

    K_D_cur = float(ss['K_D']); K_F_cur = float(ss['K_F'])
    Y_D_cur = float(ss['Y_D']); Y_F_cur = float(ss['Y_F'])

    delta_D_cal = calibration_start['alpha_D'] * Y_D_cur / K_D_cur - rk_D_target
    delta_F_cal = calibration_start['alpha_F'] * Y_F_cur / K_F_cur - rk_F_target
    calibration_start.update({'delta_D': delta_D_cal, 'delta_F': delta_F_cal})

    print(f"Depreciation calibration:  delta_D = {delta_D_cal:.6f}  delta_F = {delta_F_cal:.6f}")
    print("Final SS re-solve with calibrated delta...")

    ss = ha.solve_steady_state(calibration_start, unknowns_ss, targets_ss, solver='broyden_custom')
    _apply_ss_anchors(ss, calibration_start)

    print(f"Verified rk_D = {float(ss['rk_D']):.6f}  (target {rk_D_target:.4f})")
    print(f"Verified rk_F = {float(ss['rk_F']):.6f}  (target {rk_F_target:.4f})")
    print(f"Final beta_D  = {float(ss['beta_D']):.10f}")
    print(f"Final beta_F  = {float(ss['beta_F']):.10f}")
    print(f"\nbeta_D={ss['beta_D']:.10f}  beta_F={ss['beta_F']:.10f}  p={ss['p']:.6f}")
    print(f"rb_D={ss['rb_D']:.6f}  rb_F={ss['rb_F']:.6f}  rdep_D={ss['rdep_D']:.6f}  rdep_F={ss['rdep_F']:.6f}")
    print(f"q_b_D={float(ss['q_b_D']):.6f}  q_b_F={float(ss['q_b_F']):.6f}")
    print("SS goods residuals:")
    print("  goods_mkt_D =", ss['goods_mkt_D'])
    print("  goods_mkt_F =", ss['goods_mkt_F'])
    print("  ca_res_D    =", ss['ca_res_D'])

    cali_D = cali_F = ss
    ss_final = copy.deepcopy(ss)

    _run_ss_residual_diagnostic(ss, calibration_start)

    return {
        **ss_results,
        'ss':       ss,
        'ss_final': ss_final,
        'cali_D':   cali_D,
        'cali_F':   cali_F,
    }


def _run_ss_residual_diagnostic(ss, calibration_start):
    def _get(k):
        return float(ss[k])

    diag = {}

    for c in ['D', 'F']:
        pdiv    = _get('p') if c == 'F' else 1.0
        eta_c   = _get(f'eta_{c}')
        lam     = _get(f'lambda_gk_{c}')
        theta_c = _get(f'theta_{c}')
        Q_c     = _get(f'Q_{c}')
        K_c     = _get(f'K_{c}')
        n_c     = _get(f'n_inter_{c}')
        kappa_c = Q_c * K_c / n_c
        if c == 'D':
            nu_K, nu_bD, nu_bF = _get('nu_K_D'), _get('nu_bD_D'), _get('nu_bF_D')
            q_h, q_x = _get('q_b_D'), _get('q_b_F')
            b_h, b_x = _get('b_D_D'), _get('b_F_D')
            Dh, Dx   = _get('Delta_bD_D'), _get('Delta_bF_D')
        else:
            nu_K, nu_bD, nu_bF = _get('nu_K_F'), _get('nu_bD_F'), _get('nu_bF_F')
            q_h, q_x = _get('q_b_F'), _get('q_b_D')
            b_h, b_x = _get('b_F_F'), _get('b_D_F')
            Dh, Dx   = _get('Delta_bF_F'), _get('Delta_bD_F')
        phi_h = q_h * b_h / (pdiv * n_c)
        phi_x = q_x * b_x / (pdiv * n_c)
        value_c = (nu_K * kappa_c
                   + (nu_bD if c == 'D' else nu_bF) * phi_h
                   + (nu_bF if c == 'D' else nu_bD) * phi_x
                   + eta_c)
        theta_tgt = value_c / lam + (1 - Dh) * phi_h + (1 - Dx) * phi_x
        diag[f'IC_{c}: θ − θ_tgt'] = theta_c - theta_tgt

    for c in ['D', 'F']:
        f_c     = _get(f'f_{c}')
        lam     = _get(f'lambda_gk_{c}')
        beta_c  = _get(f'beta_inter_{c}')
        rk_c    = _get(f'rk_{c}')
        rdep_c  = _get(f'rdep_{c}')
        eta_c   = _get(f'eta_{c}')
        theta_c = _get(f'theta_{c}')
        Omega_p1 = f_c + (1 - f_c) * lam * theta_c
        rb_h = _get('rb_actual_D' if c == 'D' else 'rb_actual_F')
        rb_x = _get('rb_actual_F' if c == 'D' else 'rb_actual_D')
        nu_K_c  = _get(f'nu_K_{c}')
        nu_bh_c = _get('nu_bD_D' if c == 'D' else 'nu_bF_F')
        nu_bx_c = _get('nu_bF_D' if c == 'D' else 'nu_bD_F')
        diag[f'P1_{c}: nu_K_res']  = nu_K_c  - beta_c * Omega_p1 * (rk_c - rdep_c)
        diag[f'P1_{c}: nu_bh_res'] = nu_bh_c - beta_c * Omega_p1 * (rb_h - rdep_c)
        diag[f'P1_{c}: nu_bx_res'] = nu_bx_c - beta_c * Omega_p1 * (rb_x - rdep_c)
        diag[f'P1_{c}: eta_res']   = eta_c   - beta_c * Omega_p1 * (1 + rdep_c)

    diag['ca_res_D'] = _get('ca_res_D')

    TOL = 1e-8
    print(f"\n{'Block residual':<55} {'Value':>14}  Status")
    print("-" * 85)
    FLAGGED = []
    for name, val in diag.items():
        ok = abs(val) <= TOL
        if not ok:
            FLAGGED.append(name)
        print(f"  {name:<53} {val:>14.6e}  {'OK' if ok else '*** FAIL'}")
    print("-" * 85)
    print("All residuals < 1e-8  ✓" if not FLAGGED else f"FLAGGED: {FLAGGED}")
