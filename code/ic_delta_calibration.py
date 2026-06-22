"""
Back-solve divertable fraction (Delta) from the binding IC constraint.

Takes ss_results from solve_steady_state, updates calibration_start['Delta_*']
in-place, and returns the same ss_results dict.
"""


def _ic_delta(phi_own, phi_cross, nu_K, nu_b_own, nu_b_cross, eta, lam, theta, ratio):
    kappa     = theta - phi_own - phi_cross
    value     = nu_K * kappa + nu_b_own * phi_own + nu_b_cross * phi_cross + eta
    denom     = phi_own + ratio * phi_cross
    delta_own = (phi_own + phi_cross - (theta - value / lam)) / denom
    return float(delta_own), float(ratio * delta_own), float(value)


def calibrate_ic_delta(ss_results):
    ss                = ss_results['ss']
    calibration_start = ss_results['calibration_start']

    ratio_D = ratio_F = 2.0

    # Country D
    phi_bD_D_ss = float(ss['q_b_D']) * float(ss['b_D_D']) / float(ss['n_inter_D'])
    phi_bF_D_ss = float(ss['q_b_F']) * float(ss['b_F_D']) / float(ss['n_inter_D'])
    D_bD_D, D_bF_D, val_D = _ic_delta(
        phi_bD_D_ss, phi_bF_D_ss,
        float(ss['nu_K_D']), float(ss['nu_bD_D']), float(ss['nu_bF_D']), float(ss['eta_D']),
        float(ss['lambda_gk_D']), float(ss['theta_D']), ratio_D,
    )

    # Country F
    n_F_ss      = float(ss['n_inter_F']) * float(ss['p'])
    phi_bF_F_ss = float(ss['q_b_F']) * float(ss['b_F_F']) / n_F_ss
    phi_bD_F_ss = float(ss['q_b_D']) * float(ss['b_D_F']) / n_F_ss
    D_bF_F, D_bD_F, val_F = _ic_delta(
        phi_bF_F_ss, phi_bD_F_ss,
        float(ss['nu_K_F']), float(ss['nu_bF_F']), float(ss['nu_bD_F']), float(ss['eta_F']),
        float(ss['lambda_gk_F']), float(ss['theta_F']), ratio_F,
    )

    calibration_start.update({
        'Delta_bD_D': D_bD_D, 'Delta_bF_D': D_bF_D,
        'Delta_bF_F': D_bF_F, 'Delta_bD_F': D_bD_F,
    })

    print("IC Delta calibration:")
    print(f"  D-bank:  Delta_bD_D = {D_bD_D:.4f}  Delta_bF_D = {D_bF_D:.4f}  (value={val_D:.6f})")
    print(f"  F-bank:  Delta_bF_F = {D_bF_F:.4f}  Delta_bD_F = {D_bD_F:.4f}  (value={val_F:.6f})")
    if D_bD_D > 1 or D_bF_D > 1 or D_bF_F > 1 or D_bD_F > 1:
        print("  WARNING (C-1): back-solved Delta > 1 — IC constraint is degenerate.")
        print("  Consider hardcoding Delta_D=0.2, Delta_F=0.4 per bank-cal branch.")

    return ss_results
