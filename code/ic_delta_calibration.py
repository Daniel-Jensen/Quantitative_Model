"""
Verify the divertable fraction (Delta) implied by the binding IC constraint
against the hardcoded calibration_start['Delta_*'] values.

C-1 (Delta_cross > 1, IC degenerate): the old single-asset lambda_gk formula in
steady_auxilliary_D/F implicitly assumed Delta=1 (full divertability) for every
asset, which forced this back-solve to return a degenerate (near- or above-1)
Delta at any realistic portfolio concentration -- see docs/eba_calibration.md
"Why C-1 forces Delta->1". steady_auxilliary_D/F now solve lambda_gk from the
multi-asset IC directly, taking Delta_bD_D/Delta_bF_D (D) and Delta_bF_F/
Delta_bD_F (F) as hardcoded calibration inputs (0.2 own / 0.4 cross, per
bank-cal). Given that fix, this back-solve is now a redundant re-derivation of
the same equation and should recover the hardcoded inputs almost exactly; a
mismatch here means the SS is not actually satisfying the multi-asset IC and is
a regression, not a calibration choice to "fix" by overwriting.

Takes ss_results from solve_steady_state, asserts back-solved Delta matches the
calibration inputs to within tolerance, and returns the same ss_results dict
unchanged (Delta_bD_D etc. are NOT overwritten -- they are the calibration).
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

    Delta_in = {
        'Delta_bD_D': calibration_start['Delta_bD_D'], 'Delta_bF_D': calibration_start['Delta_bF_D'],
        'Delta_bF_F': calibration_start['Delta_bF_F'], 'Delta_bD_F': calibration_start['Delta_bD_F'],
    }
    Delta_out = {'Delta_bD_D': D_bD_D, 'Delta_bF_D': D_bF_D, 'Delta_bF_F': D_bF_F, 'Delta_bD_F': D_bD_F}

    print("IC Delta consistency check (back-solved vs. hardcoded calibration input):")
    print(f"  D-bank:  Delta_bD_D = {D_bD_D:.4f} (input {Delta_in['Delta_bD_D']:.4f})"
          f"  Delta_bF_D = {D_bF_D:.4f} (input {Delta_in['Delta_bF_D']:.4f})  (value={val_D:.6f})")
    print(f"  F-bank:  Delta_bF_F = {D_bF_F:.4f} (input {Delta_in['Delta_bF_F']:.4f})"
          f"  Delta_bD_F = {D_bD_F:.4f} (input {Delta_in['Delta_bD_F']:.4f})  (value={val_F:.6f})")

    TOL = 1e-6
    mismatches = {k: (Delta_in[k], Delta_out[k]) for k in Delta_in
                  if abs(Delta_in[k] - Delta_out[k]) > TOL}
    if mismatches:
        print(f"  *** REGRESSION: back-solved Delta does not match calibration input "
              f"(tol={TOL}): {mismatches}")
        print("  The SS is not satisfying the multi-asset IC with the hardcoded Delta --")
        print("  check steady_auxilliary_D/F's lambda_gk formula, not this back-solve.")
    if D_bD_D > 1 or D_bF_D > 1 or D_bF_F > 1 or D_bD_F > 1:
        print("  *** WARNING (C-1): back-solved Delta > 1 -- IC constraint is degenerate.")

    return ss_results
