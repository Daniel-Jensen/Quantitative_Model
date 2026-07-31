"""
Verify that the solved steady state satisfies the multi-asset incentive
constraint at the hardcoded calibration_start['Delta_*'] values.

History. C-1 (Delta_cross > 1, IC degenerate): the old single-asset lambda_gk
formula in steady_auxilliary_D/F implicitly assumed Delta=1 (full divertability)
for every asset, which forced a back-solve to return a degenerate (near- or
above-1) Delta at any realistic portfolio concentration. steady_auxilliary_D/F
now solve lambda_gk from the multi-asset IC directly, taking Delta_bD_D/
Delta_bF_D (D) and Delta_bF_F/Delta_bD_F (F) as genuine calibration inputs.

**Changed 2026-07-31 (collateral mapping).** This module used to *back-solve*
Delta and compare it to the input. That back-solve had one equation and two
unknowns, so it closed the system with a hardcoded

    ratio = Delta_cross / Delta_own = 2.0

which was an undocumented structural convention masquerading as a check: it is
exactly why the inherited 0.2/0.4 pair "passed". The convention is also
infeasible at measured EBA concentration, where GK well-posedness needs
Delta_own > ~0.73 while ratio=2.0 with Delta_cross <= 1 caps Delta_own at 0.5.

Delta_own and Delta_cross are now **free structural parameters**, and this module
checks the thing that actually matters: that the IC binds at the solved SS, i.e.

    value  ==  lambda_gk * [theta - (1-Delta_own)*phi_own - (1-Delta_cross)*phi_cross]

with value = nu_K*kappa + nu_b_own*phi_own + nu_b_cross*phi_cross + eta. This is
a genuine one-equation residual with no free closure and no ratio assumption.

Takes ss_results from solve_steady_state, asserts the IC residual is zero to
tolerance, and returns ss_results unchanged (Delta_* are NOT overwritten -- they
are the calibration).
"""


def ic_residual(phi_own, phi_cross, nu_K, nu_b_own, nu_b_cross, eta, lam, theta,
                Delta_own, Delta_cross):
    """Residual of the binding multi-asset IC at the solved SS.

    Returns (residual, value, divertable_leverage). Zero residual means the
    steady state satisfies the IC with exactly the calibrated Delta pair.
    """
    kappa = theta - phi_own - phi_cross
    value = nu_K * kappa + nu_b_own * phi_own + nu_b_cross * phi_cross + eta
    divertable = theta - (1.0 - Delta_own) * phi_own - (1.0 - Delta_cross) * phi_cross
    return float(value - lam * divertable), float(value), float(divertable)


def calibrate_ic_delta(ss_results):
    ss                = ss_results['ss']
    cal               = ss_results['calibration_start']

    # Country D
    phi_bD_D_ss = float(ss['q_b_D']) * float(ss['b_D_D']) / float(ss['n_inter_D'])
    phi_bF_D_ss = float(ss['q_b_F']) * float(ss['b_F_D']) / float(ss['n_inter_D'])
    res_D, val_D, div_D = ic_residual(
        phi_bD_D_ss, phi_bF_D_ss,
        float(ss['nu_K_D']), float(ss['nu_bD_D']), float(ss['nu_bF_D']), float(ss['eta_D']),
        float(ss['lambda_gk_D']), float(ss['theta_D']),
        float(cal['Delta_bD_D']), float(cal['Delta_bF_D']),
    )

    # Country F
    n_F_ss      = float(ss['n_inter_F']) * float(ss['p'])
    phi_bF_F_ss = float(ss['q_b_F']) * float(ss['b_F_F']) / n_F_ss
    phi_bD_F_ss = float(ss['q_b_D']) * float(ss['b_D_F']) / n_F_ss
    res_F, val_F, div_F = ic_residual(
        phi_bF_F_ss, phi_bD_F_ss,
        float(ss['nu_K_F']), float(ss['nu_bF_F']), float(ss['nu_bD_F']), float(ss['eta_F']),
        float(ss['lambda_gk_F']), float(ss['theta_F']),
        float(cal['Delta_bF_F']), float(cal['Delta_bD_F']),
    )

    print("IC residual check (SS satisfies the multi-asset IC at the calibrated Delta):")
    print(f"  D-bank:  Delta_own = {float(cal['Delta_bD_D']):.4f}  Delta_cross = {float(cal['Delta_bF_D']):.4f}"
          f"   value = {val_D:.6f}  divertable_lev = {div_D:.6f}  residual = {res_D:+.3e}")
    print(f"  F-bank:  Delta_own = {float(cal['Delta_bF_F']):.4f}  Delta_cross = {float(cal['Delta_bD_F']):.4f}"
          f"   value = {val_F:.6f}  divertable_lev = {div_F:.6f}  residual = {res_F:+.3e}")

    TOL = 1e-8
    bad = {c: r for c, r in (('D', res_D), ('F', res_F)) if abs(r) > TOL}
    if bad:
        raise ValueError(
            f"IC residual non-zero (tol={TOL}): {bad}\n"
            "The steady state does not satisfy the multi-asset IC at the calibrated "
            "Delta. This is a regression in steady_auxilliary_D/F's lambda_gk "
            "formula, not a calibration choice to 'fix' by overwriting Delta.")

    # Divertable leverage must be positive, else the IC is economically vacuous
    # (the banker can divert a negative amount).
    for c, d in (('D', div_D), ('F', div_F)):
        if d <= 0.0:
            raise ValueError(
                f"{c}: divertable leverage {d:.6f} <= 0 — IC is vacuous. "
                "Delta is too low for this portfolio concentration; see "
                "steady_state.assert_gk_well_posed.")

    return ss_results
