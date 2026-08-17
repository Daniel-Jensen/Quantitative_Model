import copy
import numpy as np
import sequence_jacobian as sj
from sequence_jacobian import simple, combine, create_model

from equations_D import (
    hh_init_D, hh_D, make_grids_D, income_D, hh_extended_D,
    smart_steady_D, market_clearing_D, steady_auxilliary_D,
    banker_div_D, sdf_D, sdf_ss_D, sdf_banker_ss_D, government_ss_D, labor_ss_D,
    government_default_D, bond_price_ss_D, bond_return_D, collateral_quality_D,
    ces_price_D, import_demand_D, deposit_rates_D, deposit_return_D,
    firm_profit_D, price_nkpc_D,
)
from equations_F import (
    hh_init_F, hh_F, make_grids_F, income_F, hh_extended_F,
    smart_steady_F, market_clearing_F, steady_auxilliary_F,
    banker_div_F, sdf_F, sdf_ss_F, sdf_banker_ss_F, government_ss_F, labor_ss_F,
    government_default_F, bond_price_ss_F, bond_return_F, collateral_quality_F,
    ces_price_F, import_demand_F, deposit_rates_F, deposit_return_F,
    firm_profit_F, price_nkpc_F,
)
from equations_global import (
    trade_balance, domestic_bond_clearing,
    portfolio_level_anchors, portfolio_adj_cost, bond_yield,
    global_goods_mkt, external_account_D,
    terms_of_trade, union_inflation,
)
# NB: import the MODULE, not the flag. `from calibration import EBA_CALIBRATION`
# binds the value at import time, so a sweep that flips the switch afterwards
# would silently keep the old portfolio targets — the same stale-binding trap as
# the regimes cache key and PSILAM_MAIN.
import calibration
from calibration import load_eba_targets, load_eba_foreign_shares


def assert_gk_well_posed(ss_in):
    """Hard guard: the Gertler-Karadi block must be economically well-posed.

    `lambda_gk` is the IC multiplier and `Omega` the banker's marginal value of
    net worth; both must be strictly positive. When they are not, the solver
    still converges, every Walras residual is machine-zero, and the IRFs look
    plausible — but the banker's continuation value is negative and nothing
    computed from the block means anything. This is the C-1 failure mode:
    silent degeneracy that passes every check the pipeline previously ran.

    From `steady_auxilliary_D/F`,
        lambda_gk = f / (D_target/(beta_inter*(1+rn)) - (1-f)*theta),
        D_target  = theta - (1-Delta_own)*phi_own - (1-Delta_cross)*phi_cross,
    so lambda_gk > 0 requires a positive denominator, i.e. approximately

        f * theta  >  (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross.

    Read: the banker's franchise value (left) must cover the non-divertable
    ("good collateral") part of the sovereign book (right). High measured
    concentration `phi_own` therefore puts a LOWER BOUND on the divertability
    `Delta_own` — which is how the EBA data partially identifies a parameter
    that has no direct empirical counterpart. See docs/eba_calibration.md.
    """
    vals, bad = [], []
    for c in ("D", "F"):
        lam, om = float(ss_in[f"lambda_gk_{c}"]), float(ss_in[f"Omega_{c}"])
        vals.append(f"{c}: lambda_gk={lam:+.6f}, Omega={om:+.6f}")
        if not (lam > 0.0 and om > 0.0):
            bad.append(c)
    if bad:
        raise ValueError(
            f"GK block is NOT well-posed — negative IC multiplier / franchise "
            f"value in {'/'.join(bad)}:\n  " + "\n  ".join(vals)
            + "\n\nThe steady state will still 'solve' and every Walras residual will be "
              "machine-zero, but the banker's continuation value is negative and all "
              "IRFs are uninterpretable. Required (approximately):\n"
              "    f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross\n"
              "Raise Delta_own (bond divertability), raise f, or lower the sovereign "
              "concentration. See assert_gk_well_posed.__doc__ and "
              "docs/eba_calibration.md 'GK feasibility'.")


def gk_feasibility_margin(theta, f, phi_own, phi_cross, Delta_own, Delta_cross):
    """f*theta - [(1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross]; must be > 0."""
    return f * theta - ((1.0 - Delta_own) * phi_own + (1.0 - Delta_cross) * phi_cross)


def min_Delta_own(theta, f, phi_own, phi_cross, Delta_cross, margin=0.0):
    """Smallest own-bond divertability consistent with a well-posed GK block."""
    return 1.0 - (f * theta - (1.0 - Delta_cross) * phi_cross - margin) / phi_own


def _apply_ss_anchors(ss_in, cal):
    # Fires on every solved SS (this is the common path for steady_state.py and
    # depreciation_calibration.py). Note psi_spread below divides by Omega, so a
    # negative Omega would silently flip the sign of the collateral channel.
    assert_gk_well_posed(ss_in)
    anchors = {
        'phi_bD_D_ss':           float(ss_in['q_b_D']) * float(ss_in['b_D_D']) / float(ss_in['n_inter_D']),
        'phi_bF_F_ss':           float(ss_in['q_b_F']) * float(ss_in['b_F_F']) / (float(ss_in['p']) * float(ss_in['n_inter_F'])),
        'b_F_D_anchor':          float(ss_in['b_F_D']),
        'b_D_F_anchor':          float(ss_in['b_D_F']),
        'excess_return_bD_D_ss': float(ss_in['rb_actual_D']) - float(ss_in['rdep_D']) - cal['T0_D'],
        'excess_return_bF_F_ss': float(ss_in['rb_actual_F']) - float(ss_in['rdep_F']) - cal['T0_F'],
        'excess_return_F_D_ss':  float(ss_in['rb_actual_F']) - float(ss_in['rdep_D']) - cal['T0_D'],
        'excess_return_D_F_ss':  float(ss_in['rb_actual_D']) - float(ss_in['rdep_F']) - cal['T0_F'],
        'psi_spread_F': float(ss_in['lambda_gk_F']) * cal['psi_lambda_B_F']
                        / (float(ss_in['beta_inter_F']) * float(ss_in['Omega_F'])),
        'psi_spread_D': float(ss_in['lambda_gk_D']) * cal['psi_lambda_B_D']
                        / (float(ss_in['beta_inter_D']) * float(ss_in['Omega_D'])),
        # macro-pru-fix: fundamental expected-loss loading per unit default probability,
        # priced by bondholders independent of the psi_lambda_B collateral friction (so it
        # survives psi_lambda_B=0). EL = (1-recovery)*[delta_b + zeta*(1-delta_b)*q_b]/q_b
        # (SS q_b(-1)=q_b). Enters only the bond FOCs (∝ def_rate(+1)=0 at SS) → SS-neutral.
        'EL_price_D': (1.0 - cal['recovery_rate_D'])
                      * (cal['delta_b_D'] + cal['zeta_writeoff_D'] * (1.0 - cal['delta_b_D']) * float(ss_in['q_b_D']))
                      / float(ss_in['q_b_D']),
        'EL_price_F': (1.0 - cal['recovery_rate_F'])
                      * (cal['delta_b_F'] + cal['zeta_writeoff_F'] * (1.0 - cal['delta_b_F']) * float(ss_in['q_b_F']))
                      / float(ss_in['q_b_F']),
        'q_b_D':   float(ss_in['q_b_D']),
        'q_b_F':   float(ss_in['q_b_F']),
        'p':       float(ss_in['p']),
        'C_D_ss':  float(ss_in['C_D']),
        'C_F_ss':  float(ss_in['C_F']),
    }
    cal.update(anchors)
    for k, v in anchors.items():
        ss_in.toplevel[k] = v
    ss_in.toplevel['b_F_D_ss']  = float(ss_in['b_F_D'])
    ss_in.toplevel['b_D_F_ss']  = float(ss_in['b_D_F'])
    ss_in.toplevel['Rgross_D']  = float(1 + ss_in['rdep_D'])
    ss_in.toplevel['Rgross_F']  = float(1 + ss_in['rdep_F'])
    _fr_D = float(ss_in['frisch_D']); _fr_F = float(ss_in['frisch_F'])
    ss_in.toplevel['X_D']   = (float(ss_in['C_D'])
                                - float(ss_in['vphi_D']) * float(ss_in['N_D'])**(1+1/_fr_D) / (1+1/_fr_D))
    ss_in.toplevel['X_F']   = (float(ss_in['C_F'])
                                - float(ss_in['vphi_F']) * float(ss_in['N_F'])**(1+1/_fr_F) / (1+1/_fr_F))
    ss_in.toplevel['U_D']   = ss_in.toplevel['X_D'] / float(ss_in['C_D'])
    ss_in.toplevel['U_F']   = ss_in.toplevel['X_F'] / float(ss_in['C_F'])
    ss_in.toplevel['Phi_D'] = float(ss_in['Phi_D'])
    ss_in.toplevel['Phi_F'] = float(ss_in['Phi_F'])
    ss_in.toplevel['value_D'] = (float(ss_in['beta_inter_D'])
                                  * float(ss_in['Omega_D']) * (1 + float(ss_in['rn_D'])))
    ss_in.toplevel['value_F'] = (float(ss_in['beta_inter_F'])
                                  * float(ss_in['Omega_F']) * (1 + float(ss_in['rn_F'])))
    for k, v in {
        'tau_mp_D': 0.0, 'tau_mp_F': 0.0,
        'T_D': 0.0,  'T_F': 0.0,
        'T_ls_D':   0.0, 'T_ls_F':   0.0,
        'b_F_D_res': 0.0, 'b_D_F_res': 0.0,
        'rb_D_res': 0.0,  'rb_F_res': 0.0,
        'labor_mkt_res_D': 0.0, 'labor_mkt_res_F': 0.0,
        'w_res_D': 0.0, 'w_res_F': 0.0,
    }.items():
        ss_in.toplevel[k] = v
    return anchors


def solve_steady_state(calibration_start):
    ha = sj.create_model([
        sdf_ss_D, sdf_banker_ss_D, government_default_D, bond_price_ss_D, bond_return_D,
        sdf_ss_F, sdf_banker_ss_F, government_default_F, bond_price_ss_F, bond_return_F,
        # Pledgeability map. Nothing in the SS consumes Delta_*_eff_* (the SS uses
        # steady_auxilliary_D/F, not intermediation_IC_D/F), but the DYNAMIC model does,
        # so ss_final must carry them. Listed here rather than injected in
        # _apply_ss_anchors so there is exactly one copy of the algebra -- a second copy
        # is how audit_artifacts/ drifted (CLAUDE.md). SS-neutral: def_rate_ss = 0.
        collateral_quality_D, collateral_quality_F,
        hh_extended_D, smart_steady_D, market_clearing_D, steady_auxilliary_D,
        banker_div_D, government_ss_D, labor_ss_D, firm_profit_D, price_nkpc_D,
        hh_extended_F, smart_steady_F, market_clearing_F, steady_auxilliary_F,
        banker_div_F, government_ss_F, labor_ss_F, firm_profit_F, price_nkpc_F,
        ces_price_D, import_demand_D, ces_price_F, import_demand_F,
        deposit_rates_D, deposit_rates_F,
        deposit_return_D, deposit_return_F,
        bond_yield,
        trade_balance, external_account_D, global_goods_mkt,
        terms_of_trade, union_inflation,
    ], name='MU HA Model 2 Country')

    unknowns_ss = {'beta_D': 0.9850, 'beta_F': 0.9850, 'p': 0.99}
    targets_ss  = ['deposit_mkt_D', 'deposit_mkt_F', 'ca_res_D']

    # ── Initial SS solve ──────────────────────────────────────────────────────
    print("Solving initial steady state...")
    ss = ha.solve_steady_state(calibration_start, unknowns_ss, targets_ss, solver='broyden_custom')

    anchors = {
        'phi_bD_D_ss':           float(ss['q_b_D']) * float(ss['b_D_D']) / float(ss['n_inter_D']),
        'phi_bF_F_ss':           float(ss['q_b_F']) * float(ss['b_F_F']) / (float(ss['p']) * float(ss['n_inter_F'])),
        'b_F_D_anchor':          float(ss['b_F_D']),
        'b_D_F_anchor':          float(ss['b_D_F']),
        'excess_return_bD_D_ss': float(ss['rb_actual_D']) - float(ss['rdep_D']) - calibration_start['T0_D'],
        'excess_return_bF_F_ss': float(ss['rb_actual_F']) - float(ss['rdep_F']) - calibration_start['T0_F'],
        'excess_return_F_D_ss':  float(ss['rb_actual_F']) - float(ss['rdep_D']) - calibration_start['T0_D'],
        'excess_return_D_F_ss':  float(ss['rb_actual_D']) - float(ss['rdep_F']) - calibration_start['T0_F'],
        'q_b_D':                 float(ss['q_b_D']),
        'q_b_F':                 float(ss['q_b_F']),
        'p':                     float(ss['p']),
        'C_D_ss':                float(ss['C_D']),
        'C_F_ss':                float(ss['C_F']),
    }
    calibration_start.update(anchors)
    for k, v in anchors.items():
        ss.toplevel[k] = v

    ss.toplevel['b_F_D_ss'] = float(ss['b_F_D'])
    ss.toplevel['b_D_F_ss'] = float(ss['b_D_F'])
    ss.toplevel['Rgross_D'] = float(1 + ss['rdep_D'])
    ss.toplevel['Rgross_F'] = float(1 + ss['rdep_F'])
    _fr_D = float(ss['frisch_D']); _fr_F = float(ss['frisch_F'])
    ss.toplevel['X_D'] = float(ss['C_D']) - float(ss['vphi_D']) * float(ss['N_D']) ** (1 + 1/_fr_D) / (1 + 1/_fr_D)
    ss.toplevel['X_F'] = float(ss['C_F']) - float(ss['vphi_F']) * float(ss['N_F']) ** (1 + 1/_fr_F) / (1 + 1/_fr_F)
    ss.toplevel['U_D'] = ss.toplevel['X_D'] / float(ss['C_D'])
    ss.toplevel['U_F'] = ss.toplevel['X_F'] / float(ss['C_F'])
    ss.toplevel['Phi_D'] = float(ss['Phi_D'])
    ss.toplevel['Phi_F'] = float(ss['Phi_F'])
    ss.toplevel['value_D'] = float(ss['beta_inter_D']) * float(ss['Omega_D']) * (1 + float(ss['rn_D']))
    ss.toplevel['value_F'] = float(ss['beta_inter_F']) * float(ss['Omega_F']) * (1 + float(ss['rn_F']))
    for k, v in {
        'tau_mp_D': 0.0, 'tau_mp_F': 0.0,
        'T_D': 0.0,  'T_F': 0.0,
        'T_ls_D': 0.0, 'T_ls_F': 0.0,
        'b_F_D_res': 0.0, 'b_D_F_res': 0.0,
        'rb_D_res': 0.0,  'rb_F_res': 0.0,
        'labor_mkt_res_D': 0.0, 'labor_mkt_res_F': 0.0,
        'w_res_D': 0.0, 'w_res_F': 0.0,
    }.items():
        ss.toplevel[k] = v

    # ── Portfolio share targeting ─────────────────────────────────────────────
    # EBA 2011 REBUILD (2026-07-31): bank-sovereign concentration, 31 Dec 2010,
    # read from data/eba_moments.json (single source of truth — do NOT hardcode
    # a second copy here). Regenerate with `python code/eba_calibration.py`.
    #   phi_bD_D = 2.390   GR banks' Greek book / capital  (own, doom loop)
    #   phi_bF_D = 0.018   GR banks' Bund  / capital       (cross)
    #   phi_bD_F = 0.069   DE banks' Greek book / capital  (cross, contagion)
    #   phi_bF_F = 2.758   DE banks' Bund  / capital       (own)
    # Gated by calibration.EBA_CALIBRATION — see the switch note there. At the
    # measured concentration the GK block has no feasible Delta (empty set), so
    # the default is the pre-EBA symmetric placeholder, which solves.
    if calibration.EBA_CALIBRATION:
        print("Targeting portfolio shares (EBA 2011, 31 Dec 2010)...")
        _eba = load_eba_targets()
        target_phi_bD_D = _eba['phi_bD_D_ss']
        target_phi_bF_D = _eba['phi_bF_D_ss']
        target_phi_bD_F = _eba['phi_bD_F_ss']
        target_phi_bF_F = _eba['phi_bF_F_ss']
    else:
        print("Targeting portfolio shares (pre-EBA symmetric placeholder)...")
        target_phi_bD_D = 0.25
        target_phi_bF_D = 0.15
        target_phi_bD_F = 0.15
        target_phi_bF_F = 0.25

    n_D  = float(ss['n_inter_D'])
    n_F  = float(ss['n_inter_F']) * float(ss['p'])
    q_D  = float(ss['q_b_D'])
    q_F  = float(ss['q_b_F'])

    # COUNTRY SIZE (2026-08-07). Each phi is a ratio to its HOLDER's net worth, so
    # every stock below is in its holder's own per-capita units -- which is exactly
    # the convention the model now uses. The aggregation happens in
    # domestic_bond_clearing via size_F, not here.
    size_F = float(calibration_start['size_F'])

    b_D_D_new = target_phi_bD_D * n_D / q_D     # D aggregate
    b_F_D_new = target_phi_bF_D * n_D / q_F     # D aggregate
    b_D_F_new = target_phi_bD_F * n_F / q_D     # per F capita
    b_F_F_new = target_phi_bF_F * n_F / q_F     # per F capita

    # Government stocks, each in ITS OWN country's units:
    #   D debt is a D aggregate;  F debt is per F capita.
    B_D_new   = b_D_D_new + size_F * b_D_F_new
    B_F_new   = b_F_F_new + b_F_D_new / size_F

    # Over-identifying checks. The two EBA moments that the pre-size-asymmetry
    # calibration could not satisfy jointly -- portfolio composition (the phi's,
    # matched by construction above) and market structure (the foreign shares
    # below). Both should now hold. B_supply_*_qgdp and the foreign shares are
    # measured directly and were never used by the targeting, which is why the
    # inconsistency went unnoticed for so long.
    _eba_all = load_eba_targets()
    _fs_eba  = load_eba_foreign_shares()
    fs_D_new = size_F * b_D_F_new / B_D_new
    fs_F_new = (b_F_D_new / size_F) / B_F_new
    print(f"  size_F = {size_F:.4f}  (F/D GDP; every F variable is per F capita)")
    print(f"  foreign-held share of the bank-held sovereign stock:"
          f"  D = {fs_D_new:.4f} (EBA {_fs_eba['D']:.4f})"
          f"   F = {fs_F_new:.6f} (EBA {_fs_eba['F']:.6f})")
    print(f"  bank-held stock in own-country quarterly GDP:"
          f"  B_D = {B_D_new:.4f} (EBA {_eba_all['B_supply_D_qgdp']:.4f})"
          f"   B_F = {B_F_new:.4f} (EBA {_eba_all['B_supply_F_qgdp']:.4f})")

    # omega_K is MEASURED (corporate+CRE EAD / K), not back-solved. It stays at
    # its calibration value; K is then an OUTPUT of the balance sheet,
    #   K = (theta·N - q_b·bonds) / (omega_K·Q),
    # and comparing it to the conventional K/Y_annual=2.7 target (K=10.8 in
    # quarterly-GDP units) is a genuine over-identifying check on the EBA
    # balance sheet. The 2026-07-22 build had this backwards: it ASSUMED
    # theta=4.0 and solved omega_K to force K=10.8, so the check was vacuous and
    # omega_K was a free plug absorbing the theta assumption.
    omega_K_D_new = float(calibration_start['omega_K_D'])
    omega_K_F_new = float(calibration_start['omega_K_F'])

    _fr = float(calibration_start['fund_rule_D'])
    # Built from the REALISED stocks, not `target_phi * n`. Post-units-fix the two
    # differ for the cross-border legs, and using the targets would compute this
    # over-identifying check against the pre-fix bond book — overstating K_implied_F
    # by ~0.13 and quietly making a broken balance sheet look like it validated.
    # Each bank's own balance sheet, in its own per-capita units -- no size_F here:
    # the F bank holds b_D_F per F capita and funds it per F capita.
    _bank_D = (float(calibration_start['theta_D']) * n_D
               - (q_D * b_D_D_new + q_F * b_F_D_new))
    _bank_F = (float(calibration_start['theta_F']) * n_F
               - (q_F * b_F_F_new + q_D * b_D_F_new))
    K_implied_D = ((1 - _fr) * _bank_D / omega_K_D_new
                   + _fr * (_bank_D + float(calibration_start['K_fund_D'])))
    K_implied_F = ((1 - _fr) * _bank_F / omega_K_F_new
                   + _fr * (_bank_F + float(calibration_start['K_fund_F'])))

    print(f"  D-bank: phi_bD_D = {target_phi_bD_D:.3f}  phi_bF_D = {target_phi_bF_D:.3f}"
          f"  omega_K_D = {omega_K_D_new:.4f}")
    print(f"  F-bank: phi_bD_F = {target_phi_bD_F:.3f}  phi_bF_F = {target_phi_bF_F:.3f}"
          f"  omega_K_F = {omega_K_F_new:.4f}")
    print(f"  over-identifying check (K from the measured balance sheet, target 10.8):"
          f"  K_D = {K_implied_D:.3f}   K_F = {K_implied_F:.3f}")

    calibration_start.update({
        'b_D_D': b_D_D_new, 'b_F_D': b_F_D_new,
        'b_D_F': b_D_F_new, 'b_F_F': b_F_F_new,
        'b_F_D_anchor': b_F_D_new, 'b_D_F_anchor': b_D_F_new,
        'phi_bF_D_ss': target_phi_bF_D,
        'omega_K_D': omega_K_D_new, 'omega_K_F': omega_K_F_new,
        'B_supply_D': B_D_new, 'b_gov_D': B_D_new, 'b_gov_ss_D': B_D_new,
        'B_supply_F': B_F_new, 'b_gov_F': B_F_new, 'b_gov_ss_F': B_F_new,
    })

    print("Re-solving SS with new portfolio allocation...")
    _unknowns_warm = {'beta_D': float(ss['beta_D']), 'beta_F': float(ss['beta_F']), 'p': float(ss['p'])}
    ss = ha.solve_steady_state(calibration_start, _unknowns_warm, targets_ss, solver='broyden_custom')
    _apply_ss_anchors(ss, calibration_start)
    print(f"SS re-solved. beta_D={float(ss['beta_D']):.8f}  p={float(ss['p']):.6f}")

    return {
        'ss':                ss,
        'ha':                ha,
        'calibration_start': calibration_start,
        'unknowns_ss':       unknowns_ss,
        'targets_ss':        targets_ss,
    }
