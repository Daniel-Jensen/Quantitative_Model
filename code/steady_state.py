import copy
import numpy as np
import sequence_jacobian as sj
from sequence_jacobian import simple, combine, create_model

from equations_D import (
    hh_init_D, hh_D, make_grids_D, income_D, hh_extended_D,
    smart_steady_D, market_clearing_D, steady_auxilliary_D,
    banker_div_D, sdf_D, sdf_ss_D, sdf_banker_ss_D, government_ss_D, labor_ss_D,
    government_default_D, bond_price_ss_D, bond_return_D,
    ces_price_D, import_demand_D, deposit_return_D,
)
from equations_F import (
    hh_init_F, hh_F, make_grids_F, income_F, hh_extended_F,
    smart_steady_F, market_clearing_F, steady_auxilliary_F,
    banker_div_F, sdf_F, sdf_ss_F, sdf_banker_ss_F, government_ss_F, labor_ss_F,
    government_default_F, bond_price_ss_F, bond_return_F,
    ces_price_F, import_demand_F, deposit_return_F,
)
from equations_global import (
    trade_balance, domestic_bond_clearing,
    portfolio_level_anchors, portfolio_adj_cost, bond_yield,
    global_goods_mkt, external_account_D,
)


def _apply_ss_anchors(ss_in, cal):
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
        hh_extended_D, smart_steady_D, market_clearing_D, steady_auxilliary_D,
        banker_div_D, government_ss_D, labor_ss_D,
        hh_extended_F, smart_steady_F, market_clearing_F, steady_auxilliary_F,
        banker_div_F, government_ss_F, labor_ss_F,
        ces_price_D, import_demand_D, ces_price_F, import_demand_F,
        deposit_return_D, deposit_return_F,
        bond_yield,
        trade_balance, external_account_D, global_goods_mkt,
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
    print("Targeting portfolio shares...")
    target_phi_bD_D = 0.25
    target_phi_bF_D = 0.15
    target_phi_bD_F = 0.15
    target_phi_bF_F = 0.25

    n_D  = float(ss['n_inter_D'])
    n_F  = float(ss['n_inter_F']) * float(ss['p'])
    q_D  = float(ss['q_b_D'])
    q_F  = float(ss['q_b_F'])

    b_D_D_new = target_phi_bD_D * n_D / q_D
    b_F_D_new = target_phi_bF_D * n_D / q_F
    b_D_F_new = target_phi_bD_F * n_F / q_D
    b_F_F_new = target_phi_bF_F * n_F / q_F
    B_D_new   = b_D_D_new + b_D_F_new
    B_F_new   = b_F_D_new + b_F_F_new

    print(f"  D-bank: phi_bD_D = {target_phi_bD_D:.3f}  phi_bF_D = {target_phi_bF_D:.3f}")
    print(f"  F-bank: phi_bD_F = {target_phi_bD_F:.3f}  phi_bF_F = {target_phi_bF_F:.3f}")

    calibration_start.update({
        'b_D_D': b_D_D_new, 'b_F_D': b_F_D_new,
        'b_D_F': b_D_F_new, 'b_F_F': b_F_F_new,
        'b_F_D_anchor': b_F_D_new, 'b_D_F_anchor': b_D_F_new,
        'phi_bF_D_ss': target_phi_bF_D,
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
