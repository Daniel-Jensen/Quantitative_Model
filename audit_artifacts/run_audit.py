"""Forensic audit runner: replicates model_v12.ipynb pipeline, then runs diagnostics.

Replicates cells 2,3,10,11,12,13,18,20,21 exactly (same order, same values),
then computes:
  A. SS residuals incl. untargeted ones (goods_mkt_D/F, ca_res_D, global_goods_res)
  B. IC / Bellman residuals at the linearization point ss_final
  C. IRF paths of untargeted residuals (Walras leak test) for TFP and default shocks
  D. Leak decomposition: factor-income timing gap and F-bank ToT revaluation gap
  E. Jacobian conditioning of the 23x23 target/unknown system
Results -> audit_artifacts/run1_results.json, audit_artifacts/irfs.npz
"""
import json, sys, copy
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import os
os.chdir(ROOT)

import sequence_jacobian as sj
from sequence_jacobian import combine

RES = {}

def jlog(key, val):
    RES[key] = val
    print(f"[audit] {key}: {val}")

# ═══════════════ CELL 2-3: calibration ═══════════════
calibration_start = {
    'frisch_D': 0.50, 'frisch_F': 0.50, 'eis_D': 0.5, 'eis_F': 0.5,
    'rdep_D': 0.000, 'rdep_F': 0.000, 'q_b_D': 0.83, 'q_b_F': 0.83,
    'Q_D': 1.0, 'Q_F': 1.0, 'alpha_D': 0.35, 'alpha_F': 0.35,
    'delta_D': 0.025, 'delta_F': 0.025, 'ksi_D': 0.50, 'ksi_F': 0.50,
    'delta_b_D': 0.10, 'delta_b_F': 0.10,
    'Y_D': 1.00, 'Y_F': 1.00, 'Y_ss_D': 1.0, 'Y_ss_F': 1.0,
    'N_D': 1.00, 'N_F': 1.00, 'w_D': 0.65, 'w_F': 0.65,
    'f_D': 0.12, 'f_F': 0.12, 'lambda_gk_D': 0.2, 'lambda_gk_F': 0.2,
    'beta_inter_D': 0.9975155088, 'beta_inter_F': 0.9975155088,
    'Delta_bD_D': 0.2, 'Delta_bF_F': 0.2, 'Delta_bF_D': 0.4, 'Delta_bD_F': 0.4,
    'lambda_BD_D': 0.06, 'lambda_BF_F': 0.06, 'lambda_BF_D': 0.06, 'lambda_BD_F': 0.06,
    'psi_lambda_B_D': 3.0, 'psi_lambda_B_F': 3.0,
    'n_inter_D': 0.75*4, 'n_inter_F': 0.75*4, 'theta_D': 4, 'theta_F': 4,
    'psi_nu_bD_D': 0.0, 'psi_nu_bD_F': 0.0, 'psi_nu_bF_D': 0.0, 'psi_nu_bF_F': 0.0,
    'B_supply_D': 0.6*4, 'B_supply_F': 0.6*4,
    'b_gov_D': 0.6*4, 'b_gov_F': 0.6*4, 'b_gov_ss_D': 0.6*4, 'b_gov_ss_F': 0.6*4,
    'tau_D': 0.181, 'tau_F': 0.181, 'lamb_D': 0.85, 'lamb_F': 0.85,
    'lamb_ss_D': 0.85, 'lamb_ss_F': 0.85, 'phi_lamb_D': 0.15, 'phi_lamb_F': 0.15,
    'shock_def_D': 0.000, 'shock_def_F': 0.0, 'T_ls_D': 0.000, 'T_ls_F': 0.000,
    'def_rate_D': 0.000, 'def_rate_F': 0.0, 'def_scale_D': 0.25, 'def_scale_F': 0.25,
    'def_curvature_D': 0.5, 'def_curvature_F': 0.5, 'def_offset_D': 0.05, 'def_offset_F': 0.05,
    'recovery_rate_D': 0.00, 'recovery_rate_F': 0.00,
    'zeta_writeoff_D': 0.0, 'zeta_writeoff_F': 0.0,
    'writeoff_enabled_D': 0.0, 'writeoff_enabled_F': 0.0,
    'chi0_D': 0.00, 'chi0_F': 0.00, 'chi1_D': 0.00, 'chi1_F': 0.00,
    'chi2_D': 2.0, 'chi2_F': 2.0,
    'T0_D': 0.000, 'T0_F': 0.000, 'T1_D': 0.0, 'T1_F': 0.0,
    'omega': 0.85, 'epsilon_trade': 1.5, 'p': 0.50,
    'phi_bF_D_ss': 0.25, 'phi_bD_F_ss': 0.25, 'psi_bF_D': 0.5, 'psi_bD_F': 0.5,
    'mu_w_D': 1.0, 'mu_w_F': 1.0, 'mc_D': 1.0, 'mc_F': 1.0,
    'rho_z_D': 0.90, 'rho_z_F': 0.90, 'sigma_z_D': 0.3, 'sigma_z_F': 0.3,
    'nZ_D': 15, 'nZ_F': 15, 'nDep_D': 500, 'nDep_F': 500,
    'Depmax_D': 150, 'Depmax_F': 150,
}
_n_D, _n_F = calibration_start['n_inter_D'], calibration_start['n_inter_F']
_B_D, _B_F = calibration_start['B_supply_D'], calibration_start['B_supply_F']
b_F_D = calibration_start['phi_bF_D_ss'] * _n_D / calibration_start['q_b_F']
b_D_F = calibration_start['phi_bD_F_ss'] * _n_F / calibration_start['q_b_D']
calibration_start.update({
    'b_F_D': b_F_D, 'b_D_F': b_D_F,
    'b_D_D': _B_D - b_D_F, 'b_F_F': _B_F - b_F_D,
    'b_F_D_anchor': b_F_D, 'b_D_F_anchor': b_D_F,
    'psi_bD_D': 0.0, 'psi_bF_F': 0.0,
})

# ═══════════════ imports of blocks ═══════════════
from equations_D import (hh_extended_D, smart_steady_D, market_clearing_D,
    steady_auxilliary_D, banker_div_D, sdf_D, sdf_ss_D, sdf_banker_ss_D,
    government_ss_D, labor_ss_D, government_default_D, bond_price_ss_D,
    bond_return_D, ces_price_D, import_demand_D, deposit_return_D,
    capital_adj_D, labor_D, labor_market_D, labor_demand_D,
    intermediation_IC_D, bank_return_D, intermediation_P1_D, k_balance_sheet_D,
    cap_adj_cost_inter_D, macro_pru_tax_D, intermediation_P2_D, banker_div_res_D,
    intermediation_P3_D, divert_bond_foc_D, tax_rule_D, capital_producer_profit_D,
    budget_residual_D, sdf_banker_D, ghh_composite_D, welfare_agg_D)
from equations_F import (hh_extended_F, smart_steady_F, market_clearing_F,
    steady_auxilliary_F, banker_div_F, sdf_F, sdf_ss_F, sdf_banker_ss_F,
    government_ss_F, labor_ss_F, government_default_F, bond_price_ss_F,
    bond_return_F, ces_price_F, import_demand_F, deposit_return_F,
    capital_adj_F, labor_F, labor_market_F, labor_demand_F,
    intermediation_IC_F, bank_return_F, intermediation_P1_F, k_balance_sheet_F,
    cap_adj_cost_inter_F, macro_pru_tax_F, intermediation_P2_F, banker_div_res_F,
    intermediation_P3_F, divert_bond_foc_F, tax_rule_F, capital_producer_profit_F,
    budget_residual_F, sdf_banker_F, ghh_composite_F, welfare_agg_F)
from equations_global import (trade_balance, domestic_bond_clearing,
    portfolio_level_anchors, divert_portfolio_adj, bond_yield,
    global_goods_mkt, external_account_D)

# ═══════════════ CELL 10: SS model + solve ═══════════════
ha = sj.create_model([
    sdf_ss_D, sdf_banker_ss_D, government_default_D, bond_price_ss_D, bond_return_D,
    sdf_ss_F, sdf_banker_ss_F, government_default_F, bond_price_ss_F, bond_return_F,
    hh_extended_D, smart_steady_D, market_clearing_D, steady_auxilliary_D,
    banker_div_D, government_ss_D, labor_ss_D,
    hh_extended_F, smart_steady_F, market_clearing_F, steady_auxilliary_F,
    banker_div_F, government_ss_F, labor_ss_F,
    ces_price_D, import_demand_D, ces_price_F, import_demand_F,
    deposit_return_D, deposit_return_F,
    bond_yield, trade_balance, external_account_D, global_goods_mkt,
], name='MU HA Model 2 Country')

unknowns_ss = {'beta_D': 0.9850, 'beta_F': 0.9850, 'p': 0.99}
targets_ss  = ['deposit_mkt_D', 'deposit_mkt_F', 'ca_res_D']
print("Solving SS (cold start)...")
ss = ha.solve_steady_state(calibration_start, unknowns_ss, targets_ss, solver='broyden_custom')

def _apply_ss_anchors(ss_in, cal):
    anchors = {
        'phi_bD_D_ss': float(ss_in['q_b_D'])*float(ss_in['b_D_D'])/float(ss_in['n_inter_D']),
        'phi_bF_F_ss': float(ss_in['q_b_F'])*float(ss_in['b_F_F'])/(float(ss_in['p'])*float(ss_in['n_inter_F'])),
        'b_F_D_anchor': float(ss_in['b_F_D']), 'b_D_F_anchor': float(ss_in['b_D_F']),
        'excess_return_bD_D_ss': float(ss_in['rb_actual_D'])-float(ss_in['rdep_D'])-cal['T0_D'],
        'excess_return_bF_F_ss': float(ss_in['rb_actual_F'])-float(ss_in['rdep_F'])-cal['T0_F'],
        'excess_return_F_D_ss': float(ss_in['rb_actual_F'])-float(ss_in['rdep_D'])-cal['T0_D'],
        'excess_return_D_F_ss': float(ss_in['rb_actual_D'])-float(ss_in['rdep_F'])-cal['T0_F'],
        'psi_spread_F': float(ss_in['lambda_gk_F'])*cal['psi_lambda_B_F']/(float(ss_in['beta_inter_F'])*float(ss_in['Omega_F'])),
        'psi_spread_D': float(ss_in['lambda_gk_D'])*cal['psi_lambda_B_D']/(float(ss_in['beta_inter_D'])*float(ss_in['Omega_D'])),
        'q_b_D': float(ss_in['q_b_D']), 'q_b_F': float(ss_in['q_b_F']),
        'p': float(ss_in['p']), 'C_D_ss': float(ss_in['C_D']), 'C_F_ss': float(ss_in['C_F']),
    }
    cal.update(anchors)
    for k, v in anchors.items():
        ss_in.toplevel[k] = v
    ss_in.toplevel['b_F_D_ss'] = float(ss_in['b_F_D'])
    ss_in.toplevel['b_D_F_ss'] = float(ss_in['b_D_F'])
    ss_in.toplevel['Rgross_D'] = float(1+ss_in['rdep_D'])
    ss_in.toplevel['Rgross_F'] = float(1+ss_in['rdep_F'])
    _fr_D = float(ss_in['frisch_D']); _fr_F = float(ss_in['frisch_F'])
    ss_in.toplevel['X_D'] = (float(ss_in['C_D']) - float(ss_in['vphi_D'])*float(ss_in['N_D'])**(1+1/_fr_D)/(1+1/_fr_D))
    ss_in.toplevel['X_F'] = (float(ss_in['C_F']) - float(ss_in['vphi_F'])*float(ss_in['N_F'])**(1+1/_fr_F)/(1+1/_fr_F))
    ss_in.toplevel['U_D'] = ss_in.toplevel['X_D']/float(ss_in['C_D'])
    ss_in.toplevel['U_F'] = ss_in.toplevel['X_F']/float(ss_in['C_F'])
    ss_in.toplevel['Phi_D'] = float(ss_in['Phi_D'])
    ss_in.toplevel['Phi_F'] = float(ss_in['Phi_F'])
    ss_in.toplevel['value_D'] = (float(ss_in['beta_inter_D'])*float(ss_in['Omega_D'])*(1+float(ss_in['rn_D'])))
    ss_in.toplevel['value_F'] = (float(ss_in['beta_inter_F'])*float(ss_in['Omega_F'])*(1+float(ss_in['rn_F'])))
    for k, v in {'tau_mp_D':0.0,'tau_mp_F':0.0,'T_D':0.0,'T_F':0.0,'T_ls_D':0.0,'T_ls_F':0.0,
                 'b_F_D_res':0.0,'b_D_F_res':0.0,'rb_D_res':0.0,'rb_F_res':0.0,
                 'labor_mkt_res_D':0.0,'labor_mkt_res_F':0.0,'w_res_D':0.0,'w_res_F':0.0}.items():
        ss_in.toplevel[k] = v
    return anchors

# replicate cell 10 anchor application (psi_spread not in cell-10 anchors, but
# _apply_ss_anchors is re-run in cells 12/13 which adds it; final state identical)
_apply_ss_anchors(ss, calibration_start)

# ═══════════════ CELL 12: portfolio share targeting + Delta calibration ═══════════════
target_phi_bD_D, target_phi_bF_D = 0.25, 0.15
target_phi_bD_F, target_phi_bF_F = 0.15, 0.25
n_D = float(ss['n_inter_D']); n_F = float(ss['n_inter_F'])*float(ss['p'])
q_D = float(ss['q_b_D']); q_F = float(ss['q_b_F']); p_ss = float(ss['p'])
b_D_D_new = target_phi_bD_D*n_D/q_D
b_F_D_new = target_phi_bF_D*n_D/q_F
b_D_F_new = target_phi_bD_F*n_F/q_D
b_F_F_new = target_phi_bF_F*n_F/q_F
B_D_new = b_D_D_new + b_D_F_new
B_F_new = b_F_D_new + b_F_F_new
calibration_start.update({
    'b_D_D': b_D_D_new, 'b_F_D': b_F_D_new, 'b_D_F': b_D_F_new, 'b_F_F': b_F_F_new,
    'b_F_D_anchor': b_F_D_new, 'b_D_F_anchor': b_D_F_new,
    'phi_bF_D_ss': target_phi_bF_D,
    'B_supply_D': B_D_new, 'b_gov_D': B_D_new, 'b_gov_ss_D': B_D_new,
    'B_supply_F': B_F_new, 'b_gov_F': B_F_new, 'b_gov_ss_F': B_F_new,
})
print("Re-solving SS with new portfolio allocation...")
_unknowns_warm = {'beta_D': float(ss['beta_D']), 'beta_F': float(ss['beta_F']), 'p': float(ss['p'])}
ss = ha.solve_steady_state(calibration_start, _unknowns_warm, targets_ss, solver='broyden_custom')
_apply_ss_anchors(ss, calibration_start)

ratio_D = ratio_F = 2.0
def _ic_delta(phi_own, phi_cross, nu_K, nu_b_own, nu_b_cross, eta, lam, theta, ratio):
    kappa = theta - phi_own - phi_cross
    value = nu_K*kappa + nu_b_own*phi_own + nu_b_cross*phi_cross + eta
    denom = phi_own + ratio*phi_cross
    delta_own = (phi_own + phi_cross - (theta - value/lam))/denom
    return float(delta_own), float(ratio*delta_own), float(value)

phi_bD_D_ss = float(ss['q_b_D'])*float(ss['b_D_D'])/float(ss['n_inter_D'])
phi_bF_D_ss = float(ss['q_b_F'])*float(ss['b_F_D'])/float(ss['n_inter_D'])
D_bD_D, D_bF_D, val_D = _ic_delta(phi_bD_D_ss, phi_bF_D_ss,
    float(ss['nu_K_D']), float(ss['nu_bD_D']), float(ss['nu_bF_D']), float(ss['eta_D']),
    float(ss['lambda_gk_D']), float(ss['theta_D']), ratio_D)
n_F_ss = float(ss['n_inter_F'])*float(ss['p'])
phi_bF_F_ss = float(ss['q_b_F'])*float(ss['b_F_F'])/n_F_ss
phi_bD_F_ss = float(ss['q_b_D'])*float(ss['b_D_F'])/n_F_ss
D_bF_F, D_bD_F, val_F = _ic_delta(phi_bF_F_ss, phi_bD_F_ss,
    float(ss['nu_K_F']), float(ss['nu_bF_F']), float(ss['nu_bD_F']), float(ss['eta_F']),
    float(ss['lambda_gk_F']), float(ss['theta_F']), ratio_F)
calibration_start.update({'Delta_bD_D': D_bD_D, 'Delta_bF_D': D_bF_D,
                          'Delta_bF_F': D_bF_F, 'Delta_bD_F': D_bD_F})
jlog('Delta_calibrated', {'Delta_bD_D': D_bD_D, 'Delta_bF_D': D_bF_D,
                          'Delta_bF_F': D_bF_F, 'Delta_bD_F': D_bD_F})

# ═══════════════ CELL 13: depreciation recalibration ═══════════════
rk_D_target = rk_F_target = 0.01
alpha_D_cal = calibration_start['alpha_D']; alpha_F_cal = calibration_start['alpha_F']
K_D_cur = float(ss['K_D']); K_F_cur = float(ss['K_F'])
Y_D_cur = float(ss['Y_D']); Y_F_cur = float(ss['Y_F'])
delta_D_cal = alpha_D_cal*Y_D_cur/K_D_cur - rk_D_target
delta_F_cal = alpha_F_cal*Y_F_cur/K_F_cur - rk_F_target
calibration_start.update({'delta_D': delta_D_cal, 'delta_F': delta_F_cal})
jlog('delta_recalibrated', {'delta_D': delta_D_cal, 'delta_F': delta_F_cal})
print("Final SS re-solve with calibrated delta...")
ss = ha.solve_steady_state(calibration_start, unknowns_ss, targets_ss, solver='broyden_custom')
_apply_ss_anchors(ss, calibration_start)
cali_D = cali_F = ss
ss_final = copy.deepcopy(ss)

# ═══════════════ A. SS residual dump ═══════════════
ss_keys = ['beta_D','beta_F','p','q_b_D','q_b_F','rb_D','rb_F','rb_actual_D','rb_actual_F',
           'rdep_D','rdep_F','rk_D','rk_F','K_D','K_F','Y_D','Y_F','C_D','C_F','I_D','I_F',
           'G_D','G_F','N_D','N_F','w_D','w_F','NX_D','NX_F','IM_D','IM_F','nfa_D',
           'n_inter_D','n_inter_F','theta_D','theta_F','b_gov_D','b_gov_F',
           'b_D_D','b_F_D','b_D_F','b_F_F','div_D','div_F','DEP_D','DEP_F',
           'D_supply_D','D_supply_F','TAX_D','TAX_F','P_CES_D','P_CES_F',
           'lambda_gk_D','lambda_gk_F','Omega_D','Omega_F','m_D','m_F','vphi_D','vphi_F',
           'goods_mkt_D','goods_mkt_F','global_goods_res','ca_res_D',
           'deposit_mkt_D','deposit_mkt_F','X_D','X_F','spread_rb',
           'nu_K_D','nu_bD_D','nu_bF_D','eta_D','nu_K_F','nu_bD_F','nu_bF_F','eta_F',
           'Z_D','Z_F','def_rate_D','def_rate_F','Phi_D','Phi_F','cap_profit_D','cap_profit_F']
RES['ss'] = {k: float(ss[k]) for k in ss_keys}
print(json.dumps({k: RES['ss'][k] for k in ['goods_mkt_D','goods_mkt_F','ca_res_D','global_goods_res','G_D','G_F','nfa_D','p']}, indent=2))

# ═══════════════ B. IC/Bellman residual at linearization point ═══════════════
def ic_resid(c):
    pdiv = float(ss['p']) if c=='F' else 1.0
    lam = float(ss[f'lambda_gk_{c}']); theta_c = float(ss[f'theta_{c}'])
    Q_c = float(ss[f'Q_{c}']); K_c = float(ss[f'K_{c}']); n_c = float(ss[f'n_inter_{c}'])
    kappa_c = Q_c*K_c/n_c; eta_c = float(ss[f'eta_{c}'])
    if c=='D':
        nu_h, nu_x = float(ss['nu_bD_D']), float(ss['nu_bF_D'])
        q_h, q_x = float(ss['q_b_D']), float(ss['q_b_F'])
        b_h, b_x = float(ss['b_D_D']), float(ss['b_F_D'])
        Dh, Dx = calibration_start['Delta_bD_D'], calibration_start['Delta_bF_D']
    else:
        nu_h, nu_x = float(ss['nu_bF_F']), float(ss['nu_bD_F'])
        q_h, q_x = float(ss['q_b_F']), float(ss['q_b_D'])
        b_h, b_x = float(ss['b_F_F']), float(ss['b_D_F'])
        Dh, Dx = calibration_start['Delta_bF_F'], calibration_start['Delta_bD_F']
    nu_K = float(ss[f'nu_K_{c}'])
    phi_h = q_h*b_h/(pdiv*n_c); phi_x = q_x*b_x/(pdiv*n_c)
    value_c = nu_K*kappa_c + nu_h*phi_h + nu_x*phi_x + eta_c
    theta_tgt = value_c/lam + (1-Dh)*phi_h + (1-Dx)*phi_x
    return theta_c - theta_tgt
jlog('ic_res_at_ss_final', {'D': ic_resid('D'), 'F': ic_resid('F')})

# k_balance_sheet residual at ss_final
K_res_D_ss = (float(ss['Q_D'])*float(ss['K_D']) + float(ss['q_b_D'])*float(ss['b_D_D'])
              + float(ss['q_b_F'])*float(ss['b_F_D']) - float(ss['theta_D'])*float(ss['n_inter_D']))
K_res_F_ss = (float(ss['Q_F'])*float(ss['K_F'])
              + (float(ss['q_b_F'])*float(ss['b_F_F'])+float(ss['q_b_D'])*float(ss['b_D_F']))/float(ss['p'])
              - float(ss['theta_F'])*float(ss['n_inter_F']))
jlog('K_res_at_ss_final', {'D': K_res_D_ss, 'F': K_res_F_ss})
# P2 residual at SS (should be 0; nonzero if m bookkeeping inconsistent)
P2_D_ss = (1-float(ss['f_D']))*((1+float(ss['rn_D']))*float(ss['n_inter_D'])+float(ss['cap_profit_D'])) + float(ss['m_D']) - float(ss['n_inter_D'])
P2_F_ss = (1-float(ss['f_F']))*((1+float(ss['rn_F']))*float(ss['n_inter_F'])+float(ss['cap_profit_F'])) + float(ss['m_F']) - float(ss['n_inter_F'])
jlog('P2_res_at_ss_final', {'D': P2_D_ss, 'F': P2_F_ss})

# ═══════════════ CELL 20: full model ═══════════════
sys.setrecursionlimit(5000)
financial_solved_D = combine([intermediation_P1_D, intermediation_IC_D]).solved(
    unknowns={'nu_K_D': float(cali_D['nu_K_D']), 'nu_bD_D': float(cali_D['nu_bD_D']),
              'nu_bF_D': float(cali_D['nu_bF_D']), 'eta_D': float(cali_D['eta_D']),
              'theta_D': float(cali_D['theta_D'])},
    targets=['nu_K_res_D','nu_bD_res_D','nu_bF_res_D','eta_res_D','ic_res_D'],
    solver='broyden_custom')
financial_solved_F = combine([intermediation_P1_F, intermediation_IC_F]).solved(
    unknowns={'nu_K_F': float(cali_F['nu_K_F']), 'nu_bF_F': float(cali_F['nu_bF_F']),
              'nu_bD_F': float(cali_F['nu_bD_F']), 'eta_F': float(cali_F['eta_F']),
              'theta_F': float(cali_F['theta_F'])},
    targets=['nu_K_res_F','nu_bF_res_F','nu_bD_res_F','eta_res_F','ic_res_F'],
    solver='broyden_custom')

ha_full = sj.create_model([
    deposit_return_D, tax_rule_D, hh_extended_D, ghh_composite_D, sdf_D, sdf_banker_D,
    government_default_D, financial_solved_D, bond_return_D, bank_return_D,
    cap_adj_cost_inter_D, macro_pru_tax_D, intermediation_P2_D, intermediation_P3_D,
    k_balance_sheet_D, capital_adj_D, capital_producer_profit_D, budget_residual_D,
    labor_D, labor_market_D, labor_demand_D, banker_div_res_D, market_clearing_D, welfare_agg_D,
    deposit_return_F, tax_rule_F, hh_extended_F, ghh_composite_F, sdf_F, sdf_banker_F,
    government_default_F, financial_solved_F, bond_return_F, bank_return_F,
    cap_adj_cost_inter_F, macro_pru_tax_F, intermediation_P2_F, intermediation_P3_F,
    k_balance_sheet_F, capital_adj_F, capital_producer_profit_F, budget_residual_F,
    labor_F, labor_market_F, labor_demand_F, banker_div_res_F, market_clearing_F, welfare_agg_F,
    ces_price_D, import_demand_D, ces_price_F, import_demand_F,
    trade_balance, external_account_D, domestic_bond_clearing, bond_yield,
    portfolio_level_anchors, divert_portfolio_adj, divert_bond_foc_D, divert_bond_foc_F,
    global_goods_mkt,
], name="Full 2-Country MU HANK")

unknowns_tp = ['K_D','n_inter_D','div_D','I_D','Q_D','b_gov_D','N_D','b_F_D','w_D','rdep_D',
               'K_F','n_inter_F','div_F','I_F','Q_F','b_gov_F','N_F','b_D_F','w_F','rdep_F',
               'p','q_b_D','q_b_F']
targets_tp = ['deposit_mkt_D','K_res_D','n_inter_val_D','div_res_D','capital_res_D','q_res_D',
              'b_gov_res_D','b_F_D_res','labor_mkt_res_D','w_res_D',
              'deposit_mkt_F','K_res_F','n_inter_val_F','div_res_F','capital_res_F','q_res_F',
              'b_gov_res_F','b_D_F_res','labor_mkt_res_F','w_res_F',
              'goods_mkt_D','rb_D_res','rb_F_res']
T = 500
exogenous = ['Z_D','shock_def_D','Z_F','shock_def_F']

print("Computing full GE Jacobian G (T=500)...")
G = ha_full.solve_jacobian(ss_final, unknowns=unknowns_tp, targets=targets_tp, inputs=exogenous, T=T)
print("G done.")

rho = 0.8
dZ = 0.01*rho**np.arange(T)
zero = np.zeros(T)
irfs_Z_D   = G @ {'Z_D': dZ,  'Z_F': zero, 'shock_def_D': zero, 'shock_def_F': zero}
irfs_def_D = G @ {'Z_D': zero,'Z_F': zero, 'shock_def_D': dZ,  'shock_def_F': zero}

# ═══════════════ C. Walras leak test ═══════════════
leakvars = ['ca_res_D','goods_mkt_F','goods_mkt_D','global_goods_res','deposit_mkt_D','deposit_mkt_F']
for nm, irf in [('Z_D', irfs_Z_D), ('def_D', irfs_def_D)]:
    out = {}
    for v in leakvars:
        try:
            path = irf[v]
            out[v] = {'max_abs': float(np.max(np.abs(path))),
                      't0': float(path[0]), 't1': float(path[1]), 't5': float(path[5])}
        except (KeyError, TypeError):
            out[v] = 'MISSING from G outputs'
    jlog(f'leak_{nm}', out)

# ═══════════════ D. Leak decomposition ═══════════════
def dec(irf, label):
    # factor gap D: d(Y - wN - mpk*K(-1)) ; SS: mpk*K = alpha*Y, w*N=(1-alpha)Y
    Yss=float(ss['Y_D']); Kss=float(ss['K_D']); Nss=float(ss['N_D']); wss=float(ss['w_D']); mpkss=float(ss['mpk_D'])
    dY=irf['Y_D']; dN=irf['N_D']; dw=irf['w_D']; dK=irf['K_D']; dmpk=irf['mpk_D']
    dKlag = np.concatenate([[0.0], dK[:-1]])
    gapD = dY - wss*dN - Nss*dw - mpkss*dKlag - Kss*dmpk
    # same for F
    YssF=float(ss['Y_F']); KssF=float(ss['K_F']); NssF=float(ss['N_F']); wssF=float(ss['w_F']); mpkssF=float(ss['mpk_F'])
    dYF=irf['Y_F']; dNF=irf['N_F']; dwF=irf['w_F']; dKF=irf['K_F']; dmpkF=irf['mpk_F']
    dKlagF = np.concatenate([[0.0], dKF[:-1]])
    gapF = dYF - wssF*dNF - NssF*dwF - mpkssF*dKlagF - KssF*dmpkF
    # F-bank ToT revaluation gap (F-goods): true bond payoff converts at 1/p_t, code uses 1/p(-1)
    # payoff_ss (D-goods) = (1+rb)*[q_F*b_FF + q_D*b_DF]; gap = payoff_ss * d(1/p_t - 1/p_{t-1})
    pss=float(ss['p'])
    payoff = ((1+float(ss['rb_actual_F']))*float(ss['q_b_F'])*float(ss['b_F_F'])
              + (1+float(ss['rb_actual_D']))*float(ss['q_b_D'])*float(ss['b_D_F']))
    dp = irf['p']; dplag = np.concatenate([[0.0], dp[:-1]])
    rev_gapF = payoff * (-(dp - dplag)/pss**2)   # d(1/p) = -dp/p^2
    out = {
        'gapD_max': float(np.max(np.abs(gapD))),
        'gapF_max': float(np.max(np.abs(gapF))),
        'rev_gapF_max': float(np.max(np.abs(rev_gapF))),
        'ca_res_D_minus_gapD_max': float(np.max(np.abs(irf['ca_res_D'] - gapD))),
        'corr_ca_gapD': float(np.corrcoef(irf['ca_res_D'][:200], gapD[:200])[0,1]),
    }
    if True:
        gmF = irf['goods_mkt_F']
        out['goodsF_minus_gapF_minus_rev_max'] = float(np.max(np.abs(gmF - gapF - rev_gapF)))
        out['goodsF_minus_gapF_max'] = float(np.max(np.abs(gmF - gapF)))
        out['corr_goodsF_(gapF+rev)'] = float(np.corrcoef(gmF[:200], (gapF+rev_gapF)[:200])[0,1])
    jlog(f'decomp_{label}', out)

try:
    dec(irfs_Z_D, 'Z_D')
    dec(irfs_def_D, 'def_D')
except Exception as e:
    jlog('decomp_error', repr(e))

# key IRF stats for mechanism validation
def stats(irf, label):
    out = {}
    for v in ['Y_D','Y_F','C_D','C_F','n_inter_D','n_inter_F','q_b_D','q_b_F','rb_D','spread_rb',
              'def_rate_D','b_gov_D','p','NX_D','b_F_D','b_D_F','theta_D','rdep_D','K_D','I_D','U_D','U_F']:
        try:
            path = irf[v]
            am = int(np.argmax(np.abs(path[:100])))
            out[v] = {'peak_t': am, 'peak': float(path[am]), 't0': float(path[0])}
        except (KeyError, TypeError):
            out[v] = 'missing'
    jlog(f'irfstats_{label}', out)
stats(irfs_def_D, 'def_D')
stats(irfs_Z_D, 'Z_D')

np.savez(ROOT/'audit_artifacts'/'irfs.npz',
         **{f'Z_D__{k}': v for k, v in irfs_Z_D.toplevel.items() if np.ndim(v)==1},
         **{f'def_D__{k}': v for k, v in irfs_def_D.toplevel.items() if np.ndim(v)==1})

# ═══════════════ E. Jacobian conditioning of target system ═══════════════
try:
    Tj = 120
    Jfull = ha_full.jacobian(ss_final, inputs=unknowns_tp, outputs=targets_tp, T=Tj)
    nU = len(unknowns_tp); nT = len(targets_tp)
    H = np.zeros((nT*Tj, nU*Tj))
    for i, tg in enumerate(targets_tp):
        for j, un in enumerate(unknowns_tp):
            try:
                blk = Jfull[tg][un]
                H[i*Tj:(i+1)*Tj, j*Tj:(j+1)*Tj] = blk @ np.eye(Tj)
            except (KeyError, TypeError):
                pass
    sv = np.linalg.svd(H, compute_uv=False)
    jlog('H_U_singvals', {'max': float(sv[0]), 'min': float(sv[-1]),
                          'cond': float(sv[0]/sv[-1]),
                          'smallest_5': [float(x) for x in sv[-5:]]})
except Exception as e:
    jlog('H_U_error', repr(e))

with open(ROOT/'audit_artifacts'/'run1_results.json','w') as fh:
    json.dump(RES, fh, indent=2, default=str)
print("AUDIT RUN COMPLETE")
