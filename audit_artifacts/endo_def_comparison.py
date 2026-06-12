"""endo_def_comparison.py
Compare IRFs: AB-audit (def_rate = f(b/Y_ss) + eps) vs endo-def (def_rate = f(b/Y) + eps).

Both models share the same steady state (f(.)=0 at SS regardless of which Y is used).
Two Jacobian solves; IRF comparison printed to stdout and saved to endo_def_results.json.
"""
import json, sys, copy
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CODE = ROOT / 'code'
sys.path.insert(0, str(CODE))   # equations_*.py live in code/
sys.path.insert(0, str(ROOT))
import os; os.chdir(ROOT)

import sequence_jacobian as sj
from sequence_jacobian import simple, combine

# ── calibration (identical to run_audit.py) ──────────────────────────────────
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
    'def_rate_D': 0.000, 'def_rate_F': 0.0,
    'def_scale_D': 0.25, 'def_scale_F': 0.25,
    'def_curvature_D': 0.5, 'def_curvature_F': 0.5,
    'def_offset_D': 0.05, 'def_offset_F': 0.05,
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
    'nZ_D': 15, 'nZ_F': 15, 'nDep_D': 500, 'nDep_F': 500, 'Depmax_D': 150, 'Depmax_F': 150,
}
_n_D = calibration_start['n_inter_D']; _n_F = calibration_start['n_inter_F']
b_F_D = calibration_start['phi_bF_D_ss'] * _n_D / calibration_start['q_b_F']
b_D_F = calibration_start['phi_bD_F_ss'] * _n_F / calibration_start['q_b_D']
calibration_start.update({
    'b_F_D': b_F_D, 'b_D_F': b_D_F,
    'b_D_D': calibration_start['B_supply_D'] - b_D_F,
    'b_F_F': calibration_start['B_supply_F'] - b_F_D,
    'b_F_D_anchor': b_F_D, 'b_D_F_anchor': b_D_F,
    'psi_bD_D': 0.0, 'psi_bF_F': 0.0,
})

# ── block imports ─────────────────────────────────────────────────────────────
from equations_D import (
    hh_extended_D, smart_steady_D, market_clearing_D, steady_auxilliary_D,
    banker_div_D, sdf_D, sdf_ss_D, sdf_banker_ss_D, government_ss_D, labor_ss_D,
    government_default_D,   # NEW: uses Y_D
    bond_price_ss_D, bond_return_D, ces_price_D, import_demand_D, deposit_return_D,
    capital_adj_D, labor_D, labor_market_D, labor_demand_D,
    intermediation_IC_D, bank_return_D, intermediation_P1_D, k_balance_sheet_D,
    cap_adj_cost_inter_D, macro_pru_tax_D, intermediation_P2_D, banker_div_res_D,
    intermediation_P3_D, divert_bond_foc_D, tax_rule_D, capital_producer_profit_D,
    budget_residual_D, sdf_banker_D, ghh_composite_D, welfare_agg_D)
from equations_F import (
    hh_extended_F, smart_steady_F, market_clearing_F, steady_auxilliary_F,
    banker_div_F, sdf_F, sdf_ss_F, sdf_banker_ss_F, government_ss_F, labor_ss_F,
    government_default_F,   # NEW: uses Y_F
    bond_price_ss_F, bond_return_F, ces_price_F, import_demand_F, deposit_return_F,
    capital_adj_F, labor_F, labor_market_F, labor_demand_F,
    intermediation_IC_F, bank_return_F, intermediation_P1_F, k_balance_sheet_F,
    cap_adj_cost_inter_F, macro_pru_tax_F, intermediation_P2_F, banker_div_res_F,
    intermediation_P3_F, divert_bond_foc_F, tax_rule_F, capital_producer_profit_F,
    budget_residual_F, sdf_banker_F, ghh_composite_F, welfare_agg_F)
from equations_global import (
    trade_balance, domestic_bond_clearing, portfolio_level_anchors,
    divert_portfolio_adj, bond_yield, global_goods_mkt, external_account_D)

# ── OLD-style blocks: use Y_ss (constant) in denominator ─────────────────────
@simple
def government_default_D_old(shock_def_D, b_gov_D, Y_ss_D, b_gov_ss_D,
                              def_scale_D, def_curvature_D, def_offset_D):
    debt_ratio_D = b_gov_D(-1) / Y_ss_D
    ss_ratio_D   = b_gov_ss_D  / Y_ss_D
    def_rate_D   = shock_def_D + def_scale_D * (
        (debt_ratio_D + def_offset_D) ** def_curvature_D
      - (ss_ratio_D  + def_offset_D) ** def_curvature_D
    )
    return def_rate_D

@simple
def government_default_F_old(shock_def_F, b_gov_F, Y_ss_F, b_gov_ss_F,
                              def_scale_F, def_curvature_F, def_offset_F):
    debt_ratio_F = b_gov_F(-1) / Y_ss_F
    ss_ratio_F   = b_gov_ss_F  / Y_ss_F
    def_rate_F   = shock_def_F + def_scale_F * (
        (debt_ratio_F + def_offset_F) ** def_curvature_F
      - (ss_ratio_F  + def_offset_F) ** def_curvature_F
    )
    return def_rate_F

# ── helpers ───────────────────────────────────────────────────────────────────
def _apply_anchors(ss_in, cal):
    anchors = {
        'phi_bD_D_ss': float(ss_in['q_b_D'])*float(ss_in['b_D_D'])/float(ss_in['n_inter_D']),
        'phi_bF_F_ss': float(ss_in['q_b_F'])*float(ss_in['b_F_F'])/(float(ss_in['p'])*float(ss_in['n_inter_F'])),
        'b_F_D_anchor': float(ss_in['b_F_D']), 'b_D_F_anchor': float(ss_in['b_D_F']),
        'excess_return_bD_D_ss': float(ss_in['rb_actual_D'])-float(ss_in['rdep_D'])-cal['T0_D'],
        'excess_return_bF_F_ss': float(ss_in['rb_actual_F'])-float(ss_in['rdep_F'])-cal['T0_F'],
        'excess_return_F_D_ss':  float(ss_in['rb_actual_F'])-float(ss_in['rdep_D'])-cal['T0_D'],
        'excess_return_D_F_ss':  float(ss_in['rb_actual_D'])-float(ss_in['rdep_F'])-cal['T0_F'],
        'psi_spread_D': float(ss_in['lambda_gk_D'])*cal['psi_lambda_B_D']/(float(ss_in['beta_inter_D'])*float(ss_in['Omega_D'])),
        'psi_spread_F': float(ss_in['lambda_gk_F'])*cal['psi_lambda_B_F']/(float(ss_in['beta_inter_F'])*float(ss_in['Omega_F'])),
        'q_b_D': float(ss_in['q_b_D']), 'q_b_F': float(ss_in['q_b_F']),
        'p': float(ss_in['p']), 'C_D_ss': float(ss_in['C_D']), 'C_F_ss': float(ss_in['C_F']),
    }
    cal.update(anchors)
    for k, v in anchors.items():
        ss_in.toplevel[k] = v
    ss_in.toplevel.update({'b_F_D_ss': float(ss_in['b_F_D']), 'b_D_F_ss': float(ss_in['b_D_F']),
                           'Rgross_D': float(1+ss_in['rdep_D']), 'Rgross_F': float(1+ss_in['rdep_F'])})
    _fr_D = float(ss_in['frisch_D']); _fr_F = float(ss_in['frisch_F'])
    ss_in.toplevel['X_D'] = float(ss_in['C_D']) - float(ss_in['vphi_D'])*float(ss_in['N_D'])**(1+1/_fr_D)/(1+1/_fr_D)
    ss_in.toplevel['X_F'] = float(ss_in['C_F']) - float(ss_in['vphi_F'])*float(ss_in['N_F'])**(1+1/_fr_F)/(1+1/_fr_F)
    ss_in.toplevel['U_D'] = ss_in.toplevel['X_D']/float(ss_in['C_D'])
    ss_in.toplevel['U_F'] = ss_in.toplevel['X_F']/float(ss_in['C_F'])
    ss_in.toplevel['Phi_D'] = float(ss_in['Phi_D']); ss_in.toplevel['Phi_F'] = float(ss_in['Phi_F'])
    ss_in.toplevel['value_D'] = float(ss_in['beta_inter_D'])*float(ss_in['Omega_D'])*(1+float(ss_in['rn_D']))
    ss_in.toplevel['value_F'] = float(ss_in['beta_inter_F'])*float(ss_in['Omega_F'])*(1+float(ss_in['rn_F']))
    for k in ('tau_mp_D','tau_mp_F','T_D','T_F','T_ls_D','T_ls_F',
              'b_F_D_res','b_D_F_res','rb_D_res','rb_F_res',
              'labor_mkt_res_D','labor_mkt_res_F','w_res_D','w_res_F'):
        ss_in.toplevel[k] = 0.0

# ── SS solve (shared) ─────────────────────────────────────────────────────────
# Use the NEW government_default for SS — identical at SS, but required by model graph
ha_ss = sj.create_model([
    sdf_ss_D, sdf_banker_ss_D, government_default_D, bond_price_ss_D, bond_return_D,
    sdf_ss_F, sdf_banker_ss_F, government_default_F, bond_price_ss_F, bond_return_F,
    hh_extended_D, smart_steady_D, market_clearing_D, steady_auxilliary_D,
    banker_div_D, government_ss_D, labor_ss_D,
    hh_extended_F, smart_steady_F, market_clearing_F, steady_auxilliary_F,
    banker_div_F, government_ss_F, labor_ss_F,
    ces_price_D, import_demand_D, ces_price_F, import_demand_F,
    deposit_return_D, deposit_return_F, bond_yield,
    trade_balance, external_account_D, global_goods_mkt,
], name='SS model')

unknowns_ss = {'beta_D': 0.9850, 'beta_F': 0.9850, 'p': 0.99}
targets_ss  = ['deposit_mkt_D', 'deposit_mkt_F', 'ca_res_D']
print("Solving SS (cold start)...")
ss = ha_ss.solve_steady_state(calibration_start, unknowns_ss, targets_ss, solver='broyden_custom')
_apply_anchors(ss, calibration_start)

# portfolio + delta recalibration (mirror cells 12-13)
target_phi_bD_D, target_phi_bF_D = 0.25, 0.15
target_phi_bD_F, target_phi_bF_F = 0.15, 0.25
n_D = float(ss['n_inter_D']); n_F = float(ss['n_inter_F'])*float(ss['p'])
q_D = float(ss['q_b_D']); q_F = float(ss['q_b_F'])
b_D_D_new = target_phi_bD_D*n_D/q_D; b_F_D_new = target_phi_bF_D*n_D/q_F
b_D_F_new = target_phi_bD_F*n_F/q_D; b_F_F_new = target_phi_bF_F*n_F/q_F
calibration_start.update({
    'b_D_D': b_D_D_new, 'b_F_D': b_F_D_new, 'b_D_F': b_D_F_new, 'b_F_F': b_F_F_new,
    'b_F_D_anchor': b_F_D_new, 'b_D_F_anchor': b_D_F_new, 'phi_bF_D_ss': target_phi_bF_D,
    'B_supply_D': b_D_D_new+b_D_F_new, 'b_gov_D': b_D_D_new+b_D_F_new, 'b_gov_ss_D': b_D_D_new+b_D_F_new,
    'B_supply_F': b_F_D_new+b_F_F_new, 'b_gov_F': b_F_D_new+b_F_F_new, 'b_gov_ss_F': b_F_D_new+b_F_F_new,
})
print("Re-solving SS with portfolio targets...")
ss = ha_ss.solve_steady_state(calibration_start,
     {'beta_D': float(ss['beta_D']), 'beta_F': float(ss['beta_F']), 'p': float(ss['p'])},
     targets_ss, solver='broyden_custom')
_apply_anchors(ss, calibration_start)

# delta recalibration
rk_target = 0.01
alpha_D = calibration_start['alpha_D']; alpha_F = calibration_start['alpha_F']
calibration_start.update({
    'delta_D': alpha_D*float(ss['Y_D'])/float(ss['K_D']) - rk_target,
    'delta_F': alpha_F*float(ss['Y_F'])/float(ss['K_F']) - rk_target,
})

# IC delta calibration
ratio = 2.0
def _ic_delta(phi_own, phi_cross, nu_K, nu_b_own, nu_b_cross, eta, lam, theta, ratio):
    kappa = theta - phi_own - phi_cross
    value = nu_K*kappa + nu_b_own*phi_own + nu_b_cross*phi_cross + eta
    denom = phi_own + ratio*phi_cross
    delta_own = (phi_own + phi_cross - (theta - value/lam)) / denom
    return float(delta_own), float(ratio*delta_own)

phi_bD_D = float(ss['q_b_D'])*float(ss['b_D_D'])/float(ss['n_inter_D'])
phi_bF_D = float(ss['q_b_F'])*float(ss['b_F_D'])/float(ss['n_inter_D'])
D_bD_D, D_bF_D = _ic_delta(phi_bD_D, phi_bF_D, float(ss['nu_K_D']), float(ss['nu_bD_D']),
                            float(ss['nu_bF_D']), float(ss['eta_D']), float(ss['lambda_gk_D']),
                            float(ss['theta_D']), ratio)
n_F_p = float(ss['n_inter_F'])*float(ss['p'])
phi_bF_F = float(ss['q_b_F'])*float(ss['b_F_F'])/n_F_p
phi_bD_F = float(ss['q_b_D'])*float(ss['b_D_F'])/n_F_p
D_bF_F, D_bD_F = _ic_delta(phi_bF_F, phi_bD_F, float(ss['nu_K_F']), float(ss['nu_bF_F']),
                            float(ss['nu_bD_F']), float(ss['eta_F']), float(ss['lambda_gk_F']),
                            float(ss['theta_F']), ratio)
calibration_start.update({'Delta_bD_D': D_bD_D, 'Delta_bF_D': D_bF_D,
                          'Delta_bF_F': D_bF_F, 'Delta_bD_F': D_bD_F})

print("Final SS solve...")
ss = ha_ss.solve_steady_state(calibration_start, unknowns_ss, targets_ss, solver='broyden_custom')
_apply_anchors(ss, calibration_start)
ss_final = copy.deepcopy(ss)
cali = ss_final

print(f"SS solved. Y_D={float(ss['Y_D']):.4f}, def_rate_D={float(ss['def_rate_D']):.6f}, "
      f"spread_rb={float(ss['spread_rb'])*400:.1f}bps")

# ── build full transition models ──────────────────────────────────────────────
sys.setrecursionlimit(5000)
T = 100

def _fin_solved(country, cali_c):
    c = country
    if c == 'D':
        return combine([intermediation_P1_D, intermediation_IC_D]).solved(
            unknowns={'nu_K_D': float(cali_c['nu_K_D']), 'nu_bD_D': float(cali_c['nu_bD_D']),
                      'nu_bF_D': float(cali_c['nu_bF_D']), 'eta_D': float(cali_c['eta_D']),
                      'theta_D': float(cali_c['theta_D'])},
            targets=['nu_K_res_D','nu_bD_res_D','nu_bF_res_D','eta_res_D','ic_res_D'],
            solver='broyden_custom')
    else:
        return combine([intermediation_P1_F, intermediation_IC_F]).solved(
            unknowns={'nu_K_F': float(cali_c['nu_K_F']), 'nu_bF_F': float(cali_c['nu_bF_F']),
                      'nu_bD_F': float(cali_c['nu_bD_F']), 'eta_F': float(cali_c['eta_F']),
                      'theta_F': float(cali_c['theta_F'])},
            targets=['nu_K_res_F','nu_bF_res_F','nu_bD_res_F','eta_res_F','ic_res_F'],
            solver='broyden_custom')

fin_D = _fin_solved('D', cali); fin_F = _fin_solved('F', cali)

COMMON_BLOCKS_D = [
    deposit_return_D, tax_rule_D, hh_extended_D, ghh_composite_D, sdf_D, sdf_banker_D,
    fin_D, bond_return_D, bank_return_D,
    cap_adj_cost_inter_D, macro_pru_tax_D, intermediation_P2_D, intermediation_P3_D,
    k_balance_sheet_D, capital_adj_D, capital_producer_profit_D, budget_residual_D,
    labor_D, labor_market_D, labor_demand_D, banker_div_res_D, market_clearing_D, welfare_agg_D,
]
COMMON_BLOCKS_F = [
    deposit_return_F, tax_rule_F, hh_extended_F, ghh_composite_F, sdf_F, sdf_banker_F,
    fin_F, bond_return_F, bank_return_F,
    cap_adj_cost_inter_F, macro_pru_tax_F, intermediation_P2_F, intermediation_P3_F,
    k_balance_sheet_F, capital_adj_F, capital_producer_profit_F, budget_residual_F,
    labor_F, labor_market_F, labor_demand_F, banker_div_res_F, market_clearing_F, welfare_agg_F,
]
GLOBAL_BLOCKS = [
    ces_price_D, import_demand_D, ces_price_F, import_demand_F,
    trade_balance, external_account_D, domestic_bond_clearing, bond_yield,
    portfolio_level_anchors, divert_portfolio_adj, divert_bond_foc_D, divert_bond_foc_F,
    global_goods_mkt,
]

# NEW model (Y_D denominator)
ha_new = sj.create_model(
    [government_default_D] + COMMON_BLOCKS_D +
    [government_default_F] + COMMON_BLOCKS_F +
    GLOBAL_BLOCKS, name='endo-def (Y_D)')

# OLD model (Y_ss denominator)
ha_old = sj.create_model(
    [government_default_D_old] + COMMON_BLOCKS_D +
    [government_default_F_old] + COMMON_BLOCKS_F +
    GLOBAL_BLOCKS, name='AB-audit (Y_ss)')

unknowns_tp = [
    'K_D', 'n_inter_D', 'div_D', 'I_D', 'Q_D', 'b_gov_D', 'N_D', 'b_F_D', 'w_D', 'rdep_D',
    'K_F', 'n_inter_F', 'div_F', 'I_F', 'Q_F', 'b_gov_F', 'N_F', 'b_D_F', 'w_F', 'rdep_F',
    'p', 'q_b_D', 'q_b_F',
]
targets_tp = [
    'deposit_mkt_D', 'K_res_D', 'n_inter_val_D', 'div_res_D',
    'capital_res_D', 'q_res_D', 'b_gov_res_D', 'b_F_D_res',
    'labor_mkt_res_D', 'w_res_D',
    'deposit_mkt_F', 'K_res_F', 'n_inter_val_F', 'div_res_F',
    'capital_res_F', 'q_res_F', 'b_gov_res_F', 'b_D_F_res',
    'labor_mkt_res_F', 'w_res_F',
    'goods_mkt_D',
    'rb_D_res', 'rb_F_res',
]
exogenous = ['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F']

print("\nComputing Jacobian — NEW model (endo-def, Y_D denominator)...")
G_new = ha_new.solve_jacobian(ss_final, unknowns=unknowns_tp, targets=targets_tp,
                               inputs=exogenous, T=T)
print("Jacobian NEW done.")

print("\nComputing Jacobian — OLD model (AB-audit, Y_ss denominator)...")
G_old = ha_old.solve_jacobian(ss_final, unknowns=unknowns_tp, targets=targets_tp,
                               inputs=exogenous, T=T)
print("Jacobian OLD done.")

# ── shocks ────────────────────────────────────────────────────────────────────
rho_def = 0.8; rho_Z = 0.8
shock_def = 0.01 * rho_def ** np.arange(T)   # 1pp default shock to D
shock_Z   = 0.01 * rho_Z   ** np.arange(T)   # 1% TFP shock to D

irfs_def_new = G_new.apply({'shock_def_D': shock_def, 'Z_D': np.zeros(T),
                             'Z_F': np.zeros(T), 'shock_def_F': np.zeros(T)})
irfs_def_old = G_old.apply({'shock_def_D': shock_def, 'Z_D': np.zeros(T),
                             'Z_F': np.zeros(T), 'shock_def_F': np.zeros(T)})
irfs_Z_new   = G_new.apply({'Z_D': shock_Z, 'shock_def_D': np.zeros(T),
                             'Z_F': np.zeros(T), 'shock_def_F': np.zeros(T)})
irfs_Z_old   = G_old.apply({'Z_D': shock_Z, 'shock_def_D': np.zeros(T),
                             'Z_F': np.zeros(T), 'shock_def_F': np.zeros(T)})

# ── comparison table ──────────────────────────────────────────────────────────
VARS = ['def_rate_D', 'Y_D', 'n_inter_D', 'spread_rb', 'b_gov_D', 'q_b_D',
        'Y_F', 'n_inter_F', 'C_D', 'C_F']
HORIZONS = [0, 1, 2, 4, 8, 16]

def _get(irf, v, t):
    try:    return float(irf[v][t])
    except: return float('nan')

def print_table(title, irf_new, irf_old, shock_path=None):
    print(f"\n{'='*72}")
    print(f"  {title}")
    print(f"{'='*72}")
    hdr = f"{'Variable':20s}  {'t':>3s}  {'OLD (Y_ss)':>12s}  {'NEW (Y_D)':>12s}  {'diff %':>8s}"
    print(hdr); print('-'*len(hdr))
    for v in VARS:
        for t in HORIZONS:
            old_v = _get(irf_old, v, t)
            new_v = _get(irf_new, v, t)
            if abs(old_v) > 1e-12:
                pct = (new_v - old_v) / abs(old_v) * 100
                pct_s = f"{pct:+8.1f}%"
            else:
                pct_s = "    n/a"
            label = v if t == HORIZONS[0] else ''
            print(f"{label:20s}  {t:>3d}  {old_v:>12.6f}  {new_v:>12.6f}  {pct_s}")
        print()

print_table("DEFAULT SHOCK to D (1pp, rho=0.8) — all variables in log-deviations from SS",
            irfs_def_new, irfs_def_old)
print_table("TFP SHOCK to D (1%, rho=0.8)  — new channel: Y_D fall raises def_rate_D",
            irfs_Z_new, irfs_Z_old)

# ── save results ──────────────────────────────────────────────────────────────
results = {}
for shock_name, irf_new, irf_old in [('default', irfs_def_new, irfs_def_old),
                                       ('tfp',     irfs_Z_new,   irfs_Z_old)]:
    results[shock_name] = {}
    for v in VARS:
        try:
            results[shock_name][v] = {
                'old': list(np.array(irf_old[v])[:30].round(8)),
                'new': list(np.array(irf_new[v])[:30].round(8)),
                'diff': list((np.array(irf_new[v]) - np.array(irf_old[v]))[:30].round(8)),
            }
        except Exception:
            pass

out_path = Path(__file__).parent / 'endo_def_results.json'
with open(out_path, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to {out_path}")
