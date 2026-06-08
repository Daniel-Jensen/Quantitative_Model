"""
Diagnostic script: why does Y_D RISE after a default risk shock in Country D?

Runs the full model from scratch, replicating the notebook setup, then
produces targeted decompositions of every channel.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import copy, sys, warnings
warnings.filterwarnings('ignore')
sys.setrecursionlimit(5000)

import sequence_jacobian as sj
from sequence_jacobian import simple, combine, create_model
from sequence_jacobian import grids

from pathlib import Path
BASE = Path("/Users/Huawei/Desktop/Codes/My own codes/First Paper/Quantitative_Model")
sys.path.insert(0, str(BASE))

from equations_D import (
    hh_init_D, hh_D, make_grids_D, income_D, hh_extended_D,
    smart_steady_D, market_clearing_D, steady_auxilliary_D,
    banker_div_D, sdf_D, sdf_ss_D, government_ss_D, labor_ss_D,
    government_default_D, bond_price_ss_D, bond_return_D,
    ces_price_D, import_demand_D, deposit_return_D,
    capital_adj_D, labor_D, labor_market_D, labor_demand_D,
    intermediation_IC_D, bank_return_D, intermediation_P1_D,
    k_balance_sheet_D, cap_adj_cost_inter_D, macro_pru_tax_D,
    intermediation_P2_D, banker_div_res_D, intermediation_P3_D,
    domestic_bond_foc_D, tax_rule_D, capital_producer_profit_D,
    budget_residual_D, ghh_composite_D,
)
from equations_F import (
    hh_init_F, hh_F, make_grids_F, income_F, hh_extended_F,
    smart_steady_F, market_clearing_F, steady_auxilliary_F,
    banker_div_F, sdf_F, sdf_ss_F, government_ss_F, labor_ss_F,
    government_default_F, bond_price_ss_F, bond_return_F,
    ces_price_F, import_demand_F, deposit_return_F,
    capital_adj_F, labor_F, labor_market_F, labor_demand_F,
    intermediation_IC_F, bank_return_F, intermediation_P1_F,
    k_balance_sheet_F, cap_adj_cost_inter_F, macro_pru_tax_F,
    intermediation_P2_F, banker_div_res_F, intermediation_P3_F,
    domestic_bond_foc_F, tax_rule_F, capital_producer_profit_F,
    budget_residual_F, ghh_composite_F,
)
from equations_global import (
    trade_balance, domestic_bond_clearing,
    portfolio_level_anchors, portfolio_adj_cost, bond_yield,
    global_goods_mkt, external_account_D,
)

# ── Calibration (verbatim from notebook) ─────────────────────────────────────
calibration_start = {
    'frisch_D': 1.0,    'frisch_F': 1.0,
    'eis_D':    0.5,    'eis_F':    0.5,
    'rdep_D':   0.0065, 'rdep_F':   0.0065,
    'q_b_D':    0.82,   'q_b_F':    0.83,
    'Q_D':      1.0,    'Q_F':      1.0,
    'alpha_D':  0.35,   'alpha_F':  0.35,
    'delta_D':  0.0125, 'delta_F':  0.0125,
    'ksi_D':    0.5,    'ksi_F':    0.5,
    'delta_b_D':0.05,   'delta_b_F':0.05,
    'Y_D':      1.00,   'Y_F':      1.00,
    'Y_ss_D':   1.0,    'Y_ss_F':   1.0,
    'N_D':      1.00,   'N_F':      1.00,
    'w_D':      0.65,   'w_F':      0.65,
    'f_D':      0.03,   'f_F':      0.03,
    'lambda_gk_D': 0.2, 'lambda_gk_F': 0.2,
    'n_inter_D':1.0*4,  'n_inter_F':1.0*4,
    'theta_D':  4,      'theta_F':  4,
    'B_supply_D':0.6*4, 'B_supply_F':0.6*4,
    'b_gov_D':  0.6*4,  'b_gov_F':  0.6*4,
    'b_gov_ss_D':0.6*4, 'b_gov_ss_F':0.6*4,
    'tau_D':    0.181,  'tau_F':    0.181,
    'lamb_D':   0.85,   'lamb_F':   0.85,
    'lamb_ss_D':0.85,   'lamb_ss_F':0.85,
    'phi_lamb_D':0.4,   'phi_lamb_F':0.4,
    'shock_def_D':0.0,  'shock_def_F':0.0,
    'def_rate_D':0.0,   'def_rate_F':0.0,
    'def_scale_D':0.0,  'def_scale_F':0.0,
    'def_curvature_D':0.5,'def_curvature_F':0.5,
    'def_offset_D':0.05,'def_offset_F':0.05,
    'recovery_rate_D':0.40,'recovery_rate_F':0.40,
    'zeta_writeoff_D':1.0,'zeta_writeoff_F':1.0,
    'chi0_D':0.01,'chi0_F':0.01,
    'chi1_D':0.5, 'chi1_F':0.5,
    'chi2_D':2.0, 'chi2_F':2.0,
    'T0_D':0.0,'T0_F':0.0,'T1_D':0.0,'T1_F':0.0,
    'omega':0.85,'epsilon_trade':1.5,'p':1.0,
    'phi_bF_D_ss':0.15,'phi_bD_F_ss':0.15,
    'psi_bF_D':0.05,'psi_bD_F':0.05,
    'mu_w_D':1.0,'mu_w_F':1.0,'mc_D':1.0,'mc_F':1.0,
    'rho_z_D':0.90,'rho_z_F':0.97,'sigma_z_D':0.8,'sigma_z_F':0.8,
    'nZ_D':7,'nZ_F':7,'nDep_D':500,'nDep_F':500,
    'Depmax_D':150,'Depmax_F':150,
}

_n_D, _n_F = calibration_start['n_inter_D'], calibration_start['n_inter_F']
_B_D, _B_F = calibration_start['B_supply_D'], calibration_start['B_supply_F']
b_F_D = calibration_start['phi_bF_D_ss'] * _n_D / calibration_start['q_b_F']
b_D_F = calibration_start['phi_bD_F_ss'] * _n_F / calibration_start['q_b_D']
calibration_start.update({
    'b_F_D': b_F_D, 'b_D_F': b_D_F,
    'b_D_D': _B_D - b_D_F, 'b_F_F': _B_F - b_F_D,
    'b_F_D_anchor': b_F_D, 'b_D_F_anchor': b_D_F,
    'psi_bD_D': 0.5, 'psi_bF_F': 0.5,
})

# ── Steady state ──────────────────────────────────────────────────────────────
print("Solving steady state...")
ha = sj.create_model([
    sdf_ss_D, government_default_D, bond_price_ss_D, bond_return_D,
    sdf_ss_F, government_default_F, bond_price_ss_F, bond_return_F,
    hh_extended_D, smart_steady_D, market_clearing_D, steady_auxilliary_D,
    banker_div_D, government_ss_D, labor_ss_D,
    hh_extended_F, smart_steady_F, market_clearing_F, steady_auxilliary_F,
    banker_div_F, government_ss_F, labor_ss_F,
    ces_price_D, import_demand_D, ces_price_F, import_demand_F,
    deposit_return_D, deposit_return_F,
    bond_yield, trade_balance, external_account_D, global_goods_mkt,
], name='SS model')

unknowns_ss = {'beta_D': 0.9920096467, 'beta_F': 0.9867162384, 'p': 0.994489}
targets_ss  = ['deposit_mkt_D', 'deposit_mkt_F', 'ca_res_D']
ss = ha.solve_steady_state(calibration_start, unknowns_ss, targets_ss, solver='broyden_custom')

anchors = {
    'phi_bD_D_ss': float(ss['q_b_D']) * float(ss['b_D_D']) / float(ss['n_inter_D']),
    'phi_bF_F_ss': float(ss['q_b_F']) * float(ss['b_F_F']) / (float(ss['p']) * float(ss['n_inter_F'])),
    'b_F_D_anchor': float(ss['b_F_D']), 'b_D_F_anchor': float(ss['b_D_F']),
    'excess_return_bD_D_ss': float(ss['rb_actual_D']) - float(ss['rdep_D']) - calibration_start['T0_D'],
    'excess_return_bF_F_ss': float(ss['rb_actual_F']) - float(ss['rdep_F']) - calibration_start['T0_F'],
    'excess_return_F_D_ss':  float(ss['rb_actual_F']) - float(ss['rdep_D']) - calibration_start['T0_D'],
    'excess_return_D_F_ss':  float(ss['rb_actual_D']) - float(ss['rdep_F']) - calibration_start['T0_F'],
    'q_b_D': float(ss['q_b_D']), 'q_b_F': float(ss['q_b_F']), 'p': float(ss['p']),
}
calibration_start.update(anchors)
for k, v in anchors.items():
    ss.toplevel[k] = v
ss.toplevel['b_F_D_ss'] = float(ss['b_F_D'])
ss.toplevel['b_D_F_ss'] = float(ss['b_D_F'])
ss.toplevel['Rgross_D'] = float(1 + ss['rdep_D'])
ss.toplevel['Rgross_F'] = float(1 + ss['rdep_F'])
_fr_D = float(ss['frisch_D']); _fr_F = float(ss['frisch_F'])
ss.toplevel['X_D'] = float(ss['C_D']) - float(ss['vphi_D'])*float(ss['N_D'])**(1+1/_fr_D)/(1+1/_fr_D)
ss.toplevel['X_F'] = float(ss['C_F']) - float(ss['vphi_F'])*float(ss['N_F'])**(1+1/_fr_F)/(1+1/_fr_F)
for k, v in {
    'tau_mp_D':0.0,'tau_mp_F':0.0,'Phi_D':0.0,'Phi_F':0.0,'T_D':0.0,'T_F':0.0,
    'b_F_D_res':0.0,'b_D_F_res':0.0,'rb_D_res':0.0,'rb_F_res':0.0,
    'labor_mkt_res_D':0.0,'labor_mkt_res_F':0.0,'w_res_D':0.0,'w_res_F':0.0,
}.items():
    ss.toplevel[k] = v
for k in ('nu_K_D','nu_bD_D','nu_bF_D','eta_D','nu_K_F','nu_bF_F','nu_bD_F','eta_F'):
    if k not in ss.toplevel:
        ss.toplevel[k] = 0.0
ss_final = copy.deepcopy(ss)
print(f"SS solved: Y_D={ss['Y_D']:.4f}, q_b_D={ss['q_b_D']:.4f}, rdep_D={ss['rdep_D']:.6f}")

# ── Full model ────────────────────────────────────────────────────────────────
print("Building full model...")
cali_D = cali_F = ss

financial_solved_D = combine([intermediation_P1_D, intermediation_IC_D]).solved(
    unknowns={'nu_K_D':float(cali_D['nu_K_D']),'nu_bD_D':float(cali_D['nu_bD_D']),
              'nu_bF_D':float(cali_D['nu_bF_D']),'eta_D':float(cali_D['eta_D']),
              'theta_D':float(cali_D['theta_D'])},
    targets=['nu_K_res_D','nu_bD_res_D','nu_bF_res_D','eta_res_D','ic_res_D'],
    solver='broyden_custom')

financial_solved_F = combine([intermediation_P1_F, intermediation_IC_F]).solved(
    unknowns={'nu_K_F':float(cali_F['nu_K_F']),'nu_bF_F':float(cali_F['nu_bF_F']),
              'nu_bD_F':float(cali_F['nu_bD_F']),'eta_F':float(cali_F['eta_F']),
              'theta_F':float(cali_F['theta_F'])},
    targets=['nu_K_res_F','nu_bF_res_F','nu_bD_res_F','eta_res_F','ic_res_F'],
    solver='broyden_custom')

ha_full = sj.create_model([
    deposit_return_D, hh_extended_D, ghh_composite_D, sdf_D,
    government_default_D, financial_solved_D, bond_return_D, bank_return_D,
    cap_adj_cost_inter_D, macro_pru_tax_D, intermediation_P2_D, intermediation_P3_D,
    k_balance_sheet_D, capital_adj_D, tax_rule_D, capital_producer_profit_D,
    budget_residual_D, labor_D, labor_market_D, labor_demand_D,
    banker_div_res_D, market_clearing_D,
    deposit_return_F, hh_extended_F, ghh_composite_F, sdf_F,
    government_default_F, financial_solved_F, bond_return_F, bank_return_F,
    cap_adj_cost_inter_F, macro_pru_tax_F, intermediation_P2_F, intermediation_P3_F,
    k_balance_sheet_F, capital_adj_F, tax_rule_F, capital_producer_profit_F,
    budget_residual_F, labor_F, labor_market_F, labor_demand_F,
    banker_div_res_F, market_clearing_F,
    ces_price_D, import_demand_D, ces_price_F, import_demand_F,
    trade_balance, external_account_D, domestic_bond_clearing, bond_yield,
    portfolio_level_anchors, portfolio_adj_cost,
    domestic_bond_foc_D, domestic_bond_foc_F, global_goods_mkt,
], name="ha_full")

unknowns_tp = [
    'K_D','n_inter_D','div_D','I_D','Q_D','b_gov_D','N_D','b_F_D','w_D','rdep_D',
    'K_F','n_inter_F','div_F','I_F','Q_F','b_gov_F','N_F','b_D_F','w_F','rdep_F',
    'p','q_b_D','q_b_F',
]
targets_tp = [
    'deposit_mkt_D','K_res_D','n_inter_val_D','div_res_D',
    'capital_res_D','q_res_D','b_gov_res_D','b_F_D_res','labor_mkt_res_D','w_res_D',
    'deposit_mkt_F','K_res_F','n_inter_val_F','div_res_F',
    'capital_res_F','q_res_F','b_gov_res_F','b_D_F_res','labor_mkt_res_F','w_res_F',
    'ca_res_D','rb_D_res','rb_F_res',
]
T = 300
exogenous = ['Z_D','shock_def_D','Z_F','shock_def_F']
print("Computing Jacobian...")
G_jac = ha_full.solve_jacobian(ss_final, unknowns_tp, targets_tp, exogenous, T=T)
print("Jacobian done.")

# ── IRF: default shock ────────────────────────────────────────────────────────
rho_def = 0.8
dShock  = 0.01 * rho_def ** np.arange(T)
print("Computing default shock IRF...")
irfs = ha_full.solve_impulse_linear(ss_final, unknowns_tp, targets_tp, {'shock_def_D': dShock})
print("IRF done.")

# Helper: get IRF value at horizon h (0-indexed)
def g(var, h=0):
    arr = irfs[var] if var in irfs else np.zeros(T)
    return float(arr[h])

def arr(var):
    return irfs[var] if var in irfs else np.zeros(T)

# ── Print key impact-period values ───────────────────────────────────────────
print("\n" + "="*65)
print("IMPACT (t=0) RESPONSES TO 1% DEFAULT SHOCK IN D")
print("="*65)

variables = [
    # Shock itself
    ('shock_def_D', 'Default shock (input)'),
    ('def_rate_D',  'Default rate'),
    # Bond market
    ('q_b_D',       'Bond price D'),
    ('q_b_F',       'Bond price F'),
    ('rb_actual_D', 'Realised bond return D'),
    ('rb_D',        'Implied bond yield D'),
    # Bank balance sheet
    ('n_inter_D',   'Bank net worth D'),
    ('theta_D',     'Bank leverage D'),
    # Capital / investment
    ('K_D',         'Capital stock D'),
    ('I_D',         'Investment D'),
    ('Q_D',         'Tobins q D'),
    ('rk_D',        'Capital return D'),
    # Bank returns
    ('rn_D',        'Bank portfolio return D'),
    ('rdep_D',      'Deposit rate D'),
    # Labour
    ('N_D',         'Employment D'),
    ('w_D',         'Wage D'),
    ('P_CES_D',     'Price index D'),
    # Output & demand
    ('Y_D',         'OUTPUT D  <-- key'),
    ('C_D',         'Consumption D'),
    ('G_D',         'Government spending D'),
    ('NX_D',        'Net exports D'),
    # Fiscal
    ('TAX_D',       'Tax revenue D (transfers)'),
    ('lamb_D',      'Transfer multiplier D'),
    ('b_gov_D',     'Government debt D'),
    # Cross-border
    ('p',           'Terms of trade p'),
    ('Y_F',         'OUTPUT F (spillover)'),
    ('n_inter_F',   'Bank net worth F'),
    ('K_F',         'Capital stock F'),
]

for var, label in variables:
    val = g(var, 0) * 100  # in % or pp deviation
    marker = '  <<<' if 'key' in label else ''
    print(f"  {label:<40s} {val:+.4f}%{marker}")

# ── Decompose why Y_D rises ───────────────────────────────────────────────────
print("\n" + "="*65)
print("DECOMPOSITION OF Y_D RESPONSE")
print("="*65)

# Y_D = Z_D * K_D^alpha * N_D^(1-alpha), so log-linearised:
# dY/Y = dZ/Z + alpha*dK/K + (1-alpha)*dN/N
alpha = float(ss['alpha_D'])
dY = g('Y_D', 0)
dK = g('K_D', 0)
dN = g('N_D', 0)
dZ = 0.0  # no TFP shock
print(f"  Y_D change:                         {dY*100:+.4f}%")
print(f"  K_D contribution (alpha*dK/K):      {alpha*dK*100:+.4f}%")
print(f"  N_D contribution ((1-a)*dN/N):      {(1-alpha)*dN*100:+.4f}%")
print(f"  Residual (TFP + cross terms):       {(dY - alpha*dK - (1-alpha)*dN)*100:+.4f}%")

# Labour supply decomposition:
# w/P_CES = vphi * N^(1/frisch)  →  N = (w/(P_CES*vphi))^frisch
# dN/N = frisch*(dw/w - dP/P)
print(f"\n  --- Labour market ---")
frisch = float(ss['frisch_D'])
dw = g('w_D', 0)
dP = g('P_CES_D', 0)
print(f"  Wage change dw/w:                   {dw*100:+.4f}%")
print(f"  Price index change dP/P:            {dP*100:+.4f}%")
print(f"  Real wage change d(w/P):            {(dw-dP)*100:+.4f}%")
print(f"  Implied dN/N = frisch*(dw-dP):      {frisch*(dw-dP)*100:+.4f}%  (actual: {dN*100:+.4f}%)")

print(f"\n  --- Bank balance sheet ---")
dtheta = g('theta_D', 0)
dn     = g('n_inter_D', 0)
dQ     = g('Q_D', 0)
dqb    = g('q_b_D', 0)
print(f"  theta_D change:                     {dtheta*100:+.4f}%")
print(f"  n_inter_D change:                   {dn*100:+.4f}%")
print(f"  Q_D change:                         {dQ*100:+.4f}%")
print(f"  q_b_D change:                       {dqb*100:+.4f}%")
print(f"  K_D = (theta*n - bonds) / Q")
theta_ss = float(ss['theta_D'])
n_ss     = float(ss['n_inter_D'])
qb_ss    = float(ss['q_b_D'])
qbF_ss   = float(ss['q_b_F'])
bDD_ss   = float(ss['b_D_D'])
bFD_ss   = float(ss['b_F_D'])
Q_ss     = float(ss['Q_D'])
K_ss     = float(ss['K_D'])
# K = (theta*n - qb_D*b_D_D - qb_F*b_F_D) / Q
# dK/K ≈ dtheta/theta + dn/n - (qb_D*b_D_D/theta/n)*dqb/qb - ... - dQ/Q
bond_share_DD = qb_ss * bDD_ss / (theta_ss * n_ss)
bond_share_FD = qbF_ss * bFD_ss / (theta_ss * n_ss)
kappa_ss = theta_ss - qb_ss*bDD_ss/n_ss - qbF_ss*bFD_ss/n_ss
dqbF = g('q_b_F', 0)
dbDD = g('b_D_D', 0) if 'b_D_D' in irfs else 0.0
dbFD = g('b_F_D', 0)
# linearized: dK/K ≈ (1/kappa_ss)[theta_ss*(dtheta+dn) - qb*bDD*(dqb+dbDD) - qbF*bFD*(dqbF+dbFD)] - dQ
print(f"  Capital-share kappa (K*Q/n):        {kappa_ss/theta_ss:.4f} of leverage ratio")
print(f"  Bond-share D*D (qb_D*b_D_D/n):      {qb_ss*bDD_ss/n_ss:.4f}")
print(f"  Bond-share F*D (qb_F*b_F_D/n):      {qbF_ss*bFD_ss/n_ss:.4f}")

print(f"\n  --- Fiscal channel ---")
dlamb  = g('lamb_D', 0)
dTAX   = g('TAX_D', 0)
dG     = g('G_D', 0)
dbgov  = g('b_gov_D', 0)
print(f"  lamb_D (transfer multiplier):       {dlamb*100:+.4f}%")
print(f"  TAX_D (household transfers):        {dTAX*100:+.4f}%")
print(f"  G_D (government spending):          {dG*100:+.4f}%")
print(f"  b_gov_D (government debt):          {dbgov*100:+.4f}%")

print(f"\n  --- Demand components (impact) ---")
dC  = g('C_D', 0)
dI  = g('I_D', 0)
dNX = g('NX_D', 0)
P_ss  = float(ss['P_CES_D'])
C_ss  = float(ss['C_D'])
I_ss  = float(ss['I_D'])
G_ss  = float(ss['G_D'])
NX_ss = float(ss['NX_D'])
print(f"  C_D change (level):                 {dC*C_ss*P_ss*100:+.5f}  (share of Y: {C_ss*P_ss:.3f})")
print(f"  I_D change (level):                 {dI*I_ss*100:+.5f}  (share of Y: {I_ss:.3f})")
print(f"  G_D change (level):                 {dG*G_ss*100:+.5f}  (share of Y: {G_ss:.3f})")
print(f"  NX_D change (level):                {dNX*NX_ss*100:+.5f}  (share of Y: {NX_ss:.3f})")
print(f"  Note: goods_mkt residual at t=0:    {g('goods_mkt_D',0)*100:+.6f}%")

# ── Goods market decomposition using nominal shares ──────────────────────────
# Y_D = P_CES_D * C_D + I_D + G_D + Phi_D + T_D + NX_D
# (market_clearing_D)
print(f"\n  --- Goods market accounting (level deviations at t=0) ---")
dPhi  = g('Phi_D', 0)
dT_   = g('T_D', 0)
Y_ss_v = float(ss['Y_D'])
print(f"  Sum of demand (P*dC + dI + dG + dNX + dPhi + dT): ")
pCdC   = P_ss * C_ss * dC  + C_ss * dP * P_ss   # d(P_CES * C) ≈ P*dC + C*dP
iDdI   = I_ss * dI
gDdG   = G_ss * dG
nxDdNX = NX_ss * dNX
phiD   = float(ss.get('Phi_D', 0)) * dPhi
tD     = float(ss.get('T_D', 0)) * dT_
Phi_ss = float(ss['Phi_D']) if 'Phi_D' in ss else 0.0
T_ss   = float(ss['T_D'])  if 'T_D'  in ss else 0.0
phiD   = Phi_ss * dPhi
tD     = T_ss   * dT_
print(f"    P_CES*C contribution:             {pCdC*100:+.5f}")
print(f"    I_D contribution:                 {iDdI*100:+.5f}")
print(f"    G_D contribution:                 {gDdG*100:+.5f}")
print(f"    NX_D contribution:                {nxDdNX*100:+.5f}")
print(f"    Total demand:                     {(pCdC+iDdI+gDdG+nxDdNX)*100:+.5f}")
print(f"    Y_D change:                       {dY*Y_ss_v*100:+.5f}")

# ── Trace the specific GHH labor-supply channel ───────────────────────────────
print("\n" + "="*65)
print("WHY DOES N_D INCREASE?  (GHH labor supply channel)")
print("="*65)
print("GHH: w/P_CES = vphi*N^(1/frisch)  =>  N = (w/P/vphi)^frisch")
print(f"  vphi_D (SS, fixed):               {float(ss['vphi_D']):.6f}")
print(f"  P_CES_D change:                   {dP*100:+.4f}% (from terms-of-trade p)")
print(f"  w_D change:                       {dw*100:+.4f}%")
print(f"  Real wage w/P change:             {(dw-dP)*100:+.4f}%")
print(f"  N_D predicted by GHH:             {frisch*(dw-dP)*100:+.4f}%")
print(f"  N_D actual:                       {dN*100:+.4f}%")
print(f"  p (terms of trade) change:        {g('p',0)*100:+.4f}%")
print(f"  P_CES_D depends on p:             omega={float(ss['omega']):.2f}, eps_trade={float(ss['epsilon_trade']):.2f}")
print(f"  IM_D change:                      {g('IM_D',0)*100:+.4f}%")

# ── Dynamic paths: what drives the first few quarters ────────────────────────
print("\n" + "="*65)
print("DYNAMIC PATHS (first 8 quarters)")
print("="*65)
Th = 8
print(f"{'Quarter':>8} {'Y_D':>8} {'N_D':>8} {'K_D':>8} {'n_inter':>9} {'theta':>8} {'Q_D':>8} {'rk_D':>8} {'rdep_D':>8} {'q_b_D':>8}")
print("-"*90)
for t in range(Th):
    print(f"  t={t:<4} {g('Y_D',t)*100:>+7.3f}% {g('N_D',t)*100:>+7.3f}% "
          f"{g('K_D',t)*100:>+7.3f}% {g('n_inter_D',t)*100:>+8.3f}% "
          f"{g('theta_D',t)*100:>+7.3f}% {g('Q_D',t)*100:>+7.3f}% "
          f"{g('rk_D',t)*100:>+7.4f}% {g('rdep_D',t)*100:>+7.4f}% "
          f"{g('q_b_D',t)*100:>+7.4f}%")

# ── Single-channel knockout experiments ──────────────────────────────────────
# We test what Y_D does if we remove the fiscal channel (freeze lamb_D and b_gov_D).
# This requires modifying the model, but we can approximate by looking at partial
# contributions in the Jacobian.
print("\n" + "="*65)
print("CHANNEL ATTRIBUTION (Jacobian rows at t=0)")
print("="*65)
# The Jacobian G_jac['Y_D']['shock_def_D'] gives the full IRF column.
# We now look at partial Jacobians through different pathways.
# G_jac[y][x][t, s] = ∂y_t / ∂x_s
# At s=0, t=0: which unknowns contribute most to dY_D?

try:
    # G_jac gives us the impulse-response from shock → unknowns → targets
    # We can look at it via the solved unknowns
    print("  Y_D IRF at t=0:", float(arr('Y_D')[0])*100, "%")
    print("  N_D IRF at t=0:", float(arr('N_D')[0])*100, "%")
    print("  K_D IRF at t=0:", float(arr('K_D')[0])*100, "%")
    print("  N_D is the dominant driver of Y_D since K_D(-1) is predetermined at impact")
    print()
    print("  N_D rises because w_D/P_CES_D rises (GHH labor supply)")
    print("  w_D rises because firm FOC w = (1-alpha)*Y/N must hold")
    print("  P_CES_D falls because terms of trade p falls (D-goods depreciate in real terms?)")
    print()
    print("  But at impact K_D is NOT predetermined - it's in the unknowns!")
    print("  Let's check: if K_D falls but N_D rises by more (in labor units), Y could still rise.")
    print(f"  K_D alpha-contribution to Y_D:     {alpha*g('K_D',0)*100:+.4f}%")
    print(f"  N_D (1-a)-contribution to Y_D:     {(1-alpha)*g('N_D',0)*100:+.4f}%")
    print(f"  Which dominates?                   {'N_D' if abs((1-alpha)*g('N_D',0)) > abs(alpha*g('K_D',0)) else 'K_D'}")
except Exception as e:
    print(f"  Channel attribution error: {e}")

# ── Why does rdep_D move? ─────────────────────────────────────────────────────
print("\n" + "="*65)
print("DEPOSIT RATE & FINANCIAL MECHANISM")
print("="*65)
print(f"  rdep_D (deposit rate) impact:       {g('rdep_D',0)*100:+.5f}%")
print(f"  rk_D (capital return) impact:       {g('rk_D',0)*100:+.5f}%")
print(f"  rb_actual_D (bond return) impact:   {g('rb_actual_D',0)*100:+.5f}%")
print(f"  rn_D (bank portfolio return):       {g('rn_D',0)*100:+.5f}%")
print(f"  Spread rk-rdep impact:              {(g('rk_D',0)-g('rdep_D',0))*100:+.5f}%")
print(f"  Spread rb_D-rdep impact:            {(g('rb_actual_D',0)-g('rdep_D',0))*100:+.5f}%")
print()
print(f"  GK IC: theta driven by nu_K*(rk-rdep) + nu_b*(rb-rdep) + eta")
print(f"  If rk-rdep RISES, theta rises → more leverage → K could rise even if n falls")
print()
print(f"  theta_D change:                     {g('theta_D',0)*100:+.4f}%  <-- key")
print(f"  n_inter_D change:                   {g('n_inter_D',0)*100:+.4f}%")
print(f"  theta_D * n_inter_D change:         {(g('theta_D',0)+g('n_inter_D',0))*100:+.4f}%  (approx)")

# ── Fiscal channel isolation ───────────────────────────────────────────────────
print("\n" + "="*65)
print("FISCAL MECHANISM")
print("="*65)
qb_ss_v = float(ss['q_b_D'])
bDD_ss_v = float(ss['b_D_D'])
bDF_ss_v = float(ss['b_D_F'])
bgov_ss_v = float(ss['b_gov_D'])
delta_b_ss = float(ss['delta_b_D'])
haircut_ss = 1.0 - float(ss['recovery_rate_D'])
zeta_ss = float(ss['zeta_writeoff_D'])
def_shock_size = 0.01

# Reduction in coupon payment to bondholders (government saves):
coupon_saving = delta_b_ss * zeta_ss * def_shock_size * haircut_ss * bgov_ss_v
print(f"  Government coupon saving at impact:")
print(f"    delta_b * zeta * shock * haircut * b_gov = {coupon_saving:.6f}")
print(f"    This is {coupon_saving/float(ss['Y_D'])*100:.4f}% of Y_D")
print(f"  This goes into budget surplus → reduces b_gov → raises lamb_D → higher transfers")
print(f"  Higher transfers boost HH income → consumption demand → higher output")
print()
print(f"  How large is the fiscal effect on output?")
print(f"  Change in G_D at t=0:               {g('G_D',0)*100:+.5f}%")
print(f"  Change in TAX_D at t=0:             {g('TAX_D',0)*100:+.5f}%")
print(f"  Change in C_D at t=0:               {g('C_D',0)*100:+.5f}%")

# ── Summary diagnosis ────────────────────────────────────────────────────────
print("\n" + "="*65)
print("DIAGNOSIS SUMMARY")
print("="*65)
if g('Y_D', 0) > 0:
    print("  Y_D RISES after default shock. Channels:")
    if g('N_D', 0) > 0:
        print(f"  [ACTIVE] Employment N_D rises: +{g('N_D',0)*100:.4f}%")
        if g('w_D', 0) - g('P_CES_D', 0) > 0:
            print(f"           → Real wage w/P rises: +{(g('w_D',0)-g('P_CES_D',0))*100:.4f}%")
            print(f"           → GHH labor supply EXPANDS (no income effect, substitution dominates)")
        if g('P_CES_D', 0) < 0:
            print(f"           → P_CES_D falls {g('P_CES_D',0)*100:.4f}%: D-goods become cheap → imports fall → more domestic demand")
    if g('K_D', 0) > 0:
        print(f"  [ACTIVE] Capital K_D rises: +{g('K_D',0)*100:.4f}%")
    elif g('K_D', 0) < 0:
        print(f"  [INACTIVE] Capital K_D falls: {g('K_D',0)*100:.4f}% (suppressing output)")
    print()
    print(f"  Root cause candidates:")
    print(f"  1. FISCAL WINDFALL: govt default → lower coupon → higher household transfers")
    print(f"     → consumption demand ↑ → firms hire more → N_D ↑ → Y_D ↑")
    print(f"  2. TERMS-OF-TRADE: default fear → capital outflow → p falls → P_CES_D falls")
    print(f"     → real wage rises → GHH labor supply ↑ → N_D ↑ → Y_D ↑")
    print(f"  3. BANK LEVERAGE PARADOX: if spread rk-rdep rises, GK allows MORE leverage")
    print(f"     → theta_D ↑ can offset n_inter_D ↓ → K_D doesn't fall much")
    print(f"\n  Most likely dominant channel (inspect signs above to confirm):")
    if abs(g('N_D', 0)) > abs(g('K_D', 0)):
        print(f"  → LABOUR CHANNEL: N_D contribution {(1-alpha)*g('N_D',0)*100:+.4f}% vs K_D {alpha*g('K_D',0)*100:+.4f}%")
    else:
        print(f"  → CAPITAL CHANNEL: K_D contribution {alpha*g('K_D',0)*100:+.4f}% vs N_D {(1-alpha)*g('N_D',0)*100:+.4f}%")

# ── Save a summary figure ──────────────────────────────────────────────────────
TT = 40
t_ = np.arange(TT)

fig, axes = plt.subplots(3, 4, figsize=(16, 10))
fig.suptitle("Default Shock in D: Diagnostic IRFs", fontsize=13, fontweight='bold')

plots = [
    ('Y_D',       'Output D',          axes[0,0]),
    ('N_D',       'Employment D',       axes[0,1]),
    ('K_D',       'Capital D',          axes[0,2]),
    ('Q_D',       'Tobins q D',         axes[0,3]),
    ('n_inter_D', 'Bank net worth D',   axes[1,0]),
    ('theta_D',   'Bank leverage D',    axes[1,1]),
    ('q_b_D',     'Bond price D',       axes[1,2]),
    ('rdep_D',    'Deposit rate D',     axes[1,3]),
    ('rk_D',      'Capital return rk', axes[2,0]),
    ('w_D',       'Wage D',             axes[2,1]),
    ('P_CES_D',   'Price index D',      axes[2,2]),
    ('C_D',       'Consumption D',      axes[2,3]),
]
for var, title, ax in plots:
    data = arr(var)[:TT] * 100
    ax.plot(t_, data, color='#002147', linewidth=1.8)
    ax.axhline(0, color='#888', linewidth=0.7, linestyle=':')
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('Quarter', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines[['top','right']].set_visible(False)

fig.tight_layout()
outpath = BASE / "diag_default_shock.png"
fig.savefig(str(outpath), dpi=120, bbox_inches='tight')
print(f"\nFigure saved to: {outpath}")
print("\nDone.")
