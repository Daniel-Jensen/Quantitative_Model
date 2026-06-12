"""TPI audit: does cb_buy_D inject unbacked resources (missing CB budget constraint)?

Prediction (walras_forensics.md §2 extended): with b_D_D = b_gov_D - b_D_F - cb_buy_D
and no CB budget anywhere, the D aggregate identity becomes
  goods_mkt_D + ca_res_D = factor_gap_D + [q_D*cb_t - (1+rb_D,t)*q_D(-1)*cb_{t-1}]
so ca_res_D under TPI = baseline-def-shock leak + CB hole. Test at gamma=10.
Also: compare discounted welfare gain dW_D vs discounted CB injection.
"""
import json, sys, copy
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import os
os.chdir(ROOT)
import sequence_jacobian as sj
from sequence_jacobian import simple, combine

RES = {}
def jlog(k, v):
    RES[k] = v
    print(f"[tpi] {k}: {v}")

src = open(ROOT/'audit_artifacts'/'run_audit.py').read()
ss_part = src.split("# ═══════════════ CELL 20")[0]
ns = {'__file__': str(ROOT/'audit_artifacts'/'run_audit.py')}
exec(compile(ss_part, 'ss_pipeline', 'exec'), ns)
ss_final = ns['ss_final']
cali = ss_final
calibration_start = ns['calibration_start']

from equations_D import (hh_extended_D, capital_adj_D, labor_D, labor_market_D, labor_demand_D,
    intermediation_IC_D, bank_return_D, intermediation_P1_D, k_balance_sheet_D,
    cap_adj_cost_inter_D, macro_pru_tax_D, intermediation_P2_D, banker_div_res_D,
    intermediation_P3_D, divert_bond_foc_D, tax_rule_D, capital_producer_profit_D,
    budget_residual_D, ces_price_D, import_demand_D, deposit_return_D,
    bond_return_D, sdf_D, sdf_banker_D, ghh_composite_D, welfare_agg_D,
    government_default_D, market_clearing_D)
from equations_F import (hh_extended_F, capital_adj_F, labor_F, labor_market_F, labor_demand_F,
    intermediation_IC_F, bank_return_F, intermediation_P1_F, k_balance_sheet_F,
    cap_adj_cost_inter_F, macro_pru_tax_F, intermediation_P2_F, banker_div_res_F,
    intermediation_P3_F, divert_bond_foc_F, tax_rule_F, capital_producer_profit_F,
    budget_residual_F, ces_price_F, import_demand_F, deposit_return_F,
    bond_return_F, sdf_F, sdf_banker_F, ghh_composite_F, welfare_agg_F,
    government_default_F, market_clearing_F)
from equations_global import (trade_balance, portfolio_level_anchors, divert_portfolio_adj,
    bond_yield, global_goods_mkt, external_account_D)

@simple
def domestic_bond_clearing_tpi(b_gov_D, b_gov_F, b_D_F, b_F_D, cb_buy_D):
    b_D_D = b_gov_D - b_D_F - cb_buy_D
    b_F_F = b_gov_F - b_F_D
    return b_D_D, b_F_F

# Mirrors the TPI-1 audit fix now in model_v12.ipynb cell TPI-1.
@simple
def budget_residual_D_tpi(b_gov_D, G_D, TAX_D, q_b_D, def_rate_D, recovery_rate_D,
                          zeta_writeoff_D, P_CES_D, delta_b_D, writeoff_enabled_D,
                          cb_buy_D):
    haircut_D      = 1.0 - recovery_rate_D
    haircut_mult_D = writeoff_enabled_D
    surv_cont_D    = 1.0 - zeta_writeoff_D * def_rate_D * haircut_D * haircut_mult_D
    coupon_D       = delta_b_D * (1.0 - def_rate_D * haircut_D * haircut_mult_D) * b_gov_D(-1)
    net_issuance_D = q_b_D * (b_gov_D - surv_cont_D * (1.0 - delta_b_D) * b_gov_D(-1))
    rem_cb_D       = (delta_b_D * (1.0 - def_rate_D * haircut_D * haircut_mult_D) * cb_buy_D(-1)
                      + q_b_D * surv_cont_D * (1.0 - delta_b_D) * cb_buy_D(-1)
                      - q_b_D * cb_buy_D)
    b_gov_res_D    = coupon_D + G_D - P_CES_D * TAX_D - net_issuance_D - rem_cb_D
    return b_gov_res_D, rem_cb_D

sys.setrecursionlimit(5000)
financial_solved_D = combine([intermediation_P1_D, intermediation_IC_D]).solved(
    unknowns={'nu_K_D': float(cali['nu_K_D']), 'nu_bD_D': float(cali['nu_bD_D']),
              'nu_bF_D': float(cali['nu_bF_D']), 'eta_D': float(cali['eta_D']),
              'theta_D': float(cali['theta_D'])},
    targets=['nu_K_res_D','nu_bD_res_D','nu_bF_res_D','eta_res_D','ic_res_D'],
    solver='broyden_custom')
financial_solved_F = combine([intermediation_P1_F, intermediation_IC_F]).solved(
    unknowns={'nu_K_F': float(cali['nu_K_F']), 'nu_bF_F': float(cali['nu_bF_F']),
              'nu_bD_F': float(cali['nu_bD_F']), 'eta_F': float(cali['eta_F']),
              'theta_F': float(cali['theta_F'])},
    targets=['nu_K_res_F','nu_bF_res_F','nu_bD_res_F','eta_res_F','ic_res_F'],
    solver='broyden_custom')

ha_full_tpi = sj.create_model([
    deposit_return_D, tax_rule_D, hh_extended_D, ghh_composite_D, sdf_D, sdf_banker_D,
    government_default_D, financial_solved_D, bond_return_D, bank_return_D,
    cap_adj_cost_inter_D, macro_pru_tax_D, intermediation_P2_D, intermediation_P3_D,
    k_balance_sheet_D, capital_adj_D, capital_producer_profit_D, budget_residual_D_tpi,
    labor_D, labor_market_D, labor_demand_D, banker_div_res_D, market_clearing_D, welfare_agg_D,
    deposit_return_F, tax_rule_F, hh_extended_F, ghh_composite_F, sdf_F, sdf_banker_F,
    government_default_F, financial_solved_F, bond_return_F, bank_return_F,
    cap_adj_cost_inter_F, macro_pru_tax_F, intermediation_P2_F, intermediation_P3_F,
    k_balance_sheet_F, capital_adj_F, capital_producer_profit_F, budget_residual_F,
    labor_F, labor_market_F, labor_demand_F, banker_div_res_F, market_clearing_F, welfare_agg_F,
    ces_price_D, import_demand_D, ces_price_F, import_demand_F,
    trade_balance, external_account_D, domestic_bond_clearing_tpi, bond_yield,
    portfolio_level_anchors, divert_portfolio_adj, divert_bond_foc_D, divert_bond_foc_F,
    global_goods_mkt,
], name="TPI")

ss_tpi = copy.deepcopy(ss_final)
ss_tpi.toplevel['cb_buy_D'] = 0.0

unknowns_tp = ['K_D','n_inter_D','div_D','I_D','Q_D','b_gov_D','N_D','b_F_D','w_D','rdep_D',
               'K_F','n_inter_F','div_F','I_F','Q_F','b_gov_F','N_F','b_D_F','w_F','rdep_F',
               'p','q_b_D','q_b_F']
targets_tp = ['deposit_mkt_D','K_res_D','n_inter_val_D','div_res_D','capital_res_D','q_res_D',
              'b_gov_res_D','b_F_D_res','labor_mkt_res_D','w_res_D',
              'deposit_mkt_F','K_res_F','n_inter_val_F','div_res_F','capital_res_F','q_res_F',
              'b_gov_res_F','b_D_F_res','labor_mkt_res_F','w_res_F',
              'goods_mkt_D','rb_D_res','rb_F_res']
T = 500
print("Computing G_tpi...")
G_tpi = ha_full_tpi.solve_jacobian(ss_tpi, unknowns=unknowns_tp, targets=targets_tp,
                                   inputs=['Z_D','shock_def_D','Z_F','shock_def_F','cb_buy_D'], T=T)
print("done")

rho = 0.8
dShock = 0.01*rho**np.arange(T)
zero = np.zeros(T)

A_def = np.array(G_tpi['spread_rb']['shock_def_D'])
A_cb  = np.array(G_tpi['spread_rb']['cb_buy_D'])
qss = float(ss_final['q_b_D']); rbss = float(ss_final['rb_actual_D'])
Css = float(ss_final['C_D']); betaD = float(ss_final['beta_D']); betaF = float(ss_final['beta_F'])

summary = {}
for gam in [0, 2, 5, 10]:
    if gam == 0:
        cb_path = zero
    else:
        spread_cl = np.linalg.solve(np.eye(T) - gam*A_cb, A_def @ dShock)
        cb_path = gam*spread_cl
    irf = G_tpi @ {'Z_D': zero, 'Z_F': zero, 'shock_def_F': zero,
                   'shock_def_D': dShock, 'cb_buy_D': cb_path}
    ca = irf['ca_res_D']; gmF = irf['goods_mkt_F']
    # predicted CB hole: q*dcb_t - (1+rb)*q*dcb_{t-1}  (linearized at SS, dq second-order)
    cb_lag = np.concatenate([[0.0], cb_path[:-1]])
    hole = qss*cb_path - (1+rbss)*qss*cb_lag
    Tw = 100
    disc_D = betaD**np.arange(Tw); disc_F = betaF**np.arange(Tw)
    W_D = float((irf['U_D'][:Tw]*disc_D*100).sum())
    W_F = float((irf['U_F'][:Tw]*disc_F*100).sum())
    inj = float((hole[:Tw]*disc_D*100/Css).sum())   # discounted injection, % of SS C units
    summary[gam] = {
        'max_ca_res': float(np.max(np.abs(ca))),
        'max_goods_mkt_F': float(np.max(np.abs(gmF))),
        'max_hole_pred': float(np.max(np.abs(hole))),
        'max_ca_minus_hole_minus_basegap': None,
        'corr_ca_hole': float(np.corrcoef(ca[:200], hole[:200])[0,1]) if gam>0 else None,
        'W_D': W_D, 'W_F': W_F,
        'disc_injection_pctC': inj,
        'cb_peak': float(np.max(np.abs(cb_path))),
        'spread_peak': float(np.max(irf['spread_rb'][:100])),
        'U_D0': float(irf['U_D'][0]),
    }
    # subtract baseline (gamma=0) ca leak to isolate CB hole
    if gam == 0:
        ca0 = ca.copy()
    else:
        summary[gam]['max_ca_minus_hole_minus_basegap'] = float(np.max(np.abs(ca - ca0 - hole)))
    jlog(f'gamma_{gam}', summary[gam])

for gam in [2, 5, 10]:
    jlog(f'dW_D_gamma_{gam}', summary[gam]['W_D'] - summary[0]['W_D'])
    jlog(f'dW_F_gamma_{gam}', summary[gam]['W_F'] - summary[0]['W_F'])
    jlog(f'injection_gamma_{gam}', summary[gam]['disc_injection_pctC'])

with open(ROOT/'audit_artifacts'/'tpi_results.json','w') as fh:
    json.dump(RES, fh, indent=2, default=str)
print("TPI TEST COMPLETE")
