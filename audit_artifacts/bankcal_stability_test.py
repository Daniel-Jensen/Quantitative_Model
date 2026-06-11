"""Post-T-2 stability: find minimal phi_lamb that restores stationarity.
T-2 removed the accidental stabilizer (state-contingent rdep windfall), exposing
the debt->spread->debt spiral; phi_lamb=0.02 is now too weak.
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
def jlog(k, v):
    RES[k] = v
    print(f"[phi] {k}: {v}", flush=True)

src = open(ROOT/'audit_artifacts'/'run_audit.py').read()
ss_part = src.split("# ═══════════════ CELL 20")[0]
ns = {'__file__': str(ROOT/'audit_artifacts'/'run_audit.py')}
exec(compile(ss_part, 'ss_pipeline', 'exec'), ns)
ss_final = ns['ss_final']
cali = ss_final

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
from equations_global import (trade_balance, domestic_bond_clearing,
    portfolio_level_anchors, divert_portfolio_adj, bond_yield,
    global_goods_mkt, external_account_D)

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
], name="full")

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
rho = 0.8
dZ = 0.01*rho**np.arange(T)
zero = np.zeros(T)

best = None
# bank-cal-style low amplification: psi_lambda_B=0, def_scale=0.05, f=0.03, delta_b longer
ss_final.toplevel['psi_lambda_B_D']=0.0; ss_final.toplevel['psi_lambda_B_F']=0.0
ss_final.toplevel['def_scale_D']=0.05;   ss_final.toplevel['def_scale_F']=0.05
ss_final.toplevel['psi_spread_D']=0.0;   ss_final.toplevel['psi_spread_F']=0.0   # follows from psi_lambda_B=0
for phl in [0.03, 0.05, 0.08, 0.12]:
    ss_try = copy.deepcopy(ss_final)
    ss_try.toplevel['phi_lamb_D'] = phl
    ss_try.toplevel['phi_lamb_F'] = phl
    print(f"--- phi_lamb = {phl} : solving Jacobian ---", flush=True)
    try:
        G = ha_full.solve_jacobian(ss_try, unknowns=unknowns_tp, targets=targets_tp,
                                   inputs=exogenous, T=T)
        irf_d = G @ {'Z_D': zero, 'Z_F': zero, 'shock_def_D': dZ, 'shock_def_F': zero}
        irf_z = G @ {'Z_D': dZ, 'Z_F': zero, 'shock_def_D': zero, 'shock_def_F': zero}
        out = {
            'b_gov_D_t499_def': float(irf_d['b_gov_D'][499]),
            'b_gov_D_t499_Z': float(irf_z['b_gov_D'][499]),
            'Y_D_t499_def': float(irf_d['Y_D'][499]),
            'max_ca_res_def': float(np.max(np.abs(irf_d['ca_res_D']))),
            'max_goodsF_def': float(np.max(np.abs(irf_d['goods_mkt_F']))),
            'max_ca_res_Z': float(np.max(np.abs(irf_z['ca_res_D']))),
            'max_goodsF_Z': float(np.max(np.abs(irf_z['goods_mkt_F']))),
            'Y_D0_def': float(irf_d['Y_D'][0]),
            'n_inter_D0_def': float(irf_d['n_inter_D'][0]),
            'n_inter_D_peak_def': float(irf_d['n_inter_D'][np.argmax(np.abs(irf_d['n_inter_D'][:60]))]),
            'spread0_def': float(irf_d['spread_rb'][0]),
            'U_D0_def': float(irf_d['U_D'][0]),
            'n_inter_D0_Z': float(irf_z['n_inter_D'][0]),
            'rdep_D_peak_Z': float(irf_z['rdep_D'][np.argmax(np.abs(irf_z['rdep_D'][:60]))]),
        }
        stable = abs(out['b_gov_D_t499_def']) < 1e-3 and abs(out['b_gov_D_t499_Z']) < 1e-3
        out['stable'] = bool(stable)
        jlog(f'phi_lamb_{phl}', out)
        if stable and best is None:
            best = phl
            np.savez(ROOT/'audit_artifacts'/'irfs_stable.npz',
                     **{f'def_D__{k}': v for k, v in irf_d.toplevel.items() if np.ndim(v)==1},
                     **{f'Z_D__{k}': v for k, v in irf_z.toplevel.items() if np.ndim(v)==1})
    except Exception as e:
        jlog(f'phi_lamb_{phl}', f'FAILED: {e!r}')

jlog('minimal_stable_phi_lamb', best)
with open(ROOT/'audit_artifacts'/'philamb_results.json','w') as fh:
    json.dump(RES, fh, indent=2, default=str)
print("PHILAMB TEST COMPLETE")
