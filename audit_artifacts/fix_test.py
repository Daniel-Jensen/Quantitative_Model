"""Fix verification: patch W-1 (K-timing), W-2 (F-bank ToT revaluation), W-3 (FOC p-conversion)
in COPIES of the blocks (model files untouched) and test whether the untargeted
residuals ca_res_D / goods_mkt_F collapse to machine zero.

Variants:
  base      : as in notebook (rerun for IRF comparison baselines)  -- skipped, loaded from irfs.npz
  fix12     : W-1 (both countries) + W-2
  fix123    : W-1 + W-2 + W-3
Outputs -> audit_artifacts/fix_results.json, audit_artifacts/irfs_fixed.npz
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
def jlog(key, val):
    RES[key] = val
    print(f"[fix] {key}: {val}")

# ── reuse the exact SS pipeline from run_audit by exec'ing its SS section ──
src = open(ROOT/'audit_artifacts'/'run_audit.py').read()
ss_part = src.split("# ═══════════════ CELL 20")[0]
# drop the json dump line that needs RES['ss'] context but keep everything else
ns = {'__file__': str(ROOT/'audit_artifacts'/'run_audit.py')}
exec(compile(ss_part, 'ss_pipeline', 'exec'), ns)
ss_final = ns['ss_final']
cali_D = cali_F = ss_final
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
from equations_global import (trade_balance, domestic_bond_clearing,
    portfolio_level_anchors, divert_portfolio_adj, bond_yield,
    global_goods_mkt, external_account_D)

# ═══════════ PATCHED BLOCKS ═══════════
# W-1: production and mpk use K(-1)  (standard SSJ timing)
@simple
def labor_D_fx(N_D, Z_D, K_D, alpha_D):
    Y_D = Z_D * K_D(-1) ** alpha_D * N_D ** (1 - alpha_D)
    return Y_D

@simple
def capital_adj_D_fx(K_D, Q_D, I_D, Z_D, N_D, alpha_D, delta_D, gamma0_D, gamma1_D, ksi_D):
    iota_D        = I_D / K_D(-1)
    mpk_D         = alpha_D * Z_D * K_D(-1) ** (alpha_D - 1) * N_D ** (1 - alpha_D)
    rk_D          = (mpk_D + (1 - delta_D) * Q_D) / Q_D(-1) - 1
    q_res_D       = Q_D - 1 / (gamma0_D * (1 - ksi_D) * iota_D ** (-ksi_D))
    capital_res_D = K_D - (1 - delta_D) * K_D(-1) - (gamma0_D * iota_D ** (1 - ksi_D) + gamma1_D) * K_D(-1)
    return iota_D, mpk_D, rk_D, q_res_D, capital_res_D

@simple
def labor_F_fx(N_F, Z_F, K_F, alpha_F):
    Y_F = Z_F * K_F(-1) ** alpha_F * N_F ** (1 - alpha_F)
    return Y_F

@simple
def capital_adj_F_fx(K_F, Q_F, I_F, Z_F, N_F, alpha_F, delta_F, gamma0_F, gamma1_F, ksi_F):
    iota_F        = I_F / K_F(-1)
    mpk_F         = alpha_F * Z_F * K_F(-1) ** (alpha_F - 1) * N_F ** (1 - alpha_F)
    rk_F          = (mpk_F + (1 - delta_F) * Q_F) / Q_F(-1) - 1
    q_res_F       = Q_F - 1 / (gamma0_F * (1 - ksi_F) * iota_F ** (-ksi_F))
    capital_res_F = K_F - (1 - delta_F) * K_F(-1) - (gamma0_F * iota_F ** (1 - ksi_F) + gamma1_F) * K_F(-1)
    return iota_F, mpk_F, rk_F, q_res_F, capital_res_F

# W-2: F-bank realized bond returns converted to F-goods with p(-1)/p
@simple
def bank_return_F_fx(theta_F, rk_F, rdep_F, b_F_F, b_D_F, n_inter_F,
                     rb_actual_F, rb_actual_D, q_b_F, q_b_D, p):
    phi_bF_lag_F = q_b_F(-1) * b_F_F(-1) / (p(-1) * n_inter_F(-1))
    phi_bD_lag_F = q_b_D(-1) * b_D_F(-1) / (p(-1) * n_inter_F(-1))
    kappa_lag_F  = theta_F(-1) - phi_bF_lag_F - phi_bD_lag_F
    rb_F_fg = (1 + rb_actual_F) * p(-1) / p - 1
    rb_D_fg = (1 + rb_actual_D) * p(-1) / p - 1
    rn_F = (kappa_lag_F  * (rk_F    - rdep_F)
            + phi_bF_lag_F * (rb_F_fg - rdep_F)
            + phi_bD_lag_F * (rb_D_fg - rdep_F)
            + rdep_F)
    return rn_F

# W-3: F-bank own-bond FOC with p-conversion (mirrors dead domestic_bond_foc_F)
@simple
def divert_bond_foc_F_fx(rb_actual_F, rdep_F, b_F_F, n_inter_F, q_b_F,
                         phi_bF_F_ss, psi_bF_F, excess_return_bF_F_ss, tau_mp_F, p,
                         psi_spread_F, def_rate_F):
    phi_bF_F   = q_b_F * b_F_F / (p * n_inter_F)
    rb_F_fg_next = (1 + rb_actual_F(+1)) * p / p(+1) - 1
    req_spread = excess_return_bF_F_ss + psi_spread_F * def_rate_F(+1)
    rb_F_res   = (rb_F_fg_next - rdep_F(+1)) - req_spread \
                 - psi_bF_F * (phi_bF_F - phi_bF_F_ss) \
                 - tau_mp_F
    return rb_F_res

sys.setrecursionlimit(5000)
def build(blocks_D_labor, blocks_D_capadj, blocks_F_labor, blocks_F_capadj, bank_ret_F, foc_F):
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
    return sj.create_model([
        deposit_return_D, tax_rule_D, hh_extended_D, ghh_composite_D, sdf_D, sdf_banker_D,
        government_default_D, financial_solved_D, bond_return_D, bank_return_D,
        cap_adj_cost_inter_D, macro_pru_tax_D, intermediation_P2_D, intermediation_P3_D,
        k_balance_sheet_D, blocks_D_capadj, capital_producer_profit_D, budget_residual_D,
        blocks_D_labor, labor_market_D, labor_demand_D, banker_div_res_D, market_clearing_D, welfare_agg_D,
        deposit_return_F, tax_rule_F, hh_extended_F, ghh_composite_F, sdf_F, sdf_banker_F,
        government_default_F, financial_solved_F, bond_return_F, bank_ret_F,
        cap_adj_cost_inter_F, macro_pru_tax_F, intermediation_P2_F, intermediation_P3_F,
        k_balance_sheet_F, blocks_F_capadj, capital_producer_profit_F, budget_residual_F,
        blocks_F_labor, labor_market_F, labor_demand_F, banker_div_res_F, market_clearing_F, welfare_agg_F,
        ces_price_D, import_demand_D, ces_price_F, import_demand_F,
        trade_balance, external_account_D, domestic_bond_clearing, bond_yield,
        portfolio_level_anchors, divert_portfolio_adj, divert_bond_foc_D, foc_F,
        global_goods_mkt,
    ], name="fixed")

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
shocks = {'Z_D': {'Z_D': dZ,'Z_F': zero,'shock_def_D': zero,'shock_def_F': zero},
          'def_D': {'Z_D': zero,'Z_F': zero,'shock_def_D': dZ,'shock_def_F': zero}}

variants = {
    'fix12':  build(labor_D_fx, capital_adj_D_fx, labor_F_fx, capital_adj_F_fx, bank_return_F_fx, divert_bond_foc_F),
    'fix123': build(labor_D_fx, capital_adj_D_fx, labor_F_fx, capital_adj_F_fx, bank_return_F_fx, divert_bond_foc_F_fx),
}

base = np.load(ROOT/'audit_artifacts'/'irfs.npz')
save = {}
for vname, mdl in variants.items():
    print(f"Solving Jacobian for {vname} ...")
    G = mdl.solve_jacobian(ss_final, unknowns=unknowns_tp, targets=targets_tp, inputs=exogenous, T=T)
    for shname, sh in shocks.items():
        irf = G @ sh
        leaks = {}
        for v in ['ca_res_D','goods_mkt_F','goods_mkt_D','global_goods_res']:
            try:
                leaks[v] = float(np.max(np.abs(irf[v])))
            except (KeyError, TypeError):
                leaks[v] = 'missing'
        jlog(f'{vname}_{shname}_leaks', leaks)
        comp = {}
        for v in ['Y_D','Y_F','C_F','n_inter_F','p','spread_rb','q_b_F','rdep_D','U_F','NX_D','n_inter_D']:
            try:
                b = base[f'{shname}__{v}']; f_ = irf[v]
                t = min(len(b), len(f_))
                comp[v] = {'base_peak': float(b[np.argmax(np.abs(b[:100]))]),
                           'fix_peak': float(f_[np.argmax(np.abs(f_[:100]))]),
                           'maxdiff': float(np.max(np.abs(b[:t]-f_[:t])))}
            except (KeyError, TypeError) as e:
                comp[v] = repr(e)
        jlog(f'{vname}_{shname}_irfcomp', comp)
        for v in ['Y_D','Y_F','C_F','C_D','n_inter_F','n_inter_D','p','spread_rb','q_b_D','q_b_F',
                  'rdep_D','rdep_F','U_D','U_F','NX_D','ca_res_D','goods_mkt_F','def_rate_D','b_gov_D',
                  'K_D','I_D','N_D','b_F_D','b_D_F','theta_D','rb_D']:
            try:
                save[f'{vname}__{shname}__{v}'] = irf[v]
            except (KeyError, TypeError):
                pass

np.savez(ROOT/'audit_artifacts'/'irfs_fixed.npz', **save)
with open(ROOT/'audit_artifacts'/'fix_results.json','w') as fh:
    json.dump(RES, fh, indent=2, default=str)
print("FIX TEST COMPLETE")
