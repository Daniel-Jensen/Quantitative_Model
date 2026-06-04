import json, numpy as np
import sequence_jacobian as sj

# ── reproduce calibration cell (cell 2) ──
nb = json.load(open('model_v12.ipynb'))
ns = {}
exec(''.join(nb['cells'][2]['source']), {'np': np}, ns)
cal = ns['calibration_start']

# ── cell 3: seed bond holdings ──
_n_D, _n_F = cal['n_inter_D'], cal['n_inter_F']
_B_D, _B_F = cal['B_supply_D'], cal['B_supply_F']
b_F_D = cal['phi_bF_D_ss'] * _n_D / cal['q_b_F']
b_D_F = cal['phi_bD_F_ss'] * _n_F / cal['q_b_D']
cal.update({
    'b_F_D': b_F_D, 'b_D_F': b_D_F,
    'b_D_D': _B_D - b_D_F, 'b_F_F': _B_F - b_F_D,
    'b_F_D_anchor': b_F_D, 'b_D_F_anchor': b_D_F,
    'psi_bD_D': 0.5, 'psi_bF_F': 0.5,
})

from equations_D import (
    smart_steady_D, market_clearing_D, steady_auxilliary_D, banker_div_D,
    sdf_ss_D, government_ss_D, labor_ss_D, government_default_D,
    bond_price_ss_D, bond_return_D, ces_price_D, import_demand_D,
    deposit_return_D, hh_extended_D, income_D, make_grids_D,
)
from equations_F import (
    smart_steady_F, market_clearing_F, steady_auxilliary_F, banker_div_F,
    sdf_ss_F, government_ss_F, labor_ss_F, government_default_F,
    bond_price_ss_F, bond_return_F, ces_price_F, import_demand_F,
    deposit_return_F, hh_extended_F,
)
from equations_global import (
    trade_balance, domestic_bond_clearing, portfolio_level_anchors,
    portfolio_adj_cost, bond_yield, global_goods_mkt, external_account_D,
)

ha = sj.create_model([
    sdf_ss_D, government_default_D, bond_price_ss_D, bond_return_D,
    sdf_ss_F, government_default_F, bond_price_ss_F, bond_return_F,
    hh_extended_D, smart_steady_D, market_clearing_D, steady_auxilliary_D,
    banker_div_D, government_ss_D, labor_ss_D,
    hh_extended_F, smart_steady_F, market_clearing_F, steady_auxilliary_F,
    banker_div_F, government_ss_F, labor_ss_F,
    ces_price_D, import_demand_D, ces_price_F, import_demand_F,
    deposit_return_D, deposit_return_F, bond_yield,
    trade_balance, external_account_D, global_goods_mkt,
], name='MU HA SS')

unknowns_ss = {'beta_D': 0.9678948189, 'beta_F': 0.9654079110, 'p': 0.997170}
targets_ss  = ['deposit_mkt_D', 'deposit_mkt_F', 'ca_res_D']

print("Solving SS ...", flush=True)
ss = ha.solve_steady_state(cal, unknowns_ss, targets_ss, solver='broyden_custom')
print("SS solved.\n", flush=True)

# ── recover per-state tax from household hetinputs at SS ──
w, N, div, tau, lamb, P = (float(ss[k]) for k in
    ['w_D','N_D','div_D','tau_D','lamb_D','P_CES_D'])
intern = ss.internals['hh_D']
e = intern['e_grid_D']
D = intern['D']                      # distribution over (z, a)
Dz = D.sum(axis=tuple(range(1, D.ndim)))   # marginal income dist
Dz = Dz / Dz.sum()

y_pre = (w * N * e + div) / P
z_post = lamb * y_pre ** (1 - tau)
t_paid = y_pre - z_post
rate_state = t_paid / y_pre

Ey  = float((Dz * y_pre).sum())
Et  = float((Dz * t_paid).sum())
avg_rate = Et / Ey

print(f"beta_D={float(ss['beta_D']):.6f}  P_CES_D={P:.4f}  w_D={w:.4f}  div_D={div:.4f}")
print(f"mean pre-tax income  E[y_pre] = {Ey:.4f}")
print(f"mean tax paid        E[t]     = {Et:.4f}")
print(f"DIST-WEIGHTED AVG TAX RATE    = {avg_rate*100:.2f}%")
print(f"  per-state rate range: {rate_state.min()*100:.1f}% (lowest e) .. {rate_state.max()*100:.1f}% (highest e)")
print()

# ── fiscal aggregates ──
Y   = float(ss['Y_D'])
TAX = float(ss['TAX_D'])
G   = float(ss['G_D'])
print(f"Y_D (quarterly) = {Y:.4f}   (annual GDP = {4*Y:.2f})")
print(f"TAX_D (revenue) = {TAX:.4f}   ->  revenue/Y = {TAX/Y*100:.2f}%   revenue/annualGDP = {TAX/(4*Y)*100:.2f}%")
print(f"G_D  (spending) = {G:.4f}   ->  G/Y       = {G/Y*100:.2f}%       G/annualGDP       = {G/(4*Y)*100:.2f}%")
print(f"P_CES*TAX = {P*TAX:.4f}")

# ── explicit budget_residual_D at SS ──
b_gov   = float(ss['b_gov_D']); q_b = float(ss['q_b_D'])
def_rate= float(ss['def_rate_D']); rec = float(ss['recovery_rate_D'])
zeta    = float(ss['zeta_writeoff_D']); db = float(ss['delta_b_D'])
haircut = 1 - rec
surv    = 1 - zeta*def_rate*haircut
coupon  = db*(1 - def_rate*haircut)*b_gov
net_iss = q_b*(b_gov - surv*(1-db)*b_gov)
budget_res = coupon + G - P*TAX - net_iss
print(f"\nbudget_residual_D at SS = coupon + G - P*TAX - net_issuance = {budget_res:.3e}")
print(f"  coupon={coupon:.4f}  net_issuance={net_iss:.4f}")
print(f"\ndeposit_mkt_D residual = {float(ss['deposit_mkt_D']):.3e}")
print(f"ca_res_D residual      = {float(ss['ca_res_D']):.3e}")
