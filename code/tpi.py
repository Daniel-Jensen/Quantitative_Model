"""
TPI (Transmission Protection Instrument) experiment.

Builds the TPI-extended model, computes closed-loop IRFs for
gamma_values = [0, 2, 5, 10], and pre-computes all welfare and
effectiveness statistics needed by tpi_plots.py.
"""
import copy
import numpy as np
import sequence_jacobian as sj
from sequence_jacobian import simple

from equations_D import (
    deposit_return_D, tax_rule_D, hh_extended_D, ghh_composite_D,
    sdf_D, sdf_banker_D, government_default_D,
    bond_return_D, bank_return_D, cap_adj_cost_inter_D, macro_pru_tax_D,
    intermediation_P2_D, intermediation_P3_D, k_balance_sheet_D,
    capital_adj_D, capital_producer_profit_D,
    labor_D, labor_market_D, labor_demand_D, banker_div_res_D,
    market_clearing_D, welfare_agg_D, ces_price_D, import_demand_D,
    divert_bond_foc_D,
)
from equations_F import (
    deposit_return_F, tax_rule_F, hh_extended_F, ghh_composite_F,
    sdf_F, sdf_banker_F, government_default_F,
    bond_return_F, bank_return_F, cap_adj_cost_inter_F, macro_pru_tax_F,
    intermediation_P2_F, intermediation_P3_F, k_balance_sheet_F,
    capital_adj_F, capital_producer_profit_F, budget_residual_F,
    labor_F, labor_market_F, labor_demand_F, banker_div_res_F,
    market_clearing_F, welfare_agg_F, ces_price_F, import_demand_F,
    divert_bond_foc_F,
)
from equations_global import (
    trade_balance, external_account_D, bond_yield,
    portfolio_level_anchors, divert_portfolio_adj, global_goods_mkt,
)

BLUE       = '#002147'
RED        = '#8C1515'
BLUE_MUTED = '#4a6f8a'
RED_MUTED  = '#c0624a'


# ── TPI-1: CB bond clearing + budget constraint (audit fix) ──────────────────
@simple
def domestic_bond_clearing_tpi(b_gov_D, b_gov_F, b_D_F, b_F_D, cb_buy_D):
    b_D_D = b_gov_D - b_D_F - cb_buy_D
    b_F_F = b_gov_F - b_F_D
    return b_D_D, b_F_F


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


def compute_tpi_irfs(G_tpi, shock_def, gamma_tpi, T):
    _has_spread = 'spread_rb' in G_tpi.outputs
    if _has_spread:
        A_def = np.array(G_tpi['spread_rb']['shock_def_D'])
        A_cb  = np.array(G_tpi['spread_rb']['cb_buy_D'])
    else:
        A_def = (np.array(G_tpi['rb_actual_D']['shock_def_D'])
                 - np.array(G_tpi['rb_actual_F']['shock_def_D']))
        A_cb  = (np.array(G_tpi['rb_actual_D']['cb_buy_D'])
                 - np.array(G_tpi['rb_actual_F']['cb_buy_D']))

    I_T = np.eye(T)
    system_matrix = I_T - gamma_tpi * A_cb
    cond = np.linalg.cond(system_matrix)
    if cond > 1e10:
        print(f"  WARNING: system matrix cond = {cond:.2e} for gamma={gamma_tpi:.1f}")

    spread_cl   = np.linalg.solve(system_matrix, A_def @ shock_def)
    cb_buy_path = gamma_tpi * spread_cl

    irfs = G_tpi @ {
        'Z_D': np.zeros(T), 'Z_F': np.zeros(T),
        'shock_def_D': shock_def, 'shock_def_F': np.zeros(T),
        'cb_buy_D': cb_buy_path,
    }
    irfs['cb_buy_D'] = cb_buy_path
    return irfs


def run_tpi(model_results):
    ha_full            = model_results['ha_full']
    financial_solved_D = model_results['financial_solved_D']
    financial_solved_F = model_results['financial_solved_F']
    ss_final           = model_results['ss_final']
    unknowns_tp        = model_results['unknowns_tp']
    targets_tp         = model_results['targets_tp']
    T                  = model_results['T']
    dShock_def_D       = model_results['dShock_def_D']
    irfs_def_D         = model_results['irfs_def_D']

    # ── Build TPI model ───────────────────────────────────────────────────────
    ha_full_tpi = sj.create_model([
        deposit_return_D, tax_rule_D, hh_extended_D, ghh_composite_D,
        sdf_D, sdf_banker_D, government_default_D, financial_solved_D,
        bond_return_D, bank_return_D, cap_adj_cost_inter_D, macro_pru_tax_D,
        intermediation_P2_D, intermediation_P3_D, k_balance_sheet_D,
        capital_adj_D, capital_producer_profit_D, budget_residual_D_tpi,
        labor_D, labor_market_D, labor_demand_D, banker_div_res_D,
        market_clearing_D, welfare_agg_D,
        deposit_return_F, tax_rule_F, hh_extended_F, ghh_composite_F,
        sdf_F, sdf_banker_F, government_default_F, financial_solved_F,
        bond_return_F, bank_return_F, cap_adj_cost_inter_F, macro_pru_tax_F,
        intermediation_P2_F, intermediation_P3_F, k_balance_sheet_F,
        capital_adj_F, capital_producer_profit_F, budget_residual_F,
        labor_F, labor_market_F, labor_demand_F, banker_div_res_F,
        market_clearing_F, welfare_agg_F,
        ces_price_D, import_demand_D, ces_price_F, import_demand_F,
        trade_balance, external_account_D, domestic_bond_clearing_tpi,
        bond_yield, portfolio_level_anchors, divert_portfolio_adj,
        divert_bond_foc_D, divert_bond_foc_F, global_goods_mkt,
    ], name="Full 2-Country MU HANK — TPI Extension")

    ss_tpi = copy.deepcopy(ss_final)
    ss_tpi.toplevel['cb_buy_D'] = 0.0

    # ── Jacobian ──────────────────────────────────────────────────────────────
    exogenous_tpi = ['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F', 'cb_buy_D']
    print(f"Computing G_tpi (T={T}, {len(exogenous_tpi)} exogenous inputs)...")
    G_tpi = ha_full_tpi.solve_jacobian(
        ss_tpi, unknowns=unknowns_tp, targets=targets_tp,
        inputs=exogenous_tpi, T=T,
    )
    print("G_tpi computed.")

    _chk = G_tpi @ {
        'Z_D': np.zeros(T), 'Z_F': np.zeros(T),
        'shock_def_D': dShock_def_D, 'shock_def_F': np.zeros(T),
        'cb_buy_D': np.zeros(T),
    }
    _err = np.max(np.abs(_chk['spread_rb'][:50] - irfs_def_D['spread_rb'][:50]))
    print(f"Sanity check G_tpi[cb=0] vs baseline G: max |err| = {_err:.2e}  (expect < 1e-8)")

    # ── Closed-loop IRFs ──────────────────────────────────────────────────────
    gamma_values = [0, 2, 5, 10]
    gamma_labels = ['γ = 0 (No TPI)', 'γ = 2  (Weak)', 'γ = 5  (Medium)', 'γ = 10 (Strong)']
    TPI_COLORS  = [BLUE, '#1a6e3a', '#c87941', RED]
    TPI_LSTYLES = ['-', '-', '--', '-.']
    TPI_MARKERS = ['', 'o', '', 's']

    irfs_tpi = {}
    for g in gamma_values:
        print(f"  gamma = {g:2d} ...", end='  ', flush=True)
        if g == 0:
            _s = {'Z_D': np.zeros(T), 'Z_F': np.zeros(T),
                  'shock_def_D': dShock_def_D, 'shock_def_F': np.zeros(T),
                  'cb_buy_D': np.zeros(T)}
            irfs_tpi[g] = G_tpi @ _s
            irfs_tpi[g]['cb_buy_D'] = np.zeros(T)
        else:
            irfs_tpi[g] = compute_tpi_irfs(G_tpi, dShock_def_D, g, T)
        _spread = irfs_tpi[g]['spread_rb'] if 'spread_rb' in irfs_tpi[g] \
                  else irfs_tpi[g]['rb_actual_D'] - irfs_tpi[g]['rb_actual_F']
        print(f"peak spread = {_spread[:100].max()*100:+.3f} pp")

    _err0 = np.max(np.abs(irfs_tpi[0]['spread_rb'][:50] - irfs_def_D['spread_rb'][:50]))
    print(f"Sanity check gamma=0 vs irfs_def_D: max |err| = {_err0:.2e}  (expect < 1e-8)")

    # ── Welfare gains ─────────────────────────────────────────────────────────
    T_disc = 100
    beta_D = float(ss_final['beta_D']); beta_F = float(ss_final['beta_F'])
    disc_D = beta_D ** np.arange(T_disc); disc_F = beta_F ** np.arange(T_disc)
    W_D  = np.array([(irfs_tpi[g]['U_D'][:T_disc] * disc_D * 100).sum() for g in gamma_values])
    W_F  = np.array([(irfs_tpi[g]['U_F'][:T_disc] * disc_F * 100).sum() for g in gamma_values])
    dW_D = W_D - W_D[0]
    dW_F = W_F - W_F[0]

    print(f"\n{'γ':>5}  {'W_D':>10}  {'W_F':>10}  {'ΔW_D vs γ=0':>13}  {'ΔW_F vs γ=0':>13}")
    print("─" * 60)
    for i, g in enumerate(gamma_values):
        print(f"{g:>5}  {W_D[i]:>10.4f}  {W_F[i]:>10.4f}  {dW_D[i]:>+13.4f}  {dW_F[i]:>+13.4f}")
    print("─" * 60)
    print("(Units: % of quarterly SS consumption, discounted over 100 quarters)")

    # ── Effectiveness curve over fine gamma grid ──────────────────────────────
    def _peak_spread(irf):
        sp = irf['spread_rb'] if 'spread_rb' in irf \
             else irf['rb_actual_D'] - irf['rb_actual_F']
        return sp[:100].max()

    peak_no_tpi  = _peak_spread(irfs_tpi[0])
    gammas_fine  = np.concatenate([np.linspace(0, 5, 25), np.linspace(5, 30, 26)[1:]])
    q_b_D_ss     = float(ss_final['q_b_D'])
    peak_arr     = np.empty(len(gammas_fine))
    cost_arr     = np.empty(len(gammas_fine))
    for i, g in enumerate(gammas_fine):
        irf_g = irfs_tpi[0] if g == 0 else compute_tpi_irfs(G_tpi, dShock_def_D, g, T)
        peak_arr[i] = _peak_spread(irf_g)
        cost_arr[i] = (irf_g['cb_buy_D'][:100] * q_b_D_ss).sum()
    frac_closed = np.clip(100.0 * (1.0 - peak_arr / peak_no_tpi), 0, 100)

    return {
        'irfs_tpi':    irfs_tpi,
        'gamma_values': gamma_values,
        'gamma_labels': gamma_labels,
        'TPI_COLORS':  TPI_COLORS,
        'TPI_LSTYLES': TPI_LSTYLES,
        'TPI_MARKERS': TPI_MARKERS,
        'dW_D':        dW_D,
        'dW_F':        dW_F,
        'W_D':         W_D,
        'W_F':         W_F,
        'ss_final':    ss_final,
        'T':           T,
        'dShock_def_D': dShock_def_D,
        'G_tpi':       G_tpi,
        'peak_no_tpi': peak_no_tpi,
        'gammas_fine': gammas_fine,
        'peak_arr':    peak_arr,
        'cost_arr':    cost_arr,
        'frac_closed': frac_closed,
        'q_b_D_ss':    q_b_D_ss,
    }
