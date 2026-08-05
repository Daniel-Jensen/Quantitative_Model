import sys
import copy
import numpy as np
import sequence_jacobian as sj
from sequence_jacobian import simple, combine

from equations_D import (
    capital_adj_D, labor_D, labor_market_D, labor_demand_D,
    intermediation_IC_D, bank_return_D, capital_fund_D, intermediation_P1_D,
    k_balance_sheet_D, cap_adj_cost_inter_D, macro_pru_tax_D,
    intermediation_P2_D, banker_div_res_D, intermediation_P3_D,
    government_default_D, divert_bond_foc_D,
    tax_rule_D, capital_producer_profit_D, budget_residual_D,
    ces_price_D, import_demand_D, deposit_return_D,
    bond_return_D, sdf_D, sdf_banker_ss_D, sdf_banker_D, ghh_composite_D,
    welfare_agg_D, market_clearing_D, hh_extended_D,
)
from equations_F import (
    capital_adj_F, labor_F, labor_market_F, labor_demand_F,
    intermediation_IC_F, bank_return_F, capital_fund_F, intermediation_P1_F,
    k_balance_sheet_F, cap_adj_cost_inter_F, macro_pru_tax_F,
    intermediation_P2_F, banker_div_res_F, intermediation_P3_F,
    government_default_F, divert_bond_foc_F,
    tax_rule_F, capital_producer_profit_F, budget_residual_F,
    ces_price_F, import_demand_F, deposit_return_F,
    bond_return_F, sdf_F, sdf_banker_ss_F, sdf_banker_F, ghh_composite_F,
    welfare_agg_F, market_clearing_F, hh_extended_F,
)
from equations_global import (
    trade_balance, domestic_bond_clearing,
    portfolio_level_anchors, divert_portfolio_adj, bond_yield,
    global_goods_mkt, external_account_D,
)


def build_block_list(financial_solved_D, financial_solved_F, *,
                     hh_D=None, hh_F=None, overrides=None):
    """The single definition of the model's block list.

    Every consumer (full_model, tpi, diagnostics/regimes) calls this. A second
    copy of the list is how the retired audit_artifacts/ harness drifted into
    silently testing a different model — see CLAUDE.md.

    financial_solved_D/F : the runtime-constructed GK solved blocks.
    hh_D/hh_F            : optionally REPLACE the household blocks with versions
                           carrying extra hetoutputs (experiments/e4_distribution
                           adds per-decile consumption).
    overrides            : {name: block} used by the TPI layer to swap in its
                           _tpi variants without keeping a second list.
    """
    o = overrides or {}
    hh_D = hh_extended_D if hh_D is None else hh_D
    hh_F = hh_extended_F if hh_F is None else hh_F

    def pick(name, default):
        return o.get(name, default)

    return [
        # Country D
        deposit_return_D, tax_rule_D, hh_D, ghh_composite_D,
        sdf_D, sdf_banker_D, government_default_D, financial_solved_D,
        bond_return_D, bank_return_D, capital_fund_D, cap_adj_cost_inter_D, macro_pru_tax_D,
        intermediation_P2_D, intermediation_P3_D, k_balance_sheet_D,
        capital_adj_D, capital_producer_profit_D,
        pick('budget_residual_D', budget_residual_D),
        labor_D, labor_market_D, labor_demand_D, banker_div_res_D,
        market_clearing_D, welfare_agg_D,
        # Country F
        deposit_return_F, tax_rule_F, hh_F, ghh_composite_F,
        sdf_F, sdf_banker_F, government_default_F, financial_solved_F,
        bond_return_F, bank_return_F, capital_fund_F, cap_adj_cost_inter_F, macro_pru_tax_F,
        intermediation_P2_F, intermediation_P3_F, k_balance_sheet_F,
        capital_adj_F, capital_producer_profit_F,
        pick('budget_residual_F', budget_residual_F),
        labor_F, labor_market_F, labor_demand_F, banker_div_res_F,
        market_clearing_F, welfare_agg_F,
        # Global
        ces_price_D, import_demand_D, ces_price_F, import_demand_F,
        trade_balance,
        pick('external_account_D', external_account_D),
        pick('domestic_bond_clearing', domestic_bond_clearing),
        bond_yield, portfolio_level_anchors, divert_portfolio_adj,
        divert_bond_foc_D, divert_bond_foc_F, global_goods_mkt,
    ]


def build_and_solve(ss_results):
    sys.setrecursionlimit(5000)

    ss_final    = ss_results['ss_final']
    cali_D      = ss_results['cali_D']
    cali_F      = ss_results['cali_F']
    calibration_start = ss_results['calibration_start']

    # ── Inner solved blocks for GK Bellman + IC ───────────────────────────────
    financial_solved_D = combine([
        intermediation_P1_D, intermediation_IC_D,
    ]).solved(
        unknowns={'nu_K_D':  float(cali_D['nu_K_D']),
                  'nu_bD_D': float(cali_D['nu_bD_D']),
                  'nu_bF_D': float(cali_D['nu_bF_D']),
                  'eta_D':   float(cali_D['eta_D']),
                  'theta_D': float(cali_D['theta_D'])},
        targets=['nu_K_res_D', 'nu_bD_res_D', 'nu_bF_res_D', 'eta_res_D', 'ic_res_D'],
        solver='broyden_custom'
    )
    financial_solved_F = combine([
        intermediation_P1_F, intermediation_IC_F,
    ]).solved(
        unknowns={'nu_K_F':  float(cali_F['nu_K_F']),
                  'nu_bF_F': float(cali_F['nu_bF_F']),
                  'nu_bD_F': float(cali_F['nu_bD_F']),
                  'eta_F':   float(cali_F['eta_F']),
                  'theta_F': float(cali_F['theta_F'])},
        targets=['nu_K_res_F', 'nu_bF_res_F', 'nu_bD_res_F', 'eta_res_F', 'ic_res_F'],
        solver='broyden_custom'
    )

    # ── Full dynamic model ────────────────────────────────────────────────────
    ha_full = sj.create_model(
        build_block_list(financial_solved_D, financial_solved_F),
        name="Full 2-Country MU HANK — GHH Preferences, Flex Price & Wage, No CB",
    )

    # ── 23×23 system ──────────────────────────────────────────────────────────
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
        'goods_mkt_D', 'rb_D_res', 'rb_F_res',
    ]
    T = 500

    # ── Market-value fiscal-rule SS reference (mv_rule=1 path) ────────────────
    # mv_gov_ss = q_b_ss · b_gov_ss so the market-value debt gap is zero at SS.
    # Set exactly from the solved SS; harmless when mv_rule=0 (term is ×0).
    ss_final.toplevel['mv_gov_ss_D'] = float(ss_final['q_b_D']) * float(ss_final['b_gov_ss_D'])
    ss_final.toplevel['mv_gov_ss_F'] = float(ss_final['q_b_F']) * float(ss_final['b_gov_ss_F'])

    # ── Jacobian ──────────────────────────────────────────────────────────────
    exogenous = ['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F']
    print(f"Computing Jacobian G (T={T}, {len(exogenous)} exogenous inputs)...")
    G = ha_full.solve_jacobian(ss_final, unknowns=unknowns_tp, targets=targets_tp,
                               inputs=exogenous, T=T)
    print("G computed successfully.")

    # ── Shocks ────────────────────────────────────────────────────────────────
    rho_Z_D      = 0.8
    dZ_D         = 0.01 * rho_Z_D ** np.arange(T)
    rho_def_D    = 0.8
    dShock_def_D = 0.01 * rho_def_D ** np.arange(T)

    irfs_Z_D = G @ {
        'Z_D': dZ_D, 'Z_F': np.zeros(T),
        'shock_def_D': np.zeros(T), 'shock_def_F': np.zeros(T)
    }
    irfs_def_D = G @ {
        'Z_D': np.zeros(T), 'Z_F': np.zeros(T),
        'shock_def_D': dShock_def_D, 'shock_def_F': np.zeros(T)
    }

    # ── Stability check ───────────────────────────────────────────────────────
    print("\n=== Stability check: debt level at t=499 (should be near 0) ===")
    print(f"  irfs_Z_D  ['b_gov_D'][499] = {irfs_Z_D['b_gov_D'][499]:.6f}")
    print(f"  irfs_def_D['b_gov_D'][499] = {irfs_def_D['b_gov_D'][499]:.6f}")
    phi_lamb = calibration_start['phi_lamb_D']
    print(f"  ρ_b (partial-eq.) = {round((0.953 * 0.95 + 0.05 - phi_lamb) / 0.953, 4)}"
          "  [target < 0.95]")
    # SSJ returns LEVEL deviations. Y_D_ss ~ 1 so Y_D*100 happens to read as a
    # percent, but n_inter_D_ss is not 1 (2.138 under BANK_SCOPE="broad", 3.0
    # pre-EBA, 0.408 under CT1), so the raw x100 was NOT a percentage and was not
    # comparable across calibrations — the "-7.227%" widely quoted in the docs is
    # the level deviation; the true impact is -3.38% of SS net worth. Both are
    # printed now: the level for continuity with the historical logs, the
    # percent-of-SS as the number to quote.
    n0, Y0 = irfs_def_D['n_inter_D'][0], irfs_def_D['Y_D'][0]
    n_ss, Y_ss = float(ss_final['n_inter_D']), float(ss_final['Y_D'])
    print(f"  n_inter_D[0] on default shock = {n0/n_ss*100:+.4f}% of SS"
          f"  (level dev {n0*100:+.4f}; negative = doom loop correct sign)")
    print(f"  Y_D[0]       on default shock = {Y0/Y_ss*100:+.4f}% of SS"
          f"  (level dev {Y0*100:+.4f}; negative = correct sign)")

    return {
        'ha_full':           ha_full,
        'financial_solved_D': financial_solved_D,
        'financial_solved_F': financial_solved_F,
        'G':                 G,
        'ss_final':          ss_final,
        'calibration_start': calibration_start,
        'unknowns_tp':       unknowns_tp,
        'targets_tp':        targets_tp,
        'T':                 T,
        'dZ_D':              dZ_D,
        'dShock_def_D':      dShock_def_D,
        'irfs_Z_D':          irfs_Z_D,
        'irfs_def_D':        irfs_def_D,
    }

