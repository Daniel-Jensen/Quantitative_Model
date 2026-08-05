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
    bond_return_D, bank_return_D, capital_fund_D, cap_adj_cost_inter_D, macro_pru_tax_D,
    intermediation_P2_D, intermediation_P3_D, k_balance_sheet_D,
    capital_adj_D, capital_producer_profit_D,
    labor_D, labor_market_D, labor_demand_D, banker_div_res_D,
    market_clearing_D, welfare_agg_D, ces_price_D, import_demand_D,
    divert_bond_foc_D,
)
from equations_F import (
    deposit_return_F, tax_rule_F, hh_extended_F, ghh_composite_F,
    sdf_F, sdf_banker_F, government_default_F,
    bond_return_F, bank_return_F, capital_fund_F, cap_adj_cost_inter_F, macro_pru_tax_F,
    intermediation_P2_F, intermediation_P3_F, k_balance_sheet_F,
    capital_adj_F, capital_producer_profit_F,
    labor_F, labor_market_F, labor_demand_F, banker_div_res_F,
    market_clearing_F, welfare_agg_F, ces_price_F, import_demand_F,
    divert_bond_foc_F,
)
from equations_global import (
    trade_balance, bond_yield,
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


# ── ECB balance sheet: capital-key conduit ────────────────────────────────────
# The CB holds cb_buy_D and passes its entire net cash flow through to the two
# treasuries by capital key: kappa_cb_F to F, (1-kappa_cb_F) to D. Each treasury
# finances its capital calls through its own fiscal rule at its own sovereign
# terms, so the funding (carry) leg of the position is endogenous. cb_flow_D is
# denominated in D-goods; the F remittance converts with /p like all bond-market
# cash flows in budget_residual_F.
@simple
def budget_residual_D_tpi(b_gov_D, G_D, TAX_D, q_b_D, def_rate_D, recovery_rate_D,
                           zeta_writeoff_D, P_CES_D, delta_b_D, writeoff_enabled_D,
                           cb_buy_D, kappa_cb_F):
    haircut_D      = 1.0 - recovery_rate_D
    haircut_mult_D = writeoff_enabled_D
    surv_cont_D    = 1.0 - zeta_writeoff_D * def_rate_D * haircut_D * haircut_mult_D
    coupon_D       = delta_b_D * (1.0 - def_rate_D * haircut_D * haircut_mult_D) * b_gov_D(-1)
    net_issuance_D = q_b_D * (b_gov_D - surv_cont_D * (1.0 - delta_b_D) * b_gov_D(-1))
    # CB net cash flow = (1+rb_actual_D)*q_b_D(-1)*cb_buy_D(-1) - q_b_D*cb_buy_D,
    # written out in coupon/survival form to match the budget's own accounting.
    cb_flow_D      = (delta_b_D * (1.0 - def_rate_D * haircut_D * haircut_mult_D) * cb_buy_D(-1)
                      + q_b_D * surv_cont_D * (1.0 - delta_b_D) * cb_buy_D(-1)
                      - q_b_D * cb_buy_D)
    rem_cb_D       = (1.0 - kappa_cb_F) * cb_flow_D
    b_gov_res_D    = coupon_D + G_D - P_CES_D * TAX_D - net_issuance_D - rem_cb_D
    return b_gov_res_D, rem_cb_D, cb_flow_D


@simple
def budget_residual_F_tpi(b_gov_F, G_F, TAX_F, q_b_F, def_rate_F, recovery_rate_F,
                           zeta_writeoff_F, p, P_CES_F, delta_b_F, writeoff_enabled_F,
                           cb_flow_D, kappa_cb_F):
    haircut_F      = 1.0 - recovery_rate_F
    haircut_mult_F = writeoff_enabled_F
    surv_cont_F    = 1.0 - zeta_writeoff_F * def_rate_F * haircut_F * haircut_mult_F
    coupon_F       = delta_b_F * (1.0 - def_rate_F * haircut_F * haircut_mult_F) * b_gov_F(-1)
    net_issuance_F = q_b_F * (b_gov_F - surv_cont_F * (1.0 - delta_b_F) * b_gov_F(-1))
    rem_cb_F       = kappa_cb_F * cb_flow_D / p
    b_gov_res_F    = (coupon_F - net_issuance_F) / p + G_F - P_CES_F * TAX_F - rem_cb_F
    return b_gov_res_F, rem_cb_F


@simple
def external_account_D_tpi(NX_D, q_b_D, q_b_F, b_F_D, b_D_F, rb_actual_F, rb_actual_D,
                           cb_buy_D, kappa_cb_F):
    # The F share of the CB book is an F claim on D: it enters D's external
    # account exactly like b_D_F. The D share stays domestic (like b_D_D).
    receipts_from_F_bonds = (1 + rb_actual_F) * q_b_F(-1) * b_F_D(-1)
    payments_on_D_bonds   = (1 + rb_actual_D) * q_b_D(-1) * (b_D_F(-1) + kappa_cb_F * cb_buy_D(-1))
    nfa_D = q_b_F * b_F_D - q_b_D * (b_D_F + kappa_cb_F * cb_buy_D)
    ca_res_D = (NX_D + receipts_from_F_bonds - payments_on_D_bonds - nfa_D)
    return nfa_D, ca_res_D


def tpi_overrides():
    """The four blocks the TPI layer swaps into the shared block list."""
    return {
        'budget_residual_D':     budget_residual_D_tpi,
        'budget_residual_F':     budget_residual_F_tpi,
        'external_account_D':    external_account_D_tpi,
        'domestic_bond_clearing': domestic_bond_clearing_tpi,
    }


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
    from full_model import build_block_list
    ha_full_tpi = sj.create_model(
        build_block_list(financial_solved_D, financial_solved_F,
                         overrides=tpi_overrides()),
        name="Full 2-Country MU HANK — TPI Extension",
    )

    ss_tpi = copy.deepcopy(ss_final)
    ss_tpi.toplevel['cb_buy_D'] = 0.0
    # inter-block CB flow: zero at SS (cb_buy_ss = 0); SSJ needs the symbol to
    # evaluate budget_residual_F_tpi's partial Jacobian
    ss_tpi.toplevel['cb_flow_D'] = 0.0
    _kap = ss_final.toplevel.get('kappa_cb_F')
    if _kap is None:
        from calibration import get_calibration
        _kap = get_calibration()['kappa_cb_F']
    kappa_cb_F = float(_kap)
    ss_tpi.toplevel['kappa_cb_F'] = kappa_cb_F

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
        _ca  = np.max(np.abs(irfs_tpi[g]['ca_res_D']))
        _gmF = np.max(np.abs(irfs_tpi[g]['goods_mkt_F']))
        print(f"peak spread = {_spread[:100].max()*100:+.3f} pp   "
              f"max|ca_res_D| = {_ca:.2e}   max|goods_mkt_F| = {_gmF:.2e}")
        if _ca > 1e-6 or _gmF > 1e-6:
            print(f"  WARNING: gamma={g}: external/goods residual > 1e-6 — CB conduit accounting leak")

    _err0 = np.max(np.abs(irfs_tpi[0]['spread_rb'][:50] - irfs_def_D['spread_rb'][:50]))
    print(f"Sanity check gamma=0 vs irfs_def_D: max |err| = {_err0:.2e}  (expect < 1e-8)")
    print(f"Sign check (gamma=0): n_inter_D[0] = {irfs_tpi[0]['n_inter_D'][0]:+.3e}, "
          f"Y_D[0] = {irfs_tpi[0]['Y_D'][0]:+.3e}  (both must be negative)")

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

    # ── ECB balance-sheet P&L: two-leg (carry + credit) decomposition ─────────
    # Off-path accounting per FRAMING_HANDOFF §2: products of first-order paths
    # (def_rate·cb_buy, spread·cb_buy) are second-order objects — hand-computed
    # here, never read off the linear DAG, which would drop them and
    # mechanically show the CB profiting with no offsetting expected loss.
    # Discounting at beta_F (creditor-side rate, handoff §7 A5).
    q_b_D_ss     = float(ss_final['q_b_D'])
    q_b_F_ss     = float(ss_final['q_b_F'])
    delta_b_D_v  = float(ss_final['delta_b_D'])
    delta_b_F_v  = float(ss_final['delta_b_F'])
    EL_price_D_v = float(ss_final['EL_price_D'])
    Y_D_ss       = float(ss_final['Y_D'])
    rb_D_ss      = delta_b_D_v * (1.0 / q_b_D_ss - 1.0)
    rb_F_ss      = delta_b_F_v * (1.0 / q_b_F_ss - 1.0)
    spread_ss    = rb_D_ss - rb_F_ss

    def cb_pnl(irf, T_pnl=100):
        """PV decomposition of the CB's D-bond position (D-goods units).

        el_pv       — credit leg: priced expected loss borne on holdings
        prem_pv     — default-risk premium income (spread deviation × position)
        carry_ss_pv — steady carry: SS yield differential × position (Lehment
                      funding leg: each treasury funds at its own bond rate)
        mtm_pv      — on-path revaluation of carried holdings
        """
        disc  = beta_F ** np.arange(T_pnl)
        cb    = irf['cb_buy_D'][:T_pnl]
        cb_l  = np.concatenate([[0.0], cb[:-1]])
        dq    = irf['q_b_D'][:T_pnl]
        dq_l  = np.concatenate([[0.0], dq[:-1]])
        defr  = irf['def_rate_D'][:T_pnl]
        dspr  = irf['spread_rb'][:T_pnl]
        purchases = cb - (1.0 - delta_b_D_v) * cb_l
        return {
            'peak_exposure': float(np.max(q_b_D_ss * cb)),
            'purchases_pv':  float((disc * q_b_D_ss * purchases).sum()),
            'el_pv':         float((disc * EL_price_D_v * defr * q_b_D_ss * cb).sum()),
            'prem_pv':       float((disc * dspr * q_b_D_ss * cb_l).sum()),
            'carry_ss_pv':   float((disc * spread_ss * q_b_D_ss * cb_l).sum()),
            'mtm_pv':        float((disc * (1.0 - delta_b_D_v) * cb_l * (dq - dq_l)).sum()),
        }

    pnl_by_gamma = {g: cb_pnl(irfs_tpi[g]) for g in gamma_values}
    pct = lambda v: 100.0 * v / Y_D_ss
    print(f"\nECB balance-sheet P&L (PV at beta_F over 100q; % of quarterly SS Y_D)")
    print(f"Capital-key conduit: F share {kappa_cb_F:.3f}, D share {1-kappa_cb_F:.3f} "
          f"(two-country renormalised euro-area key)")
    print(f"{'γ':>4} {'peak expos.':>12} {'purch. PV':>10} {'EL PV':>10} {'prem PV':>10} "
          f"{'SS-carry PV':>12} {'MTM PV':>10} {'loading':>8}")
    print("─" * 82)
    for g in gamma_values:
        d = pnl_by_gamma[g]
        loading = d['prem_pv'] / d['el_pv'] if d['el_pv'] > 1e-16 else float('nan')
        print(f"{g:>4} {pct(d['peak_exposure']):>11.3f}% {pct(d['purchases_pv']):>9.3f}% "
              f"{pct(d['el_pv']):>9.4f}% {pct(d['prem_pv']):>9.4f}% "
              f"{pct(d['carry_ss_pv']):>11.4f}% {pct(d['mtm_pv']):>9.4f}% {loading:>8.2f}")
    print("─" * 82)
    print("loading = prem PV / EL PV  (theory: ≈ 1 + psi_spread/EL_price at small γ, "
          "declining in γ — the self-extinguishing premium)")
    for g in gamma_values[1:]:
        d = pnl_by_gamma[g]
        print(f"  γ={g:<2}: F bears EL PV = {pct(kappa_cb_F*d['el_pv']):.4f}% Y_D, "
              f"receives prem PV = {pct(kappa_cb_F*d['prem_pv']):.4f}% "
              f"(memo at full euro-area key 26.1%: EL {pct(0.261*d['el_pv']):.4f}%)")

    # ── Effectiveness curve over fine gamma grid ──────────────────────────────
    def _peak_spread(irf):
        sp = irf['spread_rb'] if 'spread_rb' in irf \
             else irf['rb_actual_D'] - irf['rb_actual_F']
        return sp[:100].max()

    peak_no_tpi  = _peak_spread(irfs_tpi[0])
    gammas_fine  = np.concatenate([np.linspace(0, 5, 25), np.linspace(5, 30, 26)[1:]])
    peak_arr     = np.empty(len(gammas_fine))
    cost_arr     = np.empty(len(gammas_fine))
    el_pv_arr    = np.empty(len(gammas_fine))
    prem_pv_arr  = np.empty(len(gammas_fine))
    mtm_pv_arr   = np.empty(len(gammas_fine))
    expos_arr    = np.empty(len(gammas_fine))
    loading_arr  = np.full(len(gammas_fine), np.nan)
    for i, g in enumerate(gammas_fine):
        irf_g = irfs_tpi[0] if g == 0 else compute_tpi_irfs(G_tpi, dShock_def_D, g, T)
        peak_arr[i] = _peak_spread(irf_g)
        cost_arr[i] = (irf_g['cb_buy_D'][:100] * q_b_D_ss).sum()
        d = cb_pnl(irf_g)
        el_pv_arr[i]   = d['el_pv']
        prem_pv_arr[i] = d['prem_pv']
        mtm_pv_arr[i]  = d['mtm_pv']
        expos_arr[i]   = d['peak_exposure']
        if d['el_pv'] > 1e-16:
            loading_arr[i] = d['prem_pv'] / d['el_pv']
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
        # ECB balance-sheet block (capital-key conduit)
        'kappa_cb_F':   kappa_cb_F,
        'pnl_by_gamma': pnl_by_gamma,
        'el_pv_arr':    el_pv_arr,
        'prem_pv_arr':  prem_pv_arr,
        'mtm_pv_arr':   mtm_pv_arr,
        'expos_arr':    expos_arr,
        'loading_arr':  loading_arr,
        'Y_D_ss':       Y_D_ss,
        'spread_ss':    spread_ss,
    }
