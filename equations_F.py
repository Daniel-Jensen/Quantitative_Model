import numpy as np
import scipy.linalg
import sequence_jacobian as sj
from sequence_jacobian import simple
from sequence_jacobian import grids
from pathlib import Path

try:
    BASE_DIR_F = Path(__file__).resolve().parent
except NameError:
    BASE_DIR_F = Path.cwd()

DATA_DIR_F = BASE_DIR_F / "Discretisation" / "Outputs"

# ── Household ─────────────────────────────────────────────────────────────────

def hh_init_F(dep_F_grid, z_F, rdep_F, eis_F):
    coh_F = (1 + rdep_F) * dep_F_grid[np.newaxis, :] + z_F[:, np.newaxis]
    Vdep_F = (1 + rdep_F) * coh_F ** (-1 / eis_F)
    return Vdep_F

@sj.het(exogenous='Pi_F', policy='dep_F', backward='Vdep_F', backward_init=hh_init_F)
def hh_F(Vdep_F_p, dep_F_grid, z_F, t_paid_F, rdep_F, beta_F, eis_F):
    uc_nextgrid_F = beta_F * Vdep_F_p
    c_nextgrid_F = uc_nextgrid_F ** (-eis_F)
    coh_F = (1 + rdep_F) * dep_F_grid[np.newaxis, :] + z_F[:, np.newaxis]

    dep_F = sj.interpolate.interpolate_y(c_nextgrid_F + dep_F_grid, coh_F, dep_F_grid)
    sj.misc.setmin(dep_F, dep_F_grid[0])

    c_F = coh_F - dep_F
    uce_F = c_F ** (-1 / eis_F)
    Vdep_F = (1 + rdep_F) * uce_F

    tax_F = t_paid_F[:, np.newaxis] + np.zeros_like(dep_F_grid[np.newaxis, :])

    return Vdep_F, dep_F, c_F, uce_F, tax_F

def make_grids_F(Depmax_F, nDep_F, nZ_F, rho_z_F, sigma_z_F):
    dep_F_grid = grids.agrid(amax=Depmax_F, n=nDep_F)

    if nZ_F == 19:
        px_path_F = DATA_DIR_F / "Px_GMAR.txt"
        x_path_F  = DATA_DIR_F / "x_vec.txt"
        markov_ctstime_F = np.loadtxt(px_path_F)
        e_grid_F         = np.loadtxt(x_path_F).flatten()
        markov_Fistime_F = scipy.linalg.expm(markov_ctstime_F)
        row_sums_F       = markov_Fistime_F.sum(axis=1)
        Pi_F             = markov_Fistime_F / row_sums_F[:, None]
    else:
        e_grid_F, _, Pi_F = grids.markov_rouwenhorst(rho=rho_z_F, sigma=sigma_z_F, N=nZ_F)

    return dep_F_grid, e_grid_F, Pi_F

def income_F(e_grid_F, w_F, N_F, div_F, tau_F, lamb_F, P_CES_F):
    y_pre_F  = (w_F * N_F * e_grid_F + div_F) / P_CES_F   # real income in bundle units
    z_F      = lamb_F * (y_pre_F ** (1 - tau_F))
    t_paid_F = y_pre_F - z_F
    return z_F, t_paid_F

hh_extended_F = hh_F.add_hetinputs([make_grids_F, income_F])

# ── Steady-state blocks ───────────────────────────────────────────────────────
# CHECK STEADY STATE OWNING OF THE foreign bonds in DOMESTIC COUNTRY
@simple
def smart_steady_F(theta_F, Y_F, n_inter_F, rdep_F, alpha_F, delta_F, f_F, N_F,
                   rb_F, rb_actual_D, b_F_F, b_D_F, Q_F):
    # MU: bonds in common currency — no 1/p conversion of D-bond portfolio or balance sheet.
    K_F          = theta_F * n_inter_F
    phi_bF_F     = b_F_F / n_inter_F
    phi_bD_F     = b_D_F / n_inter_F              # no 1/p: same currency
    rk_F         = alpha_F * Y_F / K_F - delta_F
    rn_F         = theta_F * (rk_F - rdep_F) + phi_bF_F * (rb_F - rdep_F) + phi_bD_F * (rb_actual_D - rdep_F) + rdep_F
    m_F          = n_inter_F * (1 - (1 - f_F) * (1 + rn_F))
    k_inter_F    = K_F
    I_F          = K_F * delta_F
    D_supply_F   = (theta_F - 1) * n_inter_F + b_F_F + b_D_F  # no 1/p: same currency
    Z_F          = Y_F / ((K_F ** alpha_F) * (N_F ** (1 - alpha_F)))
    rdep_ante_F  = rdep_F
    cap_profit_F = Q_F * (K_F - (1 - delta_F) * K_F(-1)) - I_F
    return K_F, rk_F, rn_F, m_F, k_inter_F, I_F, D_supply_F, Z_F, rdep_ante_F, cap_profit_F

@simple
def market_clearing_F(Y_F, C_F, I_F, G_F, NX_F, DEP_F, D_supply_F):
    goods_mkt_F   = Y_F - C_F - I_F - G_F - NX_F
    deposit_mkt_F = DEP_F - D_supply_F
    return goods_mkt_F, deposit_mkt_F

@simple
def ces_price_F(omega, epsilon_trade, p):
    # F's domestic good price = 1 (in F-units); D-good import price = 1/p (in F-units).
    # So F's CES price index uses (1/p)^(1-eps), NOT p^(1-eps).
    P_CES_F = (omega + (1 - omega) * (1 / p) ** (1 - epsilon_trade)) ** (1 / (1 - epsilon_trade))
    return P_CES_F

@simple
def import_demand_F(C_F, I_F, G_F, omega, epsilon_trade, p, P_CES_F):
    A_F  = C_F + I_F + G_F
    # D-good price in F-units = 1/p; relative price for F's CES demand = P_CES_F / (1/p) = P_CES_F * p.
    # When p rises, D goods are cheaper in F-units → F imports more D goods.  (IM_F is in D-units.)
    IM_F = (1 - omega) * (P_CES_F * p) ** epsilon_trade * A_F
    return A_F, IM_F

@simple
def steady_auxilliary_F(theta_F, rk_F, rdep_F, delta_F, alpha_F, Y_F, K_F, N_F, lambda_gk_F, beta_F, ksi_F, rn_F):
    iota_F   = delta_F
    mpk_F    = alpha_F * (Y_F / K_F)
    w_F      = (1 - alpha_F) * Y_F / N_F
    Omega_F  = theta_F * lambda_gk_F / (beta_F * (1 + rn_F))
    nu_F     = beta_F * Omega_F * (rk_F - rdep_F)
    eta_F    = beta_F * Omega_F * (1 + rdep_F)
    gamma0_F = delta_F ** ksi_F / (1 - ksi_F)
    gamma1_F = -delta_F * ksi_F / (1 - ksi_F)
    return iota_F, mpk_F, w_F, Omega_F, nu_F, eta_F, gamma0_F, gamma1_F

@simple
def banker_div_F(rn_F, n_inter_F):
    div_F = rn_F * n_inter_F
    return div_F

@simple
def sdf_ss_F(beta_F):
    # SS-only block: SDF = beta (C constant at SS). Breaks cycle in ha.
    SDF_F = beta_F
    return SDF_F

@simple
def sdf_F(beta_F, C_F, eis_F):
    SDF_F = beta_F * (C_F(+1) / C_F) ** (-1 / eis_F)
    return SDF_F

@simple
def government_ss_F(TAX_F, rb_F, b_gov_F):
    # Use promised yield rb_F — consistent with budget_residual_F (interest rate channel).
    G_F = TAX_F - rb_F * b_gov_F
    return G_F

@simple
def labor_ss_F(w_F, N_F, UCE_F, frisch_F, mu_w_F):
    vphi_F = (1 / mu_w_F) * w_F * UCE_F / (N_F ** (1 / frisch_F))
    return vphi_F

@simple
def bond_return_F(rb_F, def_rate_F, recovery_rate_F):
    haircut_F   = 1.0 - recovery_rate_F
    rb_actual_F = (1 - def_rate_F * haircut_F) * (1 + rb_F(-1)) - 1
    return rb_actual_F

# ── Off-steady-state blocks ───────────────────────────────────────────────────

@simple
def capital_adj_F(Y_F, K_F, Q_F, I_F, alpha_F, delta_F, gamma0_F, gamma1_F, ksi_F, mc_F):
    iota_F        = I_F / K_F(-1)
    mpk_F         = alpha_F * Y_F / K_F(-1)
    rk_F          = (mc_F * mpk_F + (1 - delta_F) * Q_F) / Q_F(-1) - 1
    q_res_F       = Q_F - 1 / (gamma0_F * (1 - ksi_F) * iota_F ** (-ksi_F))
    capital_res_F = K_F - (1 - delta_F) * K_F(-1) - (gamma0_F * iota_F ** (1 - ksi_F) + gamma1_F) * K_F(-1)
    return iota_F, mpk_F, rk_F, q_res_F, capital_res_F

@simple
def labor_F(Y_F, Z_F, K_F, alpha_F, w_F):
    N_F  = (Y_F / (Z_F * K_F(-1) ** alpha_F)) ** (1 / (1 - alpha_F))
    mc_F = w_F * N_F / ((1 - alpha_F) * Y_F)
    return N_F, mc_F

@simple
def labor_market_F(w_F, UCE_F, N_F, vphi_F, frisch_F):
    mrs_F           = vphi_F * N_F ** (1 / frisch_F) / UCE_F
    labor_mkt_res_F = w_F - mrs_F
    return labor_mkt_res_F



@simple
def portfolio_foc_bD_F(rb_actual_D, rdep_F, b_D_F, n_inter_F, p,
                       phi_bD_F_ss, psi_bD_F, excess_return_D_F_ss, mp_wedge_D):
    # D-bank-style pricing of D-bonds held by F-banks: the wedge follows the
    # sovereign being priced (D), not the bank doing the pricing.
    # See domestic_bond_foc_D for the sign convention on mp_wedge_D.
    phi_bD_F            = (1 / p) * b_D_F / n_inter_F
    rb_actual_D_in_F_p1 = (p / p(+1)) * (1 + rb_actual_D(+1)) - 1
    foc_bD_res_F        = (rb_actual_D_in_F_p1 - rdep_F(+1)) - excess_return_D_F_ss \
                          - psi_bD_F * (phi_bD_F - phi_bD_F_ss) \
                          - mp_wedge_D
    return foc_bD_res_F


@simple
def intermediation_IC_F(nu_F, eta_F, lambda_gk_F):
    theta_F = eta_F / (lambda_gk_F - nu_F)
    return theta_F

@simple
def bank_return_F(theta_F, rk_F, rdep_F, b_F_F, b_D_F, n_inter_F, rb_actual_F, rb_actual_D,
                  phi_bD_F_ss, psi_bD_F):
    phi_bF_lag_F = b_F_F(-1) / n_inter_F(-1)
    phi_bD_lag_F = b_D_F(-1) / n_inter_F(-1)

    rn_F = (theta_F(-1) * (rk_F - rdep_F)
            + phi_bF_lag_F * (rb_actual_F - rdep_F)
            + phi_bD_lag_F * (rb_actual_D - rdep_F)
            + rdep_F)

    rn_F = rn_F - (psi_bD_F / 2) * (phi_bD_lag_F - phi_bD_F_ss) ** 2
    return rn_F

@simple
def intermediation_P1_F(rk_F, rdep_F, nu_F, lambda_gk_F, eta_F, theta_F, SDF_F, f_F):
    Omega_p1_F = f_F + (1 - f_F) * lambda_gk_F * theta_F(+1)
    nu_res_F   = nu_F  - SDF_F * Omega_p1_F * (rk_F(+1) - rdep_F(+1))
    eta_res_F  = eta_F - SDF_F * Omega_p1_F * (1 + rdep_F(+1))
    return nu_res_F, eta_res_F

@simple
def k_balance_sheet_F(Q_F, theta_F, n_inter_F, K_F):
    K_res_F = Q_F * K_F - theta_F * n_inter_F
    return K_res_F

@simple
def firm_profit_F(mc_F, Y_F):
    # Monopoly profit from sticky-price markup: (1 - mc_F) * Y_F.
    # Zero at SS (mc_F = 1), first-order off-SS. Routed to the financial
    # intermediary (banks own goods-producing firms in this GK setup).
    firm_profit_F = (1 - mc_F) * Y_F
    return firm_profit_F


@simple
def intermediation_P2_F(rn_F, n_inter_F, m_F, f_F, cap_profit_F, firm_profit_F):
    gross_income_F = (1 + rn_F) * n_inter_F(-1) + cap_profit_F + firm_profit_F
    n_inter_val_F  = (1 - f_F) * gross_income_F + m_F - n_inter_F
    return n_inter_val_F

@simple
def banker_div_res_F(rn_F, n_inter_F, div_F, m_F, f_F, cap_profit_F, firm_profit_F):
    gross_income_F = (1 + rn_F) * n_inter_F(-1) + cap_profit_F + firm_profit_F
    net_div_F      = f_F * gross_income_F - m_F
    div_res_F      = div_F - net_div_F
    return div_res_F

@simple
def intermediation_P3_F(Q_F, K_F, n_inter_F, b_F_F, b_D_F):
    # MU: b_D_F already in common currency — no 1/p conversion.
    D_supply_F = Q_F * K_F + b_F_F + b_D_F - n_inter_F
    return D_supply_F

@simple
def interest_rates_F(def_rate_F, recovery_rate_F, SDF_F):
    # Used in SS model only. In ha_full, rb_F is an unknown pinned by domestic_bond_foc_F.
    haircut_F = 1 - def_rate_F(+1) * (1.0 - recovery_rate_F)
    rb_F      = 1 / (SDF_F * haircut_F) - 1
    return rb_F


@simple
def domestic_bond_foc_F(rb_actual_F, rdep_F, b_F_F, n_inter_F,
                         phi_bF_F_ss, psi_bF_F, excess_return_bF_F_ss, mp_wedge_F):
    # GK-consistent bond pricing for F-country: rb_F adjusts until F-banks
    # willingly hold b_F_F (residual after D-banks take b_F_D).
    # At SS: phi_bF_F = phi_bF_F_ss AND mp_wedge_F = 0 → rb_F_res = 0.
    # See domestic_bond_foc_D for the sign convention on mp_wedge_F.
    phi_bF_F = b_F_F / n_inter_F
    rb_F_res = (rb_actual_F(+1) - rdep_F(+1)) - excess_return_bF_F_ss \
               - psi_bF_F * (phi_bF_F - phi_bF_F_ss) \
               - mp_wedge_F
    return rb_F_res


@simple
def government_default_F(shock_def_F, b_gov_F, Y_F, b_gov_ss_F, Y_ss_F,
                          def_scale_F, def_curvature_F, def_offset_F):
    # Symmetric to government_default_D — see comments there.
    debt_gap_F = b_gov_F(-1) / Y_F(-1) - b_gov_ss_F / Y_ss_F
    endog_F    = ((debt_gap_F + def_offset_F) ** def_curvature_F
                  - def_offset_F ** def_curvature_F)
    def_rate_F = shock_def_F + def_scale_F * endog_F
    return def_rate_F

@simple
def tax_rule_F(b_gov_F, lamb_ss_F, b_gov_ss_F, phi_lamb_F):
    lamb_F = lamb_ss_F - phi_lamb_F * (b_gov_F(-1) - b_gov_ss_F)
    return lamb_F



@simple
def capital_producer_profit_F(Q_F, K_F, I_F, delta_F):
    cap_profit_F = Q_F * (K_F - (1 - delta_F) * K_F(-1)) - I_F
    return cap_profit_F

@simple
def budget_residual_F(b_gov_F, G_F, TAX_F, rb_F, def_rate_F, recovery_rate_F, zeta_writeoff_F):
    # Symmetric to budget_residual_D — see comments there.
    haircut_F             = 1.0 - recovery_rate_F
    rb_actual_F_implied   = (1 - def_rate_F * haircut_F) * (1 + rb_F(-1)) - 1
    effective_repayment_F = (
        zeta_writeoff_F         * (1 + rb_actual_F_implied) * b_gov_F(-1)
        + (1 - zeta_writeoff_F) * (1 + rb_F(-1))            * b_gov_F(-1)
    )
    b_gov_res_F           = effective_repayment_F + G_F - TAX_F - b_gov_F
    return b_gov_res_F


# ── NK frictions ─────────────────────────────────────────────────────────────

@simple
def fisher_F(i_union, pi_F):
    rdep_ante_F = i_union - pi_F(+1)
    rdep_F      = i_union(-1) - pi_F
    return rdep_ante_F, rdep_F

@simple
def pricing_F(mc_F, pi_F, kappa_p_F, mu_p_F, SDF_F):
    nkpc_p_res_F = pi_F - kappa_p_F * (mc_F - 1 / mu_p_F) - SDF_F * pi_F(+1)
    return nkpc_p_res_F

@simple
def wage_setting_F(w_F, pi_w_F, pi_F, N_F, UCE_F, vphi_F, frisch_F, kappa_w_F, mu_w_F, beta_F):
    nkpc_w_res_F = (pi_w_F
                    - kappa_w_F * (vphi_F * N_F ** (1 + 1 / frisch_F)
                                   - (1 / mu_w_F) * w_F * N_F * UCE_F)
                    - beta_F * pi_w_F(+1))
    w_res_F = w_F - w_F(-1) * (1 + pi_w_F) / (1 + pi_F)
    return nkpc_w_res_F, w_res_F


