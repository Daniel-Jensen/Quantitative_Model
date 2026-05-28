import numpy as np
import scipy.linalg
import sequence_jacobian as sj
from sequence_jacobian import simple
from sequence_jacobian import grids
from pathlib import Path

try:
    BASE_DIR_D = Path(__file__).resolve().parent
except NameError:
    BASE_DIR_D = Path.cwd()

DATA_DIR_D = BASE_DIR_D / "Discretisation" / "Outputs"


# ── HOUSEHOLD ─── #############################################################################################

def hh_init_D(dep_D_grid, z_D, rdep_D, eis_D):
    coh_D = (1 + rdep_D) * dep_D_grid[np.newaxis, :] + z_D[:, np.newaxis]
    Vdep_D = (1 + rdep_D) * coh_D ** (-1 / eis_D)
    return Vdep_D


@sj.het(exogenous='Pi_D', policy='dep_D', backward='Vdep_D', backward_init=hh_init_D)
def hh_D(Vdep_D_p, dep_D_grid, z_D, t_paid_D, rdep_D, beta_D, eis_D):
    uc_nextgrid_D = beta_D * Vdep_D_p
    c_nextgrid_D = uc_nextgrid_D ** (-eis_D)
    coh_D = (1 + rdep_D) * dep_D_grid[np.newaxis, :] + z_D[:, np.newaxis]

    dep_D = sj.interpolate.interpolate_y(c_nextgrid_D + dep_D_grid, coh_D, dep_D_grid)
    sj.misc.setmin(dep_D, dep_D_grid[0])

    c_D = coh_D - dep_D
    uce_D = c_D ** (-1 / eis_D)
    Vdep_D = (1 + rdep_D) * uce_D

    tax_D = t_paid_D[:, np.newaxis] + np.zeros_like(dep_D_grid[np.newaxis, :])

    return Vdep_D, dep_D, c_D, uce_D, tax_D


def make_grids_D(Depmax_D, nDep_D, nZ_D, rho_z_D, sigma_z_D):
    dep_D_grid = grids.agrid(amax=Depmax_D, n=nDep_D)

    if nZ_D == 19:
        px_path_D = DATA_DIR_D / "Px_GMAR.txt"
        x_path_D  = DATA_DIR_D / "x_vec.txt"
        markov_ctstime_D = np.loadtxt(px_path_D)
        e_grid_D         = np.loadtxt(x_path_D).flatten()
        markov_distime_D = scipy.linalg.expm(markov_ctstime_D)
        row_sums_D       = markov_distime_D.sum(axis=1)
        Pi_D             = markov_distime_D / row_sums_D[:, None]
    else:
        e_grid_D, _, Pi_D = grids.markov_rouwenhorst(rho=rho_z_D, sigma=sigma_z_D, N=nZ_D)

    return dep_D_grid, e_grid_D, Pi_D


def income_D(e_grid_D, w_D, N_D, div_D, tau_D, lamb_D, P_CES_D):
    y_pre_D  = (w_D * N_D * e_grid_D + div_D) / P_CES_D   # real income in bundle units
    z_D      = lamb_D * (y_pre_D ** (1 - tau_D))
    t_paid_D = y_pre_D - z_D
    return z_D, t_paid_D


hh_extended_D = hh_D.add_hetinputs([make_grids_D, income_D])


# ── STEADY STATE EQUATIONS ─── #############################################################################################

@simple
def smart_steady_D(theta_D, Y_D, n_inter_D, rdep_D, alpha_D, delta_D, f_D, N_D,
                   rb_actual_D, rb_actual_F, b_D_D, b_F_D, Q_D, q_b_D, q_b_F,
                   chi0_D, chi1_D, chi2_D, T0_D, T1_D, def_rate_D):
    K_D          = (theta_D * n_inter_D - q_b_D * b_D_D - q_b_F * b_F_D) / Q_D
    phi_bD_D     = q_b_D * b_D_D / n_inter_D   # market-value share of D-bonds
    phi_bF_D     = q_b_F * b_F_D / n_inter_D   # market-value share of F-bonds
    kappa_D      = theta_D - phi_bD_D - phi_bF_D
    rk_D         = alpha_D * Y_D / K_D - delta_D
    arg_D        = -rk_D * K_D / (K_D + chi0_D)
    Phi_D        = (chi1_D / chi2_D) * (arg_D ** 2) ** (chi2_D / 2) * (K_D + chi0_D)
    # Macroprudential bond tax at SS
    T_D          = (T0_D + T1_D * def_rate_D) * (b_D_D + b_F_D)
    # rn = PURE portfolio return.  Phi (intermediary cap-adj cost) and T (bond
    # tax) hit net worth directly in intermediation_P2_D, not via the return
    # rate.  Keeping them out of rn ensures Omega/ν/η computed here are
    # consistent with the transition's bank_return_D.
    rn_D         = (kappa_D * (rk_D - rdep_D)
                    + phi_bD_D * (rb_actual_D - rdep_D)
                    + phi_bF_D * (rb_actual_F - rdep_D)
                    + rdep_D)
    # SS net worth identity (intermediation_P2_D = 0 at SS):
    #   n = (1-f)·(1+rn)·n + m - Phi - T   ⟹   m = n·(1 - (1-f)(1+rn)) + Phi + T
    m_D          = n_inter_D * (1 - (1 - f_D) * (1 + rn_D)) + Phi_D + T_D
    k_inter_D    = K_D
    I_D          = K_D * delta_D
    D_supply_D   = (theta_D - 1) * n_inter_D
    Z_D          = Y_D / ((K_D ** alpha_D) * (N_D ** (1 - alpha_D)))
    rdep_ante_D  = rdep_D
    cap_profit_D = Q_D * (K_D - (1 - delta_D) * K_D(-1)) - I_D
    return K_D, rk_D, rn_D, m_D, k_inter_D, I_D, D_supply_D, Z_D, rdep_ante_D, cap_profit_D, Phi_D, T_D

@simple
def market_clearing_D(Y_D, C_D, I_D, G_D, NX_D, DEP_D, D_supply_D, P_CES_D):
    goods_mkt_D   = Y_D - P_CES_D * (C_D + I_D + G_D) - NX_D
    deposit_mkt_D = P_CES_D * DEP_D - D_supply_D
    return goods_mkt_D, deposit_mkt_D


@simple
def ces_price_D(omega, epsilon_trade, p):
    P_CES_D = (omega + (1 - omega) * p ** (1 - epsilon_trade)) ** (1 / (1 - epsilon_trade))
    return P_CES_D


@simple
def import_demand_D(C_D, I_D, G_D, omega, epsilon_trade, p, P_CES_D):
    A_D  = C_D + I_D + G_D
    IM_D = (1 - omega) * (P_CES_D / p) ** epsilon_trade * A_D
    return A_D, IM_D


@simple
def steady_auxilliary_D(theta_D, rk_D, rdep_D, delta_D, alpha_D, Y_D, K_D, N_D,
                        beta_D, ksi_D, rn_D, f_D,
                        rb_actual_D, rb_actual_F):
    # Multi-asset GK: separate Bellman ν for each risky asset.
    # λ_gk is DERIVED so that the transition's Bellman Ω = f + (1-f)·λ·θ is
    # consistent with the SS IC binding β·Ω·(1+rn) = λ·θ.  Joint solution:
    #     λ = f / (θ·[1/(β(1+rn)) − (1−f)])
    # This eliminates the SS-vs-transition inconsistency that produced
    # nonzero P1 Bellman residuals when λ was fixed at calibration.
    iota_D       = delta_D
    mpk_D        = alpha_D * (Y_D / K_D)
    w_D          = (1 - alpha_D) * Y_D / N_D
    lambda_gk_D  = f_D / (theta_D * (1 / (beta_D * (1 + rn_D)) - (1 - f_D)))
    Omega_D      = f_D + (1 - f_D) * lambda_gk_D * theta_D
    nu_K_D       = beta_D * Omega_D * (rk_D        - rdep_D)
    nu_bD_D      = beta_D * Omega_D * (rb_actual_D - rdep_D)
    nu_bF_D      = beta_D * Omega_D * (rb_actual_F - rdep_D)
    eta_D        = beta_D * Omega_D * (1 + rdep_D)
    gamma0_D     = delta_D ** ksi_D / (1 - ksi_D)
    gamma1_D     = -delta_D * ksi_D / (1 - ksi_D)
    return iota_D, mpk_D, w_D, Omega_D, lambda_gk_D, nu_K_D, nu_bD_D, nu_bF_D, eta_D, gamma0_D, gamma1_D


@simple
def banker_div_D(rn_D, n_inter_D):
    div_D = rn_D * n_inter_D
    return div_D


@simple
def sdf_ss_D(beta_D):
    SDF_D = beta_D
    return SDF_D

@simple
def sdf_D(beta_D, C_D, eis_D):
    SDF_D = beta_D * (C_D(+1) / C_D) ** (-1 / eis_D)
    return SDF_D


@simple
def government_ss_D(TAX_D, q_b_D, b_gov_D):
    # At SS: gov raises q_b_D per face-value unit issued; interest cost = (1−q_b_D)·b_gov.
    G_D = TAX_D - (1.0 - q_b_D) * b_gov_D
    return G_D


@simple
def labor_ss_D(w_D, N_D, UCE_D, frisch_D, mu_w_D):
    vphi_D = (1 / mu_w_D) * w_D * UCE_D / (N_D ** (1 / frisch_D))
    return vphi_D

@simple
def bond_return_D(def_rate_D, recovery_rate_D, q_b_D):
    haircut_D   = 1.0 - recovery_rate_D
    rb_actual_D = (1 - def_rate_D * haircut_D) / q_b_D(-1) - 1
    return rb_actual_D


# ── OFF STEADY STATE EQUATIONS ─── #############################################################################################

@simple
def capital_adj_D(Y_D, K_D, Q_D, I_D, alpha_D, delta_D, gamma0_D, gamma1_D, ksi_D, mc_D):
    iota_D        = I_D / K_D(-1)
    mpk_D         = alpha_D * Y_D / K_D(-1)
    # Under monopolistic competition firms rent capital at mc * MPK (not just MPK).
    rk_D          = (mc_D * mpk_D + (1 - delta_D) * Q_D) / Q_D(-1) - 1
    q_res_D       = Q_D - 1 / (gamma0_D * (1 - ksi_D) * iota_D ** (-ksi_D))
    capital_res_D = K_D - (1 - delta_D) * K_D(-1) - (gamma0_D * iota_D ** (1 - ksi_D) + gamma1_D) * K_D(-1)
    return iota_D, mpk_D, rk_D, q_res_D, capital_res_D

@simple
def capital_producer_profit_D(Q_D, K_D, I_D, delta_D):
    cap_profit_D = Q_D * (K_D - (1 - delta_D) * K_D(-1)) - I_D
    return cap_profit_D


@simple
def labor_D(Y_D, Z_D, K_D, alpha_D, w_D):
    N_D  = (Y_D / (Z_D * K_D(-1) ** alpha_D)) ** (1 / (1 - alpha_D))
    mc_D = w_D * N_D / ((1 - alpha_D) * Y_D)
    return N_D, mc_D


@simple
def labor_market_D(w_D, UCE_D, N_D, vphi_D, frisch_D):
    mrs_D           = vphi_D * N_D ** (1 / frisch_D) / UCE_D
    labor_mkt_res_D = w_D - mrs_D
    return labor_mkt_res_D



@simple
def intermediation_IC_D(nu_K_D, nu_bD_D, nu_bF_D, eta_D,
                        Q_D, K_D, q_b_D, q_b_F, b_D_D, b_F_D, n_inter_D,
                        lambda_gk_D, theta_D):
    kappa_D     = Q_D   * K_D   / n_inter_D
    phi_bD_D    = q_b_D * b_D_D / n_inter_D
    phi_bF_D    = q_b_F * b_F_D / n_inter_D
    theta_tgt_D = (nu_K_D * kappa_D + nu_bD_D * phi_bD_D + nu_bF_D * phi_bF_D + eta_D) / lambda_gk_D
    ic_res_D    = theta_D - theta_tgt_D
    return ic_res_D


@simple
def bank_return_D(theta_D, rk_D, rdep_D, b_D_D, b_F_D, n_inter_D,
                  rb_actual_D, rb_actual_F, q_b_D, q_b_F):
    # Portfolio shares use MARKET VALUE of bond positions at t-1.
    phi_bD_lag_D = q_b_D(-1) * b_D_D(-1) / n_inter_D(-1)
    phi_bF_lag_D = q_b_F(-1) * b_F_D(-1) / n_inter_D(-1)
    kappa_lag_D  = theta_D(-1) - phi_bD_lag_D - phi_bF_lag_D
    rn_D = (kappa_lag_D  * (rk_D - rdep_D)
            + phi_bD_lag_D * (rb_actual_D - rdep_D)
            + phi_bF_lag_D * (rb_actual_F - rdep_D)
            + rdep_D)
    return rn_D


@simple
def intermediation_P1_D(rk_D, rb_actual_D, rb_actual_F, rdep_D,
                        nu_K_D, nu_bD_D, nu_bF_D, eta_D,
                        lambda_gk_D, theta_D, SDF_D, f_D):
    # Bellman: marginal franchise value of each asset.  Ω is the augmented SDF
    # weighting (continuation value of leverage).  Three ν's let default risk
    # bite ν_bD without contaminating ν_K — sovereign shocks tighten the IC
    # directly, not just via writedown.
    Omega_p1_D    = f_D + (1 - f_D) * lambda_gk_D * theta_D(+1)
    nu_K_res_D    = nu_K_D  - SDF_D * Omega_p1_D * (rk_D(+1)        - rdep_D(+1))
    nu_bD_res_D   = nu_bD_D - SDF_D * Omega_p1_D * (rb_actual_D(+1) - rdep_D(+1))
    nu_bF_res_D   = nu_bF_D - SDF_D * Omega_p1_D * (rb_actual_F(+1) - rdep_D(+1))
    eta_res_D     = eta_D   - SDF_D * Omega_p1_D * (1 + rdep_D(+1))
    return nu_K_res_D, nu_bD_res_D, nu_bF_res_D, eta_res_D


@simple
def k_balance_sheet_D(Q_D, theta_D, n_inter_D, K_D, b_D_D, b_F_D, q_b_D, q_b_F):
    # Balance sheet at MARKET VALUES: Q*K + q_b_D*b_D_D + q_b_F*b_F_D = theta*n.
    K_res_D = Q_D * K_D + q_b_D * b_D_D + q_b_F * b_F_D - theta_D * n_inter_D
    return K_res_D


@simple
def firm_profit_D(mc_D, Y_D):
    firm_profit_D = (1 - mc_D) * Y_D
    return firm_profit_D


@simple
def cap_adj_cost_inter_D(K_D, rk_D, chi0_D, chi1_D, chi2_D):
    # Auclert (2019) intermediary capital adjustment cost.
    arg_D = (K_D - (1.0 + rk_D) * K_D(-1)) / (K_D(-1) + chi0_D)
    Phi_D = (chi1_D / chi2_D) * (arg_D ** 2) ** (chi2_D / 2) * (K_D(-1) + chi0_D)
    return Phi_D


@simple
def macro_pru_tax_D(b_D_D, b_F_D, def_rate_D, T0_D, T1_D):
    # Macroprudential bond tax: T = (T0 + T1·ProbDefault) · total bond holdings
    # (both D-bonds and F-bonds held by the D-bank).  Matches smart_steady_D.
    tau_mp_D = T0_D + T1_D * def_rate_D
    T_D      = tau_mp_D * (b_D_D + b_F_D)
    return tau_mp_D, T_D


@simple
def intermediation_P2_D(rn_D, n_inter_D, m_D, f_D, cap_profit_D, firm_profit_D,
                        b_D_D, def_rate_D, recovery_rate_D,
                        b_F_D, def_rate_F, recovery_rate_F,
                        Phi_D, T_D):
    haircut_D      = 1.0 - recovery_rate_D
    haircut_F      = 1.0 - recovery_rate_F
    writedown_D    = def_rate_D * haircut_D * b_D_D(-1)
    writedown_F    = def_rate_F * haircut_F * b_F_D(-1)
    gross_income_D = (1 + rn_D) * n_inter_D(-1) + cap_profit_D + firm_profit_D
    n_inter_val_D  = ((1 - f_D) * gross_income_D + m_D
                      - writedown_D - writedown_F - Phi_D - T_D - n_inter_D)
    return n_inter_val_D


@simple
def banker_div_res_D(rn_D, n_inter_D, div_D, m_D, f_D, cap_profit_D, firm_profit_D):
    gross_income_D = (1 + rn_D) * n_inter_D(-1) + cap_profit_D + firm_profit_D
    net_div_D      = f_D * gross_income_D - m_D
    div_res_D      = div_D - net_div_D
    return div_res_D


@simple
def intermediation_P3_D(Q_D, K_D, n_inter_D, b_D_D, b_F_D, q_b_D, q_b_F):
    D_supply_D = Q_D * K_D + q_b_D * b_D_D + q_b_F * b_F_D - n_inter_D
    return D_supply_D


@simple
def bond_price_ss_D(SDF_D, def_rate_D, recovery_rate_D):
    # No-arbitrage SS bond price: pay q_b_D today, receive (1−default·haircut) at maturity.
    # In ha_full q_b_D is an unknown pinned by domestic_bond_foc_D; this block is SS-only.
    haircut_D = 1.0 - recovery_rate_D
    q_b_D     = SDF_D * (1.0 - def_rate_D(+1) * haircut_D)
    return q_b_D


@simple
def domestic_bond_foc_D(rb_actual_D, rdep_D, b_D_D, n_inter_D, q_b_D,
                         phi_bD_D_ss, psi_bD_D, excess_return_bD_D_ss, tau_mp_D):
    # Portfolio share at market value; tau_mp_D is the marginal macroprudential tax.
    phi_bD_D = q_b_D * b_D_D / n_inter_D
    rb_D_res = (rb_actual_D(+1) - rdep_D(+1)) - excess_return_bD_D_ss \
               - psi_bD_D * (phi_bD_D - phi_bD_D_ss) \
               - tau_mp_D
    return rb_D_res


# ==> GOVERMENT EQUATIONS
@simple
def government_default_D(shock_def_D, b_gov_D, Y_D, b_gov_ss_D, Y_ss_D,
                          def_scale_D, def_curvature_D, def_offset_D):
    debt_gap_D = b_gov_D(-1) / Y_D(-1) - b_gov_ss_D / Y_ss_D
    endog_D    = ((debt_gap_D + def_offset_D) ** def_curvature_D
                  - def_offset_D ** def_curvature_D)
    def_rate_D = shock_def_D + def_scale_D * endog_D
    return def_rate_D


@simple
def tax_rule_D(b_gov_D, lamb_ss_D, b_gov_ss_D, phi_lamb_D):
    lamb_D = lamb_ss_D - phi_lamb_D * (b_gov_D(-1) - b_gov_ss_D)
    return lamb_D


@simple
def budget_residual_D(b_gov_D, G_D, TAX_D, q_b_D, def_rate_D, recovery_rate_D, zeta_writeoff_D):
    # Budget: q_b_D·b_gov_D + TAX_D = G_D + repayment_of_old_bonds
    # b_gov_D is FACE VALUE; gov raises q_b_D per unit issued.
    # Repayment: zeta=1 → write off default losses; zeta=0 → pay full face value.
    haircut_D             = 1.0 - recovery_rate_D
    effective_repayment_D = (1.0 - zeta_writeoff_D * def_rate_D * haircut_D) * b_gov_D(-1)
    b_gov_res_D           = effective_repayment_D + G_D - TAX_D - q_b_D * b_gov_D
    return b_gov_res_D

