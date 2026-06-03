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

def hh_init_D(dep_D_grid, z_D, Rgross_D, eis_D, vphi_D, N_D, frisch_D):
    coh_D  = Rgross_D * dep_D_grid[np.newaxis, :] + z_D[:, np.newaxis]
    v_D    = vphi_D * N_D ** (1 + 1/frisch_D) / (1 + 1/frisch_D)
    Vdep_D = Rgross_D * (coh_D - v_D) ** (-1 / eis_D)
    return Vdep_D


@sj.het(exogenous='Pi_D', policy='dep_D', backward='Vdep_D', backward_init=hh_init_D)
def hh_D(Vdep_D_p, dep_D_grid, z_D, t_paid_D, Rgross_D, beta_D, eis_D, vphi_D, N_D, frisch_D):
    # GHH composite: x = c - v(N), v(N) = vphi*N^(1+1/frisch)/(1+1/frisch)
    v_D           = vphi_D * N_D ** (1 + 1/frisch_D) / (1 + 1/frisch_D)
    uc_nextgrid_D = beta_D * Vdep_D_p
    x_nextgrid_D  = uc_nextgrid_D ** (-eis_D)
    coh_D         = Rgross_D * dep_D_grid[np.newaxis, :] + z_D[:, np.newaxis]

    # EGM: x + a' = coh - v  →  endogenous grid x_nextgrid + a' matched to coh - v
    dep_D = sj.interpolate.interpolate_y(x_nextgrid_D + dep_D_grid, coh_D - v_D, dep_D_grid)
    sj.misc.setmin(dep_D, dep_D_grid[0])

    x_D    = coh_D - v_D - dep_D   # GHH composite
    c_D    = x_D + v_D              # total consumption = x + v(N)
    uce_D  = x_D ** (-1 / eis_D)   # marginal utility of composite
    Vdep_D = Rgross_D * uce_D

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


@simple
def deposit_return_D(rdep_D, P_CES_D):
    # Bundle-real gross deposit return: corrects for P_CES revaluation between t-1 and t.
    # At SS P_CES_D(-1)/P_CES_D = 1, so Rgross_D = 1 + rdep_D identically.
    Rgross_D = (1 + rdep_D) * P_CES_D(-1) / P_CES_D
    return Rgross_D


# ── STEADY STATE EQUATIONS ─── #############################################################################################

@simple
def smart_steady_D(theta_D, Y_D, n_inter_D, rdep_D, alpha_D, delta_D, f_D, N_D,
                   rb_actual_D, rb_actual_F, b_D_D, b_F_D, Q_D, q_b_D, q_b_F,
                   chi0_D, chi1_D, chi2_D, T0_D, T1_D, def_rate_D):
    K_D          = (theta_D * n_inter_D - q_b_D * b_D_D - q_b_F * b_F_D) / Q_D
    phi_bD_D     = q_b_D * b_D_D / n_inter_D
    phi_bF_D     = q_b_F * b_F_D / n_inter_D
    kappa_D      = theta_D - phi_bD_D - phi_bF_D
    rk_D         = alpha_D * Y_D / K_D - delta_D
    arg_D        = -rk_D * K_D / (K_D + chi0_D)
    Phi_D        = (chi1_D / chi2_D) * (arg_D ** 2) ** (chi2_D / 2) * (K_D + chi0_D)
    T_D          = (T0_D + T1_D * def_rate_D) * (b_D_D + b_F_D)
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
    cap_profit_D = Q_D * (K_D - (1 - delta_D) * K_D(-1)) - I_D
    return K_D, rk_D, rn_D, m_D, k_inter_D, I_D, D_supply_D, Z_D, cap_profit_D, Phi_D, T_D

@simple
def market_clearing_D(Y_D, C_D, I_D, G_D, NX_D, DEP_D, D_supply_D, P_CES_D, Phi_D, T_D, cap_profit_D):
    # cap_profit_D = Q·ΔK_net − I is added here to match its injection into bank net worth
    # via intermediation_P2_D. Without it the resource constraint debits only I while the
    # balance-sheet records Q·ΔK_net, leaving an unbacked (Q−1)·ΔK_net wedge.
    goods_mkt_D   = Y_D - (P_CES_D * C_D + I_D + G_D + Phi_D + T_D + cap_profit_D) - NX_D
    deposit_mkt_D = P_CES_D * DEP_D - D_supply_D
    return goods_mkt_D, deposit_mkt_D


@simple
def ces_price_D(omega, epsilon_trade, p):
    P_CES_D = (omega + (1 - omega) * p ** (1 - epsilon_trade)) ** (1 / (1 - epsilon_trade))
    return P_CES_D


@simple
def import_demand_D(C_D, I_D, G_D, omega, epsilon_trade, p, P_CES_D):
    A_D  = C_D + I_D + G_D
    IM_D = (1 - omega) * (P_CES_D / p) ** epsilon_trade * C_D
    return A_D, IM_D


@simple
def steady_auxilliary_D(theta_D, rk_D, rdep_D, delta_D, alpha_D, Y_D, K_D, N_D,
                        beta_D, ksi_D, rn_D, f_D,
                        rb_actual_D, rb_actual_F):
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
def banker_div_D(rn_D, n_inter_D, Phi_D, T_D):
    # Consistent with banker_div_res_D: div = f·gross − m = rn·n − Phi − T at SS.
    div_D = rn_D * n_inter_D - Phi_D - T_D
    return div_D


@simple
def sdf_ss_D(beta_D):
    SDF_D = beta_D
    return SDF_D

@simple
def sdf_D(beta_D, X_D, eis_D):
    # GHH: SDF uses composite x = c - v(N) instead of c
    SDF_D = beta_D * (X_D(+1) / X_D) ** (-1 / eis_D)
    return SDF_D


@simple
def ghh_composite_D(C_D, vphi_D, N_D, frisch_D):
    # Aggregate GHH composite X = C - v(N); homogeneous v(N) → X = C - v(N) exactly
    X_D = C_D - vphi_D * N_D ** (1 + 1/frisch_D) / (1 + 1/frisch_D)
    return X_D


@simple
def government_ss_D(TAX_D, q_b_D, b_gov_D, P_CES_D, delta_b_D,
                    def_rate_D, recovery_rate_D, zeta_writeoff_D):
    # At SS (b_gov constant): G = TAX*P + net_issuance - coupon.
    # Consistent with budget_residual_D at steady state.
    haircut_D   = 1.0 - recovery_rate_D
    surv_cont_D = 1.0 - zeta_writeoff_D * def_rate_D * haircut_D
    coupon_D    = delta_b_D * (1.0 - def_rate_D * haircut_D) * b_gov_D
    net_iss_D   = q_b_D * (1.0 - surv_cont_D * (1.0 - delta_b_D)) * b_gov_D
    G_D         = P_CES_D * TAX_D + net_iss_D - coupon_D
    return G_D


@simple
def labor_ss_D(w_D, N_D, frisch_D, mu_w_D, P_CES_D):
    # GHH: UCE = x^(-1/eis) cancels from intratemporal FOC → vphi independent of UCE
    vphi_D = (1 / mu_w_D) * (w_D / P_CES_D) / (N_D ** (1 / frisch_D))
    return vphi_D

@simple
def bond_return_D(def_rate_D, recovery_rate_D, q_b_D, delta_b_D, zeta_writeoff_D):
    haircut_D        = 1.0 - recovery_rate_D
    # Maturing coupon: always haircut by def_rate * haircut.
    # Continuation claim: haircut only when zeta_writeoff_D > 0 (partial/full write-down of outstanding stock).
    # zeta_writeoff_D = 0 → current model (only coupon haircut); = 1 → continuation also written down.
    current_payoff_D = delta_b_D * (1.0 - def_rate_D * haircut_D)
    continuation_D   = (1.0 - delta_b_D) * q_b_D * (1.0 - zeta_writeoff_D * def_rate_D * haircut_D)
    rb_actual_D      = (current_payoff_D + continuation_D) / q_b_D(-1) - 1.0
    return rb_actual_D


# ── OFF STEADY STATE EQUATIONS ─── #############################################################################################

@simple
def capital_adj_D(K_D, Q_D, I_D, Z_D, N_D, alpha_D, delta_D, gamma0_D, gamma1_D, ksi_D):
    iota_D        = I_D / K_D(-1)
    mpk_D         = alpha_D * Z_D * K_D ** (alpha_D - 1) * N_D ** (1 - alpha_D)
    rk_D          = (mpk_D + (1 - delta_D) * Q_D) / Q_D(-1) - 1
    q_res_D       = Q_D - 1 / (gamma0_D * (1 - ksi_D) * iota_D ** (-ksi_D))
    capital_res_D = K_D - (1 - delta_D) * K_D(-1) - (gamma0_D * iota_D ** (1 - ksi_D) + gamma1_D) * K_D(-1)
    return iota_D, mpk_D, rk_D, q_res_D, capital_res_D

@simple
def capital_producer_profit_D(Q_D, K_D, I_D, delta_D):
    cap_profit_D = Q_D * (K_D - (1 - delta_D) * K_D(-1)) - I_D
    return cap_profit_D


@simple
def labor_D(N_D, Z_D, K_D, alpha_D):
    Y_D = Z_D * K_D ** alpha_D * N_D ** (1 - alpha_D)
    return Y_D


@simple
def labor_market_D(w_D, N_D, vphi_D, frisch_D, P_CES_D):
    labor_mkt_res_D = w_D / P_CES_D - vphi_D * N_D ** (1 / frisch_D)
    return labor_mkt_res_D


@simple
def labor_demand_D(w_D, Y_D, N_D, alpha_D):
    w_res_D = w_D - (1 - alpha_D) * Y_D / N_D
    return w_res_D



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
    Omega_p1_D    = f_D + (1 - f_D) * lambda_gk_D * theta_D(+1)
    nu_K_res_D    = nu_K_D  - SDF_D * Omega_p1_D * (rk_D(+1)        - rdep_D(+1))
    nu_bD_res_D   = nu_bD_D - SDF_D * Omega_p1_D * (rb_actual_D(+1) - rdep_D(+1))
    nu_bF_res_D   = nu_bF_D - SDF_D * Omega_p1_D * (rb_actual_F(+1) - rdep_D(+1))
    eta_res_D     = eta_D   - SDF_D * Omega_p1_D * (1 + rdep_D(+1))
    return nu_K_res_D, nu_bD_res_D, nu_bF_res_D, eta_res_D


@simple
def k_balance_sheet_D(Q_D, theta_D, n_inter_D, K_D, b_D_D, b_F_D, q_b_D, q_b_F):
    K_res_D = Q_D * K_D + q_b_D * b_D_D + q_b_F * b_F_D - theta_D * n_inter_D
    return K_res_D



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
def intermediation_P2_D(rn_D, n_inter_D, m_D, f_D, cap_profit_D, Phi_D, T_D):
    # Writedown terms removed: rb_actual already embeds the default haircut via
    # rb_actual = (1 − def·haircut)/q_b(-1) − 1, so deducting them again double-counts.
    gross_income_D = (1 + rn_D) * n_inter_D(-1) + cap_profit_D
    n_inter_val_D  = (1 - f_D) * gross_income_D + m_D - Phi_D - T_D - n_inter_D
    return n_inter_val_D


@simple
def banker_div_res_D(rn_D, n_inter_D, div_D, m_D, f_D, cap_profit_D):
    gross_income_D = (1 + rn_D) * n_inter_D(-1) + cap_profit_D
    net_div_D      = f_D * gross_income_D - m_D
    div_res_D      = div_D - net_div_D
    return div_res_D


@simple
def intermediation_P3_D(Q_D, K_D, n_inter_D, b_D_D, b_F_D, q_b_D, q_b_F):
    D_supply_D = Q_D * K_D + q_b_D * b_D_D + q_b_F * b_F_D - n_inter_D
    return D_supply_D


@simple
def bond_price_ss_D(SDF_D, def_rate_D, recovery_rate_D, delta_b_D, zeta_writeoff_D):
    # SS fixed point: q = SDF*[delta_b*(1-d*h) + (1-delta_b)*q*(1-zeta*d*h)]
    # → q = SDF*delta_b*(1-d*h) / (1 - SDF*(1-delta_b)*(1-zeta*d*h))
    haircut_D   = 1.0 - recovery_rate_D
    surv_cont_D = 1.0 - zeta_writeoff_D * def_rate_D * haircut_D
    q_b_D       = (
        SDF_D * delta_b_D * (1.0 - def_rate_D * haircut_D)
        / (1.0 - SDF_D * (1.0 - delta_b_D) * surv_cont_D)
    )
    return q_b_D


@simple
def domestic_bond_foc_D(rb_actual_D, rdep_D, b_D_D, n_inter_D, q_b_D,
                         phi_bD_D_ss, psi_bD_D, excess_return_bD_D_ss, tau_mp_D):
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
def budget_residual_D(b_gov_D, G_D, TAX_D, q_b_D, def_rate_D, recovery_rate_D, zeta_writeoff_D, P_CES_D, delta_b_D):
    haircut_D      = 1.0 - recovery_rate_D
    surv_cont_D    = 1.0 - zeta_writeoff_D * def_rate_D * haircut_D
    # Coupon: maturing share always fully haircut (consistent with current_payoff in bond_return_D).
    coupon_D       = delta_b_D * (1.0 - def_rate_D * haircut_D) * b_gov_D(-1)
    # Net issuance: unmatured legacy stock written down by surv_cont when zeta > 0.
    net_issuance_D = q_b_D * (b_gov_D - surv_cont_D * (1.0 - delta_b_D) * b_gov_D(-1))
    b_gov_res_D    = coupon_D + G_D - P_CES_D * TAX_D - net_issuance_D
    return b_gov_res_D

