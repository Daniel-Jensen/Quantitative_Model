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

def hh_init_F(dep_F_grid, z_F, Rgross_F, eis_F, vphi_F, N_F, frisch_F):
    coh_F  = Rgross_F * dep_F_grid[np.newaxis, :] + z_F[:, np.newaxis]
    v_F    = vphi_F * N_F ** (1 + 1/frisch_F) / (1 + 1/frisch_F)
    Vdep_F = Rgross_F * (coh_F - v_F) ** (-1 / eis_F)
    return Vdep_F

@sj.het(exogenous='Pi_F', policy='dep_F', backward='Vdep_F', backward_init=hh_init_F)
def hh_F(Vdep_F_p, dep_F_grid, z_F, t_paid_F, Rgross_F, beta_F, eis_F, vphi_F, N_F, frisch_F):
    # GHH composite: x = c - v(N), v(N) = vphi*N^(1+1/frisch)/(1+1/frisch)
    v_F           = vphi_F * N_F ** (1 + 1/frisch_F) / (1 + 1/frisch_F)
    uc_nextgrid_F = beta_F * Vdep_F_p
    x_nextgrid_F  = uc_nextgrid_F ** (-eis_F)
    coh_F         = Rgross_F * dep_F_grid[np.newaxis, :] + z_F[:, np.newaxis]

    # EGM: x + a' = coh - v  →  endogenous grid x_nextgrid + a' matched to coh - v
    dep_F = sj.interpolate.interpolate_y(x_nextgrid_F + dep_F_grid, coh_F - v_F, dep_F_grid)
    sj.misc.setmin(dep_F, dep_F_grid[0])

    x_F    = coh_F - v_F - dep_F   # GHH composite
    c_F    = x_F + v_F              # total consumption = x + v(N)
    uce_F  = x_F ** (-1 / eis_F)   # marginal utility of composite
    Vdep_F = Rgross_F * uce_F

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

def income_F(e_grid_F, w_F, N_F, div_F, tau_F, lamb_F, P_CES_F, T_ls_F):
    y_pre_F  = (w_F * N_F * e_grid_F + div_F) / P_CES_F
    z_F      = lamb_F * (y_pre_F ** (1 - tau_F)) - T_ls_F
    t_paid_F = y_pre_F - z_F
    return z_F, t_paid_F

hh_extended_F = hh_F.add_hetinputs([make_grids_F, income_F])

@simple
def deposit_return_F(rdep_F, P_CES_F):
    # Bundle-real gross deposit return: corrects for P_CES revaluation between t-1 and t.
    # At SS P_CES_F(-1)/P_CES_F = 1, so Rgross_F = 1 + rdep_F identically.
    Rgross_F = (1 + rdep_F) * P_CES_F(-1) / P_CES_F
    return Rgross_F


# ── Steady-state blocks ───────────────────────────────────────────────────────
# CHECK STEADY STATE OWNING OF THE foreign bonds in DOMESTIC COUNTRY
@simple
def smart_steady_F(theta_F, Y_F, n_inter_F, rdep_F, alpha_F, delta_F, f_F, N_F,
                   rb_actual_F, rb_actual_D, b_F_F, b_D_F, Q_F, q_b_F, q_b_D,
                   chi0_F, chi1_F, chi2_F, T0_F, T1_F, def_rate_F, p):
    # Bonds are D-good (numeraire) claims; F-bank assets/NW are in F-goods → divide by p.
    K_F          = (theta_F * n_inter_F - (q_b_F * b_F_F + q_b_D * b_D_F) / p) / Q_F
    phi_bF_F     = q_b_F * b_F_F / (p * n_inter_F)
    phi_bD_F     = q_b_D * b_D_F / (p * n_inter_F)
    kappa_F      = theta_F - phi_bF_F - phi_bD_F
    rk_F         = alpha_F * Y_F / K_F - delta_F
    arg_F        = -rk_F * K_F / (K_F + chi0_F)
    Phi_F        = (chi1_F / chi2_F) * (arg_F ** 2) ** (chi2_F / 2) * (K_F + chi0_F)
    T_F          = (T0_F + T1_F * def_rate_F) * (b_F_F + b_D_F) / p
    # rn = PURE portfolio return.  See smart_steady_D.
    rn_F         = (kappa_F * (rk_F - rdep_F)
                    + phi_bF_F * (rb_actual_F - rdep_F)
                    + phi_bD_F * (rb_actual_D - rdep_F)
                    + rdep_F)
    m_F          = n_inter_F * (1 - (1 - f_F) * (1 + rn_F)) + Phi_F + T_F
    k_inter_F    = K_F
    I_F          = K_F * delta_F
    D_supply_F   = (theta_F - 1) * n_inter_F
    Z_F          = Y_F / ((K_F ** alpha_F) * (N_F ** (1 - alpha_F)))
    cap_profit_F = Q_F * (K_F - (1 - delta_F) * K_F(-1)) - I_F
    return K_F, rk_F, rn_F, m_F, k_inter_F, I_F, D_supply_F, Z_F, cap_profit_F, Phi_F, T_F

@simple
def market_clearing_F(Y_F, C_F, I_F, G_F, NX_F, DEP_F, D_supply_F, P_CES_F, Phi_F, T_F):
    # F (P_F = 1): only C is in bundle units; I, G, Phi, T are in domestic goods.
    goods_mkt_F = Y_F - (P_CES_F * C_F + I_F + G_F + Phi_F + T_F) - NX_F
    deposit_mkt_F = P_CES_F * DEP_F - D_supply_F
    return goods_mkt_F, deposit_mkt_F

@simple
def ces_price_F(omega, epsilon_trade, p):
    # F's domestic good price = 1 (in F-units); D-good import price = 1/p (in F-units).
    # So F's CES price index uses (1/p)^(1-eps), NOT p^(1-eps).
    P_CES_F = (omega + (1 - omega) * (1 / p) ** (1 - epsilon_trade)) ** (1 / (1 - epsilon_trade))
    return P_CES_F

@simple
def import_demand_F(C_F, omega, epsilon_trade, p, P_CES_F):
    IM_F = (1 - omega) * (P_CES_F * p) ** epsilon_trade * C_F
    return IM_F

@simple
def steady_auxilliary_F(theta_F, rk_F, rdep_F, delta_F, alpha_F, Y_F, K_F, N_F,
                        beta_inter_F, ksi_F, rn_F, f_F,
                        rb_actual_F, rb_actual_D):
    iota_F       = delta_F
    mpk_F        = alpha_F * (Y_F / K_F)
    w_F          = (1 - alpha_F) * Y_F / N_F
    lambda_gk_F  = f_F / (theta_F * (1 / (beta_inter_F * (1 + rn_F)) - (1 - f_F)))
    Omega_F      = f_F + (1 - f_F) * lambda_gk_F * theta_F
    nu_K_F       = beta_inter_F * Omega_F * (rk_F        - rdep_F)
    nu_bF_F      = beta_inter_F * Omega_F * (rb_actual_F - rdep_F)
    nu_bD_F      = beta_inter_F * Omega_F * (rb_actual_D - rdep_F)
    eta_F        = beta_inter_F * Omega_F * (1 + rdep_F)
    gamma0_F     = delta_F ** ksi_F / (1 - ksi_F)
    gamma1_F     = -delta_F * ksi_F / (1 - ksi_F)
    return iota_F, mpk_F, w_F, Omega_F, lambda_gk_F, nu_K_F, nu_bF_F, nu_bD_F, eta_F, gamma0_F, gamma1_F

@simple
def banker_div_F(rn_F, n_inter_F, Phi_F, T_F):
    # Consistent with banker_div_res_F: div = f·gross − m = rn·n − Phi − T at SS.
    div_F = rn_F * n_inter_F - Phi_F - T_F
    return div_F

@simple
def sdf_ss_F(beta_F):
    # SS-only block: SDF = beta (C constant at SS). Breaks cycle in ha.
    SDF_F = beta_F
    return SDF_F

@simple
def sdf_F(beta_F, X_F, eis_F):
    # GHH: SDF uses composite x = c - v(N) instead of c
    SDF_F = beta_F * (X_F(+1) / X_F) ** (-1 / eis_F)
    return SDF_F

@simple
def sdf_banker_ss_F(beta_inter_F):
    SDF_banker_F = beta_inter_F
    return SDF_banker_F

@simple
def sdf_banker_F(beta_inter_F, X_F, eis_F):
    SDF_banker_F = beta_inter_F * (X_F(+1) / X_F) ** (-1 / eis_F)
    return SDF_banker_F


@simple
def ghh_composite_F(C_F, vphi_F, N_F, frisch_F):
    # Aggregate GHH composite X = C - v(N); homogeneous v(N) → X = C - v(N) exactly
    X_F = C_F - vphi_F * N_F ** (1 + 1/frisch_F) / (1 + 1/frisch_F)
    return X_F

@simple
def government_ss_F(TAX_F, q_b_F, b_gov_F, p, P_CES_F, delta_b_F,
                    def_rate_F, recovery_rate_F, zeta_writeoff_F, writeoff_enabled_F):
    haircut_F      = 1.0 - recovery_rate_F
    haircut_mult_F = writeoff_enabled_F
    surv_cont_F    = 1.0 - zeta_writeoff_F * def_rate_F * haircut_F * haircut_mult_F
    coupon_F       = delta_b_F * (1.0 - def_rate_F * haircut_F * haircut_mult_F) * b_gov_F
    net_iss_F      = q_b_F * (1.0 - surv_cont_F * (1.0 - delta_b_F)) * b_gov_F
    G_F            = P_CES_F * TAX_F + (net_iss_F - coupon_F) / p
    return G_F

@simple
def labor_ss_F(w_F, N_F, frisch_F, mu_w_F, P_CES_F):
    # GHH: UCE cancels from intratemporal FOC → vphi independent of UCE
    vphi_F = (1 / mu_w_F) * (w_F / P_CES_F) / (N_F ** (1 / frisch_F))
    return vphi_F

@simple
def bond_return_F(def_rate_F, recovery_rate_F, q_b_F, delta_b_F, zeta_writeoff_F, writeoff_enabled_F):
    # writeoff_enabled_F = 0: pure sovereign risk shock, no haircuts on cash flows.
    haircut_F        = 1.0 - recovery_rate_F
    haircut_mult_F   = writeoff_enabled_F
    current_payoff_F = delta_b_F * (1.0 - def_rate_F * haircut_F * haircut_mult_F)
    continuation_F   = (1.0 - delta_b_F) * q_b_F * (1.0 - zeta_writeoff_F * def_rate_F * haircut_F * haircut_mult_F)
    rb_actual_F      = (current_payoff_F + continuation_F) / q_b_F(-1) - 1.0
    return rb_actual_F

# ── Off-steady-state blocks ───────────────────────────────────────────────────

@simple
def capital_adj_F(K_F, Q_F, I_F, Z_F, N_F, alpha_F, delta_F, gamma0_F, gamma1_F, ksi_F):
    iota_F        = I_F / K_F(-1)
    mpk_F         = alpha_F * Z_F * K_F ** (alpha_F - 1) * N_F ** (1 - alpha_F)
    rk_F          = (mpk_F + (1 - delta_F) * Q_F) / Q_F(-1) - 1
    q_res_F       = Q_F - 1 / (gamma0_F * (1 - ksi_F) * iota_F ** (-ksi_F))
    capital_res_F = K_F - (1 - delta_F) * K_F(-1) - (gamma0_F * iota_F ** (1 - ksi_F) + gamma1_F) * K_F(-1)
    return iota_F, mpk_F, rk_F, q_res_F, capital_res_F

@simple
def labor_F(N_F, Z_F, K_F, alpha_F):
    Y_F = Z_F * K_F ** alpha_F * N_F ** (1 - alpha_F)
    return Y_F

@simple
def labor_market_F(w_F, N_F, vphi_F, frisch_F, P_CES_F):
    # GHH: UCE = x^(-1/eis) divides out of both sides → w/P = vphi*N^(1/frisch)
    labor_mkt_res_F = w_F / P_CES_F - vphi_F * N_F ** (1 / frisch_F)
    return labor_mkt_res_F


@simple
def labor_demand_F(w_F, Y_F, N_F, alpha_F):
    # Firm FOC: w = (1−α)·Y/N. Pins the wage in ha_full (drop labor_mkt_res_F there).
    w_res_F = w_F - (1 - alpha_F) * Y_F / N_F
    return w_res_F




@simple
def intermediation_IC_F(nu_K_F, nu_bF_F, nu_bD_F, eta_F,
                        Q_F, K_F, q_b_F, q_b_D, b_F_F, b_D_F, n_inter_F,
                        lambda_gk_F, Delta_bF_F, Delta_bD_F, theta_F, p,
                        def_rate_F, def_rate_D, psi_lambda_B_F):
    kappa_F      = Q_F   * K_F   / n_inter_F
    phi_bF_F     = q_b_F * b_F_F / (p * n_inter_F)
    phi_bD_F     = q_b_D * b_D_F / (p * n_inter_F)
    # GK multi-asset IC — see intermediation_IC_D for derivation.
    Delta_bF_eff = Delta_bF_F + psi_lambda_B_F * def_rate_F(+1)
    Delta_bD_eff = Delta_bD_F + psi_lambda_B_F * def_rate_D(+1)
    value_F      = nu_K_F * kappa_F + nu_bF_F * phi_bF_F + nu_bD_F * phi_bD_F + eta_F
    theta_tgt_F  = (value_F / lambda_gk_F
                    + (1 - Delta_bF_eff) * phi_bF_F
                    + (1 - Delta_bD_eff) * phi_bD_F)
    ic_res_F     = theta_F - theta_tgt_F
    return ic_res_F

@simple
def bank_return_F(theta_F, rk_F, rdep_F, b_F_F, b_D_F, n_inter_F,
                  rb_actual_F, rb_actual_D, q_b_F, q_b_D, p):
    phi_bF_lag_F = q_b_F(-1) * b_F_F(-1) / (p(-1) * n_inter_F(-1))
    phi_bD_lag_F = q_b_D(-1) * b_D_F(-1) / (p(-1) * n_inter_F(-1))
    kappa_lag_F  = theta_F(-1) - phi_bF_lag_F - phi_bD_lag_F
    rn_F = (kappa_lag_F  * (rk_F        - rdep_F)
            + phi_bF_lag_F * (rb_actual_F - rdep_F)
            + phi_bD_lag_F * (rb_actual_D - rdep_F)
            + rdep_F)
    return rn_F

@simple
def intermediation_P1_F(rk_F, rb_actual_F, rb_actual_D, rdep_F,
                        nu_K_F, nu_bF_F, nu_bD_F, eta_F,
                        lambda_gk_F, theta_F, SDF_banker_F, f_F):
    Omega_p1_F    = f_F + (1 - f_F) * lambda_gk_F * theta_F(+1)
    nu_K_res_F    = nu_K_F  - SDF_banker_F * Omega_p1_F * (rk_F(+1)        - rdep_F(+1))
    nu_bF_res_F   = nu_bF_F - SDF_banker_F * Omega_p1_F * (rb_actual_F(+1) - rdep_F(+1))
    nu_bD_res_F   = nu_bD_F - SDF_banker_F * Omega_p1_F * (rb_actual_D(+1) - rdep_F(+1))
    eta_res_F     = eta_F   - SDF_banker_F * Omega_p1_F * (1 + rdep_F(+1))
    return nu_K_res_F, nu_bF_res_F, nu_bD_res_F, eta_res_F

@simple
def k_balance_sheet_F(Q_F, theta_F, n_inter_F, K_F, b_F_F, b_D_F, q_b_F, q_b_D, p):
    K_res_F = Q_F * K_F + (q_b_F * b_F_F + q_b_D * b_D_F) / p - theta_F * n_inter_F
    return K_res_F


@simple
def cap_adj_cost_inter_F(K_F, rk_F, chi0_F, chi1_F, chi2_F):
    arg_F = (K_F - (1.0 + rk_F) * K_F(-1)) / (K_F(-1) + chi0_F)
    Phi_F = (chi1_F / chi2_F) * (arg_F ** 2) ** (chi2_F / 2) * (K_F(-1) + chi0_F)
    return Phi_F


@simple
def macro_pru_tax_F(b_F_F, b_D_F, def_rate_F, T0_F, T1_F, p):
    # Tax base = total bond holdings (face value in D-goods); convert to F-goods via /p.
    tau_mp_F = T0_F + T1_F * def_rate_F
    T_F      = tau_mp_F * (b_F_F + b_D_F) / p
    return tau_mp_F, T_F


@simple
def intermediation_P2_F(rn_F, n_inter_F, m_F, f_F, cap_profit_F):
    # Writedown terms removed: rb_actual already embeds the default haircut via
    # rb_actual = (1 − def·haircut)/q_b(-1) − 1, so deducting them again double-counts.
    gross_income_F = (1 + rn_F) * n_inter_F(-1) + cap_profit_F
    n_inter_val_F  = (1 - f_F) * gross_income_F + m_F - n_inter_F
    return n_inter_val_F

@simple
def banker_div_res_F(rn_F, n_inter_F, div_F, m_F, f_F, cap_profit_F, Phi_F, T_F):
    gross_income_F = (1 + rn_F) * n_inter_F(-1) + cap_profit_F
    net_div_F      = f_F * gross_income_F - m_F - Phi_F - T_F
    div_res_F      = div_F - net_div_F
    return div_res_F

@simple
def intermediation_P3_F(Q_F, K_F, n_inter_F, b_F_F, b_D_F, q_b_F, q_b_D, p):
    D_supply_F = Q_F * K_F + (q_b_F * b_F_F + q_b_D * b_D_F) / p - n_inter_F
    return D_supply_F

@simple
def bond_price_ss_F(SDF_banker_F, def_rate_F, recovery_rate_F, delta_b_F, zeta_writeoff_F, writeoff_enabled_F):
    haircut_F      = 1.0 - recovery_rate_F
    haircut_mult_F = writeoff_enabled_F
    surv_cont_F    = 1.0 - zeta_writeoff_F * def_rate_F * haircut_F * haircut_mult_F
    q_b_F          = (
        SDF_banker_F * delta_b_F * (1.0 - def_rate_F * haircut_F * haircut_mult_F)
        / (1.0 - SDF_banker_F * (1.0 - delta_b_F) * surv_cont_F)
    )
    return q_b_F


@simple
def domestic_bond_foc_F(rb_actual_F, rdep_F, b_F_F, n_inter_F, q_b_F,
                         phi_bF_F_ss, psi_bF_F, excess_return_bF_F_ss, tau_mp_F, p):
    phi_bF_F     = q_b_F * b_F_F / (p * n_inter_F)
    # Expected F-good return on D-good bond: (1+rb)·p/p(+1) − 1
    rb_F_fg_next = (1 + rb_actual_F(+1)) * p / p(+1) - 1
    rb_F_res     = (rb_F_fg_next - rdep_F(+1)) - excess_return_bF_F_ss \
                   - psi_bF_F * (phi_bF_F - phi_bF_F_ss) \
                   - tau_mp_F
    return rb_F_res


@simple
def government_default_F(shock_def_F, b_gov_F, Y_ss_F, b_gov_ss_F,
                          def_scale_F, def_curvature_F, def_offset_F):
    debt_ratio_F = b_gov_F(-1) / Y_ss_F
    ss_ratio_F   = b_gov_ss_F  / Y_ss_F
    def_rate_F   = shock_def_F + def_scale_F * (
        (debt_ratio_F + def_offset_F) ** def_curvature_F
      - (ss_ratio_F  + def_offset_F) ** def_curvature_F
    )
    return def_rate_F

@simple
def tax_rule_F(b_gov_F, b_gov_ss_F, phi_lamb_F):
    T_ls_F = phi_lamb_F * (b_gov_F(-1) - b_gov_ss_F)
    return T_ls_F


@simple
def capital_producer_profit_F(Q_F, K_F, I_F, delta_F):
    cap_profit_F = Q_F * (K_F - (1 - delta_F) * K_F(-1)) - I_F
    return cap_profit_F

@simple
def budget_residual_F(b_gov_F, G_F, TAX_F, q_b_F, def_rate_F, recovery_rate_F, zeta_writeoff_F, p, P_CES_F, delta_b_F, writeoff_enabled_F):
    haircut_F      = 1.0 - recovery_rate_F
    haircut_mult_F = writeoff_enabled_F
    surv_cont_F    = 1.0 - zeta_writeoff_F * def_rate_F * haircut_F * haircut_mult_F
    coupon_F       = delta_b_F * (1.0 - def_rate_F * haircut_F * haircut_mult_F) * b_gov_F(-1)
    net_issuance_F = q_b_F * (b_gov_F - surv_cont_F * (1.0 - delta_b_F) * b_gov_F(-1))
    b_gov_res_F    = (coupon_F - net_issuance_F) / p + G_F - P_CES_F * TAX_F
    return b_gov_res_F


@simple
def divert_bond_foc_F(rb_actual_F, rdep_F, b_F_F, n_inter_F, q_b_F,
                      phi_bF_F_ss, psi_bF_F, excess_return_bF_F_ss, tau_mp_F, p,
                      psi_spread_F, def_rate_F):
    phi_bF_F   = q_b_F * b_F_F / (p * n_inter_F)
    # IC-theory derived required spread: additive default loading independent of SS excess return.
    # psi_spread_F = lambda_gk_F * psi_lambda_B_F / (beta_inter_F * Omega_F), computed in _apply_ss_anchors.
    req_spread = excess_return_bF_F_ss + psi_spread_F * def_rate_F(+1)
    rb_F_res   = (rb_actual_F(+1) - rdep_F(+1)) - req_spread \
                 - psi_bF_F * (phi_bF_F - phi_bF_F_ss) \
                 - tau_mp_F
    return rb_F_res


@simple
def welfare_agg_F(X_F, C_F_ss):
    # GHH utility composite normalised by SS consumption.
    # In IRFs the deviation ΔX_F/C_F_ss gives welfare change as a fraction of SS consumption.
    U_F = X_F / C_F_ss
    return U_F
