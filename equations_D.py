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
                   rb_actual_D, rb_actual_F, b_D_D, b_F_D, Q_D):
    K_D          = theta_D * n_inter_D
    phi_bD_D     = b_D_D / n_inter_D
    phi_bF_D     = b_F_D / n_inter_D          
    rk_D         = alpha_D * Y_D / K_D - delta_D
    rn_D         = theta_D * (rk_D - rdep_D) + phi_bD_D * (rb_actual_D - rdep_D) + phi_bF_D * (rb_actual_F - rdep_D) + rdep_D
    m_D          = n_inter_D * (1 - (1 - f_D) * (1 + rn_D))
    k_inter_D    = K_D
    I_D          = K_D * delta_D
    D_supply_D   = (theta_D - 1) * n_inter_D + b_D_D + b_F_D 
    Z_D          = Y_D / ((K_D ** alpha_D) * (N_D ** (1 - alpha_D)))
    rdep_ante_D  = rdep_D
    cap_profit_D = Q_D * (K_D - (1 - delta_D) * K_D(-1)) - I_D
    return K_D, rk_D, rn_D, m_D, k_inter_D, I_D, D_supply_D, Z_D, rdep_ante_D, cap_profit_D

@simple
def market_clearing_D(Y_D, C_D, I_D, G_D, NX_D, DEP_D, D_supply_D):
    goods_mkt_D   = Y_D - C_D - I_D - G_D - NX_D
    deposit_mkt_D = DEP_D - D_supply_D
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
def steady_auxilliary_D(theta_D, rk_D, rdep_D, delta_D, alpha_D, Y_D, K_D, N_D, lambda_gk_D, beta_D, ksi_D, rn_D):
    iota_D   = delta_D
    mpk_D    = alpha_D * (Y_D / K_D)
    w_D      = (1 - alpha_D) * Y_D / N_D
    Omega_D  = theta_D * lambda_gk_D / (beta_D * (1 + rn_D))
    nu_D     = beta_D * Omega_D * (rk_D - rdep_D)
    eta_D    = beta_D * Omega_D * (1 + rdep_D)
    gamma0_D = delta_D ** ksi_D / (1 - ksi_D)
    gamma1_D = -delta_D * ksi_D / (1 - ksi_D)
    return iota_D, mpk_D, w_D, Omega_D, nu_D, eta_D, gamma0_D, gamma1_D


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
def government_ss_D(TAX_D, rb_D, b_gov_D):
    # Use promised yield rb_D — consistent with budget_residual_D (interest rate channel).
    G_D = TAX_D - rb_D * b_gov_D
    return G_D


@simple
def labor_ss_D(w_D, N_D, UCE_D, frisch_D, mu_w_D):
    # Calibrate vphi_D so wage NKPC residual = 0 at SS:
    # vphi_D * N^(1/φ) = (1/μ_w) * w * UCE   → with μ_w=1 this is the frictionless MRS.
    vphi_D = (1 / mu_w_D) * w_D * UCE_D / (N_D ** (1 / frisch_D))
    return vphi_D

@simple
def bond_return_D(rb_D, def_rate_D, recovery_rate_D):
    haircut_D   = 1.0 - recovery_rate_D
    rb_actual_D = (1 - def_rate_D * haircut_D) * (1 + rb_D(-1)) - 1
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
def portfolio_foc_bF_D(rb_actual_F, rdep_D, b_F_D, n_inter_D, p,
                       phi_bF_D_ss, psi_bF_D, excess_return_F_D_ss):
    phi_bF_D            =  b_F_D / n_inter_D
    rb_actual_F_in_D_p1 = (1 + rb_actual_F) - 1
    foc_bF_res_D        = (rb_actual_F_in_D_p1 - rdep_D(+1)) - excess_return_F_D_ss \
                          - psi_bF_D * (phi_bF_D - phi_bF_D_ss)
    return foc_bF_res_D


@simple
def macroprudential_D(def_rate_D, lambda_gk_D, phi_macro_D):
    # Countercyclical tightening of the GK IC constraint via lambda_eff_D.
    # Higher lambda_eff → regulator treats more bank assets as divertable →
    # IC binds more tightly → theta = eta/(lambda_eff - nu) falls →
    # banks hold less capital per unit of equity → K_D↓ → I_D↓ → Y_D↓.
    # At SS (def_rate_D = 0): lambda_eff_D = lambda_gk_D — no SS distortion.
    lambda_eff_D = lambda_gk_D + phi_macro_D * def_rate_D
    return lambda_eff_D


@simple
def intermediation_IC_D(nu_D, eta_D, lambda_eff_D):
    theta_D = eta_D / (lambda_eff_D - nu_D)
    return theta_D


@simple
def bank_return_D(theta_D, rk_D, rdep_D, b_D_D, b_F_D, n_inter_D, rb_actual_D, rb_actual_F,
                  phi_bF_D_ss, psi_bF_D):
    phi_bD_lag_D = b_D_D(-1) / n_inter_D(-1)
    phi_bF_lag_D = b_F_D(-1) / n_inter_D(-1)

    rn_D = (theta_D(-1) * (rk_D - rdep_D)
            + phi_bD_lag_D * (rb_actual_D - rdep_D)
            + phi_bF_lag_D * (rb_actual_F - rdep_D)
            + rdep_D)

    rn_D = rn_D - (psi_bF_D / 2) * (phi_bF_lag_D - phi_bF_D_ss) ** 2
    return rn_D


@simple
def intermediation_P1_D(rk_D, rdep_D, nu_D, lambda_eff_D, eta_D, theta_D, SDF_D, f_D):
    # lambda_eff_D(+1): use expected future regulatory lambda in bank value function,
    # consistent with GK timing (IC constraint applies at t+1).
    Omega_p1_D = f_D + (1 - f_D) * lambda_eff_D(+1) * theta_D(+1)
    nu_res_D   = nu_D  - SDF_D * Omega_p1_D * (rk_D(+1) - rdep_D(+1))
    eta_res_D  = eta_D - SDF_D * Omega_p1_D * (1 + rdep_D(+1))
    return nu_res_D, eta_res_D


@simple
def k_balance_sheet_D(Q_D, theta_D, n_inter_D, K_D):
    # GK balance sheet constraint restored: Q*K = theta*n_inter.
    # Macroprudential tightening enters via lambda_eff_D → theta_D (see macroprudential_D).
    K_res_D = Q_D * K_D - theta_D * n_inter_D
    return K_res_D


@simple
def firm_profit_D(mc_D, Y_D):
    # Monopoly profit from sticky-price markup: (1 - mc_D) * Y_D.
    # Zero at SS (mc_D = 1), first-order off-SS. Routed to the financial
    # intermediary (banks own goods-producing firms in this GK setup).
    firm_profit_D = (1 - mc_D) * Y_D
    return firm_profit_D


@simple
def intermediation_P2_D(rn_D, n_inter_D, m_D, f_D, cap_profit_D, firm_profit_D):
    gross_income_D = (1 + rn_D) * n_inter_D(-1) + cap_profit_D + firm_profit_D
    n_inter_val_D  = (1 - f_D) * gross_income_D + m_D - n_inter_D
    return n_inter_val_D


@simple
def banker_div_res_D(rn_D, n_inter_D, div_D, m_D, f_D, cap_profit_D, firm_profit_D):
    gross_income_D = (1 + rn_D) * n_inter_D(-1) + cap_profit_D + firm_profit_D
    net_div_D      = f_D * gross_income_D - m_D
    div_res_D      = div_D - net_div_D
    return div_res_D


@simple
def intermediation_P3_D(Q_D, K_D, n_inter_D, b_D_D, b_F_D):
    # MU: b_F_D already in common currency — no p conversion.
    D_supply_D = Q_D * K_D + b_D_D + b_F_D - n_inter_D
    return D_supply_D


@simple
def interest_rates_D(def_rate_D, recovery_rate_D, SDF_D):
    # Used in SS model only. In ha_full, rb_D is an unknown pinned by domestic_bond_foc_D.
    haircut_D = 1 - def_rate_D(+1) * (1.0 - recovery_rate_D)
    rb_D      = 1 / (SDF_D * haircut_D) - 1
    return rb_D


@simple
def domestic_bond_foc_D(rb_actual_D, rdep_D, b_D_D, n_inter_D,
                         phi_bD_D_ss, psi_bD_D, excess_return_bD_D_ss):
    # GK-consistent bond pricing: rb_D adjusts until D-banks willingly hold
    # b_D_D (the residual after foreign banks take b_D_F).
    # At SS: phi_bD_D = phi_bD_D_ss → rb_D_res = 0 by construction.
    # After default shock: b_D_D rises (foreign flight) → phi_bD_D rises → rb_D rises.
    phi_bD_D = b_D_D / n_inter_D
    rb_D_res = (rb_actual_D(+1) - rdep_D(+1)) - excess_return_bD_D_ss \
               - psi_bD_D * (phi_bD_D - phi_bD_D_ss)
    return rb_D_res


# ==> GOVERMENT EQUATIONS
@simple
def government_default_D(shock_def_D):
    def_rate_D = shock_def_D
    return def_rate_D


@simple
def tax_rule_D(b_gov_D, lamb_ss_D, b_gov_ss_D, phi_lamb_D):
    lamb_D = lamb_ss_D - phi_lamb_D * (b_gov_D(-1) - b_gov_ss_D)
    return lamb_D


@simple
def budget_residual_D(b_gov_D, G_D, TAX_D, rb_D, def_rate_D, recovery_rate_D, zeta_writeoff_D):
    # zeta_writeoff_D ∈ [0, 1] dials between two default regimes:
    #   ζ = 1  (default) Treasury repays only the post-haircut amount
    #          rb_actual_D_implied; the haircut is a fiscal transfer from
    #          bondholders to the Treasury. Matches the v11 design and the
    #          empirical Eurozone-restructuring narrative (sovereign captures
    #          real fiscal relief; trade flows shift through internal devaluation).
    #   ζ = 0  Treasury repays the full promised yield rb_D(-1); the haircut is
    #          pure deadweight (no agent captures it). Bondholders still take
    #          the loss via rb_actual_D in bond_return_D. Matches PR #2 design.
    # At SS def_rate_D = 0 ⇒ rb_actual_D_implied = rb_D, so the residual is
    # invariant to ζ at SS — no recalibration needed when ζ moves.
    haircut_D             = 1.0 - recovery_rate_D
    rb_actual_D_implied   = (1 - def_rate_D * haircut_D) * (1 + rb_D(-1)) - 1
    effective_repayment_D = (
        zeta_writeoff_D         * (1 + rb_actual_D_implied) * b_gov_D(-1)
        + (1 - zeta_writeoff_D) * (1 + rb_D(-1))            * b_gov_D(-1)
    )
    b_gov_res_D           = effective_repayment_D + G_D - TAX_D - b_gov_D
    return b_gov_res_D


# ── NK frictions (price NKPC, wage NKPC, Taylor rule, Fisher) ────────────────

@simple
def fisher_D(i_union, pi_D):
    # MU: both countries share the single nominal rate i_union set by the union CB.
    # Real rates differ across countries via inflation differentials.
    rdep_ante_D = i_union - pi_D(+1)
    rdep_D      = i_union(-1) - pi_D
    return rdep_ante_D, rdep_D

@simple
def taylor_rule_union(i_union, pi_D, pi_F, Y_D, Y_F, Y_ss_D, Y_ss_F,
                      omega_union_U, rho_i_U, phi_pi_U, phi_y_U, r_star_U, eps_m_D):
    pi_union         = omega_union_U * pi_D + (1 - omega_union_U) * pi_F
    Y_union          = omega_union_U * Y_D  + (1 - omega_union_U) * Y_F
    Y_ss_union       = omega_union_U * Y_ss_D + (1 - omega_union_U) * Y_ss_F
    taylor_res_union = (i_union
                        - rho_i_U * i_union(-1)
                        - (1 - rho_i_U) * (r_star_U + phi_pi_U * pi_union
                                           + phi_y_U * (Y_union / Y_ss_union - 1))
                        - eps_m_D)
    return taylor_res_union

@simple
def pricing_D(mc_D, pi_D, kappa_p_D, mu_p_D, SDF_D):
    nkpc_p_res_D = pi_D - kappa_p_D * (mc_D - 1 / mu_p_D) - SDF_D * pi_D(+1)
    return nkpc_p_res_D

@simple
def wage_setting_D(w_D, pi_w_D, pi_D, N_D, UCE_D, vphi_D, frisch_D, kappa_w_D, mu_w_D, beta_D):
    # Rotemberg wage NKPC: π_w = κ_w * (MRS - w/μ_w) + β E[π_w(+1)].
    # vphi_D calibrated so MRS = w/μ_w at SS → residual = 0.
    nkpc_w_res_D = (pi_w_D
                    - kappa_w_D * (vphi_D * N_D ** (1 + 1 / frisch_D)
                                   - (1 / mu_w_D) * w_D * N_D * UCE_D)
                    - beta_D * pi_w_D(+1))
    # Real wage law of motion: w_t = w_{t-1} * (1 + π^w_t) / (1 + π_t).
    w_res_D = w_D - w_D(-1) * (1 + pi_w_D) / (1 + pi_D)
    return nkpc_w_res_D, w_res_D
