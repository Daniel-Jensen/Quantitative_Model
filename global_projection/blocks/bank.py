# TWO-COUNTRY GERTLER-KARADI / BOCOLA (2016) FINANCIAL-INTERMEDIARY BLOCK.
# Each bank holds capital + domestic + foreign bonds. Cross-border legs convert
# via p (D-goods per F-good), so p up => F-goods more expensive.
# Kernel convention (Bocola Prop. 1): f = exit/payout share per period, so
# Omega = beta*[f + (1-f)*alpha'] puts weight (1-f) on the franchise value.
# The IC multiplier is occasionally binding; the dynamic complementarity is Bocola's
# closed form in solver_recursive/point_map.py. It binds with equality at the SS.
import numpy as np
from scipy.optimize import brentq

def _least_root(resid, what, v_lo=1e-6, v_hi=1e6, n_scan=300):
    # LEAST ROOT OF A SCALAR RESIDUAL, FOUND BY LOG-GRID SIGN SCAN + BRENTQ.
    # Selecting the LEAST root matters: the franchise fixed point can fold, and
    # value iteration from below would land on this root too (see calibration.py).
    grid = np.geomspace(v_lo, v_hi, n_scan)
    vals = np.array([resid(v) for v in grid])
    fin  = np.isfinite(vals)
    sc   = np.where(np.diff(np.sign(vals[fin])) != 0)[0]
    if len(sc) == 0:
        raise RuntimeError(f"No sign change in the alpha fixed point ({what}).")
    gf = grid[fin]
    i  = sc[0]
    return brentq(resid, gf[i], gf[i + 1], xtol=1e-13, rtol=1e-13)


def _alpha_ss_fixed_point(beta_inter, f, lambda_K, rk_ss, rdep_ss):
    # SCALAR FIXED POINT FOR THE FRANCHISE VALUE alpha (AND mu, Omega) AT THE SS.
    def resid(a):
        # BELLMAN GAP AT A CANDIDATE FRANCHISE VALUE a.
        Omega = beta_inter * (f + (1 - f) * a)
        mu    = Omega * (rk_ss - rdep_ss) / lambda_K
        if mu >= 1.0:
            return np.inf
        return Omega * (1 + rdep_ss) / (1 - mu) - a

    alpha_ss = _least_root(resid, f"rk_ss={rk_ss:.6f}, rdep_ss={rdep_ss:.6f}")
    Omega_ss = beta_inter * (f + (1 - f) * alpha_ss)
    mu_ss    = Omega_ss * (rk_ss - rdep_ss) / lambda_K
    return alpha_ss, mu_ss, Omega_ss


def calibrate_bank_targets(beta_inter, f, rdep, theta_target, spread_target):
    # SOLVE THE SINGLE lambda AND ENTRANT TRANSFER omega_ent FROM LEVERAGE + SPREAD TARGETS.
    s = spread_target

    def resid(a):
        # BELLMAN GAP AT A CANDIDATE FRANCHISE VALUE a, WITH lambda = alpha/theta
        # FOLDED IN SO alpha = Omega(1+rdep)/(1-mu) STAYS A SCALAR FIXED POINT.
        Omega = beta_inter * (f + (1 - f) * a)
        mu    = Omega * s * theta_target / a
        if mu >= 1.0:
            return np.inf
        return Omega * (1 + rdep) / (1 - mu) - a

    alpha = _least_root(resid, f"theta={theta_target}, spread={s}")
    Omega = beta_inter * (f + (1 - f) * alpha)
    mu    = Omega * s * theta_target / alpha
    lambda_single = alpha / theta_target

    D_val     = 1.0 - (1 - f) * (1 + rdep)   # net-worth accumulation discount
    omega_ent = D_val / theta_target - (1 - f) * s
    if omega_ent <= 0.0:
        raise RuntimeError(
            f"Infeasible targets: omega_ent={omega_ent:.4e} <= 0 "
            f"(theta={theta_target}, spread={s}, f={f}, rdep={rdep})."
        )
    return lambda_single, omega_ent, alpha, mu, Omega


def steady_state_bank(cal, rk_ss, Kap_ss, Q_bD_ss, Q_bF_ss,
                      b_dom_ss, b_for_ss, p_ss, country="D", L_wc_ss=0.0):
    # STEADY-STATE BANK BLOCK: PRICES, MULTIPLIERS, NET WORTH, LEVERAGE, DEPOSITS.
    # L_wc_ss is the WORKING-CAPITAL LOAN (Bocola 2016 SV.C): "firms need to borrow a
    # fraction phi of the wage bill before production takes place. These loans are
    # obtained from the bankers ... and they pay the gross return R_W = R* + lambda*
    # mu/E[Lambda]". So the loan is a BANK ASSET, deposit-funded, inside the divertable
    # base at the SAME lambda (his residual_model_open.m: mu = N/(lambda*(Q*K + q*B +
    # psi*(1-alpha)*gdp))), and the interest accrues to the BANK -- not, as before, to
    # households as a dividend, which made the credit spread a pure intra-period transfer
    # and turned the risk channel expansionary.
    # It earns the SAME excess return s = lambda*mu/Omega as the other classes, so the
    # leverage identity theta = D_val/((1-f)s + omega_ent) is unchanged and the
    # calibration still hits its targets with the enlarged asset base.
    f          = cal[f"f_{country}"]
    rdep_ss    = cal[f"r_dep_{country}_target"]
    lambda_K   = cal[f"lambda_K_{country}"]
    lambda_bD  = cal[f"lambda_bD_{country}"]
    lambda_bF  = cal[f"lambda_bF_{country}"]
    omega_ent  = cal[f"omega_ent_{country}"]

    alpha_ss, mu_ss, Omega_ss = _alpha_ss_fixed_point(
        cal[f"beta_inter_{country}"], f, lambda_K, rk_ss, rdep_ss
    )
    px = _ss_bond_prices(cal, country, rdep_ss, mu_ss, Omega_ss,
                         lambda_K, lambda_bD, lambda_bF)
    D_val = 1.0 - (1 - f) * (1 + rdep_ss)
    if D_val <= 0:
        raise ValueError(f"[{country}] D={D_val} <= 0: no stationary net-worth rest point.")
    n_ss_IC, total_assets, n_ss_ACCUM = _ss_net_worth(
        country, f, rk_ss, rdep_ss, omega_ent, lambda_K, lambda_bD, lambda_bF, alpha_ss,
        D_val, Kap_ss, b_dom_ss, b_for_ss, p_ss, L_wc_ss, px)
    n_ss = n_ss_ACCUM
    if n_ss <= 0:
        raise ValueError(f"[{country}] n_ss={n_ss:.4f} <= 0 at rk_ss={rk_ss:.6f}.")
    pf = _ss_portfolio(country, f, rk_ss, rdep_ss, omega_ent, n_ss, total_assets,
                       Kap_ss, b_dom_ss, b_for_ss, p_ss, L_wc_ss, px)

    return dict(
        alpha_ss=alpha_ss, mu_ss=mu_ss, Omega_ss=Omega_ss,
        n_ss=n_ss, n_ss_IC=n_ss_IC, n_ss_ACCUM=n_ss_ACCUM,
        kappa_ss=pf["kappa_ss"], phi_bdom_ss=pf["phi_bdom_ss"],
        phi_bfor_ss=pf["phi_bfor_ss"], theta_ss=pf["theta_ss"], div_ss=pf["div_ss"],
        phi_wc_ss=pf["phi_wc_ss"], L_wc_ss=L_wc_ss, IC_spread_wc=px["IC_spread_wc"],
        Dep_supply_ss=(pf["theta_ss"] - 1) * n_ss,
        # the P STATE is the bank's obligation NET of the working-capital receivable
        # (Bocola: P' = R*(assets - N') - R_W*L), which is what makes ng = X - P hold
        # without carrying L as a separate state
        P_state_ss=((1 + rdep_ss) * ((pf["theta_ss"] - 1) * n_ss)
                    - (1 + rdep_ss + px["IC_spread_wc"]) * L_wc_ss),
        rb_dom_ss=px["rb_dom_ss"], rb_for_ss=px["rb_for_ss"],
        Q_bdom_IC=px["Q_bdom_ss"],
        IC_spread_dom=px["IC_spread_dom"], IC_spread_for=px["IC_spread_for"],
        lambda_K=lambda_K, lambda_bD=lambda_bD, lambda_bF=lambda_bF,
    )


def _ss_bond_prices(cal, country, rdep_ss, mu_ss, Omega_ss, lambda_K, lambda_bD, lambda_bF):
    # IC-REQUIRED EXCESS RETURNS lambda*mu/Omega AND THE HM PERPETUITY PRICES THEY IMPLY.
    IC_spread_dom = lambda_bD * mu_ss / Omega_ss
    IC_spread_for = lambda_bF * mu_ss / Omega_ss
    IC_spread_wc  = lambda_K * mu_ss / Omega_ss      # r_wc - rdep, Bocola's R_W - R*
    # "dom"/"for" are relative to the country: D's home bond is the D-bond
    delta_b_D, delta_b_F = cal["delta_b_D"], cal["delta_b_F"]
    db_dom, db_for = (delta_b_D, delta_b_F) if country == "D" else (delta_b_F, delta_b_D)
    return dict(IC_spread_dom=IC_spread_dom, IC_spread_for=IC_spread_for,
                IC_spread_wc=IC_spread_wc,
                Q_bdom_ss=db_dom / (rdep_ss + db_dom + IC_spread_dom),
                Q_bfor_ss=db_for / (rdep_ss + db_for + IC_spread_for),
                rb_dom_ss=rdep_ss + IC_spread_dom, rb_for_ss=rdep_ss + IC_spread_for)


def _ss_net_worth(country, f, rk_ss, rdep_ss, omega_ent, lambda_K, lambda_bD, lambda_bF,
                  alpha_ss, D_val, Kap_ss, b_dom_ss, b_for_ss, p_ss, L_wc_ss, px):
    # NET WORTH FROM THE BINDING IC (n_IC) AND FROM ACCUMULATION (n_ACCUM), AND ASSETS.
    # Foreign leg valued in home goods: D holds F-bonds (xp), F holds D-bonds (/p).
    # The operand order below is load-bearing, not stylistic — the stage-1 SS solve
    # runs to xtol 1e-12 and the TPI Q-floor solves sit close enough to a Newton
    # knife edge that even a one-ulp reassociation here can flip one of them from
    # converged to stalled. Keep the conversion inside each product.
    Q_bdom_ss, Q_bfor_ss = px["Q_bdom_ss"], px["Q_bfor_ss"]
    IC_spread_dom, IC_spread_for = px["IC_spread_dom"], px["IC_spread_for"]
    IC_spread_wc = px["IC_spread_wc"]
    bdom_val = Q_bdom_ss * b_dom_ss
    if country == "D":
        n_ss_IC = (lambda_K * Kap_ss
                   + lambda_bD * Q_bdom_ss * b_dom_ss
                   + lambda_bF * p_ss * Q_bfor_ss * b_for_ss
                   + lambda_K * L_wc_ss) / alpha_ss
        total_assets = Kap_ss + bdom_val + p_ss * Q_bfor_ss * b_for_ss + L_wc_ss
        n_ss_ACCUM = (
            ((1 - f) * (rk_ss - rdep_ss) + omega_ent) * Kap_ss
            + ((1 - f) * IC_spread_dom + omega_ent) * Q_bdom_ss * b_dom_ss
            + ((1 - f) * IC_spread_for + omega_ent) * p_ss * Q_bfor_ss * b_for_ss
            + ((1 - f) * IC_spread_wc + omega_ent) * L_wc_ss
        ) / D_val
    else:
        bfor_val = Q_bfor_ss * b_for_ss / p_ss
        n_ss_IC = (lambda_K * Kap_ss
                   + lambda_bD * Q_bdom_ss * b_dom_ss
                   + lambda_bF * Q_bfor_ss * b_for_ss / p_ss
                   + lambda_K * L_wc_ss) / alpha_ss
        total_assets = Kap_ss + bdom_val + bfor_val + L_wc_ss
        n_ss_ACCUM = (
            ((1 - f) * (rk_ss - rdep_ss) + omega_ent) * Kap_ss
            + ((1 - f) * IC_spread_dom + omega_ent) * bdom_val
            + ((1 - f) * IC_spread_for + omega_ent) * bfor_val
            + ((1 - f) * IC_spread_wc + omega_ent) * L_wc_ss
        ) / D_val
    return n_ss_IC, total_assets, n_ss_ACCUM


def _ss_portfolio(country, f, rk_ss, rdep_ss, omega_ent, n_ss, total_assets,
                  Kap_ss, b_dom_ss, b_for_ss, p_ss, L_wc_ss, px):
    # PORTFOLIO SHARES OF NET WORTH, LEVERAGE, THE NET RETURN AND BANK DIVIDENDS.
    Q_bdom_ss, Q_bfor_ss = px["Q_bdom_ss"], px["Q_bfor_ss"]
    kappa_ss    = Kap_ss / n_ss
    phi_bdom_ss = Q_bdom_ss * b_dom_ss / n_ss
    phi_bfor_ss = (p_ss * Q_bfor_ss * b_for_ss / n_ss if country == "D"
                   else Q_bfor_ss * b_for_ss / (p_ss * n_ss))
    phi_wc_ss   = L_wc_ss / n_ss
    theta_ss    = kappa_ss + phi_bdom_ss + phi_bfor_ss + phi_wc_ss
    rn_ss = (kappa_ss * (rk_ss - rdep_ss)
             + phi_bdom_ss * (px["rb_dom_ss"] - rdep_ss)
             + phi_bfor_ss * (px["rb_for_ss"] - rdep_ss)
             + phi_wc_ss * px["IC_spread_wc"]
             + rdep_ss)
    div_ss = f * (1 + rn_ss) * n_ss - omega_ent * total_assets
    return dict(kappa_ss=kappa_ss, phi_bdom_ss=phi_bdom_ss, phi_bfor_ss=phi_bfor_ss,
                phi_wc_ss=phi_wc_ss, theta_ss=theta_ss, div_ss=div_ss)


