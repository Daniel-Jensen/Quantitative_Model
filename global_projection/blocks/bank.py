# TWO-COUNTRY GERTLER-KARADI / BOCOLA (2016) FINANCIAL-INTERMEDIARY BLOCK.
# Each bank holds capital + domestic + foreign bonds. Cross-border legs convert
# via p (D-goods per F-good), so p up => F-goods more expensive.
# Kernel convention (Bocola Prop. 1): f = exit/payout share per period, so
# Omega = beta*[f + (1-f)*alpha'] puts weight (1-f) on the franchise value.
# The IC multiplier is occasionally binding: mu = Omega*E[rk'-rdep]/lambda from
# the capital FOC is valid in both regimes, and the complementarity itself is
# imposed in transition.py (Fischer-Burmeister). It binds with equality at the SS.
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
    beta_inter = cal[f"beta_inter_{country}"]
    lambda_K   = cal[f"lambda_K_{country}"]
    lambda_bD  = cal[f"lambda_bD_{country}"]
    lambda_bF  = cal[f"lambda_bF_{country}"]
    omega_ent  = cal[f"omega_ent_{country}"]
    delta_b_D  = cal["delta_b_D"]
    delta_b_F  = cal["delta_b_F"]

    alpha_ss, mu_ss, Omega_ss = _alpha_ss_fixed_point(
        beta_inter, f, lambda_K, rk_ss, rdep_ss
    )
    IC_spread_dom = lambda_bD * mu_ss / Omega_ss
    IC_spread_for = lambda_bF * mu_ss / Omega_ss
    IC_spread_wc  = lambda_K * mu_ss / Omega_ss      # r_wc - rdep, Bocola's R_W - R*

    # "dom"/"for" are relative to the country: D's home bond is the D-bond
    db_dom, db_for = (delta_b_D, delta_b_F) if country == "D" else (delta_b_F, delta_b_D)
    Q_bdom_ss = db_dom / (rdep_ss + db_dom + IC_spread_dom)
    Q_bfor_ss = db_for / (rdep_ss + db_for + IC_spread_for)

    rb_dom_ss = rdep_ss + IC_spread_dom
    rb_for_ss = rdep_ss + IC_spread_for

    D_val = 1.0 - (1 - f) * (1 + rdep_ss)
    if D_val <= 0:
        raise ValueError(f"[{country}] D={D_val} <= 0: no stationary net-worth rest point.")

    # Foreign leg valued in home goods: D holds F-bonds (xp), F holds D-bonds (/p).
    # The operand order below is load-bearing, not stylistic — the stage-1 SS solve
    # runs to xtol 1e-12 and the TPI Q-floor solves sit close enough to a Newton
    # knife edge that even a one-ulp reassociation here can flip one of them from
    # converged to stalled. Keep the conversion inside each product.
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

    n_ss = n_ss_ACCUM
    if n_ss <= 0:
        raise ValueError(f"[{country}] n_ss={n_ss:.4f} <= 0 at rk_ss={rk_ss:.6f}.")

    kappa_ss    = Kap_ss / n_ss
    phi_bdom_ss = Q_bdom_ss * b_dom_ss / n_ss
    phi_bfor_ss = (p_ss * Q_bfor_ss * b_for_ss / n_ss if country == "D"
                   else Q_bfor_ss * b_for_ss / (p_ss * n_ss))
    phi_wc_ss   = L_wc_ss / n_ss
    theta_ss    = kappa_ss + phi_bdom_ss + phi_bfor_ss + phi_wc_ss

    rn_ss = (kappa_ss * (rk_ss - rdep_ss)
             + phi_bdom_ss * (rb_dom_ss - rdep_ss)
             + phi_bfor_ss * (rb_for_ss - rdep_ss)
             + phi_wc_ss * IC_spread_wc
             + rdep_ss)
    div_ss = f * (1 + rn_ss) * n_ss - omega_ent * total_assets

    return dict(
        alpha_ss=alpha_ss, mu_ss=mu_ss, Omega_ss=Omega_ss,
        n_ss=n_ss, n_ss_IC=n_ss_IC, n_ss_ACCUM=n_ss_ACCUM,
        kappa_ss=kappa_ss, phi_bdom_ss=phi_bdom_ss, phi_bfor_ss=phi_bfor_ss,
        theta_ss=theta_ss, div_ss=div_ss, phi_wc_ss=phi_wc_ss,
        L_wc_ss=L_wc_ss, IC_spread_wc=IC_spread_wc,
        Dep_supply_ss=(theta_ss - 1) * n_ss,
        # the P STATE is the bank's obligation NET of the working-capital receivable
        # (Bocola: P' = R*(assets - N') - R_W*L), which is what makes ng = X - P hold
        # without carrying L as a separate state
        P_state_ss=((1 + rdep_ss) * ((theta_ss - 1) * n_ss)
                    - (1 + rdep_ss + IC_spread_wc) * L_wc_ss),
        rb_dom_ss=rb_dom_ss, rb_for_ss=rb_for_ss,
        Q_bdom_IC=Q_bdom_ss,
        IC_spread_dom=IC_spread_dom, IC_spread_for=IC_spread_for,
        lambda_K=lambda_K, lambda_bD=lambda_bD, lambda_bF=lambda_bF,
    )


def bank_backward(rk_D, rk_F, rdep_D, rdep_F, p_path,
                  cal, ss_bk_D, ss_bk_F,
                  def_price_D=None, risk_D=None,
                  tpi_D=None, s_tpi_D=None):
    # BACKWARD PASS: VALUE SLOPES, BOND PRICES, CROSS-BORDER FOC HOLDINGS.
    # Prices come from marginal conditions only (no stocks) — that is what lets
    # debt be forward-integrated afterwards. Only D is default-risky.
    #   risk_D=None    -> risk-neutral pricing of def_price_D
    #   risk_D=dict    -> Bocola two-branch expectations (pi replaces def_price_D)
    #   tpi_D=dict     -> adds a THIRD branch, the CB backstop reneging; requires
    #                     risk_D, since TPI pricing rides on the two-branch kernel
    #   s_tpi_D        -> REALIZED purchase regime, independent of tpi_D's priced
    #                     channel (the mechanical Q-floor below reads only this)
    T = len(rk_D)

    if def_price_D is None:
        def_price_D = np.zeros(T)

    f_D        = cal["f_D"];          f_F        = cal["f_F"]
    bi_D       = cal["beta_inter_D"]; bi_F       = cal["beta_inter_F"]
    lK_D       = cal["lambda_K_D"];   lK_F       = cal["lambda_K_F"]
    lbD_D      = cal["lambda_bD_D"];  lbD_F      = cal["lambda_bD_F"]
    lbF_D      = cal["lambda_bF_D"];  lbF_F      = cal["lambda_bF_F"]
    db_D       = cal["delta_b_D"];    db_F       = cal["delta_b_F"]
    psi_bFD    = cal["psi_bF_D"]
    psi_bDF    = cal["psi_bD_F"]
    b_F_D_ss   = cal["b_F_D_ss"]
    b_D_F_ss   = cal["b_D_F_ss"]
    exc_FD_ss  = cal["excess_return_F_D_ss"]
    exc_DF_ss  = cal["excess_return_D_F_ss"]
    rec_D      = cal["recovery_rate_D"]

    alpha_D_path = np.empty(T);  mu_D_path = np.empty(T)
    alpha_F_path = np.empty(T);  mu_F_path = np.empty(T)
    Omega_D_path = np.empty(T);  Omega_F_path = np.empty(T)
    Q_bD_path    = np.empty(T);  Q_bF_path = np.empty(T)
    b_F_D_path   = np.empty(T);  b_D_F_path = np.empty(T)
    ic_bD_D_path = np.empty(T)
    ic_bF_F_path = np.empty(T)

    # terminal conditions (SS values as the t=T continuation)
    alpha_D_next = ss_bk_D["alpha_ss"]
    alpha_F_next = ss_bk_F["alpha_ss"]
    ic_spread_bD_ss = ss_bk_D["lambda_bD"] * ss_bk_D["mu_ss"] / ss_bk_D["Omega_ss"]
    ic_spread_bF_ss = ss_bk_F["lambda_bF"] * ss_bk_F["mu_ss"] / ss_bk_F["Omega_ss"]
    Q_bD_ss_val = db_D / (cal["r_dep_D_target"] + db_D + ic_spread_bD_ss)
    Q_bF_ss_val = db_F / (cal["r_dep_F_target"] + db_F + ic_spread_bF_ss)
    Q_bD_next   = Q_bD_ss_val
    Q_bF_next   = Q_bF_ss_val

    risk_mode = risk_D is not None
    tpi_mode  = tpi_D is not None
    if tpi_mode and not risk_mode:
        raise ValueError("tpi_D requires risk_D to also be provided (TPI pricing "
                         "nests inside the two-branch default-risk kernel).")
    if risk_mode:
        pi_path = np.asarray(risk_D["pi"])
        Om_d_D  = np.broadcast_to(risk_D["Omega_d_D"], T)
        Om_d_F  = np.broadcast_to(risk_D["Omega_d_F"], T)
        rk_d_D  = risk_D["rk_d_D"];  rk_d_F = risk_D["rk_d_F"]
        Q_bD_d  = risk_D["Q_bD_d"];  Q_bF_d = risk_D["Q_bF_d"]
        p_d     = risk_D["p_d"]
        surv_d  = float(np.asarray(risk_D.get("surv_d", rec_D)))   # priced-event survival

        if tpi_mode:
            pi_tpi_path = np.asarray(tpi_D["pi"])
            Om_tpi_D    = np.broadcast_to(tpi_D["Omega_tpi_D"], T)
            Om_tpi_F    = np.broadcast_to(tpi_D["Omega_tpi_F"], T)
            rk_tpi_d_D  = tpi_D["rk_tpi_d_D"];  rk_tpi_d_F = tpi_D["rk_tpi_d_F"]
            Q_bD_tpi_d  = tpi_D["Q_bD_tpi_d"];  Q_bF_tpi_d = tpi_D["Q_bF_tpi_d"]
            p_tpi_d     = tpi_D["p_tpi_d"]
            surv_tpi_d  = float(np.asarray(tpi_D.get("surv_tpi_d", 1.0)))
        else:
            # inert placeholders: pi_tpi = 1 forces w_tpi = 0 exactly below, so
            # these values can never contribute
            pi_tpi_path = np.ones(T)
            Om_tpi_D = np.zeros(T);  Om_tpi_F = np.zeros(T)
            rk_tpi_d_D = 0.0;        rk_tpi_d_F = 0.0
            Q_bD_tpi_d = 1.0;        Q_bF_tpi_d = 1.0
            p_tpi_d = 1.0
            surv_tpi_d = 1.0

    for t in range(T - 1, -1, -1):
        rk_D_next = rk_D[t + 1] if t + 1 < T else rk_D[T - 1]   # rk[T] ~ rk[T-1]
        rk_F_next = rk_F[t + 1] if t + 1 < T else rk_F[T - 1]
        p_next    = p_path[t + 1] if t + 1 < T else p_path[t]

        if not risk_mode:
            Omega_D = bi_D * (f_D + (1 - f_D) * alpha_D_next)
            mu_D    = Omega_D * (rk_D_next - rdep_D[t]) / lK_D
            if mu_D >= 1.0:
                raise RuntimeError(f"D-bank mu_D={mu_D:.4f} >= 1 at t={t}; IC infeasible.")
            alpha_D = Omega_D * (1 + rdep_D[t]) / (1 - mu_D)

            defp_D_next    = def_price_D[t + 1] if t + 1 < T else 0.0
            surv_D_price   = 1.0 - defp_D_next * (1.0 - rec_D)   # priced default prob
            ic_spread_bD_D = lbD_D * mu_D / Omega_D
            payoff_D_nd = db_D + (1 - db_D) * Q_bD_next
            Q_bD = surv_D_price * payoff_D_nd / (1 + rdep_D[t] + ic_spread_bD_D)

            Omega_F = bi_F * (f_F + (1 - f_F) * alpha_F_next)
            mu_F    = Omega_F * (rk_F_next - rdep_F[t]) / lK_F
            if mu_F >= 1.0:
                raise RuntimeError(f"F-bank mu_F={mu_F:.4f} >= 1 at t={t}; IC infeasible.")
            alpha_F = Omega_F * (1 + rdep_F[t]) / (1 - mu_F)

            ic_spread_bF_F = lbF_F * mu_F / Omega_F
            Q_bF = (db_F + (1 - db_F) * Q_bF_next) / (1 + rdep_F[t] + ic_spread_bF_F)

            rb_F_in_D = (((db_F + (1 - db_F) * Q_bF_next) / Q_bF)
                         * p_next / p_path[t] - 1)
            rb_D_in_F = ((surv_D_price * payoff_D_nd / Q_bD)
                         * p_path[t] / p_next - 1)
            ic_required_bF_D = lbF_D * mu_D / Omega_D
            ic_required_bD_F = lbD_F * mu_F / Omega_F
        else:
            # Bocola two-branch (or +TPI three-branch) step: bankers weight the
            # no-event continuation against a default branch (prob pi_def1) and,
            # in tpi_mode, a "backstop reneged" branch. pi_tpi1 = 1 forces
            # w_tpi = 0 in exact floating point, collapsing this to two branches.
            pi_def1 = pi_path[t + 1] if t + 1 < T else 0.0
            pi_tpi1 = pi_tpi_path[t + 1] if t + 1 < T else 1.0

            w_nd  = (1.0 - pi_def1) * pi_tpi1
            w_def = pi_def1
            w_tpi = (1.0 - pi_def1) * (1.0 - pi_tpi1)

            Omega_nd_D  = bi_D * (f_D + (1 - f_D) * alpha_D_next)
            Omega_til_D = w_nd * Omega_nd_D + w_def * Om_d_D[t] + w_tpi * Om_tpi_D[t]
            mu_D = (w_nd * Omega_nd_D * (rk_D_next - rdep_D[t])
                    + w_def * Om_d_D[t] * (rk_d_D - rdep_D[t])
                    + w_tpi * Om_tpi_D[t] * (rk_tpi_d_D - rdep_D[t])) / lK_D
            if mu_D >= 1.0:
                raise RuntimeError(f"D-bank mu_D={mu_D:.4f} >= 1 at t={t}; IC infeasible.")
            alpha_D = Omega_til_D * (1 + rdep_D[t]) / (1 - mu_D)

            payoff_D_nd  = db_D + (1 - db_D) * Q_bD_next
            payoff_D_d   = surv_d * (db_D + (1 - db_D) * Q_bD_d)          # haircut
            payoff_D_tpi = surv_tpi_d * (db_D + (1 - db_D) * Q_bD_tpi_d)  # repriced only
            ic_spread_bD_D = lbD_D * mu_D / Omega_til_D
            Q_bD_free = ((w_nd * Omega_nd_D * payoff_D_nd
                         + w_def * Om_d_D[t] * payoff_D_d
                         + w_tpi * Om_tpi_D[t] * payoff_D_tpi)
                        / (Omega_til_D * (1 + rdep_D[t]) + lbD_D * mu_D))


            Omega_nd_F  = bi_F * (f_F + (1 - f_F) * alpha_F_next)
            Omega_til_F = w_nd * Omega_nd_F + w_def * Om_d_F[t] + w_tpi * Om_tpi_F[t]
            mu_F = (w_nd * Omega_nd_F * (rk_F_next - rdep_F[t])
                    + w_def * Om_d_F[t] * (rk_d_F - rdep_F[t])
                    + w_tpi * Om_tpi_F[t] * (rk_tpi_d_F - rdep_F[t])) / lK_F
            if mu_F >= 1.0:
                raise RuntimeError(f"F-bank mu_F={mu_F:.4f} >= 1 at t={t}; IC infeasible.")
            alpha_F = Omega_til_F * (1 + rdep_F[t]) / (1 - mu_F)

            # F-bonds take no haircut in either branch, but the price jumps to the
            # branch level (safe-haven repricing)
            payoff_F_nd  = db_F + (1 - db_F) * Q_bF_next
            payoff_F_d   = db_F + (1 - db_F) * Q_bF_d
            payoff_F_tpi = db_F + (1 - db_F) * Q_bF_tpi_d
            ic_spread_bF_F = lbF_F * mu_F / Omega_til_F
            Q_bF = ((w_nd * Omega_nd_F * payoff_F_nd
                    + w_def * Om_d_F[t] * payoff_F_d
                    + w_tpi * Om_tpi_F[t] * payoff_F_tpi)
                   / (Omega_til_F * (1 + rdep_F[t]) + lbF_F * mu_F))

            # cross-border FOCs: certainty-equivalent returns under Omega-tilde,
            # branch legs valued at branch prices
            gross_F_in_D_nd  = payoff_F_nd / Q_bF * p_next / p_path[t]
            gross_F_in_D_d   = payoff_F_d / Q_bF * p_d / p_path[t]
            gross_F_in_D_tpi = payoff_F_tpi / Q_bF * p_tpi_d / p_path[t]
            rb_F_in_D = ((w_nd * Omega_nd_D * gross_F_in_D_nd
                          + w_def * Om_d_D[t] * gross_F_in_D_d
                          + w_tpi * Om_tpi_D[t] * gross_F_in_D_tpi) / Omega_til_D) - 1

            # D-bonds are valued at the ACTUAL traded price (post-floor when TPI is
            # mechanically active): the F-bank buys at the same CB-supported price
            gross_D_in_F_nd  = payoff_D_nd / Q_bD * p_path[t] / p_next
            gross_D_in_F_d   = payoff_D_d / Q_bD * p_path[t] / p_d
            gross_D_in_F_tpi = payoff_D_tpi / Q_bD * p_path[t] / p_tpi_d
            rb_D_in_F = ((w_nd * Omega_nd_F * gross_D_in_F_nd
                          + w_def * Om_d_F[t] * gross_D_in_F_d
                          + w_tpi * Om_tpi_F[t] * gross_D_in_F_tpi) / Omega_til_F) - 1

            ic_required_bF_D = lbF_D * mu_D / Omega_til_D
            ic_required_bD_F = lbD_F * mu_F / Omega_til_F
            Omega_D = Omega_til_D
            Omega_F = Omega_til_F

        # cross-border FOC holdings (end-of-period)
        b_F_D_t = (b_F_D_ss
                   + (rb_F_in_D - rdep_D[t] - exc_FD_ss - ic_required_bF_D)
                   / psi_bFD)
        b_D_F_t = (b_D_F_ss
                   + (rb_D_in_F - rdep_F[t] - exc_DF_ss - ic_required_bD_F)
                   / psi_bDF)

        alpha_D_path[t] = alpha_D;  mu_D_path[t] = mu_D;  Omega_D_path[t] = Omega_D
        alpha_F_path[t] = alpha_F;  mu_F_path[t] = mu_F;  Omega_F_path[t] = Omega_F
        Q_bD_path[t]    = Q_bD;     Q_bF_path[t] = Q_bF
        b_F_D_path[t]   = b_F_D_t;  b_D_F_path[t] = b_D_F_t
        ic_bD_D_path[t] = ic_spread_bD_D
        ic_bF_F_path[t] = ic_spread_bF_F

        alpha_D_next = alpha_D;  alpha_F_next = alpha_F
        Q_bD_next    = Q_bD;     Q_bF_next    = Q_bF

    return dict(
        alpha_D=alpha_D_path, mu_D=mu_D_path, Omega_D=Omega_D_path,
        alpha_F=alpha_F_path, mu_F=mu_F_path, Omega_F=Omega_F_path,
        Q_bD=Q_bD_path, Q_bF=Q_bF_path,
        b_F_D=b_F_D_path, b_D_F=b_D_F_path,
        ic_spread_bD_D=ic_bD_D_path, ic_spread_bF_F=ic_bF_F_path,
        Q_bD_ss_val=Q_bD_ss_val, Q_bF_ss_val=Q_bF_ss_val,
    )


def bank_forward(Kap_D, Kap_F, Q_D, Q_F, rk_D, rk_F, rdep_D, rdep_F, p_path,
                 b_D_D_path, b_F_F_path, bwd, cal, ss_bk_D, ss_bk_F,
                 def_real_D=None,
                 init_D=None, init_F=None,
                 Q_bD_lag0=None, Q_bF_lag0=None, p_lag0=None,
                 L_wc_D=None, L_wc_F=None):
    # FORWARD PASS: NET WORTH, DIVIDENDS, DEPOSIT SUPPLY FROM REALIZED RETURNS.
    # L_wc_* are the WORKING-CAPITAL LOANS the bank holds (Bocola SV.C): assets and
    # divertable base like any other class, earning lambda*mu/Omega over the deposit
    # rate. They default to the SS level so a constant-SS input path still reproduces
    # the SS exactly. (Legacy path-based block: the projection solver uses point_map.)
    T = len(Kap_D)

    if def_real_D is None:
        def_real_D = np.zeros(T)

    f_D   = cal["f_D"];           f_F   = cal["f_F"]
    lK_D  = cal["lambda_K_D"];    lK_F  = cal["lambda_K_F"]
    lbD_D = cal["lambda_bD_D"];   lbD_F = cal["lambda_bD_F"]
    lbF_D = cal["lambda_bF_D"];   lbF_F = cal["lambda_bF_F"]
    rec_D = cal["recovery_rate_D"]
    db_D  = cal["delta_b_D"];     db_F  = cal["delta_b_F"]

    L_wc_D = (np.full(T, ss_bk_D["L_wc_ss"]) if L_wc_D is None
              else np.asarray(L_wc_D, dtype=float))
    L_wc_F = (np.full(T, ss_bk_F["L_wc_ss"]) if L_wc_F is None
              else np.asarray(L_wc_F, dtype=float))

    Q_bD_path = bwd["Q_bD"];  Q_bF_path = bwd["Q_bF"]
    b_F_D_path = bwd["b_F_D"];  b_D_F_path = bwd["b_D_F"]
    alpha_D_path = bwd["alpha_D"];  alpha_F_path = bwd["alpha_F"]

    # realized bond returns on positions bought at t-1 (REALIZED survival only)
    Q_bD_l0 = bwd["Q_bD_ss_val"] if Q_bD_lag0 is None else Q_bD_lag0
    Q_bF_l0 = bwd["Q_bF_ss_val"] if Q_bF_lag0 is None else Q_bF_lag0
    Q_bD_lag = np.concatenate(([Q_bD_l0], Q_bD_path[:-1]))
    Q_bF_lag = np.concatenate(([Q_bF_l0], Q_bF_path[:-1]))

    surv_D_real = 1.0 - np.asarray(def_real_D) * (1.0 - rec_D)
    rb_D_path = (db_D * surv_D_real + (1 - db_D) * Q_bD_path * surv_D_real) / Q_bD_lag - 1
    rb_F_path = (db_F + (1 - db_F) * Q_bF_path) / Q_bF_lag - 1

    n_IC_D = np.empty(T);  n_ACCUM_D = np.empty(T)
    rn_D   = np.empty(T);  div_D     = np.empty(T)
    n_IC_F = np.empty(T);  n_ACCUM_F = np.empty(T)
    rn_F   = np.empty(T);  div_F     = np.empty(T)

    if init_D is None:
        init_D = dict(n_prev=ss_bk_D["n_ss"], kappa_prev=ss_bk_D["kappa_ss"],
                      phi_bdom_prev=ss_bk_D["phi_bdom_ss"],
                      phi_bfor_prev=ss_bk_D["phi_bfor_ss"],
                      rdep_prev=cal["r_dep_D_target"])
    n_D_prev        = init_D["n_prev"]
    phi_wc_D_prev   = init_D.get("phi_wc_prev", ss_bk_D["phi_wc_ss"])
    kappa_D_prev    = init_D["kappa_prev"]
    phi_bdom_D_prev = init_D["phi_bdom_prev"]
    phi_bfor_D_prev = init_D["phi_bfor_prev"]
    rdep_D_prev     = init_D["rdep_prev"]

    if init_F is None:
        init_F = dict(n_prev=ss_bk_F["n_ss"], kappa_prev=ss_bk_F["kappa_ss"],
                      phi_bdom_prev=ss_bk_F["phi_bdom_ss"],
                      phi_bfor_prev=ss_bk_F["phi_bfor_ss"],
                      rdep_prev=cal["r_dep_F_target"])
    n_F_prev        = init_F["n_prev"]
    phi_wc_F_prev   = init_F.get("phi_wc_prev", ss_bk_F["phi_wc_ss"])
    kappa_F_prev    = init_F["kappa_prev"]
    phi_bdom_F_prev = init_F["phi_bdom_prev"]
    phi_bfor_F_prev = init_F["phi_bfor_prev"]
    rdep_F_prev     = init_F["rdep_prev"]

    p_l0 = cal.get("p_ss", 1.0) if p_lag0 is None else p_lag0

    for t in range(T):
        p_t   = p_path[t]
        p_lag = p_path[t - 1] if t > 0 else p_l0

        # D-bank: rb_D on D-bonds, rb_F on F-bonds converted F -> D goods
        rb_D_t = rb_D_path[t]
        rb_F_t = (1.0 + rb_F_path[t]) * p_t / p_lag - 1

        rwc_D_t = lK_D * bwd["mu_D"][t] / bwd["Omega_D"][t]      # r_wc - rdep
        rn_D_t = (kappa_D_prev * (rk_D[t] - rdep_D_prev)
                  + phi_bdom_D_prev * (rb_D_t - rdep_D_prev)
                  + phi_bfor_D_prev * (rb_F_t - rdep_D_prev)
                  + phi_wc_D_prev * rwc_D_t
                  + rdep_D_prev)
        gross_D = (1 + rn_D_t) * n_D_prev
        total_assets_D = (Q_D[t] * Kap_D[t] + Q_bD_path[t] * b_D_D_path[t]
                          + p_t * Q_bF_path[t] * b_F_D_path[t] + L_wc_D[t])
        entrant_D = cal["omega_ent_D"] * total_assets_D
        n_ACCUM_D_t = (1 - f_D) * gross_D + entrant_D
        n_ACCUM_D[t] = n_ACCUM_D_t
        rn_D[t]  = rn_D_t
        div_D[t] = f_D * gross_D - entrant_D

        # n_IC is the net worth at which the IC binds exactly (F-bonds enter as
        # p x Q_bF); slack = alpha*(n - n_IC) >= 0
        n_IC_D[t] = (lK_D * Q_D[t] * Kap_D[t]
                     + lbD_D * Q_bD_path[t] * b_D_D_path[t]
                     + lbF_D * p_t * Q_bF_path[t] * b_F_D_path[t]
                     + lK_D * L_wc_D[t]) / alpha_D_path[t]
        # portfolio shares on ACTUAL net worth (equal to n_IC only when it binds)
        kappa_D_prev    = Q_D[t] * Kap_D[t] / n_ACCUM_D_t
        phi_bdom_D_prev = Q_bD_path[t] * b_D_D_path[t] / n_ACCUM_D_t
        phi_bfor_D_prev = p_t * Q_bF_path[t] * b_F_D_path[t] / n_ACCUM_D_t
        phi_wc_D_prev   = L_wc_D[t] / n_ACCUM_D_t
        n_D_prev        = n_ACCUM_D_t
        rdep_D_prev     = rdep_D[t]

        # F-bank: rb_F on F-bonds, rb_D on D-bonds converted D -> F goods
        rb_F_fg_t = rb_F_path[t]
        rb_D_fg_t = (1 + rb_D_path[t]) * p_lag / p_t - 1

        rwc_F_t = lK_F * bwd["mu_F"][t] / bwd["Omega_F"][t]      # r_wc - rdep
        rn_F_t = (kappa_F_prev * (rk_F[t] - rdep_F_prev)
                  + phi_bdom_F_prev * (rb_F_fg_t - rdep_F_prev)
                  + phi_bfor_F_prev * (rb_D_fg_t - rdep_F_prev)
                  + phi_wc_F_prev * rwc_F_t
                  + rdep_F_prev)
        gross_F = (1 + rn_F_t) * n_F_prev
        total_assets_F = (Q_F[t] * Kap_F[t]
                          + Q_bF_path[t] * b_F_F_path[t]
                          + Q_bD_path[t] * b_D_F_path[t] / p_t + L_wc_F[t])
        entrant_F = cal["omega_ent_F"] * total_assets_F
        n_ACCUM_F_t = (1 - f_F) * gross_F + entrant_F
        n_ACCUM_F[t] = n_ACCUM_F_t
        rn_F[t]  = rn_F_t
        div_F[t] = f_F * gross_F - entrant_F

        n_IC_F[t] = (lK_F * Q_F[t] * Kap_F[t]
                     + lbF_F * Q_bF_path[t] * b_F_F_path[t]
                     + lbD_F * Q_bD_path[t] * b_D_F_path[t] / p_t
                     + lK_F * L_wc_F[t]) / alpha_F_path[t]
        kappa_F_prev    = Q_F[t] * Kap_F[t] / n_ACCUM_F_t
        phi_bdom_F_prev = Q_bF_path[t] * b_F_F_path[t] / n_ACCUM_F_t
        phi_bfor_F_prev = Q_bD_path[t] * b_D_F_path[t] / (p_t * n_ACCUM_F_t)
        phi_wc_F_prev   = L_wc_F[t] / n_ACCUM_F_t
        n_F_prev        = n_ACCUM_F_t
        rdep_F_prev     = rdep_F[t]

    theta_D = ((Q_D * Kap_D + Q_bD_path * b_D_D_path + p_path * Q_bF_path * b_F_D_path
                + L_wc_D) / n_ACCUM_D)
    theta_F = ((Q_F * Kap_F + Q_bF_path * b_F_F_path + Q_bD_path * b_D_F_path / p_path
                + L_wc_F) / n_ACCUM_F)

    return dict(
        n_IC_D=n_IC_D, n_D=n_ACCUM_D, rn_D=rn_D, div_D=div_D,
        theta_D=theta_D, Dep_supply_D=(theta_D - 1) * n_ACCUM_D,
        n_IC_F=n_IC_F, n_F=n_ACCUM_F, rn_F=rn_F, div_F=div_F,
        theta_F=theta_F, Dep_supply_F=(theta_F - 1) * n_ACCUM_F,
        rb_D=rb_D_path, rb_F=rb_F_path,
        b_D_D=b_D_D_path, b_F_F=b_F_F_path,
    )
