"""Two-country Gertler-Karadi financial intermediary block (Bocola 2016 variant).

Each country has a bank that holds three assets:
  D-bank: domestic capital K_D, domestic D-bonds b_D_D, foreign F-bonds b_F_D
  F-bank: domestic capital K_F, domestic F-bonds b_F_F, foreign D-bonds b_D_F

Bond denomination convention (consistent throughout):
  D-bonds (Q_bD): D-good claims — priced by the D-bank's backward pass using rdep_D.
  F-bonds (Q_bF): F-good claims — priced by the F-bank's backward pass using rdep_F.
  Cross-border positions:
    D-bank's F-bond leg:   D-good value = p · Q_bF · b_F_D
    F-bank's D-bond leg:   F-good value = Q_bD · b_D_F / p

NOTE on p convention: p = price of F-goods in D-good units (D-goods per
F-good). An increase in p means F-goods are more expensive.
"""
import numpy as np
from scipy.optimize import brentq


# ─────────────────────────────────────────────────────────────────────────────
# Steady-state helpers
# ─────────────────────────────────────────────────────────────────────────────

def _alpha_ss_fixed_point(beta_inter, f, lambda_K, rk_ss, rdep_ss,
                           v_lo=1e-6, v_hi=1e6, n_scan=300):
    # Calculates the marginal value of networth for a banker V(n)=alpha_ss*n_t.  Omega = beta_inter*((1-f)+f*alpha_ss).  And IC multiplier
    
    def resid(a):
        Omega = beta_inter * ((1 - f) + f * a)
        mu    = Omega * (rk_ss - rdep_ss) / lambda_K
        if mu >= 1.0:
            return np.inf
        return Omega * (1 + rdep_ss) / (1 - mu) - a

    grid = np.geomspace(v_lo, v_hi, n_scan)
    vals = np.array([resid(v) for v in grid])
    fin  = np.isfinite(vals)
    sc   = np.where(np.diff(np.sign(vals[fin])) != 0)[0]
    if len(sc) == 0:
        raise RuntimeError(
            f"No sign change in alpha_ss fixed point (rk_ss={rk_ss:.6f}, "
            f"rdep_ss={rdep_ss:.6f}); check lambda_K / beta_inter calibration."
        )
    gf  = grid[fin]
    i   = sc[0]
    alpha_ss = brentq(resid, gf[i], gf[i + 1], xtol=1e-13, rtol=1e-13)
    Omega_ss = beta_inter * ((1 - f) + f * alpha_ss)
    mu_ss    = Omega_ss * (rk_ss - rdep_ss) / lambda_K
    return alpha_ss, mu_ss, Omega_ss


def calibrate_bank_targets(beta_inter, f, rdep, theta_target, spread_target,
                           v_lo=1e-6, v_hi=1e6, n_scan=300):
    """Single-λ Gertler-Karadi SS calibration.

    Solves the divertability λ (common to all asset classes — Bocola IC) and the
    entrant transfer ω_ent that make the deterministic steady state hit:
      • leverage       θ = (total assets)/n = theta_target   (binding IC ⇒ θ = α/λ)
      • credit spread  rk − rdep = spread_target             (equal across assets
                                                              under a single λ)

    Inverts `_alpha_ss_fixed_point`: with λ = α/θ folded in, μ = Ω·s·θ/α and
    α = Ω(1+rdep)/(1−μ) is a scalar fixed point in α; then λ = α/θ and, from the
    net-worth accumulation SS, ω_ent = [1−(1−f)(1+rdep)]/θ − (1−f)·s.

    Returns (lambda_single, omega_ent, alpha, mu, Omega).
    """
    s = spread_target

    def resid(a):
        Omega = beta_inter * ((1 - f) + f * a)
        mu    = Omega * s * theta_target / a
        if mu >= 1.0:
            return np.inf
        return Omega * (1 + rdep) / (1 - mu) - a

    grid = np.geomspace(v_lo, v_hi, n_scan)
    vals = np.array([resid(v) for v in grid])
    fin  = np.isfinite(vals)
    sc   = np.where(np.diff(np.sign(vals[fin])) != 0)[0]
    if len(sc) == 0:
        raise RuntimeError(
            f"No sign change in bank-target α fixed point (theta={theta_target}, "
            f"spread={s}); check beta_inter / f / targets."
        )
    gf = grid[fin]
    i  = sc[0]
    alpha = brentq(resid, gf[i], gf[i + 1], xtol=1e-13, rtol=1e-13)

    Omega = beta_inter * ((1 - f) + f * alpha)
    mu    = Omega * s * theta_target / alpha
    lambda_single = alpha / theta_target

    D_val     = 1.0 - (1 - f) * (1 + rdep)
    omega_ent = D_val / theta_target - (1 - f) * s
    if omega_ent <= 0.0:
        raise RuntimeError(
            f"Infeasible targets: omega_ent={omega_ent:.4e} ≤ 0 "
            f"(theta={theta_target}, spread={s}, f={f}, rdep={rdep})."
        )
    return lambda_single, omega_ent, alpha, mu, Omega


def steady_state_bank(cal, rk_ss, Kap_ss, Q_bD_ss, Q_bF_ss,
                      b_dom_ss, b_for_ss, p_ss, country="D"):
    
    #Steady-state bank block for one country.
    f           = cal[f"f_{country}"]
    rdep_ss     = cal[f"r_dep_{country}_target"]
    beta_inter  = cal[f"beta_inter_{country}"]
    lambda_K    = cal[f"lambda_K_{country}"]
    
    # Divertability of each bond type as seen by this bank
    lambda_bD   = cal[f"lambda_bD_{country}"]   # D-bond divertability for this bank (Equal)
    lambda_bF   = cal[f"lambda_bF_{country}"]   # F-bond divertability for this bank (Equal)
    omega_ent   = cal[f"omega_ent_{country}"]   # Entrant transfer fraction of total assets
    delta_b_D   = cal["delta_b_D"]
    delta_b_F   = cal["delta_b_F"]

    #Getting to see the multipliers on IC-consistent bond prices
    alpha_ss, mu_ss, Omega_ss = _alpha_ss_fixed_point(
        beta_inter, f, lambda_K, rk_ss, rdep_ss
    )
    #Getting the spreads on the bonds from the IC constraint
    IC_spread_dom = lambda_bD * mu_ss / Omega_ss  # IC spread on domestic bond
    IC_spread_for = lambda_bF * mu_ss / Omega_ss  # IC spread on foreign bond

    if country == "D":
        Q_bdom_ss = delta_b_D / (rdep_ss + delta_b_D + IC_spread_dom)  # D-bond market price
        Q_bfor_ss = delta_b_F / (rdep_ss + delta_b_F + IC_spread_for)  # F-bond market price
    else:

        Q_bdom_ss = delta_b_F / (rdep_ss + delta_b_F + IC_spread_dom)  # F-bond price
        Q_bfor_ss = delta_b_D / (rdep_ss + delta_b_D + IC_spread_for)  # D-bond price

    # Bond excess returns at SS
    rb_dom_ss = rdep_ss + IC_spread_dom
    rb_for_ss = rdep_ss + IC_spread_for

    # Net worth from IC-consistent bond prices:  n = (λ_K·K + λ_bD·Q_bD·b_D + λ_bF·p·Q_bF·b_F) 
    if country == "D":
        ic_numerator = (lambda_K * Kap_ss
                        + lambda_bD * Q_bdom_ss * b_dom_ss
                        + lambda_bF * p_ss * Q_bfor_ss * b_for_ss)   # F-bonds: p×Q_bF (F→D goods)
    else:
        ic_numerator = (lambda_K * Kap_ss
                        + lambda_bD * Q_bdom_ss * b_dom_ss            # F-bonds already F-goods
                        + lambda_bF * Q_bfor_ss * b_for_ss / p_ss)   # D-bonds ÷p → F-goods

    n_ss_IC = ic_numerator / alpha_ss

    # Net worth from forward accumulation
    D_val = 1.0 - (1 - f) * (1 + rdep_ss) #Internmediary discount factor
    if D_val <= 0:
        raise ValueError(f"[{country}] D={D_val} ≤ 0: no stationary net-worth rest point.")

    if country == "D":
        total_assets = Kap_ss + Q_bdom_ss * b_dom_ss + p_ss * Q_bfor_ss * b_for_ss
        n_ss_ACCUM = (
            ((1 - f) * (rk_ss     - rdep_ss) + omega_ent) * Kap_ss
            + ((1 - f) * IC_spread_dom + omega_ent) * Q_bdom_ss * b_dom_ss
            + ((1 - f) * IC_spread_for + omega_ent) * p_ss * Q_bfor_ss * b_for_ss
        ) / D_val
    else:
        Kap_val   = Kap_ss
        bdom_val  = Q_bdom_ss * b_dom_ss           
        bfor_val  = Q_bfor_ss * b_for_ss / p_ss 
        total_assets = Kap_val + bdom_val + bfor_val
        n_ss_ACCUM = (
            ((1 - f) * (rk_ss     - rdep_ss) + omega_ent) * Kap_val
            + ((1 - f) * IC_spread_dom + omega_ent) * bdom_val
            + ((1 - f) * IC_spread_for + omega_ent) * bfor_val
        ) / D_val

    n_ss = n_ss_ACCUM
    if n_ss <= 0:
        raise ValueError(f"[{country}] n_ss={n_ss:.4f} ≤ 0 at rk_ss={rk_ss:.6f}.")

    #Leverage ratios
    if country == "D":
        kappa_ss    = Kap_ss / n_ss
        phi_bdom_ss = Q_bdom_ss * b_dom_ss / n_ss
        phi_bfor_ss = p_ss * Q_bfor_ss * b_for_ss / n_ss         
    else:
        kappa_ss    = Kap_ss / n_ss
        phi_bdom_ss = Q_bdom_ss * b_dom_ss / n_ss                  
        phi_bfor_ss = Q_bfor_ss * b_for_ss / (p_ss * n_ss)        

    #Total leverge ratio
    theta_ss = kappa_ss + phi_bdom_ss + phi_bfor_ss

    #Returns on net worth and gross income
    rn_ss = (kappa_ss * (rk_ss - rdep_ss)
             + phi_bdom_ss * (rb_dom_ss - rdep_ss)
             + phi_bfor_ss * (rb_for_ss - rdep_ss)
             + rdep_ss)
    gross_income_ss = (1 + rn_ss) * n_ss
    entrant_ss      = omega_ent * total_assets
    div_ss          = f * gross_income_ss - entrant_ss
    Dep_supply_ss   = (theta_ss - 1) * n_ss

    return dict(
        alpha_ss=alpha_ss, mu_ss=mu_ss, Omega_ss=Omega_ss,
        n_ss=n_ss, n_ss_IC=n_ss_IC, n_ss_ACCUM=n_ss_ACCUM,
        kappa_ss=kappa_ss, phi_bdom_ss=phi_bdom_ss, phi_bfor_ss=phi_bfor_ss,
        theta_ss=theta_ss, rn_ss=rn_ss, div_ss=div_ss, entrant_ss=entrant_ss,
        Dep_supply_ss=Dep_supply_ss,
        rb_dom_ss=rb_dom_ss, rb_for_ss=rb_for_ss,
        Q_bdom_IC=Q_bdom_ss, Q_bfor_IC=Q_bfor_ss,
        IC_spread_dom=IC_spread_dom, IC_spread_for=IC_spread_for,
        lambda_K=lambda_K, lambda_bD=lambda_bD, lambda_bF=lambda_bF,
        Kap_ss=Kap_ss, total_assets_ss=total_assets,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Transition-path bank block — backward pass (prices, multipliers, FOC holdings)
# ─────────────────────────────────────────────────────────────────────────────

def bank_backward(rk_D, rk_F, rdep_D, rdep_F, p_path,
                  cal, ss_bk_D, ss_bk_F,
                  def_price_D=None, def_price_F=None, risk_D=None):
    """Backward pass for both banks: value-function slopes, bond prices, and
    cross-border FOC holdings.  Needs NO holdings or debt stocks — prices come
    from marginal conditions only, which is what allows the debt stock to be
    forward-integrated afterwards (see transition._inner_economy).

    def_price_D/F : (T,) PRICED default probability paths (expected haircut
                    enters Q recursions and expected-return FOCs).  None → 0.

    risk_D : None → risk-neutral pricing (expected-haircut surv form, current
             behaviour).  Otherwise a dict implementing the Bocola (2016)
             RISK CHANNEL via two-branch expectations over the D-default
             event (see risk_branch.py):
               pi        : (T,) prob. default occurs at t (replaces the
                           def_price_D surv factor in D pricing — do not
                           double count; def_price_D itself is UNUSED in
                           risk mode.  def_price_F still applies: F default
                           is assumed independent of the D-event and priced
                           risk-neutrally — its survival factor multiplies
                           both F-bond branch payoffs; there is no F default
                           branch)
               Omega_d_D : (T,) D-bank branch discount Λ^d·[(1−f)+f·α^d(0)]
               Omega_d_F : (T,) F-bank analogue
               rk_d_D, rk_d_F : scalars — branch h=0 capital returns
               Q_bD_d, Q_bF_d : scalars — branch h=0 bond prices
               p_d       : scalar — branch h=0 real exchange rate
             With Ω^d > Ω^nd multiplying the low default-branch payoffs, both
             bonds and capital carry an endogenous risk premium beyond the
             expected loss (precautionary deleveraging).  pi ≡ 0 or
             Ω^d = Ω^nd with branch prices = base prices reproduce the
             risk-neutral formulas exactly (nesting).

    Returns dict (all length T unless noted):
      alpha_D, mu_D, Omega_D (=Ω̃_{t+1} used at t), ic_spread_bD_D,
      alpha_F, mu_F, Omega_F, ic_spread_bF_F,
      Q_bD, Q_bF          : bond price paths
      b_F_D, b_D_F        : cross-border holdings from adjustment-cost FOCs
      Q_bD_ss_val, Q_bF_ss_val : scalars, IC-consistent SS prices (lag anchors)
    """
    T = len(rk_D)

    # checks if there's anything that changes the risk of default
    if def_price_D is None:
        def_price_D = np.zeros(T)
    if def_price_F is None:
        def_price_F = np.zeros(T)

    # unpack some calibration params
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
    rec_D      = cal["recovery_rate_D"]; rec_F = cal["recovery_rate_F"]

    # Initialize paths for backward recursion -> margin slopes, bond prices, and FOC holdings.
    alpha_D_path = np.empty(T);  mu_D_path = np.empty(T)
    alpha_F_path = np.empty(T);  mu_F_path = np.empty(T)
    Omega_D_path = np.empty(T);  Omega_F_path = np.empty(T)
    Q_bD_path    = np.empty(T);  Q_bF_path = np.empty(T)
    b_F_D_path   = np.empty(T);  b_D_F_path = np.empty(T)
    ic_bD_D_path = np.empty(T)   
    ic_bF_F_path = np.empty(T)   

    # Terminal period conditions 
    alpha_D_next = ss_bk_D["alpha_ss"]
    alpha_F_next = ss_bk_F["alpha_ss"]
    ic_spread_bD_ss = ss_bk_D["lambda_bD"] * ss_bk_D["mu_ss"] / ss_bk_D["Omega_ss"]
    ic_spread_bF_ss = ss_bk_F["lambda_bF"] * ss_bk_F["mu_ss"] / ss_bk_F["Omega_ss"]
    Q_bD_ss_val = cal["delta_b_D"] / (cal["r_dep_D_target"] + cal["delta_b_D"] + ic_spread_bD_ss)
    Q_bF_ss_val = cal["delta_b_F"] / (cal["r_dep_F_target"] + cal["delta_b_F"] + ic_spread_bF_ss)
    Q_bD_next   = Q_bD_ss_val
    Q_bF_next   = Q_bF_ss_val

    # Risk-channel inputs - in the case of no change in the risk of default it doesn't matter
    risk_mode = risk_D is not None
    if risk_mode:
        pi_path   = np.asarray(risk_D["pi"])
        Om_d_D    = np.broadcast_to(risk_D["Omega_d_D"], T)
        Om_d_F    = np.broadcast_to(risk_D["Omega_d_F"], T)
        rk_d_D    = risk_D["rk_d_D"];  rk_d_F = risk_D["rk_d_F"]
        Q_bD_d    = risk_D["Q_bD_d"];  Q_bF_d = risk_D["Q_bF_d"]
        p_d       = risk_D["p_d"]
        # Survival factor of the priced event (partial restructuring allowed;
        # must match the default branch's realized haircut)
        surv_d    = float(np.asarray(risk_D.get("surv_d", rec_D)))

    # Backward recursion over time
    for t in range(T - 1, -1, -1):
        # terminal conditions rk[T] is unknown; use rk[T-1] as the SS approximation.
        rk_D_next = rk_D[t + 1] if t + 1 < T else rk_D[T - 1]
        rk_F_next = rk_F[t + 1] if t + 1 < T else rk_F[T - 1]

        defp_D_next = def_price_D[t + 1] if t + 1 < T else 0.0
        defp_F_next = def_price_F[t + 1] if t + 1 < T else 0.0
        p_next = p_path[t + 1] if t + 1 < T else p_path[t]  # terminal: p constant

        if not risk_mode:
            # Risk-neutral backward step
            #continuation value of the bankers 
            Omega_D = bi_D * ((1 - f_D) + f_D * alpha_D_next)
            #the ic multiplier
            mu_D    = Omega_D * (rk_D_next - rdep_D[t]) / lK_D
            if mu_D >= 1.0:
                raise RuntimeError(f"D-bank mu_D={mu_D:.4f} ≥ 1 at t={t}; IC infeasible.")
            # the franchise value (marginal value of net worth)
            alpha_D = Omega_D * (1 + rdep_D[t]) / (1 - mu_D)

            # HM pricing: Q_bD = surv^e_{t+1}·(db + (1-db)·Q_next) / (1 + rdep + IC_spread)
            ic_spread_bD_D = lbD_D * mu_D / Omega_D
            # surv = 1 when no priced risk (TFP run); < 1 in the CK sunspot
            # experiment — this is exactly where the sunspot enters pricing.
            surv_D_price   = 1.0 - defp_D_next * (1.0 - rec_D)
            Q_bD = surv_D_price * (db_D + (1 - db_D) * Q_bD_next) / (1 + rdep_D[t] + ic_spread_bD_D)  # numerator: expected payoff; denominator: funding + IC cost

            # F-bank repetition 
            Omega_F = bi_F * ((1 - f_F) + f_F * alpha_F_next)
            mu_F    = Omega_F * (rk_F_next - rdep_F[t]) / lK_F
            if mu_F >= 1.0:
                raise RuntimeError(f"F-bank mu_F={mu_F:.4f} ≥ 1 at t={t}; IC infeasible.")
            alpha_F = Omega_F * (1 + rdep_F[t]) / (1 - mu_F)

            ic_spread_bF_F = lbF_F * mu_F / Omega_F
            surv_F_price   = 1.0 - defp_F_next * (1.0 - rec_F)
            Q_bF = surv_F_price * (db_F + (1 - db_F) * Q_bF_next) / (1 + rdep_F[t] + ic_spread_bF_F)

            # Cross-border FOCs: expected returns with priced survival and their pricing parameter
            rb_F_in_D = ((surv_F_price * (db_F + (1 - db_F) * Q_bF_next) / Q_bF)
                         * p_next / p_path[t] - 1)
            ic_required_bF_D = lbF_D * mu_D / Omega_D
            rb_D_in_F = ((surv_D_price * (db_D + (1 - db_D) * Q_bD_next) / Q_bD)
                         * p_path[t] / p_next - 1)
            ic_required_bD_F = lbD_F * mu_F / Omega_F
        else:
            # ── Bocola risk-channel step: two-branch expectations ──
            # With prob pi1 the D government defaults at t+1: bond stock is
            # haircut to rec_D, prices/returns jump to the default-branch
            # values, and the banker's discount is Ω^d (high marginal value).
            pi1 = pi_path[t + 1] if t + 1 < T else 0.0

            # D-bank
            Omega_nd_D  = bi_D * ((1 - f_D) + f_D * alpha_D_next)
            Omega_til_D = (1 - pi1) * Omega_nd_D + pi1 * Om_d_D[t]
            mu_D = ((1 - pi1) * Omega_nd_D * (rk_D_next - rdep_D[t])
                    + pi1 * Om_d_D[t] * (rk_d_D - rdep_D[t])) / lK_D
            if mu_D >= 1.0:
                raise RuntimeError(f"D-bank mu_D={mu_D:.4f} ≥ 1 at t={t}; IC infeasible.")
            alpha_D = Omega_til_D * (1 + rdep_D[t]) / (1 - mu_D)

            payoff_D_nd = db_D + (1 - db_D) * Q_bD_next
            payoff_D_d  = surv_d * (db_D + (1 - db_D) * Q_bD_d)   # event haircut on the whole claim
            ic_spread_bD_D = lbD_D * mu_D / Omega_til_D
            Q_bD = (((1 - pi1) * Omega_nd_D * payoff_D_nd + pi1 * Om_d_D[t] * payoff_D_d)
                    / (Omega_til_D * (1 + rdep_D[t]) + lbD_D * mu_D))

            # F-bank (same aggregate D-default event drives its expectations)
            Omega_nd_F  = bi_F * ((1 - f_F) + f_F * alpha_F_next)
            Omega_til_F = (1 - pi1) * Omega_nd_F + pi1 * Om_d_F[t]
            mu_F = ((1 - pi1) * Omega_nd_F * (rk_F_next - rdep_F[t])
                    + pi1 * Om_d_F[t] * (rk_d_F - rdep_F[t])) / lK_F
            if mu_F >= 1.0:
                raise RuntimeError(f"F-bank mu_F={mu_F:.4f} ≥ 1 at t={t}; IC infeasible.")
            alpha_F = Omega_til_F * (1 + rdep_F[t]) / (1 - mu_F)

            # F-bonds carry no haircut from the D-event but their price jumps
            # to Q_bF_d in it (safe-haven repricing → negative risk premium).
            # F's OWN default risk (independent of the D-event, priced
            # risk-neutrally — no F default branch) enters as the same
            # survival factor as in risk-neutral mode, on both branches.
            surv_F_price = 1.0 - defp_F_next * (1.0 - rec_F)
            payoff_F_nd = surv_F_price * (db_F + (1 - db_F) * Q_bF_next)
            payoff_F_d  = surv_F_price * (db_F + (1 - db_F) * Q_bF_d)
            ic_spread_bF_F = lbF_F * mu_F / Omega_til_F
            Q_bF = (((1 - pi1) * Omega_nd_F * payoff_F_nd + pi1 * Om_d_F[t] * payoff_F_d)
                    / (Omega_til_F * (1 + rdep_F[t]) + lbF_F * mu_F))

            # Cross-border FOCs: certainty-equivalent returns under Ω̃
            # (E[Ω'(1+rb')]/Ω̃ − 1), branch legs converted with branch prices.
            gross_F_in_D_nd = payoff_F_nd / Q_bF * p_next / p_path[t]
            gross_F_in_D_d  = payoff_F_d / Q_bF * p_d / p_path[t]
            rb_F_in_D = (((1 - pi1) * Omega_nd_D * gross_F_in_D_nd
                          + pi1 * Om_d_D[t] * gross_F_in_D_d) / Omega_til_D) - 1
            ic_required_bF_D = lbF_D * mu_D / Omega_til_D

            gross_D_in_F_nd = payoff_D_nd / Q_bD * p_path[t] / p_next
            gross_D_in_F_d  = payoff_D_d / Q_bD * p_path[t] / p_d
            rb_D_in_F = (((1 - pi1) * Omega_nd_F * gross_D_in_F_nd
                          + pi1 * Om_d_F[t] * gross_D_in_F_d) / Omega_til_F) - 1
            ic_required_bD_F = lbD_F * mu_F / Omega_til_F

            Omega_D = Omega_til_D   # stored below (Ω̃ used at t)
            Omega_F = Omega_til_F

        # Cross-border FOC holdings (end-of-period, after realized returns)
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


# ─────────────────────────────────────────────────────────────────────────────
# Transition-path bank block — forward pass (net worth, dividends, deposits)
# ─────────────────────────────────────────────────────────────────────────────

def bank_forward(Kap_D, Kap_F, Q_D, Q_F, rk_D, rk_F, rdep_D, rdep_F, p_path,
                 b_D_D_path, b_F_F_path, bwd, cal, ss_bk_D, ss_bk_F,
                 def_real_D=None, def_real_F=None,
                 init_D=None, init_F=None,
                 Q_bD_lag0=None, Q_bF_lag0=None, p_lag0=None):
    # GIVEN THE PRICES AND CROSS BORDER HOLDINGS FROM THE BACKWARD PASS, LET US CALCULATE THE FORWARD PATHS OF NET WORTH AND DIVIDENDS
    T = len(Kap_D)

    if def_real_D is None:
        def_real_D = np.zeros(T)
    if def_real_F is None:
        def_real_F = np.zeros(T)

    # unpack calibration
    f_D   = cal["f_D"];           f_F   = cal["f_F"]
    lK_D  = cal["lambda_K_D"];    lK_F  = cal["lambda_K_F"]
    lbD_D = cal["lambda_bD_D"];   lbD_F = cal["lambda_bD_F"]
    lbF_D = cal["lambda_bF_D"];   lbF_F = cal["lambda_bF_F"]
    rec_D = cal["recovery_rate_D"]; rec_F = cal["recovery_rate_F"]
    db_D  = cal["delta_b_D"];     db_F  = cal["delta_b_F"]

    # Entrant equity: "proportional" (GK, entrant = ω·assets_t) leaves n with
    # an exact unit root conditional on returns; "anchored" fixes the injection
    # at the SS level (optionally scaled by franchise value α^phi_entry), which
    # lowers the conditional eigenvalue to 1 − ω_ent·θ_ss (see calibration.py).
    entrant_mode = cal.get("entrant_mode", "proportional")
    phi_entry    = cal.get("phi_entry", 0.0)
    assets_ss_D  = ss_bk_D["total_assets_ss"];  alpha_ss_D = ss_bk_D["alpha_ss"]
    assets_ss_F  = ss_bk_F["total_assets_ss"];  alpha_ss_F = ss_bk_F["alpha_ss"]

    # unpack from the backward path 
    Q_bD_path = bwd["Q_bD"];  Q_bF_path = bwd["Q_bF"]
    b_F_D_path = bwd["b_F_D"];  b_D_F_path = bwd["b_D_F"]
    alpha_D_path = bwd["alpha_D"];  alpha_F_path = bwd["alpha_F"]

    # Realized returns at t on bonds bought at t-1 — REALIZED survival only.
    Q_bD_l0 = bwd["Q_bD_ss_val"] if Q_bD_lag0 is None else Q_bD_lag0
    Q_bF_l0 = bwd["Q_bF_ss_val"] if Q_bF_lag0 is None else Q_bF_lag0
    Q_bD_lag = np.concatenate(([Q_bD_l0], Q_bD_path[:-1]))
    Q_bF_lag = np.concatenate(([Q_bF_l0], Q_bF_path[:-1]))

    surv_D_real = 1.0 - np.asarray(def_real_D) * (1.0 - rec_D)
    surv_F_real = 1.0 - np.asarray(def_real_F) * (1.0 - rec_F)
    rb_D_path = (db_D * surv_D_real + (1 - db_D) * Q_bD_path * surv_D_real) / Q_bD_lag - 1
    rb_F_path = (db_F * surv_F_real + (1 - db_F) * Q_bF_path * surv_F_real) / Q_bF_lag - 1


    # initialize forward paths
    n_IC_D = np.empty(T);  n_ACCUM_D = np.empty(T)
    rn_D   = np.empty(T);  div_D     = np.empty(T)
    n_IC_F = np.empty(T);  n_ACCUM_F = np.empty(T)
    rn_F   = np.empty(T);  div_F     = np.empty(T)

    # D-bank forward state (SS by default
    if init_D is None:
        init_D = dict(n_prev=ss_bk_D["n_ss"], kappa_prev=ss_bk_D["kappa_ss"],
                      phi_bdom_prev=ss_bk_D["phi_bdom_ss"],
                      phi_bfor_prev=ss_bk_D["phi_bfor_ss"],
                      rdep_prev=cal["r_dep_D_target"])
    n_D_prev        = init_D["n_prev"]
    kappa_D_prev    = init_D["kappa_prev"]
    phi_bdom_D_prev = init_D["phi_bdom_prev"]
    phi_bfor_D_prev = init_D["phi_bfor_prev"]
    rdep_D_prev     = init_D["rdep_prev"]

    # F-bank forward state
    if init_F is None:
        init_F = dict(n_prev=ss_bk_F["n_ss"], kappa_prev=ss_bk_F["kappa_ss"],
                      phi_bdom_prev=ss_bk_F["phi_bdom_ss"],
                      phi_bfor_prev=ss_bk_F["phi_bfor_ss"],
                      rdep_prev=cal["r_dep_F_target"])
    n_F_prev        = init_F["n_prev"]
    kappa_F_prev    = init_F["kappa_prev"]
    phi_bdom_F_prev = init_F["phi_bdom_prev"]
    phi_bfor_F_prev = init_F["phi_bfor_prev"]
    rdep_F_prev     = init_F["rdep_prev"]

    p_l0 = cal.get("p_ss", 1.0) if p_lag0 is None else p_lag0

    for t in range(T):
        #price conversion
        p_t   = p_path[t]
        p_lag = p_path[t - 1] if t > 0 else p_l0

        # D-bank earns rb_D on D-bonds (D-goods) and rb_F on F-bonds (convert F→D goods).
        rb_D_t = rb_D_path[t]
        rb_F_t = (1.0 + rb_F_path[t]) * p_t / p_lag - 1   # F-goods → D-goods

        #return on net worth 
        rn_D_t = (kappa_D_prev * (rk_D[t] - rdep_D_prev)
                  + phi_bdom_D_prev * (rb_D_t - rdep_D_prev)
                  + phi_bfor_D_prev * (rb_F_t - rdep_D_prev)
                  + rdep_D_prev)
        gross_D = (1 + rn_D_t) * n_D_prev
        total_assets_D = (Q_D[t] * Kap_D[t] + Q_bD_path[t] * b_D_D_path[t]
                          + p_t * Q_bF_path[t] * b_F_D_path[t])
        if entrant_mode == "anchored":
            entrant_D = (cal["omega_ent_D"] * assets_ss_D
                         * (alpha_D_path[t] / alpha_ss_D) ** phi_entry)
        else:
            entrant_D = cal["omega_ent_D"] * total_assets_D
        #assets that remain in the bank after dividend payout and new entrants
        n_ACCUM_D_t = (1 - f_D) * gross_D + entrant_D
        #net dividends paid to shareholders (after new entrants)
        div_D_t     = f_D * gross_D - entrant_D
        n_ACCUM_D[t] = n_ACCUM_D_t
        rn_D[t]  = rn_D_t
        div_D[t] = div_D_t

        # D-bank IC in D-goods (F-bonds enter as p×Q_bF)
        n_IC_D_t = (lK_D * Q_D[t] * Kap_D[t]
                    + lbD_D * Q_bD_path[t] * b_D_D_path[t]
                    + lbF_D * p_t * Q_bF_path[t] * b_F_D_path[t]) / alpha_D_path[t]
        n_IC_D[t] = n_IC_D_t
        #calculate the leverage ratios for the next period
        kappa_D_prev    = Q_D[t] * Kap_D[t] / n_IC_D_t
        phi_bdom_D_prev = Q_bD_path[t] * b_D_D_path[t] / n_IC_D_t
        phi_bfor_D_prev = p_t * Q_bF_path[t] * b_F_D_path[t] / n_IC_D_t
        n_D_prev        = n_ACCUM_D_t
        rdep_D_prev     = rdep_D[t]

        # ── F-bank ─────────────────────────────────────────────────────────
        rb_F_fg_t = rb_F_path[t]                            # F-bonds: F-goods ✓
        rb_D_fg_t = (1 + rb_D_path[t]) * p_lag / p_t - 1    # D-bonds: D-goods → F-goods

        rn_F_t = (kappa_F_prev * (rk_F[t] - rdep_F_prev)
                  + phi_bdom_F_prev * (rb_F_fg_t - rdep_F_prev)
                  + phi_bfor_F_prev * (rb_D_fg_t - rdep_F_prev)
                  + rdep_F_prev)
        gross_F = (1 + rn_F_t) * n_F_prev
        total_assets_F = (Q_F[t] * Kap_F[t]
                          + Q_bF_path[t] * b_F_F_path[t]
                          + Q_bD_path[t] * b_D_F_path[t] / p_t)
        if entrant_mode == "anchored":
            entrant_F = (cal["omega_ent_F"] * assets_ss_F
                         * (alpha_F_path[t] / alpha_ss_F) ** phi_entry)
        else:
            entrant_F = cal["omega_ent_F"] * total_assets_F
        n_ACCUM_F_t = (1 - f_F) * gross_F + entrant_F
        div_F_t     = f_F * gross_F - entrant_F
        n_ACCUM_F[t] = n_ACCUM_F_t
        rn_F[t]  = rn_F_t
        div_F[t] = div_F_t

        # F-bank IC in F-goods (F-bonds are F-good claims; D-bonds ÷p)
        n_IC_F_t = (lK_F * Q_F[t] * Kap_F[t]
                    + lbF_F * Q_bF_path[t] * b_F_F_path[t]
                    + lbD_F * Q_bD_path[t] * b_D_F_path[t] / p_t) / alpha_F_path[t]
        n_IC_F[t] = n_IC_F_t

        kappa_F_prev    = Q_F[t] * Kap_F[t] / n_IC_F_t
        phi_bdom_F_prev = Q_bF_path[t] * b_F_F_path[t] / n_IC_F_t
        phi_bfor_F_prev = Q_bD_path[t] * b_D_F_path[t] / (p_t * n_IC_F_t)
        n_F_prev        = n_ACCUM_F_t
        rdep_F_prev     = rdep_F[t]

    theta_D = ((Q_D * Kap_D + Q_bD_path * b_D_D_path + p_path * Q_bF_path * b_F_D_path)
               / n_ACCUM_D)
    theta_F = ((Q_F * Kap_F + Q_bF_path * b_F_F_path + Q_bD_path * b_D_F_path / p_path)
               / n_ACCUM_F)

    Dep_supply_D = (theta_D - 1) * n_ACCUM_D
    Dep_supply_F = (theta_F - 1) * n_ACCUM_F

    return dict(
        n_IC_D=n_IC_D, n_D=n_ACCUM_D, rn_D=rn_D, div_D=div_D,
        theta_D=theta_D, Dep_supply_D=Dep_supply_D,
        n_IC_F=n_IC_F, n_F=n_ACCUM_F, rn_F=rn_F, div_F=div_F,
        theta_F=theta_F, Dep_supply_F=Dep_supply_F,
        rb_D=rb_D_path, rb_F=rb_F_path,
        b_D_D=b_D_D_path, b_F_F=b_F_F_path,
    )
