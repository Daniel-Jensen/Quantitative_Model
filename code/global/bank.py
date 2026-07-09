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

Multi-asset incentive constraint (IC), Bocola (2016) eq. (3):
  A banker can divert a fraction lambda of TOTAL assets, so with the linear
  value function V_t(n) = alpha_t · n the constraint is
    D-bank (D-goods):
      lambda_K · Q_D·K_D + lambda_bD · Q_bD·b_D_D + lambda_bF · p·Q_bF·b_F_D ≤ alpha_D · n_D
    F-bank (F-goods):
      lambda_K · Q_F·K_F + lambda_bF · Q_bF·b_F_F + lambda_bD · Q_bD·b_D_F/p ≤ alpha_F · n_F
  Following Bocola the baseline calibration sets a SINGLE lambda per bank
  (lambda_K = lambda_bD = lambda_bF), so the constraint is on total assets and
  banks cannot relax it by substituting between bonds and capital.  The code
  keeps the three slots separate for robustness exercises.

When the IC binds (imposed throughout, standard in perfect foresight), the
closed-form backward pass is (Bocola 2016 eq. (1)-(2), GK 2011):

  Omega_{t+1}   = beta_inter · [(1−f) + f · alpha_{t+1}]
  mu_t          = Omega_{t+1} · (rk_{t+1} − rdep_t) / lambda_K    [capital FOC]
  alpha_t       = Omega_{t+1} · (1 + rdep_t) / (1 − mu_t)          [Bellman]
  Q_bX_t        = surv^e_{t+1} · (delta_bX + (1−delta_bX)·Q_bX_{t+1})
                    / (1 + rdep_t + lambda_bX·mu_t/Omega_{t+1})

PRICED vs REALIZED default (Bocola 2016 experiment design):
  surv^e_{t+1} = 1 − def_price_{t+1}·(1−recovery)   [expected — enters PRICES
                  and all expected-return FOCs; the Bocola pass-through shock]
  surv^r_t     = 1 − def_real_t·(1−recovery)        [realized — enters realized
                  returns rb and the government's coupon/stock flows]
  The baseline Cole-Kehoe experiment sets def_price = sunspot (in the crisis
  zone) and def_real = 0: pure news of future default lowers Q on impact,
  inflicting a mark-to-market loss on legacy bond holders — net worth falls
  although no default ever happens (Bocola's "pass-through of sovereign risk").
  Ex post, banks earn above-required bond returns while beliefs persist
  (bought cheap, repaid in full), so net worth recovers as the sunspot decays.

Cross-border positions from portfolio adjustment-cost FOC (expected returns
use def_price):
  b_F_D_t = b_F_D_ss + [E_t rb_F_in_D_{t+1} − rdep_D_t
                         − excess_return_F_D_ss − lambda_bF·mu_D_t/Omega_D_{t+1}] / psi_bF_D
  (and symmetrically for F-bank's D-bond holding b_D_F_t)

Bond market clearing (in transition.py): banks jointly hold the END-of-period
outstanding stock from the government's budget identity:
  b_D_D_t + b_D_F_t = b_gov_D_eop_t
  b_F_D_t + b_F_F_t = b_gov_F_eop_t

Net worth has two characterisations that must agree (outer residual):
  n_IC    = IC-binding allocation size (desired by bank given alpha)
  n_ACCUM = forward accumulation (true state, carried from last period)
Their difference (n_IC − n_ACCUM) / n_ss is the capital-market residual
fed to the outer Newton solver to pin Kap_path for each country.
"""
import numpy as np
from scipy.optimize import brentq


# ─────────────────────────────────────────────────────────────────────────────
# Steady-state helpers
# ─────────────────────────────────────────────────────────────────────────────

def _alpha_ss_fixed_point(beta_inter, f, lambda_K, rk_ss, rdep_ss,
                           v_lo=1e-6, v_hi=1e6, n_scan=300):
    """Solve the self-referential fixed point for alpha_ss.

    At SS: alpha = Omega(alpha) · (1+rdep) / (1−mu(alpha))
    where Omega = beta_inter · [(1−f) + f·alpha]
          mu    = Omega · (rk − rdep) / lambda_K

    Uses scan-then-brentq (identical pattern to old bank.py).
    """
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


def steady_state_bank(cal, rk_ss, Kap_ss, Q_bD_ss, Q_bF_ss,
                      b_dom_ss, b_for_ss, p_ss, country="D"):
    """Steady-state bank block for one country.

    Arguments
    ---------
    rk_ss   : steady-state capital return (from outer solve)
    Kap_ss  : capital stock (from Cobb-Douglas demand)
    Q_bD_ss : D-bond price at SS (from government.py)
    Q_bF_ss : F-bond price at SS
    b_dom_ss: domestic-bond holding (D-good units; = B_gov_D - b_D_F_ss or B_gov_F - b_F_D_ss)
    b_for_ss: foreign-bond holding (D-good units; e.g. b_F_D_ss for D-bank)
    p_ss    : real exchange rate at SS (= 1 at symmetric SS)
    country : "D" or "F"

    Returns
    -------
    dict with keys: alpha_ss, mu_ss, Omega_ss, n_ss, n_ss_IC, n_ss_ACCUM,
                    kappa_ss, phi_bdom_ss, phi_bfor_ss, theta_ss,
                    rn_ss, div_ss, entrant_ss, Dep_supply_ss,
                    rb_dom_ss, rb_for_ss,
                    lambda_K, lambda_bD, lambda_bF, Kap_ss
    """
    f           = cal[f"f_{country}"]
    rdep_ss     = cal[f"r_dep_{country}_target"]
    beta_inter  = cal[f"beta_inter_{country}"]
    lambda_K    = cal[f"lambda_K_{country}"]
    # Divertability of each bond type as seen by this bank
    lambda_bD   = cal[f"lambda_bD_{country}"]   # D-bond divertability for this bank
    lambda_bF   = cal[f"lambda_bF_{country}"]   # F-bond divertability for this bank
    omega_ent   = cal[f"omega_ent_{country}"]
    delta_b_D   = cal["delta_b_D"]
    delta_b_F   = cal["delta_b_F"]

    alpha_ss, mu_ss, Omega_ss = _alpha_ss_fixed_point(
        beta_inter, f, lambda_K, rk_ss, rdep_ss
    )

    # IC-consistent SS bond prices: GK excess-return on each bond = lambda_b*mu/Omega.
    # Bond FOC at SS gives the MARKET price Q_bX = delta_bX/(rdep+delta_bX+IC_spread).
    # This is the fixed point of the backward pricing recurrence used in bank_backward.
    IC_spread_dom = lambda_bD * mu_ss / Omega_ss  # IC spread on domestic bond
    IC_spread_for = lambda_bF * mu_ss / Omega_ss  # IC spread on foreign bond

    if country == "D":
        Q_bdom_ss = delta_b_D / (rdep_ss + delta_b_D + IC_spread_dom)  # D-bond market price
        Q_bfor_ss = delta_b_F / (rdep_ss + delta_b_F + IC_spread_for)  # F-bond market price
    else:
        # F-bank: domestic = F-bonds (divertability lambda_bD is lambda_bF_F in cal)
        #         foreign  = D-bonds (divertability lambda_bF is lambda_bD_F in cal)
        Q_bdom_ss = delta_b_F / (rdep_ss + delta_b_F + IC_spread_dom)  # F-bond price
        Q_bfor_ss = delta_b_D / (rdep_ss + delta_b_D + IC_spread_for)  # D-bond price

    # Bond excess returns at SS (= IC spread above deposit rate)
    rb_dom_ss = rdep_ss + IC_spread_dom
    rb_for_ss = rdep_ss + IC_spread_for

    # n from IC binding: n_IC = (lambda_K·Q_K·K + lambda_bD·Q_bD·b_dom + lambda_bF·Q_bF·b_for) / alpha
    # For F-bank, D-bonds are D-good claims; divide by p to get F-good IC value
    if country == "D":
        ic_numerator = (lambda_K * Kap_ss
                        + lambda_bD * Q_bdom_ss * b_dom_ss
                        + lambda_bF * p_ss * Q_bfor_ss * b_for_ss)   # F-bonds: p×Q_bF (F→D goods)
    else:
        ic_numerator = (lambda_K * Kap_ss
                        + lambda_bD * Q_bdom_ss * b_dom_ss            # F-bonds already F-goods
                        + lambda_bF * Q_bfor_ss * b_for_ss / p_ss)   # D-bonds ÷p → F-goods

    n_ss_IC = ic_numerator / alpha_ss

    # n from forward accumulation: D·n = excess returns + entrant transfer
    # D = 1 − (1−f)·(1+rdep)
    D_val = 1.0 - (1 - f) * (1 + rdep_ss)
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
        # F-bank: assets in F-goods. F-bonds (Q_bF) are F-good claims → no conversion.
        # D-bonds (Q_bD) are D-good claims → divide by p_ss to get F-goods.
        Kap_val   = Kap_ss
        bdom_val  = Q_bdom_ss * b_dom_ss            # F-bonds: Q_bF × b_F_F already F-goods
        bfor_val  = Q_bfor_ss * b_for_ss / p_ss     # D-bonds: Q_bD × b_D_F ÷ p → F-goods
        total_assets = Kap_val + bdom_val + bfor_val
        n_ss_ACCUM = (
            ((1 - f) * (rk_ss     - rdep_ss) + omega_ent) * Kap_val
            + ((1 - f) * IC_spread_dom + omega_ent) * bdom_val
            + ((1 - f) * IC_spread_for + omega_ent) * bfor_val
        ) / D_val

    n_ss = n_ss_ACCUM
    if n_ss <= 0:
        raise ValueError(f"[{country}] n_ss={n_ss:.4f} ≤ 0 at rk_ss={rk_ss:.6f}.")

    if country == "D":
        kappa_ss    = Kap_ss / n_ss
        phi_bdom_ss = Q_bdom_ss * b_dom_ss / n_ss
        phi_bfor_ss = p_ss * Q_bfor_ss * b_for_ss / n_ss          # F-bonds valued in D-goods
    else:
        kappa_ss    = Kap_ss / n_ss
        phi_bdom_ss = Q_bdom_ss * b_dom_ss / n_ss                  # F-bonds already F-goods
        phi_bfor_ss = Q_bfor_ss * b_for_ss / (p_ss * n_ss)        # D-bonds ÷p

    theta_ss = kappa_ss + phi_bdom_ss + phi_bfor_ss

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
        # IC-consistent bond prices (used by steady_state.py to update Tax_ss)
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

    if def_price_D is None:
        def_price_D = np.zeros(T)
    if def_price_F is None:
        def_price_F = np.zeros(T)

    f_D        = cal["f_D"];          f_F        = cal["f_F"]
    bi_D       = cal["beta_inter_D"]; bi_F       = cal["beta_inter_F"]
    lK_D       = cal["lambda_K_D"];   lK_F       = cal["lambda_K_F"]
    lbD_D      = cal["lambda_bD_D"];  lbD_F      = cal["lambda_bD_F"]
    lbF_D      = cal["lambda_bF_D"];  lbF_F      = cal["lambda_bF_F"]
    db_D       = cal["delta_b_D"];    db_F       = cal["delta_b_F"]
    psi_bFD    = cal["psi_bF_D"]      # D-bank adj cost for F-bond deviation
    psi_bDF    = cal["psi_bD_F"]      # F-bank adj cost for D-bond deviation
    b_F_D_ss   = cal["b_F_D_ss"]
    b_D_F_ss   = cal["b_D_F_ss"]
    exc_FD_ss  = cal["excess_return_F_D_ss"]
    exc_DF_ss  = cal["excess_return_D_F_ss"]
    rec_D      = cal["recovery_rate_D"]; rec_F = cal["recovery_rate_F"]

    alpha_D_path = np.empty(T);  mu_D_path = np.empty(T)
    alpha_F_path = np.empty(T);  mu_F_path = np.empty(T)
    Omega_D_path = np.empty(T);  Omega_F_path = np.empty(T)
    Q_bD_path    = np.empty(T);  Q_bF_path = np.empty(T)
    b_F_D_path   = np.empty(T);  b_D_F_path = np.empty(T)
    ic_bD_D_path = np.empty(T)   # lambda_bD_D·mu_D/Omega_D (liquidity premium on D-bonds)
    ic_bF_F_path = np.empty(T)   # lambda_bF_F·mu_F/Omega_F (liquidity premium on F-bonds)

    alpha_D_next = ss_bk_D["alpha_ss"]
    alpha_F_next = ss_bk_F["alpha_ss"]

    # IC-consistent SS bond prices (terminal condition for backward pass).
    ic_spread_bD_ss = ss_bk_D["lambda_bD"] * ss_bk_D["mu_ss"] / ss_bk_D["Omega_ss"]
    ic_spread_bF_ss = ss_bk_F["lambda_bF"] * ss_bk_F["mu_ss"] / ss_bk_F["Omega_ss"]
    Q_bD_ss_val = cal["delta_b_D"] / (cal["r_dep_D_target"] + cal["delta_b_D"] + ic_spread_bD_ss)
    Q_bF_ss_val = cal["delta_b_F"] / (cal["r_dep_F_target"] + cal["delta_b_F"] + ic_spread_bF_ss)
    Q_bD_next   = Q_bD_ss_val
    Q_bF_next   = Q_bF_ss_val

    # Risk-channel inputs (two-branch expectations over the D-default event)
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

    for t in range(T - 1, -1, -1):
        # At the terminal period, rk[T] is unknown; use rk[T-1] as the SS approximation.
        rk_D_next = rk_D[t + 1] if t + 1 < T else rk_D[T - 1]
        rk_F_next = rk_F[t + 1] if t + 1 < T else rk_F[T - 1]
        # Next-period PRICED default probability (expected haircut in prices)
        defp_D_next = def_price_D[t + 1] if t + 1 < T else 0.0
        defp_F_next = def_price_F[t + 1] if t + 1 < T else 0.0
        p_next = p_path[t + 1] if t + 1 < T else p_path[t]  # terminal: p constant

        if not risk_mode:
            # ── Risk-neutral backward step (expected-haircut surv form) ──
            # D-bank
            Omega_D = bi_D * ((1 - f_D) + f_D * alpha_D_next)
            mu_D    = Omega_D * (rk_D_next - rdep_D[t]) / lK_D
            if mu_D >= 1.0:
                raise RuntimeError(f"D-bank mu_D={mu_D:.4f} ≥ 1 at t={t}; IC infeasible.")
            alpha_D = Omega_D * (1 + rdep_D[t]) / (1 - mu_D)

            # HM pricing: Q_bD = surv^e_{t+1}·(db + (1-db)·Q_next) / (1 + rdep + IC_spread)
            ic_spread_bD_D = lbD_D * mu_D / Omega_D
            surv_D_price   = 1.0 - defp_D_next * (1.0 - rec_D)
            Q_bD = surv_D_price * (db_D + (1 - db_D) * Q_bD_next) / (1 + rdep_D[t] + ic_spread_bD_D)

            # F-bank
            Omega_F = bi_F * ((1 - f_F) + f_F * alpha_F_next)
            mu_F    = Omega_F * (rk_F_next - rdep_F[t]) / lK_F
            if mu_F >= 1.0:
                raise RuntimeError(f"F-bank mu_F={mu_F:.4f} ≥ 1 at t={t}; IC infeasible.")
            alpha_F = Omega_F * (1 + rdep_F[t]) / (1 - mu_F)

            ic_spread_bF_F = lbF_F * mu_F / Omega_F
            surv_F_price   = 1.0 - defp_F_next * (1.0 - rec_F)
            Q_bF = surv_F_price * (db_F + (1 - db_F) * Q_bF_next) / (1 + rdep_F[t] + ic_spread_bF_F)

            # Cross-border FOCs: expected returns with priced survival
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
    """Forward pass for both banks: realized returns, net-worth accumulation,
    IC-implied net worth, dividends and deposit supply.

    b_D_D_path, b_F_F_path : domestic bond holdings from market clearing
                             (end-of-period government stock minus the foreign
                             bank's FOC holding — computed in transition.py)
    bwd                    : output dict of bank_backward()
    def_real_D/F           : (T,) REALIZED default (haircut) paths — enter
                             realized bond returns only.  None → 0 (Bocola
                             risk-only experiment: priced but never realized).
    init_D/init_F          : optional dicts with keys (n_prev, kappa_prev,
                             phi_bdom_prev, phi_bfor_prev, rdep_prev) — the
                             bank state carried into period 0 when the path
                             starts mid-crisis (default branch / policy runs).
                             None → steady-state values (current behaviour).
    Q_bD_lag0, Q_bF_lag0   : period −1 bond prices for realized returns at
                             t=0 (None → IC-consistent SS prices).
    p_lag0                 : period −1 real exchange rate (None → p_ss).

    Returns dict of paths: n_IC_D, n_D, rn_D, div_D, theta_D, Dep_supply_D,
    rb_D, and F analogues, plus pass-through of holdings b_D_D, b_F_F.
    """
    T = len(Kap_D)

    if def_real_D is None:
        def_real_D = np.zeros(T)
    if def_real_F is None:
        def_real_F = np.zeros(T)

    f_D   = cal["f_D"];           f_F   = cal["f_F"]
    lK_D  = cal["lambda_K_D"];    lK_F  = cal["lambda_K_F"]
    lbD_D = cal["lambda_bD_D"];   lbD_F = cal["lambda_bD_F"]
    lbF_D = cal["lambda_bF_D"];   lbF_F = cal["lambda_bF_F"]
    rec_D = cal["recovery_rate_D"]; rec_F = cal["recovery_rate_F"]
    db_D  = cal["delta_b_D"];     db_F  = cal["delta_b_F"]

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

    n_IC_D = np.empty(T);  n_ACCUM_D = np.empty(T)
    rn_D   = np.empty(T);  div_D     = np.empty(T)
    n_IC_F = np.empty(T);  n_ACCUM_F = np.empty(T)
    rn_F   = np.empty(T);  div_F     = np.empty(T)

    # D-bank forward state (SS by default; overridable for mid-path starts)
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
        p_t   = p_path[t]
        p_lag = p_path[t - 1] if t > 0 else p_l0

        # D-bank earns rb_D on D-bonds (D-goods) and rb_F on F-bonds (convert F→D goods).
        rb_D_t = rb_D_path[t]
        rb_F_t = (1.0 + rb_F_path[t]) * p_t / p_lag - 1   # F-goods → D-goods

        rn_D_t = (kappa_D_prev * (rk_D[t] - rdep_D_prev)
                  + phi_bdom_D_prev * (rb_D_t - rdep_D_prev)
                  + phi_bfor_D_prev * (rb_F_t - rdep_D_prev)
                  + rdep_D_prev)
        gross_D = (1 + rn_D_t) * n_D_prev
        total_assets_D = (Q_D[t] * Kap_D[t] + Q_bD_path[t] * b_D_D_path[t]
                          + p_t * Q_bF_path[t] * b_F_D_path[t])
        entrant_D   = cal["omega_ent_D"] * total_assets_D
        n_ACCUM_D_t = (1 - f_D) * gross_D + entrant_D
        div_D_t     = f_D * gross_D - entrant_D
        n_ACCUM_D[t] = n_ACCUM_D_t
        rn_D[t]  = rn_D_t
        div_D[t] = div_D_t

        # D-bank IC in D-goods (F-bonds enter as p×Q_bF)
        n_IC_D_t = (lK_D * Q_D[t] * Kap_D[t]
                    + lbD_D * Q_bD_path[t] * b_D_D_path[t]
                    + lbF_D * p_t * Q_bF_path[t] * b_F_D_path[t]) / alpha_D_path[t]
        n_IC_D[t] = n_IC_D_t

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
        entrant_F   = cal["omega_ent_F"] * total_assets_F
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
