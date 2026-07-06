"""Two-country Gertler-Karadi financial intermediary block.

Each country has a bank that holds three assets:
  D-bank: domestic capital K_D, domestic D-bonds b_D_D, foreign F-bonds b_F_D
  F-bank: domestic capital K_F, domestic F-bonds b_F_F, foreign D-bonds b_D_F

All bond quantities and prices are expressed in D-good units (D is the
monetary-union numeraire). F-bank balance-sheet quantities convert to
F-goods by dividing by the real exchange rate p.

NOTE on p convention: p = price of F-goods in D-good units (D-goods per
F-good). An increase in p means F-goods are more expensive. This is the
OPPOSITE of the real exchange rate convention used in some macro texts
(where q = P*/P = cost of domestic in foreign). The trade.py formulas
and all FX conversions in this file are consistent with this definition.

Multi-asset incentive constraint (IC):
  D-bank (D-goods):
    lambda_K · Q_D · K_D + lambda_bD · Q_bD · b_D_D + lambda_bF · Q_bF · b_F_D ≤ alpha_D · n_D

  F-bank (F-goods, bonds divided by p):
    lambda_K · Q_F · K_F + lambda_bF · Q_bF · b_F_F/p + lambda_bD · Q_bD · b_D_F/p ≤ alpha_F · n_F

When the IC binds (it always does in deterministic perfect-foresight), the
linear value function V_t(n) = alpha_t · n (Bocola 2016 Result 1, extended
to three assets) gives closed-form backward-pass equations:

  Omega_{t+1}   = beta_inter · [(1−f) + f · alpha_{t+1}]
  mu_t          = Omega_{t+1} · (rk_{t+1} − rdep_t) / lambda_K    [capital FOC]
  alpha_t       = Omega_{t+1} · (1 + rdep_t) / (1 − mu_t)          [Bellman]
  Q_bD_t = surv_D_{t+1}·(delta_b_D + (1−delta_b_D)·Q_bD_{t+1}) / (1 + rdep_t + lambda_bD·mu_t/Omega_{t+1})
  Q_bF_t = surv_F_{t+1}·(delta_b_F + (1−delta_b_F)·Q_bF_{t+1}) / (1 + rdep_F_t + lambda_bF·mu_F_t/Omega_F_{t+1})
  surv_{t+1} = 1 − def_{t+1}·(1−recovery)·writeoff_enabled
  At SS or writeoff_enabled=0: surv=1, collapsing to Q_ss = delta_b/(rdep+delta_b+IC_spread).

The HM bond pricing here uses the IC-derived no-arbitrage relationship:
the excess return on each bond over the deposit rate must equal the IC
spread lambda_b · mu / Omega.  The survival factor surv_{t+1} prices in
expected write-off losses under perfect foresight; with writeoff_enabled=0
surv=1 and the formula reduces to the risk-free HM recursion.

Cross-border positions from portfolio adjustment-cost FOC:
  D-bank holds F-bonds: b_F_D_t from FOC
    rb_F_in_D_{t+1} − rdep_D_t = excess_return_F_D_ss + psi_bF_D · (b_F_D_t − b_F_D_ss)
                                  + lambda_bF_D · mu_D_t / Omega_D_{t+1}
  → b_F_D_t = b_F_D_ss + [rb_F_in_D_{t+1} − rdep_D_t
                            − excess_return_F_D_ss − lambda_bF_D · mu_D_t / Omega_D_{t+1}]
                           / psi_bF_D

  (and symmetrically for F-bank's D-bond holding b_D_F_t)

Bond market clearing (both banks hold government supply jointly):
  b_D_D_t + b_D_F_t = B_gov_D   (D-bonds: D-bank + F-bank = supply)
  b_F_D_t + b_F_F_t = B_gov_F   (F-bonds: D-bank + F-bank = supply)

Bocola-Dovis (2019) rollover-risk spread:
  lbD_D_eff = lbD_D + psi_bd_D · xi_{t+1}   (D-bank IC for D-bonds)
  lbD_F_eff = lbD_F + psi_bd_D · xi_{t+1}   (F-bank IC for D-bonds — contagion)
The sunspot xi_{t+1} ∈ [0,1] is the next-period run probability, supplied
exogenously from solve_transition_bd().  It tightens the IC when rollover
risk is high, lowering Q_bD.  With psi_bd_D=0 this collapses to base GK.

Net worth has two characterisations that must agree (outer residual):
  n_IC   = IC-binding allocation size (desired by bank given alpha)
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
    # This is the fixed point of the backward pricing recurrence used in solve_bank_paths.
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

    # Scale factor for F-bank: F-bank balance sheet in F-goods (divide by p)
    p_scale = 1.0 if country == "D" else p_ss

    # n from IC binding: n_IC = (lambda_K·Q_K·K + lambda_bD·Q_bD·b_dom + lambda_bF·Q_bF·b_for) / alpha
    # For F-bank, bonds are in D-goods; divide by p to get F-good IC value
    if country == "D":
        ic_numerator = (lambda_K * Kap_ss
                        + lambda_bD * Q_bdom_ss * b_dom_ss
                        + lambda_bF * Q_bfor_ss * b_for_ss)
    else:
        ic_numerator = (lambda_K * Kap_ss
                        + (lambda_bD * Q_bdom_ss * b_dom_ss
                           + lambda_bF * Q_bfor_ss * b_for_ss) / p_ss)

    n_ss_IC = ic_numerator / alpha_ss

    # n from forward accumulation: D·n = excess returns + entrant transfer
    # D = 1 − (1−f)·(1+rdep)
    D_val = 1.0 - (1 - f) * (1 + rdep_ss)
    if D_val <= 0:
        raise ValueError(f"[{country}] D={D_val} ≤ 0: no stationary net-worth rest point.")

    if country == "D":
        total_assets = Kap_ss + Q_bdom_ss * b_dom_ss + Q_bfor_ss * b_for_ss
        n_ss_ACCUM = (
            ((1 - f) * (rk_ss     - rdep_ss) + omega_ent) * Kap_ss
            + ((1 - f) * IC_spread_dom + omega_ent) * Q_bdom_ss * b_dom_ss
            + ((1 - f) * IC_spread_for + omega_ent) * Q_bfor_ss * b_for_ss
        ) / D_val
    else:
        # F-bank: assets in F-goods (divide bond values by p_ss)
        Kap_val   = Kap_ss
        bdom_val  = Q_bdom_ss * b_dom_ss / p_ss   # F-bonds in F-goods
        bfor_val  = Q_bfor_ss * b_for_ss / p_ss   # D-bonds in F-goods
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
        phi_bfor_ss = Q_bfor_ss * b_for_ss / n_ss
    else:
        kappa_ss    = Kap_ss / n_ss
        phi_bdom_ss = Q_bdom_ss * b_dom_ss / (p_ss * n_ss)
        phi_bfor_ss = Q_bfor_ss * b_for_ss / (p_ss * n_ss)

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
# Transition-path bank block
# ─────────────────────────────────────────────────────────────────────────────

def solve_bank_paths(Kap_D, Kap_F, Q_D, Q_F, rk_D, rk_F,
                     rdep_D, rdep_F, p_path, B_gov_D, B_gov_F,
                     cal, ss_bk_D, ss_bk_F,
                     def_D_path=None, def_F_path=None,
                     sunspot_D_path=None, sunspot_F_path=None):
    """Transition-path bank block for both countries simultaneously.

    Arguments
    ---------
    Kap_D, Kap_F   : capital paths (T,) — outer Newton unknowns
    Q_D,  Q_F     : Tobin's Q paths (T,) — from capital.py
    rk_D, rk_F   : capital return paths (T,) — from capital.py
    rdep_D, rdep_F: deposit rate paths (T,) — outer Newton unknowns
    p_path        : real exchange rate path (T,) — outer Newton unknown
    B_gov_D/F     : total government bond supply (scalars, fixed)
    cal           : calibration dict
    ss_bk_D/F    : steady-state bank dicts from steady_state_bank()
    def_D_path, def_F_path: default rate paths (T,) or None (→ zeros)

    Returns
    -------
    dict with keys (all length T):
      For D: alpha_D, mu_D, Q_bD, b_F_D, b_D_D, n_IC_D, n_D, rn_D, div_D,
             theta_D, Dep_supply_D, rb_D (realised return on D-bonds)
      For F: alpha_F, mu_F, Q_bF, b_D_F, b_F_F, n_IC_F, n_F, rn_F, div_F,
             theta_F, Dep_supply_F, rb_F
      Shared: rb_D (D-bond realised return to holders), rb_F (F-bond return)
    """
    T = len(Kap_D)

    if def_D_path is None:
        def_D_path = np.zeros(T)
    if def_F_path is None:
        def_F_path = np.zeros(T)

    # Pull calibration for each bank
    f_D        = cal["f_D"];          f_F        = cal["f_F"]
    bi_D       = cal["beta_inter_D"]; bi_F       = cal["beta_inter_F"]
    lK_D       = cal["lambda_K_D"];   lK_F       = cal["lambda_K_F"]
    lbD_D      = cal["lambda_bD_D"];  lbD_F      = cal["lambda_bD_F"]
    lbF_D      = cal["lambda_bF_D"];  lbF_F      = cal["lambda_bF_F"]
    om_D       = cal["omega_ent_D"];  om_F       = cal["omega_ent_F"]
    db_D       = cal["delta_b_D"];    db_F       = cal["delta_b_F"]
    psi_bFD    = cal["psi_bF_D"]      # D-bank adj cost for F-bond deviation
    psi_bDF    = cal["psi_bD_F"]      # F-bank adj cost for D-bond deviation
    b_F_D_ss   = cal["b_F_D_ss"]
    b_D_F_ss   = cal["b_D_F_ss"]
    exc_FD_ss  = cal["excess_return_F_D_ss"]
    exc_DF_ss  = cal["excess_return_D_F_ss"]
    rec_D      = cal["recovery_rate_D"]; rec_F = cal["recovery_rate_F"]
    psi_bd_D   = cal.get("psi_bd_D", 0.0)   # BD sunspot spread sensitivity (D-bonds)
    psi_bd_F   = cal.get("psi_bd_F", 0.0)   # BD sunspot spread sensitivity (F-bonds)

    # ── Backward pass ─────────────────────────────────────────────────────────
    # At each t, given alpha_{t+1} for both banks, solve closed-form:
    #   Omega_{t+1} = beta_inter·[(1−f) + f·alpha_{t+1}]
    #   mu_t        = Omega_{t+1}·(rk_{t+1} − rdep_t)/lambda_K
    #   alpha_t     = Omega_{t+1}·(1+rdep_t)/(1−mu_t)
    #   Q_bD_t      = 1/(rdep_D_t + delta_b_D + lambda_bD_D·mu_D_t/Omega_D_{t+1})
    #   Q_bF_t      = 1/(rdep_F_t + delta_b_F + lambda_bF_F·mu_F_t/Omega_F_{t+1})
    #   b_F_D_t     from portfolio-adj-cost FOC (linear, closed-form)
    #   b_D_F_t     same

    alpha_D_path = np.empty(T)
    mu_D_path    = np.empty(T)
    alpha_F_path = np.empty(T)
    mu_F_path    = np.empty(T)
    Q_bD_path    = np.empty(T)  # D-bond price path
    Q_bF_path    = np.empty(T)  # F-bond price path
    b_F_D_path   = np.empty(T)  # D-bank's F-bond holding
    b_D_F_path   = np.empty(T)  # F-bank's D-bond holding

    # Effective divertability paths (base + BD sunspot tightening)
    lbD_D_eff_path = np.empty(T)   # D-bank IC for D-bonds
    lbF_D_eff_path = np.empty(T)   # D-bank IC for F-bonds
    lbD_F_eff_path = np.empty(T)   # F-bank IC for D-bonds (contagion channel)
    lbF_F_eff_path = np.empty(T)   # F-bank IC for F-bonds

    alpha_D_next = ss_bk_D["alpha_ss"]
    alpha_F_next = ss_bk_F["alpha_ss"]

    # IC-consistent SS bond prices (terminal condition for backward pass).
    # Bond pricing FOC at SS: Q = delta_b / (rdep + delta_b + lambda_b*mu/Omega)
    # Use ss_bk mu_ss/Omega_ss for the IC spread (already solved in SS bank block).
    ic_spread_bD_D = ss_bk_D["lambda_bD"] * ss_bk_D["mu_ss"] / ss_bk_D["Omega_ss"]
    ic_spread_bF_F = ss_bk_F["lambda_bF"] * ss_bk_F["mu_ss"] / ss_bk_F["Omega_ss"]
    Q_bD_ss_val = cal["delta_b_D"] / (cal["r_dep_D_target"] + cal["delta_b_D"] + ic_spread_bD_D)
    Q_bF_ss_val = cal["delta_b_F"] / (cal["r_dep_F_target"] + cal["delta_b_F"] + ic_spread_bF_F)
    Q_bD_next   = Q_bD_ss_val
    Q_bF_next   = Q_bF_ss_val

    for t in range(T - 1, -1, -1):
        # At the terminal period, rk[T] is unknown; use rk[T-1] as the SS approximation.
        # Using cal["rk_D_guess"] instead would introduce an alpha error at t=T-1
        # because rk_guess ≠ rk_ss in general.
        rk_D_next = rk_D[t + 1] if t + 1 < T else rk_D[T - 1]
        rk_F_next = rk_F[t + 1] if t + 1 < T else rk_F[T - 1]
        # Next-period default rates (survival in bond pricing)
        def_D_next = def_D_path[t + 1] if t + 1 < T else 0.0
        def_F_next = def_F_path[t + 1] if t + 1 < T else 0.0
        # BD: next-period sunspot (coordination failure probability → IC tightening)
        xi_D_next = (sunspot_D_path[t + 1] if sunspot_D_path is not None and t + 1 < T else 0.0)
        xi_F_next = (sunspot_F_path[t + 1] if sunspot_F_path is not None and t + 1 < T else 0.0)

        # ── D-bank backward step ──
        Omega_D     = bi_D * ((1 - f_D) + f_D * alpha_D_next)
        mu_D        = Omega_D * (rk_D_next - rdep_D[t]) / lK_D
        if mu_D >= 1.0:
            raise RuntimeError(f"D-bank mu_D={mu_D:.4f} ≥ 1 at t={t}; IC infeasible.")
        alpha_D     = Omega_D * (1 + rdep_D[t]) / (1 - mu_D)

        # BD: sunspot tightens IC for sovereign bonds (rollover risk premium)
        lbD_D_eff = lbD_D + psi_bd_D * xi_D_next   # D-bank IC for D-bonds
        lbF_D_eff = lbF_D + psi_bd_F * xi_F_next   # D-bank IC for F-bonds

        # HM pricing: Q_bD = surv_{t+1}·(db + (1-db)·Q_next) / (1 + rdep + IC_spread)
        # Survival from fundamental default only; BD sunspot enters denominator.
        ic_spread_bD_D   = lbD_D_eff * mu_D / Omega_D
        surv_D_for_price = 1.0 - def_D_next * (1.0 - rec_D)
        Q_bD = surv_D_for_price * (db_D + (1 - db_D) * Q_bD_next) / (1 + rdep_D[t] + ic_spread_bD_D)

        # ── F-bank backward step ──
        Omega_F     = bi_F * ((1 - f_F) + f_F * alpha_F_next)
        mu_F        = Omega_F * (rk_F_next - rdep_F[t]) / lK_F
        if mu_F >= 1.0:
            raise RuntimeError(f"F-bank mu_F={mu_F:.4f} ≥ 1 at t={t}; IC infeasible.")
        alpha_F     = Omega_F * (1 + rdep_F[t]) / (1 - mu_F)

        # BD: F-bank's IC for D-bonds also tightens with D-sunspot (contagion channel)
        lbD_F_eff = lbD_F + psi_bd_D * xi_D_next   # F-bank IC for D-bonds
        lbF_F_eff = lbF_F + psi_bd_F * xi_F_next   # F-bank IC for F-bonds

        # F-bank FOC for F-bonds (symmetric)
        ic_spread_bF_F   = lbF_F_eff * mu_F / Omega_F
        surv_F_for_price = 1.0 - def_F_next * (1.0 - rec_F)
        Q_bF = surv_F_for_price * (db_F + (1 - db_F) * Q_bF_next) / (1 + rdep_F[t] + ic_spread_bF_F)

        # ── Cross-border positions from portfolio adj-cost FOC ──
        # Survival for expected return at t+1 uses next-period default rate (B-2 fix).
        def_F_t1 = def_F_path[t + 1] if t + 1 < T else 0.0
        def_D_t1 = def_D_path[t + 1] if t + 1 < T else 0.0
        surv_F = 1.0 - def_F_t1 * (1.0 - rec_F)
        surv_D = 1.0 - def_D_t1 * (1.0 - rec_D)

        rb_F_raw = (db_F * surv_F + (1 - db_F) * Q_bF_next * surv_F) / Q_bF - 1
        # Convert F-bond return to D-goods using p_{t+1}/p_t
        p_next = p_path[t + 1] if t + 1 < T else p_path[t]  # terminal: p constant
        rb_F_in_D = (1 + rb_F_raw) * p_next / p_path[t] - 1

        ic_required_bF_D = lbF_D_eff * mu_D / Omega_D
        b_F_D_t = (b_F_D_ss
                   + (rb_F_in_D - rdep_D[t] - exc_FD_ss - ic_required_bF_D)
                   / psi_bFD)

        # F-bank's D-bond FOC:
        rb_D_raw = (db_D * surv_D + (1 - db_D) * Q_bD_next * surv_D) / Q_bD - 1
        # Convert D-bond return to F-goods using p_t/p_{t+1}
        rb_D_in_F = (1 + rb_D_raw) * p_path[t] / p_next - 1

        ic_required_bD_F = lbD_F_eff * mu_F / Omega_F
        b_D_F_t = (b_D_F_ss
                   + (rb_D_in_F - rdep_F[t] - exc_DF_ss - ic_required_bD_F)
                   / psi_bDF)

        alpha_D_path[t] = alpha_D;  mu_D_path[t] = mu_D
        alpha_F_path[t] = alpha_F;  mu_F_path[t] = mu_F
        Q_bD_path[t]    = Q_bD;     Q_bF_path[t] = Q_bF
        b_F_D_path[t]   = b_F_D_t
        b_D_F_path[t]   = b_D_F_t
        lbD_D_eff_path[t] = lbD_D_eff
        lbF_D_eff_path[t] = lbF_D_eff
        lbD_F_eff_path[t] = lbD_F_eff
        lbF_F_eff_path[t] = lbF_F_eff

        alpha_D_next = alpha_D;  alpha_F_next = alpha_F
        Q_bD_next    = Q_bD;     Q_bF_next    = Q_bF

    # Bond market clearing — B_gov may be scalar (SS) or time-varying path (BD/CK)
    B_gov_D_arr = np.broadcast_to(B_gov_D, T).copy() if np.ndim(B_gov_D) == 0 else np.asarray(B_gov_D)
    B_gov_F_arr = np.broadcast_to(B_gov_F, T).copy() if np.ndim(B_gov_F) == 0 else np.asarray(B_gov_F)
    b_D_D_path = B_gov_D_arr - b_D_F_path   # D-bonds held by D-bank
    b_F_F_path = B_gov_F_arr - b_F_D_path   # F-bonds held by F-bank

    # Realised returns at t on bonds bought at t-1
    Q_bD_lag = np.concatenate(([Q_bD_ss_val], Q_bD_path[:-1]))
    Q_bF_lag = np.concatenate(([Q_bF_ss_val], Q_bF_path[:-1]))

    surv_D_path = 1.0 - def_D_path * (1.0 - rec_D)
    surv_F_path = 1.0 - def_F_path * (1.0 - rec_F)
    rb_D_path   = (db_D * surv_D_path + (1 - db_D) * Q_bD_path * surv_D_path) / Q_bD_lag - 1
    rb_F_path   = (db_F * surv_F_path + (1 - db_F) * Q_bF_path * surv_F_path) / Q_bF_lag - 1

    # ── Forward pass: accumulate n_ACCUM and compute n_IC ─────────────────────
    n_IC_D   = np.empty(T)
    n_ACCUM_D = np.empty(T)
    rn_D     = np.empty(T)
    div_D    = np.empty(T)

    n_IC_F   = np.empty(T)
    n_ACCUM_F = np.empty(T)
    rn_F     = np.empty(T)
    div_F    = np.empty(T)

    # D-bank forward state
    n_D_prev       = ss_bk_D["n_ss"]
    kappa_D_prev   = ss_bk_D["kappa_ss"]
    phi_bdom_D_prev = ss_bk_D["phi_bdom_ss"]
    phi_bfor_D_prev = ss_bk_D["phi_bfor_ss"]
    rdep_D_prev    = cal["r_dep_D_target"]

    # F-bank forward state
    n_F_prev       = ss_bk_F["n_ss"]
    kappa_F_prev   = ss_bk_F["kappa_ss"]
    phi_bdom_F_prev = ss_bk_F["phi_bdom_ss"]
    phi_bfor_F_prev = ss_bk_F["phi_bfor_ss"]
    rdep_F_prev    = cal["r_dep_F_target"]

    for t in range(T):
        p_t   = p_path[t]
        # p_lag: exchange rate at t-1, needed for FX conversions of realised bond returns.
        # Q_bF is denominated in F-goods (F-bank prices F-bonds against rdep_F in F-goods);
        # Q_bD is denominated in D-goods (D-bank prices D-bonds against rdep_D in D-goods).
        p_lag = p_path[t - 1] if t > 0 else cal.get("p_ss", 1.0)

        # D-bank earns rb_D on D-bonds (D-goods, no conversion) and rb_F on F-bonds.
        # rb_F_path is a F-good return (Q_bF in F-goods); convert to D-goods:
        #   D-good return = (1 + rb_F_path) * p_t / p_lag  (bought at p_lag per F-good, sold at p_t)
        rb_D_t = rb_D_path[t]
        rb_F_t = (1.0 + rb_F_path[t]) * p_t / p_lag - 1   # F-goods → D-goods

        # D-bank net worth return (all terms in D-goods)
        rn_D_t = (kappa_D_prev * (rk_D[t] - rdep_D_prev)
                  + phi_bdom_D_prev * (rb_D_t - rdep_D_prev)
                  + phi_bfor_D_prev * (rb_F_t - rdep_D_prev)   # rb_F_t already in D-goods
                  + rdep_D_prev)
        gross_D = (1 + rn_D_t) * n_D_prev
        # Entrant transfer: proportional to total assets in D-goods
        total_assets_D = Q_D[t] * Kap_D[t] + Q_bD_path[t] * b_D_D_path[t] + Q_bF_path[t] * b_F_D_path[t]
        entrant_D  = cal["omega_ent_D"] * total_assets_D
        n_ACCUM_D_t = (1 - f_D) * gross_D + entrant_D
        div_D_t    = f_D * gross_D - entrant_D
        n_ACCUM_D[t] = n_ACCUM_D_t
        rn_D[t]  = rn_D_t
        div_D[t] = div_D_t

        # D-bank IC: uses effective lambda (base + BD sunspot tightening)
        n_IC_D_t = (lK_D * Q_D[t] * Kap_D[t]
                    + lbD_D_eff_path[t] * Q_bD_path[t] * b_D_D_path[t]
                    + lbF_D_eff_path[t] * Q_bF_path[t] * b_F_D_path[t]) / alpha_D_path[t]
        n_IC_D[t] = n_IC_D_t

        # Update D-bank portfolio ratios for next forward step
        kappa_D_prev    = Q_D[t] * Kap_D[t] / n_IC_D_t
        phi_bdom_D_prev = Q_bD_path[t] * b_D_D_path[t] / n_IC_D_t
        phi_bfor_D_prev = Q_bF_path[t] * b_F_D_path[t] / n_IC_D_t
        n_D_prev        = n_ACCUM_D_t
        rdep_D_prev     = rdep_D[t]

        # ── F-bank ─────────────────────────────────────────────────────────
        # Q_bF denominated in F-goods: rb_F_path is already a F-good return → use directly.
        # Q_bD denominated in D-goods: convert D-bond return to F-goods via p_lag/p_t.
        rb_F_fg_t = rb_F_path[t]                                # F-bonds: F-goods ✓
        rb_D_fg_t = (1 + rb_D_path[t]) * p_lag / p_t - 1      # D-bonds: D-goods → F-goods

        rn_F_t = (kappa_F_prev * (rk_F[t] - rdep_F_prev)
                  + phi_bdom_F_prev * (rb_F_fg_t - rdep_F_prev)   # F-bonds in F-goods
                  + phi_bfor_F_prev * (rb_D_fg_t - rdep_F_prev)   # D-bonds in F-goods
                  + rdep_F_prev)
        gross_F = (1 + rn_F_t) * n_F_prev
        # Total F-bank assets in F-goods
        total_assets_F = (Q_F[t] * Kap_F[t]
                          + Q_bF_path[t] * b_F_F_path[t] / p_t
                          + Q_bD_path[t] * b_D_F_path[t] / p_t)
        entrant_F   = cal["omega_ent_F"] * total_assets_F
        n_ACCUM_F_t = (1 - f_F) * gross_F + entrant_F
        div_F_t     = f_F * gross_F - entrant_F
        n_ACCUM_F[t] = n_ACCUM_F_t
        rn_F[t]  = rn_F_t
        div_F[t] = div_F_t

        # F-bank IC in F-goods: effective lambda (D-bond term uses contagion channel)
        n_IC_F_t = (lK_F * Q_F[t] * Kap_F[t]
                    + lbF_F_eff_path[t] * Q_bF_path[t] * b_F_F_path[t] / p_t
                    + lbD_F_eff_path[t] * Q_bD_path[t] * b_D_F_path[t] / p_t) / alpha_F_path[t]
        n_IC_F[t] = n_IC_F_t

        kappa_F_prev    = Q_F[t] * Kap_F[t] / n_IC_F_t
        phi_bdom_F_prev = Q_bF_path[t] * b_F_F_path[t] / (p_t * n_IC_F_t)
        phi_bfor_F_prev = Q_bD_path[t] * b_D_F_path[t] / (p_t * n_IC_F_t)
        n_F_prev        = n_ACCUM_F_t
        rdep_F_prev     = rdep_F[t]

    theta_D = ((Q_D * Kap_D + Q_bD_path * b_D_D_path + Q_bF_path * b_F_D_path)
               / n_ACCUM_D)
    theta_F = ((Q_F * Kap_F + Q_bF_path * b_F_F_path / p_path + Q_bD_path * b_D_F_path / p_path)
               / n_ACCUM_F)

    Dep_supply_D = (theta_D - 1) * n_ACCUM_D
    Dep_supply_F = (theta_F - 1) * n_ACCUM_F

    return dict(
        # D-bank
        alpha_D=alpha_D_path, mu_D=mu_D_path,
        n_IC_D=n_IC_D, n_D=n_ACCUM_D, rn_D=rn_D, div_D=div_D,
        theta_D=theta_D, Dep_supply_D=Dep_supply_D,
        b_D_D=b_D_D_path, b_F_D=b_F_D_path,
        # F-bank
        alpha_F=alpha_F_path, mu_F=mu_F_path,
        n_IC_F=n_IC_F, n_F=n_ACCUM_F, rn_F=rn_F, div_F=div_F,
        theta_F=theta_F, Dep_supply_F=Dep_supply_F,
        b_D_F=b_D_F_path, b_F_F=b_F_F_path,
        # Shared bond prices and returns
        Q_bD=Q_bD_path, Q_bF=Q_bF_path,
        rb_D=rb_D_path, rb_F=rb_F_path,
    )
