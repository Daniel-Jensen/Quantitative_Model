"Government block: Hatchondo-Martinez (2009) geometric-decay perpetuity bonds, Bohn (1998) fiscal rule, Cole-Kehoe (2000) self-fulfilling default crisis zones."
import numpy as np


# ── Steady-state helpers ──────────────────────────────────────────────────────

def hm_bond_price_ss(rdep_ss, delta_b):
    """HM perpetuity price at a risk-free steady state.

    Bond pays coupon delta_b per period, decays at rate (1-delta_b).
    No-arbitrage SS price: Q_B_ss = delta_b / (rdep_ss + delta_b).
    """
    return delta_b / (rdep_ss + delta_b)


def hm_bond_return_ss(Q_B_ss, delta_b):
    """Realised return at SS on HM perpetuity (= rdep_ss by no-arbitrage)."""
    return (delta_b + (1 - delta_b) * Q_B_ss) / Q_B_ss - 1


def govt_steady_state(cal, rdep_ss, country):
    """Steady-state government block.

    At SS: def_rate = 0, bond stock = B_gov_ss (exogenous), no net issuance.

    Returns dict: Q_B_ss, rb_ss, Tax_ss, b_gov_ss, coupon_ss.
    """
    delta_b  = cal[f"delta_b_{country}"]
    B_gov_ss = cal[f"B_gov_{country}_ss"]
    G        = cal[f"G_{country}"]

    Q_B_ss    = hm_bond_price_ss(rdep_ss, delta_b)
    rb_ss     = hm_bond_return_ss(Q_B_ss, delta_b)
    coupon_ss = delta_b * B_gov_ss
    # Government budget at SS: G + coupon = Tax + issuance_proceeds
    # Issuance = delta_b*B_gov new bonds at price Q_B_ss
    Tax_ss    = G + coupon_ss * (1.0 - Q_B_ss)   # = G + delta_b*B_gov*rdep/(rdep+delta_b)

    return dict(Q_B_ss=Q_B_ss, rb_ss=rb_ss,
                Tax_ss=Tax_ss, b_gov_ss=B_gov_ss, coupon_ss=coupon_ss)


# ── Cole-Kehoe default probability ───────────────────────────────────────────

def ck_default_prob(b_gov, Y_ss, cal, sunspot, country):
    """Cole-Kehoe (2000) default probability for a single period.

    Three zones based on debt-to-output ratio b_gov/Y_ss:
      Safe zone   (b/Y < b_ck_low):           returns 0.
      Crisis zone (b_ck_low ≤ b/Y < b_ck_high): returns sunspot ∈ [0,1].
      Certain-default (b/Y ≥ b_ck_high):      returns 1.

    `sunspot` is the exogenous probability of default conditional on the crisis
    zone being active.  In a perfect-foresight MIT-shock setting, sunspot is an
    AR(1) path that starts at some peak and decays; it replaces the old AR(1)
    def_D_path.
    """
    b_low  = cal[f"b_ck_low_{country}"]
    b_high = cal[f"b_ck_high_{country}"]
    b_y    = b_gov / Y_ss
    return float(np.where(b_y < b_low, 0.0,
                          np.where(b_y >= b_high, 1.0, sunspot)))


# ── Bocola-Dovis fundamental default probability ─────────────────────────────

def bd_fundamental_default_prob(b_gov, Y_ss, cal, country):
    """BD fundamental default: only when debt exceeds the hard upper bound b_ck_high.

    Unlike Cole-Kehoe, there is no crisis zone — the sunspot enters through the
    IC spread in bank.py.  This function purely captures the solvency threshold
    (Bocola-Dovis B̄): certain default once debt is unsustainable regardless of
    lender beliefs.
    """
    b_high = cal.get(f"b_ck_high_{country}", 99.0)
    return 1.0 if b_gov / Y_ss >= b_high else 0.0


# ── Endogenous debt accumulation ─────────────────────────────────────────────

def integrate_b_gov(b_gov0, Q_bD_path, def_rate_path, Tax_path, cal, country):
    """Forward-integrate government debt stock from the HM budget identity.

    At each t the government:
      - Has outstanding bonds b_gov[t] (beginning of period, = b_gov0 at t=0).
      - Fraction delta_b matures; outstanding stock falls by factor (1-delta_b).
      - Default: fraction def_rate[t] of bonds partially haircut; effective
        survival surv[t] = 1 - def_rate[t]*(1 - recovery_rate).
      - Budget: Tax[t] - G = coupon[t] - issuance_proceeds[t]
          coupon[t]    = delta_b * b_gov[t] * surv[t]
          new_bonds[t] = (G + coupon[t] - Tax[t]) / Q_bD[t]
          b_gov[t+1]   = (1-delta_b)*b_gov[t]*surv[t] + new_bonds[t]

    Verification at SS: Tax_ss = G + delta_b*B_gov_ss*(1-Q_bD_ss)
      → b_gov[1] = (1-db)*B_ss + (G + db*B_ss - Tax_ss)/Q_ss = B_ss. ✓

    Parameters
    ----------
    b_gov0        : scalar, initial stock (= B_gov_ss at period 0)
    Q_bD_path     : (T,) bond price path from bank backward pass
    def_rate_path : (T,) default probability path (from CK outer loop)
    Tax_path      : (T,) tax revenue path (from govt_transition)
    cal           : calibration dict
    country       : "D" or "F"

    Returns
    -------
    b_gov_path : (T,) array — stocks at beginning of periods 1..T
    """
    delta_b       = cal[f"delta_b_{country}"]
    recovery_rate = cal[f"recovery_rate_{country}"]
    G             = cal[f"G_{country}"]
    T             = len(Q_bD_path)

    b_gov_path = np.empty(T)
    b_gov      = float(b_gov0)

    for t in range(T):
        surv_t    = 1.0 - def_rate_path[t] * (1.0 - recovery_rate)
        coupon_t  = delta_b * b_gov * surv_t
        new_bonds = (G + coupon_t - Tax_path[t]) / Q_bD_path[t]
        b_gov     = (1.0 - delta_b) * b_gov * surv_t + new_bonds
        b_gov_path[t] = b_gov

    return b_gov_path


# ── Transition-path government block ─────────────────────────────────────────

def bohn_tax(b_gov, b_gov_ss, Tax_ss, phi_lamb):
    """Bohn fiscal rule: lump-sum tax rises when debt exceeds steady state."""
    return Tax_ss + phi_lamb * (b_gov - b_gov_ss)


def govt_transition(cal, gs, Q_B_path, def_rate_path, b_gov_path, country, p_path=None):
    """Government flows along a transition path.

    Parameters
    ----------
    gs           : steady-state government dict (from govt_steady_state).
    Q_B_path     : (T,) bond price path from bank backward pass.
    def_rate_path: (T,) default probability path from CK outer loop.
    b_gov_path   : (T,) beginning-of-period debt stock path.
                   Entry t is the stock carried into period t (= b_gov_ss initially).
    p_path       : (T,) real exchange rate (required for country="F").

    The Bohn fiscal rule adjusts taxes in response to lagged debt:
      Tax[t] = Tax_ss + phi_lamb * (b_gov_path[t] - b_gov_ss)

    Survival factor always active (no writeoff gate):
      surv[t] = 1 - def_rate[t] * (1 - recovery_rate)

    Currency convention:
      All bonds are D-good claims; coupons and issuance proceeds are D-good flows.
      For country="F" the D-good flows are converted to F-goods via p_path so
      that Tax_F is in F-goods, consistent with F-household income accounting.

    Returns dict: Tax, coupon, net_issuance, b_gov — all shape (T,), own-good units.
    """
    delta_b       = cal[f"delta_b_{country}"]
    recovery_rate = cal[f"recovery_rate_{country}"]
    phi_lamb      = cal[f"phi_lamb_{country}"]
    Tax_ss        = gs["Tax_ss"]
    b_gov_ss      = gs["b_gov_ss"]

    T = len(Q_B_path)
    if p_path is None:
        p_path = np.ones(T)

    # Bohn rule: tax responds to beginning-of-period (lagged) debt
    Tax_D     = Tax_ss + phi_lamb * (b_gov_path - b_gov_ss)    # D-goods, shape (T,)

    # Survival factor (always active in CK framework)
    surv      = 1.0 - def_rate_path * (1.0 - recovery_rate)

    # D-good coupon and roll-over issuance
    coupon_D  = delta_b * b_gov_path * surv
    net_iss_D = delta_b * b_gov_path * Q_B_path

    if country == "F":
        # Convert D-good bond flows to F-goods. p = price of F-goods in D-goods.
        Tax_out      = Tax_D     / p_path
        coupon_out   = coupon_D  / p_path
        net_iss_out  = net_iss_D / p_path
    else:
        Tax_out     = Tax_D
        coupon_out  = coupon_D
        net_iss_out = net_iss_D

    return dict(Tax=Tax_out, coupon=coupon_out, net_issuance=net_iss_out,
                b_gov=b_gov_path)
