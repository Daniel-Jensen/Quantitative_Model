"""Trade block: CES consumption aggregator and bilateral trade flows.

Real exchange rate: p = price of F-goods in D-good units (D-goods per F-good).
  - D is the monetary-union numeraire for bonds and asset prices.
  - p > 1 means F-goods are expensive relative to D-goods.
  - When p rises, D-country imports become more costly; F-country terms-of-trade improve.
  - At the symmetric steady state: p = 1.

CES price indices (Dixit-Stiglitz):
  P_CES_D = [omega_home + (1 − omega_home) · p^(1−eta)]^[1/(1−eta)]
  P_CES_F = [omega_home + (1 − omega_home) · (1/p)^(1−eta)]^[1/(1−eta)]

where omega_home = home-good weight, eta = epsilon_trade = trade elasticity.

Import demand (in units of the imported good):
  IM_D = (1 − omega_home) · (P_CES_D / p)^eta · C_D   (D imports F-goods)
  IM_F = (1 − omega_home) · (P_CES_F · p)^eta · C_F   (F imports D-goods)

Trade balance (in D-goods):
  NX_D = IM_F − p · IM_D   (receipts from F imports minus D's import bill)
  NX_F = IM_D − IM_F / p   (symmetric: Walras identity NX_D + p·NX_F = 0)

All functions accept scalar or NumPy arrays for p, C, etc. to support
both steady-state evaluation and vectorised transition paths.
"""
import numpy as np


def ces_price(p, cal, country="D"):
    """CES price index for country's consumption basket.

    p: real exchange rate (price of F-good in D-goods), scalar or array.
    """
    omega = cal["omega_home"]
    eta   = cal["epsilon_trade"]
    exp   = 1.0 - eta

    if country == "D":
        # D uses D-goods (price=1) and F-goods (price=p in D-goods, so 1 per F-good)
        # CES over D-good and F-good in D-good units:
        #   P_D = [omega·1^(1-eta) + (1-omega)·p^(1-eta)]^(1/(1-eta))
        inside = omega + (1.0 - omega) * p ** exp
    else:
        # F uses F-goods (price=1) and D-goods (price=1/p in F-goods)
        #   P_F = [omega·1^(1-eta) + (1-omega)·(1/p)^(1-eta)]^(1/(1-eta))
        inside = omega + (1.0 - omega) * (1.0 / p) ** exp

    return inside ** (1.0 / exp)


def import_demand(p, C, P_CES, cal, country="D"):
    """Volume of imports demanded by `country` (in units of the foreign good).

    D imports F-goods: IM_D = (1−omega)·(P_CES_D/p)^eta · C_D
    F imports D-goods: IM_F = (1−omega)·(P_CES_F·p)^eta · C_F

    The relative-price term adjusts for how expensive the import is in
    terms of the domestic price index.
    """
    omega = cal["omega_home"]
    eta   = cal["epsilon_trade"]

    if country == "D":
        return (1.0 - omega) * (P_CES / p) ** eta * C
    else:
        return (1.0 - omega) * (P_CES * p) ** eta * C


def trade_balance(p, IM_D, IM_F):
    """Net exports in D-good units for each country.

    NX_D = value of F-goods imported by F (= IM_F D-goods sold to F)
           minus D's import bill (IM_D F-goods × p D-goods/F-good).
    NX_F is the symmetric expression in F-good units.

    Walras identity: NX_D + p·NX_F = 0.
    """
    NX_D = IM_F - p * IM_D
    NX_F = IM_D - IM_F / p
    return NX_D, NX_F


def external_account(NX_D, Q_bF, b_F_D, Q_bD, b_D_F,
                     Q_bF_lag, b_F_D_lag, Q_bD_lag, b_D_F_lag,
                     rb_F, rb_D):
    """Current account residual for country D (Walras-redundant diagnostic).

    CA_D = NX_D + interest receipts on F-bonds - interest payments on D-bonds
           - change in net foreign assets (NFA).

    This should equal zero to machine precision after the solver converges;
    it is NOT imposed as a residual (Walras-redundant).
    """
    receipts = (1.0 + rb_F) * Q_bF_lag * b_F_D_lag
    payments  = (1.0 + rb_D) * Q_bD_lag * b_D_F_lag
    nfa_now   = Q_bF * b_F_D - Q_bD * b_D_F
    nfa_lag   = Q_bF_lag * b_F_D_lag - Q_bD_lag * b_D_F_lag
    ca_resid  = NX_D + receipts - payments - (nfa_now - nfa_lag)
    return ca_resid
