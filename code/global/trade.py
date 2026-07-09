"Trade block: CES consumption aggregator and bilateral trade flows." 


import numpy as np


def ces_price(p, cal, country="D"):
    #CES price index for country's consummptions basket

    omega = cal["omega_home"]
    eta   = cal["epsilon_trade"]
    exp   = 1.0 - eta

    if country == "D":
        inside = omega + (1.0 - omega) * p ** exp
    else:
        # F uses F-goods (price=1) and D-goods (price=1/p in F-goods)
        #   P_F = [omega·1^(1-eta) + (1-omega)·(1/p)^(1-eta)]^(1/(1-eta))
        inside = omega + (1.0 - omega) * (1.0 / p) ** exp

    return inside ** (1.0 / exp)


def import_demand(p, C, P_CES, cal, country="D"):
  #Volume of imports
    omega = cal["omega_home"]
    eta   = cal["epsilon_trade"]

    if country == "D":
        return (1.0 - omega) * (P_CES / p) ** eta * C
    else:
        return (1.0 - omega) * (P_CES * p) ** eta * C


def trade_balance(p, IM_D, IM_F):
    #Net exports 
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
