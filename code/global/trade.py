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
