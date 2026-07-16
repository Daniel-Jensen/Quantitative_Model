# **Trade block: CES consumption aggregator and bilateral trade flows.**
import numpy as np


def ces_price(p, cal, country="D"):
    # **CES price index of a country's consumption basket.**
    omega = cal["omega_home"]
    eta   = cal["epsilon_trade"]
    exp   = 1.0 - eta
    if country == "D":
        inside = omega + (1.0 - omega) * p ** exp
    else:
        inside = omega + (1.0 - omega) * (1.0 / p) ** exp   # F imports D-goods at 1/p
    return inside ** (1.0 / exp)


def import_demand(p, C, P_CES, cal, country="D"):
    # **Import volume from the CES demand curve.**
    omega = cal["omega_home"]
    eta   = cal["epsilon_trade"]
    if country == "D":
        return (1.0 - omega) * (P_CES / p) ** eta * C
    else:
        return (1.0 - omega) * (P_CES * p) ** eta * C


def trade_balance(p, IM_D, IM_F):
    # **Net exports of each country (in own-good units).**
    NX_D = IM_F - p * IM_D
    NX_F = IM_D - IM_F / p
    return NX_D, NX_F
