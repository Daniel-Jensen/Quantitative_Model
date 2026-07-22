# TRADE BLOCK: CES CONSUMPTION AGGREGATOR AND BILATERAL TRADE FLOWS.
# p = D-goods per F-good.


def ces_price(p, cal, country="D"):
    # CES PRICE INDEX OF A COUNTRY'S CONSUMPTION BASKET.
    omega = cal["omega_home"]
    eta   = cal["epsilon_trade"]
    exp   = 1.0 - eta
    if country == "D":
        inside = omega + (1.0 - omega) * p ** exp
    else:
        inside = omega + (1.0 - omega) * (1.0 / p) ** exp   # F imports D-goods at 1/p
    return inside ** (1.0 / exp)


def import_demand(p, C, P_CES, cal, country="D"):
    # IMPORT VOLUME FROM THE CES DEMAND CURVE.
    omega = cal["omega_home"]
    eta   = cal["epsilon_trade"]
    if country == "D":
        return (1.0 - omega) * (P_CES / p) ** eta * C
    return (1.0 - omega) * (P_CES * p) ** eta * C


def trade_balance(p, IM_D, IM_F):
    # NET EXPORTS OF EACH COUNTRY, IN OWN-GOOD UNITS.
    NX_D = IM_F - p * IM_D
    NX_F = IM_D - IM_F / p
    return NX_D, NX_F
