# TRADE BLOCK: CES CONSUMPTION AGGREGATOR AND BILATERAL TRADE FLOWS.
# p = D-goods per F-good. Every quantity is PER CAPITA of its own country; the
# country masses size_D / size_F enter only where the two countries' flows are
# added together (trade_balance). Home bias is per-country and size-consistent:
# balanced trade at p = 1 needs size_D*(1-omega_D) = size_F*(1-omega_F), which
# calibration.py imposes when it derives omega_home_F. size_F = size_D reproduces
# the symmetric block exactly.


def _omega(cal, country):
    # HOME-GOODS WEIGHT FOR ONE COUNTRY, falling back to the legacy scalar key.
    return cal.get(f"omega_home_{country}", cal["omega_home"])


def ces_price(p, cal, country="D"):
    # CES PRICE INDEX OF A COUNTRY'S CONSUMPTION BASKET.
    omega = _omega(cal, country)
    eta   = cal["epsilon_trade"]
    exp   = 1.0 - eta
    if country == "D":
        inside = omega + (1.0 - omega) * p ** exp
    else:
        inside = omega + (1.0 - omega) * (1.0 / p) ** exp   # F imports D-goods at 1/p
    return inside ** (1.0 / exp)


def import_demand(p, C, P_CES, cal, country="D"):
    # IMPORT VOLUME FROM THE CES DEMAND CURVE, PER CAPITA OF THE IMPORTING COUNTRY.
    omega = _omega(cal, country)
    eta   = cal["epsilon_trade"]
    if country == "D":
        return (1.0 - omega) * (P_CES / p) ** eta * C
    return (1.0 - omega) * (P_CES * p) ** eta * C


def size_ratio(cal):
    # F's MASS RELATIVE TO D. One place, so no caller writes size_F/size_D by hand.
    return cal.get("size_F", 1.0) / cal.get("size_D", 1.0)


def trade_balance(p, IM_D, IM_F, cal=None):
    # NET EXPORTS OF EACH COUNTRY, IN OWN-GOOD UNITS, PER CAPITA OF THAT COUNTRY.
    # IM_D is D's per-capita import of F-goods (F units); IM_F is F's per-capita
    # import of D-goods (D units). D's exports are consumed by size_F F-agents per
    # size_D D-agents, so the cross-country leg carries the mass ratio. cal=None
    # keeps the old symmetric signature working (ratio 1).
    sz = 1.0 if cal is None else size_ratio(cal)
    NX_D = sz * IM_F - p * IM_D
    NX_F = IM_D / sz - IM_F / p
    return NX_D, NX_F
