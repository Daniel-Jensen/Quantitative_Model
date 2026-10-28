# THE REST POINT UNDER EXACT CLEARING: THE ZERO-SHOCK PATH FROM THE FITTED REST POINT, CLEARED EXACTLY EACH QUARTER.
#   WT=dff4 python3 exact_rest.py <tag> [max quarters]
import sys
ARGV = list(sys.argv)
sys.argv = ["x"]
import ha_path as hp
if __name__ == "__main__":
    import json
    import numpy as np
    from global_projection.solver_recursive.recursive_experiment import (stochastic_rest_point, read_exact,
                                                                         read_at, advance)
    import tpi_national as tn
    tag, Tmax = ARGV[1], int(ARGV[2]) if len(ARGV) > 2 else 400
    c = hp.context(tag)
    cal0 = c["cal"]
    keys = ("Y_D", "C_D", "I_D", "N_D", "n_D", "mu_D", "rdep_D", "r_wc_D", "Q_bD", "Q_bF", "p", "Tax_D",
            "B_D", "A_D", "m_cb", "M_cb_new", "Y_F", "C_F")
    res = {}
    for econ in ("off", "on"):
        cal = dict(cal0, tpi_on=(econ == "on"))
        r = c["rules"][econ]
        S = stochastic_rest_point(r, cal, c["ss"], c["sproc"], verbose=False)
        fit = read_at(r, cal, c["ss"], c["sproc"], S.copy())
        o, x0 = None, None
        for t in range(Tmax):
            o = read_exact(r, cal, c["ss"], c["sproc"], S, x0)
            x0 = o["_x"]
            Sn = advance(o, S, c["sproc"], r.grid)
            d = float(np.max(np.abs(Sn - S)))
            S = Sn
            if d < 1e-10:
                break
        legs = tn.price_legs(r, cal, c["ss"], c["sproc"], S)
        db, dbF = cal["delta_b_D"], cal["delta_b_F"]
        sov = lambda q, qF: 4e4 * (db * (1 - q) / q - dbF * (1 - qF) / qF)
        keys = tuple(k for k in keys if k in o and k in fit)
        res[econ] = dict(quarters=t + 1, dS=d, exact={k: float(np.squeeze(o[k])) for k in keys},
                         fitted={k: float(np.squeeze(fit[k])) for k in keys}, legs=legs,
                         sov_exact=sov(float(o["Q_bD"]), float(o["Q_bF"])),
                         sov_fitted=sov(float(fit["Q_bD"]), float(fit["Q_bF"])))
        print(f"[{tag} {econ}] exact rest after {t+1} quarters (|dS| {d:.1e}); sov spread fitted@S0 "
              f"{res[econ]['sov_fitted']:.1f} -> exact rest {res[econ]['sov_exact']:.1f} bp; legs EL "
              f"{legs['exp_loss_bp']:.1f} RP {legs['risk_prem_bp']:+.1f} LP {legs['liq_prem_bp']:.1f} mu {legs['mu']:.4f}",
              flush=True)
    a, b = res["off"]["exact"], res["on"]["exact"]
    for k in keys:
        rate = k.startswith(("rdep", "r_wc"))
        dv = 4e4 * (b[k] - a[k]) if rate else (100 * (b[k] / a[k] - 1) if abs(a[k]) > 1e-9 else b[k] - a[k])
        print(f"  {k:9s} off {a[k]:+.6f}  on {b[k]:+.6f}   {'bp/yr' if rate else '%'} {dv:+.3f}")
    json.dump(res, open(f"tpi_runs/{tag}_exact_rest.json", "w"), indent=1, default=float)
