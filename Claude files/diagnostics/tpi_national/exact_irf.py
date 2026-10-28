# THE HEADLINE-SHOCK IRF STARTED FROM THE EXACT-CLEARING REST POINT, AGAINST THE ONE STARTED FROM THE FITTED REST POINT.
#   WT=dff4 python3 exact_irf.py <tag> <econ> [T] [rest quarters]
import sys
ARGV = list(sys.argv)
sys.argv = ["x"]
import ha_path as hp
if __name__ == "__main__":
    import json
    import numpy as np
    from global_projection.solver_recursive.recursive_experiment import (stochastic_rest_point, read_exact,
                                                                         advance, s_from_pd)
    from global_projection.solver_recursive.state_grid import IS
    tag, econ = ARGV[1], ARGV[2]
    T = int(ARGV[3]) if len(ARGV) > 3 else 25
    Tr = int(ARGV[4]) if len(ARGV) > 4 else 600
    c = hp.context(tag)
    cal = dict(c["cal"], tpi_on=(econ == "on"))
    r, ss, sp = c["rules"][econ], c["ss"], c["sproc"]
    S0 = stochastic_rest_point(r, cal, ss, sp, verbose=False)

    def path(S, T, x0=None):
        outs = []
        for _ in range(T):
            o = read_exact(r, cal, ss, sp, S, x0)
            x0 = o["_x"]
            outs.append(o)
            S = advance(o, S, sp, r.grid)
        return outs, S

    _, Sx = path(S0.copy(), Tr)
    keys = ("Y_D", "C_D", "I_D", "n_D", "Q_bD", "rdep_D", "r_wc_D", "mu_D", "m_cb", "p")
    db, dbF = cal["delta_b_D"], cal["delta_b_F"]

    def irf(Sstart):
        Ss = Sstart.copy(); Ss[IS] = s_from_pd(0.02)
        sh, _ = path(Ss, T); rf, _ = path(Sstart.copy(), T)
        g = lambda outs, k: np.array([float(o[k]) for o in outs])
        out = {}
        for k in ("Y_D", "C_D", "I_D", "n_D", "Q_bD"):
            out[k] = 100 * (g(sh, k) / g(rf, k) - 1)
        out["spread"] = 4e4 * ((g(sh, "r_wc_D") - g(sh, "rdep_D")) - (g(rf, "r_wc_D") - g(rf, "rdep_D")))
        sov = lambda outs: 4e4 * (db * (1 - g(outs, "Q_bD")) / g(outs, "Q_bD") - dbF * (1 - g(outs, "Q_bF")) / g(outs, "Q_bF"))
        out["sov_level"] = sov(sh); out["sov_ref"] = sov(rf)
        out["m_cb"] = 100 * g(sh, "m_cb") / cal["B_gov_D_ss"]
        out["mu"] = g(sh, "mu_D"); out["mu_ref"] = g(rf, "mu_D")
        return out

    res = {"fitted_start": irf(S0), "exact_start": irf(Sx)}
    for lab, d in res.items():
        print(f"[{tag} {econ}] {lab}: Y " + " ".join(f"{x:+.3f}" for x in d["Y_D"][:6])
              + f" | C0 {d['C_D'][0]:+.3f} cumC {d['C_D'].sum():+.2f} cumY {d['Y_D'].sum():+.2f} n0 {d['n_D'][0]:+.2f}"
              f" sov0 {d['sov_level'][0]:.1f} (ref {d['sov_ref'][0]:.1f}) dspread0 {d['spread'][0]:+.1f} m0 {d['m_cb'][0]:.1f}%"
              f" mu0 {d['mu'][0]:.4f} (ref {d['mu_ref'][0]:.4f})", flush=True)
    json.dump({k: {kk: np.asarray(vv).tolist() for kk, vv in v.items()} for k, v in res.items()},
              open(f"tpi_runs/{tag}_{econ}_exact_irf.json", "w"), indent=1)
