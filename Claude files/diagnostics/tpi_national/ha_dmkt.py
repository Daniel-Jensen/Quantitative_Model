# D'S DEPOSIT MARKET CLEARED THROUGH THE HETEROGENEOUS HOUSEHOLDS, AT THE GLOBAL SOLUTION'S OTHER PRICES.
#   WT=dff4 python3 ha_dmkt.py <tag> <econ off|on> [T]
# The global solution's paths (wages, hours, prices, transfers, the banks' and the Eurosystem's choices) are
# held fixed. D's deposit rate is the price that clears: by deposit-UIP, extra household deposits are lent
# abroad and lower D's rate by kappa_nfa per unit of the position relative to GDP,
#   r_D,t = r_D,t(global) - kappa * (A_HA,t - A_global,t) * P_t / Y_ss,
# and the households' deposits A_HA respond to that rate path. Newton on the rate path; each evaluation is one
# EGM backward pass and one distribution forward pass.
import sys
ARGV = list(sys.argv)
sys.argv = ["x"]
import ha_path as hp
if __name__ == "__main__":
    import json
    import numpy as np
    from global_projection.blocks.household import solve_backward_transition
    from global_projection.blocks.distribution import forward_paths, aggregate_assets
    from global_projection.solver_recursive.recursive_experiment import s_from_pd
    tag, econ = ARGV[1], ARGV[2]
    T = int(ARGV[3]) if len(ARGV) > 3 else 60
    hp.context(tag); hp._setup_rest(econ); c = hp._CTX
    cal, ss, rest = c["cal"], c["ss"], c["rest"]["D"]
    kappa, Yss = cal["kappa_nfa"], ss["ss_firm_D"]["Y_ss"]
    a_grid, Pi, e = ss["a_grid_D"], ss["Pi_D"], ss["e_D"]
    fr, chi = cal["frisch_D"], cal["chi_D"]

    def prep(path):
        # THE GLOBAL PATH'S PRICES AND THE HOUSEHOLDS' OTHER RESOURCES (AT THE GLOBAL DEPOSIT RATE).
        o = path["outs"]
        P = np.array([x["P_CES_D"] for x in o]); R = 1.0 + np.array([x["rdep_D"] for x in o])
        N = np.array([x["_x"][0] for x in o]); w = np.array([x["w_D"] for x in o])
        W = np.array([x["W_D"] for x in o]); inc = np.array([x["inc_D"] for x in o])
        A = np.array([float(x["A_D"]) for x in o])
        rc = np.empty(T); rc[0] = rest["R"] * rest["P"] / P[0] - 1.0; rc[1:] = R[:-1] * P[:-1] / P[1:] - 1.0
        A_lag = np.concatenate([[aggregate_assets(rest["D0"], a_grid)], A[:-1]])
        labour = (w / P) * N
        tau = W / P + inc - labour - (1.0 + rc) * A_lag
        return dict(P=P, R=R, A=A, rc0=rc[0], y=labour[:, None] * e[None, :] + tau[:, None],
                    vN=chi * N ** (1 + 1 / fr) / (1 + 1 / fr))

    def households(g, dr):
        # DEPOSITS AND CONSUMPTION WHEN D'S DEPOSIT RATE IS SHIFTED BY dr (OWN DEPOSITS EARN IT).
        R1 = g["R"] + dr
        rc = np.empty(T); rc[0] = g["rc0"]; rc[1:] = R1[:-1] * g["P"][:-1] / g["P"][1:] - 1.0
        r_path = np.concatenate([rc, [R1[-1] * g["P"][-1] / rest["P"] - 1.0]])
        cp, ap = solve_backward_transition(a_grid, Pi, r_path, g["y"], rest["c_term"], rest["beta"],
                                           cal["sigma_D"], cal["a_min_D"], vN_path=g["vN"])
        A, C, _ = forward_paths(rest["D0"], ap, cp, a_grid, Pi)
        return A, C, cp

    def clear(g):
        # NEWTON ON THE RATE PATH: G(dr) = dr + kappa (A_HA(dr) - A_global) P / Y_ss = 0.
        G = lambda dr: dr + kappa * (households(g, dr)[0] - g["A"]) * g["P"] / Yss
        dr = np.zeros(T)
        f = G(dr)
        h = 1e-6
        J = np.empty((T, T))
        for j in range(T):
            d = dr.copy(); d[j] += h
            J[:, j] = (G(d) - f) / h
        hist = [float(np.max(np.abs(f)))]
        for it in range(30):
            step = -np.linalg.solve(J, f)
            a = 1.0
            while a > 1e-3:
                fn = G(dr + a * step)
                if np.linalg.norm(fn) < np.linalg.norm(f):
                    break
                a /= 2
            dr, f = dr + a * step, fn
            hist.append(float(np.max(np.abs(f))))
            if hist[-1] < 1e-11:
                break
        A, C, cp = households(g, dr)
        return dict(dr=dr, A=A, C=C, c0=cp[0], err=hist[-1], hist=hist, cond=float(np.linalg.cond(J)))

    s1 = s_from_pd(0.02)
    paths = {k: hp.ra_path(tag, econ, s, T) for k, s in (("ref", None), ("shock", s1))}
    g = {k: prep(p) for k, p in paths.items()}
    pe = {k: households(g[k], np.zeros(T)) for k in g}
    ge = {k: clear(g[k]) for k in g}
    pct = lambda a, b: 100 * (a / b - 1)
    C_ra = pct(np.array([float(x["C_D"]) for x in paths["shock"]["outs"]]),
               np.array([float(x["C_D"]) for x in paths["ref"]["outs"]]))
    out = dict(C_ra=C_ra, C_pe=pct(pe["shock"][1], pe["ref"][1]), C_ge=pct(ge["shock"]["C"], ge["ref"]["C"]),
               A_ra=pct(g["shock"]["A"], g["ref"]["A"]), A_pe=pct(pe["shock"][0], pe["ref"][0]),
               A_ge=pct(ge["shock"]["A"], ge["ref"]["A"]),
               dr_ge_bp=4e4 * (ge["shock"]["dr"] - ge["ref"]["dr"]),
               dr_level_bp={k: (4e4 * ge[k]["dr"]).tolist() for k in ge},
               err={k: ge[k]["err"] for k in ge}, cond={k: ge[k]["cond"] for k in ge})
    print(f"[{tag} {econ}] deposit-market clearing: |G| ref {ge['ref']['err']:.1e} shock {ge['shock']['err']:.1e}; "
          f"cond(J) {ge['shock']['cond']:.1f}; Newton history {['%.0e' % x for x in ge['shock']['hist']]}")
    for q in (0, 1, 2, 4, 8, 12, 24):
        print(f"  q{q:2d}: C  RA {out['C_ra'][q]:+7.3f}%  HA-PE {out['C_pe'][q]:+7.3f}%  HA-cleared {out['C_ge'][q]:+7.3f}%"
              f"   deposits RA {out['A_ra'][q]:+6.3f}% PE {out['A_pe'][q]:+6.3f}% cleared {out['A_ge'][q]:+6.3f}%"
              f"   clearing rate shift {out['dr_ge_bp'][q]:+6.2f} bp/yr")
    print(f"  25q sums: C RA {out['C_ra'][:25].sum():+.2f}  HA-PE {out['C_pe'][:25].sum():+.2f}  HA-cleared {out['C_ge'][:25].sum():+.2f}")
    # who cuts consumption on impact, cleared
    D0 = rest["D0"]
    if D0.shape != ge["ref"]["c0"].shape:
        D0 = D0.T
    ma = D0.sum(axis=1); cum = np.cumsum(ma) - 0.5 * ma; qi = np.minimum((cum * 5).astype(int), 4)
    rows = []
    for k in range(5):
        m = qi == k
        base = float((D0[m] * ge["ref"]["c0"][m]).sum())
        d = float((D0[m] * (ge["shock"]["c0"][m] - ge["ref"]["c0"][m])).sum())
        rows.append(dict(quintile=k + 1, dC_pct=100 * d / base))
    print("  impact C by wealth quintile (cleared): " + "  ".join(f"Q{r['quintile']} {r['dC_pct']:+.3f}%" for r in rows))
    out["quintiles_ge"] = rows
    json.dump({k: (np.asarray(v).tolist() if isinstance(v, np.ndarray) else v) for k, v in out.items()},
              open(f"tpi_runs/{tag}_{econ}_ha_dmkt.json", "w"), indent=1)
