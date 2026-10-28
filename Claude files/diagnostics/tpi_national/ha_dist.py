# WHO CUTS CONSUMPTION: THE IMPACT RESPONSE OF D CONSUMPTION BY WEALTH QUINTILE (HETEROGENEOUS HOUSEHOLDS).
#   WT=dff4 python3 ha_dist.py <tag> <econ> <pe|ge> [T]
# pe: the households face the representative-household paths; ge: the cleared paths saved by ha_homotopy.py.
import sys
ARGV = list(sys.argv)
sys.argv = ["x"]
import ha_path as hp
if __name__ == "__main__":
    import json
    import numpy as np
    import global_projection.blocks.household as hh
    from global_projection.solver_recursive.recursive_experiment import s_from_pd
    tag, econ, mode = ARGV[1], ARGV[2], ARGV[3]
    T = int(ARGV[4]) if len(ARGV) > 4 else 60
    hp.context(tag); hp._setup_rest(econ); c = hp._CTX
    s1 = s_from_pd(0.02)
    # capture the consumption policies households() computes
    grab, orig = {}, hh.solve_backward_transition
    def wrap(*a, **k):
        cp, ap = orig(*a, **k)
        grab["c"] = cp
        return cp, ap
    hh.solve_backward_transition = wrap
    pol = {}
    if mode == "pe":
        for lab, s in (("ref", None), ("shock", s1)):
            p = hp.ra_path(tag, econ, s, T)
            hp.households(p, "D", np.array([float(o["A_D"]) for o in p["outs"]]))
            pol[lab] = grab["c"][0]
    else:
        z = np.load(f"tpi_runs/{tag}_{econ}_ha_T{T}_u.npz")
        for lab, s, u in (("ref", None, z["u_ref"]), ("shock", s1, z["u_sho"])):
            p = hp.simulate(econ, u[:T], u[T:], s)
            hp.households(p, "D", u[T:])
            pol[lab] = grab["c"][0]
    D0 = c["rest"]["D"]["D0"]
    a_grid = c["ss"]["a_grid_D"]
    if D0.shape != pol["ref"].shape:
        D0 = D0.T
    # wealth quintiles of the rest-point distribution, by asset grid point
    ma = D0.sum(axis=1)
    cum = np.cumsum(ma) - 0.5 * ma
    q = np.minimum((cum * 5).astype(int), 4)
    rows = []
    for k in range(5):
        m = q == k
        base = float((D0[m] * pol["ref"][m]).sum())
        d = float((D0[m] * (pol["shock"][m] - pol["ref"][m])).sum())
        rows.append(dict(quintile=k + 1, mass=float(ma[m].sum()), a_lo=float(a_grid[m].min()), a_hi=float(a_grid[m].max()),
                         dC_pct=100 * d / base, share_of_total_dC=d))
    tot = sum(r["share_of_total_dC"] for r in rows)
    base_all = float((D0 * pol["ref"]).sum())
    at_limit = float(ma[0])
    print(f"[{tag} {econ} {mode}] impact dC_D {100 * tot / base_all:+.3f}%; mass at the borrowing limit {100 * at_limit:.1f}%")
    for r in rows:
        r["share_of_total_dC"] = r["share_of_total_dC"] / tot if tot else float("nan")
        print(f"  Q{r['quintile']}: assets {r['a_lo']:.2f}-{r['a_hi']:.2f}  mass {100 * r['mass']:.1f}%  dC {r['dC_pct']:+.3f}%  "
              f"share of the aggregate change {100 * r['share_of_total_dC']:.0f}%")
    json.dump(dict(rows=rows, total_pct=100 * tot / base_all, at_limit=at_limit),
              open(f"tpi_runs/{tag}_{econ}_ha_dist_{mode}.json", "w"), indent=1)
