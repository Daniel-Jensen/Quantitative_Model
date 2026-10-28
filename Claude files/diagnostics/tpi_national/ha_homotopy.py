# HA CLEARING ALONG THE HEADLINE-SHOCK PATH: REFERENCE AND SHOCKED PATHS, RESIDUAL BY DECADE, PATHS SAVED.
#   WT=dff4 python3 ha_homotopy.py <tag> <econ> [T] [workers]
import sys
ARGV = list(sys.argv)
sys.argv = ["x"]
import ha_path as hp
if __name__ == "__main__":
    import json
    import numpy as np
    from global_projection.solver_recursive.recursive_experiment import s_from_pd
    from global_projection.solver_recursive.state_grid import IS
    tag, econ = ARGV[1], ARGV[2]
    T = int(ARGV[3]) if len(ARGV) > 3 else 60
    W = int(ARGV[4]) if len(ARGV) > 4 else 3
    hp.context(tag); hp._setup_rest(econ); c = hp._CTX
    s0, s1 = float(c["S0"][econ][IS]), s_from_pd(0.02)
    u_ref, ref, e0 = hp.solve_path(tag, econ, None, T, W)
    # the shocked path from its own representative-household path, with a fresh Jacobian
    u, sho, e = hp.solve_path(tag, econ, s1, T, W)
    errs = [e]
    f, _ = hp.residual(econ, u, s1, sho)
    print("  |F_HA| by decade:", " ".join(f"{np.max(np.abs(f[T:][k:k+10])):.0e}" for k in range(0, T, 10)),
          "| |F_UIP| by decade:", " ".join(f"{np.max(np.abs(f[:T][k:k+10])):.0e}" for k in range(0, T, 10)), flush=True)
    ha = hp.irf(hp.table(sho), hp.table(ref))
    ra = hp.irf(hp.table(hp.ra_path(tag, econ, s1, T)), hp.table(hp.ra_path(tag, econ, None, T)))
    for lab, d in (("HA", ha), ("RA", ra)):
        print(f"   {lab}: Y " + " ".join(f"{x:+.3f}" for x in d["Y_D"][:8])
              + f" | cumY25 {d['Y_D'][:25].sum():+.2f}  C0 {d['C_D'][0]:+.3f}  cumC25 {d['C_D'][:25].sum():+.2f}"
              f"  I0 {d['I_D'][0]:+.3f}  n0 {d['n_D'][0]:+.2f}  rdep0 {d['rdep_D'][0]:+.1f}bp  sov0 {d['sov_bp_level'][0]:.1f}bp", flush=True)
    out = {econ: {"err": [e0] + errs, "map_resid": max(max(o["_resid"] for o in p["outs"]) for p in (ref, sho)),
                  "ha": {k: v.tolist() for k, v in ha.items()}, "ra": {k: v.tolist() for k, v in ra.items()},
                  "resid": f.tolist()}}
    json.dump(out, open(f"tpi_runs/{tag}_ha_T{T}.json", "w"), indent=1)
    np.savez(f"tpi_runs/{tag}_{econ}_ha_T{T}_u.npz", u_ref=u_ref, u_sho=u)
