# RESUME THE SHOCKED HA PATH FROM A SAVED ITERATE WITH A LARGER FINITE-DIFFERENCE STEP; REPORT THE HOUSEHOLDS' OWN C.
#   WT=dff4 python3 ha_ge2.py <tag> <econ> <h> [T] [workers] [maxit]
import sys
ARGV = list(sys.argv)
sys.argv = ["x"]
import ha_path as hp
if __name__ == "__main__":
    import json
    import numpy as np
    from global_projection.solver_recursive.recursive_experiment import s_from_pd
    tag, econ, h = ARGV[1], ARGV[2], float(ARGV[3])
    T = int(ARGV[4]) if len(ARGV) > 4 else 60
    W = int(ARGV[5]) if len(ARGV) > 5 else 3
    mx = int(ARGV[6]) if len(ARGV) > 6 else 60
    hp.FD_H = h
    hp.context(tag); hp._setup_rest(econ); c = hp._CTX
    s1 = s_from_pd(0.02)
    z = np.load(f"tpi_runs/{tag}_{econ}_ha_T{T}_u.npz")
    u_ref = z["u_ref"]
    u, sho, e = hp.solve_path(tag, econ, s1, T, W, u0=z["u_sho"], maxit=mx)
    ref = hp.simulate(econ, u_ref[:T], u_ref[T:], None)
    f, _ = hp.residual(econ, u, s1, sho)
    print("  |F_HA| by decade:", " ".join(f"{np.max(np.abs(f[T:][k:k+10])):.0e}" for k in range(0, T, 10)),
          "| |F_UIP| by decade:", " ".join(f"{np.max(np.abs(f[:T][k:k+10])):.0e}" for k in range(0, T, 10)), flush=True)
    hA1, hC1 = hp.households(sho, "D", u[T:]); hA0, hC0 = hp.households(ref, "D", u_ref[T:])
    ha = hp.irf(hp.table(sho), hp.table(ref))
    ha["C_HA"] = 100 * (hC1 / hC0 - 1)
    ha["A_HA"] = 100 * (hA1 / hA0 - 1)
    ra = hp.irf(hp.table(hp.ra_path(tag, econ, s1, T)), hp.table(hp.ra_path(tag, econ, None, T)))
    print(f"   HA households' own C: " + " ".join(f"{x:+.3f}" for x in ha["C_HA"][:8]) + f" | cum25 {ha['C_HA'][:25].sum():+.2f}")
    print(f"   map C at the iterate: " + " ".join(f"{x:+.3f}" for x in ha["C_D"][:8]) + f" | cum25 {ha['C_D'][:25].sum():+.2f}")
    print(f"   RA C:                 " + " ".join(f"{x:+.3f}" for x in ra["C_D"][:8]) + f" | cum25 {ra['C_D'][:25].sum():+.2f}")
    print(f"   rdep_D HA-GE vs RA (bp/yr): " + " ".join(f"{a - b:+.2f}" for a, b in zip(ha["rdep_D"][:25:4], ra["rdep_D"][:25:4])))
    out = {econ: {"err": [float(e)], "h": h, "ha": {k: np.asarray(v).tolist() for k, v in ha.items()},
                  "ra": {k: np.asarray(v).tolist() for k, v in ra.items()}, "resid": f.tolist()}}
    json.dump(out, open(f"tpi_runs/{tag}_ha_T{T}_h{h:g}.json", "w"), indent=1)
    np.savez(f"tpi_runs/{tag}_{econ}_ha_T{T}_h{h:g}_u.npz", u_ref=u_ref, u_sho=u)
