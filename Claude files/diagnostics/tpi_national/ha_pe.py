# PARTIAL EQUILIBRIUM: HETEROGENEOUS D HOUSEHOLDS FACING THE GLOBAL EQUILIBRIUM'S PRICE AND INCOME PATHS.
#   WT=dff4 python3 ha_pe.py <tag> <econ off|on> [T]
import sys, json
from pathlib import Path
tag, econ = sys.argv[1], sys.argv[2]
T = int(sys.argv[3]) if len(sys.argv) > 3 else 60
sys.argv = ["x"]
import ha_path as hp
if __name__ == "__main__":
    import numpy as np
    from global_projection.solver_recursive.recursive_experiment import s_from_pd
    hp.context(tag); hp._setup_rest(econ); c = hp._CTX
    out = {}
    for lab, s in (("ref", None), ("shock", s_from_pd(0.02))):
        ra = hp.ra_path(tag, econ, s, T)
        A_ra = np.array([float(o["A_D"]) for o in ra["outs"]]); C_ra = np.array([float(o["C_D"]) for o in ra["outs"]])
        hA, hC = hp.households(ra, "D", A_ra)
        out[lab] = dict(A_ra=A_ra, C_ra=C_ra, A_ha=hA, C_ha=hC, P=np.array([o["P_CES_D"] for o in ra["outs"]]),
                        Y=np.array([o["Y_D"] for o in ra["outs"]]))
    r, s_ = out["ref"], out["shock"]
    d = lambda k: 100 * (s_[k] / r[k] - 1)
    ra_C, ha_C, ra_A, ha_A = d("C_ra"), d("C_ha"), d("A_ra"), d("A_ha")
    # the cross-border position the HA saving gap implies, and the deposit-rate feedback through UIP
    gap = (s_["A_ha"] - s_["A_ra"]) - (r["A_ha"] - r["A_ra"])
    dV = gap * s_["P"]
    kappa = c["cal"]["kappa_nfa"]
    rate_fb = 4e4 * kappa * dV / c["ss"]["ss_firm_D"]["Y_ss"]
    print(f"[{tag} {econ}] T={T}  rest beta {c['rest']['D']['beta']:.6f}")
    for q in (0, 1, 2, 4, 8, 12, 20, 24):
        print(f"  q{q:2d}: C  RA {ra_C[q]:+7.3f}%  HA {ha_C[q]:+7.3f}%   deposits RA {ra_A[q]:+6.3f}%  HA {ha_A[q]:+6.3f}%"
              f"   HA-RA deposit gap {gap[q]:+.4f} (={100*dV[q]:+.2f}% of qtr GDP), UIP rate feedback {rate_fb[q]:+.2f} bp/yr")
    print(f"  25q cumulative: C RA {ra_C[:25].sum():+.2f}  HA {ha_C[:25].sum():+.2f};  max |rate feedback| {np.max(np.abs(rate_fb[:25])):.2f} bp/yr")
    json.dump({k: v.tolist() for k, v in dict(ra_C=ra_C, ha_C=ha_C, ra_A=ra_A, ha_A=ha_A, gap=gap, rate_fb=rate_fb).items()},
              open(f"tpi_runs/{tag}_{econ}_ha_pe.json", "w"), indent=1)
