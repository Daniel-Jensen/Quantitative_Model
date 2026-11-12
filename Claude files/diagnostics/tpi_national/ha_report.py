# RA, HA PARTIAL EQUILIBRIUM AND HA GENERAL EQUILIBRIUM, SIDE BY SIDE, FOR ONE RUN AND ECONOMY.
#   python3 ha_report.py <tag> <econ> [T]
import json, sys
import numpy as np
tag, econ = sys.argv[1], sys.argv[2]
T = int(sys.argv[3]) if len(sys.argv) > 3 else 60
ge = json.load(open(f"tpi_runs/{tag}_ha_T{T}.json"))[econ]
pe = json.load(open(f"tpi_runs/{tag}_{econ}_ha_pe.json"))
ha, ra = ({k: np.array(v) for k, v in ge[w].items()} for w in ("ha", "ra"))
pC = np.array(pe["ha_C"])
print(f"[{tag} {econ}] solve error (relative deposits) ref {ge['err'][0]:.1e}, shock {ge['err'][1]:.1e}; map clearing {ge['map_resid']:.1e}")
qs = (0, 1, 2, 4, 8, 12, 24)
rows = [("C_D (%)", ra["C_D"], pC, ha["C_D"]), ("Y_D (%)", ra["Y_D"], None, ha["Y_D"]),
        ("deposits A_D (%)", ra["A_D"], np.array(pe["ha_A"]), ha["A_D"]),
        ("rdep_D (bp/yr)", ra["rdep_D"], None, ha["rdep_D"]), ("lending spread (bp/yr)", ra["spread_bp"], None, ha["spread_bp"]),
        ("sov spread (bp/yr)", ra["sov_bp"], None, ha["sov_bp"]), ("n_D (%)", ra["n_D"], None, ha["n_D"]),
        ("I_D (%)", ra["I_D"], None, ha["I_D"]), ("C_F (%)", ra["C_F"], None, ha["C_F"])]
print("  variable                 model    " + "  ".join(f"q{q:<5d}" for q in qs) + "  sum q0-24")
for name, r, p, h in rows:
    for lab, x in (("RA", r), ("HA-PE", p), ("HA-GE", h)):
        if x is None:
            continue
        print(f"  {name:24s} {lab:6s} " + "  ".join(f"{x[q]:+6.3f}" for q in qs) + f"  {x[:25].sum():+7.2f}")
