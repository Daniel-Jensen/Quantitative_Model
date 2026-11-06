# THE TWO STOCHASTIC REST POINTS (NO TPI, TPI) STATE BY STATE, AND THE SPREAD AT EACH UNDER BOTH RULE SETS: IS THE TPI's
# REST-POINT GAP A RULE EFFECT OR A DIFFERENT REST STATE?
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
from rung_anatomy import economy
from global_projection.solver_recursive.recursive_experiment import stochastic_rest_point, read_exact, _sov_bp
from global_projection.solver_recursive.state_grid import STATE_NAMES

off = pickle.load(open(ROOT / "results" / "GLOBAL" / "data" / "rules_risk.pkl", "rb"))
on = pickle.load(open(ROOT / "results" / "GLOBAL" / "data" / "rules_tpi.pkl", "rb"))
cal, ss, sproc = economy()
cal["tpi_cap_bp"] = 200.0
for k in ("mu_eps", "book_eps", "floor_eps"):
    cal[k] = 0.0
c_off = dict(cal, tpi_on=False)
c_on = dict(cal, tpi_on=True, tpi_eps=getattr(on, "tpi_eps_solved", 1e-4), tpi_gz=True)
S0 = stochastic_rest_point(off, c_off, ss, sproc, verbose=False)
S1 = stochastic_rest_point(on, c_on, ss, sproc, verbose=False)
print("state      no-TPI rest    TPI rest      diff %")
for j, n in enumerate(STATE_NAMES):
    d = 100 * (S1[j] / S0[j] - 1) if abs(S0[j]) > 1e-9 else np.nan
    print(f"{n:8s} {S0[j]:12.5f} {S1[j]:12.5f}   {d:+8.3f}")
for lab, S in (("no-TPI rest", S0), ("TPI rest", S1)):
    a, b = read_exact(off, c_off, ss, sproc, S.copy()), read_exact(on, c_on, ss, sproc, S.copy())
    print(f"at the {lab}: spread no-TPI rules {_sov_bp(a, cal):6.1f} bp, TPI rules {_sov_bp(b, cal):6.1f} bp;"
          f" rdep_D {4e4 * a['rdep_D']:6.1f} / {4e4 * b['rdep_D']:6.1f} bp; n_D {a['n_D']:.4f} / {b['n_D']:.4f}")
