# THE COST OF ONE FULL COLLOCATION RESIDUAL F (150 POINT MAPS) AND OF ONE COMPLEX-STEP JACOBIAN COLUMN, SERIAL.
import pickle
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))
from rung_anatomy import economy
from global_projection.solver_recursive.collocation import make_residual, pack, _column

ROOT = HERE.parents[2]
cal, ss, sproc = economy()
sv = pickle.load(open(ROOT / "results" / "GLOBAL_tpi_walk" / "tpi_rung_00060260.pkl", "rb"))
rules = sv["rules"].copy()
cal["tpi_cap_bp"] = sv["cap"]
F = make_residual(rules, cal, ss, sproc, regimes=(0, 1), n_gh=rules.n_gh or 5)
x = pack(rules, (0, 1))
t = time.perf_counter(); f = F(x); t1 = time.perf_counter() - t
t = time.perf_counter(); c = _column(F, x, f, 1e-6, 17, "complex"); t2 = time.perf_counter() - t
print(f"unknowns {x.size}, residuals {f.size}, max|F| {np.max(np.abs(f)):.1e}")
print(f"one F (150 point maps, real): {t1:.2f} s;  one complex-step column: {t2:.2f} s")
print(f"=> one full Jacobian = {x.size} columns ~ {x.size * t2 / 60:.0f} core-min, ~{x.size * t2 / 60 / 8:.1f} min on 8 cores")
nz = np.abs(c) > 1e-14
print(f"one column touches {nz.sum()} of {f.size} residuals (the perturbed rule's refit reaches every point's continuation)")
