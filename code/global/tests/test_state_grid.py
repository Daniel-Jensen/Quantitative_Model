# STATE-GRID GATE: SMOLYAK COUNTS, ON-GRID EXACTNESS, QUADRATIC EXACTNESS,
# THE AGREED 8-STATE BOX, AND THE s-PROCESS EXPERIMENT MAPPING.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from solver_recursive.state_grid import (SmolyakGrid, build_state_box, default_prob,
                        s_process_params, STATE_NAMES)
from common import get_ss


def test_counts_and_exactness():
    # KNOWN POINT COUNTS (d=8 mu=2 -> 145) AND COLLOCATION EXACTNESS.
    g = SmolyakGrid(-np.ones(8), np.ones(8), mu=2)
    assert g.n == 145, f"expected 145 Smolyak points in 8-D, got {g.n}"
    rng = np.random.default_rng(0)
    vals = rng.standard_normal(g.n)
    assert np.max(np.abs(g.eval(g.fit(vals), g.points) - vals)) < 1e-9


def test_quadratic_exactness():
    # LEVEL 2 REPRODUCES COMPLETE QUADRATICS INCL. CROSS TERMS (KK04 PROPERTY).
    rng = np.random.default_rng(1)
    d = 8
    g = SmolyakGrid(-np.ones(d), np.ones(d), mu=2)
    A = rng.standard_normal((d, d)); A = 0.5 * (A + A.T)
    b = rng.standard_normal(d); c0 = rng.standard_normal()

    def f(x):
        # RANDOM COMPLETE QUADRATIC WITH ALL CROSS TERMS.
        return np.einsum("ij,jk,ik->i", x, A, x) + x @ b + c0

    coef = g.fit(f(g.points))
    x = rng.uniform(-1, 1, size=(800, d))
    assert np.max(np.abs(g.eval(coef, x) - f(x))) < 1e-9


def test_state_box_and_s_process():
    # THE AGREED BOX CONTAINS THE SS POINT; s-MAPPING MATCHES THE EXPERIMENT.
    cal, ss = get_ss()
    g = build_state_box(ss, cal)
    assert g.d == 8 and g.n == 145
    bk_D, bk_F = ss["ss_bank_D"], ss["ss_bank_F"]
    ss_pt = np.array([ss["Kap_D_ss"], ss["Kap_F_ss"],
                      (1 + cal["r_dep_D_target"]) * bk_D["Dep_supply_ss"],
                      (1 + cal["r_dep_F_target"]) * bk_F["Dep_supply_ss"],
                      cal["B_gov_D_ss"],
                      (1 + cal["r_dep_D_target"]) * ss["A_D_ss"],
                      (1 + cal["r_dep_F_target"]) * ss["A_F_ss"],
                      s_process_params(cal)["s_star"]])
    assert np.all(ss_pt > g.lo) and np.all(ss_pt < g.hi), "SS not interior"
    assert len(STATE_NAMES) == 8

    sp = s_process_params(cal)
    assert abs(default_prob(sp["s_star"]) - 0.001) < 1e-10   # rest p^d = 0.1%
    s_2sd = sp["s_star"] + 2.0 * sp["sigma_s"]
    assert abs(default_prob(s_2sd) - 0.02) < 1e-10           # +2sd -> p^d = 2%
    assert sp["rho_s"] == 0.95                                # experiment rho


if __name__ == "__main__":
    test_counts_and_exactness()
    test_quadratic_exactness()
    test_state_box_and_s_process()
    print("test_state_grid: ALL PASSED")
