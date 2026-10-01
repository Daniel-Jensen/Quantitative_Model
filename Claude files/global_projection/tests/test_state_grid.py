# STATE-GRID GATE: SMOLYAK COUNTS, ON-GRID EXACTNESS, QUADRATIC EXACTNESS,
# THE STATE BOX, AND THE s-PROCESS EXPERIMENT MAPPING.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from solver_recursive.state_grid import (SmolyakGrid, build_state_box, default_prob,
                        s_process_params, STATE_NAMES, IS)
from solver_recursive.recursive_main import ss_state
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
    # THE STATE BOX CONTAINS THE SS POINT; s-PROCESS IS BOCOLA'S.
    # The SS point comes from ss_state, NOT a literal rebuilt here. The old hand-built
    # 7-element array was a duplicate of the convention, so when the state vector grew
    # to 9 it kept asserting against the wrong thing instead of failing -- the same
    # failure mode as the duplicated 7-entry x_ss literal in recursive_main._sweep.
    cal, ss = get_ss()
    g = build_state_box(ss, cal, mu=1)
    sp = s_process_params(cal)
    ss_pt = ss_state(ss, cal, sp)
    assert g.d == len(STATE_NAMES) == ss_pt.size, (
        f"grid dim {g.d}, STATE_NAMES {len(STATE_NAMES)}, ss_state {ss_pt.size}")
    assert np.all(ss_pt > g.lo) and np.all(ss_pt < g.hi), "SS not interior"

    assert abs(default_prob(sp["s_star"]) - 0.001) < 1e-10   # rest p^d = 0.1%
    assert sp["rho_s"] == 0.95                                # Bocola param(22)
    # 0.63/sqrt(2): Bocola's reported param(23) is 0.63, but his GaussHermite.m is the
    # PHYSICISTS' rule under a sigma*node map, so his solved model's innovation sd is
    # 0.63/sqrt(2). gh_nodes here is the probabilists' rule, so the EFFECTIVE sd is what
    # has to match. See calibration.py.
    assert abs(sp["sigma_s"] - 0.63 / np.sqrt(2.0)) < 1e-4


def test_s_box_covers_the_process():
    # THE s BOX MUST BE WIDE RELATIVE TO THE INNOVATION, NOT JUST TO s* .
    # The failure this guards against is silent: a box narrow in CONDITIONAL sd makes
    # grid.clip truncate the quadrature, which mean-reverts s far harder than rho_s and
    # kills the long bond's repricing. Bocola's own box is +-9.8 conditional sd.
    cal, ss = get_ss()
    sp = s_process_params(cal)
    g = build_state_box(ss, cal, mu=1)
    lo, hi = g.lo[IS], g.hi[IS]
    assert lo < sp["s_star"] < hi
    for half in ((sp["s_star"] - lo), (hi - sp["s_star"])):
        assert half / sp["sigma_s"] > 5.0, (
            f"s box only {half / sp['sigma_s']:.2f} conditional sd wide")
    unc = sp["sigma_s"] / np.sqrt(1 - sp["rho_s"] ** 2)
    assert (hi - sp["s_star"]) / unc > 2.0, "s box under 2 unconditional sd"

    # and the clip must not distort the process in the ergodic core
    x, w = np.polynomial.hermite_e.hermegauss(5)
    w = w / w.sum()
    grid_s = np.linspace(sp["s_star"] - unc, sp["s_star"] + unc, 21)
    Ec = np.array([w @ np.clip((1 - sp["rho_s"]) * sp["s_star"] + sp["rho_s"] * s
                               + sp["sigma_s"] * x, lo, hi) for s in grid_s])
    rho_eff = np.polyfit(grid_s, Ec, 1)[0]
    assert rho_eff > 0.93, f"box clip cuts effective rho_s to {rho_eff:.3f}"


def test_rotated_box_round_trip_and_identity():
    # THE ROTATED BOX IS A CHANGE OF COORDINATES ONLY: identity nests exactly, the
    # transform round-trips, collocation stays exact, and clip lands inside.
    rng = np.random.default_rng(3)
    lo, hi = -np.ones(4), np.ones(4)
    g0 = SmolyakGrid(lo, hi, mu=2)
    assert g0.rot is None and np.allclose(g0.points, g0.from_unit(g0.points_unit))

    R = np.eye(4)
    R[np.ix_((1, 2), (1, 2))] = np.array([[1.0, -0.35], [2.0, -1.2]])
    c = np.array([0.3, -0.2, 0.1, 0.0])
    g = SmolyakGrid(lo, hi, mu=2, rot=R, centre=c)
    assert g.n == g0.n
    # round trip natural <-> unit
    u = rng.uniform(-1, 1, size=(200, 4))
    assert np.max(np.abs(g.to_unit(g.from_unit(u)) - u)) < 1e-10
    # collocation still exact at the (natural-coordinate) points
    vals = rng.standard_normal(g.n)
    assert np.max(np.abs(g.eval(g.fit(vals), g.points) - vals)) < 1e-9
    # clip projects onto the box IN ROTATED COORDS and reports zero violation after
    far = g.from_unit(rng.uniform(-3, 3, size=(200, 4)))
    assert np.max(g.outside(g.clip(far))) < 1e-10
    # a rotation with rot=I and centre=0 must reproduce the axis-aligned grid exactly
    gi = SmolyakGrid(lo, hi, mu=2, rot=np.eye(4), centre=np.zeros(4))
    assert np.max(np.abs(gi.points - g0.points)) < 1e-12


if __name__ == "__main__":
    test_counts_and_exactness()
    test_quadratic_exactness()
    test_state_box_and_s_process()
    test_s_box_covers_the_process()
    test_rotated_box_round_trip_and_identity()
    print("test_state_grid: ALL PASSED")
