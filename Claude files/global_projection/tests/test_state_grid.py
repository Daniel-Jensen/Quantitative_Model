# STATE-GRID GATE: SMOLYAK COUNTS, ON-GRID EXACTNESS, QUADRATIC EXACTNESS,
# THE STATE BOX, THE s-PROCESS EXPERIMENT MAPPING, THE C1 BOX EDGE AND THE AXIS REFINEMENT.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import numpy as np

from global_projection.solver_recursive.state_grid import (SmolyakGrid, build_state_box, default_prob,
                        s_process_params, STATE_NAMES, IS, CLIP_EPS, build_aligned_box,
                        align_weights, _box_centre, chebyshev_basis_1d)
from global_projection.solver_recursive.recursive_main import ss_state
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
    # THE STATE BOX CONTAINS THE SS POINT; THE s-PROCESS IS THE CALIBRATION'S.
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

    assert abs(default_prob(sp["s_star"]) - cal["pd_rest"]) < 1e-10   # rest p^d = cal["pd_rest"]
    # appendix B (2026-10-26): rho_s 0.9, sigma_s 0.41, read from the calibration (rho_s was
    # hard-coded at Bocola's 0.95 inside s_process_params); gh_nodes is the probabilists' rule,
    # so sigma_s is the innovation sd as is
    assert sp["rho_s"] == cal["rho_s"] == 0.9
    assert sp["sigma_s"] == cal["sigma_s"] == 0.41


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
    assert rho_eff > sp["rho_s"] - 0.02, f"box clip cuts effective rho_s to {rho_eff:.3f}"


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


def test_c1_box_edge():
    # THE MODEL BOX READS OUTSIDE STATES THROUGH A C1 SATURATION: EXACT INSIDE, BOUNDED, NO CORNER.
    cal, ss = get_ss()
    g = build_state_box(ss, cal, mu=1)
    assert g.clip_eps == CLIP_EPS > 0.0
    rng = np.random.default_rng(3)
    inside = g.from_unit(rng.uniform(-1, 1, size=(200, g.d)))
    assert np.array_equal(g.clip(inside), inside), "a state inside the box must be read unchanged"
    far = g.from_unit(rng.uniform(-8, 8, size=(200, g.d)))
    assert np.max(np.abs(g.to_unit(g.clip(far)))) <= 1.0 + CLIP_EPS + 1e-9
    # one-sided slopes across the upper edge of every dimension agree (C1), and match 1 inside
    h = 1e-7
    for q in range(g.d):
        e = np.zeros(g.d); e[q] = 1.0
        up = lambda t: g.to_unit(g.clip(g.from_unit(t * e)))[0, q]
        s_in, s_out = (up(1.0) - up(1.0 - h)) / h, (up(1.0 + h) - up(1.0)) / h
        assert abs(s_in - 1.0) < 1e-5 and abs(s_out - 1.0) < 1e-5, (STATE_NAMES[q], s_in, s_out)
    # a raw grid keeps the hard clip
    raw = SmolyakGrid(-np.ones(2), np.ones(2), mu=1)
    assert np.array_equal(raw.clip(np.array([[3.0, -0.5]])), np.array([[1.0, -0.5]]))


def test_axis_refine():
    # THE AXIS REFINEMENT: THE COARSE GRID PLUS m - 3 NODES ON THE s AXIS, SQUARE, NESTED, EXACT ON ITS SPAN.
    d = 12
    g0 = SmolyakGrid(-np.ones(d), np.ones(d), mu=1)
    X = np.random.default_rng(5).uniform(-1, 1, size=(300, d))
    for m in (3, 5, 9):
        g = SmolyakGrid(-np.ones(d), np.ones(d), mu=1, refine=(IS, m, "axis"))
        assert g.n == g0.n + m - 3 and g.max_deg == max(2, m - 1)
        assert all(np.any(np.all(np.isclose(g.points_unit, p), axis=1)) for p in g0.points_unit)
        # the new nodes lie on the s axis and add degree in s alone: no new s x state cross term
        new = ~np.array([np.any(np.all(np.isclose(g0.points_unit, p), axis=1)) for p in g.points_unit])
        assert np.all(np.delete(g.points_unit[new], IS, axis=1) == 0.0)
        assert np.all(np.delete(g.degrees[g0.n:], IS, axis=1) == 0)
        f = lambda U: (chebyshev_basis_1d(U[:, IS], m - 1)[:, m - 1] + 0.3 * U[:, 0] ** 2
                       + 0.5 * U[:, 4] - 0.2)
        err = np.max(np.abs(g.eval(g.fit(f(g.points_unit)), X) - f(X)))
        assert err < 1e-12, (m, err)
        assert np.linalg.cond(g._Phi) < 1e3
    try:
        SmolyakGrid(-np.ones(d), np.ones(d), mu=1, refine=(IS, 4, "axis"))
        raise AssertionError("4 nodes are not nested with the coarse s nodes")
    except AssertionError as e:
        assert "2^k + 1" in str(e)
    # the tensor default is unchanged
    assert SmolyakGrid(-np.ones(d), np.ones(d), mu=1, refine=(IS, 5)).n == 115
    print(f"  axis refinement: 25 + (m - 3) points, nested, exact on its span (m = 3, 5, 9)")


def test_cross_refine():
    # THE CROSS REFINEMENT: THE AXIS GRID PLUS FOUR NODES PER s x STATE PAIR, SQUARE, NESTED, EXACT ON
    # s x STATE PRODUCTS UP TO DEGREE 2 IN EACH, AND THE s BOX FIXED WHEREVER THE REST POINT IS.
    d = 12
    X = np.random.default_rng(7).uniform(-1, 1, size=(300, d))
    for m in (3, 5, 9):
        ga = SmolyakGrid(-np.ones(d), np.ones(d), mu=1, refine=(IS, m, "axis"))
        g = SmolyakGrid(-np.ones(d), np.ones(d), mu=1, refine=(IS, m, "cross"))
        assert g.n == ga.n + 4 * (d - 1)
        assert all(np.any(np.all(np.isclose(g.points_unit, p), axis=1)) for p in ga.points_unit)
        f = lambda U: (chebyshev_basis_1d(U[:, IS], m - 1)[:, m - 1] + 0.7 * (2 * U[:, IS] ** 2 - 1) * U[:, 4]
                       - 0.4 * U[:, IS] * (2 * U[:, 0] ** 2 - 1) + 0.3 * U[:, 2] ** 2 - 0.2)
        err = np.max(np.abs(g.eval(g.fit(f(g.points_unit)), X) - f(X)))
        assert err < 1e-12, (m, err)
        assert np.linalg.cond(g._Phi) < 1e3
        # the additive axis grid cannot carry the cross term
        assert np.max(np.abs(ga.eval(ga.fit(f(ga.points_unit)), X) - f(X))) > 1e-2
    cal, ss = get_ss()
    lo0 = build_state_box(ss, cal, mu=1).lo[IS]
    sp = s_process_params(dict(cal, pd_rest=0.0002))
    assert abs(sp["s_star"] - np.log(0.0002 / 0.9998)) < 1e-14
    assert build_state_box(ss, dict(cal, pd_rest=0.0002), mu=1).lo[IS] == lo0
    print(f"  cross refinement: axis + 4 (d - 1) points, exact on s x state products; s box fixed in pd_rest")


def test_box_at():
    # A MOVED CENTRE FOR THE BANK STOCKS SHIFTS ONLY THEM, KEEPING THIS SS's BAND WIDTHS; AT THE
    # SS's OWN CENTRE THE BOX IS UNCHANGED.
    from global_projection.solver_recursive.state_grid import BOX_STOCKS
    cal, ss = get_ss()
    c = _box_centre(ss, cal)
    g0 = build_state_box(ss, cal, mu=1)
    g1 = build_state_box(ss, cal, mu=1, at={k: c[k] for k in BOX_STOCKS})
    assert np.array_equal(g0.points, g1.points)
    shift = {k: 0.01 * abs(c[k]) + 1e-3 for k in BOX_STOCKS}
    g2 = build_state_box(ss, cal, mu=1, at={k: c[k] + v for k, v in shift.items()})
    moved = [STATE_NAMES.index(k) for k in BOX_STOCKS]
    d_lo, d_hi = g2.lo - g0.lo, g2.hi - g0.hi
    assert np.allclose(d_lo[moved], list(shift.values()), rtol=0, atol=1e-12)
    assert np.allclose(d_hi[moved], list(shift.values()), rtol=0, atol=1e-12)
    rest = [i for i in range(g0.d) if i not in moved]
    assert np.array_equal(g2.lo[rest], g0.lo[rest]) and np.array_equal(g2.hi[rest], g0.hi[rest])
    print("  box at: own centre unchanged; a moved centre shifts the bank stocks only")


def test_aligned_box_identity_axes():
    # THE ALIGNED BUILDER WITH IDENTITY AXES AND THE OLD V5 BANDS REBUILDS THE V5 BOX EXACTLY.
    cal, ss = get_ss()
    c, W = _box_centre(ss, cal), align_weights(ss, cal)
    half = np.array([0.02 * c["K_D"], 0.0226 * c["P_D"], 0.12 * c["b_DD"], 0.02 * c["P_D"]]) * W
    spec = dict(axes=np.eye(4).tolist(), lo=(-half).tolist(), hi=half.tolist(),
                bands=dict(K_F=[-0.02, 0.02], P_F=[-0.0226, 0.0226], b_DF=[-0.12, 0.12],
                           b_FD=[-0.12, 0.12], Z_D=[-0.03, 0.03], M_cb=0.32, O_cb=0.25))
    v5 = dict(k_band=0.02, p_band_D=0.0226, p_band_F=0.0226, b_band=0.12, w_band=0.02, m_band=0.32)
    for refine in (None, (IS, 5)):
        g0 = build_state_box(ss, cal, mu=1, refine=refine, **v5)
        g1 = build_aligned_box(ss, cal, spec, mu=1, refine=refine)
        assert g0.n == g1.n
        assert np.max(np.abs(g1.points - g0.points) / np.maximum(np.abs(g0.points), 1e-12)) < 1e-12
        x = g0.from_unit(np.random.default_rng(1).uniform(-3, 3, size=(50, g0.d)))
        assert np.allclose(g1.clip(x), g0.clip(x), rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    test_counts_and_exactness()
    test_quadratic_exactness()
    test_state_box_and_s_process()
    test_s_box_covers_the_process()
    test_rotated_box_round_trip_and_identity()
    test_c1_box_edge()
    test_axis_refine()
    test_cross_refine()
    test_box_at()
    test_aligned_box_identity_axes()
    print("test_state_grid: ALL PASSED")
