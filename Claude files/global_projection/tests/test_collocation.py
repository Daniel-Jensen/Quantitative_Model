# GLOBAL COLLOCATION SOLVER: PACKING, THE IDENTITY RESIDUALS, AND A REAL SOLVE.
# The solver is the one thing between the period map and every reported number, so it
# needs its own gate. Three checks, cheapest first.
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from common import get_ss                                             # noqa: E402
from global_projection.solver_recursive.state_grid import (build_state_box, s_process_params,  # noqa: E402
                                         IS)
from global_projection.solver_recursive.decision_rules import RuleSet, STORE_RULES, SOLVE  # noqa: E402
from global_projection.solver_recursive.recursive_main import (calibrate_household_anchors,  # noqa: E402
                                             ss_state, ss_x, time_iteration)
from global_projection.solver_recursive.point_map import point_residuals                # noqa: E402
from global_projection.solver_recursive import collocation as C                         # noqa: E402

BOX = dict(k_band=0.02, p_band_D=0.04, p_band_F=0.04, b_band=0.12, w_band=0.04)


def _setup(refine=None, n_regimes=2):
    cal, ss = get_ss()
    cal["nw_floor_frac"] = 0.15
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    grid = build_state_box(ss, cal, mu=1, refine=refine, **BOX)
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=n_regimes)
    rules.n_gh = 5
    return cal, ss, sproc, rules


def test_pack_roundtrip():
    # THE FLAT UNKNOWN VECTOR MUST BE A LOSSLESS VIEW OF THE RULE VALUES.
    # The stacking convention is the one place a silent reordering would corrupt every
    # residual without raising, so it is checked rather than trusted.
    _, _, _, rules = _setup()
    regimes = tuple(range(rules.n_regimes))
    theta = C.pack(rules, regimes)
    assert theta.size == len(STORE_RULES) * len(regimes) * rules.grid.n
    back = C.unpack(theta, rules, regimes)
    for k in STORE_RULES:
        for d in regimes:
            assert np.allclose(back[k][d], rules.vals[k][d], rtol=0, atol=1e-12), k
    print("  pack/unpack round-trip exact over "
          f"{theta.size} unknowns: PASSED")


def test_identity_residuals_are_the_readoffs():
    # THE SIX ADDED RESIDUALS MUST BE log(guess / point_map's implied value).
    # This is Bocola's residual 2 (log(alp/alp_impl)) generalised to every object the
    # old code READ OFF a frozen continuation. If the wiring were wrong the system
    # would still "converge" -- to the wrong fixed point.
    cal, ss, sproc, rules = _setup()
    F = C.make_residual(rules, cal, ss, sproc, regimes=(0,), no_default=True, n_gh=5)
    r = F(C.pack(rules, (0,))).reshape(rules.grid.n, C.N_RES)
    from global_projection.solver_recursive.decision_rules import DERIVED, to_fit
    for i in (0, rules.grid.n // 2, rules.grid.n - 1):
        S = rules.grid.points[i]
        x = np.array([rules.vals[k][0][i] for k in SOLVE])
        res, out = point_residuals(S, 0, x, rules, cal, ss, sproc, n_gh=5,
                                   no_default=True)
        assert np.allclose(r[i, :C.N_RES_POINT], res, rtol=0, atol=1e-12)
        for q, k in enumerate(DERIVED):
            want = to_fit(k, rules.vals[k][0][i]) - to_fit(k, out[k])
            assert abs(r[i, C.N_RES_POINT + q] - want) < 1e-12, k
    print(f"  {C.N_RES} residuals per point = {C.N_RES_POINT} from point_map "
          f"+ {len(DERIVED)} identities: PASSED")


def test_solve_reaches_the_floor():
    # A REAL SOLVE MUST ROOT THE WHOLE SYSTEM, not merely settle a damped iteration.
    # d = 0 at pi = 0, which is the stage whose answer is pinned independently: the
    # steady state must come back out of it exactly.
    cal, ss, sproc, rules = _setup()
    time_iteration(rules, cal, ss, sproc, regimes=(0,), no_default=True, damp=0.5,
                   tol=1e-4, max_it=12, n_gh=5, verbose=False)
    ok, its, worst = C.solve_collocation(rules, cal, ss, sproc, regimes=(0,),
                                         no_default=True, n_gh=5, backend="parsolve",
                                         maxit=20, verbose=False, label=" test")
    assert ok, f"collocation solve did not converge (max|F| = {worst:.2e})"
    assert worst <= 10 * C.TOL_MAXF, f"max|F| = {worst:.2e}"
    # the SS is a collocation node, so the solved rules must reproduce it there
    S0 = ss_state(ss, cal, sproc)
    x = np.array([float(rules.eval(k, 0, np.atleast_2d(S0))[0]) for k in SOLVE])
    assert np.allclose(x, ss_x(ss, cal), rtol=2e-6, atol=2e-6), \
        f"SS policy not reproduced: max dev {np.max(np.abs(x - ss_x(ss, cal))):.2e}"
    print(f"  d=0 collocation solve: converged in {its} Newton steps, "
          f"max|F| = {worst:.1e}, SS reproduced: PASSED")


def test_refined_grid_is_square_and_exact():
    # THE TENSOR-REFINED GRID MUST STILL BE AN EXACT INTERPOLATION OPERATOR.
    _, _, _, rules = _setup(refine=(IS, 5))
    g = rules.grid
    assert g._Phi.shape == (g.n, g.n), g._Phi.shape
    # sparse(mu=1) over the other 11 states = 23 points, times the 5 dense s nodes
    assert g.n == 5 * (2 * (g.d - 1) + 1) == 115, (g.n, g.d)
    rng = np.random.default_rng(0)
    y = rng.normal(size=g.n)
    assert np.max(np.abs(g.eval(g.fit(y), g.points) - y)) < 1e-10
    assert np.linalg.cond(g._Phi) < 500
    print(f"  refined grid {g.n} points, degree {g.max_deg}, "
          f"cond {np.linalg.cond(g._Phi):.0f}, on-node exact: PASSED")


if __name__ == "__main__":
    test_pack_roundtrip()
    test_identity_residuals_are_the_readoffs()
    test_refined_grid_is_square_and_exact()
    test_solve_reaches_the_floor()
    print("test_collocation: ALL PASSED")
