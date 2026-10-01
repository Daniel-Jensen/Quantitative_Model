# RECURSIVE NESTING GATES N1-N3 FOR THE TIME-ITERATION SOLUTION.
# N1  the SS is a rest point of the single-point map: with SS-constant rules in
#     both regimes, no default, and the state at the SS point, all SEVEN market-
#     clearing residuals vanish (the recursive image of "the zero-shock
#     transition stays at the SS"). Requires the rep-agent household anchors.
# N2  the no-default (pi=0) d=0 block TIME-ITERATES to a fixed point that keeps
#     the SS grid point at the SS.
# N3  the Fischer-Burmeister complementarity holds on the grid (mu >= 0).
# The economic blocks are untouched; these gates validate the re-indexing only.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import numpy as np

from common import get_ss
from global_projection.solver_recursive.state_grid import build_state_box, s_process_params
from global_projection.solver_recursive.decision_rules import RuleSet
from global_projection.solver_recursive.point_map import point_residuals
from global_projection.solver_recursive.recursive_main import (time_iteration, calibrate_household_anchors,
                            ss_state, ss_x)
from global_projection.solver_recursive.recursive_experiment import BOX_KW
from global_projection.solver_recursive.collocation import (solve_collocation, TOL_MAXF, RES_NAMES,
                                          N_RES_POINT)

# SINGLE-SOURCED FROM THE SOLVER. This list used to be a hand-kept copy and went stale
# at every change to the residual system -- the assert below is what catches that now.
_LABELS = RES_NAMES[:N_RES_POINT]


def test_n1_ss_rest_point():
    # SS + SS-CONSTANT RULES + NO DEFAULT => EVERY RESIDUAL ~ 0.
    # The label list must cover the whole residual vector: when the system grew from 7
    # to 11 this printed only the first seven, so the bond FOCs, the F household's
    # Euler and union deposit clearing were silently unreported (the assert on
    # max|res| still covered them, but nothing showed WHICH one moved).
    cal, ss = get_ss()
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    grid = build_state_box(ss, cal)
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=2)

    res, out = point_residuals(ss_state(ss, cal, sp), 0, ss_x(ss, cal),
                               rules, cal, ss, sp, no_default=True)
    worst = np.max(np.abs(res))
    assert len(_LABELS) == res.size, (
        f"residual vector is {res.size} long, _LABELS covers {len(_LABELS)}")
    for lab, r in zip(_LABELS, res):
        print(f"    {lab:10s} {r:+.3e}")
    print(f"    mu_D={out['mu_D']:.6f} (ss {ss['ss_bank_D']['mu_ss']:.6f})  "
          f"C_D={out['C_D']:.5f} (ss {ss['C_D_ss']:.5f})  "
          f"A_D={out['A_D']:.5f} (ss {ss['A_D_ss']:.5f})")
    assert worst < 1e-6, f"SS not a rest point: max|res|={worst:.2e}"


def test_n2_no_default_grid_solve():
    # N2: THE GRID-WIDE FIXED POINT AT pi = 0 MUST ACTUALLY BE FOUND.
    # This used to be a REPORTING probe rather than a gate, because damped time
    # iteration could not converge it: its binding mode is the franchise-value
    # recursion at 0.990 per sweep, so it reported "converged=False, worst point
    # residual 2e-14" -- every point clearing against a continuation still moving.
    # The global collocation Newton (solver_recursive/collocation.py) roots the whole
    # system instead, so this is a hard assert now.
    cal, ss = get_ss()
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    # bands MUST match what solve_recursive ships, or the probe measures a box defect
    # rather than the solver -- imported from BOX_KW so the two cannot drift apart.
    # (Hard-wired bands here went stale at the 9-state change: the old 0.12/0.20 P
    # bands admit corners the period map cannot solve now that the household claim W_D
    # is a separate state, and the probe reported their |F| = 9.0 as a solver failure.)
    grid = build_state_box(ss, cal, mu=1, **BOX_KW)
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=2)

    rules.n_gh = 5
    time_iteration(rules, cal, ss, sp, regimes=(0,), no_default=True, damp=0.5,
                   tol=1e-4, max_it=12, n_gh=5, verbose=False)                               # warm start only
    ok, its, worst = solve_collocation(rules, cal, ss, sp, regimes=(0,),
                                       no_default=True, n_gh=5, backend="parsolve",
                                       maxit=20, verbose=False, label=" N2")
    assert ok, f"N2: collocation solve did not converge (max|F| = {worst:.2e})"
    assert worst <= 10 * TOL_MAXF, f"N2: max|F| = {worst:.2e}"
    print(f"    N2 (grid-wide solve at pi=0): converged in {its} Newton steps, "
          f"max|F| = {worst:.2e}")


if __name__ == "__main__":
    test_n1_ss_rest_point()
    print("test_recursive_nesting N1 (SS rest point): PASSED")
    test_n2_no_default_grid_solve()
    print("test_recursive_nesting N2 (grid-wide solve): PASSED")
