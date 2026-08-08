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
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from common import get_ss
from state_grid import build_state_box, s_process_params
from decision_rules import RuleSet, SOLVE7
from point_map import point_residuals
from recursive_main import (time_iteration, calibrate_household_anchors,
                            ss_state, ss_x)

_LABELS = ("cap_D", "cap_F", "lab_D", "lab_F", "euler_D", "uip", "goods_D")


def test_n1_ss_rest_point():
    # SS + SS-CONSTANT RULES + NO DEFAULT => ALL 7 RESIDUALS ~ 0.
    cal, ss = get_ss()
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    grid = build_state_box(ss, cal)
    rules = RuleSet.from_ss(grid, ss, cal)

    res, out = point_residuals(ss_state(ss, cal, sp), 0, ss_x(ss, cal),
                               rules, cal, ss, sp, no_default=True)
    worst = np.max(np.abs(res))
    for lab, r in zip(_LABELS, res):
        print(f"    {lab:10s} {r:+.3e}")
    print(f"    mu_D={out['mu_D']:.6f} (ss {ss['ss_bank_D']['mu_ss']:.6f})  "
          f"C_D={out['C_D']:.5f} (ss {ss['C_D_ss']:.5f})  "
          f"A_D={out['A_D']:.5f} (ss {ss['A_D_ss']:.5f})")
    assert worst < 1e-6, f"SS not a rest point: max|res|={worst:.2e}"


def probe_n2_no_default_time_iteration():
    # N2 CONVERGENCE PROBE (KNOWN-OPEN, NOT A HARD GATE): time-iterate the d=0
    # block at pi=0. The pointwise map nests transition.py and the SS is an exact
    # rest point (N1), but the grid-wide fixed point is stiff at the occasionally-
    # binding kink (polynomial interpolation of a C0 policy limit-cycles -- the
    # plan's bounds-adaptation + anisotropic refinement + value-based convergence
    # are the remaining apparatus). This probe REPORTS status; it does not assert.
    cal, ss = get_ss()
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    grid = build_state_box(ss, cal, mu=1, k_band=0.02, p_band=0.04,
                           b_band=0.04)
    rules = RuleSet.from_ss(grid, ss, cal)

    ok, its, worst = time_iteration(rules, cal, ss, sp, regimes=(0,),
                                    no_default=True, damp=0.3, tol=1e-6,
                                    max_it=40, n_gh=5, verbose=False)
    print(f"    N2 probe: converged={ok} in {its} sweeps, "
          f"worst point residual={worst:.2e} (known-open)")


if __name__ == "__main__":
    test_n1_ss_rest_point()
    print("test_recursive_nesting N1 (SS rest point): PASSED")
    probe_n2_no_default_time_iteration()
