# CROSS-BORDER BOND FOCS: THE p'/p VALUATION AND THE SS FRICTION REFERENCES.
# X1  each bank values its foreign bond in its own good: the D bank's F-bond FOC carries p'/p
#     and the F bank's D-bond FOC p/p'. Raising every continuation p by (1+g), with nothing else
#     moved, shifts residual 13 by exactly g*payF and residual 9 by (1/(1+g) - 1)*payD at the SS.
# X2  with UNEQUAL bank spreads (what the rest-point calibration delivers) the deterministic SS is
#     still an exact rest point of the period map: the friction references absorb the cross-border
#     return gap, and the foreign leg earns its issuer's spread in the bank identity.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import numpy as np

from common import get_ss
from calibration.global_projection import get_calibration
from global_projection.steady_state import solve_steady_state
from global_projection.solver_recursive.state_grid import s_process_params
from global_projection.solver_recursive.decision_rules import RuleSet
from global_projection.solver_recursive.point_map import point_residuals
from global_projection.solver_recursive.recursive_main import (calibrate_household_anchors,
                                                               ss_state, ss_x)
from global_projection.solver_recursive.recursive_experiment import normal_box
from global_projection.solver_recursive.collocation import RES_NAMES

I9, I13 = RES_NAMES.index("bondD_F"), RES_NAMES.index("bondF_D")


def _ss_rest(cal, ss):
    # THE SS STATE, THE SS UNKNOWNS AND SS-CONSTANT RULES, WITH THE HOUSEHOLD ANCHORS SET.
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    rules = RuleSet.from_ss(normal_box(ss, cal, mu=1), ss, cal)
    return sp, rules, ss_state(ss, cal, sp), ss_x(ss, cal)


def test_x1_foreign_bond_valued_at_p_ratio():
    # X1: THE CONTINUATION p ENTERS THE TWO CROSS-BORDER FOCS EXACTLY AS p/p' AND p'/p.
    # The raw FOC gaps, not the residuals: the no-short-sale pairs transform the residual.
    cal, ss = get_ss()
    sp, rules, S0, x0 = _ss_rest(cal, ss)
    _, o0 = point_residuals(S0, 0, x0, rules, cal, ss, sp, no_default=True)
    g = 0.01
    for d in (0, 1):
        rules.set_values("p", d, np.full(rules.grid.n, ss["p_ss"] * (1.0 + g)))
    _, o1 = point_residuals(S0, 0, x0, rules, cal, ss, sp, no_default=True)
    payF = cal["delta_b_F"] + (1.0 - cal["delta_b_F"]) * ss["Q_bF_ss"]
    payD = cal["delta_b_D"] + (1.0 - cal["delta_b_D"]) * ss["Q_bD_ss"]
    d13 = o1["foc_FD_gap"] - o0["foc_FD_gap"]
    d9 = o1["foc_DF_gap"] - o0["foc_DF_gap"]
    assert abs(d13 - g * payF) < 1e-12, f"D bank's F-bond FOC: {d13:.6e} vs {g * payF:.6e}"
    assert abs(d9 - (1.0 / (1.0 + g) - 1.0) * payD) < 1e-12, f"F bank's D-bond FOC: {d9:.6e}"
    print(f"    X1: p' up 1% moves res13 by {d13:+.6f} (= g*payF), res9 by {d9:+.6f} (= -g/(1+g)*payD)")


def test_x2_unequal_spreads_rest_point():
    # X2: s_D != s_F, AND THE DETERMINISTIC SS IS STILL A REST POINT OF THE PERIOD MAP.
    cal = get_calibration()
    cal["credit_spread_target_D"] = cal["credit_spread_target_D"] + 10.0 / 4e4   # +10 bp/yr
    ss = solve_steady_state(cal, verbose=False)
    sp, rules, S0, x0 = _ss_rest(cal, ss)
    res, out = point_residuals(S0, 0, x0, rules, cal, ss, sp, no_default=True)
    worst = float(np.max(np.abs(res)))
    # the no-short-sale pairs' smoothing leaves ~book_eps^2/(2*book) ~ 1e-7 at the SS, as in N1
    assert worst < 1e-6, f"SS not a rest point with unequal spreads: max|res| = {worst:.2e}"
    # without the references the cross-border FOCs carry the whole spread gap
    c2 = dict(cal, b_D_F_ref=cal["b_D_F_ss"], b_F_D_ref=cal["b_F_D_ss"])
    res2, _ = point_residuals(S0, 0, x0, rules, c2, ss, sp, no_default=True)
    assert min(abs(res2[I9]), abs(res2[I13])) > 1e-6, "the references did nothing"
    print(f"    X2: spreads D {4e4*cal['credit_spread_target_D']:.1f} / F "
          f"{4e4*cal['credit_spread_target_F']:.1f} bp: max|res| = {worst:.1e} "
          f"(without the references: res9 {res2[I9]:+.1e}, res13 {res2[I13]:+.1e})")


if __name__ == "__main__":
    test_x1_foreign_bond_valued_at_p_ratio()
    print("test_cross_border X1 (p'/p valuation): PASSED")
    test_x2_unequal_spreads_rest_point()
    print("test_cross_border X2 (unequal spreads rest point): PASSED")
