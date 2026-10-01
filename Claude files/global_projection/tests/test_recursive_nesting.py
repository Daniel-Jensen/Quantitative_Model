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
from solver_recursive.state_grid import build_state_box, s_process_params
from solver_recursive.decision_rules import RuleSet, SOLVE7
from solver_recursive.point_map import point_residuals
from solver_recursive.recursive_main import (time_iteration, calibrate_household_anchors,
                            ss_state, ss_x)
from solver_recursive.recursive_experiment import BOX_KW
from solver_recursive.collocation import (solve_collocation, TOL_MAXF, RES_NAMES,
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
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=4)

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
    # THE BACKSTOP MUST BE ASLEEP AT THE STEADY STATE. It is inactive in normal times by
    # construction (the peg is the rest-point price, which the SS price sits above), so a
    # non-zero CB position here would mean the backstop is intervening in a state it has
    # no business in -- and every reported IRF is then measured against a polluted base.
    assert out["m_ltro_D"] == 0.0, (
        f"facility drawn in a no-backstop regime at the SS: {out['m_ltro_D']:.3e}")


def test_n3_phi_zero_nests_the_no_cb_model():
    # N3: phi = 0 MUST REPRODUCE THE NO-BACKSTOP MODEL EXACTLY, NOT TO TOLERANCE.
    # This is the same doctrine as pi = 0 nesting the risk-free model and
    # size_F = size_D the symmetric one: a policy switch at its off value must leave
    # ZERO trace, so any measured effect of the backstop is the backstop.
    cal, ss = get_ss()
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    grid = build_state_box(ss, cal, mu=1, **BOX_KW)
    x0 = ss_x(ss, cal)
    r2 = RuleSet.from_ss(grid, ss, cal, n_regimes=2)
    r4 = RuleSet.from_ss(grid, ss, cal, n_regimes=4)
    r2.n_gh = r4.n_gh = 5
    # TWO INDEPENDENT OFF SWITCHES. phi = 0 removes the facility from the EXPECTATION;
    # ltro = 0 removes it from the CONSTRAINT. Testing both separately is what makes any
    # measured effect attributable to the backstop rather than to the extra regimes or
    # the wider basis.
    #
    # THEY HOLD TO DIFFERENT TOLERANCES, AND THE DIFFERENCE IS NOT SLOPPINESS.
    # At phi = 0 the facility regimes carry weight vectors that are identically zero, so
    # _regime_weights aliases them to regime 0, np.dot contributes an exact 0.0 and the
    # sum is BIT-IDENTICAL to the two-regime model. At phi > 0 with ltro = 0 the economy
    # is the same but the ARITHMETIC is not: the same expectation is accumulated as
    # (1-phi)*a + phi*a instead of a, which in floating point differs by an ULP. Demanding
    # bit-identity there would be demanding that addition be associative. One ULP of the
    # residual is 14 orders below the acceptance floor.
    m0 = cal["ltro_D"]
    for label, tol, kw in (("phi_ltro = 0", 0.0, dict(phi_ltro=0.0, ltro_D=m0)),
                           ("ltro = 0", 1e-14, dict(phi_ltro=0.75, ltro_D=0.0))):
        cal.update(kw); cal["ltro_F"] = 0.0
        worst = 0.0
        for i in range(grid.n):
            S = grid.points[i]
            a, _ = point_residuals(S, 0, x0, r2, cal, ss, sp, n_gh=5)
            b, _ = point_residuals(S, 0, x0, r4, cal, ss, sp, n_gh=5)
            worst = max(worst, float(np.max(np.abs(a - b))))
        assert worst <= tol, f"{label} does not nest: max gap {worst:.3e} > {tol:.0e}"
        how = "EXACT (bit-for-bit)" if tol == 0.0 else f"max gap {worst:.1e} (~1 ULP)"
        print(f"    N3 ({label}): {how} over {grid.n} points x "
              f"{N_RES_POINT} equations")
    cal["phi_ltro"], cal["ltro_D"] = 0.0, m0


def test_n4_the_facility_touches_only_the_constraint():
    # N4: THE LTRO MUST MOVE THE INCENTIVE CONSTRAINT AND NOTHING ELSE.
    # The design claim is that lending at the deposit rate changes the COMPOSITION of the
    # bank's funding, not its size or its cost, so every budget identity is untouched and
    # the whole effect is in mu. That claim is cheap to assert and expensive to get wrong:
    # if the facility leaked into the funding side, resources would appear from nowhere
    # and Walras would leak -- the failure mode that dogged the earlier purchase design.
    # Comparing the SAME state and the SAME policy vector with the facility off and on,
    # the deposit-market residual and the carried obligation must be BIT-IDENTICAL, while
    # the multiplier must weakly fall and strictly fall somewhere.
    cal, ss = get_ss()
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    grid = build_state_box(ss, cal, mu=1, **BOX_KW)
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=4)
    rules.n_gh = 5
    x0 = ss_x(ss, cal)
    cal["phi_ltro"], cal["ltro_F"] = 0.5, 0.0
    m = cal["ltro_D"]
    i_dep = _LABELS.index("dep_clear")
    n_relieved = 0
    for i in range(grid.n):
        S = grid.points[i]
        cal["ltro_D"] = 0.0
        r_off, o_off = point_residuals(S, 1, x0, rules, cal, ss, sp, n_gh=5)
        cal["ltro_D"] = m
        r_on, o_on = point_residuals(S, 1, x0, rules, cal, ss, sp, n_gh=5)
        assert r_off[i_dep] == r_on[i_dep], (
            f"the facility moved deposit clearing at point {i}: "
            f"{r_off[i_dep]:.17e} vs {r_on[i_dep]:.17e}")
        for k in ("dep_D", "Pp_D", "Vp_dep", "nfa_dep_D", "n_D"):
            assert o_off[k] == o_on[k], f"the facility moved {k} at point {i}"
        assert o_on["m_ltro_D"] == m, "facility not drawn in the CB regime"
        assert o_on["mu_D"] <= o_off["mu_D"] + 1e-15, (
            f"the facility TIGHTENED the constraint at point {i}: "
            f"{o_off['mu_D']:.6e} -> {o_on['mu_D']:.6e}")
        assert o_on["slack_D"] >= o_off["slack_D"] - 1e-12
        n_relieved += o_on["mu_D"] < o_off["mu_D"]
    cal["phi_ltro"] = 0.0
    assert n_relieved > 0, "the facility never relieved the constraint anywhere"
    print(f"    N4: deposit clearing, dep_D, P', V' and n_D bit-identical with the "
          f"facility on; mu strictly lower at {n_relieved}/{grid.n} points")


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
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=4)

    rules.n_gh = 5
    time_iteration(rules, cal, ss, sp, regimes=(0,), no_default=True, damp=0.5,
                   tol=1e-4, max_it=12, n_gh=5, verbose=False,
                   no_cb=True)                                      # warm start only
    ok, its, worst = solve_collocation(rules, cal, ss, sp, regimes=(0,),
                                       no_default=True, n_gh=5, backend="parsolve",
                                       maxit=20, verbose=False, label=" N2",
                                       no_cb=True)
    assert ok, f"N2: collocation solve did not converge (max|F| = {worst:.2e})"
    assert worst <= 10 * TOL_MAXF, f"N2: max|F| = {worst:.2e}"
    print(f"    N2 (grid-wide solve at pi=0): converged in {its} Newton steps, "
          f"max|F| = {worst:.2e}")


if __name__ == "__main__":
    test_n1_ss_rest_point()
    print("test_recursive_nesting N1 (SS rest point): PASSED")
    test_n3_phi_zero_nests_the_no_cb_model()
    print("test_recursive_nesting N3 (nesting): PASSED")
    test_n4_the_facility_touches_only_the_constraint()
    print("test_recursive_nesting N4 (facility touches only the constraint): PASSED")
    test_n2_no_default_grid_solve()
    print("test_recursive_nesting N2 (grid-wide solve): PASSED")
