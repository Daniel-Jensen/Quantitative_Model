# DECISION-RULE LAYER GATE: THE LOG PARAMETERISATION IS A CHANGE OF VARIABLE ONLY.
# Values are stored in levels and only the Chebyshev fit/eval go through log (or
# log-gross for rates), so: node values must round-trip exactly, constants must stay
# exact, the interpolant must be positive everywhere by construction, and a rule with
# a saturated corner must no longer drag the fit at the ergodic centre negative.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from solver_recursive.state_grid import SmolyakGrid
from solver_recursive.decision_rules import (RuleSet, ALL_RULES, LOG_RULES,
                                             GROSS_RULES, to_fit, from_fit)
from common import get_ss


def test_transform_round_trip():
    # to_fit / from_fit are exact inverses on the admissible range.
    v = np.array([0.05, 0.5, 1.0, 3.0, 40.0])
    for name in ("alpha_D", "Q_bD", "C_D", "p", "N_D"):
        assert np.max(np.abs(from_fit(name, to_fit(name, v)) - v)) < 1e-12
    r = np.array([-0.02, 0.0, 0.003, 0.05])            # rates may be negative
    for name in GROSS_RULES:
        assert np.max(np.abs(from_fit(name, to_fit(name, r)) - r)) < 1e-12
    # untransformed rules pass through untouched
    assert np.max(np.abs(from_fit("__none__", to_fit("__none__", v)) - v)) < 1e-15


def test_node_exactness_and_positivity():
    # The fit still reproduces its own node values, and evaluation is positive.
    rng = np.random.default_rng(0)
    g = SmolyakGrid(-np.ones(4), np.ones(4), mu=2)
    rs = RuleSet(g)
    for name in ("alpha_D", "Q_bD", "rdep_D", "p"):
        vals = (rng.uniform(0.2, 5.0, g.n) if name != "rdep_D"
                else rng.uniform(-0.01, 0.05, g.n))
        rs.set_values(name, 0, vals)
        assert np.max(np.abs(rs.eval(name, 0, g.points) - vals)) < 1e-9, name
    x = rng.uniform(-1, 1, size=(500, 4))
    for name in ("alpha_D", "Q_bD", "p"):
        assert np.all(rs.eval(name, 0, x) > 0.0), f"{name} went non-positive off-node"
    assert np.all(rs.eval("rdep_D", 0, x) > -1.0), "gross deposit rate went negative"


def test_constants_are_exact():
    # A constant rule must be reproduced to machine precision (the SS cold start).
    cal, ss = get_ss()
    g = SmolyakGrid(-np.ones(5), np.ones(5), mu=2)
    rs = RuleSet(g)
    rng = np.random.default_rng(1)
    x = rng.uniform(-1, 1, size=(300, 5))
    for name in ALL_RULES:
        c = 0.003 if name in GROSS_RULES else 1.37
        rs.set_values(name, 0, np.full(g.n, c))
        assert np.max(np.abs(rs.eval(name, 0, x) - c)) < 1e-11, name


def test_saturated_corner_does_not_go_negative():
    # The Table-7 failure mode: one node driven to a guard used to move the fit at the
    # centre by ~26% AND make it negative. In logs the interpolant cannot go negative,
    # which is the property the point-map guards used to have to enforce by clipping.
    g = SmolyakGrid(np.array([-1.0]), np.array([1.0]), mu=3)
    rs = RuleSet(g)
    u = g.points[:, 0]
    base = 1.10 + 0.02 * np.exp((u + 1.0) * 1.375)
    vals = base.copy()
    vals[np.argmax(u)] = 40.0                      # corner saturates the alpha cap
    rs.set_values("alpha_D", 0, vals)
    sweep = np.linspace(-1, 1, 601)[:, None]
    got = rs.eval("alpha_D", 0, sweep)
    assert np.all(got > 0.0), f"log fit went non-positive: min {got.min():.4f}"
    assert np.max(np.abs(rs.eval("alpha_D", 0, g.points) - vals)) < 1e-8


if __name__ == "__main__":
    test_transform_round_trip()
    test_node_exactness_and_positivity()
    test_constants_are_exact()
    test_saturated_corner_does_not_go_negative()
    print("test_decision_rules: ALL PASSED")
