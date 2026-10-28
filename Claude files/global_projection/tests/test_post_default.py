# THE POST-DEFAULT BOX: ITS RULE SET, THE BLENDED CONTINUATION AND THE TWO-SET UNKNOWN VECTOR.
# P1  the blend weight is 1 at the post-default cluster and 0 at the SS, monotone in between, with
#     a continuous slope at both ends of the handover (a quintic step along the recovery line).
# P2  a post-default set that copies the normal rules on the normal box changes no continuation.
# P3  pack/unpack round-trip with both sets; the residual stacks set by set.
# P4  copy() copies the post-default set and re-links it to the copy.
# P5  freshly seeded, the post-default rules reproduce the old edge read wherever they hold alone.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import numpy as np

from common import get_ss
from global_projection.solver_recursive.state_grid import s_process_params
from global_projection.solver_recursive.decision_rules import RuleSet, STORE_RULES
from global_projection.solver_recursive.recursive_main import calibrate_household_anchors, ss_state
from global_projection.solver_recursive.recursive_experiment import (normal_box, make_blend,
                                                                     _attach_post, BOX_SPECS)
from global_projection.solver_recursive import collocation as C


def _rules():
    # THE PRODUCTION NORMAL BOX WITH STEADY-STATE RULES.
    cal, ss = get_ss()
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    rules = RuleSet.from_ss(normal_box(ss, cal, mu=1), ss, cal)
    rules.n_gh = 5
    return cal, ss, sp, rules


def test_p1_blend_weight():
    # P1: 1 AT THE CLUSTER, 0 AT THE SS, MONOTONE, C1 AT THE HANDOVER ENDS.
    cal, ss, sp, rules = _rules()
    _attach_post(rules, ss, cal)
    b, S0 = rules.blend, ss_state(ss, cal, sp)

    def w(p):
        # THE WEIGHT AT THE STATE p NET-WORTH UNITS ALONG THE RECOVERY LINE FROM THE SS
        S = S0.copy()
        S[b["block"]] = b["centre"] + p * b["direction"] / b["weights"]
        return float(rules.blend_weight(S)[0])

    assert w(BOX_SPECS["blend"]["cluster"]) == 1.0 and w(0.0) == 0.0
    ps = np.linspace(BOX_SPECS["blend"]["cluster"] - 0.1, 0.1, 401)
    ws = np.array([w(p) for p in ps])
    assert np.all((ws >= 0.0) & (ws <= 1.0)) and np.all(np.diff(ws) <= 1e-15)
    h = 1e-6
    for p0 in (b["post_until"], b["normal_from"]):
        assert abs(w(p0 + h) - w(p0 - h)) / (2 * h) < 1e-3, p0
    print(f"    P1: weight 1 up to {b['post_until']:+.3f}, 0 from {b['normal_from']:+.3f} "
          f"(net-worth units along the recovery line), flat at both ends")


def test_p2_identical_sets_change_nothing():
    # P2: BLENDING A SET WITH ITS OWN COPY IS THE PLAIN READ.
    cal, ss, sp, rules = _rules()
    rules.attach_post(rules.copy(), make_blend(ss, cal))
    x = rules.grid.from_unit(np.random.default_rng(3).uniform(-1.5, 1.5, size=(64, rules.grid.d)))
    for d in (0, 1):
        r, r0 = rules.eval_cont(d, x), rules.eval_all(d, rules.grid.clip(x))
        for k in STORE_RULES:
            assert np.allclose(r[k], r0[k], rtol=1e-13, atol=1e-15), k
    print("    P2: a post-default copy of the normal rules leaves every continuation unchanged")


def test_p3_two_set_packing():
    # P3: THE UNKNOWN VECTOR AND THE RESIDUAL WITH BOTH SETS.
    cal, ss, sp, rules = _rules()
    _attach_post(rules, ss, cal)
    n0, n2 = rules.grid.n, rules.post.grid.n
    theta = C.pack(rules, (0, 1))
    assert theta.size == len(STORE_RULES) * 2 * (n0 + n2)
    back = C.unpack(theta, rules, (0, 1))
    assert len(back) == 2
    for rs, vals in zip((rules, rules.post), back):
        for k in STORE_RULES:
            for d in (0, 1):
                assert np.allclose(vals[k][d], rs.vals[k][d], rtol=0, atol=1e-12), k
    assert C.pack(rules, (0, 1), sets=("post",)).size == len(STORE_RULES) * 2 * n2
    F = C.make_residual(rules, cal, ss, sp, regimes=(0, 1), n_gh=5)
    assert F(theta).size == (n0 + n2) * 2 * C.N_RES
    print(f"    P3: {theta.size} unknowns over {n0} + {n2} points round-trip exactly")


def test_p4_copy_relinks():
    # P4: A COPY OWNS ITS POST-DEFAULT SET.
    cal, ss, sp, rules = _rules()
    _attach_post(rules, ss, cal)
    c = rules.copy()
    assert c.post is not rules.post and c.post.home is c and c.blend is rules.blend
    for k in STORE_RULES:
        assert np.array_equal(c.post.vals[k][1], rules.post.vals[k][1]), k
    print("    P4: copy() copies and re-links the post-default set")


def test_p5_seed_reproduces_the_edge_read():
    # P5: WHERE THE POST-DEFAULT RULES HOLD ALONE, THE SEEDED CONTINUATION IS THE OLD ONE.
    cal, ss, sp, rules = _rules()
    _attach_post(rules, ss, cal)
    pts = rules.post.grid.points
    alone = rules.blend_weight(pts) == 1.0
    assert alone.sum() >= 5, "too few post-default nodes beyond the handover"
    for d in (0, 1):
        r = rules.eval_cont(d, pts[alone])
        r0 = rules.eval_all(d, rules.grid.clip(pts[alone]))
        for k in STORE_RULES:
            assert np.allclose(r[k], r0[k], rtol=1e-10, atol=1e-12), k
    print(f"    P5: at the {int(alone.sum())} post-default nodes beyond the handover the seeded "
          f"continuation equals the old edge read")


if __name__ == "__main__":
    test_p1_blend_weight(); print("test_post_default P1 (blend weight): PASSED")
    test_p2_identical_sets_change_nothing(); print("test_post_default P2 (identical sets): PASSED")
    test_p3_two_set_packing(); print("test_post_default P3 (two-set packing): PASSED")
    test_p4_copy_relinks(); print("test_post_default P4 (copy): PASSED")
    test_p5_seed_reproduces_the_edge_read(); print("test_post_default P5 (seed): PASSED")
