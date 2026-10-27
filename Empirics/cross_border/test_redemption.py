"""Tests for the maturity-ladder roll-forward."""

import math

import pytest

from Empirics.cross_border import config
from Empirics.cross_border.redemption import redeemed_share, passive_stock


def test_bucket_entirely_beyond_horizon_is_untouched():
    # 10Y+ bucket, three-quarter horizon: nothing matures.
    assert redeemed_share(10.0, 30.0, 0.75, "uniform") == 0.0


def test_bucket_entirely_inside_horizon_is_fully_redeemed():
    # 0-3M bucket, three-quarter horizon: all of it has matured.
    assert redeemed_share(0.0, 0.25, 0.75, "uniform") == 1.0


def test_uniform_placement_prorates_a_straddling_bucket():
    # 3M-1Y bucket is 0.75 years wide; a horizon of 0.625 is halfway through it.
    assert redeemed_share(0.25, 1.0, 0.625, "uniform") == pytest.approx(0.5)


def test_near_placement_redeems_all_of_a_straddling_bucket():
    # All mass sits at the near edge (0.25), which is inside a 0.625 horizon.
    assert redeemed_share(0.25, 1.0, 0.625, "near") == 1.0


def test_far_placement_redeems_none_of_a_straddling_bucket():
    # All mass sits at the far edge (1.0), which is beyond a 0.625 horizon.
    assert redeemed_share(0.25, 1.0, 0.625, "far") == 0.0


def test_zero_horizon_leaves_the_whole_ladder_standing():
    ladder = {name: 10.0 for name in config.BUCKET_NAMES}
    assert passive_stock(ladder, 0.0, "uniform") == pytest.approx(70.0)


def test_all_short_ladder_is_fully_run_off():
    ladder = {name: 0.0 for name in config.BUCKET_NAMES}
    ladder["0_3M"] = 42.0
    assert passive_stock(ladder, 0.75, "uniform") == pytest.approx(0.0)


def test_all_long_ladder_survives_intact():
    ladder = {name: 0.0 for name in config.BUCKET_NAMES}
    ladder["10Y_PLUS"] = 42.0
    assert passive_stock(ladder, 0.75, "uniform") == pytest.approx(42.0)


def test_placements_bound_the_uniform_point_estimate():
    # "far" postpones maturity, so the MOST survives; "near" the least.
    ladder = {name: 10.0 for name in config.BUCKET_NAMES}
    near = passive_stock(ladder, 0.625, "near")
    unif = passive_stock(ladder, 0.625, "uniform")
    far = passive_stock(ladder, 0.625, "far")
    assert near <= unif <= far
    assert near < far  # the bounds must be informative, not degenerate


def test_unknown_placement_is_rejected():
    with pytest.raises(ValueError, match="placement"):
        redeemed_share(0.25, 1.0, 0.5, "middle")


def test_missing_bucket_is_rejected_rather_than_silently_zero():
    ladder = {"0_3M": 1.0}  # six buckets absent
    with pytest.raises(KeyError):
        passive_stock(ladder, 0.75, "uniform")


def test_missing_last_bucket_is_rejected_not_just_missing_first():
    # Same guard, but omitting the LAST bucket (10Y_PLUS) rather than an early
    # one, so the test actually exercises name-keyed lookup rather than
    # position-keyed lookup that happens to stop early.
    ladder = {name: 1.0 for name in config.BUCKET_NAMES if name != "10Y_PLUS"}
    with pytest.raises(KeyError):
        passive_stock(ladder, 0.75, "uniform")


@pytest.mark.parametrize("placement", ["uniform", "near", "far"])
def test_horizon_exactly_on_lower_edge_redeems_nothing(placement):
    # horizon == lo: the bucket has not started maturing yet, for any
    # placement assumption. Pins the `<=` boundary against a future flip.
    assert redeemed_share(0.25, 1.0, 0.25, placement) == 0.0


@pytest.mark.parametrize("placement", ["uniform", "near", "far"])
def test_horizon_exactly_on_upper_edge_redeems_everything(placement):
    # horizon == hi: the bucket has fully matured, for any placement
    # assumption. Pins the `>=` boundary against a future flip.
    assert redeemed_share(0.25, 1.0, 1.0, placement) == 1.0


def test_negative_face_value_is_rejected_rather_than_inverting_bounds():
    # A negative entry (e.g. a mis-parsed "net direct position" instead of
    # gross exposure) silently flips near/far if left unchecked: negating a
    # weight reverses which bound is larger. The module's whole guarantee is
    # near <= uniform <= far, which holds only for non-negative face values.
    ladder = {name: 10.0 for name in config.BUCKET_NAMES}
    ladder["3M_1Y"] = -50.0
    with pytest.raises(ValueError, match="3M_1Y"):
        passive_stock(ladder, 0.625, "uniform")


def test_nan_horizon_is_rejected_for_every_placement():
    # NaN fails both edge comparisons (`<=` and `>=` are False for NaN), so
    # execution falls through to the interior branch. "uniform" would at
    # least produce a loud NaN, but "near"/"far" would silently return
    # 1.0/0.0 without ever reading `horizon`. Must be rejected up front.
    for placement in config.PLACEMENTS:
        with pytest.raises(ValueError, match="horizon"):
            redeemed_share(0.25, 1.0, math.nan, placement)
