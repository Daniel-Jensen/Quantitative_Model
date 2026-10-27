"""Tests for the active/passive holdings decomposition."""

import pytest

from Empirics.cross_border import config
from Empirics.cross_border.reallocation import decompose


def _flat_ladder(value_per_bucket):
    return {name: value_per_bucket for name in config.BUCKET_NAMES}


def test_identity_closes_exactly():
    d = decompose(stock_start=70.0, stock_end=55.0,
                  ladder_start=_flat_ladder(10.0), horizon=0.75)
    assert d["passive"] + d["active"] == pytest.approx(d["stock_end"], abs=1e-12)
    assert d["stock_start"] - d["passive"] == pytest.approx(d["redeemed"], abs=1e-12)


def test_holder_that_only_let_bonds_mature_shows_zero_active():
    ladder = _flat_ladder(10.0)
    passive = 10.0 * 5 + 10.0 * (1.0 - (0.75 - 0.25) / 0.75)  # 1Y+ buckets, then 3M-1Y
    d = decompose(stock_start=70.0, stock_end=passive,
                  ladder_start=ladder, horizon=0.75)
    assert d["active"] == pytest.approx(0.0, abs=1e-9)


def test_selling_beyond_run_off_is_negative_active():
    d = decompose(stock_start=70.0, stock_end=10.0,
                  ladder_start=_flat_ladder(10.0), horizon=0.75)
    assert d["active"] < 0.0


def test_buying_shows_positive_active_even_when_the_stock_falls():
    # The whole point of the exercise: a falling stock can be net BUYING once
    # scheduled redemptions are stripped out.
    ladder = {name: 0.0 for name in config.BUCKET_NAMES}
    ladder["0_3M"] = 100.0          # all of it matures inside the horizon
    d = decompose(stock_start=100.0, stock_end=30.0,
                  ladder_start=ladder, horizon=0.75)
    assert d["passive"] == pytest.approx(0.0)
    assert d["active"] == pytest.approx(30.0)
    assert d["active_pct_of_start"] == pytest.approx(30.0)


def test_active_pct_is_relative_to_the_starting_stock():
    d = decompose(stock_start=200.0, stock_end=0.0,
                  ladder_start={**{n: 0.0 for n in config.BUCKET_NAMES},
                                "10Y_PLUS": 200.0},
                  horizon=0.75)
    # Nothing matures, so the entire fall is active selling: -100% of the start.
    assert d["active_pct_of_start"] == pytest.approx(-100.0)


def test_zero_starting_stock_gives_a_null_percentage_not_a_divide_by_zero():
    d = decompose(stock_start=0.0, stock_end=0.0,
                  ladder_start={n: 0.0 for n in config.BUCKET_NAMES},
                  horizon=0.75)
    assert d["active_pct_of_start"] is None


def test_bounds_are_returned_and_ordered():
    d = decompose(stock_start=70.0, stock_end=40.0,
                  ladder_start=_flat_ladder(10.0), horizon=0.75)
    # "far" postpones redemption, so more survives passively, so the same
    # observed end-stock implies MORE active selling. Hence active_far is the
    # lower bound on active reallocation.
    assert d["active_far"] <= d["active"] <= d["active_near"]


def test_ladder_inconsistent_with_the_starting_stock_is_rejected():
    with pytest.raises(ValueError, match="ladder"):
        decompose(stock_start=100.0, stock_end=50.0,
                  ladder_start=_flat_ladder(10.0), horizon=0.75)  # ladder sums to 70
