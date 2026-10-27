"""Maturity-ladder roll-forward: the passive-holdings counterfactual.

Purpose (paper role)
--------------------
Greek debt matured and was redeemed at par with official money. A foreign bank
whose Greek holdings fell 40% may have made no portfolio decision at all - it
held to maturity while the Troika refinanced. Raw holdings changes therefore
measure the maturity ladder, not retrenchment.

This module answers one question: how much face value survives from t to t+h if
the holder buys nothing? Everything else is active portfolio choice, and only
active portfolio choice is comparable to the model's cross-border position,
which is a quantity chosen each period rather than a stock that runs off.

The EBA residual-maturity buckets are coarse, so the within-bucket distribution
is unobserved. We report a uniform-placement point estimate bracketed by two
bounds: "near" (all mass at the near edge, fastest possible run-off) and "far"
(all mass at the far edge, slowest). No I/O, no network.

Precondition: face values must be non-negative (gross exposure). The
near <= uniform <= far bound guarantee only holds for non-negative weights;
a negative entry (e.g. a mis-parsed "net direct position") silently inverts
which bound is larger, so `passive_stock` rejects it rather than returning a
backwards answer.
"""

from __future__ import annotations

import math

from Empirics.cross_border import config


def redeemed_share(lo: float, hi: float, horizon: float, placement: str) -> float:
    """Share of one bucket's face value that has matured by `horizon` years.

    Parameters
    ----------
    lo, hi
        Bucket edges in years, e.g. (0.25, 1.0) for the 3M-1Y bucket.
    horizon
        Years forward from the snapshot date. Must not be NaN: a NaN fails
        both edge comparisons below, so it would otherwise fall through to
        the interior branch and produce a silent, wrong "near"/"far" answer
        instead of an error.
    placement
        Within-bucket distribution assumption: "uniform", "near" or "far".
    """
    if placement not in config.PLACEMENTS:
        raise ValueError(
            f"unknown placement {placement!r}; expected one of {config.PLACEMENTS}"
        )
    if math.isnan(horizon):
        raise ValueError("horizon is NaN")
    if horizon <= lo:
        return 0.0
    if horizon >= hi:
        return 1.0
    # The horizon falls strictly inside the bucket, so the assumption bites.
    if placement == "uniform":
        return (horizon - lo) / (hi - lo)
    if placement == "near":
        return 1.0  # all mass at `lo`, which the horizon has already passed
    return 0.0      # "far": all mass at `hi`, which the horizon has not reached


def passive_stock(ladder: dict, horizon: float, placement: str = "uniform") -> float:
    """Face value surviving to `horizon` with zero purchases.

    `ladder` maps every name in `config.BUCKET_NAMES` to a face value, and
    every value must be non-negative (gross exposure, not a net position). A
    missing bucket raises rather than defaulting to zero: a silently absent
    bucket would understate the surviving stock and overstate retrenchment,
    which is the exact error this module exists to prevent. A negative entry
    raises for the same reason from the other direction: it silently inverts
    which of "near"/"far" is the larger bound, breaking the near <= uniform
    <= far guarantee this module promises.
    """
    total = 0.0
    for name, lo, hi in config.BUCKETS:
        value = ladder[name]  # KeyError on a missing bucket is deliberate
        if value < 0.0:
            raise ValueError(
                f"negative face value {value!r} in bucket {name!r}; "
                "passive_stock requires gross (non-negative) exposures"
            )
        total += value * (1.0 - redeemed_share(lo, hi, horizon, placement))
    return total
