"""Active/passive decomposition of a holder's sovereign position.

Purpose (paper role)
--------------------
    stock_end = passive_survival(ladder_start, horizon) + active

`active` is portfolio CHOICE - net purchases over the interval - and it is the
only term comparable to the structural model's cross-border position, which is a
quantity chosen each period rather than a stock that runs off. A holder whose
Greek book halved while every bond in it matured did not retrench.

The identity closes by construction (`active` is defined residually), and
`decompose` asserts that closure at 1e-9. That is deliberate: it mirrors the
self-verifying style of experiments/e2_dy_decomposition.py, so a units or
alignment error is caught here rather than in a figure caption.

No I/O, no network.
"""

from __future__ import annotations

from Empirics.cross_border import config
from Empirics.cross_border.redemption import passive_stock

_CLOSURE_TOL = 1e-9


def decompose(stock_start: float, stock_end: float, ladder_start: dict,
              horizon: float) -> dict:
    """Split the change in a holding into scheduled run-off and active choice.

    All quantities are NOMINAL (face value). See plan Decision 1: the EBA
    exercises report nominal exposures with a maturity ladder, which is what the
    roll-forward needs; BIS book-value claims cannot be decomposed this way.

    Returns a dict with the point estimate (`uniform` placement) plus the
    `near`/`far` bounds the coarse ladder forces us to report.
    """
    ladder_total = sum(ladder_start[name] for name in config.BUCKET_NAMES)
    if abs(ladder_total - stock_start) > 1e-6 * max(1.0, abs(stock_start)):
        raise ValueError(
            f"ladder sums to {ladder_total:.6f} but stock_start is "
            f"{stock_start:.6f}; the ladder must exhaust the starting stock"
        )

    passive = passive_stock(ladder_start, horizon, "uniform")
    active = stock_end - passive

    passive_near = passive_stock(ladder_start, horizon, "near")
    passive_far = passive_stock(ladder_start, horizon, "far")

    # More passive survival implies more active selling for the same end-stock,
    # so the "far" placement gives the LOWER bound on active reallocation.
    active_near = stock_end - passive_near
    active_far = stock_end - passive_far

    if abs((passive + active) - stock_end) > _CLOSURE_TOL:
        raise AssertionError(
            f"decomposition failed to close: passive {passive:.12f} + active "
            f"{active:.12f} != stock_end {stock_end:.12f}"
        )

    pct = None if stock_start == 0.0 else 100.0 * active / stock_start

    return {
        "stock_start": stock_start,
        "stock_end": stock_end,
        "passive": passive,
        "redeemed": stock_start - passive,
        "active": active,
        "active_near": active_near,
        "active_far": active_far,
        "active_pct_of_start": pct,
    }
