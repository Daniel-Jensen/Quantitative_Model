"""Normalised EBA sovereign-exposure panel.

Purpose (paper role)
--------------------
The four EBA exercises are the only public source giving bank x sovereign x
residual-maturity exposures through the Greek crisis. They use different
workbook layouts, different column names and different exposure concepts, so
they are normalised here into one long panel with a single documented schema
and a validator that fails loudly.

One row = one bank x one snapshot x one sovereign issuer.

Concept discipline (plan Decision 2): the BASELINE uses gross direct long
exposure in the banking book, the one concept reported on a comparable basis at
both 2010Q4 and 2011Q3. The 2010Q1 exercise predates the net-direct-position and
banking/trading-book split, so it is robustness, not baseline.

Run:
    /opt/anaconda3/envs/ssj/bin/python -m Empirics.cross_border.eba_panel

Outputs (Empirics/outputs/):
    cb_eba_panel.csv    the normalised long panel
"""

from __future__ import annotations

import os
import re

import pandas as pd
from pandas.api.types import is_numeric_dtype

from Empirics.cross_border import config

PANEL_COLUMNS = (
    "bank_id", "bank_name", "bank_country", "snapshot", "issuer",
    "exposure_nominal", "concept", "book",
    # Pre-crisis conditioning variables for the cross-sectional test. Carried on
    # every row of a bank (redundantly) so the regression frame needs no second
    # merge. "_pre" means the EARLIEST available vintage, never contemporaneous:
    # Greek bank capital is endogenous to the event once HFSF recapitalisation
    # lands in 2012-13.
    "ct1_capital_pre", "ct1_ratio_pre",
) + config.BUCKET_NAMES

PANEL_CSV = os.path.join(config.OUT, "cb_eba_panel.csv")

# Columns that must be numeric before any arithmetic (ladder sums, sign checks,
# ratio bounds) touches them. A string-typed column raises a raw TypeError deep
# inside pandas that points nowhere near the real defect; this check raises
# early, naming the offending column.
_NUMERIC_COLUMNS = (
    ("exposure_nominal", "ct1_capital_pre", "ct1_ratio_pre") + config.BUCKET_NAMES
)

# Exposure "concept" and accounting "book" vocabulary actually reported by the
# EBA exercises (see data_raw/README.md): gross vs net direct long exposure,
# banking vs trading book. A value outside this set is a parse error, not a
# reporting choice, and must not silently fall out of any downstream mask.
_VALID_CONCEPTS = frozenset({"gross_direct_long", "net_direct_long"})
_VALID_BOOKS = frozenset({"banking", "trading"})

# issuer / bank_country are ISO-2 country codes. The EBA panel legitimately
# carries ~30 sovereign issuers (including non-euro, non-EU ones), so this
# checks only the SHAPE of the code, not membership in a small whitelist -
# config.DOMESTIC/CORE/PERIPHERY_EX_GR/STRESS_ISSUERS name only the handful of
# countries the structural model cares about, not the issuer universe.
_ISO2_RE = re.compile(r"^[A-Z]{2}$")

# ct1_ratio_pre is a percentage point (e.g. 8.5 means 8.5%), matching the
# convention fixed by the test fixtures. Bounds are generous (no real EBA bank
# is anywhere near 100%) but a ratio must be positive and cannot exceed 100.
_CT1_RATIO_MAX_PCT = 100.0

# Ladder-vs-total tolerance. EBA workbooks round to whole millions, so an exact
# match is not attainable. The tolerance is 0.5% of the exposure, floored at
# 1.0 unit to absorb whole-million rounding on large banks - but the floor
# itself is capped at a modest fraction of the exposure (20%) so it cannot
# swamp the check for small exposures (e.g. a 1.0 floor against a 0.4 exposure
# would tolerate a 250% relative error).
_LADDER_TOL_FRAC = 0.005
_LADDER_TOL_FLOOR = 1.0
_LADDER_TOL_FLOOR_MAX_FRAC = 0.2


def validate_panel(df: pd.DataFrame) -> None:
    """Raise if `df` is not a well-formed exposure panel."""
    missing = [c for c in PANEL_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"panel is missing columns: {missing}")

    for col in _NUMERIC_COLUMNS:
        if not is_numeric_dtype(df[col]):
            raise ValueError(
                f"column {col!r} must be numeric, got dtype {df[col].dtype}"
            )

    bad_snap = sorted(set(df["snapshot"]) - set(config.EBA_SNAPSHOTS))
    if bad_snap:
        raise ValueError(
            f"unknown snapshot labels {bad_snap}; expected "
            f"{sorted(config.EBA_SNAPSHOTS)}"
        )

    for col in ("issuer", "bank_country"):
        bad_shape = ~df[col].astype(str).str.match(_ISO2_RE)
        if bad_shape.any():
            bad_values = sorted(set(df.loc[bad_shape, col].astype(str)))
            raise ValueError(
                f"{col} must be a two-letter uppercase ISO-style code; found "
                f"{bad_values[:5]}"
            )

    bad_concept = ~df["concept"].isin(_VALID_CONCEPTS)
    if bad_concept.any():
        bad_values = sorted(set(df.loc[bad_concept, "concept"].astype(str)))
        raise ValueError(
            f"unknown concept values {bad_values}; expected one of "
            f"{sorted(_VALID_CONCEPTS)}"
        )

    bad_book = ~df["book"].isin(_VALID_BOOKS)
    if bad_book.any():
        bad_values = sorted(set(df.loc[bad_book, "book"].astype(str)))
        raise ValueError(
            f"unknown book values {bad_values}; expected one of "
            f"{sorted(_VALID_BOOKS)}"
        )

    dupes = df.duplicated(subset=["bank_id", "snapshot", "issuer"])
    if dupes.any():
        rows = df.loc[dupes, ["bank_id", "snapshot", "issuer"]].to_dict("records")
        raise ValueError(f"duplicate bank/snapshot/issuer rows: {rows[:5]}")

    # NaN must be caught before any arithmetic: a NaN in exposure_nominal makes
    # `< 0` evaluate False, and a NaN bucket is silently treated as zero by
    # `.sum()`, so a shifted-header-row parse (everything NaN) would otherwise
    # pass both the sign check and the ladder-sum check.
    if df["exposure_nominal"].isna().any():
        bad_ids = df.loc[df["exposure_nominal"].isna(), "bank_id"].tolist()
        raise ValueError(f"exposure_nominal contains NaN for bank_id(s) {bad_ids[:5]}")

    bucket_cols = list(config.BUCKET_NAMES)
    bucket_nan = df[bucket_cols].isna().any(axis=1)
    if bucket_nan.any():
        bad_ids = df.loc[bucket_nan, "bank_id"].tolist()
        raise ValueError(
            f"bucket columns {bucket_cols} contain NaN for bank_id(s) {bad_ids[:5]}"
        )

    gross = df["concept"] == "gross_direct_long"
    if (df.loc[gross, "exposure_nominal"] < 0).any():
        raise ValueError("negative exposure_nominal on a gross_direct_long row")

    bad_bucket_sign = (df.loc[gross, bucket_cols] < 0).any(axis=1)
    if bad_bucket_sign.any():
        bad_ids = df.loc[gross].loc[bad_bucket_sign, "bank_id"].tolist()
        raise ValueError(
            f"negative bucket value(s) on a gross_direct_long row for "
            f"bank_id(s) {bad_ids[:5]}"
        )

    # cross_section divides by ct1_capital_pre; a zero or missing value there
    # would silently produce inf rather than an error. No upper plausibility
    # bound is imposed here: the panel's unit convention (millions vs billions)
    # is set by the Task 5 adapters, not fixed yet, so a bound picked now could
    # reject a genuine EBA bank on a unit guess rather than a real defect.
    if df["ct1_capital_pre"].isna().any() or (df["ct1_capital_pre"] <= 0).any():
        raise ValueError(
            "ct1_capital_pre must be present and strictly positive on every row"
        )

    # ct1_ratio_pre mirrors ct1_capital_pre's NaN/positivity check, plus an
    # upper bound: it is unambiguously a percentage point (see convention
    # above), so 100 is a safe generous ceiling no real EBA bank approaches.
    bad_ratio = (
        df["ct1_ratio_pre"].isna()
        | (df["ct1_ratio_pre"] <= 0)
        | (df["ct1_ratio_pre"] > _CT1_RATIO_MAX_PCT)
    )
    if bad_ratio.any():
        raise ValueError(
            "ct1_ratio_pre must be present, strictly positive and no greater "
            f"than {_CT1_RATIO_MAX_PCT} (percentage points) on every row"
        )

    ladder_sum = df[bucket_cols].sum(axis=1)
    floor = (df["exposure_nominal"].abs() * _LADDER_TOL_FLOOR_MAX_FRAC).clip(
        upper=_LADDER_TOL_FLOOR
    )
    tol = pd.concat(
        [df["exposure_nominal"].abs() * _LADDER_TOL_FRAC, floor], axis=1
    ).max(axis=1)
    off = (ladder_sum - df["exposure_nominal"]).abs() > tol
    if off.any():
        first = df.loc[off].iloc[0]
        raise ValueError(
            f"ladder does not sum to exposure_nominal for bank "
            f"{first['bank_id']} / {first['snapshot']} / {first['issuer']}: "
            f"ladder {ladder_sum[off].iloc[0]:.3f} vs exposure "
            f"{first['exposure_nominal']:.3f}"
        )


def ladder_of(row) -> dict:
    """Extract one row's residual-maturity ladder as a bucket-keyed dict."""
    return {name: float(row[name]) for name in config.BUCKET_NAMES}
