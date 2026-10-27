"""Tests for the normalised EBA exposure panel."""

import pandas as pd
import pytest

from Empirics.cross_border import config
from Empirics.cross_border.eba_panel import PANEL_COLUMNS, validate_panel, ladder_of


def _minimal_row(**over):
    row = {
        "bank_id": "DE001",
        "bank_name": "Test Bank AG",
        "bank_country": "DE",
        "snapshot": "2010Q4",
        "issuer": "GR",
        "exposure_nominal": 70.0,
        "concept": "gross_direct_long",
        "book": "banking",
        "ct1_capital_pre": 22.8,
        "ct1_ratio_pre": 8.5,
    }
    row.update({name: 10.0 for name in config.BUCKET_NAMES})
    row.update(over)
    return row


def test_panel_columns_are_exactly_the_documented_schema():
    assert PANEL_COLUMNS == (
        "bank_id", "bank_name", "bank_country", "snapshot", "issuer",
        "exposure_nominal", "concept", "book",
        "ct1_capital_pre", "ct1_ratio_pre",
    ) + config.BUCKET_NAMES


def test_a_wellformed_panel_validates():
    df = pd.DataFrame([_minimal_row()])
    validate_panel(df)  # must not raise


def test_missing_column_is_rejected():
    df = pd.DataFrame([_minimal_row()]).drop(columns=["issuer"])
    with pytest.raises(ValueError, match="missing columns"):
        validate_panel(df)


def test_unknown_snapshot_label_is_rejected():
    df = pd.DataFrame([_minimal_row(snapshot="2011Q1")])
    with pytest.raises(ValueError, match="snapshot"):
        validate_panel(df)


def test_ladder_that_does_not_sum_to_the_exposure_is_rejected():
    df = pd.DataFrame([_minimal_row(exposure_nominal=999.0)])
    with pytest.raises(ValueError, match="ladder"):
        validate_panel(df)


def test_duplicate_bank_snapshot_issuer_rows_are_rejected():
    df = pd.DataFrame([_minimal_row(), _minimal_row()])
    with pytest.raises(ValueError, match="duplicate"):
        validate_panel(df)


def test_negative_exposure_is_rejected_for_a_gross_long_concept():
    df = pd.DataFrame([_minimal_row(exposure_nominal=-70.0)])
    with pytest.raises(ValueError, match="negative"):
        validate_panel(df)


def test_zero_capital_is_rejected_rather_than_dividing_by_zero_downstream():
    df = pd.DataFrame([_minimal_row(ct1_capital_pre=0.0)])
    with pytest.raises(ValueError, match="ct1_capital_pre"):
        validate_panel(df)


def test_ladder_of_returns_a_dict_keyed_by_bucket_name():
    df = pd.DataFrame([_minimal_row()])
    ladder = ladder_of(df.iloc[0])
    assert set(ladder) == set(config.BUCKET_NAMES)
    assert ladder["0_3M"] == 10.0


# --------------------------------------------------------------------------- #
# Regression tests added after code review (finding numbers per the review).
# --------------------------------------------------------------------------- #

# --- 1: NaN in exposure_nominal / bucket columns must be rejected before any
# arithmetic touches them (the shifted-header-row failure mode). ---

def test_nan_exposure_nominal_is_rejected():
    df = pd.DataFrame([_minimal_row(exposure_nominal=float("nan"))])
    with pytest.raises(ValueError, match="NaN"):
        validate_panel(df)


def test_single_nan_bucket_value_is_rejected():
    df = pd.DataFrame([_minimal_row(**{"0_3M": float("nan")})])
    with pytest.raises(ValueError, match="NaN"):
        validate_panel(df)


def test_all_nan_row_is_rejected():
    over = {"exposure_nominal": float("nan")}
    over.update({name: float("nan") for name in config.BUCKET_NAMES})
    df = pd.DataFrame([_minimal_row(**over)])
    with pytest.raises(ValueError, match="NaN"):
        validate_panel(df)


# --- 2: per-bucket non-negativity for gross rows, independent of the total. ---

def test_negative_bucket_value_is_rejected_even_when_ladder_sums_correctly():
    # 0_3M=-100, 3M_1Y=170, remaining 5 buckets at the _minimal_row default of
    # 10.0 each (=50): ladder sums to 120, matched to exposure_nominal so the
    # *total* check alone would wave this through.
    over = {"0_3M": -100.0, "3M_1Y": 170.0, "exposure_nominal": 120.0}
    df = pd.DataFrame([_minimal_row(**over)])
    with pytest.raises(ValueError, match="negative"):
        validate_panel(df)


# --- 3 / 9: concept and book must be whitelisted, independent of the sign check. ---

def test_unknown_concept_is_rejected():
    df = pd.DataFrame([_minimal_row(concept="bogus_concept")])
    with pytest.raises(ValueError, match="concept"):
        validate_panel(df)


def test_unknown_concept_does_not_silently_disable_the_negative_exposure_guard():
    # Ladder kept consistent with the negative total (each bucket -10.0, summing
    # to -70.0) so the ladder-tolerance check cannot fire and mask the real bug:
    # under the old code, an unrecognised concept drops the row out of the
    # `gross` mask entirely, so NEITHER the exposure- nor bucket-level negative
    # check ever sees it, and the row validates clean.
    over = {"concept": "bogus_concept", "exposure_nominal": -70.0}
    over.update({name: -10.0 for name in config.BUCKET_NAMES})
    df = pd.DataFrame([_minimal_row(**over)])
    with pytest.raises(ValueError, match="concept"):
        validate_panel(df)


def test_unknown_book_is_rejected():
    df = pd.DataFrame([_minimal_row(book="mezzanine")])
    with pytest.raises(ValueError, match="book"):
        validate_panel(df)


# --- 4: ct1_ratio_pre validated like its sibling ct1_capital_pre. ---

def test_nan_ct1_ratio_pre_is_rejected():
    df = pd.DataFrame([_minimal_row(ct1_ratio_pre=float("nan"))])
    with pytest.raises(ValueError, match="ct1_ratio_pre"):
        validate_panel(df)


def test_negative_ct1_ratio_pre_is_rejected():
    df = pd.DataFrame([_minimal_row(ct1_ratio_pre=-5.0)])
    with pytest.raises(ValueError, match="ct1_ratio_pre"):
        validate_panel(df)


def test_implausibly_large_ct1_ratio_pre_is_rejected():
    df = pd.DataFrame([_minimal_row(ct1_ratio_pre=250.0)])
    with pytest.raises(ValueError, match="ct1_ratio_pre"):
        validate_panel(df)


# --- 5 (issuer/bank_country shape, replacing the rejected whitelist idea): ---

def test_issuer_with_wrong_shape_is_rejected():
    df = pd.DataFrame([_minimal_row(issuer="Greece")])
    with pytest.raises(ValueError, match="issuer"):
        validate_panel(df)


def test_bank_country_with_wrong_shape_is_rejected():
    df = pd.DataFrame([_minimal_row(bank_country="de")])
    with pytest.raises(ValueError, match="bank_country"):
        validate_panel(df)


# --- 6: numeric-dtype assertions, naming the offending column. ---

def test_string_typed_exposure_nominal_is_rejected_with_a_clear_message():
    df = pd.DataFrame([_minimal_row(exposure_nominal="70.0")])
    with pytest.raises(ValueError, match="exposure_nominal"):
        validate_panel(df)


def test_string_typed_bucket_column_is_rejected_with_a_clear_message():
    df = pd.DataFrame([_minimal_row(**{"0_3M": "ten"})])
    with pytest.raises(ValueError, match="0_3M"):
        validate_panel(df)


def test_string_typed_ct1_capital_pre_is_rejected_with_a_clear_message():
    df = pd.DataFrame([_minimal_row(ct1_capital_pre="22.8")])
    with pytest.raises(ValueError, match="ct1_capital_pre"):
        validate_panel(df)


def test_string_typed_ct1_ratio_pre_is_rejected_with_a_clear_message():
    df = pd.DataFrame([_minimal_row(ct1_ratio_pre="8.5")])
    with pytest.raises(ValueError, match="ct1_ratio_pre"):
        validate_panel(df)


# --- 8: the ladder tolerance floor must not swamp the relative check at small
# exposures, while still absorbing whole-million rounding at large ones. ---

def test_large_exposure_still_tolerates_whole_million_rounding():
    # exposure=100, ladder off by 0.8 (rounding-sized) -> must NOT raise.
    over = {"exposure_nominal": 100.0, "0_3M": 40.8}
    df = pd.DataFrame([_minimal_row(**over)])
    validate_panel(df)  # must not raise


def test_small_exposure_catches_a_large_relative_ladder_error():
    # exposure=0.4, ladder off by 0.9999 (250% relative error) -> must raise.
    over = {"exposure_nominal": 0.4}
    over.update({name: 0.0 for name in config.BUCKET_NAMES})
    over["0_3M"] = 1.3999
    df = pd.DataFrame([_minimal_row(**over)])
    with pytest.raises(ValueError, match="ladder"):
        validate_panel(df)
