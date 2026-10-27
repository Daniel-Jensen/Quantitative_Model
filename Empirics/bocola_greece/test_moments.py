"""Gates on the Bocola reconstruction and the Greek counterparts.

The point of these tests is that the *reconstruction* claims in
docs/greek_bocola_calibration.md are checkable. They run offline: the Bocola
numbers come from the archived Data.xls, the Greek ones from the cached series
CSV. Nothing here touches the model.
"""
import os

import numpy as np
import pandas as pd
import pytest

from . import moments as M


# ── the reconstruction of Bocola's own construction ──────────────────────────
@pytest.fixture(scope="module")
def italy():
    return M.bocola_italy()


def test_leverage_identity_is_exact(italy):
    """lev = TA / (TA - TL + F51_liab + F52_liab), Appendix B.C, to machine zero."""
    assert np.abs(italy.lev - italy.lev_xls).max() == 0.0


def test_multiplier_equation_is_exact_given_the_spread(italy):
    """eq. (B2) in the archive's convention, to machine zero."""
    mu = M.mu_of(italy.spread_xls, italy.lev, divide_by_R=False)
    assert np.abs(mu - italy.mu_xls).dropna().max() == 0.0


def test_spread_from_stored_rate_levels_is_a_rounding_artefact(italy):
    """Re-deriving the spread from the 6-dp rate levels costs ~1e-5, not more."""
    err = np.abs(italy.spread - italy.spread_xls).dropna().max()
    assert err < 1e-4, err
    assert err > 0.0, "if this is exactly zero the levels are not rounded"


def test_published_and_archive_multiplier_conventions_agree_closely(italy):
    """The paper divides the spread by R_t; the spreadsheet does not."""
    a = M.mu_of(italy.spread_xls, italy.lev, divide_by_R=False)
    b = M.mu_of(italy.spread_xls, italy.lev, divide_by_R=True)
    rel = ((a - b).abs() / a.abs()).dropna()
    assert rel.max() < 0.005, rel.max()


def test_eight_basis_points_is_the_figure1_window_annualised(italy):
    """The paper's '8 basis points' = mean Euribor-Eurepo, 2002Q2-2012Q4, x4."""
    lo, hi = M.WINDOWS["bocola_figure1_2002Q2_2012Q4"]
    sl = italy.loc[(italy.index >= lo) & (italy.index <= hi)]
    bp = sl.spread_xls.mean() * 4e4
    assert 7.5 < bp < 9.0, bp


def test_estimation_window_spread_is_materially_smaller(italy):
    """2002Q2-2009Q4, the window actually estimated, averages ~5bp not ~8bp."""
    lo, hi = M.WINDOWS["bocola_estimation_2002Q2_2009Q4"]
    sl = italy.loc[(italy.index >= lo) & (italy.index <= hi)]
    assert 4.5 < sl.spread_xls.mean() * 4e4 < 5.6


def test_table1_leverage_is_not_the_measured_italian_leverage(italy):
    """lev_bg = 5 is a literature convention; the data say ~8."""
    assert M.BOCOLA_TABLE1["lev_bg"][0] == 5.0
    assert italy.lev.mean() > 7.5


def test_b2_at_table1_leverage_and_the_8bp_spread_returns_the_estimate():
    """(lev=5, 8.2bp/yr) -> mu = 0.001, which is the Table 2 posterior mean."""
    mu = M.mu_of(8.19e-4 / 4, M.BOCOLA_TABLE1["lev_bg"][0])
    assert abs(mu - M.BOCOLA_MU_POSTERIOR["mean"]) < 5e-5, mu


def test_wedge_inverts_the_multiplier():
    for mu in (1e-4, 1e-3, 1.2309e-2):
        s = M.wedge_of(mu, 5.0)
        assert abs(M.mu_of(s, 5.0) - mu) < 1e-12


def test_italian_top5_exposure_moment_reproduces_table1():
    e = M.eba_exposure_moments()["Italy, Bocola top-5 (archive)"]
    assert abs(e["net_over_ta"] - M.BOCOLA_TABLE1["exp_bg"][0]) < 5e-4
    assert abs(1 / e["avg_maturity_q"] - M.BOCOLA_TABLE1["pi"][0]) < 2e-3


# ── the Greek counterparts ───────────────────────────────────────────────────
@pytest.fixture(scope="module")
def series():
    if not os.path.exists(M.SERIES_CSV):
        pytest.skip("run Empirics/bocola_greece/moments.py once to build the cache")
    return pd.read_csv(M.SERIES_CSV, index_col=0)


def test_eurostat_reproduces_bocola_italian_leverage(series, italy):
    """The comparability test: same definition, different provider, same series.

    If this fails the Greek leverage number is not a Bocola counterpart.
    """
    cmp = pd.concat({"b": italy.lev,
                     "e": pd.to_numeric(series.lev_eurostat_IT, errors="coerce")},
                    axis=1).dropna()
    assert len(cmp) > 50
    assert cmp.b.corr(cmp.e) > 0.95
    assert abs(cmp.b.mean() / cmp.e.mean() - 1) < 0.05


def test_greek_precrisis_leverage_is_close_to_italian(series):
    """~8 both sides: the Italian leverage moment needs no Greek correction."""
    lo, hi = M.WINDOWS["greek_precrisis_2002Q2_2007Q2"]
    sl = series.loc[(series.index >= lo) & (series.index <= hi)]
    gr = pd.to_numeric(sl.lev_eurostat_GR, errors="coerce").mean()
    it = pd.to_numeric(sl.lev_eurostat_IT, errors="coerce").mean()
    assert 6.5 < gr < 9.5, gr
    assert abs(gr - it) < 1.5, (gr, it)


def test_greek_crisis_leverage_is_not_a_steady_state(series):
    """Greek bank equity was all but gone by 2012; >20 is a solvency reading."""
    lo, hi = M.WINDOWS["greek_crisis_2010Q1_2012Q4"]
    sl = series.loc[(series.index >= lo) & (series.index <= hi)]
    assert pd.to_numeric(sl.lev_eurostat_GR, errors="coerce").mean() > 20


def test_the_deposit_premium_proxy_fails_the_sign_test_in_italy(series):
    """Italian retail time-deposit rates sit BELOW the money market pre-crisis.

    A negative funding wedge implies mu < 0, which the model cannot represent.
    This is why the liability-side proxy is classified not-comparable for the
    steady state rather than merely imperfect.
    """
    lo, hi = M.WINDOWS["greek_precrisis_2002Q2_2007Q2"]
    sl = series.loc[(series.index >= lo) & (series.index <= hi)]
    assert pd.to_numeric(sl.wedge_depIT_over_eonia, errors="coerce").mean() < 0
