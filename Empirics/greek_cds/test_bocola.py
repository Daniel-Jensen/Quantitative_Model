"""Gates on the Bocola (2016) CDS -> physical default probability reconstruction.

/opt/anaconda3/envs/ssj/bin/python -m pytest Empirics/greek_cds/ -q
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from greek_cds import bocola as B  # noqa: E402


@pytest.fixture(scope="module")
def dp():
    return B.archive_default_probability().dropna(subset=["p"])


@pytest.fixture(scope="module")
def lm():
    return B.archive_lagrange_multiplier()


def test_equation_19_reproduces_his_default_probability(dp):
    """p^d = CDS * lambda[lev(1-mu)+mu] / (Delta * E[Lhat|d=1]), machine zero."""
    p = B.physical_default_prob(dp.cds, dp.lev, dp.mu, dp.esdf)
    assert np.abs(p - dp.p).max() < 1e-15


def test_correction_column_is_the_change_of_measure(dp):
    c = B.correction(dp.lev, dp.mu, dp.esdf)
    assert np.abs(c - dp.correction).max() == 0.0


def test_lambda_is_exactly_0p21_and_haircut_0p55(dp):
    """Both are over-identified by the archive: 39 quarters, two unknowns."""
    implied = dp.p * B.HAIRCUT * dp.esdf / (dp.cds * (dp.lev * (1 - dp.mu) + dp.mu))
    assert np.abs(implied - 0.21).max() < 5e-6


def test_equation_21_reproduces_his_sdf(lm):
    """Lhat = b exp(-dc)[(1-omega) + omega*lambda*lev], b=.9982 omega=.95.

    1999Q1 is excluded: his own first cell is stale by 0.27%, implying a
    leverage of 6.5131 against the 6.4949 stored one row away.
    """
    d = lm.loc["1999Q2":].dropna(subset=["sdf", "dc", "lev"])
    assert np.abs(B.sdf(d.dc, d.lev) - d.sdf).max() < 1e-9


def test_sdf_discount_parameter_is_exactly_0p9982(lm):
    d = lm.dropna(subset=["sdf_nolev", "dc"])
    b = d.sdf_nolev * np.exp(d.dc)
    assert np.abs(b - B.SDF_B).max() < 1e-12


def test_physical_probability_lies_below_the_risk_adjusted_one(dp):
    """The correction must be below one.

    Default is a state in which the banker's marginal value of wealth is high,
    so the risk-adjusted measure over-weights it and the physical probability
    is the smaller of the two.  A correction above one is not a large number;
    it is a signal that the leverage input is not a steady-state reading.
    """
    assert (dp.correction < 1.0).all()
    assert 0.85 < dp.correction.min() and dp.correction.max() < 0.95


def test_esdf_projection_tracks_his_own_series_within_half_a_percent():
    """The Lhat -> E[Lhat|d=1] step is the only one the archive does not pin.

    It is an affine projection, not his exact calculation, so it is gated as a
    bound rather than asserted exact.  Half a percent on E[Lhat|d=1] is half a
    percent on the default probability - two orders below the question being
    asked of it.
    """
    a0, a1, err = B.fit_esdf_projection()
    assert err < 0.5e-2 * 2.0
    assert 1.07 < (a0 + a1 * 8.0) / ((1 - B.SDF_OMEGA) + B.SDF_OMEGA * B.LAMBDA * 8.0) < 1.20


def test_units_a_one_basis_point_cds_is_a_tiny_quarterly_probability():
    """Guards the two places units go wrong: bp->decimal, and annual->quarterly."""
    cds_q = 1.0 / 1e4 / 4.0                      # 1bp per annum
    p = B.physical_default_prob(cds_q, lev=8.0, mu=0.0007,
                                esdf=B.esdf_projection(0.0, 8.0))
    assert 3e-5 < p < 5e-5                       # ~0.004% per quarter


def test_italian_precrisis_probability_is_far_below_the_model_calibration(dp):
    """Bocola's own Italy: 0.014% per quarter, a seventh of def_rate_ss = 0.10%."""
    pre = dp.loc["2002Q2":"2007Q2"].p.mean()
    assert 1.0e-4 < pre < 2.0e-4
    assert pre / 0.0010 < 0.25


def test_the_one_year_cds_is_the_contract_and_it_is_not_the_five_year():
    """Do not silently substitute tenors: they differ by a factor of ~3."""
    m = B.archive_italian_cds_monthly().loc["2002-04":"2007-06"]
    assert m["1y"].mean() < m["5y"].mean()
    assert 2.0 < m["5y"].mean() / m["1y"].mean() < 5.0


# ---------------------------------------------------------------------------
# The Greek CDS extract. Licensed vendor data, so it is git-ignored and these
# checks skip where it has not been pulled.
# ---------------------------------------------------------------------------
from greek_cds import cds_input, fetch_lseg, greek  # noqa: E402

pytestmark_needs_extract = pytest.mark.skipif(
    fetch_lseg.load() is None, reason="no LSEG extract; run fetch_lseg.py")


@pytestmark_needs_extract
def test_greek_cds_has_no_pre_crisis_observation():
    """The availability finding, asserted so it cannot be forgotten.

    Greek sovereign CDS on the 2003 protocol opens 2008-07-21 and ends at the
    PSI credit event. Bocola's pre-crisis window has no Greek counterpart, and
    any claim about Greek default risk before 2008 is an inference from bonds,
    never a CDS observation.
    """
    st = cds_input.staleness_report()["GR_1Y"]
    assert st["first"] == cds_input.FIRST_QUOTE
    assert st["last"] == cds_input.LAST_QUOTE
    assert st["share_unchanged"] < 0.25       # a price, not a monthly mark


@pytestmark_needs_extract
def test_this_composite_reproduces_bocolas_markit_series():
    """Same contract, different composite: correlated but not identical.

    The level ratio is not noise - it is a vendor difference that passes
    one-for-one into the default probability, so it is asserted as a bound and
    quoted wherever the Greek numbers are.
    """
    v = cds_input.validate_against_bocola(B.archive_default_probability().cds)
    assert v["quarters"] >= 12
    assert v["corr"] > 0.98
    assert 1.05 < v["level_ratio"] < 1.25


@pytestmark_needs_extract
def test_bocolas_linear_cds_formula_breaks_at_greek_spreads():
    """A numerical-approximation issue, not an economic result.

    Equation (18) prices a one-period bill with the premium paid regardless of
    survival, so the implied probability is linear in the spread and crosses
    one. At Italian spreads the two forms agree to better than a per cent; at
    the Greek 2011Q4 quote the linear form returns 55% against 34%, and by
    2012Q1 it returns a probability above one.
    """
    it = cds_input.load_italian_cds("1Y").dropna()
    lin = greek.risk_adjusted_prob(it, B.HAIRCUT, survival=False)
    sur = greek.risk_adjusted_prob(it, B.HAIRCUT, survival=True)
    assert (np.abs(lin - sur) / sur).max() < 0.05

    gr = cds_input.load_greek_cds("1Y").dropna()
    lin_g = greek.risk_adjusted_prob(gr, B.HAIRCUT, survival=False)
    sur_g = greek.risk_adjusted_prob(gr, B.HAIRCUT, survival=True)
    assert lin_g.max() > 1.0          # the linear form is unusable here
    assert sur_g.max() < 1.0          # the survival form is not


@pytestmark_needs_extract
def test_the_greek_leverage_correction_is_dead_before_the_cds_series_starts():
    """The two halves of eq. (19) do not overlap for Greece.

    Greek MFI leverage passes 20 in 2008Q4 - before any sovereign event - so the
    correction exceeds one from that quarter on, which is impossible when
    default is a high-marginal-value state. Only 2008Q3 has both a CDS quote and
    a leverage reading that is still a solvency ratio.
    """
    df = greek.build()
    g = greek.greek_conversion(df)
    live = g.loc["2008Q3":"2012Q1"]
    usable = live.lev.apply(lambda x: np.isfinite(x) and 0 < x < 25)
    assert usable.sum() <= 8
    assert (g.correction.loc["2008Q4":"2011Q4"] > 1.0).all()


@pytestmark_needs_extract
def test_the_first_observable_greek_quarter_sits_at_the_model_calibration():
    """2008Q3: Greek 1-year CDS 22bp, physical p^d about 0.09% per quarter.

    Close enough to the model's 0.10% that the calibration cannot be rejected
    on it - and, separately, close to the 24bp break-even derived from the
    transformation before the data were obtained.
    """
    df = greek.build()
    ans = greek.greek_answer(df, greek.greek_conversion(df))
    r = ans.loc["2008Q3"]
    assert 20.0 < r.cds_1y_bp < 25.0
    assert 0.0007 < r.p_D55 < 0.0011
