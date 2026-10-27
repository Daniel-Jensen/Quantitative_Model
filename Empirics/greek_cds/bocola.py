"""Bocola (2016) CDS -> physical sovereign default probability: exact reconstruction.

Source of record: Luigi Bocola, "The Pass-Through of Sovereign Risk," JPE 124(4),
2016, 879-926, together with his replication archive (``data/bocola2016/``).

The chain, in his notation, is

    (18)   CDS_t = phat^d_t * Delta / R_interbank,t
    (19)   p^d_t = CDS_t * lambda [ lev_t (1 - mu_t) + mu_t ] / ( Delta * E_t[Lhat_{t+1} | d_{t+1}=1] )
    (21)   Lhat_{t,t+1} = b * exp(-dc_{t+1}) * [ (1 - omega) + omega * lambda * lev_{t+1} ]

phat(.) is the probability measure defined in his eq. (17); it is the banker's
risk-adjusted measure, NOT the physical one. Footnote 15 gives

    R_interbank = { a (1 - mu) + lambda mu } / E[Lhat],      a = lambda * lev

so the bracket in (19) is exactly R_interbank * E[Lhat] and (19) is the ordinary
change of measure

    p^d = phat^d * E[Lhat] / E[Lhat | d=1].

We therefore name the ratio

    CORRECTION_t = lambda [ lev_t (1-mu_t) + mu_t ] / E_t[Lhat_{t+1} | d_{t+1}=1]
                 = R_interbank,t * E_t[Lhat_{t+1}] / E_t[Lhat_{t+1} | d_{t+1}=1]

which is the column Bocola stores in ``Data.xls``. It is below one because the
banker's marginal value of wealth is high in a default state, so the physical
probability lies below the risk-adjusted one.

Everything here reproduces the archive to machine zero except the single step
from Lhat to E[Lhat | d=1], which the published paper defers to an online
appendix that is not in the archive. That step is handled by
``esdf_projection`` and is accurate to 0.3% on his own Italian series.
"""

import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ARCHIVE = os.path.join(ROOT, "data", "bocola2016", "Data.xls")

# ---------------------------------------------------------------------------
# Parameters recovered from the archive, all to machine zero.
# ---------------------------------------------------------------------------
LAMBDA = 0.21     # eq. (19)/(21); = alpha_bg / lev_bg with lev_bg = 5
HAIRCUT = 0.55    # Delta, "consistent with the Greek debt restructuring of 2012"
SDF_B = 0.99820   # b in eq. (21)
SDF_OMEGA = 0.95  # omega in eq. (21)
R_GROSS = 1.003   # R_bg, gross quarterly risk-free rate (Table 1)


def sdf(dc, lev, b=SDF_B, omega=SDF_OMEGA, lam=LAMBDA):
    """Equation (21): the banker's stochastic discount factor Lhat_{t,t+1}."""
    return b * np.exp(-np.asarray(dc, float)) * (
        (1.0 - omega) + omega * lam * np.asarray(lev, float)
    )


def multiplier(spread_q, lev):
    """Equation (B2) with the archive's convention (no division by R).

    ``spread_q`` is the 1-week EURIBOR minus EUREPO expressed as a difference of
    GROSS QUARTERLY rates, exactly as the archive stores it.
    """
    x = np.asarray(spread_q, float) * np.asarray(lev, float)
    return x / (1.0 + x)


def correction(lev, mu, esdf, lam=LAMBDA):
    """R_interbank * E[Lhat] / E[Lhat | d=1]: the measure change in eq. (19)."""
    lev = np.asarray(lev, float)
    mu = np.asarray(mu, float)
    return lam * (lev * (1.0 - mu) + mu) / np.asarray(esdf, float)


def physical_default_prob(cds_q, lev, mu, esdf, haircut=HAIRCUT, lam=LAMBDA):
    """Equation (19). ``cds_q`` is the 1-year CDS spread as a QUARTERLY rate.

    A CDS quoted at ``s`` basis points per annum enters as ``s / 1e4 / 4``.
    """
    return np.asarray(cds_q, float) * correction(lev, mu, esdf, lam) / haircut


# ---------------------------------------------------------------------------
# Archive readers
# ---------------------------------------------------------------------------
def _sheet(name, cols, names):
    df = pd.read_excel(ARCHIVE, sheet_name=name, header=None).iloc[11:, cols]
    df.columns = ["per"] + names
    df = df[df.per.notna()].copy()
    per = df.pop("per")
    try:
        idx = pd.PeriodIndex(per, freq="Q")
    except Exception:
        idx = pd.PeriodIndex(pd.to_datetime(per), freq="Q")
    return df.set_index(idx).apply(pd.to_numeric, errors="coerce")


def archive_default_probability():
    """Bocola's ``Default Probability`` sheet: his own inputs and his own answer."""
    return _sheet("Default Probability", [0, 2, 4, 6, 8, 10, 12],
                  ["cds", "mu", "lev", "esdf", "correction", "p"])


def archive_lagrange_multiplier():
    """His ``Lagrange_Multiplier`` sheet, which carries dc and the SDF itself."""
    return _sheet("Lagrange_Multiplier", [0, 2, 4, 6, 8, 10, 14, 16],
                  ["spread", "lev", "mu", "dy", "dc", "sdf", "sdf_nolev"])


def archive_italian_cds_monthly():
    """Monthly Markit CDS term structure, Italy. Columns are tenors in years."""
    raw = pd.read_excel(ARCHIVE, sheet_name="Italian CDS", header=None)
    tenors = {4: "1y", 6: "2y", 8: "3y", 12: "5y", 14: "7y", 16: "10y",
              18: "15y", 20: "20y", 22: "30y", 2: "6m"}
    out = raw.iloc[11:, [0] + list(tenors)]
    out.columns = ["per"] + list(tenors.values())
    out = out[out.per.notna()].copy()
    per = out.pop("per")
    return out.set_index(pd.PeriodIndex(per, freq="M")).apply(pd.to_numeric, errors="coerce")


# ---------------------------------------------------------------------------
# The one step the archive does not pin down.
# ---------------------------------------------------------------------------
def fit_esdf_projection():
    """Project E[Lhat | d=1] on the same state Lhat uses.

    The paper computes E_t[Lhat_{t+1} | d_{t+1}=1] with a model-based step that
    only the (unavailable) online appendix describes. On his own data the object
    is affine in leverage once the consumption term is divided out:

        E[Lhat | d=1] = b exp(-dc) * (a0 + a1 * lev)

    Returns (a0, a1, max_abs_error). Economically the fitted coefficients say
    E[Lhat | d=1] sits about 13% above Lhat, i.e. the banker's marginal value of
    wealth is a seventh higher in a default state than on average.
    """
    dp = archive_default_probability()
    lm = archive_lagrange_multiplier()
    d = pd.DataFrame({"esdf": dp.esdf, "lev": lm.lev, "dc": lm.dc}).dropna()
    y = (d.esdf / (SDF_B * np.exp(-d.dc))).values
    A = np.c_[np.ones(len(d)), d.lev.values]
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    err = np.abs(y - A @ coef).max() * SDF_B
    return float(coef[0]), float(coef[1]), float(err)


def esdf_projection(dc, lev, coef=None):
    """E[Lhat | d=1] evaluated out of sample, e.g. on Greek data."""
    a0, a1, _ = coef if coef is not None else fit_esdf_projection()
    return SDF_B * np.exp(-np.asarray(dc, float)) * (a0 + a1 * np.asarray(lev, float))
