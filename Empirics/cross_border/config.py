"""Shared configuration for the cross-border sovereign ownership project.

Purpose (paper role)
--------------------
Empirical component #1 of the TPI paper. Measures who actively bought and sold
Greek government bonds over 2009Q4-2012Q1, separating portfolio choice from
passive maturity run-off, and supplies the moment the structural model's
cross-border retrenchment is held to.

Design spec: docs/superpowers/specs/2026-09-08-cross-border-ownership-design.md

This module holds constants only. It performs no I/O.
"""

from __future__ import annotations

import os

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))
EMPIRICS = os.path.dirname(HERE)
OUT = os.path.join(EMPIRICS, "outputs")
RAW = os.path.join(HERE, "data_raw")

# --------------------------------------------------------------------------- #
# Sample
# --------------------------------------------------------------------------- #
# The March 2012 PSI mechanically converts every private holding, so it is the
# terminal event, not an observation. Post-PSI data is descriptive only.
SAMPLE_START = "2009-10-01"
SAMPLE_END = "2012-03-31"

# EBA exposure snapshots: label -> (reference date, source exercise)
EBA_SNAPSHOTS = {
    "2010Q1": ("2010-03-31", "EBA 2010 EU-wide stress test"),
    "2010Q4": ("2010-12-31", "EBA 2011 EU-wide stress test"),
    "2011Q3": ("2011-09-30", "EBA December 2011 capital exercise"),
    "2012Q2": ("2012-06-30", "EBA 2012 capital exercise final update"),
}

# Baseline decomposition interval (see plan Decision 2).
BASELINE_FROM = "2010Q4"
BASELINE_TO = "2011Q3"

# Horizon in YEARS between the baseline snapshots: 2010-12-31 -> 2011-09-30.
BASELINE_HORIZON_YEARS = 0.75

# --------------------------------------------------------------------------- #
# Residual-maturity buckets as reported by the EBA exercises.
# (name, lower edge in years, upper edge in years)
# The 10Y+ bucket has no reported upper edge; 30.0 is the conventional cap and
# is immaterial at horizons below one year (nothing in it matures).
# --------------------------------------------------------------------------- #
BUCKETS = (
    ("0_3M",     0.0,  0.25),
    ("3M_1Y",    0.25, 1.0),
    ("1Y_2Y",    1.0,  2.0),
    ("2Y_3Y",    2.0,  3.0),
    ("3Y_5Y",    3.0,  5.0),
    ("5Y_10Y",   5.0,  10.0),
    ("10Y_PLUS", 10.0, 30.0),
)

BUCKET_NAMES = tuple(b[0] for b in BUCKETS)

# Within-bucket placement assumptions. "uniform" is the point estimate;
# "near" and "far" are the bounds the coarse ladder forces us to report.
PLACEMENTS = ("near", "uniform", "far")

# --------------------------------------------------------------------------- #
# Holder groups. ISO-2 country of the reporting bank's nationality.
# --------------------------------------------------------------------------- #
DOMESTIC = "GR"
CORE = ("DE", "FR", "NL", "BE", "AT")
PERIPHERY_EX_GR = ("IT", "ES", "PT", "IE")

# Sovereign issuers whose exposure enters the balance-sheet-stress control.
# Greece is EXCLUDED by construction: including it would make the control
# mechanically correlated with the dependent variable.
STRESS_ISSUERS = ("IT", "ES", "PT", "IE")
