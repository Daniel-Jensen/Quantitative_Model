"""The Greek CDS leg of Bocola's conversion, from LSEG Workspace.

Pulled by ``fetch_lseg.py``; see that module for the instrument choice and the
credential handling.  This module only converts the cached extract into the
units equation (19) expects and records what the series can and cannot support.

Availability, established 2026-09-17 against the live entitlement
-----------------------------------------------------------------
Greek sovereign CDS on the 2003 ISDA protocol runs **2008-07-21 to 2012-03-08**,
933 quoted days, ending at the PSI credit event.  Nothing earlier exists on this
entitlement at any tenor or in either currency: the EUR-denominated Greek
contracts return empty, and the 2014-protocol contracts that replaced them carry
no history before 2023.  **There is therefore no pre-crisis Greek sovereign CDS
observation to be had here.**  That is a fact about the market as much as about
the vendor - a Eurosystem study of euro-area sovereign CDS (Fontana and
Scheicher, ECB WP 1271) also could not begin before 2006, and used 10-year
quotes only.

The Italian contract on the same chain starts on the same day, which is what
makes the vendor comparison in ``validate_against_bocola`` possible: over the
2008Q3-2011Q4 overlap this composite correlates 0.993 with Bocola's Markit/WRDS
series but sits about 13% above it in level.
"""

import numpy as np
import pandas as pd

from greek_cds import fetch_lseg

FIRST_QUOTE = "2008-07-21"
LAST_QUOTE = "2012-03-08"        # Greek PSI credit-event auction 2012-03-09


def load_greek_cds(tenor="1Y"):
    """Quarterly Greek CDS as a decimal QUARTERLY rate, or None if not pulled."""
    s = fetch_lseg.quarterly("GR", tenor)
    return None if s is None else s / 1e4 / 4.0


def load_italian_cds(tenor="1Y"):
    s = fetch_lseg.quarterly("IT", tenor)
    return None if s is None else s / 1e4 / 4.0


def validate_against_bocola(archive_cds_q):
    """Compare this composite with Bocola's Markit/WRDS Italian 1-year series.

    ``archive_cds_q`` is his quarterly column as a decimal quarterly rate.
    Returns correlation, mean absolute percentage gap, and the level ratio.
    The level ratio is the number that matters: it is a vendor-composite
    difference, not a measurement error, and it propagates one-for-one into the
    default probability.
    """
    mine = load_italian_cds("1Y")
    if mine is None:
        return None
    d = pd.DataFrame({"bocola": archive_cds_q, "lseg": mine}).dropna()
    if len(d) < 4:
        return None
    return {
        "quarters": int(len(d)),
        "corr": float(d.bocola.corr(d.lseg)),
        "mean_abs_pct_gap": float((100 * (d.lseg - d.bocola).abs() / d.bocola).mean()),
        "level_ratio": float(d.lseg.mean() / d.bocola.mean()),
    }


def market_recovery():
    """The recovery the composite assumes when quoting the par spread.

    Bocola sets Delta = 0.55 off the realised 2012 PSI.  The quoting convention
    is 40% recovery, i.e. Delta = 0.60.  A par spread is only a probability
    jointly with a recovery assumption, so the two must not be mixed silently.
    """
    d = fetch_lseg.load()
    if d is None:
        return None
    return {c.split("_")[0]: float(d[c].dropna().mean() / 100.0)
            for c in d.columns if c.endswith("_5Y_RECOV_RATE")}


def staleness_report():
    d = fetch_lseg.load()
    if d is None:
        return None
    out = {}
    for c in d.columns:
        if not c.endswith("_PAR_MID1"):
            continue
        s = d[c].dropna()
        cty, ten, _ = c.split("_", 2)
        out[f"{cty}_{ten}"] = {
            "quotes": int(len(s)),
            "first": str(s.index[0].date()),
            "last": str(s.index[-1].date()),
            "share_unchanged": float((s.diff() == 0).mean()),
            "mean_bp": float(s.mean()),
        }
    return out


MISSING_DATA_NOTE = """\
No cached LSEG extract. Start Workspace, then:
  LSEG_APP_KEY=$(cat ~/.config/lseg/appkey) \\
    /opt/anaconda3/envs/ssj/bin/python Empirics/greek_cds/fetch_lseg.py
"""
