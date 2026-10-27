"""Pull the Greek and Italian sovereign CDS term structure from LSEG Workspace.

    LSEG_APP_KEY=$(cat ~/.config/lseg/appkey) \
      /opt/anaconda3/envs/ssj/bin/python Empirics/greek_cds/fetch_lseg.py

Requires the Workspace desktop application to be running on this machine: the
App Key authenticates a DESKTOP session, which proxies through Workspace's local
port. It is read from the environment, or failing that from a 0600 file outside
the repository. It is never written into a source file and the extract it
produces is licensed vendor data, so ``data/`` is git-ignored.

Instrument choice
-----------------
Bocola uses Markit 1-year Italian sovereign CDS, RED 4AB951, complete
restructuring, through WRDS. The LSEG counterpart is the composite
``<CTY>GV<TENOR><CCY>A?=R``, where the protocol letter is the thing that matters:

  * ``AB`` = Cum-Restructuring, **2003 ISDA protocol**.  This is the contract in
    force through the sovereign crisis and the one whose clause matches
    Bocola's CR.  History runs 2008-07-21 to the 2012-03-09 Greek credit event.
  * ``AC`` = Cum-Restructuring, **2014 protocol**.  A different contract, and on
    this entitlement it carries no history before 2023 for Greece.

Mixing the two would splice two different contracts across the credit event.
USD and EUR denominations are also distinct contracts and are pulled separately;
the EUR Greek series is empty on this entitlement, so USD is the only option.

``PAR_MID1`` is the par spread in basis points per annum - the field comparable
to Markit's ``parspread``. ``RECOV_RATE`` is the recovery the composite assumes
when quoting, and is worth carrying because it is the market's convention (40%)
rather than Bocola's haircut (45%).
"""

import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "data", "cds_panel.csv")

COUNTRIES = {"GR": "Hellenic Republic", "IT": "Italy"}
TENORS = ["6M", "1Y", "2Y", "3Y", "5Y", "7Y", "10Y"]
PROTOCOL = "AB"          # 2003 ISDA; see module docstring
CCY = "US"
FIELDS = ["PAR_MID1", "RECOV_RATE", "DEFLT_PROB"]
START, END = "2001-01-01", "2012-12-31"


def ric(country, tenor, protocol=PROTOCOL, ccy=CCY):
    return f"{country}GV{tenor}{ccy}{protocol}=R"


def _app_key():
    k = os.environ.get("LSEG_APP_KEY")
    if k:
        return k.strip()
    path = os.path.expanduser("~/.config/lseg/appkey")
    if os.path.exists(path):
        return open(path).read().strip()
    raise SystemExit(
        "No LSEG App Key. Set LSEG_APP_KEY, or put the key in "
        "~/.config/lseg/appkey with mode 0600. Never commit it."
    )


def fetch():
    import lseg.data as ld
    session = ld.session.desktop.Definition(app_key=_app_key()).get_session()
    ld.session.set_default(session)
    session.open()
    frames = {}
    for c in COUNTRIES:
        for t in TENORS:
            r = ric(c, t)
            try:
                h = ld.get_history(universe=r, fields=FIELDS,
                                   start=START, end=END, interval="1D")
            except Exception as e:                       # noqa: BLE001
                print(f"  {r}: {type(e).__name__} {str(e)[:70]}")
                continue
            if len(h) == 0:
                print(f"  {r}: empty")
                continue
            frames[(c, t)] = h
            p = h["PAR_MID1"].dropna()
            print(f"  {r}: {len(p)} quotes {p.index[0].date()} -> {p.index[-1].date()}")
    if not frames:
        raise SystemExit("nothing returned; is Workspace running?")
    out = pd.concat(frames, axis=1)
    out.columns = ["_".join(str(x) for x in col) for col in out.columns]
    os.makedirs(os.path.dirname(RAW), exist_ok=True)
    out.to_csv(RAW)
    print(f"wrote {RAW}  ({out.shape[0]} days x {out.shape[1]} columns)")
    return out


def load():
    """Daily panel from the cached extract, or None if it has not been pulled."""
    if not os.path.exists(RAW):
        return None
    d = pd.read_csv(RAW, index_col=0, parse_dates=True)
    return d


def quarterly(country, tenor, field="PAR_MID1"):
    """Quarterly average of the daily par spread, annualised basis points."""
    d = load()
    if d is None:
        return None
    col = f"{country}_{tenor}_{field}"
    if col not in d.columns:
        return None
    s = d[col].dropna()
    return s.groupby(pd.PeriodIndex(s.index, freq="Q")).mean()


def staleness(country, tenor):
    """Share of quoted days on which the composite does not move.

    A sovereign contract quoted by few dealers repeats yesterday's number; this
    is the cheapest available read on whether the series is a price or a mark.
    """
    d = load()
    s = d[f"{country}_{tenor}_PAR_MID1"].dropna()
    return float((s.diff() == 0).mean())


if __name__ == "__main__":
    fetch()
