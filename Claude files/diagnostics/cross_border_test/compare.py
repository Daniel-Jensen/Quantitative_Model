# CROSS-BORDER TEST (2026-10-30): DO FOREIGN (GERMAN) BANKS BUY OR SELL GREEK GOVERNMENT DEBT AS THE GREEK SPREAD RISES?
# DATA: BIS consolidated banking statistics (claims on the Greek official sector, ultimate-risk basis, USD m -> EUR at the
# ECB end-of-quarter rate), ECB BSI (Greek MFIs' holdings of euro-area government debt securities, EUR m, mostly GGBs),
# FRED 10y yields (Empirics/outputs/spreads_fred.csv). Raw downloads in raw/.
import csv
from pathlib import Path
HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
QUARTERS = ("2009-Q4", "2010-Q1", "2010-Q2", "2010-Q3", "2010-Q4", "2011-Q1", "2011-Q2", "2011-Q3", "2011-Q4")


def bis(country):
    # BIS CLAIMS ON THE GREEK OFFICIAL SECTOR BY REPORTING COUNTRY, USD m.
    return {r["TIME_PERIOD"]: float(r["OBS_VALUE"]) for r in csv.DictReader(open(RAW / "bis_official_GR.csv"))
            if r["L_REP_CTY"] == country and r["OBS_VALUE"] not in ("", "NaN")}


def main():
    fx = {r["TIME_PERIOD"]: float(r["OBS_VALUE"]) for r in csv.DictReader(open(RAW / "ecb_usd_eur_eop.csv"))}
    gr = {r["TIME_PERIOD"]: float(r["OBS_VALUE"]) for r in csv.DictReader(open(RAW / "bsi_gr_U2.csv"))}
    sp = {r["date"][:7]: float(r["GRC_spread"]) for r in csv.DictReader(open(HERE.parents[2] / "Empirics/outputs/spreads_fred.csv"))
          if r["GRC_spread"]}
    de, fo = bis("DE"), bis("5A")
    q2m = lambda q: f"{q[:4]}-{int(q[-1]) * 3:02d}"
    print(f"{'quarter':8s} {'GR spread bp':>12s} {'DE banks EURbn':>14s} {'all foreign EURbn':>17s} {'GR banks EURbn':>14s} "
          f"{'DE share %':>10s} {'foreign share %':>15s}")
    for q in QUARTERS:
        d, f, g = de[q] / fx[q] / 1e3, fo[q] / fx[q] / 1e3, gr[q2m(q)] / 1e3
        print(f"{q:8s} {100 * sp[q2m(q)]:12.0f} {d:14.1f} {f:17.1f} {g:14.1f} {100 * d / (d + g):10.1f} {100 * f / (f + g):15.1f}")


if __name__ == "__main__":
    main()
