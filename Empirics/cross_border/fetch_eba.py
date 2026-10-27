"""Acquire the four EBA exposure workbooks by scraping the publication pages.

Purpose (paper role)
--------------------
Layer 2 acquisition for empirical component #1. The EBA exercises are the only
public source giving bank x sovereign x residual-maturity exposures through the
Greek crisis, and they are published as spreadsheets on pages the EBA has
reorganised more than once.

The scrape is deliberately permissive: it collects every spreadsheet link on a
candidate page and ranks them by filename keywords, because the exact file names
have changed. It verifies nothing about CONTENT - eba_panel.validate_panel and
the adapter acceptance tests are what reject a wrong file - but it records a
SHA-256 of every download so provenance names exactly which bytes produced the
results.

Falls back to a clear instruction when a page yields nothing. These are four
one-off historical files; hand-downloading them is a perfectly good outcome.

Run:
    /opt/anaconda3/envs/ssj/bin/python -m Empirics.cross_border.fetch_eba

Outputs (Empirics/cross_border/data_raw/):
    eba_2010_st.xls*, eba_2011_st.xls*, eba_2011_ce.xls*, eba_2012_ce.xls*
    eba_downloads.json   url, filename, sha256 and retrieval date per file
"""

from __future__ import annotations

import glob
import hashlib
import json
import os
import re
import time
import urllib.parse
import urllib.request

from Empirics.cross_border import config

UA = "Mozilla/5.0 (compatible; academic-research/1.0)"
MANIFEST = os.path.join(config.RAW, "eba_downloads.json")

# stem -> (candidate landing pages, filename keywords worth points)
SOURCES = {
    "eba_2011_st": (
        ["https://www.eba.europa.eu/risk-and-data/risk-analysis/eu-wide-stress-testing/2011",
         "https://www.eba.europa.eu/risk-analysis-and-data/eu-wide-stress-testing/2011"],
        ("sovereign", "exposure"),
    ),
    "eba_2011_ce": (
        ["https://www.eba.europa.eu/risk-and-data/risk-analysis/eu-wide-recapitalisation",
         "https://www.eba.europa.eu/eu-capital-exercise-final-results"],
        ("sovereign", "exposure", "capital"),
    ),
    "eba_2010_st": (
        ["https://www.eba.europa.eu/risk-and-data/risk-analysis/eu-wide-stress-testing/2010",
         "https://www.eba.europa.eu/risk-analysis-and-data/eu-wide-stress-testing/2010"],
        ("sovereign", "exposure"),
    ),
    "eba_2012_ce": (
        ["https://www.eba.europa.eu/eu-capital-exercise-final-results",
         "https://www.eba.europa.eu/risk-and-data/risk-analysis/eu-wide-recapitalisation"],
        ("sovereign", "exposure", "capital"),
    ),
}

_LINK = re.compile(r"""href=["']([^"']+\.(?:xls|xlsx|zip))["']""", re.I)


def _get(url, binary=False):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120) as r:
        blob = r.read()
    return blob if binary else blob.decode("utf-8", errors="replace")


def _candidates(page_url):
    """Every spreadsheet link on a page, as absolute URLs."""
    try:
        html = _get(page_url)
    except Exception as exc:
        print("    page failed: %s (%s)" % (page_url, exc))
        return []
    return [urllib.parse.urljoin(page_url, h) for h in _LINK.findall(html)]


def _score(url, keywords):
    name = urllib.parse.unquote(url).lower()
    return sum(1 for k in keywords if k in name)


def fetch_one(stem, pages, keywords):
    for page in pages:
        print("  scanning %s" % page)
        links = _candidates(page)
        if not links:
            continue
        best = max(links, key=lambda u: _score(u, keywords))
        if _score(best, keywords) == 0:
            print("    %d spreadsheets found, none matching %s" % (len(links), keywords))
            continue
        ext = os.path.splitext(urllib.parse.urlparse(best).path)[1] or ".xlsx"
        dest = os.path.join(config.RAW, stem + ext)
        print("    downloading %s" % best)
        blob = _get(best, binary=True)
        with open(dest, "wb") as fh:
            fh.write(blob)
        time.sleep(1.0)   # be polite to a public regulator's web server
        return {
            "stem": stem,
            "url": best,
            "path": dest,
            "sha256": hashlib.sha256(blob).hexdigest(),
            "bytes": len(blob),
            "retrieved": time.strftime("%Y-%m-%d"),
        }
    return None


def main():
    os.makedirs(config.RAW, exist_ok=True)
    got, missing = [], []
    for stem, (pages, keywords) in SOURCES.items():
        if glob.glob(os.path.join(config.RAW, stem + ".xls*")):
            print("%s: already present, skipping" % stem)
            continue
        print("%s:" % stem)
        rec = fetch_one(stem, pages, keywords)
        if rec:
            got.append(rec)
        else:
            missing.append(stem)

    if got:
        existing = json.load(open(MANIFEST)) if os.path.exists(MANIFEST) else []
        json.dump(existing + got, open(MANIFEST, "w"), indent=2)
        print("\nwrote %s (+%d records)" % (MANIFEST, len(got)))

    if missing:
        print("\nNOT ACQUIRED: %s" % ", ".join(missing))
        print("Download these by hand from")
        print("  https://www.eba.europa.eu/risk-and-data/risk-analysis/eu-wide-stress-testing")
        print("and save them under those stems in %s" % config.RAW)
        print("(see data_raw/README.md for what to take from each).")


if __name__ == "__main__":
    main()
