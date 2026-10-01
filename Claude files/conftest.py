"""pytest set-up for the tests kept in Claude files/: the repo root (package imports such as
calibration.*, linear_ssj.*), experiments/ and diagnostics/regimes/ (their flat module
imports) go on sys.path.

The global_projection tests are scripts (python3 "Claude files/global_projection/tests/X.py")
and are skipped here: their helper module common.py would shadow experiments/common.py.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (os.path.join(ROOT, "diagnostics", "regimes"), os.path.join(ROOT, "experiments"), ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)

collect_ignore_glob = ["global_projection/*"]
