# M2 (2026-10-26): THE REST-POINT CALIBRATION LOOP AT THE APPENDIX-B TARGETS, ON THE COARSE GRID.
#   python3 "Claude files/diagnostics/tex_calibration/calibrate_srp.py" [--start runs/m2_start.json]
# The checkpoint runs/m2_state.json resumes an interrupted loop; --start only sets the first instruments.
import argparse
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default=None)
    ap.add_argument("--max-evals", type=int, default=12)
    a = ap.parse_args()

    import global_projection.solver_recursive.calibrate_stochastic as CS
    if a.start:
        start = json.loads((HERE / a.start).read_text())
        base = CS.get_calibration

        def with_start():
            # THE CALIBRATION WITH THE LOOP'S FIRST INSTRUMENTS (evaluate overrides them each step).
            cal = base()
            cal.update(start)
            return cal
        CS.get_calibration = with_start
    settings = dict(NW_FLOOR=0.15, RISK_MU_VEC=None, ROTATE_P=False, S_REFINE=0)
    x, m = CS.calibrate(settings, HERE / "runs" / "m2_state.json", max_evals=a.max_evals)
    if m is not None:
        (HERE / "runs" / "m2_instruments.json").write_text(
            json.dumps({p[0]: float(v) for p, v in zip(CS.PAIRS, x)}, indent=1))


if __name__ == "__main__":
    main()
