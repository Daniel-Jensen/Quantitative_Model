# THE REST-POINT CALIBRATION LOOP ON A SYNTHETIC MOMENT MAP (NO MODEL SOLVE).
# C1  Broyden in tolerance-scaled units reaches a known root of a coupled, mildly nonlinear map
#     within the tolerances, in a handful of evaluations.
# C2  the checkpoint resumes: stopping after two evaluations and calling again finishes from
#     the saved state, with the same answer as an uninterrupted run.
# C3  a failed solve is rejected and the step halved from the last accepted point.
import sys, os, json, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
from pathlib import Path

import numpy as np

import global_projection.solver_recursive.calibrate_stochastic as CS
from calibration.global_projection import get_calibration

SETTINGS = dict(NW_FLOOR=0.15, RISK_MU_VEC=None, ROTATE_P=False, S_REFINE=5)


def _synthetic(fail_at=None):
    # A MOMENT MAP WITH A KNOWN ROOT x*: OWN SLOPES ~1, CROSS EFFECTS, AND CURVATURE.
    cal = get_calibration()
    t = CS._targets(cal)
    x_star = t + np.array([0.31, 0.12, 4.0, 1.0, 0.004, -0.002]) * np.array(
        [1, 1, 1 / 4e4, 1 / 4e4, 1, 1])
    rng = np.random.default_rng(3)
    A = np.eye(6) + 0.08 * rng.standard_normal((6, 6))
    calls = []

    def evaluate(x, settings, verbose=True):
        calls.append(np.array(x))
        u = (np.asarray(x) - x_star) / CS.TOL
        m = t + CS.TOL * (A @ u + 0.002 * u * np.abs(u))
        ok = not (fail_at is not None and len(calls) == fail_at)
        return m, dict(solve_ok=ok, solve_worst=1e-10, read_resid=1e-6, p_ss=1.15, p_rest=1.15)
    return evaluate, x_star, calls


def test_c1_converges_to_known_root():
    # C1: CONVERGENCE WITHIN TOLERANCE, FEW EVALUATIONS.
    ev, x_star, calls = _synthetic()
    CS.evaluate = ev
    with tempfile.TemporaryDirectory() as d:
        x, m = CS.calibrate(SETTINGS, Path(d) / "state.json")
    t = CS._targets(get_calibration())
    assert m is not None and np.all(np.abs(m - t) <= CS.TOL), "did not converge"
    assert len(calls) <= 8, f"{len(calls)} evaluations"
    print(f"    C1: converged in {len(calls)} evaluations, max scaled gap "
          f"{np.max(np.abs((m - t) / CS.TOL)):.2f}")


def test_c2_resumes_from_checkpoint():
    # C2: TWO EVALUATIONS, STOP, RESUME; SAME ANSWER AS ONE RUN.
    ev, _, calls = _synthetic()
    CS.evaluate = ev
    with tempfile.TemporaryDirectory() as d:
        x_full, _ = CS.calibrate(SETTINGS, Path(d) / "a.json")
        n_full = len(calls)
        calls.clear()
        CS.calibrate(SETTINGS, Path(d) / "b.json", max_evals=2)
        assert len(json.loads((Path(d) / "b.json").read_text())["history"]) == 2
        x_res, _ = CS.calibrate(SETTINGS, Path(d) / "b.json")
    assert np.array_equal(x_full, x_res), "resumed run differs from the uninterrupted one"
    assert len(calls) == n_full, "the resumed run repeated or skipped an evaluation"
    print(f"    C2: resumed after 2 of {n_full} evaluations, identical instruments")


def test_c3_rejects_failed_solve():
    # C3: A FAILED SOLVE HALVES THE STEP FROM THE LAST ACCEPTED POINT AND IS NOT A MOMENT.
    ev, _, calls = _synthetic(fail_at=2)
    CS.evaluate = ev
    with tempfile.TemporaryDirectory() as d:
        x, m = CS.calibrate(SETTINGS, Path(d) / "c.json")
        hist = json.loads((Path(d) / "c.json").read_text())["history"]
    assert hist[1].get("rejected"), "the failed solve was accepted"
    assert np.allclose(np.array(hist[1]["x_next"]),
                       0.5 * (np.array(hist[0]["x"]) + np.array(hist[1]["x"])))
    assert m is not None, "did not converge after the rejection"
    print(f"    C3: failed solve rejected and halved; converged in {len(calls)} evaluations")


if __name__ == "__main__":
    test_c1_converges_to_known_root()
    print("test_calibrate_stochastic C1 (convergence): PASSED")
    test_c2_resumes_from_checkpoint()
    print("test_calibrate_stochastic C2 (resume): PASSED")
    test_c3_rejects_failed_solve()
    print("test_calibrate_stochastic C3 (rejection): PASSED")
