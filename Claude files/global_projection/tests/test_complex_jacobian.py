# THE COMPLEX-STEP JACOBIAN (2026-10-28): THE RESIDUAL MUST CARRY A COMPLEX PERTURBATION END TO END.
# A float() or a float array anywhere on the path drops the imaginary part with only a warning, and the
# column comes back silently zero, so the checks run with that warning turned into an error.
import sys, os, warnings
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from common import get_ss                                                                   # noqa: E402
from global_projection.solver_recursive.state_grid import s_process_params                 # noqa: E402
from global_projection.solver_recursive.decision_rules import RuleSet                       # noqa: E402
from global_projection.solver_recursive.recursive_main import calibrate_household_anchors   # noqa: E402
from global_projection.solver_recursive import collocation as C                             # noqa: E402
from global_projection.solver_recursive.recursive_experiment import normal_box              # noqa: E402


def _setup():
    # SS RULES ON THE COARSE BOX, PERTURBED OFF THE SS SO EVERY BLOCK IS ACTIVE, SMOOTHED KINKS.
    cal, ss = get_ss()
    cal = dict(cal, nw_floor_frac=0.5, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    rules = RuleSet.from_ss(normal_box(ss, cal, mu=1), ss, cal, n_regimes=2)
    rules.n_gh = 5
    F = C.make_residual(rules, cal, ss, sproc, (0, 1), False, 5)
    x = C.pack(rules, (0, 1))
    x = x + 1e-3 * np.random.default_rng(3).standard_normal(x.size)
    return F, x


def test_complex_step_matches_central_differences():
    # 30 COLUMNS: THE COMPLEX STEP AGAINST CENTRAL DIFFERENCES, AND F(x + 0j) AGAINST F(x).
    F, x = _setup()
    with warnings.catch_warnings():
        warnings.simplefilter("error", np.exceptions.ComplexWarning)
        f = F(x)
        fc = F(x.astype(complex))
        # complex exp/log/pow round differently in the last bits; only the imaginary part is used
        assert np.max(np.abs(fc.real - f)) < 1e-12 and not np.any(fc.imag), "the real part must not move"
        cols = np.random.default_rng(5).choice(x.size, 30, replace=False)
        h = 1e-6
        worst = 0.0
        for i in cols:
            cs = C._column(F, x, f, h, i, "complex")
            e = np.zeros(x.size); e[i] = h
            cd = (F(x + e) - F(x - e)) / (2 * h)
            # central differences are good to ~h^2 * F''' and round-off ~1e-16/h
            err = np.max(np.abs(cs - cd)) / max(np.max(np.abs(cd)), 1e-3)
            worst = max(worst, err)
            assert np.any(cs != 0.0) or not np.any(np.abs(cd) > 1e-8), f"column {i} came back zero"
    assert worst < 1e-6, worst
    print(f"  complex step vs central differences, 30 columns: worst relative gap {worst:.1e}: PASSED")


if __name__ == "__main__":
    test_complex_step_matches_central_differences()
