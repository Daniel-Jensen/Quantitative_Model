# ACCURACY AND ECONOMIC DIAGNOSTICS FOR A SOLVED RULE SET.
# Unit-free log10 Euler errors at collocation points and at off-grid draws;
# constraint binding frequency. The Euler residuals are already unit-free
# (deposit Euler ~ R*E[Lam]-1; the two asset Eulers minus lam*mu), so the
# log10 of their max magnitude is the natural accuracy report.
import numpy as np

from time_iteration import point_block


def euler_errors(grid, coef, cal, quad, points):
    # MAX-ABS EULER RESIDUAL AT EACH SUPPLIED STATE (SOLVING av/mu INTERNALLY).
    # Re-solves the pointwise rule values implied by the coefficients (i.e.
    # evaluates the rules, then measures how far the Eulers miss).
    errs = np.empty(points.shape[0])
    for i, s in enumerate(points):
        C = grid.eval(coef["C"], s)[0]
        R = grid.eval(coef["R"], s)[0]
        QB = grid.eval(coef["QB"], s)[0]
        res, _, _, _ = point_block([C, R, QB], s, grid, coef, cal, quad)
        errs[i] = np.max(np.abs(res))
    return errs


def _stat(e):
    # MEAN, MEDIAN, MAX, AND 99th-PCTILE LOG10 ERROR (GUARDING log10(0)).
    le = np.log10(np.maximum(e, 1e-16))
    return dict(mean=le.mean(), med=np.median(le), p99=np.quantile(le, 0.99),
                max=le.max())


def accuracy_report(grid, coef, cal, quad, n_off=2000, seed=0, path=None):
    # LOG10 EULER ERRORS ON-GRID, AT UNIFORM OFF-GRID DRAWS, AND (IF SUPPLIED)
    # ALONG THE ERGODIC PATH -- the last is the accuracy that actually matters,
    # since the model never visits the box corners a uniform draw penalizes.
    rng = np.random.default_rng(seed)
    off_pts = grid.sample_interior(n_off, rng)   # inside the grid domain
    rep = dict(on=_stat(euler_errors(grid, coef, cal, quad, grid.points)),
               off=_stat(euler_errors(grid, coef, cal, quad, off_pts)))
    if path is not None:
        idx = np.linspace(0, len(path) - 1, min(2000, len(path))).astype(int)
        rep["ergodic"] = _stat(euler_errors(grid, coef, cal, quad, path[idx]))
    return rep


def print_accuracy(rep):
    # HUMAN-READABLE ACCURACY TABLE (LOG10 EULER ERRORS).
    print("Euler-equation accuracy  (log10 max-abs residual)")
    for key in ("on", "off", "ergodic"):
        if key in rep:
            s = rep[key]
            label = {"on": "on-grid ", "off": "off-grid", "ergodic": "ergodic "}[key]
            print(f"  {label} : mean {s['mean']:+.2f}  med {s['med']:+.2f}  "
                  f"p99 {s['p99']:+.2f}  max {s['max']:+.2f}")
