# SMOLYAK LIBRARY GATE: COUNTS, ON-GRID EXACTNESS, QUADRATIC EXACTNESS,
# LEVEL-2->3 ERROR DECAY, ANISOTROPY, AND A BROCK-MIRMAN TIME-ITERATION SMOKE
# TEST (LOG-POLICY IS LINEAR IN LOG STATES, SO THE GRID MUST NAIL IT).
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from smolyak import SmolyakGrid


def test_counts():
    # KNOWN SMOLYAK POINT COUNTS: d=2 mu=2 -> 13; d=6 mu=2 -> 85; d=6 mu=3 -> 389.
    assert SmolyakGrid(-np.ones(2), np.ones(2), mu=2).n == 13
    assert SmolyakGrid(-np.ones(6), np.ones(6), mu=2).n == 85
    assert SmolyakGrid(-np.ones(6), np.ones(6), mu=3).n == 389


def test_on_grid_exactness():
    # THE INTERPOLANT REPRODUCES ANY FUNCTION EXACTLY AT ITS OWN NODES.
    rng = np.random.default_rng(0)
    g = SmolyakGrid(np.array([0.5, -2.0, 1.0]), np.array([1.5, 3.0, 4.0]), mu=2)
    vals = rng.standard_normal(g.n)
    coef = g.fit(vals)
    assert np.max(np.abs(g.eval(coef, g.points) - vals)) < 1e-10


def test_quadratic_exactness():
    # LEVEL 2 SPANS ALL COMPLETE QUADRATICS INCL. CROSS TERMS (KK04 PROPERTY).
    rng = np.random.default_rng(1)
    d = 6
    lo = rng.uniform(-2.0, -0.5, d)
    hi = rng.uniform(0.5, 2.0, d)
    g = SmolyakGrid(lo, hi, mu=2)
    A = rng.standard_normal((d, d)); A = 0.5 * (A + A.T)
    b = rng.standard_normal(d); c0 = rng.standard_normal()

    def f(x):
        # RANDOM COMPLETE QUADRATIC WITH ALL CROSS TERMS.
        return np.einsum("ij,jk,ik->i", x, A, x) + x @ b + c0

    coef = g.fit(f(g.points))
    x_test = rng.uniform(lo, hi, size=(1000, d))
    err = np.max(np.abs(g.eval(coef, x_test) - f(x_test)))
    assert err < 1e-9, f"quadratic exactness violated: {err:.2e}"


def test_error_decay():
    # SUSTAINED GEOMETRIC (SPECTRAL) DECAY OF OFF-GRID ERROR AS LEVEL RISES.
    # A smooth 6D exponential: measured max errors run ~1e-1, ~1e-2, ~5e-4 at
    # mu=2,3,4 -- each level ~10-20x better. Asserting the DECAY RATIO is the
    # honest test of convergence; a fixed absolute level just encodes the
    # function's steepness, not the grid's correctness.
    rng = np.random.default_rng(2)
    d = 6
    lo, hi = -np.ones(d), np.ones(d)
    a = rng.uniform(-0.4, 0.4, d)

    def f(x):
        # SMOOTH ANISOTROPIC EXPONENTIAL TEST FUNCTION.
        return np.exp(x @ a)

    x_test = rng.uniform(lo, hi, size=(2000, d))
    errs = []
    for mu in (2, 3, 4):
        g = SmolyakGrid(lo, hi, mu=mu)
        coef = g.fit(f(g.points))
        errs.append(np.max(np.abs(g.eval(coef, x_test) - f(x_test))))
    assert errs[1] < 0.25 * errs[0], f"weak decay 2->3: {errs[0]:.2e} -> {errs[1]:.2e}"
    assert errs[2] < 0.25 * errs[1], f"weak decay 3->4: {errs[1]:.2e} -> {errs[2]:.2e}"
    assert errs[2] < 1e-3, f"level-4 error too large: {errs[2]:.2e}"


def test_anisotropic():
    # PER-DIMENSION LEVEL CAPS GIVE A VALID SQUARE SYSTEM, ON-GRID EXACT.
    rng = np.random.default_rng(3)
    g = SmolyakGrid(-np.ones(3), np.ones(3), mu=3, mu_vec=[3, 2, 1])
    vals = rng.standard_normal(g.n)
    coef = g.fit(vals)
    assert np.max(np.abs(g.eval(coef, g.points) - vals)) < 1e-10
    iso = SmolyakGrid(-np.ones(3), np.ones(3), mu=3)
    assert g.n < iso.n   # caps must prune points


def test_brock_mirman_time_iteration():
    # FULL-LOOP SMOKE: TIME ITERATION ON THE SAME INFRASTRUCTURE RECOVERS THE
    # ANALYTIC BROCK-MIRMAN POLICY C = (1-alpha*beta)*z*k^alpha (log-linear in
    # log states, so a level-2 grid on (log k, log z) is exact up to solver tol).
    alpha, beta = 0.36, 0.96
    rho, sig = 0.9, 0.02
    # ergodic box for (log k, log z); BM steady state log k* = log(alpha*beta)/(1-alpha)
    logk_ss = np.log(alpha * beta) / (1.0 - alpha)
    lo = np.array([logk_ss - 0.5, -4.0 * sig / np.sqrt(1 - rho ** 2)])
    hi = np.array([logk_ss + 0.5, 4.0 * sig / np.sqrt(1 - rho ** 2)])
    g = SmolyakGrid(lo, hi, mu=2)

    gh_x, gh_w = np.polynomial.hermite_e.hermegauss(5)   # E[f(eps)], eps ~ N(0,1)
    gh_w = gh_w / gh_w.sum()

    def c_analytic(logk, logz):
        # CLOSED-FORM BROCK-MIRMAN CONSUMPTION.
        return (1.0 - alpha * beta) * np.exp(logz + alpha * logk)

    coef = g.fit(np.log(c_analytic(g.points[:, 0], g.points[:, 1]) * 1.35))  # start off the truth
    for it in range(600):
        logk, logz = g.points[:, 0], g.points[:, 1]
        y = np.exp(logz + alpha * logk)
        c_prev = np.exp(g.eval(coef, g.points))
        c_new = np.empty(g.n)
        for p in range(g.n):
            # solve the euler for c at this node given next period's rule
            def euler(c):
                # UNIT-FREE EULER RESIDUAL AT CANDIDATE CONSUMPTION c.
                kp = y[p] - c
                if kp <= 0:
                    return 1e6
                logkp = np.log(kp)
                logzp = rho * logz[p] + sig * gh_x
                cp = np.exp(g.eval(coef, np.column_stack(
                    [np.full(5, logkp), logzp])))
                rhs = beta * np.sum(gh_w * alpha * np.exp(logzp)
                                    * kp ** (alpha - 1.0) / cp)
                return 1.0 / c - rhs
            c_lo, c_hi = 1e-8, y[p] * (1 - 1e-8)
            from scipy.optimize import brentq
            c_new[p] = brentq(euler, c_lo, c_hi, xtol=1e-14)
        coef_new = g.fit(np.log(c_new))
        step = np.max(np.abs(coef_new - coef))
        coef = 0.5 * coef + 0.5 * coef_new
        if step < 1e-12:
            break

    x_test = np.random.default_rng(4).uniform(lo, hi, size=(500, 2))
    c_num = np.exp(g.eval(coef, x_test))
    c_true = c_analytic(x_test[:, 0], x_test[:, 1])
    rel = np.max(np.abs(c_num / c_true - 1.0))
    assert rel < 1e-6, f"Brock-Mirman policy error {rel:.2e}"


if __name__ == "__main__":
    test_counts()
    test_on_grid_exactness()
    test_quadratic_exactness()
    test_error_decay()
    test_anisotropic()
    test_brock_mirman_time_iteration()
    print("test_smolyak: ALL PASSED")
