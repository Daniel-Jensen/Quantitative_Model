# FAST-KERNEL EQUIVALENCE: THE NUMBA KERNELS MUST REPRODUCE THE PURE-NUMPY
# REFERENCE IMPLEMENTATIONS TO FLOAT PRECISION, AND THE BINCOUNT SCATTER MUST
# MATCH THE HISTORICAL np.add.at SCATTER.
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402

from blocks import fast_kernels  # noqa: E402
from blocks.household import solve_backward_transition  # noqa: E402
from blocks.distribution import (forward_iterate, forward_paths,  # noqa: E402
                          get_lottery_weights)


def _random_problem(seed=0, T=40, n_a=80, n_e=3):
    # RANDOM HOUSEHOLD PROBLEM EXERCISING BOTH THE INTERIOR AND CONSTRAINED BRANCHES.
    rng = np.random.default_rng(seed)
    a_grid = np.linspace(0.0, 1.0, n_a) ** 2 * 60.0
    Pi = rng.uniform(0.05, 1.0, (n_e, n_e))
    Pi /= Pi.sum(axis=1, keepdims=True)
    r_path = 0.005 + 0.002 * rng.standard_normal(T + 1)
    # include a low-income state so the borrowing constraint binds at the
    # bottom of the grid (exercises the constrained branch)
    y_path = 0.15 + 0.9 * rng.uniform(size=(T, n_e))
    vN_path = 0.30 + 0.05 * rng.uniform(size=T)
    beta, sigma, a_min = 0.97, 2.0, 0.0
    c_term = np.maximum(0.02 * a_grid[:, None] + y_path[-1][None, :], 1e-8)
    return dict(a_grid=a_grid, Pi=Pi, r_path=r_path, y_path=y_path,
                vN_path=vN_path, beta=beta, sigma=sigma, a_min=a_min,
                c_term=c_term)


def test_hh_backward_equivalence():
    # THE EGM BACKWARD KERNEL MUST MATCH THE NUMPY REFERENCE.
    if not fast_kernels.HAVE_NUMBA:
        print("  numba not importable — hh_backward equivalence skipped")
        return
    for seed in (0, 1, 2):
        p = _random_problem(seed)
        args = (p["a_grid"], p["Pi"], p["r_path"], p["y_path"], p["c_term"],
                p["beta"], p["sigma"], p["a_min"])
        c_np, a_np = solve_backward_transition(*args, vN_path=p["vN_path"],
                                               use_fast=False)
        c_nb, a_nb = solve_backward_transition(*args, vN_path=p["vN_path"],
                                               use_fast=True)
        assert np.allclose(c_np, c_nb, rtol=0, atol=1e-12), \
            f"c mismatch (seed {seed}): {np.max(np.abs(c_np - c_nb)):.2e}"
        assert np.allclose(a_np, a_nb, rtol=0, atol=1e-12), \
            f"a_pol mismatch (seed {seed}): {np.max(np.abs(a_np - a_nb)):.2e}"


def test_dist_forward_equivalence():
    # THE DISTRIBUTION FORWARD KERNEL MUST MATCH THE NUMPY REFERENCE.
    if not fast_kernels.HAVE_NUMBA:
        print("  numba not importable — dist_forward equivalence skipped")
        return
    for seed in (0, 1):
        p = _random_problem(seed)
        c_np, a_np = solve_backward_transition(
            p["a_grid"], p["Pi"], p["r_path"], p["y_path"], p["c_term"],
            p["beta"], p["sigma"], p["a_min"], vN_path=p["vN_path"],
            use_fast=False)
        rng = np.random.default_rng(seed + 100)
        n_a, n_e = p["a_grid"].size, p["Pi"].shape[0]
        D0 = rng.uniform(size=(n_a, n_e))
        D0 /= D0.sum()
        ref = forward_paths(D0, a_np, c_np, p["a_grid"], p["Pi"],
                            use_fast=False)
        fast = forward_paths(D0, a_np, c_np, p["a_grid"], p["Pi"],
                             use_fast=True)
        for name, x_ref, x_fast in zip(("A", "C", "D_start"), ref, fast):
            assert np.allclose(x_ref, x_fast, rtol=0, atol=1e-12), \
                (f"{name} mismatch (seed {seed}): "
                 f"{np.max(np.abs(np.asarray(x_ref) - np.asarray(x_fast))):.2e}")


def test_bincount_scatter_matches_add_at():
    # forward_iterate'S ONE-SHOT BINCOUNT MUST EQUAL THE HISTORICAL PER-e
    # np.add.at SCATTER, UP TO FLOAT ASSOCIATIVITY.
    rng = np.random.default_rng(7)
    n_a, n_e = 120, 4
    a_grid = np.linspace(0.0, 1.0, n_a) ** 2 * 40.0
    Pi = rng.uniform(0.05, 1.0, (n_e, n_e))
    Pi /= Pi.sum(axis=1, keepdims=True)
    D = rng.uniform(size=(n_a, n_e))
    D /= D.sum()
    a_pol = rng.uniform(-1.0, 45.0, size=(n_a, n_e))   # includes off-grid

    new = forward_iterate(D, a_pol, a_grid, Pi)

    idx_lo, idx_hi, w_lo, w_hi = get_lottery_weights(a_pol, a_grid)
    pre = np.zeros((n_a, n_e))
    for e in range(n_e):
        np.add.at(pre[:, e], idx_lo[:, e], D[:, e] * w_lo[:, e])
        np.add.at(pre[:, e], idx_hi[:, e], D[:, e] * w_hi[:, e])
    ref = pre @ Pi

    assert np.allclose(ref, new, rtol=0, atol=1e-15), \
        f"scatter mismatch: {np.max(np.abs(ref - new)):.2e}"
    assert abs(new.sum() - 1.0) < 1e-12, "mass not conserved"


if __name__ == "__main__":
    test_bincount_scatter_matches_add_at()
    test_hh_backward_equivalence()
    test_dist_forward_equivalence()
    print("test_fast_kernels: ALL PASSED")
