# CROSS-SECTIONAL DISTRIBUTION VIA THE YOUNG (2010) LOTTERY METHOD.
import numpy as np


def get_lottery_weights(a_pol, a_grid):
    # SPLIT EACH OFF-GRID SAVINGS CHOICE ACROSS THE TWO NEAREST GRIDPOINTS.
    a_pol_c = np.clip(a_pol, a_grid[0], a_grid[-1])

    idx_hi = np.searchsorted(a_grid, a_pol_c, side="right")
    idx_hi = np.clip(idx_hi, 1, len(a_grid) - 1)
    idx_lo = idx_hi - 1

    denom = a_grid[idx_hi] - a_grid[idx_lo]
    weight_hi = np.where(denom > 0, (a_pol_c - a_grid[idx_lo]) / denom, 0.0)
    return idx_lo, idx_hi, 1.0 - weight_hi, weight_hi


def forward_iterate(D, a_pol, a_grid, Pi):
    # ONE FORWARD STEP: LOTTERY SCATTER OVER THE ASSET GRID, THEN THE MARKOV STEP.
    n_a, n_e = D.shape
    idx_lo, idx_hi, w_lo, w_hi = get_lottery_weights(a_pol, a_grid)

    # single bincount over flattened (a, e) indices: ~4x faster than per-e add.at
    cols = np.arange(n_e)
    flat = np.concatenate([(idx_lo * n_e + cols).ravel(),
                           (idx_hi * n_e + cols).ravel()])
    wts  = np.concatenate([(D * w_lo).ravel(), (D * w_hi).ravel()])
    pre  = np.bincount(flat, weights=wts, minlength=n_a * n_e).reshape(n_a, n_e)

    return pre @ Pi


def stationary_distribution(a_pol, a_grid, Pi, pi_e_stationary, tol, maxiter=100_000):
    # ITERATE THE FORWARD OPERATOR TO THE STATIONARY DISTRIBUTION OVER (a, e).
    n_a, n_e = a_pol.shape
    D = np.zeros((n_a, n_e))
    D[0, :] = pi_e_stationary   # start as a point mass at a_min

    for _ in range(maxiter):
        D_new = forward_iterate(D, a_pol, a_grid, Pi)
        if np.max(np.abs(D_new - D)) < tol:
            return D_new
        D = D_new

    raise RuntimeError(
        f"Distribution iteration did not converge (diff={np.max(np.abs(D_new - D)):.2e})")


def forward_paths(D0, a_pol_path, c_path, a_grid, Pi, use_fast=True):
    # FORWARD-SIMULATE THE DISTRIBUTION OVER THE TRANSITION (NUMBA KERNEL OR NUMPY).
    import fast_kernels

    T = a_pol_path.shape[0]
    if use_fast and fast_kernels.HAVE_NUMBA:
        return fast_kernels.dist_forward(
            np.ascontiguousarray(D0), np.ascontiguousarray(a_pol_path),
            np.ascontiguousarray(c_path), np.ascontiguousarray(a_grid),
            np.ascontiguousarray(Pi))

    # C_t on the START-of-period dist, A_t on the END-of-period one;
    # D_start[t] is the dist entering t (D_start[0]=D0, D_start[T]=terminal)
    A_path  = np.empty(T)
    C_path  = np.empty(T)
    D_start = np.empty((T + 1,) + D0.shape)
    D = D0
    for t in range(T):
        D_start[t] = D
        C_path[t]  = aggregate_consumption(D, c_path[t])
        D = forward_iterate(D, a_pol_path[t], a_grid, Pi)
        A_path[t] = aggregate_assets(D, a_grid)
    D_start[T] = D
    return A_path, C_path, D_start


def aggregate_assets(D, a_grid):
    # DISTRIBUTION-WEIGHTED MEAN ASSETS.
    return float(np.sum(D * a_grid[:, None]))


def aggregate_consumption(D, c_pol):
    # DISTRIBUTION-WEIGHTED MEAN CONSUMPTION.
    return float(np.sum(D * c_pol))
