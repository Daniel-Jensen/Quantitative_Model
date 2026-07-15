"""Non-stochastic simulation of the cross-sectional distribution via the
'lottery' method (Young, 2010): savings choices that fall between grid
points are split across the two nearest gridpoints with probability
weights chosen to match the mean exactly.
"""
import numpy as np


def stationary_distribution(a_pol, a_grid, Pi, pi_e_stationary, tol, maxiter=100_000):
    # COMPUTE STATIONARY DISTRIBUTION OF ASSETS AND INCOME GIVEN POLICY FUNCTIONS
    n_a, n_e = a_pol.shape

    # This is the propability distribution over (a, e). 
    # D[i, j] = the fraction of the population (mass, between 0 and 1) currently holding asset level a_grid[i] and being in income state e_grid[j].
    D = np.zeros((n_a, n_e))
    D[0, :] = pi_e_stationary  # start as point mass at a_min

    for it in range(maxiter):
        #The distribution fitted to the grid
        D_new = forward_iterate(D, a_pol, a_grid, Pi)
        if np.max(np.abs(D_new - D)) < tol:
            return D_new
        D = D_new

    raise RuntimeError(f"Distribution iteration did not converge (diff={np.max(np.abs(D_new - D)):.2e})")


def forward_iterate(D, a_pol, a_grid, Pi):
    # LOTTERY-WEIGHT SCATTER OF THE DISTRIBUTION, THEN THE MARKOV STEP.
    # Both legs go through ONE np.bincount over flattened (a, e) indices —
    # same sums as the historical per-e np.add.at, ~4x faster.
    n_a, n_e = D.shape
    idx_lo, idx_hi, w_lo, w_hi = get_lottery_weights(a_pol, a_grid)

    cols = np.arange(n_e)
    flat = np.concatenate([(idx_lo * n_e + cols).ravel(),
                           (idx_hi * n_e + cols).ravel()])
    wts  = np.concatenate([(D * w_lo).ravel(), (D * w_hi).ravel()])
    pre = np.bincount(flat, weights=wts, minlength=n_a * n_e).reshape(n_a, n_e)

    return pre @ Pi


def forward_paths(D0, a_pol_path, c_path, a_grid, Pi, use_fast=True):
    """Forward-simulate the distribution over the whole transition.

    Returns (A_path, C_path, D_start) with the exact timing convention of the
    transition solver: C_t aggregates over the START-of-period distribution,
    A_t over the END-of-period one; D_start[t] is the distribution entering
    period t (D_start[0] = D0, D_start[T] = terminal).

    Dispatches to the numba kernel (fast_kernels.py) when importable and
    `use_fast` (cal["use_numba"]); the pure-numpy loop below is the reference
    implementation.
    """
    import fast_kernels
    T = a_pol_path.shape[0]
    if use_fast and fast_kernels.HAVE_NUMBA:
        return fast_kernels.dist_forward(
            np.ascontiguousarray(D0), np.ascontiguousarray(a_pol_path),
            np.ascontiguousarray(c_path), np.ascontiguousarray(a_grid),
            np.ascontiguousarray(Pi))

    A_path = np.empty(T)
    C_path = np.empty(T)
    D_start = np.empty((T + 1,) + D0.shape)
    D = D0
    for t in range(T):
        D_start[t] = D
        C_path[t] = aggregate_consumption(D, c_path[t])
        D = forward_iterate(D, a_pol_path[t], a_grid, Pi)
        A_path[t] = aggregate_assets(D, a_grid)
    D_start[T] = D
    return A_path, C_path, D_start


def get_lottery_weights(a_pol, a_grid):
    # YOUNG (2010) LOTTERY METHOD
    a_min, a_max = a_grid[0], a_grid[-1]
    a_pol_c = np.clip(a_pol, a_min, a_max)

    # If the policy function a' is between two grid points, we split the mass between the two nearest grid points.
    idx_hi = np.searchsorted(a_grid, a_pol_c, side="right")
    idx_hi = np.clip(idx_hi, 1, len(a_grid) - 1)
    idx_lo = idx_hi - 1

    # Computes those weights
    denom = a_grid[idx_hi] - a_grid[idx_lo]
    weight_hi = np.where(denom > 0, (a_pol_c - a_grid[idx_lo]) / denom, 0.0)
    weight_lo = 1.0 - weight_hi
    return idx_lo, idx_hi, weight_lo, weight_hi


# Calculate the aggregate assets 
def aggregate_assets(D, a_grid):
    return float(np.sum(D * a_grid[:, None]))

# Calculate the aggregate Consumtpion 
def aggregate_consumption(D, c_pol):
    return float(np.sum(D * c_pol))
