"""Non-stochastic simulation of the cross-sectional distribution via the
'lottery' method (Young, 2010): savings choices that fall between grid
points are split across the two nearest gridpoints with probability
weights chosen to match the mean exactly.
"""
import numpy as np


def stationary_distribution(a_pol, a_grid, Pi, pi_e_stationary, tol, maxiter=100_000):
    "Compute the stationary distribution of assets and income given policy functions"
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
    "Gets the output from the get_lottery_weights function and then uses np.add.at to update the distribution D."
    n_a, n_e = D.shape
    idx_lo, idx_hi, w_lo, w_hi = get_lottery_weights(a_pol, a_grid)

    pre = np.zeros((n_a, n_e))
    for e in range(n_e):
        np.add.at(pre[:, e], idx_lo[:, e], D[:, e] * w_lo[:, e])
        np.add.at(pre[:, e], idx_hi[:, e], D[:, e] * w_hi[:, e])

    return pre @ Pi


def get_lottery_weights(a_pol, a_grid):
    "Young (2010) lottery method"
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
