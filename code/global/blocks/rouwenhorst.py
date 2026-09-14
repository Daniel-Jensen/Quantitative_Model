# ROUWENHORST (1995) DISCRETIZATION OF AN AR(1) INCOME PROCESS.
import numpy as np


def _rouwenhorst_Pi(rho, n):
    # RECURSIVE CONSTRUCTION OF THE n-STATE ROUWENHORST TRANSITION MATRIX.
    p = (1 + rho) / 2
    if n == 2:
        return np.array([[p, 1 - p], [1 - p, p]])

    Pi_prev = _rouwenhorst_Pi(rho, n - 1)
    Pi = np.zeros((n, n))
    Pi[:-1, :-1] += p * Pi_prev
    Pi[:-1, 1:]  += (1 - p) * Pi_prev
    Pi[1:, :-1]  += (1 - p) * Pi_prev
    Pi[1:, 1:]   += p * Pi_prev
    Pi[1:-1, :] /= 2   # interior rows double-counted above
    return Pi


def stationary_distribution(Pi):
    # STATIONARY DISTRIBUTION OF A MARKOV CHAIN VIA POWER ITERATION.
    dist = np.ones(Pi.shape[0]) / Pi.shape[0]
    for _ in range(10_000):
        new_dist = dist @ Pi
        if np.max(np.abs(new_dist - dist)) < 1e-14:
            dist = new_dist
            break
        dist = new_dist
    return dist / dist.sum()


def rouwenhorst(rho, sigma, n=2):
    # DISCRETIZE AN AR(1): RETURNS (e_grid NORMALIZED TO MEAN 1, Pi, STATIONARY DIST).
    Pi = _rouwenhorst_Pi(rho, n)
    psi = sigma * np.sqrt(n - 1) / np.sqrt(1 - rho ** 2)
    e_grid = np.exp(np.linspace(-psi, psi, n))

    pi_stationary = stationary_distribution(Pi)
    e_grid = e_grid / (pi_stationary @ e_grid)

    return e_grid, Pi, pi_stationary
