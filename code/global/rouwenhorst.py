# **Rouwenhorst (1995) discretization of an AR(1) income process.**
import numpy as np


def _rouwenhorst_Pi(rho, n):
    # **Recursive construction of the n-state Rouwenhorst transition matrix.**
    p = (1 + rho) / 2
    if n == 2:
        return np.array([[p, 1 - p], [1 - p, p]])

    Pi_prev = _rouwenhorst_Pi(rho, n - 1)
    Pi = np.zeros((n, n))
    Pi[:-1, :-1] += p * Pi_prev
    Pi[:-1, 1:] += (1 - p) * Pi_prev
    Pi[1:, :-1] += (1 - p) * Pi_prev
    Pi[1:, 1:] += p * Pi_prev
    Pi[1:-1, :] /= 2   # interior rows double-counted above
    return Pi


def stationary_distribution(Pi):
    # **Stationary distribution of a Markov chain via power iteration.**
    n = Pi.shape[0]
    dist = np.ones(n) / n
    for _ in range(10_000):
        new_dist = dist @ Pi
        if np.max(np.abs(new_dist - dist)) < 1e-14:
            dist = new_dist
            break
        dist = new_dist
    return dist / dist.sum()


def rouwenhorst(rho, sigma, n=2):
    # **Discretize an AR(1): returns (e_grid normalized to mean 1, Pi, stationary dist).**
    Pi = _rouwenhorst_Pi(rho, n)
    psi = sigma * np.sqrt(n - 1) / np.sqrt(1 - rho ** 2)
    x_grid = np.linspace(-psi, psi, n)
    e_grid = np.exp(x_grid)

    pi_stationary = stationary_distribution(Pi)
    e_grid = e_grid / (pi_stationary @ e_grid)

    return e_grid, Pi, pi_stationary
