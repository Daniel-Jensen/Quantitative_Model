import numpy as np


def calculate_gini(values, weights):
    """
    Gini coefficient given values and their corresponding weights/probabilities.
    """
    values_flat  = values.flatten()
    weights_flat = weights.flatten()

    mask = weights_flat > 0
    values_flat  = values_flat[mask]
    weights_flat = weights_flat[mask]

    if len(values_flat) == 0:
        return 0.0

    sorted_indices        = np.argsort(values_flat)
    values_sorted         = values_flat[sorted_indices]
    weights_sorted        = weights_flat[sorted_indices]
    weights_sorted        = weights_sorted / weights_sorted.sum()

    cumulative_weights    = np.cumsum(weights_sorted)
    cumulative_values     = np.cumsum(values_sorted * weights_sorted)

    total_value = cumulative_values[-1]
    if total_value <= 0:
        return 0.0

    cumulative_value_share = cumulative_values / total_value
    lorenz_area            = np.trapz(cumulative_value_share, cumulative_weights)
    return 1 - 2 * lorenz_area

def hh_inner_deposit(D, dep, c):
    """Gini coefficients for the 1-asset (deposit-only) household block."""
    gini_dep = np.full_like(D, calculate_gini(dep, D))
    gini_c = np.full_like(D, calculate_gini(c, D))
    return gini_dep, gini_c


def hh_inner2(D, a, b, c):
    wealth      = a + b
    #gini_a      = np.full_like(D, calculate_gini(a,      D))
    #gini_c      = np.full_like(D, calculate_gini(c,      D))
    #gini_wealth = np.full_like(D, calculate_gini(wealth, D))
    return gini_a , gini_b, gini_c, gini_wealth




def hh_inner_deposit3(D, dep, c):
    wealth      = dep + c
    gini_dep    = np.full_like(D, calculate_gini(dep,    D))
    gini_c      = np.full_like(D, calculate_gini(c,      D))
    gini_wealth = np.full_like(D, calculate_gini(wealth, D))
    return gini_dep, gini_c, gini_wealth