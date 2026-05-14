import numpy as np
import matplotlib.pyplot as plt

def get_lorenz_curve(values, weights):
    """
    Calculates the cumulative population shares and value shares for a Lorenz curve.
    """
    values_flat  = values.flatten()
    weights_flat = weights.flatten()

    mask = weights_flat > 0
    values_flat  = values_flat[mask]
    weights_flat = weights_flat[mask]

    if len(values_flat) == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0])

    sorted_indices = np.argsort(values_flat)
    values_sorted  = values_flat[sorted_indices]
    weights_sorted = weights_flat[sorted_indices]
    weights_sorted = weights_sorted / weights_sorted.sum()

    cumulative_weights = np.cumsum(weights_sorted)
    cumulative_values  = np.cumsum(values_sorted * weights_sorted)

    total_value = cumulative_values[-1]
    if total_value <= 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0])

    cumulative_value_share = cumulative_values / total_value

    # Insert (0,0) at the beginning so the plot starts smoothly from the origin
    cumulative_weights = np.insert(cumulative_weights, 0, 0.0)
    cumulative_value_share = np.insert(cumulative_value_share, 0, 0.0)

    return cumulative_weights, cumulative_value_share


def calculate_gini(values, weights):
    """
    Gini coefficient given values and their corresponding weights/probabilities.
    """
    cum_weights, cum_value_shares = get_lorenz_curve(values, weights)
    
    # Failsafe for empty arrays or edge cases
    if len(cum_weights) <= 2: 
        return 0.0

    lorenz_area = np.trapz(cum_value_shares, cum_weights)
    return 1 - 2 * lorenz_area


def plot_consumption_and_wealth_gini(D, c, wealth, title="Lorenz Curves: Consumption vs Wealth"):
    """
    Plots the Lorenz curves for Consumption and Wealth, 
    and displays their Gini coefficients in the legend.
    """
    # Fetch coordinates
    c_weights, c_shares = get_lorenz_curve(c, D)
    w_weights, w_shares = get_lorenz_curve(wealth, D)
    
    # Calculate Gini scalars for the legend
    gini_c = calculate_gini(c, D)
    gini_w = calculate_gini(wealth, D)
    
    plt.figure(figsize=(8, 6))
    
    # Plot Lorenz curves
    plt.plot(c_weights, c_shares, label=f"Consumption (Gini = {gini_c:.3f})", color="royalblue", linewidth=2)
    plt.plot(w_weights, w_shares, label=f"Wealth (Gini = {gini_w:.3f})", color="firebrick", linewidth=2)
    
    # Line of equality
    plt.plot([0, 1], [0, 1], linestyle="--", color="black", alpha=0.7, label="Line of Equality")
    
    plt.title(title, fontsize=14)
    plt.xlabel("Cumulative Share of Population", fontsize=12)
    plt.ylabel("Cumulative Share of Value", fontsize=12)
    plt.legend(loc="upper left", fontsize=12)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


# ── Simplified Sequence-Jacobian Blocks ───────────────────────────────────────

def hh_inner_D(D, wealth, c):
    gini_c      = np.full_like(D, calculate_gini(c,      D))
    gini_wealth = np.full_like(D, calculate_gini(wealth, D))
    return gini_c, gini_wealth


def hh_inner_F(D_F, wealth, c):
    gini_c      = np.full_like(D, calculate_gini(c,      D_F))
    gini_wealth = np.full_like(D, calculate_gini(wealth, D_F))
    return gini_c, gini_wealth