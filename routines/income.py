import numpy as np

# Compute income z_i for each e_i in e_grid, given parameters w, vphi, frisch, lambda_tax, and tau
# From that calculate tax revenue T = E_D[tax paid] = E_D[x_i - z_i], where x_i is pre-tax income and z_i is after-tax income


def income_1(e_grid, w, vphi, frisch, lambda_tax, tau):

    eff_wage = e_grid * w                                    # e_i * w
    exp_L    = frisch / (1 + tau * frisch)                   # exponent for L
    L        = ((1 - tau) * lambda_tax* eff_wage ** (1 - tau) / vphi) ** exp_L    # optimal L_i
    x_grid   = eff_wage * L                                 # gross income
    z_grid   = lambda_tax * (eff_wage * L) ** (1 - tau)      # after-tax income
    return z_grid, x_grid



def income_2(e_grid, w, vphi, frisch, lambda_tax, tau):

    eff_wage = e_grid * w                                    # e_i * w
    exp_L    = frisch / (1 + tau * frisch)                   # exponent for L
    L        = ((1 - tau) * lambda_tax* eff_wage ** (1 - tau) / vphi) ** exp_L    # optimal L_i
    x_grid   = eff_wage * L                                 # gross income
    z_grid   = lambda_tax * (eff_wage * L) ** (1 - tau)      # after-tax income
    return z_grid, x_grid



def income_linear(e_grid, w, vphi, frisch, t_lin):
    L      = ((1 - t_lin) * w * e_grid / vphi) ** frisch
    x_grid = w * e_grid * L          # gross labor income
    z_grid = (1 - t_lin) * x_grid   # after-tax income
    return z_grid, x_grid


def income_hbfed(e_grid, T, Y, alpha, div):
    e_norm = e_grid / e_grid.mean()
    z_grid = ((1 - alpha) * Y - T) * e_norm + div
    return z_grid


def tax_revenue_1(D, x_grid, z_grid):
    D_z = D.sum(axis=(1, 2))                 # marginal over (b, a) → shape (nZ,)
    T   = np.dot(D_z, x_grid - z_grid)       # E_D[tax paid]
    return T



def tax_revenue_2(D, x_grid, z_grid):
    D_z = D.sum(axis=(1, 2))                 # marginal over (b, a) → shape (nZ,)
    T   = np.dot(D_z, x_grid - z_grid)       # E_D[tax paid]
    return T