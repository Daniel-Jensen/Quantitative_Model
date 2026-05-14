
import os
import numpy as np
import scipy.linalg
from sequence_jacobian import grids


_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "..", "Discretisation", "Outputs")



def make_grids_1(bmax, amax, kmax, nB, nA, nK, nZ, rho_z, sigma_z):

    b_grid = grids.agrid(amax=bmax, n=nB)
    a_grid = grids.agrid(amax=amax, n=nA)
    k_grid = grids.agrid(amax=kmax, n=nK)[::-1].copy()

    if nZ == 19:

        markov_ctstime = np.loadtxt(os.path.join(_DATA_DIR, "Px_GMAR.txt"))
        e_grid = np.loadtxt(os.path.join(_DATA_DIR, "x_vec.txt")).flatten()  # already in levels

        # Continuous-time → discrete-time
        markov_distime = scipy.linalg.expm(markov_ctstime)

        # Row normalize
        row_sums = markov_distime.sum(axis=1)
        Pi = markov_distime / row_sums[:, None]

    else:
        e_grid, _, Pi = grids.markov_rouwenhorst(rho=rho_z, sigma=sigma_z, N=nZ)

    return b_grid, a_grid, k_grid, e_grid, Pi



def make_grids_deposit(Depmax, nDep, nZ, rho_z, sigma_z):
    """Single-asset (deposit) grid — no illiquid or capital grid."""
    dep_grid = grids.agrid(amax=Depmax, n=nDep)

    if nZ == 19:
        markov_ctstime = np.loadtxt(os.path.join(_DATA_DIR, "Px_GMAR.txt"))
        e_grid = np.loadtxt(os.path.join(_DATA_DIR, "x_vec.txt")).flatten()
        markov_distime = scipy.linalg.expm(markov_ctstime)
        row_sums = markov_distime.sum(axis=1)
        Pi = markov_distime / row_sums[:, None]
    else:
        e_grid, _, Pi = grids.markov_rouwenhorst(rho=rho_z, sigma=sigma_z, N=nZ)

    return dep_grid, e_grid, Pi


def make_grids_2(bmax, amax, kmax, nB, nA, nK, nZ, rho_z, sigma_z):

    b_grid = grids.agrid(amax=bmax, n=nB)
    a_grid = grids.agrid(amax=amax, n=nA)
    k_grid = grids.agrid(amax=kmax, n=nK)[::-1].copy()

    if nZ == 19:

        markov_ctstime = np.loadtxt(os.path.join(_DATA_DIR, "Px_GMAR.txt"))
        e_grid = np.loadtxt(os.path.join(_DATA_DIR, "x_vec.txt")).flatten()  # already in levels

        # Continuous-time → discrete-time
        markov_distime = scipy.linalg.expm(markov_ctstime)

        # Row normalize
        row_sums = markov_distime.sum(axis=1)
        Pi = markov_distime / row_sums[:, None]

    else:
        e_grid, _, Pi = grids.markov_rouwenhorst(rho=rho_z, sigma=sigma_z, N=nZ)

    return b_grid, a_grid, k_grid, e_grid, Pi
