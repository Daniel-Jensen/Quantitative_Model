# SMOLYAK SPARSE GRID + CHEBYSHEV BASIS (KRUEGER-KUBLER 2004 NESTED CONSTRUCTION).
# Nested 1D Chebyshev-extrema levels; grid = union of tensor products of the
# per-level "new point" sets over multi-indices i with sum(i_j - 1) <= mu;
# basis = same index set over the per-level "new degree" sets, so the
# collocation matrix is square and the interpolant is exact at the nodes.
# Supports anisotropic per-dimension level caps (mu_vec). Model-free.
import numpy as np
from scipy.linalg import lu_factor, lu_solve


def _level_points(i):
    # FULL 1D CHEBYSHEV-EXTREMA SET OF LEVEL i (m = 1, 3, 5, 9, ... POINTS).
    if i == 1:
        return np.array([0.0])
    m = 2 ** (i - 1) + 1
    return -np.cos(np.pi * np.arange(m) / (m - 1))


def _new_points(i):
    # POINTS INTRODUCED AT LEVEL i (DISJOINT ACROSS LEVELS BY NESTEDNESS).
    if i == 1:
        return np.array([0.0])
    if i == 2:
        return np.array([-1.0, 1.0])
    return _level_points(i)[1::2]   # odd positions are absent from level i-1


def _new_degrees(i):
    # CHEBYSHEV DEGREES INTRODUCED AT LEVEL i (|new degrees| = |new points|).
    if i == 1:
        return np.array([0])
    if i == 2:
        return np.array([1, 2])
    m_prev = 2 ** (i - 2) + 1
    m = 2 ** (i - 1) + 1
    return np.arange(m_prev, m)


def _multi_indices(d, mu, mu_vec):
    # ALL LEVEL MULTI-INDICES i (EACH >= 1) WITH sum(i-1) <= mu AND i-1 <= mu_vec.
    out = []

    def rec(prefix, budget):
        # DEPTH-FIRST ENUMERATION UNDER THE REMAINING LEVEL BUDGET.
        j = len(prefix)
        if j == d:
            out.append(tuple(prefix))
            return
        for lev in range(1, min(budget, mu_vec[j]) + 2):
            rec(prefix + [lev], budget - (lev - 1))

    rec([], mu)
    return out


def chebyshev_basis_1d(x, max_deg):
    # T_0..T_max_deg AT POINTS x VIA THE RECURRENCE (SHAPE (len(x), max_deg+1)).
    x = np.asarray(x, dtype=float)
    T = np.empty((x.size, max_deg + 1))
    T[:, 0] = 1.0
    if max_deg >= 1:
        T[:, 1] = x
    for k in range(2, max_deg + 1):
        T[:, k] = 2.0 * x * T[:, k - 1] - T[:, k - 2]
    return T


class SmolyakGrid:
    # SPARSE COLLOCATION GRID ON A BOX [lo, hi]^d WITH SQUARE CHEBYSHEV BASIS.

    def __init__(self, lo, hi, mu=2, mu_vec=None):
        # BUILD POINTS, BASIS DEGREES, AND THE LU-FACTORED COLLOCATION MATRIX.
        self.lo = np.asarray(lo, dtype=float)
        self.hi = np.asarray(hi, dtype=float)
        self.d = self.lo.size
        assert self.hi.shape == (self.d,) and np.all(self.hi > self.lo)
        self.mu = int(mu)
        self.mu_vec = (np.full(self.d, self.mu, dtype=int) if mu_vec is None
                       else np.asarray(mu_vec, dtype=int))
        assert self.mu_vec.shape == (self.d,)

        idx = _multi_indices(self.d, self.mu, self.mu_vec)
        pts, degs = [], []
        for i_vec in idx:
            axes_p = [_new_points(i) for i in i_vec]
            axes_d = [_new_degrees(i) for i in i_vec]
            # tensor products of the per-level new sets (disjoint across i_vec)
            mesh_p = np.meshgrid(*axes_p, indexing="ij")
            mesh_d = np.meshgrid(*axes_d, indexing="ij")
            pts.append(np.column_stack([m.ravel() for m in mesh_p]))
            degs.append(np.column_stack([m.ravel() for m in mesh_d]))
        self.points_unit = np.vstack(pts)
        self.degrees = np.vstack(degs).astype(int)
        self.n = self.points_unit.shape[0]
        assert self.degrees.shape[0] == self.n

        self.points = self.from_unit(self.points_unit)
        self.max_deg = int(self.degrees.max())
        self._lu = lu_factor(self._basis_unit(self.points_unit))

    def to_unit(self, x):
        # MAP NATURAL COORDINATES TO [-1, 1]^d.
        return 2.0 * (np.atleast_2d(x) - self.lo) / (self.hi - self.lo) - 1.0

    def from_unit(self, u):
        # MAP [-1, 1]^d COORDINATES TO THE NATURAL BOX.
        return self.lo + 0.5 * (np.atleast_2d(u) + 1.0) * (self.hi - self.lo)

    def _basis_unit(self, u):
        # BASIS MATRIX AT UNIT-BOX POINTS: PRODUCTS OF PER-DIMENSION CHEBYSHEVS.
        u = np.atleast_2d(u)
        B = np.ones((u.shape[0], self.n))
        for j in range(self.d):
            Tj = chebyshev_basis_1d(u[:, j], self.max_deg)
            B *= Tj[:, self.degrees[:, j]]
        return B

    def basis(self, x):
        # BASIS MATRIX AT NATURAL-COORDINATE POINTS (EXTRAPOLATES OUTSIDE BOX).
        return self._basis_unit(self.to_unit(x))

    def fit(self, values):
        # COLLOCATION COEFFICIENTS FROM VALUES AT self.points ((n,) OR (n, k)).
        return lu_solve(self._lu, np.asarray(values, dtype=float))

    def eval(self, coeffs, x):
        # EVALUATE THE INTERPOLANT AT ARBITRARY NATURAL-COORDINATE POINTS.
        return self.basis(x) @ coeffs

    def clip(self, x):
        # PROJECT POINTS INTO THE BOX (FOR SIMULATION USE; DOCUMENTED, NOT SILENT).
        return np.clip(np.atleast_2d(x), self.lo, self.hi)

    def sample_interior(self, n, rng):
        # UNIFORM DRAWS INSIDE THE GRID DOMAIN (PHYSICAL == BOX HERE).
        return rng.uniform(self.lo, self.hi, size=(n, self.d))
