# SMOLYAK SPARSE GRID + CHEBYSHEV BASIS FOR THE RECURSIVE GLOBAL SOLUTION.
# Nested Chebyshev-extrema construction (Krueger-Kubler 2004): the grid is the
# union of tensor products of per-level "new point" sets over multi-indices i
# with sum(i_j - 1) <= mu; the basis uses the same index set over "new degree"
# sets, so the collocation matrix is square and the interpolant is exact at
# the nodes. Model-free; the agreed 8-state box builder lives at the bottom.
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
    return np.arange(m_prev, 2 ** (i - 1) + 1)


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
    # SPARSE COLLOCATION GRID ON A BOX [lo, hi]^d WITH A SQUARE CHEBYSHEV BASIS.

    def __init__(self, lo, hi, mu=2, mu_vec=None):
        # BUILD POINTS, BASIS DEGREES, AND THE LU-FACTORED COLLOCATION MATRIX.
        self.lo = np.asarray(lo, dtype=float)
        self.hi = np.asarray(hi, dtype=float)
        self.d = self.lo.size
        assert self.hi.shape == (self.d,) and np.all(self.hi > self.lo)
        self.mu = int(mu)
        self.mu_vec = (np.full(self.d, self.mu, dtype=int) if mu_vec is None
                       else np.asarray(mu_vec, dtype=int))

        pts, degs = [], []
        for i_vec in _multi_indices(self.d, self.mu, self.mu_vec):
            axes_p = [_new_points(i) for i in i_vec]
            axes_d = [_new_degrees(i) for i in i_vec]
            mesh_p = np.meshgrid(*axes_p, indexing="ij")
            mesh_d = np.meshgrid(*axes_d, indexing="ij")
            pts.append(np.column_stack([m.ravel() for m in mesh_p]))
            degs.append(np.column_stack([m.ravel() for m in mesh_d]))
        self.points_unit = np.vstack(pts)
        self.degrees = np.vstack(degs).astype(int)
        self.n = self.points_unit.shape[0]
        self.points = self.from_unit(self.points_unit)
        self.max_deg = int(self.degrees.max())
        self._Phi = self._basis_unit(self.points_unit)     # (n, n) collocation basis
        self._lu = lu_factor(self._Phi)

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

    def fit_weighted(self, values, w=None, ridge=0.0):
        # RIDGE-REGULARISED WEIGHTED LEAST-SQUARES FIT: the fix for global-fit corner
        # poisoning at mu=2. w in [0,1] per point (0 = point EXCLUDED from the fit, so
        # an unsolvable/frozen corner cannot leak into the coefficients); the ridge
        # penalises high-degree coefficients (damps the Gibbs wiggle at the kink) with
        # a small uniform floor for invertibility when points are masked. w=None,
        # ridge=0 falls back to the exact square solve (unchanged behaviour).
        values = np.asarray(values, dtype=float)
        if w is None and ridge == 0.0:
            return lu_solve(self._lu, values)
        w = np.ones(self.n) if w is None else np.asarray(w, dtype=float)
        A = self._Phi.T * w                                    # Phi^T @ diag(w)
        lhs = A @ self._Phi                                    # Phi^T W Phi
        diag_scale = float(np.mean(np.diag(lhs))) + 1e-12
        td = self.degrees.sum(axis=1).astype(float)            # total degree per basis fn
        reg = diag_scale * (ridge * td / max(td.max(), 1.0) + 1e-6)
        lhs[np.diag_indices_from(lhs)] += reg
        return np.linalg.solve(lhs, A @ values)

    def eval(self, coeffs, x):
        # EVALUATE THE INTERPOLANT AT ARBITRARY NATURAL-COORDINATE POINTS.
        return self.basis(x) @ coeffs

    def clip(self, x):
        # PROJECT POINTS INTO THE BOX (SIMULATION USE; EXPLICIT, NOT SILENT).
        return np.clip(np.atleast_2d(x), self.lo, self.hi)


# TIER-3 STATE ORDER (2026-07-30, reduced to 6): THE ONE PLACE THE CONVENTION LIVES.
#   0 K_D    1 K_F    2 P_D    3 P_F    4 B_D    5 s
# P_X = (1+rdep_X,t-1)*dep_X,t-1 (bank gross deposit obligation, own-good). This
# IS the household's gross deposit CLAIM (a deposit liability is a household asset),
# so the household's carried wealth is P_X/P_CES_X -- there is NO separate W state.
# Dropping the redundant W_D/W_F states removes the off-manifold points where
# household supply and bank demand could not both be set freely; deposits now clear
# by quantity (A = dep/P_CES) and the deposit RATE is pinned by the household Euler.
# d_D in {0,1} is the discrete default regime (separate coefficient sets).
STATE_NAMES = ("K_D", "K_F", "P_D", "P_F", "B_D", "s", "Z_D")


def build_state_box(ss, cal, s_lo=-9.0, s_hi=-3.5, k_band=0.03, p_band=0.25,
                    b_band=0.30, b_lo_frac=None, mu=2, mu_vec=None, z_band=0.03):
    # THE 7-STATE BOX AROUND THE (RISKY) STEADY STATE, WIDE-LOW WHERE DEFAULT
    # CUTS STOCKS. s spans p^d ~ 0.01% .. ~3% (logistic); bands sized from the
    # transition-path experience (K moves little -- kept under the ~4.5% Jermann
    # feasibility band; P/B are the movers). b_lo_frac sets the B LOWER bound as a
    # fraction of B_ss (default 1-b_band); the default regime's surviving debt
    # (~recovery*B) needs it near recovery_rate so the d=1 continuation stays
    # ON-GRID (else Bp_D ~ 0.45 B_ss extrapolates below a symmetric box). mu is
    # the isotropic Smolyak level; mu_vec (per-state levels, order K_D,K_F,P_D,
    # P_F,B_D,s) overrides it for anisotropic refinement where the constraint
    # boundary moves (s, P, B) without paying for resolution in the near-fixed K.
    bk_D, bk_F = ss["ss_bank_D"], ss["ss_bank_F"]
    K_D, K_F = ss["Kap_D_ss"], ss["Kap_F_ss"]
    P_D = (1.0 + cal["r_dep_D_target"]) * bk_D["Dep_supply_ss"]
    P_F = (1.0 + cal["r_dep_F_target"]) * bk_F["Dep_supply_ss"]
    B_D = cal["B_gov_D_ss"]
    Z_D = cal["Z_ss_D"]                        # deterministic TFP state (7th dim)
    b_lo = (1 - b_band if b_lo_frac is None else b_lo_frac) * B_D
    lo = np.array([(1 - k_band) * K_D, (1 - k_band) * K_F,
                   (1 - p_band) * P_D, (1 - p_band) * P_F,
                   b_lo, s_lo, (1 - z_band) * Z_D])
    hi = np.array([(1 + k_band) * K_D, (1 + k_band) * K_F,
                   (1 + p_band) * P_D, (1 + p_band) * P_F,
                   (1 + b_band) * B_D, s_hi, (1 + z_band) * Z_D])
    return SmolyakGrid(lo, hi, mu=(mu if mu_vec is None else int(max(mu_vec))),
                       mu_vec=mu_vec)


def default_prob(s):
    # PRICED ONE-QUARTER-AHEAD DEFAULT PROBABILITY: LOGISTIC IN THE s FACTOR.
    return 1.0 / (1.0 + np.exp(-np.asarray(s, dtype=float)))


def s_process_params(cal):
    # AR(1) FOR s MATCHED TO THE CURRENT EXPERIMENT (USER DECISION 2026-07-27):
    # p^d(s*) ~ 0.1% at rest; rho_s = the experiment's 0.95; sigma_s sized so a
    # +2sd innovation lifts p^d to ~2% (the pi_0 of the existing headline run).
    s_star = np.log(0.001 / (1.0 - 0.001))            # p^d(s*) = 0.1%
    s_high = np.log(0.02 / (1.0 - 0.02))              # p^d = 2% target
    rho_s = 0.95
    sigma_s = (s_high - s_star) / 2.0                 # 2sd one-shot innovation
    # TFP (Z_D) as a DETERMINISTIC 7th state: perfect-foresight AR(1), no innovation
    # (Z is never shocked on the ergodic set -- Z_star = Z_ss -- so its rule slice is
    # exercised only off-grid, in the TFP experiment that reads along a Z-decay path).
    return dict(s_star=s_star, rho_s=rho_s, sigma_s=sigma_s,
                z_star=cal["Z_ss_D"], rho_z=0.9)
