# SMOLYAK SPARSE GRID AND CHEBYSHEV BASIS (KRUEGER-KUBLER), PLUS THE MODEL'S STATE BOX.
import numpy as np
from scipy.linalg import lu_factor, lu_solve


def _level_points(i):
    # THE 1-D CHEBYSHEV EXTREMA OF LEVEL i.
    if i == 1:
        return np.array([0.0])
    m = 2 ** (i - 1) + 1
    return -np.cos(np.pi * np.arange(m) / (m - 1))


def _level_points_m(m):
    # m CHEBYSHEV EXTREMA ON [-1, 1].
    if m < 2:
        return np.array([0.0])
    return -np.cos(np.pi * np.arange(m) / (m - 1))


def _new_points(i):
    # THE POINTS NEW AT LEVEL i.
    if i == 1:
        return np.array([0.0])
    if i == 2:
        return np.array([-1.0, 1.0])
    return _level_points(i)[1::2]  # the odd positions are the new ones


def _new_degrees(i):
    # THE CHEBYSHEV DEGREES NEW AT LEVEL i.
    if i == 1:
        return np.array([0])
    if i == 2:
        return np.array([1, 2])
    m_prev = 2 ** (i - 2) + 1
    return np.arange(m_prev, 2 ** (i - 1) + 1)


def _multi_indices(d, mu, mu_vec):
    # EVERY LEVEL MULTI-INDEX WITHIN THE SMOLYAK BUDGET.
    out = []

    def rec(prefix, budget):
        # DEPTH-FIRST ENUMERATION.
        j = len(prefix)
        if j == d:
            out.append(tuple(prefix))
            return
        for lev in range(1, min(budget, mu_vec[j]) + 2):
            rec(prefix + [lev], budget - (lev - 1))

    rec([], mu)
    return out


def chebyshev_basis_1d(x, max_deg):
    # CHEBYSHEV POLYNOMIALS T_0..T_max_deg AT x.
    x = np.asarray(x, dtype=float)
    T = np.empty((x.size, max_deg + 1))
    T[:, 0] = 1.0
    if max_deg >= 1:
        T[:, 1] = x
    for k in range(2, max_deg + 1):
        T[:, k] = 2.0 * x * T[:, k - 1] - T[:, k - 2]
    return T


def _sparse_nodes(d, mu, mu_vec):
    # SMOLYAK NODES AND THEIR DEGREES ON [-1, 1]^d.
    pts, degs = [], []
    for i_vec in _multi_indices(d, mu, mu_vec):
        axes_p = [_new_points(i) for i in i_vec]
        axes_d = [_new_degrees(i) for i in i_vec]
        mesh_p = np.meshgrid(*axes_p, indexing="ij")
        mesh_d = np.meshgrid(*axes_d, indexing="ij")
        pts.append(np.column_stack([m.ravel() for m in mesh_p]))
        degs.append(np.column_stack([m.ravel() for m in mesh_d]))
    return np.vstack(pts), np.vstack(degs).astype(int)


def _tensor_refine(pts, degs, keep, refine, d):
    # A SPARSE GRID TIMES m DENSE NODES IN ONE DIMENSION.
    r, m_s = refine
    u_s = _level_points_m(m_s)
    d_s = np.arange(m_s)
    nb = pts.shape[0]
    P = np.empty((nb * m_s, d))
    G = np.empty((nb * m_s, d), dtype=int)
    P[:, keep] = np.repeat(pts, m_s, axis=0)
    G[:, keep] = np.repeat(degs, m_s, axis=0)
    P[:, r] = np.tile(u_s, nb)
    G[:, r] = np.tile(d_s, nb)
    return P, G


class SmolyakGrid:
    # SPARSE COLLOCATION GRID ON A (POSSIBLY ROTATED) BOX, WITH A SQUARE CHEBYSHEV BASIS.

    def __init__(self, lo, hi, mu=2, mu_vec=None, rot=None, centre=None,
                 refine=None):
        # BUILD THE NODES AND THE FACTORED COLLOCATION MATRIX.
        self.lo = np.asarray(lo, dtype=float)
        self.hi = np.asarray(hi, dtype=float)
        self.d = self.lo.size
        assert self.hi.shape == (self.d,) and np.all(self.hi > self.lo)
        self.mu = int(mu)
        self.mu_vec = (np.full(self.d, self.mu, dtype=int) if mu_vec is None
                       else np.asarray(mu_vec, dtype=int))
        self.rot = None if rot is None else np.asarray(rot, dtype=float)
        self.centre = (np.zeros(self.d) if centre is None
                       else np.asarray(centre, dtype=float))
        # rot need only be invertible, not orthogonal
        self.rot_inv = None if self.rot is None else np.linalg.inv(self.rot)

        self.refine = None if refine is None else (int(refine[0]), int(refine[1]))
        if self.refine is None:
            self.points_unit, self.degrees = _sparse_nodes(self.d, self.mu, self.mu_vec)
        else:
            # refined: sparse in every other dimension, dense in this one
            keep = [j for j in range(self.d) if j != self.refine[0]]
            pts, degs = _sparse_nodes(self.d - 1, self.mu, self.mu_vec[keep])
            self.points_unit, self.degrees = _tensor_refine(pts, degs, keep, self.refine,
                                                            self.d)
        self.n = self.points_unit.shape[0]
        self.points = self.from_unit(self.points_unit)
        self.max_deg = int(self.degrees.max())
        self._Phi = self._basis_unit(self.points_unit)
        self._lu = lu_factor(self._Phi)

    def _fwd(self, x):
        # NATURAL -> BOX COORDINATES.
        x = np.atleast_2d(x)
        return x if self.rot is None else (x - self.centre) @ self.rot.T

    def _bwd(self, z):
        # BOX COORDINATES -> NATURAL.
        z = np.atleast_2d(z)
        return z if self.rot is None else z @ self.rot_inv.T + self.centre

    def to_unit(self, x):
        # NATURAL COORDINATES -> [-1, 1]^d.
        return 2.0 * (self._fwd(x) - self.lo) / (self.hi - self.lo) - 1.0

    def from_unit(self, u):
        # [-1, 1]^d -> NATURAL COORDINATES.
        return self._bwd(self.lo + 0.5 * (np.atleast_2d(u) + 1.0) * (self.hi - self.lo))

    def _basis_unit(self, u):
        # THE BASIS MATRIX AT UNIT-BOX POINTS.
        u = np.atleast_2d(u)
        B = np.ones((u.shape[0], self.n))
        for j in range(self.d):
            Tj = chebyshev_basis_1d(u[:, j], self.max_deg)
            B *= Tj[:, self.degrees[:, j]]
        return B

    def basis(self, x):
        # THE BASIS MATRIX AT NATURAL POINTS (EXTRAPOLATES OUTSIDE THE BOX).
        return self._basis_unit(self.to_unit(x))

    def fit(self, values):
        # COLLOCATION COEFFICIENTS FROM VALUES AT THE NODES.
        return lu_solve(self._lu, np.asarray(values, dtype=float))

    def fit_weighted(self, values, w=None, ridge=0.0):
        # WEIGHTED RIDGE FIT; THE DEFAULTS GIVE THE EXACT SOLVE.
        values = np.asarray(values, dtype=float)
        if w is None and ridge == 0.0:
            return lu_solve(self._lu, values)
        w = np.ones(self.n) if w is None else np.asarray(w, dtype=float)
        A = self._Phi.T * w
        lhs = A @ self._Phi
        diag_scale = float(np.mean(np.diag(lhs))) + 1e-12
        td = self.degrees.sum(axis=1).astype(float)
        reg = diag_scale * (ridge * td / max(td.max(), 1.0) + 1e-6)
        lhs[np.diag_indices_from(lhs)] += reg
        return np.linalg.solve(lhs, A @ values)

    def eval(self, coeffs, x):
        # EVALUATE THE INTERPOLANT AT NATURAL POINTS.
        return self.basis(x) @ coeffs

    def clip(self, x):
        # PROJECT POINTS INTO THE BOX.
        return self._bwd(np.clip(self._fwd(x), self.lo, self.hi))

    def outside(self, x):
        # HOW FAR EACH POINT LIES OUTSIDE THE BOX, AS A SHARE OF ITS WIDTH.
        z = self._fwd(x)
        return np.maximum(np.maximum(self.lo - z, z - self.hi), 0.0) / (self.hi - self.lo)


# the 12 states in order; M_cb and O_cb are the TPI book, zero at the SS
STATE_NAMES = ("K_D", "K_F", "P_D", "P_F", "b_DD", "b_DF", "b_FD", "V_dep",
               "s", "Z_D", "M_cb", "O_cb")


def _band(spec, default):
    # A (LOWER, UPPER) BAND FROM A SCALAR OR A PAIR.
    if spec is None:
        spec = default
    if np.isscalar(spec):
        return float(spec), float(spec)
    lo, hi = spec
    return float(lo), float(hi)


# named state indices, so no caller writes S[5]
IK_D, IK_F, IP_D, IP_F, IBDD, IBDF, IBFD, IV, IS, IZ, IM, IO = range(12)
NSTATE = len(STATE_NAMES)

# the s box covers this many unconditional sd either side of s*
S_COVER_SD = 2.75


def build_state_box(ss, cal, s_lo=None, s_hi=None, s_halfwidth=None, k_band=0.03,
                    p_band=0.25, p_band_D=None, p_band_F=None,
                    b_band=0.30, b_lo_frac=None, mu=2, mu_vec=None, z_band=0.03, w_band=0.04,
                    m_band=1.0, o_band=0.25, rot=None, centre=None, refine=None):
    # THE STATE BOX AROUND THE STEADY STATE.
    c = _box_centre(ss, cal)
    s_lo, s_hi = _s_bounds(cal, s_lo, s_hi, s_halfwidth)
    pD_lo, pD_hi = _band(p_band_D, p_band)
    pF_lo, pF_hi = _band(p_band_F, p_band)
    b_lo_f = (1 - b_band if b_lo_frac is None else b_lo_frac)
    # V and the TPI book are zero at the SS, so their bands are absolute
    V_half = w_band * c["P_D"]
    M_half = m_band * c["b_DD"]
    O_half = o_band * M_half
    lo = np.array([(1 - k_band) * c["K_D"], (1 - k_band) * c["K_F"],
                   (1 - pD_lo) * c["P_D"], (1 - pF_lo) * c["P_F"],
                   b_lo_f * c["b_DD"], b_lo_f * c["b_DF"], b_lo_f * c["b_FD"], -V_half,
                   s_lo, (1 - z_band) * c["Z_D"], -M_half, -O_half])
    hi = np.array([(1 + k_band) * c["K_D"], (1 + k_band) * c["K_F"],
                   (1 + pD_hi) * c["P_D"], (1 + pF_hi) * c["P_F"],
                   (1 + b_band) * c["b_DD"], (1 + b_band) * c["b_DF"], (1 + b_band) * c["b_FD"],
                   +V_half,
                   s_hi, (1 + z_band) * c["Z_D"], +M_half, +O_half])
    if rot is not None:
        rot = np.asarray(rot, dtype=float)
        centre = np.asarray(centre, dtype=float)
        lo, hi = _rotate_box(lo, hi, rot, centre)
    else:
        rot, centre = np.eye(lo.size), np.zeros(lo.size)
    # the TPI shear: the box bounds b_DD + M and O - rho*M, where purchases actually go
    shear = np.eye(lo.size)
    shear[IBDD, IM] = 1.0
    shear[IO, IM] = -(1.0 + cal["r_dep_D_target"]) * ss["Q_bD_ss"]
    # refine adds dense nodes in s, where p^d(s) curves
    return SmolyakGrid(lo, hi, mu=(mu if mu_vec is None else int(max(mu_vec))),
                       mu_vec=mu_vec, rot=shear @ rot, centre=centre, refine=refine)


def _box_centre(ss, cal):
    # THE STEADY-STATE CENTRE OF THE BOX.
    b_DF = cal["b_D_F_ss"]
    return dict(K_D=ss["Kap_D_ss"], K_F=ss["Kap_F_ss"],
                P_D=ss["ss_bank_D"]["P_state_ss"],  # net of the working-capital receivable
                P_F=ss["ss_bank_F"]["P_state_ss"],
                b_DD=cal["B_gov_D_ss"] - b_DF, b_DF=b_DF,
                b_FD=cal["b_F_D_ss"],
                Z_D=cal["Z_ss_D"])


def _s_bounds(cal, s_lo, s_hi, s_halfwidth):
    # THE s BOUNDS: S_COVER_SD UNCONDITIONAL SD EITHER SIDE OF s*.
    _sp = s_process_params(cal)
    s_star = _sp["s_star"]
    if s_halfwidth is None:
        s_halfwidth = S_COVER_SD * _sp["sigma_s"] / np.sqrt(1.0 - _sp["rho_s"] ** 2)
    s_lo = s_star - s_halfwidth if s_lo is None else s_lo
    s_hi = s_star + s_halfwidth if s_hi is None else s_hi
    return s_lo, s_hi


def _rotate_box(lo, hi, rot, centre):
    # THE ROTATED BOX THAT REACHES THE NATURAL HALF-WIDTHS.
    t = 0.5 * (hi - lo)
    A = np.abs(np.linalg.inv(rot))
    b = np.linalg.solve(A, t)
    if np.any(b <= 0.0):
        b = A @ t  # fallback: a larger box
    mid = rot @ (0.5 * (lo + hi) - centre)
    return mid - b, mid + b


def default_prob(s):
    # THE PRICED DEFAULT PROBABILITY, LOGISTIC IN s.
    return 1.0 / (1.0 + np.exp(-np.asarray(s, dtype=float)))


def s_process_params(cal):
    # THE AR(1) FOR THE RISK FACTOR s.
    s_star = np.log(0.001 / (1.0 - 0.001))  # p^d = 0.1% at rest
    rho_s = 0.95  # Bocola's estimate
    sigma_s = float(cal.get("sigma_s", 0.63))  # the calibration sets Bocola's effective 0.4455
    # Z_D is deterministic: it decays at rho_z with no innovation
    return dict(s_star=s_star, rho_s=rho_s, sigma_s=sigma_s,
                z_star=cal["Z_ss_D"], rho_z=0.9)
