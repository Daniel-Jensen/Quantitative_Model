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


def _level_points_m(m):
    # m CHEBYSHEV-EXTREMA NODES ON [-1, 1] (the dense factor of a refined grid).
    if m < 2:
        return np.array([0.0])
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
    # The box may be stated in ROTATED coordinates z = rot @ (x - centre) (Bocola's
    # V-transform, model_solution_mean.m): the collocation box is drawn on z while
    # every public method still speaks NATURAL coordinates x, so RuleSet / point_map /
    # the drivers are unchanged. rot=None is the identity and reproduces the
    # axis-aligned grid bit for bit.

    def __init__(self, lo, hi, mu=2, mu_vec=None, rot=None, centre=None,
                 refine=None):
        # BUILD POINTS, BASIS DEGREES, AND THE LU-FACTORED COLLOCATION MATRIX.
        # refine=(dim, m) makes the grid a CARTESIAN PRODUCT of a sparse Smolyak grid
        # over the other dimensions and a DENSE m-node Chebyshev grid in `dim`. Raising
        # one dimension's Smolyak level instead raises the GLOBAL budget -- on the
        # 10-state box [1,...,1,2,1] costs 165 points -- whereas the tensor factor buys
        # degree m-1 in that one dimension and full interaction with the sparse basis for
        # nb*m points. Both factors are square interpolation operators, so the product is
        # square and exact at its nodes, and every method below is unchanged.
        # Measured on a test function with this model's curvature profile (near-linear in
        # the wealth states, logistic in s, with an interaction), relative RMS error:
        #   isotropic mu=1   21 pts  1.9e-1      sparse(mu=1) x cheb_s(5)   95 pts  2.5e-2
        #   isotropic mu=2  221 pts  3.9e-2      sparse(mu=1) x cheb_s(9)  171 pts  1.1e-3
        # i.e. 36x better than the isotropic mu=2 grid at fewer points, and the curvature
        # it resolves is exactly the logistic p^d(s) that the mu=1 quadratic mis-fits
        # (3.30% against a true 1.98% at the headline shock).
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
        # inverse cached once; rot must be invertible but need NOT be orthogonal
        # (the eigenbasis of a non-symmetric transition Jacobian generally is not)
        self.rot_inv = None if self.rot is None else np.linalg.inv(self.rot)

        self.refine = None if refine is None else (int(refine[0]), int(refine[1]))
        if self.refine is None:
            base_d, base_mu_vec, keep = self.d, self.mu_vec, None
        else:
            # sparse factor over every dimension EXCEPT the refined one
            keep = [j for j in range(self.d) if j != self.refine[0]]
            base_d, base_mu_vec = self.d - 1, self.mu_vec[keep]
        pts, degs = [], []
        for i_vec in _multi_indices(base_d, self.mu, base_mu_vec):
            axes_p = [_new_points(i) for i in i_vec]
            axes_d = [_new_degrees(i) for i in i_vec]
            mesh_p = np.meshgrid(*axes_p, indexing="ij")
            mesh_d = np.meshgrid(*axes_d, indexing="ij")
            pts.append(np.column_stack([m.ravel() for m in mesh_p]))
            degs.append(np.column_stack([m.ravel() for m in mesh_d]))
        pts = np.vstack(pts)
        degs = np.vstack(degs).astype(int)
        if self.refine is None:
            self.points_unit, self.degrees = pts, degs
        else:
            r, m_s = self.refine
            u_s = _level_points_m(m_s)                 # dense Chebyshev extrema
            d_s = np.arange(m_s)                       # degrees 0 .. m-1
            nb = pts.shape[0]
            P = np.empty((nb * m_s, self.d))
            G = np.empty((nb * m_s, self.d), dtype=int)
            P[:, keep] = np.repeat(pts, m_s, axis=0)
            G[:, keep] = np.repeat(degs, m_s, axis=0)
            P[:, r] = np.tile(u_s, nb)
            G[:, r] = np.tile(d_s, nb)
            self.points_unit, self.degrees = P, G
        self.n = self.points_unit.shape[0]
        self.points = self.from_unit(self.points_unit)
        self.max_deg = int(self.degrees.max())
        self._Phi = self._basis_unit(self.points_unit)     # (n, n) collocation basis
        self._lu = lu_factor(self._Phi)

    def _fwd(self, x):
        # NATURAL -> ROTATED BOX COORDINATES (identity when rot is None).
        x = np.atleast_2d(x)
        return x if self.rot is None else (x - self.centre) @ self.rot.T

    def _bwd(self, z):
        # ROTATED BOX COORDINATES -> NATURAL (identity when rot is None).
        z = np.atleast_2d(z)
        return z if self.rot is None else z @ self.rot_inv.T + self.centre

    def to_unit(self, x):
        # MAP NATURAL COORDINATES TO [-1, 1]^d.
        return 2.0 * (self._fwd(x) - self.lo) / (self.hi - self.lo) - 1.0

    def from_unit(self, u):
        # MAP [-1, 1]^d COORDINATES TO THE NATURAL BOX.
        return self._bwd(self.lo + 0.5 * (np.atleast_2d(u) + 1.0) * (self.hi - self.lo))

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
        # Clipping happens in the ROTATED coordinates the box is drawn on, so a
        # rotated box projects onto its own faces, not onto an axis-aligned hull.
        return self._bwd(np.clip(self._fwd(x), self.lo, self.hi))

    def outside(self, x):
        # PER-DIMENSION BOX VIOLATION IN ROTATED COORDS, AS A FRACTION OF BOX WIDTH.
        # Zero inside; the accuracy diagnostic uses it to report how far the ergodic
        # path leaves the collocation box instead of silently clipping.
        z = self._fwd(x)
        return np.maximum(np.maximum(self.lo - z, z - self.hi), 0.0) / (self.hi - self.lo)


# STATE ORDER (2026-08-25, 9 states): THE ONE PLACE THE CONVENTION LIVES.
#   0 K_D  1 K_F  2 P_D  3 P_F  4 b_DD  5 b_DF  6 V_dep  7 s  8 Z_D
# P_X = (1+rdep_X,t-1)*dep_X,t-1 is the BANK's gross deposit obligation. It used to
# double as the household's gross CLAIM, which is what let the state vector stop at 7 --
# but that identity only holds under NATIONAL deposit clearing, where each household is
# force-fed exactly its own bank's funding need. Under the union deposit market the two
# diverge, and what they diverge BY is the carried cross-border deposit position
#     V = W_D - P_D    =>    W_D = P_D + V,    W_F = P_F - V/p
# so V is the state and BOTH household claims come off it. Carrying W_D instead and
# deriving W_F from the union identity (P_D + p*P_F - W_D)/p was measurably worse: W_F
# then inherits the fit error of P_D, P_F and W_D at once, and because C_F = W_F/P_CES
# + inc - A_F is a ~0.79 difference of ~8-sized terms, a 0.06% error in W_F became a
# 1.1% error in euler_F -- 22x the D-side error, corr(|euler_F|, |W_F gap|) = 0.97.
# V is ZERO at the symmetric SS and small everywhere, so neither claim is a difference
# of large numbers and the two countries are symmetric in their error. Its law of
# motion is exactly V' = (1+rdep_D)*nfa_dep_D (the WC deductions cancel).
# b_DD / b_DF are the two banks' carried holdings of the D sovereign, and b_FD is the D
# bank's carried holding of the F sovereign (b_FF = B_gov_F - b_FD closes it; the F stock
# itself is fixed, so the SPLIT is the only F state needed). Each was one fixed SS share
# until the corresponding bond market was made to clear; once a bank bids through its own
# FOC, last period's split is genuinely part of the state (the HM perpetuity pays
# payD*b_lag PER HOLDER). Leaving b_FD fixed meant the D bank held F bonds with NO
# first-order condition justifying the position -- capital and the D bond had their FOC in
# the residual system and the F leg did not -- which also froze Q_bF and made a
# flight-to-safety substitution impossible by construction.
# d_D in {0,1} is the discrete default regime (separate coefficient sets).
# THE CB BACKSTOP ADDS NO STATE. An LTRO changes the COMPOSITION of the bank's funding
# -- divertable deposits for non-divertable central-bank credit -- at an unchanged rate,
# so no stock is carried across periods and no budget identity moves; the whole effect is
# in the incentive constraint. (A BOND-PURCHASE backstop does need a state for the CB's
# book. That was built and measured: purchases can only remove the liquidity premium,
# 0.2-0.7% of the price here, because they work by pushing mu down and mu is floored at
# zero. See docs/ltro_backstop_plan.md and git history for the implementation.)
# phi, the per-period activation probability, is deliberately NOT a state either: its box
# [0,1] centres at 0.5, so no collocation node would have phi = 0 with every other state
# at its own centre and the steady state would stop being a grid point. It is a
# per-experiment scalar (cal["phi_ltro"]) -- one solve per activation, each EXACT at its
# own phi rather than quadratically interpolated.
STATE_NAMES = ("K_D", "K_F", "P_D", "P_F", "b_DD", "b_DF", "b_FD", "V_dep",
               "s", "Z_D")


def _band(spec, default):
    # (LOWER, UPPER) FRACTIONAL BAND FROM A SCALAR OR A (lo, hi) PAIR.
    if spec is None:
        spec = default
    if np.isscalar(spec):
        return float(spec), float(spec)
    lo, hi = spec
    return float(lo), float(hi)


# NAMED STATE INDICES, single-sourced. Every experiment that advances the state by hand
# imports these instead of writing S[5]/S[6], which is how the 7->9 state change would
# otherwise have silently relabelled s and Z_D as b_DF and W_D.
IK_D, IK_F, IP_D, IP_F, IBDD, IBDF, IBFD, IV, IS, IZ = range(10)
NSTATE = len(STATE_NAMES)

# s-BOX COVERAGE in unconditional sd of the s process (see build_state_box).
S_COVER_SD = 2.75


def build_state_box(ss, cal, s_lo=None, s_hi=None, s_halfwidth=None, k_band=0.03,
                    p_band=0.25, p_band_D=None, p_band_F=None,
                    b_band=0.30, b_lo_frac=None, mu=2, mu_vec=None, z_band=0.03, w_band=0.04,
                    rot=None, centre=None, refine=None):
    # THE 9-STATE BOX AROUND THE (RISKY) STEADY STATE, WIDE-LOW WHERE DEFAULT
    # CUTS STOCKS. b_lo_frac sets the B LOWER bound as a fraction of B_ss (default
    # 1-b_band); the default regime's surviving debt (~recovery*B) needs it near
    # recovery_rate so the d=1 continuation stays ON-GRID. mu is the isotropic
    # Smolyak level; mu_vec (per-state levels) overrides it for anisotropic
    # refinement where the constraint boundary moves (s, P) without paying for
    # resolution in the near-fixed K.
    #
    # s BOX: s* +- s_halfwidth, and the halfwidth is Bocola's own COVERAGE -- +-2.16
    # UNCONDITIONAL sd of the s process (model_solution_mean.m bounds(6,:) =
    # [-4.35,4.35] around his s* = -7.06) -- so it is computed from the process, not
    # hard-wired. It used to be the literal 4.35, which is +-2.16 sd only at
    # sigma_s = 0.63; at the corrected 0.4455 (see calibration.py) the same number is
    # +-3.05 sd, 41% wider than Bocola covers. The extra width is pure cost: out there
    # p^d ~ 0 and Q_bD must be flat at the risk-free perpetuity price ~0.946, but the
    # mu=1 quadratic OVERSHOOTS to 0.981-0.986 and even turns non-monotone, and that
    # overshoot is what feeds the bond-FOC residual (corr(|bondFOC_D|, |s-s*|) = 0.885,
    # fitted Q_bD off by up to 0.95%). Explicit s_lo/s_hi still override.
    # The earlier hard-wired [-9,-3.5] was only +-1.39/+2.26 CONDITIONAL innovation sd
    # wide against sigma_s = 1.5075, so grid.clip truncated ~18% of the quadrature mass
    # and cut the effective persistence of s from 0.95 to 0.80 -- a 62% attenuation of
    # the long bond's repricing, the mechanism the model exists to measure. Coverage in
    # sd units is the invariant; an absolute halfwidth silently breaks on any sigma change.
    #
    # p_band_D / p_band_F take a scalar or an explicit (lower, upper) pair; both fall
    # back to p_band. NOTE (measured, docs): the P-transition Jacobian at the SS has
    # eigenvalues 0.766 / -0.942 (stable) but |J| has spectral radius 1.96, so NO
    # axis-aligned box centred on the SS is one-step invariant. Pass rot/centre (see
    # recursive_main.p_block_rotation) to draw the box on the eigenbasis instead,
    # where the map is diagonal and every box IS invariant.
    bk_D, bk_F = ss["ss_bank_D"], ss["ss_bank_F"]
    K_D, K_F = ss["Kap_D_ss"], ss["Kap_F_ss"]
    P_D = bk_D["P_state_ss"]                   # net of the WC receivable, see bank.py
    P_F = bk_F["P_state_ss"]
    # THE D-SOVEREIGN STOCK IS CARRIED AS ITS TWO HOLDINGS, NOT AS ONE TOTAL. Once the
    # split is chosen by the banks' own FOCs it is no longer a fixed share of B, so last
    # period's split is genuinely part of the state: the HM perpetuity payoff is
    # payD*b_lag PER HOLDER. B_D = b_DD + b_DF is recovered wherever the total is needed.
    b_DF = cal["b_D_F_ss"]
    b_DD = cal["B_gov_D_ss"] - b_DF
    b_FD = cal["b_F_D_ss"]                     # D bank's holding of the F sovereign
    Z_D = cal["Z_ss_D"]                        # deterministic TFP state (9th dim)
    _sp = s_process_params(cal)
    s_star = _sp["s_star"]
    if s_halfwidth is None:
        # COVERAGE IN UNCONDITIONAL sd. Bocola's own is 2.16, but the box must also
        # CONTAIN the experiment and it cannot be shifted: s* has to stay the box centre
        # or it stops being a collocation node and the exact SS rest point goes with it.
        # The headline shock (p^d 0.10% -> 1.98%) is itself 2.11 sd, so at 2.16 it sits
        # at 97% of the half-width -- on the boundary, where the fit is worst and the
        # dynamic IRF would start at the edge. 2.75 puts it at 77% with room to move,
        # and is still 10% narrower than the 3.05 sd that the literal 4.35 halfwidth
        # became once sigma_s was corrected.
        s_halfwidth = S_COVER_SD * _sp["sigma_s"] / np.sqrt(1.0 - _sp["rho_s"] ** 2)
    s_lo = s_star - s_halfwidth if s_lo is None else s_lo
    s_hi = s_star + s_halfwidth if s_hi is None else s_hi
    pD_lo, pD_hi = _band(p_band_D, p_band)
    pF_lo, pF_hi = _band(p_band_F, p_band)
    b_lo_f = (1 - b_band if b_lo_frac is None else b_lo_frac)
    # V IS ZERO AT THE SS, so its band is ABSOLUTE (a fractional band round 0 collapses
    # the dimension). w_band is read as a fraction OF P_D, which keeps the caller's units
    # comparable to the other wealth states: the measured ergodic |nfa| runs to 0.137
    # against P_D = 7.74, so w_band = 0.04 gives +-0.31 -- roughly 2x the ergodic reach.
    V_half = w_band * P_D
    lo = np.array([(1 - k_band) * K_D, (1 - k_band) * K_F,
                   (1 - pD_lo) * P_D, (1 - pF_lo) * P_F,
                   b_lo_f * b_DD, b_lo_f * b_DF, b_lo_f * b_FD, -V_half,
                   s_lo, (1 - z_band) * Z_D])
    hi = np.array([(1 + k_band) * K_D, (1 + k_band) * K_F,
                   (1 + pD_hi) * P_D, (1 + pF_hi) * P_F,
                   (1 + b_band) * b_DD, (1 + b_band) * b_DF, (1 + b_band) * b_FD,
                   +V_half,
                   s_hi, (1 + z_band) * Z_D])
    if rot is not None:
        # The bands above are NATURAL half-widths t. The collocation box lives on
        # z = rot(x - centre), and a z-box |z| <= b reaches natural half-widths
        # |rot^-1| b, so solve |rot^-1| b = t for the z half-widths. (For the
        # un-rotated dimensions rot^-1 is the identity there and b = t exactly.)
        rot = np.asarray(rot, dtype=float)
        centre = np.asarray(centre, dtype=float)
        t = 0.5 * (hi - lo)
        A = np.abs(np.linalg.inv(rot))
        b = np.linalg.solve(A, t)
        if np.any(b <= 0.0):
            b = A @ t          # fallback: a superset box, still invariant on z
        mid = rot @ (0.5 * (lo + hi) - centre)
        lo, hi = mid - b, mid + b
    # refine=(dim, m) tensors a DENSE m-node Chebyshev factor onto one dimension.
    # The risk experiment passes refine=(IS, m): the logistic p^d(s) is where all the
    # curvature is, and raising the Smolyak level to reach it would pay for resolution
    # in nine other dimensions that are near-linear. See SmolyakGrid.__init__.
    return SmolyakGrid(lo, hi, mu=(mu if mu_vec is None else int(max(mu_vec))),
                       mu_vec=mu_vec, rot=rot, centre=centre, refine=refine)


def default_prob(s):
    # PRICED ONE-QUARTER-AHEAD DEFAULT PROBABILITY: LOGISTIC IN THE s FACTOR.
    return 1.0 / (1.0 + np.exp(-np.asarray(s, dtype=float)))


def s_process_params(cal):
    # AR(1) FOR THE SOVEREIGN-RISK FACTOR s. p^d(s*) = 0.1% at rest; rho_s = 0.95 and
    # sigma_s = 0.63 are Bocola's Table 2 posterior means (param(22), param(23) of
    # model_solution_mean.m).
    # sigma_s USED to be sized so a SINGLE +2sd innovation lifted p^d from 0.1% to 2%,
    # which forces sigma_s = 1.5075 -- 2.4x Bocola's. No plausible box holds that: it
    # made the box only +-1.4/2.3 conditional sd wide, so the box clip mean-reverted s
    # far harder than rho_s does (effective rho 0.80) and the long bond stopped pricing
    # persistent risk. Bocola never asks one quarter to do that: p^d reaches 2% because
    # a small-sigma, rho = 0.95 process WANDERS there. The experiments set the level of
    # s directly, so nothing downstream depends on the one-step interpretation.
    # NB Bocola's own GaussHermite.m maps z = rho*x + sigma*node with PHYSICISTS'
    # nodes and no sqrt(2), so his solved model contains sigma_eff = 0.63/sqrt(2) =
    # 0.4455; set cal["sigma_s"] = 0.4455 to reproduce his numbers exactly rather than
    # his reported parameter.
    s_star = np.log(0.001 / (1.0 - 0.001))            # p^d(s*) = 0.1%
    rho_s = 0.95                                      # Bocola param(22)
    sigma_s = float(cal.get("sigma_s", 0.63))         # Bocola param(23)
    # TFP (Z_D) as a DETERMINISTIC 7th state: perfect-foresight AR(1), no innovation
    # (Z is never shocked on the ergodic set -- Z_star = Z_ss -- so its rule slice is
    # exercised only off-grid, in the TFP experiment that reads along a Z-decay path).
    return dict(s_star=s_star, rho_s=rho_s, sigma_s=sigma_s,
                z_star=cal["Z_ss_D"], rho_z=0.9)
