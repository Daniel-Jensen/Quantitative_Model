# PCA-ROTATED SMOLYAK GRID: DECORRELATES THE ENDOGENOUS STATES (K,B,P).
# The ergodic cloud is elongated along a diagonal (more capital -> more
# deposits -> higher obligations), so an axis-aligned box containing it has
# deeply-infeasible anti-correlated corners (low assets + high obligations)
# that the model never visits but that corrupt the collocation fit. Gridding in
# the principal-component coordinates of the ergodic covariance aligns the box
# with the cloud, so the corners sit near its true extent. The exogenous states
# (Dz, g) keep their physical axes. This class presents a PHYSICAL-coordinate
# interface (points, eval, fit, clip) so the period map / expectations /
# simulator use it unchanged; the Chebyshev basis lives in rotated space.
import numpy as np

from smolyak import SmolyakGrid

NENDOG = 3   # K, B, P are rotated; Dz, g (indices 3,4) are passed through.


class RotatedGrid:
    def __init__(self, mean, V, pc_lo, pc_hi, exog_lo, exog_hi, mu=2):
        # mean, V: ergodic mean and eigenvector matrix (columns) of (K,B,P).
        # pc_lo/pc_hi: box in principal-component coords; exog_*: Dz,g box.
        self.mean = np.asarray(mean, float)
        self.V = np.asarray(V, float)
        lo = np.concatenate([pc_lo, exog_lo])
        hi = np.concatenate([pc_hi, exog_hi])
        self._g = SmolyakGrid(lo, hi, mu=mu)
        self.n = self._g.n
        self.d = self._g.d
        self.points = self.to_phys(self._g.points)      # physical collocation pts
        # physical bounding box (for reporting / interior checks)
        corners = self.to_phys(self._g.clip(self._g.points))
        self.lo = self.points.min(axis=0)
        self.hi = self.points.max(axis=0)

    def to_grid(self, x):
        # PHYSICAL STATE -> ROTATED GRID COORDS.
        x = np.atleast_2d(np.asarray(x, float))
        pc = (x[:, :NENDOG] - self.mean) @ self.V
        return np.column_stack([pc, x[:, NENDOG:]])

    def to_phys(self, w):
        # ROTATED GRID COORDS -> PHYSICAL STATE.
        w = np.atleast_2d(np.asarray(w, float))
        endog = self.mean + w[:, :NENDOG] @ self.V.T
        return np.column_stack([endog, w[:, NENDOG:]])

    def fit(self, values):
        # COLLOCATION COEFFICIENTS FROM VALUES AT self.points (SAME ORDER).
        return self._g.fit(values)

    def eval(self, coeffs, x):
        # EVALUATE THE INTERPOLANT AT PHYSICAL POINTS.
        return self._g.eval(coeffs, self.to_grid(x))

    def clip(self, x):
        # PROJECT PHYSICAL POINTS INTO THE ROTATED BOX (RETURN PHYSICAL).
        return self.to_phys(self._g.clip(self.to_grid(x)))

    def sample_interior(self, n, rng):
        # UNIFORM DRAWS INSIDE THE ROTATED PARALLELEPIPED (RETURN PHYSICAL).
        w = rng.uniform(self._g.lo, self._g.hi, size=(n, self.d))
        return self.to_phys(w)


def rotated_from_path(path, mu=2, k_endog=4.0, k_exog=3.5, exog_lo=None,
                      exog_hi=None):
    # BUILD A ROTATED GRID FROM AN ERGODIC PATH: PCA ON (K,B,P), physical rest.
    # Handles any number of trailing exogenous dims (2 restricted, 3 full).
    # exog_lo/exog_hi override the +-k_exog*sd defaults (e.g. to keep the s box
    # within the solved range).
    endog = path[:, :NENDOG]
    mean = endog.mean(axis=0)
    cov = np.cov(endog, rowvar=False)
    evals, V = np.linalg.eigh(cov)                  # ascending; V columns orthonormal
    sd = np.sqrt(np.maximum(evals, 1e-16))
    pc_lo, pc_hi = -k_endog * sd, k_endog * sd
    ex = path[:, NENDOG:]
    if exog_lo is None:
        exog_lo = ex.mean(0) - k_exog * ex.std(0)
        exog_hi = ex.mean(0) + k_exog * ex.std(0)
    return RotatedGrid(mean, V, pc_lo, pc_hi, np.asarray(exog_lo, float),
                       np.asarray(exog_hi, float), mu=mu)
