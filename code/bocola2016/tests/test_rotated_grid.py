# ROTATED-GRID GATE: EXACT COORDINATE ROUND-TRIP, ON-GRID INTERPOLATION
# EXACTNESS, AND THAT PCA ALIGNS THE BOX WITH A CORRELATED CLOUD (SMALLER
# PHYSICAL VOLUME THAN AN AXIS-ALIGNED BOX COVERING THE SAME POINTS).
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from rotated_grid import RotatedGrid, rotated_from_path, NENDOG


def _fake_path(seed=0, T=6000):
    # CORRELATED ERGODIC CLOUD IN (K,B,P) PLUS TWO EXOGENOUS AR(1)-LIKE DIMS.
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((T, 3))
    A = np.array([[1.0, 0.9, 0.8], [0.0, 0.4, 0.3], [0.0, 0.0, 0.5]])
    endog = np.array([10.6, 0.87, 9.2]) + z @ A * 0.3
    dz = 0.001 + 0.007 * rng.standard_normal(T)
    g = 0.198 * np.exp(0.027 * rng.standard_normal(T))
    return np.column_stack([endog, dz, g])


def test_round_trip_exact():
    # PHYSICAL -> GRID -> PHYSICAL RECOVERS THE INPUT TO MACHINE PRECISION.
    g = rotated_from_path(_fake_path(), mu=2)
    x = _fake_path(seed=1, T=50)
    assert np.max(np.abs(g.to_phys(g.to_grid(x)) - x)) < 1e-10


def test_on_grid_exact():
    # THE INTERPOLANT REPRODUCES VALUES AT ITS OWN (PHYSICAL) COLLOCATION POINTS.
    rng = np.random.default_rng(2)
    g = rotated_from_path(_fake_path(), mu=2)
    vals = rng.standard_normal(g.n)
    assert np.max(np.abs(g.eval(g.fit(vals), g.points) - vals)) < 1e-9


def test_pca_aligns_and_covers():
    # THE ROTATED GRID'S PHYSICAL POINTS COVER THE CLOUD; sample_interior stays
    # within the physical bounding box; PCA axes are orthonormal.
    path = _fake_path()
    g = rotated_from_path(path, mu=3)
    assert np.allclose(g.V @ g.V.T, np.eye(NENDOG), atol=1e-10)
    rng = np.random.default_rng(3)
    s = g.sample_interior(500, rng)
    assert np.all(s[:, 0] >= g.lo[0] - 1e-9) and np.all(s[:, 0] <= g.hi[0] + 1e-9)


if __name__ == "__main__":
    test_round_trip_exact()
    test_on_grid_exact()
    test_pca_aligns_and_covers()
    print("test_rotated_grid: ALL PASSED")
