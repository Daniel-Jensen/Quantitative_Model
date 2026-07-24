# PERTURBATION GATE: EXACT BGP LINEARIZATION POINT, BLANCHARD-KAHN SATISFIED,
# STABLE TRANSITION, AND SIGN-SENSIBLE IMPULSE RESPONSES.
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings("ignore")

import numpy as np

from calibration import solve_bgp
from perturbation import solve_perturbation, w_bgp, F_det, irf, NX


def test_bgp_is_linearization_point():
    # THE DETERMINISTIC RESIDUAL VANISHES AT THE BGP (F(w_bg, w_bg) = 0).
    cal = solve_bgp(delta_mode="standard")
    w0 = w_bgp(cal)
    assert np.max(np.abs(F_det(w0, w0, cal))) < 1e-10


def test_blanchard_kahn_and_stability():
    # KLEIN SOLVES (5 STABLE ROOTS) AND THE TRANSITION MATRIX IS STABLE.
    cal = solve_bgp(delta_mode="standard")
    sol = solve_perturbation(cal)              # raises if BK fails
    ev = np.abs(np.linalg.eigvals(sol["P"]))
    assert np.all(ev < 1.0 + 1e-8), f"unstable transition: {ev}"
    assert np.max(ev) > 0.95                   # the persistent net-worth root


def test_irf_signs():
    # TFP GROWTH SHOCK RAISES CONSUMPTION ON IMPACT AND MEAN-REVERTS;
    # SPENDING SHOCK RAISES DEBT AND LOWERS THE BOND PRICE.
    cal = solve_bgp(delta_mode="standard")
    sol = solve_perturbation(cal)
    base = w_bgp(cal)

    z = irf(sol, 3, cal["sig_z"], T=12)        # +1sd Dz
    assert z[0, 5] > base[5], "C should rise on a positive TFP shock"
    assert abs(z[6, 3] - base[3]) < abs(z[0, 3] - base[3]), "Dz must revert"

    gsh = irf(sol, 4, cal["sig_g"] * cal["g_star"], T=12)   # +1sd g
    assert gsh[1, 1] > base[1], "govt debt B should rise on a spending shock"
    assert gsh[0, 7] < base[7], "bond price Q_B should fall on a spending shock"


if __name__ == "__main__":
    test_bgp_is_linearization_point()
    test_blanchard_kahn_and_stability()
    test_irf_signs()
    print("test_perturbation: ALL PASSED")
