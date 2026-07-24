# RESTRICTED-MODEL GATE: BGP REST POINT (SHOCKS OFF => EXACT FIXED POINT),
# THEN THE STOCHASTIC SOLVE'S ACCURACY AND ECONOMIC SANITY.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from calibration import solve_bgp
from smolyak import SmolyakGrid
from expectations import gauss_hermite_2d
from time_iteration import initial_coef, point_block


def build_box(cal, wide=False):
    # STATE BOX AROUND THE BGP; EXOGENOUS DIMS CENTERED AT gamma, g_star.
    # Capital moves little (high K/Y) so its band is tight; government debt B
    # is the main endogenous mover. Ergodic exogenous ranges: Dz ~ +-3.5 sd,
    # log g ~ +-3.5 sd of their AR(1) stationary distributions.
    K, B, P = cal["K_bg"], cal["B_bg"], cal["P_bg"]
    sd_z = cal["sig_z"] / np.sqrt(1 - cal["rho_z"] ** 2)
    sd_g = cal["sig_g"] / np.sqrt(1 - cal["rho_g"] ** 2)
    kw, bw, pw = (0.14, 0.45, 0.30) if wide else (0.08, 0.25, 0.18)
    lo = np.array([(1 - kw) * K, (1 - bw) * B, (1 - pw) * P,
                   cal["gamma"] - 3.5 * sd_z, np.log(cal["g_star"]) - 3.5 * sd_g])
    hi = np.array([(1 + kw) * K, (1 + bw) * B, (1 + pw) * P,
                   cal["gamma"] + 3.5 * sd_z, np.log(cal["g_star"]) + 3.5 * sd_g])
    lo[4], hi[4] = np.exp(lo[4]), np.exp(hi[4])   # g in levels
    return lo, hi


def test_bgp_rest_point():
    # DETERMINISTIC BGP IS AN EXACT FIXED POINT OF THE PERIOD MAP + EULERS.
    cal0 = solve_bgp(delta_mode="standard")     # primary calibration
    lo, hi = build_box(cal0)                    # box from stochastic ergodic sd
    cal = dict(cal0); cal["sig_z"] = 0.0; cal["sig_g"] = 0.0   # then shocks off
    grid = SmolyakGrid(lo, hi, mu=2)
    coef, _ = initial_coef(grid, cal)             # constant BGP rules
    quad = gauss_hermite_2d(5)

    state = np.array([cal["K_bg"], cal["B_bg"], cal["P_bg"],
                      cal["gamma"], cal["g_star"]])
    x = np.array([cal["C_bg"], cal["R_bg"], cal["qb_bg"]])
    res, mu, av, st = point_block(x, state, grid, coef, cal, quad)

    # period map reproduces the BGP stocks exactly
    assert abs(st["Kp"] - cal["K_bg"]) < 1e-10, f"Kp {st['Kp']} vs {cal['K_bg']}"
    assert abs(st["Bp"] - cal["B_bg"]) < 1e-10, f"Bp {st['Bp']} vs {cal['B_bg']}"
    assert abs(st["Pp"] - cal["P_bg"]) < 1e-10, f"Pp {st['Pp']} vs {cal['P_bg']}"
    assert abs(st["N"] - cal["N_bg"]) < 1e-10, f"N {st['N']} vs {cal['N_bg']}"
    assert abs(st["I"] - cal["I_bg"]) < 1e-10, f"I {st['I']} vs {cal['I_bg']}"
    assert abs(st["Y"] - cal["Y_bg"]) < 1e-10, f"Y {st['Y']} vs {cal['Y_bg']}"
    assert abs(st["Q_K"] - 1.0) < 1e-10, f"Q_K {st['Q_K']} != 1"
    # Euler residuals vanish; mu, av recover their BGP values
    assert np.max(np.abs(res)) < 1e-9, f"residuals {res}"
    assert abs(mu - cal["mu_bg"]) < 1e-9, f"mu {mu} vs {cal['mu_bg']}"
    assert abs(av - cal["av_bg"]) < 1e-9, f"av {av} vs {cal['av_bg']}"


if __name__ == "__main__":
    test_bgp_rest_point()
    print("test_restricted (BGP rest point): PASSED")
