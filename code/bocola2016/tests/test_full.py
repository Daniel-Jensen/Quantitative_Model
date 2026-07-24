# FULL-MODEL GATE: WITH pi^d=0 AND SHOCKS OFF, THE 6-STATE TWO-REGIME MACHINERY
# REDUCES TO THE RESTRICTED BGP REST POINT (period map, two-branch quadrature,
# and default plumbing all nest correctly). Also: at pi^d=0 the d=0 full-model
# residual equals the restricted residual at random states.
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings("ignore")

import numpy as np

from calibration import solve_bgp
from smolyak import SmolyakGrid
from expectations_full import gauss_hermite_3d, default_prob
from time_iteration_full import point_block_full
from time_iteration import point_block as point_block_restr
from expectations import gauss_hermite_2d


def _full_box(cal):
    # 6-DIM BOX AROUND THE BGP; s CENTERED AT s_star.
    sd = lambda sig, rho: sig / np.sqrt(1 - rho ** 2)
    K, B, P = cal["K_bg"], cal["B_bg"], cal["P_bg"]
    lo = np.array([0.85 * K, 0.75 * B, 0.8 * P, cal["gamma"] - 0.03,
                   np.log(cal["g_star"]) - 0.12, cal["s_star"] - 4.0])
    hi = np.array([1.15 * K, 1.25 * B, 1.2 * P, cal["gamma"] + 0.03,
                   np.log(cal["g_star"]) + 0.12, cal["s_star"] + 4.0])
    lo[4], hi[4] = np.exp(lo[4]), np.exp(hi[4])
    return lo, hi


def _const_full_coef(grid, cal):
    # CONSTANT BGP RULES IN BOTH REGIMES.
    n = grid.n
    base = dict(C=np.full(n, cal["C_bg"]), R=np.full(n, cal["R_bg"]),
                QB=np.full(n, cal["qb_bg"]), av=np.full(n, cal["av_bg"]))
    return {k: (grid.fit(v), grid.fit(v)) for k, v in base.items()}


def test_default_prob_logistic():
    # pi^d(s_star) MATCHES BOCOLA'S ~0.09% BGP DEFAULT PROBABILITY.
    cal = solve_bgp()
    assert abs(default_prob(cal["s_star"]) - 0.0009) < 5e-4
    assert abs(default_prob(0.0) - 0.5) < 1e-12


def test_full_bgp_rest_point_nesting():
    # pi^d=0, SHOCKS OFF: FULL d=0 REST POINT REPRODUCES THE RESTRICTED BGP.
    cal = solve_bgp(delta_mode="standard")
    cal = dict(cal); cal["no_default"] = True
    cal["sig_z"] = cal["sig_g"] = cal["sig_s"] = 0.0
    lo, hi = _full_box(cal)
    grid = SmolyakGrid(lo, hi, mu=2)
    coef = _const_full_coef(grid, cal)
    quad = gauss_hermite_3d(3)

    state = np.array([cal["K_bg"], cal["B_bg"], cal["P_bg"], cal["gamma"],
                      cal["g_star"], cal["s_star"]])
    x = np.array([cal["C_bg"], cal["R_bg"], cal["qb_bg"]])
    res, mu, av, st = point_block_full(x, state, 0, grid, coef, cal, quad)

    assert abs(st["Kp"] - cal["K_bg"]) < 1e-9
    assert abs(st["Bp"] - cal["B_bg"]) < 1e-9
    assert abs(st["Pp"] - cal["P_bg"]) < 1e-9
    assert np.max(np.abs(res)) < 1e-8, f"residuals {res}"
    assert abs(mu - cal["mu_bg"]) < 1e-8
    assert abs(av - cal["av_bg"]) < 1e-8


def test_full_matches_restricted_offbgp():
    # pi^d=0: THE d=0 FULL RESIDUAL EQUALS THE RESTRICTED RESIDUAL (rules flat
    # in s; the extra es-integration averages the same value).
    cal = solve_bgp(delta_mode="standard")
    calf = dict(cal); calf["no_default"] = True
    lo6, hi6 = _full_box(cal)
    g6 = SmolyakGrid(lo6, hi6, mu=2)
    g5 = SmolyakGrid(lo6[:5], hi6[:5], mu=2)
    coef6 = _const_full_coef(g6, cal)
    coef5 = {k: g5.fit(np.full(g5.n, cal[{"C": "C_bg", "R": "R_bg",
             "QB": "qb_bg", "av": "av_bg"}[k]])) for k in ("C", "R", "QB", "av")}
    q3, q2 = gauss_hermite_3d(3), gauss_hermite_2d(3)

    rng = np.random.default_rng(0)
    for _ in range(6):
        s5 = lo6[:5] + rng.random(5) * (hi6[:5] - lo6[:5])
        s6 = np.append(s5, cal["s_star"])
        x = np.array([cal["C_bg"], cal["R_bg"], cal["qb_bg"]])
        rf, _, _, _ = point_block_full(x, s6, 0, g6, coef6, calf, q3)
        rr, _, _, _ = point_block_restr(x, s5, g5, coef5, cal, q2)
        assert np.max(np.abs(rf - rr)) < 1e-9, f"full!=restricted: {rf-rr}"


if __name__ == "__main__":
    test_default_prob_logistic()
    test_full_bgp_rest_point_nesting()
    test_full_matches_restricted_offbgp()
    print("test_full: ALL PASSED")
