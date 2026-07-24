# CALIBRATION GATE: THE FAITHFUL (euler) BGP HITS EVERY STATED TABLE-1/2 TARGET
# EXACTLY, AND THE BACKED-OUT PARAMETERS SATISFY THE BGP IDENTITIES THEY CAME
# FROM. The implausible K/Y is a documented consequence, not a target, so it is
# asserted only to be internally consistent (i/y hits), not "reasonable".
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from calibration import base_params, solve_bgp


def test_targets_hit_exactly():
    # EVERY STATED FOOTNOTE-14 TARGET REPRODUCED TO 1e-12 IN euler MODE.
    p = base_params()
    b = solve_bgp(p, delta_mode="euler")
    # investment-output target (delta was backed out to hit this)
    assert abs(b["iy_implied"] - p["iy_bg"]) < 1e-12
    # risk-free / household Euler: beta = G/R exactly
    assert abs(b["beta"] * p["R_bg"] / b["G"] - 1.0) < 1e-14
    # leverage target: binding leverage av/lambda = lev_bg
    assert abs(b["av_bg"] / b["lam"] - p["lev_bg"]) < 1e-12
    # bond/capital from exposure target
    assert abs(b["BK"] - p["exp_bg"] / (1 - p["exp_bg"])) < 1e-14
    # q_b normalized 1 <-> coupon reproduces R_B on BGP with Q_B=1
    R_B = (p["pi_mat"] + (1 - p["pi_mat"]) * (b["iota"] + 1.0)) / 1.0
    assert abs(R_B - b["R_B_bg"]) < 1e-12


def test_marginal_value_and_spread():
    # av RECURSION AND THE CLOSED-FORM BGP CREDIT SPREAD ARE CONSISTENT.
    p = base_params()
    b = solve_bgp(p)
    av, mu, psi = b["av_bg"], p["mu_bg"], p["psi"]
    # av = [(1-psi)+psi*av]/(1-mu)  <=>  av(1-mu) = (1-psi)+psi*av
    assert abs(av * (1 - mu) - ((1 - psi) + psi * av)) < 1e-13
    # spread = lambda*mu/E[Lambda_hat'] = mu*R/(lev*(1-mu))
    spread_alt = mu * p["R_bg"] / (p["lev_bg"] * (1 - mu))
    assert abs(b["RK_minus_R"] - spread_alt) < 1e-13
    # ~8bp ANNUALIZED (per-quarter decimal x4 x1e4 = annual bps); matches
    # Bocola's own ~8bp interbank-spread anchor
    annual_bps = b["RK_minus_R"] * 4 * 1e4
    assert 5.0 < annual_bps < 12.0, f"BGP spread {annual_bps:.2f} bps ann"


def test_jermann_normalization():
    # ON BGP: Q_K = 1/Phi'(x*) = 1  AND  Phi(x*) = x*  (adj_bg = 0).
    p = base_params()
    b = solve_bgp(p)
    xi, x = p["xi"], b["xstar"]
    a1, a2 = b["a1"], b["a2"]
    phi = a1 * x ** (1 - xi) + a2
    phi_prime = a1 * (1 - xi) * x ** (-xi)
    assert abs(phi - x) < 1e-12, "Phi(x*) != x* (adj_bg=0 violated)"
    assert abs(phi_prime - 1.0) < 1e-12, "Q_K != 1 (Phi'(x*) != 1)"


def test_labor_foc():
    # chi REPRODUCES l_bg FROM THE INTRATEMPORAL CONDITION AT THE BGP SHARES.
    p = base_params()
    b = solve_bgp(p)
    # chi*l^nu * c = (1-alpha)*Y/l ; with Y=1, c=cy, w=(1-alpha)/l
    l = ((1 - p["alpha"]) / (b["chi"] * b["cy"])) ** (1.0 / (1.0 + p["nu"]))
    assert abs(l - p["l_bg"]) < 1e-12


def test_delta_modes_bracket_target():
    # euler MODE HITS i/y EXACTLY; standard MODE (delta=0.025) OVERSHOOTS IT,
    # CONFIRMING THE DOCUMENTED TENSION (tiny mu => low delta => high K/Y).
    be = solve_bgp(delta_mode="euler")
    bs = solve_bgp(delta_mode="standard")
    assert abs(be["iy_implied"] - be["iy_bg"]) < 1e-12
    assert bs["iy_implied"] > bs["iy_bg"] + 0.05   # standard delta misses high
    assert be["KY"] > bs["KY"]                     # euler delta => higher K/Y


if __name__ == "__main__":
    test_targets_hit_exactly()
    test_marginal_value_and_spread()
    test_jermann_normalization()
    test_labor_foc()
    test_delta_modes_bracket_target()
    print("test_calibration: ALL PASSED")
