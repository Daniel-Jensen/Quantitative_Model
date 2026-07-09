"""General-equilibrium regression tests:
 1. the zero-shock transition is a fixed point of the full pipeline
    (including the endogenous debt integration), and
 2. Walras-redundant residuals stay at the numerical floor even when the
    sunspot moves the debt stock (the pre-fix leak was 1.6e-2)."""
import numpy as np

from common import get_ss, transition_residuals
from transition import solve_transition, solve_transition_ck

WALRAS_TOL = 2e-6   # SS household-grid floor is ~9e-7


def test_zero_shock_fixed_point():
    cal, ss = get_ss()
    T = cal["T"]
    out = solve_transition(ss, cal, np.full(T, cal["Z_ss_D"]),
                           np.full(T, cal["Z_ss_F"]), verbose=False)
    assert np.max(np.abs(out["Y_D"] / ss["ss_firm_D"]["Y_ss"] - 1)) < 1e-5
    assert np.max(np.abs(out["n_D"] / ss["ss_bank_D"]["n_ss"] - 1)) < 1e-5
    assert np.max(np.abs(out["b_gov_D"] / cal["B_gov_D_ss"] - 1)) < 1e-6
    res = transition_residuals(out, cal)
    assert res["goods_F"] < WALRAS_TOL, f"goods_F = {res['goods_F']:.2e}"
    assert res["goods_D"] < 1e-9
    return out


def test_walras_with_moving_debt():
    """Risk-only sunspot moves b_gov by construction; the goods-F residual
    must stay at the grid floor because banks absorb the true issuance."""
    cal, ss = get_ss()
    T = cal["T"]
    sun = 0.03 * 0.9 ** np.arange(T)
    out = solve_transition_ck(ss, cal, np.full(T, cal["Z_ss_D"]),
                              np.full(T, cal["Z_ss_F"]),
                              sunspot_D_path=sun, verbose=False)
    moved = np.max(np.abs(out["b_gov_D"] / cal["B_gov_D_ss"] - 1))
    assert moved > 1e-3, f"debt did not move ({moved:.2e}) — shock inactive?"
    res = transition_residuals(out, cal)
    assert res["goods_F"] < WALRAS_TOL, f"goods_F = {res['goods_F']:.2e}"
    assert res["goods_D"] < 1e-9
    assert res["dep_D"] < 1e-6 and res["dep_F"] < 1e-6
    # Bond market clearing against the true end-of-period stock
    assert np.max(np.abs(out["b_D_D"] + out["b_D_F"] - out["b_gov_D"])) < 1e-10
    return out


if __name__ == "__main__":
    test_zero_shock_fixed_point()
    test_walras_with_moving_debt()
    print("test_transition_walras: ALL PASSED")
