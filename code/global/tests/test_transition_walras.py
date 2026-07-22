# GENERAL-EQUILIBRIUM REGRESSION TESTS: THE ZERO-SHOCK FIXED POINT AND WALRAS
# REDUNDANCY UNDER A MOVING DEBT STOCK (THE PRE-FIX LEAK WAS 1.6e-2).
import numpy as np

from common import get_ss
from transition import solve_transition, market_residuals

WALRAS_TOL = 2e-6   # SS household-grid floor is ~9e-7


def test_zero_shock_fixed_point():
    # THE ZERO-SHOCK TRANSITION MUST BE A FIXED POINT OF THE FULL PIPELINE.
    cal, ss = get_ss()
    T = cal["T"]
    out = solve_transition(ss, cal, np.full(T, cal["Z_ss_D"]),
                           np.full(T, cal["Z_ss_F"]), verbose=False)
    assert np.max(np.abs(out["Y_D"] / ss["ss_firm_D"]["Y_ss"] - 1)) < 1e-5
    assert np.max(np.abs(out["n_D"] / ss["ss_bank_D"]["n_ss"] - 1)) < 1e-5
    assert np.max(np.abs(out["b_gov_D"] / cal["B_gov_D_ss"] - 1)) < 1e-6
    res = market_residuals(out, cal)
    assert res["goods_F"] < WALRAS_TOL, f"goods_F = {res['goods_F']:.2e}"
    assert res["goods_D"] < 1e-9
    return out


def test_walras_with_moving_debt():
    # WITH b_gov MOVING, goods_F MUST STAY AT THE GRID FLOOR, THE IC
    # COMPLEMENTARITY MUST HOLD, AND THE UNION DEPOSIT MARKET MUST CLEAR AT PARITY
    # WITH THE CROSS-BORDER DEPOSIT POSITION ALIVE.
    cal, ss = get_ss()
    T = cal["T"]
    pi = 0.03 * 0.9 ** np.arange(T)
    out = solve_transition(ss, cal, np.full(T, cal["Z_ss_D"]),
                           np.full(T, cal["Z_ss_F"]),
                           def_price_D=pi, verbose=False)
    moved = np.max(np.abs(out["b_gov_D"] / cal["B_gov_D_ss"] - 1))
    assert moved > 1e-3, f"debt did not move ({moved:.2e}) — shock inactive?"
    res = market_residuals(out, cal)
    # goods_F is THE leak detector for the union integration: an unassigned
    # RER valuation profit on the interbank position would show up here.
    assert res["goods_F"] < WALRAS_TOL, f"goods_F = {res['goods_F']:.2e}"
    assert res["goods_D"] < 1e-9
    assert res["dep_union"] < 1e-6, f"union clearing = {res['dep_union']:.2e}"
    assert res["uip"] < 1e-9, f"deposit UIP = {res['uip']:.2e}"
    # the cross-border deposit position must move (margin alive)
    assert np.max(np.abs(out["nfa_dep_D"])) > 1e-3, \
        "cross-border deposit position did not move — union market inactive?"
    # Occasionally-binding IC: both complementarity legs non-negative,
    # product at the numerical floor
    assert res["mu_min_D"] > -1e-9 and res["mu_min_F"] > -1e-9
    assert res["slack_min_D"] > -1e-7 and res["slack_min_F"] > -1e-7
    assert res["cap_D"] < 1e-8 and res["cap_F"] < 1e-8
    # Bond market clearing against the true end-of-period stock
    assert np.max(np.abs(out["b_D_D"] + out["b_D_F"] - out["b_gov_D"])) < 1e-10
    return out


if __name__ == "__main__":
    test_zero_shock_fixed_point()
    test_walras_with_moving_debt()
    print("test_transition_walras: ALL PASSED")
