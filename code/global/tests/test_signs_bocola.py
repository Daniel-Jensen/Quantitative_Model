"""Acceptance criteria for the centerpiece risk-only experiment (Bocola 2016
pass-through).  An exogenous rise in the PRICED default probability π_t must:
  - reprice sovereign bonds down (MTM),
  - shrink both banks' net worth on impact (D directly, F via cross-holdings),
  - widen the lending spread rk − rdep,
  - contract output on impact,
  - raise the debt stock and Bohn taxes (beliefs worsen fundamentals),
while no default is ever realized.  These are the signs that were wrong in
the pre-rework model (sovereign risk was expansionary)."""
import numpy as np

from common import get_ss, transition_residuals
from transition import solve_transition


def test_risk_only_shock_signs():
    cal, ss = get_ss()
    T = cal["T"]
    pi = 0.03 * 0.90 ** np.arange(T)
    out = solve_transition(ss, cal, np.full(T, cal["Z_ss_D"]),
                           np.full(T, cal["Z_ss_F"]),
                           def_price_D=pi, verbose=False)

    Y_ss = ss["ss_firm_D"]["Y_ss"]
    n_ss_D = ss["ss_bank_D"]["n_ss"]
    n_ss_F = ss["ss_bank_F"]["n_ss"]

    # No default is ever realized — pure expectations experiment
    assert np.all(out["def_real_D"] == 0)

    # Bond repricing and bank losses on impact
    assert out["Q_bD"][0] < ss["Q_bD_ss"], "Q_bD must fall on impact"
    assert out["n_D"][0] < n_ss_D, "D-bank net worth must fall (MTM loss)"
    assert out["n_F"][0] < n_ss_F, "F-bank net worth must fall (contagion)"

    # Real contraction on impact (old acceptance criterion: positive = bug)
    assert out["Y_D"][0] < Y_ss, "Y_D[0] must fall"
    assert out["C_D"][0] < ss["C_D_ss"], "C_D[0] must fall (austerity + dividends)"

    # Lending spread rk − rdep widens
    rdep_lag = np.concatenate([[cal["r_dep_D_target"]], out["rdep_D"][:-1]])
    lend_spread = out["rk_D"] - rdep_lag
    assert np.max(lend_spread) > ss["rk_D_ss"] + 1e-4, "lending spread must widen"

    # Sovereign yield spread widens
    db = cal["delta_b_D"]
    y_sov = db / out["Q_bD"][0] - db
    y_ss = db / ss["Q_bD_ss"] - db
    assert y_sov > y_ss + 1e-4, "sovereign yield must rise"

    # Fiscal channel: rollover at depressed prices raises debt and taxes
    assert np.max(out["b_gov_D"]) > cal["B_gov_D_ss"] * 1.001, "debt must rise"
    assert np.max(out["Tax_D"]) > ss["Tax_D_ss"] + 1e-5, "Bohn taxes must rise"

    # Predetermined capital (Bocola eq. 6): with Z_0 = Z_ss and K_{-1} = K_ss,
    # impact output moves through hours ALONE — exact Cobb-Douglas identity.
    lhs = np.log(out["Y_D"][0] / Y_ss)
    rhs = (1 - cal["alpha_D"]) * np.log(out["N_D"][0] / 1.0)
    assert abs(lhs - rhs) < 1e-12, \
        f"impact-composition identity violated: {lhs:.3e} vs {rhs:.3e}"

    # Accounting stays closed while all of this happens
    res = transition_residuals(out, cal)
    assert res["goods_F"] < 2e-6, f"goods_F = {res['goods_F']:.2e}"


def test_headline_calibration_comovement():
    """Acceptance criterion of the comovement fix: at the HEADLINE experiment
    (pi = 1%·0.95^t, the main.py configuration) impact output and consumption
    must both fall. I_D[0] and the risk-on n ordering are reported, not
    asserted, on the first pass."""
    cal, ss = get_ss()
    T = cal["T"]
    pi = 0.01 * 0.95 ** np.arange(T)
    out = solve_transition(ss, cal, np.full(T, cal["Z_ss_D"]),
                           np.full(T, cal["Z_ss_F"]),
                           def_price_D=pi, verbose=False)
    Y_ss = ss["ss_firm_D"]["Y_ss"]
    assert out["Y_D"][0] < Y_ss, \
        f"headline Y_D[0] did not fall: {out['Y_D'][0] / Y_ss - 1:+.4%}"
    assert out["C_D"][0] < ss["C_D_ss"], \
        f"headline C_D[0] did not fall: {out['C_D'][0] / ss['C_D_ss'] - 1:+.4%}"
    i_dev = out["I_D"][0] / ss["ss_firm_D"]["I_ss"] - 1
    y_tr = np.min(out["Y_D"] / Y_ss - 1)
    print(f"  [headline 1%*0.95^t] Y_D[0] {out['Y_D'][0] / Y_ss - 1:+.4%}  "
          f"Y trough {y_tr:+.4%}  I_D[0] {i_dev:+.3%} (reported, not asserted)")


if __name__ == "__main__":
    test_risk_only_shock_signs()
    test_headline_calibration_comovement()
    print("test_signs_bocola: ALL PASSED")
