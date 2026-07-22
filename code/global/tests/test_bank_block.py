# BANK-BLOCK UNIT TESTS: SS FIXED POINT, BOND FOC IDENTITIES UNDER PRICED DEFAULT
# RISK, AND NO-ARBITRAGE AROUND A REALIZED DEFAULT EVENT.
import numpy as np

from common import get_ss, ss_input_paths
from bank import bank_backward, bank_forward


def _run_bank(cal, ss, def_price_D=None, def_real_D=None):
    # RUN THE BACKWARD AND FORWARD BANK PASSES ON CONSTANT SS INPUT PATHS.
    paths = ss_input_paths(cal, ss)
    bwd = bank_backward(paths["rk_D"], paths["rk_F"], paths["rdep_D"],
                        paths["rdep_F"], paths["p_path"], cal,
                        ss["ss_bank_D"], ss["ss_bank_F"],
                        def_price_D=def_price_D)
    b_D_D = np.full(cal["T"], cal["B_gov_D_ss"]) - bwd["b_D_F"]
    b_F_F = np.full(cal["T"], cal["B_gov_F_ss"]) - bwd["b_F_D"]
    fwd = bank_forward(paths["Kap_D"], paths["Kap_F"], paths["Q_D"], paths["Q_F"],
                       paths["rk_D"], paths["rk_F"], paths["rdep_D"], paths["rdep_F"],
                       paths["p_path"], b_D_D, b_F_F, bwd, cal,
                       ss["ss_bank_D"], ss["ss_bank_F"], def_real_D=def_real_D)
    return bwd, fwd


def test_fixed_point_at_ss():
    # AT SS INPUTS THE BANK BLOCK MUST REPRODUCE EVERY SS OBJECT.
    cal, ss = get_ss()
    bwd, fwd = _run_bank(cal, ss)
    assert np.max(np.abs(bwd["Q_bD"] - ss["Q_bD_ss"])) < 1e-12
    assert np.max(np.abs(fwd["n_D"] - ss["ss_bank_D"]["n_ss"])) < 1e-9
    assert np.max(np.abs(fwd["n_IC_D"] - fwd["n_D"])) < 1e-9
    assert np.max(np.abs(fwd["rb_D"] - ss["ss_bank_D"]["rb_dom_ss"])) < 1e-12


def test_priced_risk_foc_identity():
    # E_t[1+rb'] = 1 + rdep_t + lambda*mu_t/Omega AT EVERY t UNDER PRICED DEFAULT
    # RISK — THE EXACT GK BOND FOC WITH AN EXPECTED HAIRCUT.
    cal, ss = get_ss()
    T = cal["T"]
    defp = 0.03 * 0.9 ** np.arange(T)
    bwd, fwd = _run_bank(cal, ss, def_price_D=defp)
    db, rec = cal["delta_b_D"], cal["recovery_rate_D"]
    r = cal["r_dep_D_target"]
    lb = cal["lambda_bD_D"]
    max_foc = 0.0
    for t in range(T - 1):
        surv_e = 1 - defp[t + 1] * (1 - rec)
        exp_ret = surv_e * (db + (1 - db) * bwd["Q_bD"][t + 1]) / bwd["Q_bD"][t] - 1
        req = r + lb * bwd["mu_D"][t] / bwd["Omega_D"][t]
        max_foc = max(max_foc, abs(exp_ret - req))
    assert max_foc < 1e-12, f"bond FOC violated: {max_foc:.2e}"

    # Risk-only economics: MTM loss at t=0, no realized loss afterwards
    assert bwd["Q_bD"][0] < ss["Q_bD_ss"]
    assert fwd["n_D"][0] < ss["ss_bank_D"]["n_ss"]
    # Ex-post bond returns exceed the required return while risk persists
    assert fwd["rb_D"][1] > ss["ss_bank_D"]["rb_dom_ss"]


def test_realized_default_no_arbitrage():
    # A ONE-OFF HAIRCUT AT t=5 IS FULLY PRICED AT t=4, SO THE REALIZED RETURN IN
    # THE DEFAULT PERIOD EQUALS THE REQUIRED RETURN.
    cal, ss = get_ss()
    T = cal["T"]
    defp = np.zeros(T); defp[5] = 0.20
    bwd, fwd = _run_bank(cal, ss, def_price_D=defp, def_real_D=defp)
    db, rec = cal["delta_b_D"], cal["recovery_rate_D"]
    r = cal["r_dep_D_target"]
    lb = cal["lambda_bD_D"]
    req_4 = r + lb * bwd["mu_D"][4] / bwd["Omega_D"][4]
    assert abs(fwd["rb_D"][5] - req_4) < 1e-12
    # Anticipation dip: price at t=4 below SS, recovery after the event
    assert bwd["Q_bD"][4] < bwd["Q_bD"][6] < ss["Q_bD_ss"] + 1e-12


def test_no_sunspot_equals_no_shock():
    # def_price = 0 MUST GIVE PATHS IDENTICAL TO THE NO-SHOCK RUN, GUARDING
    # AGAINST ANY HIDDEN SUNSPOT CHANNEL RE-APPEARING.
    cal, ss = get_ss()
    bwd0, fwd0 = _run_bank(cal, ss)
    bwdz, fwdz = _run_bank(cal, ss, def_price_D=np.zeros(cal["T"]))
    assert np.array_equal(bwd0["Q_bD"], bwdz["Q_bD"])
    assert np.array_equal(fwd0["n_D"], fwdz["n_D"])


if __name__ == "__main__":
    test_fixed_point_at_ss()
    test_priced_risk_foc_identity()
    test_realized_default_no_arbitrage()
    test_no_sunspot_equals_no_shock()
    print("test_bank_block: ALL PASSED")
