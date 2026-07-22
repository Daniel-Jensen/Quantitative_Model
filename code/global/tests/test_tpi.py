# TPI TESTS: NESTING OF THE PRICED THREE-BRANCH KERNEL AND THE MECHANICAL PRICE
# FLOOR, THE CB-BUDGET-CLOSURE IDENTITY (THE DIRECT REGRESSION GUARD AGAINST THE
# HISTORICAL TPI-1 LEAK), AND THE DIRECTIONAL PREDICTIONS.
import numpy as np

from common import get_ss, ss_input_paths
from bank import bank_backward
from transition import solve_transition, market_residuals
from risk_branch import solve_transition_risk, solve_tpi_branch, bond_decomposition


def _bwd(cal, ss, **kw):
    # RUN bank_backward ON CONSTANT STEADY-STATE INPUT PATHS.
    paths = ss_input_paths(cal, ss)
    return bank_backward(paths["rk_D"], paths["rk_F"], paths["rdep_D"],
                         paths["rdep_F"], paths["p_path"], cal,
                         ss["ss_bank_D"], ss["ss_bank_F"], **kw)


def _risk_D(cal, ss, defp):
    # A NON-DEGENERATE TWO-BRANCH risk_D DICT FOR KERNEL-ALGEBRA TESTS.
    T = cal["T"]
    bk = ss["ss_bank_D"]; bkF = ss["ss_bank_F"]
    Om_nd_D = cal["beta_inter_D"] * (cal["f_D"] + (1 - cal["f_D"]) * bk["alpha_ss"])
    Om_nd_F = cal["beta_inter_F"] * (cal["f_F"] + (1 - cal["f_F"]) * bkF["alpha_ss"])
    return dict(
        pi=defp,
        Omega_d_D=np.full(T, Om_nd_D * 1.8), Omega_d_F=np.full(T, Om_nd_F * 1.4),
        rk_d_D=-0.05, rk_d_F=-0.02, Q_bD_d=0.85, Q_bF_d=1.0, p_d=1.1,
    )


def test_tpi_off_nests_baseline():
    # tpi_D=None MUST REPRODUCE THE TWO-BRANCH PASS BIT-FOR-BIT, AND s_tpi_D=None
    # MUST EQUAL s_tpi_D=zeros WHATEVER tpi_D/risk_D CONTAIN.
    # The mechanical channel is independent of the priced channel.
    cal, ss = get_ss()
    T = cal["T"]
    defp = 0.03 * 0.9 ** np.arange(T)
    risk_D = _risk_D(cal, ss, defp)

    base = _bwd(cal, ss, risk_D=risk_D)
    same = _bwd(cal, ss, risk_D=risk_D, tpi_D=None, s_tpi_D=None)
    for k in ("Q_bD", "Q_bF", "alpha_D", "mu_D", "b_D_F", "b_F_D", "cb_buy_D"):
        assert np.allclose(base[k], same[k], rtol=0, atol=1e-13), k

    # s_tpi_D=None vs s_tpi_D=zeros: equivalent, both leave cb_buy_D=0 and
    # Q_bD=Q_bD_free, whatever the (active, non-degenerate) tpi_D pricing is.
    active_tpi = dict(
        pi=np.full(T, 0.7),
        Omega_tpi_D=np.full(T, 0.9), Omega_tpi_F=np.full(T, 0.9),
        rk_tpi_d_D=-0.3, rk_tpi_d_F=-0.2, Q_bD_tpi_d=0.5, Q_bF_tpi_d=0.6, p_tpi_d=1.2,
    )
    none_s  = _bwd(cal, ss, risk_D=risk_D, tpi_D=active_tpi, s_tpi_D=None)
    zero_s  = _bwd(cal, ss, risk_D=risk_D, tpi_D=active_tpi, s_tpi_D=np.zeros(T))
    for k in ("Q_bD", "Q_bF", "cb_buy_D"):
        assert np.allclose(none_s[k], zero_s[k], rtol=0, atol=1e-14), k
    assert np.allclose(none_s["cb_buy_D"], 0.0, atol=1e-14)
    assert np.allclose(none_s["Q_bD"], none_s["Q_bD_free"], atol=1e-14)

    # s_tpi_D=None is also inert with NO tpi_D or risk_D at all (risk-neutral path)
    base_rn = _bwd(cal, ss, def_price_D=defp)
    same_rn = _bwd(cal, ss, def_price_D=defp, s_tpi_D=None)
    assert np.allclose(base_rn["Q_bD"], same_rn["Q_bD"], rtol=0, atol=1e-14)


def test_tpi_requires_risk_mode():
    # tpi_D WITHOUT risk_D MUST RAISE: TPI PRICING NESTS INSIDE THE TWO-BRANCH
    # KERNEL AND IS NOT DEFINED STANDALONE.
    cal, ss = get_ss()
    T = cal["T"]
    tpi_D = dict(pi=np.ones(T), Omega_tpi_D=np.ones(T), Omega_tpi_F=np.ones(T),
                rk_tpi_d_D=0.0, rk_tpi_d_F=0.0, Q_bD_tpi_d=1.0, Q_bF_tpi_d=1.0, p_tpi_d=1.0)
    try:
        _bwd(cal, ss, tpi_D=tpi_D)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_pi_tpi_one_nests_two_branch():
    # pi_tpi_D_path = 1 (BACKSTOP CERTAIN) COLLAPSES THE THREE-BRANCH KERNEL TO
    # THE TWO-BRANCH FORMULA EXACTLY, BECAUSE w_tpi IS ALGEBRAICALLY ZERO.
    cal, ss = get_ss()
    T = cal["T"]
    defp = 0.03 * 0.9 ** np.arange(T)
    risk_D = _risk_D(cal, ss, defp)

    base = _bwd(cal, ss, risk_D=risk_D)
    tpi_D = dict(
        pi=np.ones(T),   # backstop certain -> w_tpi = 0 exactly
        Omega_tpi_D=np.full(T, 42.0), Omega_tpi_F=np.full(T, 42.0),  # must be inert
        rk_tpi_d_D=-0.7, rk_tpi_d_F=-0.7, Q_bD_tpi_d=0.2, Q_bF_tpi_d=0.2, p_tpi_d=3.0,
    )
    three_branch = _bwd(cal, ss, risk_D=risk_D, tpi_D=tpi_D)
    for k in ("Q_bD", "Q_bF", "alpha_D", "mu_D", "b_D_F", "b_F_D"):
        assert np.allclose(base[k], three_branch[k], rtol=0, atol=1e-13), k


def test_q_floor_mechanism_direction():
    # THE FLOOR NEVER LOWERS THE PRICE, NEVER BUYS A NEGATIVE QUANTITY, AND DOES
    # TRIGGER SOMEWHERE UNDER A REALISTIC SOLVED RISK PREMIUM.
    # At the unshocked SS the true gap is zero everywhere, so any purchase there
    # is pure _CB_SMOOTH_EPS noise, orders of magnitude below a real intervention.
    cal, ss = get_ss()
    T = cal["T"]
    Z_D = np.full(T, cal["Z_ss_D"]); Z_F = np.full(T, cal["Z_ss_F"])
    pi = 0.02 * 0.85 ** np.arange(T)
    on = solve_transition_risk(ss, cal, Z_D, Z_F, pi_D_path=pi, verbose=False)
    risk_D = dict(pi=pi, **on["risk_D_inputs"])

    s_tpi = np.ones(T)
    withtpi = bank_backward(on["rk_D"], on["rk_F"], on["rdep_D"], on["rdep_F"], on["p"],
                            cal, ss["ss_bank_D"], ss["ss_bank_F"],
                            risk_D=risk_D, s_tpi_D=s_tpi)
    assert np.all(withtpi["Q_bD"] >= withtpi["Q_bD_free"] - 1e-12)
    assert np.all(withtpi["cb_buy_D"] >= -1e-12)
    stress_cb_buy = np.max(withtpi["cb_buy_D"])
    assert stress_cb_buy > 0, "floor never triggered — test setup too mild"

    off = _bwd(cal, ss, s_tpi_D=s_tpi)
    ss_noise = np.max(off["cb_buy_D"])
    assert ss_noise < 1e-3 * cal["B_gov_D_ss"], "SS smoothing noise too large"
    assert ss_noise < 0.05 * stress_cb_buy, \
        "SS smoothing noise not negligible relative to a real stress purchase"


def test_tpi_branch_solve_and_nesting():
    # solve_tpi_branch RETURNS A VALID TRANSITION, AND DEGENERATE BRANCH OBJECTS
    # COLLAPSE THE THREE-BRANCH KERNEL TO THE TWO-BRANCH VALUE.
    cal, ss = get_ss()
    T = cal["T"]
    Z_D = np.full(T, cal["Z_ss_D"]); Z_F = np.full(T, cal["Z_ss_F"])
    off = solve_transition(ss, cal, Z_D, Z_F, verbose=False)

    branch = solve_tpi_branch(off, ss, cal, tau=1, verbose=False)
    assert branch["Q_bD"][0] > 0
    assert np.all(np.isfinite(branch["Q_bD"]))

    defp = 0.03 * 0.9 ** np.arange(T)
    base = _bwd(cal, ss, def_price_D=defp)
    # degenerate branch Ω must equal the SS continuation Ω_nd (NOT beta_inter
    # alone — Ω_nd = beta_inter*(f+(1-f)*alpha_ss)), else the branch isn't
    # actually degenerate and this test proves nothing about the kernel algebra.
    Om_nd_D = cal["beta_inter_D"] * (cal["f_D"] + (1 - cal["f_D"]) * ss["ss_bank_D"]["alpha_ss"])
    Om_nd_F = cal["beta_inter_F"] * (cal["f_F"] + (1 - cal["f_F"]) * ss["ss_bank_F"]["alpha_ss"])
    risk_D = dict(
        pi=defp,
        Omega_d_D=np.full(T, Om_nd_D), Omega_d_F=np.full(T, Om_nd_F),
        rk_d_D=ss["rk_D_ss"], rk_d_F=ss["rk_F_ss"],
        Q_bD_d=base["Q_bD_ss_val"], Q_bF_d=base["Q_bF_ss_val"], p_d=ss["p_ss"],
    )
    tpi_D = dict(
        pi=1.0 - defp,   # arbitrary priced-reneging path
        Omega_tpi_D=np.full(T, Om_nd_D), Omega_tpi_F=np.full(T, Om_nd_F),
        rk_tpi_d_D=ss["rk_D_ss"], rk_tpi_d_F=ss["rk_F_ss"],
        Q_bD_tpi_d=base["Q_bD_ss_val"], Q_bF_tpi_d=base["Q_bF_ss_val"], p_tpi_d=ss["p_ss"],
        surv_tpi_d=1.0,
    )
    three = _bwd(cal, ss, def_price_D=None, risk_D=risk_D, tpi_D=tpi_D)
    # degenerate branches equal to the (zero-risk) SS continuation: exact
    # equality once risk has died out (tail of the path), close near the front
    assert np.allclose(three["Q_bD"][60:], base["Q_bD"][60:], atol=1e-9)


def test_cb_budget_closure():
    # THE LOAD-BEARING TEST: MECHANICAL TPI MUST NOT REOPEN THE HISTORICAL TPI-1
    # WALRAS LEAK. goods_D/goods_F must stay at the SAME tolerance as the no-TPI
    # baseline across a grid of psi_cb_D and activation dates, not a looser one.
    cal, ss = get_ss()
    T = cal["T"]
    Z_D = np.full(T, cal["Z_ss_D"]); Z_F = np.full(T, cal["Z_ss_F"])
    pi = 0.02 * 0.85 ** np.arange(T)
    off = solve_transition(ss, cal, Z_D, Z_F, def_price_D=pi, verbose=False)
    res_off = market_residuals(off, cal, ss=ss)
    assert res_off["goods_D"] < 5e-9
    assert res_off["goods_F"] < 2e-6

    for psi in (0.5, 1.0):
        for on_date in (0, 5):
            cal_i = dict(cal); cal_i["psi_cb_D"] = psi
            s_tpi = np.zeros(T); s_tpi[on_date:] = 1.0
            on = solve_transition(ss, cal_i, Z_D, Z_F, def_price_D=pi,
                                  s_tpi_D=s_tpi, verbose=False, y0=off["y_vec"])
            res_on = market_residuals(on, cal_i, ss=ss)
            assert res_on["goods_D"] < 5e-9, \
                f"psi_cb_D={psi} on_date={on_date}: goods_D={res_on['goods_D']:.2e}"
            assert res_on["goods_F"] < 2e-6, \
                f"psi_cb_D={psi} on_date={on_date}: goods_F={res_on['goods_F']:.2e}"
            # the CB remittance identity itself (trivial by construction, but
            # catches lag/timing typos): rem_cb_D = coupon+continuation(cb_lag) - cost(cb)
            surv = 1.0 - np.asarray(on["def_real_D"]) * (1.0 - cal_i["recovery_rate_D"])
            cb_lag = np.concatenate(([0.0], on["cb_buy_D"][:-1]))
            expected = (cal_i["delta_b_D"] * surv * cb_lag
                       + on["Q_bD"] * surv * (1 - cal_i["delta_b_D"]) * cb_lag
                       - on["Q_bD"] * on["cb_buy_D"])
            assert np.allclose(on["rem_cb_D"], expected, atol=1e-12)


def test_tpi_directional_predictions():
    # SAME SHOCK, TPI-ON VS TPI-OFF: SMALLER Q_bD DECLINE, SMALLER LENDING-SPREAD
    # PEAK, SHALLOWER Y_D TROUGH.
    # rem_cb_D must also be non-zero somewhere, which catches a silently-inert
    # implementation that would otherwise pass budget closure trivially.
    cal, ss = get_ss()
    T = cal["T"]
    Z_D = np.full(T, cal["Z_ss_D"]); Z_F = np.full(T, cal["Z_ss_F"])
    pi = 0.02 * 0.85 ** np.arange(T)
    off = solve_transition(ss, cal, Z_D, Z_F, def_price_D=pi, verbose=False)

    s_tpi = np.ones(T)
    on = solve_transition(ss, cal, Z_D, Z_F, def_price_D=pi, s_tpi_D=s_tpi,
                          verbose=False, y0=off["y_vec"])

    assert np.any(np.abs(on["rem_cb_D"]) > 1e-8), "CB rebate is uniformly zero"
    assert np.max(on["cb_buy_D"]) > 0
    dec_on  = bond_decomposition(on, ss, cal)
    dec_off = bond_decomposition(off, ss, cal)
    assert (ss["Q_bD_ss"] - on["Q_bD"][0]) < (ss["Q_bD_ss"] - off["Q_bD"][0]), \
        "TPI must shrink the peak MTM repricing relative to the no-TPI run"
    rdep_lag_on  = np.concatenate([[cal["r_dep_D_target"]], on["rdep_D"][:-1]])
    rdep_lag_off = np.concatenate([[cal["r_dep_D_target"]], off["rdep_D"][:-1]])
    lend_on  = np.max((on["rk_D"] - rdep_lag_on) - ss["rk_D_ss"])
    lend_off = np.max((off["rk_D"] - rdep_lag_off) - ss["rk_D_ss"])
    assert lend_on <= lend_off + 1e-10, "TPI must not raise the lending spread"
    assert np.min(on["Y_D"]) >= np.min(off["Y_D"]) - 1e-10, \
        "TPI must not deepen the output trough"

    res_on = market_residuals(on, cal, ss=ss)
    assert res_on["goods_D"] < 5e-9
    assert res_on["goods_F"] < 2e-6


def test_extract_init_state_carries_cb_buy_lag():
    # extract_init_state MUST PROPAGATE cb_buy_D_lag0 TO ANY BRANCH LAUNCHED
    # MID-PATH. Without it the branches silently drop the CB's legacy coupon and
    # continuation income at their own h=0, understating branch resources.
    from risk_branch import extract_init_state
    cal, ss = get_ss()
    T = cal["T"]
    Z_D = np.full(T, cal["Z_ss_D"]); Z_F = np.full(T, cal["Z_ss_F"])
    pi = 0.02 * 0.85 ** np.arange(T)
    s_tpi = np.ones(T)
    on = solve_transition(ss, cal, Z_D, Z_F, def_price_D=pi, s_tpi_D=s_tpi, verbose=False)
    assert np.max(on["cb_buy_D"]) > 0, "test setup: CB must actually hold something"

    for tau in (1, 5):
        init = extract_init_state(on, ss, cal, tau)
        assert "cb_buy_D_lag0" in init
        assert init["cb_buy_D_lag0"] == on["cb_buy_D"][tau - 1]


def test_full_three_way_budget_closure():
    # THE DEFAULT AND TPI-RENEGED BRANCHES PRICED SIMULTANEOUSLY — THE FULL OUTER
    # FIXED POINT — MUST CLOSE goods_D/goods_F AT THE SAME TOLERANCE AS EVERY
    # SIMPLER CONFIGURATION.
    cal, ss = get_ss()
    T = cal["T"]
    Z_D = np.full(T, cal["Z_ss_D"]); Z_F = np.full(T, cal["Z_ss_F"])
    pi_D = 0.01 * 0.95 ** np.arange(T)
    pi_tpi = 1.0 - 0.03 * 0.9 ** np.arange(T)
    s_tpi = np.ones(T)

    out = solve_transition_risk(ss, cal, Z_D, Z_F, pi_D_path=pi_D,
                                pi_tpi_D_path=pi_tpi, s_tpi_D_path=s_tpi,
                                verbose=False, max_rounds=12)
    assert out["branch"] is not None
    assert out["tpi_branch"] is not None
    res = market_residuals(out, cal, ss=ss)
    assert res["goods_D"] < 5e-9, f"goods_D = {res['goods_D']:.2e}"
    assert res["goods_F"] < 2e-6, f"goods_F = {res['goods_F']:.2e}"
    assert np.max(out["cb_buy_D"]) > 0
    assert np.any(np.abs(out["rem_cb_D"]) > 1e-8)


if __name__ == "__main__":
    test_tpi_off_nests_baseline()
    test_tpi_requires_risk_mode()
    test_pi_tpi_one_nests_two_branch()
    test_q_floor_mechanism_direction()
    test_tpi_branch_solve_and_nesting()
    test_cb_budget_closure()
    test_tpi_directional_predictions()
    test_extract_init_state_carries_cb_buy_lag()
    test_full_three_way_budget_closure()
    print("test_tpi: ALL PASSED")
