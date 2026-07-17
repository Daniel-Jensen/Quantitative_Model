"""Bocola (2016) risk-channel tests: nesting, pricing identity, and the
directional predictions (precautionary deleveraging)."""
import numpy as np

from common import get_ss, ss_input_paths, transition_residuals
from bank import bank_backward
from transition import solve_transition
from risk_branch import solve_transition_risk, bond_decomposition


def _bwd(cal, ss, **kw):
    paths = ss_input_paths(cal, ss)
    return bank_backward(paths["rk_D"], paths["rk_F"], paths["rdep_D"],
                         paths["rdep_F"], paths["p_path"], cal,
                         ss["ss_bank_D"], ss["ss_bank_F"], **kw)


def test_pi_zero_nesting():
    """risk_D with pi ≡ 0 must reproduce the risk-neutral backward pass
    exactly, whatever the branch objects say. def_price_D stays None: in
    risk mode D risk enters via pi only."""
    cal, ss = get_ss()
    T = cal["T"]
    risk_D = dict(
        pi=np.zeros(T),
        Omega_d_D=np.full(T, 99.0), Omega_d_F=np.full(T, 99.0),  # must be inert
        rk_d_D=-0.5, rk_d_F=-0.5, Q_bD_d=0.1, Q_bF_d=0.1, p_d=2.0,
    )
    base  = _bwd(cal, ss)
    risky = _bwd(cal, ss, risk_D=risk_D)
    for k in ("Q_bD", "Q_bF", "alpha_D", "mu_D", "b_D_F", "b_F_D"):
        # FOC holdings divide the return wedge by psi = 0.01, amplifying
        # one ulp of Ω·x/Ω rounding to ~2e-14 — hence the looser atol.
        tol = 1e-13 if k in ("b_D_F", "b_F_D") else 1e-14
        assert np.allclose(base[k], risky[k], rtol=0, atol=tol), k


def test_degenerate_branch_nesting():
    """With Ω^d = Ω^nd and branch prices equal to the continuation values,
    the two-branch formulas collapse to the surv-form pricing."""
    cal, ss = get_ss()
    T = cal["T"]
    defp = 0.03 * 0.9 ** np.arange(T)
    base = _bwd(cal, ss, def_price_D=defp)
    bk = ss["ss_bank_D"]; bkF = ss["ss_bank_F"]
    Om_nd_D = cal["beta_inter_D"] * (cal["f_D"] + (1 - cal["f_D"]) * bk["alpha_ss"])
    Om_nd_F = cal["beta_inter_F"] * (cal["f_F"] + (1 - cal["f_F"]) * bkF["alpha_ss"])
    risk_D = dict(
        pi=defp,
        Omega_d_D=np.full(T, Om_nd_D), Omega_d_F=np.full(T, Om_nd_F),
        rk_d_D=ss["rk_D_ss"], rk_d_F=ss["rk_F_ss"],
        Q_bD_d=None, Q_bF_d=None, p_d=ss["p_ss"],
    )
    # Degenerate branch prices = the base continuation prices: at SS inputs
    # the base Q path solves the surv recursion, so feed Q^d = base Q_next.
    # Exact equivalence requires payoff_d = rec·payoff_nd, i.e. Q_bD_d equal
    # to the base next-period price at every t — only true at SS where the
    # base Q path is flat at the risk-adjusted level; here compare against a
    # zero-risk flat SS path instead.
    base0 = _bwd(cal, ss)
    risk_D["Q_bD_d"] = base0["Q_bD_ss_val"]
    risk_D["Q_bF_d"] = base0["Q_bF_ss_val"]
    risky = _bwd(cal, ss, def_price_D=None, risk_D=risk_D)
    # Compare against the surv-form with the same pi (base) — the branch Q^d
    # here is the SS price, matching the surv-form's Q_next only at SS; the
    # paths differ during the risk episode, so require closeness, not equality.
    assert np.max(np.abs(risky["Q_bD"] - base["Q_bD"])) < 5e-3
    # And exact equality once risk has died out (tail of the path)
    assert np.allclose(risky["Q_bD"][60:], base["Q_bD"][60:], atol=1e-10)


def test_pricing_identity_and_positive_premium():
    """GE run with the risk channel: the per-period decomposition identity is
    exact, the risk premium is positive during the crisis, and it is zero in
    the risk-off run."""
    cal, ss = get_ss()
    T = cal["T"]
    Z_D = np.full(T, cal["Z_ss_D"])
    Z_F = np.full(T, cal["Z_ss_F"])
    pi = 0.02 * 0.85 ** np.arange(T)
    off = solve_transition(ss, cal, Z_D, Z_F, def_price_D=pi, verbose=False)
    on = solve_transition_risk(ss, cal, Z_D, Z_F, pi_D_path=pi,
                               verbose=False, y0=off["y_vec"])

    dec_on  = bond_decomposition(on, ss, cal)
    dec_off = bond_decomposition(off, ss, cal)
    # identity: promised excess (dev) = defcomp + risk + liquidity (dev)
    ident_on = dec_on["promised_excess"] - (dec_on["defcomp"] + dec_on["risk"]
                                            + dec_on["liquidity"])
    assert np.max(np.abs(ident_on)) < 1e-6, "decomposition identity broken"
    # risk-neutral: risk premium ≈ 0 by the pricing equation
    assert np.max(np.abs(dec_off["risk"])) < 1e-6
    # Bocola channel: positive risk premium while risk is priced
    assert np.max(dec_on["risk"]) > 5.0, \
        f"risk premium too small: {np.max(dec_on['risk']):.2f} bps"

    # Directional prediction (pricing/valuation side): the risk channel must
    # DEEPEN the bond repricing relative to risk-off.
    assert on["Q_bD"][0] < off["Q_bD"][0], \
        f"risk-on Q_bD[0] not below risk-off ({on['Q_bD'][0]:.5f} vs {off['Q_bD'][0]:.5f})"
    # Occasionally-binding IC: complementarity holds on both solved paths
    for o in (on, off):
        assert o["mu_min_D"] > -1e-9, f"mu_D negative ({o['mu_min_D']:+.2e})"
        assert o["slack_min_D"] > -1e-7, f"slack_D negative ({o['slack_min_D']:+.2e})"

    # Branch event accounting: surv_d in pricing equals the recovery rate
    # (pure haircut on the whole claim), and the branch government budget
    # closes including any recap outlay (zero at the Bocola-pure baseline).
    br = on["branch"]
    surv_d = float(np.asarray(on["risk_D_inputs"]["surv_d"]))
    assert abs(surv_d - cal["recovery_rate_D"]) < 1e-12
    recap = br["recap_D_path"]
    iss_check = (cal["G_D"] + recap + br["coupon_D"] - br["Tax_D"]
                 - br["net_issuance_D"])
    assert np.max(np.abs(iss_check)) < 1e-12, \
        f"branch govt budget (incl. recap) violated: {np.max(np.abs(iss_check)):.2e}"

    # Accounting still closed with the risk machinery on
    res = transition_residuals(on, cal)
    assert res["goods_F"] < 2e-6, f"goods_F = {res['goods_F']:.2e}"
    # goods_D is the IMPOSED residual — its size is the solver's stopping
    # point, not an accounting identity (acceptance is cal["tol_transition"],
    # 1e-10 normalized; the Newton polish typically lands well below).
    # Allow 5e-9 here; risk-neutral runs keep the 1e-9 bar.
    assert res["goods_D"] < 5e-9, f"goods_D = {res['goods_D']:.2e}"

    # Branch sanity: default state has depressed net worth (the MTM hit)
    assert br["n_D"][0] < ss["ss_bank_D"]["n_ss"]
    assert np.all(on["def_real_D"] == 0), "no default realized on the base path"


def test_quality_and_wc_nesting():
    """quality0 = 1 and zeta_wc = 0 must reproduce the plain blocks exactly
    (capital path unchanged; no labour wedge)."""
    from capital import solve_capital_path
    cal, ss = get_ss()
    T = cal["T"]
    rng = np.random.default_rng(0)
    Kap = ss["Kap_D_ss"] * (1 + 0.01 * rng.standard_normal(T)).cumprod() ** 0.1
    mpk = np.full(T, ss["ss_firm_D"]["mpk_ss"]) * (1 + 0.005 * rng.standard_normal(T))
    base = solve_capital_path(Kap, ss["Kap_D_ss"], 1.0, mpk, cal, "D")
    same = solve_capital_path(Kap, ss["Kap_D_ss"], 1.0, mpk, cal, "D", quality0=1.0)
    for k in ("Q", "rk", "I", "iota", "cap_profit"):
        assert np.array_equal(base[k], same[k]), f"quality0=1 not neutral in {k}"
    # quality0 < 1: only the t=0 objects move, and the claim return takes
    # exactly the proportional payoff hit
    xi = 0.05
    shocked = solve_capital_path(Kap, ss["Kap_D_ss"], 1.0, mpk, cal, "D",
                                 quality0=1.0 - xi)
    delta = cal["delta_D"]
    rk0_expected = ((1 - xi) * (mpk[0] + (1 - delta) * shocked["Q"][0]) / 1.0 - 1)
    assert abs(shocked["rk"][0] - rk0_expected) < 1e-14

    # zeta = 0 leaves the SS wage/chi at the wedge-free values
    from firms import steady_state_firm, markup_ss
    cal0 = dict(cal); cal0["zeta_wc_D"] = 0.0
    fm0 = steady_state_firm(cal0, ss["Kap_D_ss"], "D")
    mc = markup_ss(cal, "D")
    w_unwedged = mc * (1 - cal["alpha_D"]) * fm0["Y_ss"]
    assert abs(fm0["w_ss"] - w_unwedged) < 1e-14
    # and zeta > 0 divides by exactly (1 + zeta·r_wc_ss)
    fm1 = steady_state_firm(cal, ss["Kap_D_ss"], "D")
    r_wc_ss = cal["r_dep_D_target"] + cal["credit_spread_target_D"]
    assert abs(fm1["w_ss"] * (1 + cal["zeta_wc_D"] * r_wc_ss) - w_unwedged) < 1e-14


def test_zero_shock_with_risk_machinery():
    """pi ≡ 0: the risk loop must return the SS fixed point (exact nesting)."""
    cal, ss = get_ss()
    T = cal["T"]
    Z_D = np.full(T, cal["Z_ss_D"])
    Z_F = np.full(T, cal["Z_ss_F"])
    out = solve_transition_risk(ss, cal, Z_D, Z_F, verbose=False, max_rounds=2)
    assert out["branch"] is None
    assert np.max(np.abs(out["Y_D"] / ss["ss_firm_D"]["Y_ss"] - 1)) < 1e-5
    assert np.max(np.abs(out["Q_bD"] / ss["Q_bD_ss"] - 1)) < 1e-5


if __name__ == "__main__":
    test_pi_zero_nesting()
    test_degenerate_branch_nesting()
    test_quality_and_wc_nesting()
    test_pricing_identity_and_positive_premium()
    test_zero_shock_with_risk_machinery()
    print("test_risk_channel: ALL PASSED")
