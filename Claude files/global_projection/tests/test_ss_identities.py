# STEADY-STATE IDENTITIES AGAINST GK/BOCOLA THEORY, AT MACHINE PRECISION.

from common import get_ss
from global_projection.blocks.trade import ces_price, import_demand, trade_balance, size_ratio


def test_bank_bellman_and_pricing():
    # BANK BELLMAN, KERNEL WEIGHTS, SINGLE lambda, AND IC-CONSISTENT BOND PRICING.
    cal, ss = get_ss()
    for c in ("D", "F"):
        bk = ss[f"ss_bank_{c}"]
        r = cal[f"r_dep_{c}_target"]
        db = cal[f"delta_b_{c}"]
        alpha, mu, Om = bk["alpha_ss"], bk["mu_ss"], bk["Omega_ss"]

        # Bellman: alpha = Omega (1+r) / (1 - mu)
        assert abs(Om * (1 + r) / (1 - mu) - alpha) < 1e-12
        # Omega = beta_inter [f + (1-f) alpha] (Bocola kernel: f = exit share,
        # weight 1-f = survival on the franchise value)
        assert abs(cal[f"beta_inter_{c}"] * (cal[f"f_{c}"] + (1 - cal[f"f_{c}"]) * alpha) - Om) < 1e-12
        # Franchise value exceeds outside option (needed for the risk channel)
        assert alpha > 1.0, f"[{c}] alpha_ss = {alpha:.4f} ≤ 1"
        # Single lambda (Bocola eq. 3): all divertabilities equal
        assert bk["lambda_K"] == bk["lambda_bD"] == bk["lambda_bF"]
        # Bond price = delta_b / (r + delta_b + lambda mu / Omega)
        spread = bk["lambda_bD"] * mu / Om
        assert abs(bk["Q_bdom_IC"] - db / (r + db + spread)) < 1e-12
        # Excess bond return = IC spread; with single lambda = capital spread
        assert abs(bk["rb_dom_ss"] - (r + spread)) < 1e-12
        assert abs(spread - (ss[f"rk_{c}_ss"] - r)) < 1e-10
        # IC-implied and accumulated net worth agree
        assert abs(bk["n_ss_IC"] / bk["n_ss_ACCUM"] - 1) < 1e-10


def test_ss_market_clearing():
    # DEPOSITS, BOTH GOODS MARKETS AND THE CURRENT ACCOUNT CLEAR AT THE SOLVED SS, AT ANY p.
    # Recomputed here from the blocks, not read back: NX from the trade block at the HA
    # consumption, net coupon income from the books at each bond's own yield.
    cal, ss = get_ss()
    p, sz = ss["p_ss"], size_ratio(cal)
    P = {c: ces_price(p, cal, c) for c in ("D", "F")}
    for c in ("D", "F"):
        # households hold deposits in basket units, banks issue them in goods
        gap = ss[f"A_{c}_ss"] * P[c] - ss[f"ss_bank_{c}"]["Dep_supply_ss"]
        assert abs(gap) < 1e-8, f"SS deposit market {c}: {gap:.2e}"
    IM_D = import_demand(p, ss["C_D_ss"], P["D"], cal, "D")
    IM_F = import_demand(p, ss["C_F_ss"], P["F"], cal, "F")
    NX = dict(zip(("D", "F"), trade_balance(p, IM_D, IM_F, cal)))
    for c in ("D", "F"):
        fm = ss[f"ss_firm_{c}"]
        walras = fm["Y_ss"] - P[c] * ss[f"C_{c}_ss"] - fm["I_ss"] - cal[f"G_{c}"] - NX[c]
        assert abs(walras) < 1e-8, f"SS goods market {c}: {walras:.2e}"
    yD = cal["r_dep_D_target"] + cal["credit_spread_target_D"]
    yF = cal["r_dep_F_target"] + cal["credit_spread_target_F"]
    nfi_D = yF * p * ss["Q_bF_ss"] * ss["b_F_D_ss"] - yD * ss["Q_bD_ss"] * ss["b_D_F_ss"]
    assert abs(NX["D"] + nfi_D) < 1e-8, f"D current account {NX['D'] + nfi_D:.2e}"
    assert abs(NX["F"] - nfi_D / (sz * p)) < 1e-8, "F current account"
    # chi pins N_ss = 1 at the real wage w/P_CES
    for c in ("D", "F"):
        assert abs(cal[f"chi_{c}"] - ss[f"ss_firm_{c}"]["w_ss"] / P[c]) < 1e-15


def test_calibration_targets():
    # THE DETERMINISTIC-SS INSTRUMENTS IN calibration.py MUST HOLD EXACTLY AT THE SOLVED SS.
    cal, ss = get_ss()
    sz = size_ratio(cal)
    # EXPOSURE IS OWN-SOVEREIGN BONDS AS A SHARE OF OWN-BANK ASSETS (EBA 2010: GR 14.4%,
    # DE 6.5%), and the debt stocks are solved to deliver it. These are the instruments; the
    # 14.4%/6.5% TARGETS are hit at the stochastic rest point (calibrate_stochastic.py).
    for c, own in (("D", ss["Q_bD_ss"] * ss["b_D_D_ss"]), ("F", ss["Q_bF_ss"] * ss["b_F_F_ss"])):
        bk = ss[f"ss_bank_{c}"]
        exposure = own / (bk["theta_ss"] * bk["n_ss"])
        assert abs(exposure - cal[f"sov_exposure_{c}"]) < 1e-12, \
            f"{c}-sovereign exposure/assets = {exposure:.6f} != {cal[f'sov_exposure_{c}']}"
    # each stock's foreign slice: F banks hold foreign_share_D of D's, D banks foreign_share_F of F's
    assert abs(ss["b_D_F_ss"] / cal["B_gov_D_ss"] - cal["foreign_share_D"]) < 1e-14
    assert abs(ss["b_F_D_ss"] / sz / cal["B_gov_F_ss"] - cal["foreign_share_F"]) < 1e-14
    # λ and ω_ent are calibrated to hit leverage and credit-spread targets
    for c in ("D", "F"):
        bkc = ss[f"ss_bank_{c}"]
        assert abs(bkc["theta_ss"] - cal[f"leverage_target_{c}"]) < 1e-6, \
            f"[{c}] leverage {bkc['theta_ss']:.4f} ≠ target {cal[f'leverage_target_{c}']}"
        spread = ss[f"rk_{c}_ss"] - cal[f"r_dep_{c}_target"]
        assert abs(spread - cal[f"credit_spread_target_{c}"]) < 1e-6, \
            f"[{c}] credit spread {spread:.5f} ≠ target {cal[f'credit_spread_target_{c}']}"


def test_bohn_debt_root():
    # THE CALIBRATED BOHN COEFFICIENT, AND A STABLE DEBT ROOT AT THE PRICE THE PERIOD MAP ISSUES AT.
    cal, ss = get_ss()
    for c in ("D", "F"):
        gs, db = ss[f"gs_{c}"], cal[f"delta_b_{c}"]
        assert gs["Q_B_ss"] == ss[f"Q_b{c}_ss"], f"[{c}] fiscal block not at the SS bond price"
        assert abs(4.0 * gs["gamma_tau"] - cal["bohn_annual"]) < 1e-15, f"[{c}] Bohn coefficient"
        root = (1.0 - db) + (db - gs["gamma_tau"]) / gs["Q_B_ss"]
        assert root < 1.0, f"[{c}] debt root {root:.10f} >= 1: the Bohn rule does not stabilise debt"


if __name__ == "__main__":
    test_bank_bellman_and_pricing()
    test_ss_market_clearing()
    test_calibration_targets()
    test_bohn_debt_root()
    print("test_ss_identities: ALL PASSED")
