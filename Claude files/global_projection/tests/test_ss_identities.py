# STEADY-STATE IDENTITIES AGAINST GK/BOCOLA THEORY, AT MACHINE PRECISION.
import numpy as np

from common import get_ss


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
    # DEPOSIT AND GOODS MARKETS MUST CLEAR AT THE SOLVED STEADY STATE.
    cal, ss = get_ss()
    for c in ("D", "F"):
        bk = ss[f"ss_bank_{c}"]
        assert abs(ss[f"A_{c}_ss"] - bk["Dep_supply_ss"]) < 1e-6
        walras = (ss[f"ss_firm_{c}"]["Y_ss"] - ss[f"C_{c}_ss"]
                  - ss[f"ss_firm_{c}"]["I_ss"] - cal[f"G_{c}"])
        assert abs(walras) < 5e-6, f"SS goods market {c}: {walras:.2e}"


def test_government_stationary():
    # govt_transition AT SS PRICES WITH NO DEFAULT MUST KEEP DEBT CONSTANT.
    from blocks.government import govt_transition
    cal, ss = get_ss()
    T = cal["T"]
    for c in ("D", "F"):
        gov = govt_transition(cal, ss[f"gs_{c}"],
                              np.full(T, ss[f"Q_b{c}_ss"]), None, c)
        assert np.max(np.abs(gov["b_gov_eop"] - cal[f"B_gov_{c}_ss"])) < 1e-10
        assert np.max(np.abs(gov["Tax"] - ss[f"gs_{c}"]["Tax_ss"])) < 1e-10


def test_calibration_targets():
    # THE BOCOLA CALIBRATION ANCHORS DOCUMENTED IN calibration.py MUST BE HIT.
    cal, ss = get_ss()
    bk = ss["ss_bank_D"]
    # THE CALIBRATION TARGET IS EXPOSURE AS A SHARE OF BANK ASSETS: Bocola's exp^bg =
    # 7.6% (Table B1, 160/2093), which is what B_gov_D_ss is set to deliver. This used
    # to assert 0.7-1.1 on exposure/NET WORTH, a threshold left over from the 3.722
    # B_gov misreading ("93% of bank EQUITY" read as a debt/GDP ratio); it has been red
    # at 0.36 ever since B_gov was corrected -- and 0.36 is right, since 7.6% of assets
    # at leverage 5 IS 0.38 of net worth (Bocola's own q*b/n = 0.38).
    assets = bk["theta_ss"] * bk["n_ss"]
    exposure = ss["Q_bD_ss"] * ss["b_D_D_ss"] / assets
    assert 0.06 < exposure < 0.09, \
        f"D-sovereign exposure/assets = {exposure:.4f} (target 0.076, Bocola Table B1)"
    # the F bank is sz times bigger, so its OWN book carries b_D_F_ss/sz of the D bond
    sz = cal["size_F"] / cal["size_D"]
    assets_F = ss["ss_bank_F"]["theta_ss"] * ss["ss_bank_F"]["n_ss"]
    exp_F = ss["Q_bF_ss"] * ss["b_F_F_ss"] / assets_F
    assert 0.06 < exp_F < 0.11, f"F-sovereign exposure/assets = {exp_F:.4f}"
    # λ and ω_ent are calibrated to hit leverage and credit-spread targets
    for c in ("D", "F"):
        bkc = ss[f"ss_bank_{c}"]
        assert abs(bkc["theta_ss"] - cal[f"leverage_target_{c}"]) < 1e-6, \
            f"[{c}] leverage {bkc['theta_ss']:.4f} ≠ target {cal[f'leverage_target_{c}']}"
        spread = ss[f"rk_{c}_ss"] - cal[f"r_dep_{c}_target"]
        assert abs(spread - cal[f"credit_spread_target_{c}"]) < 1e-6, \
            f"[{c}] credit spread {spread:.5f} ≠ target {cal[f'credit_spread_target_{c}']}"


if __name__ == "__main__":
    test_bank_bellman_and_pricing()
    test_ss_market_clearing()
    test_government_stationary()
    test_calibration_targets()
    print("test_ss_identities: ALL PASSED")
