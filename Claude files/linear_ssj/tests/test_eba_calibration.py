"""Regression checks for eba_calibration.py — run standalone:

    /opt/anaconda3/envs/ssj/bin/python code/test_eba_calibration.py
"""
import eba_calibration as E


def approx(a, b, tol):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def test_validate_deutsche_bank():
    E.validate()  # asserts DE017 == Deutsche Bank 2010


def test_gdp_parse():
    gdp = E.load_gdp()
    approx(gdp["EL"]["2010"], 223_590.5, 1.0)   # Greek 2010 nominal GDP, EUR m
    approx(gdp["DE"]["2010"], 2_615_260.0, 1.0)  # German 2010 nominal GDP, EUR m


def test_cross_holding_matrix():
    df = E.load_eba()
    de = E._de_banks(df)
    approx(E.sovereign_book(df, "GR", E.GR_BANKS), 54_447, 5)   # b_D_D
    approx(E.sovereign_book(df, "GR", de),          7_934, 5)   # b_D_F (contagion)
    approx(E.sovereign_book(df, "DE", E.GR_BANKS),    411, 5)   # b_F_D
    approx(E.sovereign_book(df, "DE", de),        315_313, 5)   # b_F_F


def test_model_targets():
    m = E.compute_moments()
    t = m["model_targets"]
    approx(t["n_inter_D"],    0.408, 0.005)
    approx(t["n_inter_F"],    0.175, 0.005)
    approx(t["phi_bD_F_ss"],  0.069, 0.002)
    approx(t["phi_bF_D_ss"],  0.018, 0.002)
    approx(t["phi_bD_D_ss"],  2.390, 0.02)
    approx(t["phi_bF_F_ss"],  2.758, 0.02)
    approx(t["B_supply_D_qgdp"], 1.116, 0.02)
    # rebuilt 2026-07-31
    approx(t["theta_D"],   5.511, 0.02)
    approx(t["theta_F"],   6.941, 0.02)
    approx(t["omega_K_D"], 0.117, 0.002)
    approx(t["omega_K_F"], 0.067, 0.002)
    approx(t["delta_b_D"], 0.0777, 0.0005)
    approx(t["delta_b_F"], 0.0568, 0.0005)


def test_ladder_exhausts_total():
    """The seven maturity buckets must sum to the reported MATURITY_TOTAL row."""
    df = E.load_eba()
    de = E._de_banks(df)
    for issuer, banks in [("GR", E.GR_BANKS), ("DE", de), ("GR", de), ("DE", E.GR_BANKS)]:
        lad = E.sovereign_ladder(df, issuer, banks)
        approx(sum(lad.values()), E.sovereign_book(df, issuer, banks), 1.0)


def test_duration_below_maturity():
    """Sanity: modified duration must be strictly below weighted-average maturity,
    and the gap must widen with the discount rate (12% GGB vs 2.9% Bund)."""
    m = E.compute_moments()
    gr, de = m["ladder"]["b_D_D"], m["ladder"]["b_F_F"]
    assert gr["modified_duration_y"] < gr["wavg_residual_maturity_y"]
    assert de["modified_duration_y"] < de["wavg_residual_maturity_y"]
    gap_gr = gr["wavg_residual_maturity_y"] - gr["modified_duration_y"]
    gap_de = de["wavg_residual_maturity_y"] - de["modified_duration_y"]
    assert gap_gr > gap_de, f"GGB gap {gap_gr:.2f} should exceed Bund gap {gap_de:.2f}"


def test_delta_b_roundtrip():
    for d in (0.02, 0.0568, 0.0777, 0.15, 0.30):
        approx(E.delta_b_for_duration(E.model_modified_duration(d)), d, 1e-9)


def test_balance_sheet_identity():
    """theta - phi_own - phi_cross must reproduce the measured capital book, i.e.
    omega_K is measured rather than a residual absorbing a theta assumption."""
    m = E.compute_moments()
    E.validate(moments=m)     # assertion 4 lives there
    t, raw = m["model_targets"], m["raw_EURm"]
    k_implied = (t["theta_D"] - t["phi_bD_D_ss"] - t["phi_bF_D_ss"]) * raw["CT1_D"]
    approx(k_implied, raw["K_bank_D"], 1.0)


def test_mtm_mechanical_channel():
    """The measured mechanical MTM hit must be an order of magnitude larger than
    what the pre-EBA placeholder (phi_own=0.25, delta_b=0.10) could generate.
    This is the finding that psi_lambda_B had been standing in for."""
    m = E.compute_moments()
    dnn = m["mtm"]["D"]["dNW_per_100bp"]
    assert -0.08 < dnn < -0.04, f"D mechanical MTM {dnn:.2%}/100bp out of range"
    pre_eba = -0.25 * E.model_modified_duration(0.10) * 0.01
    assert abs(dnn) > 5 * abs(pre_eba), \
        f"measured {dnn:.2%} vs pre-EBA placeholder {pre_eba:.2%}"


def test_adverse_scenario_not_used():
    """Guard: the 2011 adverse scenario excluded banking-book sovereign default,
    so it must never enter the moment set."""
    m = E.compute_moments()
    assert "adverse" in m["meta"]["NOT_used"]
    assert "adverse_scenario_CT1_depletion" in m["identification"]["deliberately_rejected"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")
