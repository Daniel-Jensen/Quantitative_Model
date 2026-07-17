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


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\nAll {len(fns)} tests passed.")
