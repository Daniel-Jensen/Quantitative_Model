"""Regression checks for the country-size asymmetry — run standalone:

    /opt/anaconda3/envs/ssj/bin/python -m pytest code/test_cross_border_units.py -v

THE CONVENTION: every F-side variable is PER F CAPITA and O(1); every D-side
variable is a D aggregate (size_D == 1). `size_F` (= 11.697, Germany/Greece 2010
GDP) appears in exactly four blocks — the only places the two countries meet:
domestic_bond_clearing, trade_balance, external_account_D, global_goods_mkt.

WHY: with both countries normalised to Y_ss = 1 the model could not match the EBA
portfolio-composition moment (phi_bD_F = 0.0075) and the market-structure moment
(foreigners hold 12.7% of the bank-held Greek stock) at the same time. Joint
consistency needs n_F/n_D = 8.85 against the model's 0.761 — a gap that is exactly
the GDP ratio. Matching one forced a ~11x error in the other.

These tests are cheap (no model solve) and guard the arithmetic only. The pipeline's
printed over-identifying check (steady_state._apply_portfolio_targets) is what
confirms both moments hold in a SOLVED steady state.
"""
import json
import os

import calibration as C

_MOMENTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "data", "eba_moments.json")


def approx(a, b, tol):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def _raw():
    with open(_MOMENTS) as fh:
        return json.load(fh)["raw_EURm"]


def test_size_ratio_matches_eurostat():
    # Germany / Greece 2010 nominal GDP: 2,615,260 / 223,590.5
    approx(C.load_eba_size_ratio(), 11.697, 1e-3)
    approx(C.get_calibration()["size_F"], 11.697, 1e-3)


def test_foreign_shares_are_pure_ratios():
    """Immune to the GDP normalisation — that is what makes them the check."""
    fs = C.load_eba_foreign_shares()
    approx(fs["D"], 7_933.604 / 62_380.699, 1e-6)   # 12.72% of the Greek bank-held stock
    approx(fs["F"], 410.741 / 315_723.883, 1e-8)    # 0.13% of the Bund bank-held stock
    assert fs["D"] > 0.12, "Greek stock should be substantially foreign-held"
    assert fs["F"] < 0.01, "Bund stock should be almost entirely domestically held"


def test_both_eba_moments_hold_jointly():
    """The whole point of size_F: composition AND market structure, together.

    Pre-size-asymmetry no calibration could do both — matching phi_bD_F forced the
    foreign share to 1.25%, and matching the foreign share forced phi_bD_F to ~0.086.

    NOTE this runs on the INITIAL GUESS, where q_b_D = 0.83 rather than its solved
    0.969. Since b_D_F = phi*n_F/q, the low q inflates the face value and the
    foreign share reads ~15% here against 12.7% in the data. So the composition
    moment is asserted exactly (it is exact by construction at any q) while market
    structure is asserted only to be in the right regime. The precise match is
    checked on the SOLVED steady state by the pipeline's printed over-identifying
    check in steady_state._apply_portfolio_targets — that is the real gate.
    """
    cal = C.get_calibration()
    size_F = cal["size_F"]
    eba = C.load_eba_targets()
    fs = C.load_eba_foreign_shares()

    # Composition: b_D_F is per F capita, so phi is a per-capita ratio. Exact.
    phi_bD_F = cal["q_b_D"] * cal["b_D_F"] / cal["n_inter_F"]
    approx(phi_bD_F, eba["phi_bD_F_ss"], 1e-12)

    # Market structure: aggregate the F holding before comparing to the D stock.
    # The defect's signature was 1.25%; anything in double digits means size_F is
    # doing its job. Bracketed generously because q_b_D is a guess at this stage.
    share = size_F * cal["b_D_F"] / (cal["b_D_D"] + size_F * cal["b_D_F"])
    assert 0.08 < share < 0.20, f"foreign share {share:.4f} outside the plausible band"
    assert share > 8.0 * 0.0125, "size_F is not aggregating — this is the old defect"
    # And it must bracket the measured moment from the correct side of the old error.
    assert share > fs["D"] * 0.8


def test_bond_clearing_weights_are_consistent():
    """calibration's residual own-holdings must clear at the weights the model uses."""
    cal = C.get_calibration()
    size_F = cal["size_F"]
    # domestic_bond_clearing: b_D_D = b_gov_D - size_F*b_D_F, b_F_F = b_gov_F - b_F_D/size_F
    approx(cal["b_D_D"] + size_F * cal["b_D_F"], cal["B_supply_D"], 1e-10)
    approx(cal["b_F_F"] + cal["b_F_D"] / size_F, cal["B_supply_F"], 1e-10)


def test_home_bias_balances_bilateral_trade():
    """A single shared omega is inconsistent with size asymmetry.

    trade_balance uses NX_D = size_F*IM_F - p*IM_D. At p = 1 and C_D ~ C_F this
    is zero only if size_F*(1-omega_F) = (1-omega_D). At the old shared
    omega = 0.85 the D-good export leg came out size_F times too large.
    """
    cal = C.get_calibration()
    size_F = cal["size_F"]
    approx(size_F * (1.0 - cal["omega_F"]), 1.0 - cal["omega_D"], 1e-12)
    # The larger country is the more closed one.
    assert cal["omega_F"] > cal["omega_D"]
    approx(cal["omega_D"], 0.85, 1e-12)     # unchanged from the long-standing value


def test_pre_eba_branch_is_unweighted():
    """size_F = 1 keeps the pre-EBA placeholder calibration bit-exact."""
    live = C.EBA_CALIBRATION
    try:
        # get_calibration reads the module attribute at call time, so patching it
        # here is enough — see the calibration_override note in CLAUDE.md.
        C.EBA_CALIBRATION = False
        cal = C.get_calibration()
        assert cal["size_F"] == 1.0
        approx(cal["omega_F"], cal["omega_D"], 1e-12)
    finally:
        C.EBA_CALIBRATION = live
