"""
Numerical verification: which mechanism drives the default-risk IRF in the modular model?

Tests two hypotheses:
  H1 (psi_bd):  state-dependent IC tightening — lbD_eff = lbD + psi_bd*xi (Bocola-Dovis)
  H2 (null):    psi_bd=0, no IC tightening; bond price drop is purely from default survival

Method:
  Run BD sunspot shock (xi_path: 10% peak, rho=0.85) twice:
    (A) baseline:   psi_bd_D=3.0  (current calibration)
    (B) psi off:    psi_bd_D=0.0  (H1 switched off)

  Metrics:
    - max |Y_D - Y_D_ss| / Y_D_ss      (output response)
    - max |Q_bD - Q_bD_ss| / Q_bD_ss   (bond-price response)
    - max |IC_spread_D|                  (spread tightening)

  If H1 is the mechanism: A shows large responses, B shows ~0 responses.
  If some other channel were active: B would still show nonzero responses.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import copy
from calibration   import get_calibration
from steady_state  import solve_steady_state
from transition    import solve_transition

T_PLOT = 20

def run_default_shock(psi_override, label):
    cal = get_calibration()
    cal["psi_bd_D"] = psi_override   # BD sunspot IC-tightening sensitivity (key used by bank.py)
    # psi_bd_F stays 0 — F-country has no default risk
    cal["T"] = 100

    ss = solve_steady_state(cal)
    Y_D_ss   = ss["ss_firm_D"]["Y_ss"]
    Q_bD_ss  = ss["Q_bD_ss"]

    rho_sun, sun0 = 0.85, 0.10
    sunspot_D_path = sun0 * rho_sun ** np.arange(cal["T"])

    out = solve_transition(
        ss, cal,
        np.full(cal["T"], cal["Z_ss_D"]),
        np.full(cal["T"], cal["Z_ss_F"]),
        sunspot_D_path=sunspot_D_path,
        verbose=False,
    )

    y_dev    = 100.0 * (np.asarray(out["Y_D"])[:T_PLOT]   / Y_D_ss  - 1.0)
    qb_dev   = 100.0 * (np.asarray(out["Q_bD"])[:T_PLOT]  / Q_bD_ss - 1.0)
    rdep_dev = 10000.0 * np.asarray(out["rdep_D"])[:T_PLOT]  # bps above 0

    print(f"\n{'='*60}")
    print(f"  {label}  (psi_bd_D = {psi_override})")
    print(f"{'='*60}")
    print(f"  max |Y_D  dev|   = {np.max(np.abs(y_dev)):8.4f} %")
    print(f"  max |Q_bD dev|   = {np.max(np.abs(qb_dev)):8.4f} %")
    print(f"  peak Y_D  dev t0 = {y_dev[0]:+8.4f} %")
    print(f"  peak Q_bD dev t0 = {qb_dev[0]:+8.4f} %")
    print(f"  peak rdep_D  bps = {rdep_dev[0]:+8.2f} bps")
    print(f"  Y_D  path (first {T_PLOT}q): {np.round(y_dev[:10], 3)}")
    print(f"  Q_bD path (first {T_PLOT}q): {np.round(qb_dev[:10], 3)}")

    # Residuals to confirm convergence
    goods_D = np.max(np.abs(out["Y_D"] - out["C_D"] - out["I_D"] - out["NX_D"]))
    dep_D   = np.max(np.abs(out["A_D"] - out["Dep_supply_D"]))
    print(f"  [residuals] goods_D={goods_D:.2e}  deposit_D={dep_D:.2e}")

    return y_dev, qb_dev


if __name__ == "__main__":
    print("\nVerification: mechanism behind default-risk IRF")
    print("Comparing Bocola-Dovis IC tightening (psi_bd=3) vs null (psi_bd=0)\n")

    y_psi3, qb_psi3 = run_default_shock(psi_override=3.0, label="BASELINE  psi_bd=3.0")
    y_psi0, qb_psi0 = run_default_shock(psi_override=0.0, label="PSI=0     psi_bd=0.0")

    print("\n" + "="*60)
    print("  VERDICT")
    print("="*60)
    ratio_y  = np.max(np.abs(y_psi3))  / max(np.max(np.abs(y_psi0)),  1e-12)
    ratio_qb = np.max(np.abs(qb_psi3)) / max(np.max(np.abs(qb_psi0)), 1e-12)
    print(f"  Y_D response ratio  (psi=3 / psi=0) = {ratio_y:.1f}x")
    print(f"  Q_bD response ratio (psi=3 / psi=0) = {ratio_qb:.1f}x")
    print()
    if ratio_y > 50 and ratio_qb > 50:
        print("  CONFIRMED: psi_bd (Bocola-Dovis IC tightening) is the sole active mechanism.")
        print("  With psi_bd=0 the sunspot shock produces negligible effects (<1% of psi_bd=3 response).")
    else:
        print("  WARNING: some other mechanism may be active even with psi_bd=0.")
        print(f"  psi=0 Y_D max = {np.max(np.abs(y_psi0)):.4f}%  (should be ~0)")
