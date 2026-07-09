"""Entry point: solve the two-country HANK-GK monetary union steady state,
a TFP-shock transition, and the centerpiece Cole-Kehoe / Bocola (2016)
sovereign-risk pass-through experiment; print diagnostics and save figures."""

import os
import time
import numpy as np

from calibration import get_calibration
from steady_state import solve_steady_state
from transition import solve_transition, solve_transition_ck
from risk_branch import solve_transition_ck_risk, bond_decomposition
from plots import (OUTDIR, plot_steady_state, plot_irf, plot_default_irf,
                   plot_risk_comparison)


def main():
    t0_total = time.perf_counter()
    os.makedirs(OUTDIR, exist_ok=True)
    cal = get_calibration()

    print("=" * 65)
    print("  Two-country HANK-GK monetary union: steady-state solve")
    print("=" * 65)
    t0 = time.perf_counter()
    ss = solve_steady_state(cal)
    print(f"  [steady state]  {time.perf_counter() - t0:.1f}s")

    bk_D = ss["ss_bank_D"]
    bk_F = ss["ss_bank_F"]
    fm_D = ss["ss_firm_D"]
    fm_F = ss["ss_firm_F"]

    print("\n── Country D (domestic / Greece) ──────────────────────────────")
    print(f"  rk_D_ss      = {ss['rk_D_ss']:.4%}  beta_D = {ss['beta_D_ss']:.6f}")
    print(f"  Kap_D_ss     = {ss['Kap_D_ss']:.4f}   Y_D = {fm_D['Y_ss']:.4f}  I_D = {fm_D['I_ss']:.4f}")
    print(f"  C_D_ss       = {ss['C_D_ss']:.4f}   A_D = {ss['A_D_ss']:.4f}")
    print(f"  w_D_ss       = {fm_D['w_ss']:.4f}   chi_D = {cal['chi_D']:.4f}")
    print(f"  n_D_ss       = {bk_D['n_ss']:.4f}   theta_D = {bk_D['theta_ss']:.4f}")
    print(f"  alpha_D_ss   = {bk_D['alpha_ss']:.4f}   mu_D_ss = {bk_D['mu_ss']:.6f}")
    print(f"  kappa_D_ss   = {bk_D['kappa_ss']:.4f}   phi_bdom = {bk_D['phi_bdom_ss']:.4f}  phi_bfor = {bk_D['phi_bfor_ss']:.4f}")
    print(f"  sov exposure = {ss['Q_bD_ss'] * ss['b_D_D_ss'] / bk_D['n_ss']:.3f}"
          f" of net worth  (Bocola GIPS fact: ≈0.93)")
    print(f"  Dep_supply_D = {bk_D['Dep_supply_ss']:.4f}   rb_dom_ss = {bk_D['rb_dom_ss']:.4%}")

    print("\n── Country F (foreign / Germany) ───────────────────────────────")
    print(f"  rk_F_ss      = {ss['rk_F_ss']:.4%}  beta_F = {ss['beta_F_ss']:.6f}")
    print(f"  Kap_F_ss     = {ss['Kap_F_ss']:.4f}   Y_F = {fm_F['Y_ss']:.4f}  I_F = {fm_F['I_ss']:.4f}")
    print(f"  C_F_ss       = {ss['C_F_ss']:.4f}   A_F = {ss['A_F_ss']:.4f}")
    print(f"  n_F_ss       = {bk_F['n_ss']:.4f}   theta_F = {bk_F['theta_ss']:.4f}")
    print(f"  kappa_F_ss   = {bk_F['kappa_ss']:.4f}   phi_bdom = {bk_F['phi_bdom_ss']:.4f}  phi_bfor = {bk_F['phi_bfor_ss']:.4f}")
    print(f"  Dep_supply_F = {bk_F['Dep_supply_ss']:.4f}   rb_dom_ss = {bk_F['rb_dom_ss']:.4%}")

    print("\n── Global / cross-border ────────────────────────────────────────")
    print(f"  p_ss          = {ss['p_ss']:.6f}  (real exchange rate; 1 = symmetric)")
    print(f"  Q_bD_ss       = {ss['Q_bD_ss']:.5f}   Q_bF_ss = {ss['Q_bF_ss']:.5f}")
    print(f"  b_D_D_ss      = {ss['b_D_D_ss']:.5f}   b_F_D_ss = {ss['b_F_D_ss']:.5f}  (D-bank holdings)")
    print(f"  b_F_F_ss      = {ss['b_F_F_ss']:.5f}   b_D_F_ss = {ss['b_D_F_ss']:.5f}  (F-bank holdings)")
    print(f"  F-bank share of D-debt = {ss['b_D_F_ss'] / cal['B_gov_D_ss']:.1%}  (contagion leg)")
    print(f"  face debt/annual GDP D = {cal['B_gov_D_ss'] / (4 * fm_D['Y_ss']):.1%}")

    print("\n── Steady-state residuals ──────────────────────────────────────")
    ic_resid_D = (bk_D["n_ss_IC"] - bk_D["n_ss_ACCUM"]) / bk_D["n_ss_ACCUM"]
    ic_resid_F = (bk_F["n_ss_IC"] - bk_F["n_ss_ACCUM"]) / bk_F["n_ss_ACCUM"]
    dep_resid_D = ss["A_D_ss"] - bk_D["Dep_supply_ss"]
    dep_resid_F = ss["A_F_ss"] - bk_F["Dep_supply_ss"]
    walras_D = fm_D["Y_ss"] - ss["C_D_ss"] - fm_D["I_ss"] - cal["G_D"]
    walras_F = fm_F["Y_ss"] - ss["C_F_ss"] - fm_F["I_ss"] - cal["G_F"]
    print(f"  IC resid D (n_IC/n_ACCUM - 1)   = {ic_resid_D:.2e}")
    print(f"  IC resid F (n_IC/n_ACCUM - 1)   = {ic_resid_F:.2e}")
    print(f"  deposit resid D (A - Dep_supply) = {dep_resid_D:.2e}")
    print(f"  deposit resid F (A - Dep_supply) = {dep_resid_F:.2e}")
    print(f"  Walras D (Y - C - I - G)         = {walras_D:.2e}")
    print(f"  Walras F (Y - C - I - G)         = {walras_F:.2e}  [diagnostic, not imposed]")

    plot_steady_state(ss, cal)
    print(f"\nFigures saved to {OUTDIR}")

    print("\n" + "=" * 65)
    print("  TFP shock in country D: rho=0.8, shock=0.01")
    print("=" * 65)
    rho_z, shock0 = 0.8, 0.01
    Z_D_path = cal["Z_ss_D"] * np.exp(shock0 * rho_z ** np.arange(cal["T"]))
    Z_F_path = np.full(cal["T"], cal["Z_ss_F"])

    t0 = time.perf_counter()
    out = solve_transition(ss, cal, Z_D_path, Z_F_path, verbose=False)
    print(f"  [TFP transition]  {time.perf_counter() - t0:.1f}s")
    _print_transition_residuals(out, cal)

    plot_irf(out, ss, cal)
    print(f"\nFigures saved to {OUTDIR}")

    # ── Centerpiece: Cole-Kehoe sunspot, Bocola (2016) pass-through ──────────
    # xi_t = probability lenders coordinate on no-rollover at t (crisis zone
    # active at SS).  PRICED into Q_bD; default never REALIZED (def_real=0).
    # Transmission: xi ↑ → Q_bD ↓ (expected-haircut pricing) → MTM loss on
    # legacy bonds → n_D ↓ → single-λ IC tightens → lending spread rk−rdep ↑
    # → I_D, Y_D ↓.  Government rolls over at depressed prices → b_gov ↑ →
    # Bohn taxes ↑ (beliefs worsen fundamentals, Cole-Kehoe).
    # Persistence matters: a 7-year bond's yield averages default risk over
    # its life, so Greek-scale spreads need a persistent sunspot (Bocola's
    # estimated s_t process is highly persistent).
    print("\n" + "=" * 65)
    print("  Cole-Kehoe sunspot (risk-only, Bocola pass-through):")
    rho_sun, sun0 = 0.95, 0.07
    print(f"  peak default prob xi_0 = {sun0:.0%} q, rho = {rho_sun}")
    print("=" * 65)
    T = cal["T"]
    sunspot_D_path = sun0 * rho_sun ** np.arange(T)
    Z_flat_D = np.full(T, cal["Z_ss_D"])
    Z_flat_F = np.full(T, cal["Z_ss_F"])

    t0 = time.perf_counter()
    # Step 1 — RISK-OFF (liquidity channel only), sunspot homotopy: warm-start
    # each step with the previous solution (large MTM repricing makes a cold
    # Newton start fragile).
    out_off, y_warm = None, None
    for scale in (0.25, 0.5, 1.0):
        out_off = solve_transition_ck(
            ss, cal, Z_flat_D, Z_flat_F,
            sunspot_D_path=scale * sunspot_D_path,
            verbose=False, y0=y_warm,
        )
        y_warm = out_off["y_vec"]
    print(f"  [risk-off base, 3-step homotopy]  {time.perf_counter() - t0:.1f}s")

    # Step 2 — RISK-ON: Bocola risk channel via the representative default
    # branch (two-branch expectations in the bank backward pass).
    t0 = time.perf_counter()
    out_ck = solve_transition_ck_risk(
        ss, cal, Z_flat_D, Z_flat_F,
        sunspot_D_path=sunspot_D_path,
        verbose=True, y0=y_warm,
    )
    print(f"  [risk-on, branch fixed point]  {time.perf_counter() - t0:.1f}s")
    _print_transition_residuals(out_ck, cal)

    # Comparable liquidity-only counterfactual: if the feasible priced event
    # is a PARTIAL restructuring (haircut ladder in risk_branch), re-solve the
    # risk-off path pricing the SAME event size so the on/off gap isolates
    # the risk channel rather than the event size.
    s_star = out_ck["branch"]["haircut_scale"]
    if s_star < 1.0:
        print(f"  priced event: partial restructuring, haircut "
              f"{s_star * (1 - cal['recovery_rate_D']):.0%} of face value "
              "(full PSI haircut infeasible without bank recap)")
        out_off = solve_transition_ck(
            ss, cal, Z_flat_D, Z_flat_F,
            sunspot_D_path=s_star * sunspot_D_path,
            verbose=False, y0=y_warm,
        )

    # ── Pass-through diagnostics (Bocola 2016 style) ──────────────────────────
    dec = bond_decomposition(out_ck, ss, cal)
    sov_spread_ann = dec["total_yield"]

    def lend_spread(out):
        rdep_lag = np.concatenate([[cal["r_dep_D_target"]], out["rdep_D"][:-1]])
        return 4e4 * ((out["rk_D"] - rdep_lag) - ss["rk_D_ss"])

    lend_on  = lend_spread(out_ck)
    lend_off = lend_spread(out_off)
    i_peak = int(np.argmax(sov_spread_ann))
    ip_l   = int(np.argmax(lend_on))
    risk_share = 1.0 - lend_off[ip_l] / lend_on[ip_l] if lend_on[ip_l] != 0 else np.nan

    print("\n── Cole-Kehoe / Bocola diagnostics (RISK CHANNEL ON) ─────────")
    print(f"  Q_bD[0]    = {out_ck['Q_bD'][0]:.5f}  (ss={ss['Q_bD_ss']:.5f})"
          f"  =>  MTM repricing = {(out_ck['Q_bD'][0] / ss['Q_bD_ss'] - 1) * 100:.2f}%"
          f"   [risk-off: {(out_off['Q_bD'][0] / ss['Q_bD_ss'] - 1) * 100:.2f}%]")
    print(f"  n_D[0]     dev = {(out_ck['n_D'][0] / bk_D['n_ss'] - 1) * 100:+.2f}%"
          f"   [risk-off: {(out_off['n_D'][0] / bk_D['n_ss'] - 1) * 100:+.2f}%]  (no default)")
    print(f"  n_F[0]     dev = {(out_ck['n_F'][0] / bk_F['n_ss'] - 1) * 100:+.2f}%  (contagion)")
    print(f"  Y_D trough = {np.min(out_ck['Y_D'] / fm_D['Y_ss'] - 1) * 100:+.3f}%"
          f"   [risk-off: {np.min(out_off['Y_D'] / fm_D['Y_ss'] - 1) * 100:+.3f}%]")
    print(f"  I_D[0]     = {(out_ck['I_D'][0] / fm_D['I_ss'] - 1) * 100:+.3f}%"
          f"   [risk-off: {(out_off['I_D'][0] / fm_D['I_ss'] - 1) * 100:+.3f}%]")
    print(f"  sov spread peak = {sov_spread_ann[i_peak]:+.0f} bps ann (t={i_peak}):"
          f"  default comp {dec['defcomp'][i_peak]:+.0f}"
          f"  + RISK PREMIUM {dec['risk'][i_peak]:+.0f}"
          f"  + liquidity {dec['liquidity'][i_peak]:+.0f}")
    print(f"  lending spread peak = {lend_on[ip_l]:+.0f} bps ann"
          f"   [risk-off: {lend_off[ip_l]:+.0f}]")
    print(f"  RISK-CHANNEL SHARE of lending-spread response = {risk_share:.0%}"
          f"   (Bocola 2016 estimate: up to 45%)")
    print(f"  b_gov_D peak = {np.max(out_ck['b_gov_D']):.3f}  (ss={cal['B_gov_D_ss']:.3f})"
          f"   Tax_D peak dev = {np.max(out_ck['Tax_D']) - ss['Tax_D_ss']:+.4f}")
    print(f"  def_real ≡ 0 on base path: {np.all(out_ck['def_real_D'] == 0)}")
    br = out_ck["branch"]
    print(f"  [branch] n_D(0)/n_ss = {br['n_D'][0] / bk_D['n_ss']:.3f}"
          f"   Y_D(0) dev = {(br['Y_D'][0] / fm_D['Y_ss'] - 1) * 100:+.2f}%"
          f"   (the default state bankers fear)")

    plot_default_irf(out_ck, ss, cal)
    plot_risk_comparison(out_ck, out_off, ss, cal)
    print(f"\nFigures saved to {OUTDIR}")

    print("\n" + "=" * 65)
    print(f"  TOTAL  {time.perf_counter() - t0_total:.1f}s")
    print("=" * 65)


def _print_transition_residuals(out, cal):
    cap_resid_D = np.max(np.abs(out["n_IC_D"] - out["n_D"]))
    cap_resid_F = np.max(np.abs(out["n_IC_F"] - out["n_F"]))
    dep_resid_D = np.max(np.abs(out["P_CES_D"] * out["A_D"] - out["Dep_supply_D"]))
    dep_resid_F = np.max(np.abs(out["P_CES_F"] * out["A_F"] - out["Dep_supply_F"]))
    goods_D = np.max(np.abs(out["Y_D"] - out["P_CES_D"] * out["C_D"] - out["I_D"]
                            - out["NX_D"] - cal["G_D"]))
    goods_F = np.max(np.abs(out["Y_F"] - out["P_CES_F"] * out["C_F"] - out["I_F"]
                            - out["NX_F"] - cal["G_F"]))
    print(f"  max|capital resid D| (n_IC − n)  = {cap_resid_D:.2e}")
    print(f"  max|capital resid F| (n_IC − n)  = {cap_resid_F:.2e}")
    print(f"  max|deposit resid D|             = {dep_resid_D:.2e}")
    print(f"  max|deposit resid F|             = {dep_resid_F:.2e}")
    print(f"  max|goods mkt D|                 = {goods_D:.2e}")
    print(f"  max|goods mkt F| [diagnostic]    = {goods_F:.2e}")


if __name__ == "__main__":
    main()
