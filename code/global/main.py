"""Entry point: solve the two-country HANK-GK monetary union steady state
and a TFP-shock transition path, print diagnostics, and save figures."""

import os
import time
import numpy as np

from calibration import get_calibration
from steady_state import solve_steady_state
from transition import solve_transition, solve_transition_ck, solve_transition_bd
from plots import OUTDIR, plot_steady_state, plot_irf, plot_default_irf


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
    print(f"  Dep_supply_D = {bk_D['Dep_supply_ss']:.4f}   rb_dom_ss = {bk_D['rb_dom_ss']:.4%}")

    print("\n── Country F (foreign / Germany) ───────────────────────────────")
    print(f"  rk_F_ss      = {ss['rk_F_ss']:.4%}  beta_F = {ss['beta_F_ss']:.6f}")
    print(f"  Kap_F_ss     = {ss['Kap_F_ss']:.4f}   Y_F = {fm_F['Y_ss']:.4f}  I_F = {fm_F['I_ss']:.4f}")
    print(f"  C_F_ss       = {ss['C_F_ss']:.4f}   A_F = {ss['A_F_ss']:.4f}")
    print(f"  w_F_ss       = {fm_F['w_ss']:.4f}   chi_F = {cal['chi_F']:.4f}")
    print(f"  n_F_ss       = {bk_F['n_ss']:.4f}   theta_F = {bk_F['theta_ss']:.4f}")
    print(f"  alpha_F_ss   = {bk_F['alpha_ss']:.4f}   mu_F_ss = {bk_F['mu_ss']:.6f}")
    print(f"  kappa_F_ss   = {bk_F['kappa_ss']:.4f}   phi_bdom = {bk_F['phi_bdom_ss']:.4f}  phi_bfor = {bk_F['phi_bfor_ss']:.4f}")
    print(f"  Dep_supply_F = {bk_F['Dep_supply_ss']:.4f}   rb_dom_ss = {bk_F['rb_dom_ss']:.4%}")

    print("\n── Global / cross-border ────────────────────────────────────────")
    print(f"  p_ss          = {ss['p_ss']:.6f}  (real exchange rate; 1 = symmetric)")
    print(f"  Q_bD_ss       = {ss['Q_bD_ss']:.5f}   Q_bF_ss = {ss['Q_bF_ss']:.5f}")
    print(f"  b_D_D_ss      = {ss['b_D_D_ss']:.5f}   b_F_D_ss = {ss['b_F_D_ss']:.5f}  (D-bank holdings)")
    print(f"  b_F_F_ss      = {ss['b_F_F_ss']:.5f}   b_D_F_ss = {ss['b_D_F_ss']:.5f}  (F-bank holdings)")
    print(f"  excess_ret FD = {cal['excess_return_F_D_ss']:.4e}   excess_ret DF = {cal['excess_return_D_F_ss']:.4e}")

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

    cap_resid_D = np.max(np.abs(out["n_IC_D"] - out["n_D"]))
    cap_resid_F = np.max(np.abs(out["n_IC_F"] - out["n_F"]))
    dep_resid_D = np.max(np.abs(out["P_CES_D"] * out["A_D"] - out["Dep_supply_D"]))
    dep_resid_F = np.max(np.abs(out["P_CES_F"] * out["A_F"] - out["Dep_supply_F"]))
    goods_D = np.max(np.abs(out["Y_D"] - out["P_CES_D"] * out["C_D"] - out["I_D"] - out["NX_D"] - cal["G_D"]))
    goods_F = np.max(np.abs(out["Y_F"] - out["P_CES_F"] * out["C_F"] - out["I_F"] - out["NX_F"] - cal["G_F"]))
    print(f"  max|capital resid D| (n_IC − n)  = {cap_resid_D:.2e}")
    print(f"  max|capital resid F| (n_IC − n)  = {cap_resid_F:.2e}")
    print(f"  max|deposit resid D|             = {dep_resid_D:.2e}")
    print(f"  max|deposit resid F|             = {dep_resid_F:.2e}")
    print(f"  max|goods mkt D|                 = {goods_D:.2e}")
    print(f"  max|goods mkt F| [diagnostic]    = {goods_F:.2e}")

    plot_irf(out, ss, cal)
    print(f"\nFigures saved to {OUTDIR}")

    # ── Bocola-Dovis sunspot shock ────────────────────────────────────────────
    # xi_t = run-probability shock, AR(1): xi_t = xi_0 * rho^t
    # Transmission: xi_{t+1} → IC tightens (lbD_eff = lbD + psi_bd * xi)
    #             → Q_bD falls → government issues more bonds → b_gov rises
    #             → beliefs transform into fundamentals (Bocola-Dovis mechanism)
    print("\n" + "=" * 65)
    print("  Bocola-Dovis sunspot shock in D: rho=0.85, peak xi=0.10")
    print(f"  psi_bd_D={cal['psi_bd_D']:.1f}  =>  peak IC spread ≈"
          f" {cal['psi_bd_D'] * 0.10 * ss['ss_bank_D']['mu_ss'] / ss['ss_bank_D']['Omega_ss']:.4f}"
          f" (quarterly)")
    print("=" * 65)
    rho_sun, sun0 = 0.85, 0.10
    sunspot_D_path = sun0 * rho_sun ** np.arange(cal["T"])

    t0 = time.perf_counter()
    out_bd = solve_transition_bd(
        ss, cal,
        np.full(cal["T"], cal["Z_ss_D"]),
        np.full(cal["T"], cal["Z_ss_F"]),
        sunspot_D_path=sunspot_D_path,
        verbose=True,
    )
    print(f"  [BD transition]  {time.perf_counter() - t0:.1f}s")

    cap_resid_D = np.max(np.abs(out_bd["n_IC_D"] - out_bd["n_D"]))
    cap_resid_F = np.max(np.abs(out_bd["n_IC_F"] - out_bd["n_F"]))
    dep_resid_D = np.max(np.abs(out_bd["P_CES_D"] * out_bd["A_D"] - out_bd["Dep_supply_D"]))
    dep_resid_F = np.max(np.abs(out_bd["P_CES_F"] * out_bd["A_F"] - out_bd["Dep_supply_F"]))
    goods_D_bd  = np.max(np.abs(out_bd["Y_D"] - out_bd["P_CES_D"] * out_bd["C_D"] - out_bd["I_D"]
                                - out_bd["NX_D"] - cal["G_D"]))
    goods_F_bd  = np.max(np.abs(out_bd["Y_F"] - out_bd["P_CES_F"] * out_bd["C_F"] - out_bd["I_F"]
                                - out_bd["NX_F"] - cal["G_F"]))
    print(f"\n── Bocola-Dovis diagnostics ──────────────────────────────")
    print(f"  Q_bD[0]   = {out_bd['Q_bD'][0]:.5f}  (ss={ss['Q_bD_ss']:.5f})"
          f"  =>  price drop = {(out_bd['Q_bD'][0]/ss['Q_bD_ss']-1)*100:.2f}%")
    print(f"  b_gov_D peak = {np.max(out_bd['b_gov_D']):.4f}  (ss={cal['B_gov_D_ss']:.4f})")
    print(f"  Tax_D[0] = {out_bd['Tax_D'][0]:.5f}  (ss Tax increases via Bohn rule)")
    print(f"  n_D[0]   = {out_bd['n_D'][0]:.5f}  (ss={ss['ss_bank_D']['n_ss']:.5f})")
    print(f"  Y_D[0]   = {out_bd['Y_D'][0]:.5f}  (ss={ss['ss_firm_D']['Y_ss']:.5f})")
    print(f"  def_D peak = {np.max(out_bd['def_D']):.4f}  (1=fundamental default crossed)")
    print(f"  max|capital resid D|             = {cap_resid_D:.2e}")
    print(f"  max|capital resid F|             = {cap_resid_F:.2e}")
    print(f"  max|deposit resid D|             = {dep_resid_D:.2e}")
    print(f"  max|deposit resid F|             = {dep_resid_F:.2e}")
    print(f"  max|goods mkt D|                 = {goods_D_bd:.2e}")
    print(f"  max|goods mkt F| [diagnostic]    = {goods_F_bd:.2e}")

    plot_default_irf(out_bd, ss, cal, out_bd["sunspot_D"], out_bd["def_D"])
    print(f"\nFigures saved to {OUTDIR}")

    print("\n" + "=" * 65)
    print(f"  TOTAL  {time.perf_counter() - t0_total:.1f}s")
    print("=" * 65)

if __name__ == "__main__":
    main()
