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
from plots import (OUTDIR, plot_steady_state, plot_household_policies,
                   plot_irf, plot_default_irf, plot_risk_comparison)


def main():
    # ── Output flags ─────────────────────────────────────────────────────────
    PRINT_SS         = True   # steady-state tables
    PRINT_TRANSITION = True   # TFP transition residuals
    PRINT_CK         = True   # Cole-Kehoe / Bocola diagnostics

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

    if PRINT_SS:
        W = 26  # label column width
        def row(label, vD, vF, note=""):
            note_str = f"  {note}" if note else ""
            return f"  {label:<{W}} {vD:>10}  {vF:>10}{note_str}"

        print(f"\n{'':─<65}")
        print(f"  {'Variable':<{W}} {'D (Greece)':>10}  {'F (Germany)':>10}  Note")
        print(f"{'':─<65}")

        # Real economy
        print(row("Y_ss",      f"{fm_D['Y_ss']:.4f}",   f"{fm_F['Y_ss']:.4f}",  "normalised (Z rescaled)"))
        print(row("K_ss",      f"{ss['Kap_D_ss']:.4f}", f"{ss['Kap_F_ss']:.4f}"))
        print(row("I_ss",      f"{fm_D['I_ss']:.4f}",   f"{fm_F['I_ss']:.4f}",  "= δ·K"))
        print(row("C_ss",      f"{ss['C_D_ss']:.4f}",   f"{ss['C_F_ss']:.4f}",  "HA aggregate"))
        print(row("A_ss (HH deposits)", f"{ss['A_D_ss']:.4f}", f"{ss['A_F_ss']:.4f}"))
        print(row("w_ss",      f"{fm_D['w_ss']:.4f}",   f"{fm_F['w_ss']:.4f}"))
        print(row("rk_ss (ann %)",  f"{ss['rk_D_ss']*400:.3f}", f"{ss['rk_F_ss']*400:.3f}", "target ≈ 1.8% ann"))
        print(row("beta_ss",   f"{ss['beta_D_ss']:.6f}", f"{ss['beta_F_ss']:.6f}"))
        print(row("chi (GHH)",  f"{cal['chi_D']:.4f}",  f"{cal['chi_F']:.4f}",  "pinned to N_ss = 1"))

        # Financial intermediary
        print(f"{'':─<65}")
        print(row("n_ss (net worth)",  f"{bk_D['n_ss']:.4f}",     f"{bk_F['n_ss']:.4f}"))
        print(row("theta (leverage)",  f"{bk_D['theta_ss']:.4f}", f"{bk_F['theta_ss']:.4f}", "target 4–6 (GK11)"))
        print(row("kappa (K/n)",       f"{bk_D['kappa_ss']:.4f}", f"{bk_F['kappa_ss']:.4f}"))
        print(row("phi_bdom (dom bond/n)", f"{bk_D['phi_bdom_ss']:.4f}", f"{bk_F['phi_bdom_ss']:.4f}"))
        print(row("phi_bfor (for bond/n)", f"{bk_D['phi_bfor_ss']:.4f}", f"{bk_F['phi_bfor_ss']:.4f}"))
        sov_D = ss['Q_bD_ss'] * ss['b_D_D_ss'] / bk_D['n_ss']
        sov_F = ss['Q_bF_ss'] * ss['b_F_F_ss'] / bk_F['n_ss']
        print(row("Q·b_dom / n",  f"{sov_D:.3f}", f"{sov_F:.3f}", "Bocola GIPS target ≈ 0.93"))
        print(row("alpha (V/n)",  f"{bk_D['alpha_ss']:.4f}", f"{bk_F['alpha_ss']:.4f}", "franchise value"))
        print(row("mu (IC mult)", f"{bk_D['mu_ss']:.6f}", f"{bk_F['mu_ss']:.6f}"))
        print(row("Dep_supply",   f"{bk_D['Dep_supply_ss']:.4f}", f"{bk_F['Dep_supply_ss']:.4f}"))
        print(row("rb_dom (ann %)", f"{bk_D['rb_dom_ss']*400:.3f}", f"{bk_F['rb_dom_ss']*400:.3f}", "= rdep + IC spread"))

        # Sovereign / cross-border
        print(f"{'':─<65}")
        print(row("p_ss (RER)",    f"{ss['p_ss']:.6f}", "—", "1 = symmetric SS"))
        print(row("Q_bD_ss",       f"{ss['Q_bD_ss']:.5f}", f"{ss['Q_bF_ss']:.5f}", "IC-consistent prices"))
        print(row("b_D_D / b_F_D (D-bank)", f"{ss['b_D_D_ss']:.4f}", f"{ss['b_F_D_ss']:.4f}", "dom / for holdings"))
        print(row("b_F_F / b_D_F (F-bank)", f"{ss['b_F_F_ss']:.4f}", f"{ss['b_D_F_ss']:.4f}"))
        print(row("F-bank share of D-debt", f"{ss['b_D_F_ss']/cal['B_gov_D_ss']:.1%}", "—", "contagion leg, target 20%"))
        print(row("B_gov / 4Y (debt/GDP)", f"{cal['B_gov_D_ss']/(4*fm_D['Y_ss']):.1%}", f"{cal['B_gov_F_ss']/(4*fm_F['Y_ss']):.1%}", "target ≈ 93%"))

        # Residuals
        ic_resid_D  = (bk_D["n_ss_IC"] - bk_D["n_ss_ACCUM"]) / bk_D["n_ss_ACCUM"]
        ic_resid_F  = (bk_F["n_ss_IC"] - bk_F["n_ss_ACCUM"]) / bk_F["n_ss_ACCUM"]
        dep_resid_D = ss["A_D_ss"] - bk_D["Dep_supply_ss"]
        dep_resid_F = ss["A_F_ss"] - bk_F["Dep_supply_ss"]
        walras_D    = fm_D["Y_ss"] - ss["C_D_ss"] - fm_D["I_ss"] - cal["G_D"]
        walras_F    = fm_F["Y_ss"] - ss["C_F_ss"] - fm_F["I_ss"] - cal["G_F"]
        print(f"{'':─<65}")
        print(f"  {'Residual':<{W}} {'D':>10}  {'F':>10}  Threshold")
        print(f"{'':─<65}")
        print(row("IC  (n_IC/n_ACCUM − 1)", f"{ic_resid_D:.2e}",  f"{ic_resid_F:.2e}",  "≤ 1e-9"))
        print(row("Deposit (A − Dep)",       f"{dep_resid_D:.2e}", f"{dep_resid_F:.2e}", "≤ 1e-9"))
        print(row("Walras (Y−C−I−G)",        f"{walras_D:.2e}",   f"{walras_F:.2e}",    "F = diagnostic only"))
        print(f"{'':─<65}")

    plot_steady_state(ss, cal)
    plot_household_policies(ss, cal)
    print(f"\nFigures saved to {OUTDIR}")

    ### TFP SHOCK TRANSITION (baseline, no default) ───────────────────────────────

    print("\n" + "=" * 65)
    print("  TFP shock in country D: rho=0.8, shock=0.01")
    print("=" * 65)

    #Generating a shock path 
    rho_z, shock0 = 0.8, 0.01
    Z_D_path = cal["Z_ss_D"] * np.exp(shock0 * rho_z ** np.arange(cal["T"])) #this makes it in logs
    Z_F_path = np.full(cal["T"], cal["Z_ss_F"])

    t0 = time.perf_counter()

    out = solve_transition(ss, cal, Z_D_path, Z_F_path, verbose=False)
    print(f"  [TFP transition]  {time.perf_counter() - t0:.1f}s")
    if PRINT_TRANSITION:
        _print_transition_residuals(out, cal)

    plot_irf(out, ss, cal)
    print(f"\nFigures saved to {OUTDIR}")

    ### DEFAULT SHOCK TRANSITION (baseline, no default) ───────────────────────────────

    print("\n" + "=" * 65)
    print("  Cole-Kehoe sunspot (risk-only, Bocola pass-through):")
    rho_sun, sun0 = 0.95, 0.07
    print(f"  peak default prob xi_0 = {sun0:.0%} q, rho = {rho_sun}")
    print("=" * 65)
    T = cal["T"]

    #Generate a sunspot path for the default shock
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
    if PRINT_TRANSITION:
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

    if PRINT_CK:
        br = out_ck["branch"]
        W2 = 30
        def ck_row(label, risk_on, risk_off="", note=""):
            off_str  = f"  {risk_off:>10}" if risk_off else ""
            note_str = f"  {note}" if note else ""
            return f"  {label:<{W2}} {risk_on:>10}{off_str}{note_str}"

        print(f"\n{'':─<72}")
        print(f"  Cole-Kehoe / Bocola pass-through  (xi_0={sun0:.0%}, rho={rho_sun})")
        print(f"  {'Statistic':<{W2}} {'Risk-on':>10}  {'Risk-off':>10}  Note")
        print(f"{'':─<72}")
        print(ck_row("Q_bD[0]  (% dev from SS)",
                     f"{(out_ck['Q_bD'][0]/ss['Q_bD_ss']-1)*100:+.2f}%",
                     f"{(out_off['Q_bD'][0]/ss['Q_bD_ss']-1)*100:+.2f}%",
                     "MTM repricing"))
        print(ck_row("n_D[0]  (% dev)",
                     f"{(out_ck['n_D'][0]/bk_D['n_ss']-1)*100:+.2f}%",
                     f"{(out_off['n_D'][0]/bk_D['n_ss']-1)*100:+.2f}%",
                     "no default on base path"))
        print(ck_row("n_F[0]  (% dev)",
                     f"{(out_ck['n_F'][0]/bk_F['n_ss']-1)*100:+.2f}%",
                     note="contagion leg"))
        print(ck_row("Y_D trough (% dev)",
                     f"{np.min(out_ck['Y_D']/fm_D['Y_ss']-1)*100:+.3f}%",
                     f"{np.min(out_off['Y_D']/fm_D['Y_ss']-1)*100:+.3f}%"))
        print(ck_row("I_D[0]  (% dev)",
                     f"{(out_ck['I_D'][0]/fm_D['I_ss']-1)*100:+.3f}%",
                     f"{(out_off['I_D'][0]/fm_D['I_ss']-1)*100:+.3f}%"))
        print(ck_row("Sov spread peak (bps ann)",
                     f"{sov_spread_ann[i_peak]:+.0f}",
                     note=f"t={i_peak}: def {dec['defcomp'][i_peak]:+.0f} + risk {dec['risk'][i_peak]:+.0f} + liq {dec['liquidity'][i_peak]:+.0f}"))
        print(ck_row("Lending spread peak (bps)",
                     f"{lend_on[ip_l]:+.0f}",
                     f"{lend_off[ip_l]:+.0f}"))
        print(ck_row("Risk-channel share",
                     f"{risk_share:.0%}",
                     note="Bocola 2016 estimate: up to 45%"))
        print(ck_row("b_gov_D peak",
                     f"{np.max(out_ck['b_gov_D']):.3f}",
                     note=f"SS={cal['B_gov_D_ss']:.3f}; Tax peak dev {np.max(out_ck['Tax_D'])-ss['Tax_D_ss']:+.4f}"))
        print(ck_row("def_real ≡ 0 on base path",
                     str(np.all(out_ck['def_real_D'] == 0))))
        print(ck_row("[branch] n_D(0)/n_ss",
                     f"{br['n_D'][0]/bk_D['n_ss']:.3f}",
                     note=f"Y_D(0) dev {(br['Y_D'][0]/fm_D['Y_ss']-1)*100:+.2f}%  (feared default state)"))
        print(f"{'':─<72}")

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

    



