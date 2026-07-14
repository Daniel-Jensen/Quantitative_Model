"""Entry point: solve the two-country HANK-GK monetary union steady state,
a TFP-shock transition, and the centerpiece Cole-Kehoe / Bocola (2016)
sovereign-risk pass-through experiment; print diagnostics and save figures.

Sections are independent and toggleable:
  - flags below (RUN_TFP / RUN_SUNSPOT), or
  - command line:  python3 main.py tfp        (TFP shock only)
                   python3 main.py sunspot    (CK/Bocola experiment only)
                   python3 main.py            (everything)
The steady state always runs — both experiments need it.
"""

import os
import sys
import time
import numpy as np

from calibration import get_calibration
from steady_state import solve_steady_state
from transition import solve_transition
from risk_branch import solve_transition_ck_risk, bond_decomposition
from plots import (OUTDIR, plot_steady_state, plot_household_policies,
                   plot_irf, plot_default_irf)

# ── Run flags: flip to False to skip a section ───────────────────────────────
RUN_TFP     = True    # TFP shock in country D
RUN_SUNSPOT = False   # Cole-Kehoe sunspot + Bocola risk channel (centerpiece)

# ── Output flags ─────────────────────────────────────────────────────────────
PRINT_SS         = True   # steady-state tables
PRINT_TRANSITION = True   # transition residual diagnostics
PRINT_CK         = True   # Cole-Kehoe / Bocola diagnostics


# ═════════════════════════════════════════════════════════════════════════════
#  SECTION 0 — STEADY STATE (always runs; both experiments start from it)
# ═════════════════════════════════════════════════════════════════════════════

def run_steady_state(cal):
    print("=" * 65)
    print("  Two-country HANK-GK monetary union: steady-state solve")
    print("=" * 65)
    t0 = time.perf_counter()
    ss = solve_steady_state(cal)
    print(f"  [steady state]  {time.perf_counter() - t0:.1f}s")

    if PRINT_SS:
        _print_ss_table(ss, cal)

    plot_steady_state(ss, cal)
    plot_household_policies(ss, cal)
    print(f"\nFigures saved to {OUTDIR}")
    return ss


def _print_ss_table(ss, cal):
    bk_D = ss["ss_bank_D"];  bk_F = ss["ss_bank_F"]
    fm_D = ss["ss_firm_D"];  fm_F = ss["ss_firm_F"]

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


# ═════════════════════════════════════════════════════════════════════════════
#  SECTION 1 — TFP SHOCK (baseline real shock, no default risk)
# ═════════════════════════════════════════════════════════════════════════════

def run_tfp(ss, cal):
    print("\n" + "=" * 65)
    print("  TFP shock in country D: rho=0.8, shock=0.01")
    print("=" * 65)

    # 1% log-deviation shock, AR(1) decay
    rho_z, shock0 = 0.8, 0.01
    Z_D_path = cal["Z_ss_D"] * np.exp(shock0 * rho_z ** np.arange(cal["T"]))
    Z_F_path = np.full(cal["T"], cal["Z_ss_F"])

    t0 = time.perf_counter()
    out = solve_transition(ss, cal, Z_D_path, Z_F_path, verbose=False)
    print(f"  [TFP transition]  {time.perf_counter() - t0:.1f}s")
    if PRINT_TRANSITION:
        _print_transition_residuals(out, cal)

    plot_irf(out, ss, cal)
    print(f"\nFigures saved to {OUTDIR}")
    return out


# ═════════════════════════════════════════════════════════════════════════════
#  SECTION 2 — COLE-KEHOE SUNSPOT + BOCOLA RISK CHANNEL (centerpiece)
#  Risk is PRICED but never REALIZED (def_real ≡ 0): pure pass-through.
# ═════════════════════════════════════════════════════════════════════════════

def run_sunspot(ss, cal):
    bk_D = ss["ss_bank_D"];  bk_F = ss["ss_bank_F"]
    fm_D = ss["ss_firm_D"]

    print("\n" + "=" * 65)
    print("  Cole-Kehoe sunspot (risk-only, Bocola pass-through):")
    rho_sun, sun0 = 0.95, 0.01
    print(f"  peak default prob xi_0 = {sun0:.0%} q, rho = {rho_sun}")
    print("=" * 65)
    T = cal["T"]

    sunspot_D_path = sun0 * rho_sun ** np.arange(T)
    Z_flat_D = np.full(T, cal["Z_ss_D"])
    Z_flat_F = np.full(T, cal["Z_ss_F"])

    # Single RISK-ON run (Bocola risk channel via the representative default
    # branch).  The 1% sunspot is small enough for a cold Newton start; the
    # 3-step homotopy and the risk-off comparison run were removed 2026-07-13
    # (one model, one figure).  solve_transition_ck_risk's round 0 still
    # solves the risk-off base internally as the fixed-point starting path.
    t0 = time.perf_counter()
    out_ck = solve_transition_ck_risk(
        ss, cal, Z_flat_D, Z_flat_F,
        sunspot_D_path=sunspot_D_path,
        verbose=True,
    )
    print(f"  [risk-on, branch fixed point]  {time.perf_counter() - t0:.1f}s")
    if PRINT_TRANSITION:
        _print_transition_residuals(out_ck, cal)

    # ── Pass-through diagnostics (Bocola 2016 style) ──────────────────────────
    dec = bond_decomposition(out_ck, ss, cal)
    sov_spread_ann = dec["total_yield"]

    rdep_lag = np.concatenate([[cal["r_dep_D_target"]], out_ck["rdep_D"][:-1]])
    lend_on  = 4e4 * ((out_ck["rk_D"] - rdep_lag) - ss["rk_D_ss"])
    i_peak = int(np.argmax(sov_spread_ann))
    ip_l   = int(np.argmax(lend_on))

    if PRINT_CK:
        br = out_ck["branch"]
        W2 = 30
        def ck_row(label, val, note=""):
            note_str = f"  {note}" if note else ""
            return f"  {label:<{W2}} {val:>10}{note_str}"

        print(f"\n{'':─<72}")
        print(f"  Cole-Kehoe / Bocola pass-through  (xi_0={sun0:.0%}, rho={rho_sun})")
        print(f"  {'Statistic':<{W2}} {'Risk-on':>10}  Note")
        print(f"{'':─<72}")
        print(ck_row("Q_bD[0]  (% dev from SS)",
                     f"{(out_ck['Q_bD'][0]/ss['Q_bD_ss']-1)*100:+.2f}%",
                     "MTM repricing"))
        print(ck_row("n_D[0]  (% dev)",
                     f"{(out_ck['n_D'][0]/bk_D['n_ss']-1)*100:+.2f}%",
                     "no default on base path"))
        print(ck_row("n_F[0]  (% dev)",
                     f"{(out_ck['n_F'][0]/bk_F['n_ss']-1)*100:+.2f}%",
                     "contagion leg"))
        print(ck_row("Y_D trough (% dev)",
                     f"{np.min(out_ck['Y_D']/fm_D['Y_ss']-1)*100:+.3f}%"))
        print(ck_row("I_D[0]  (% dev)",
                     f"{(out_ck['I_D'][0]/fm_D['I_ss']-1)*100:+.3f}%"))
        print(ck_row("Sov spread peak (bps ann)",
                     f"{sov_spread_ann[i_peak]:+.0f}",
                     f"t={i_peak}: def {dec['defcomp'][i_peak]:+.0f} + risk {dec['risk'][i_peak]:+.0f} + liq {dec['liquidity'][i_peak]:+.0f}"))
        print(ck_row("Lending spread peak (bps)",
                     f"{lend_on[ip_l]:+.0f}"))
        print(ck_row("b_gov_D peak",
                     f"{np.max(out_ck['b_gov_D']):.3f}",
                     f"SS={cal['B_gov_D_ss']:.3f}; Tax peak dev {np.max(out_ck['Tax_D'])-ss['Tax_D_ss']:+.4f}"))
        print(ck_row("def_real ≡ 0 on base path",
                     str(np.all(out_ck['def_real_D'] == 0))))
        print(ck_row("[branch] n_D(0)/n_ss",
                     f"{br['n_D'][0]/bk_D['n_ss']:.3f}",
                     f"Y_D(0) dev {(br['Y_D'][0]/fm_D['Y_ss']-1)*100:+.2f}%  (feared default state)"))
        print(f"{'':─<72}")

    plot_default_irf(out_ck, ss, cal)
    print(f"\nFigures saved to {OUTDIR}")
    return out_ck


# ═════════════════════════════════════════════════════════════════════════════
#  Shared diagnostics
# ═════════════════════════════════════════════════════════════════════════════

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


def main():
    t0_total = time.perf_counter()
    os.makedirs(OUTDIR, exist_ok=True)
    cal = get_calibration()

    ss = run_steady_state(cal)

    if RUN_TFP:
        run_tfp(ss, cal)
    # if RUN_SUNSPOT:
    #     run_sunspot(ss, cal)

    print("\n" + "=" * 65)
    print(f"  TOTAL  {time.perf_counter() - t0_total:.1f}s")
    print("=" * 65)


if __name__ == "__main__":
    # CLI selection overrides the flags: `python3 main.py tfp`, `... sunspot`
    # args = {a.lower() for a in sys.argv[1:]}
    # if args:
    #     RUN_TFP     = bool(args & {"tfp", "all"})
    #     RUN_SUNSPOT = bool(args & {"sunspot", "ck", "all"})
    main()
