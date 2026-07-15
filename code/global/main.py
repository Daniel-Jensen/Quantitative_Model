"""Entry point: IRFs of the two-country HANK-GK monetary union.

Pipeline (steady state always runs; both experiments start from it):
  1. TFP shock in country D            → output/tfp_irf.png
  2. Cole-Kehoe sunspot in country D   → output/default_irf.png
     (Bocola pass-through: sovereign-default risk is PRICED but never
      REALIZED — bond prices fall, bank net worth drops, spreads rise)
     optional: output/bond_decomposition.png (spread = default
     compensation + risk premium + liquidity premium)

Toggle the sections with the flags below.  Total runtime is dominated by
the transition solves (7T-unknown Newton each, see transition.py).
"""

import os
import time
import numpy as np

from calibration import get_calibration
from steady_state import solve_steady_state
from transition import solve_transition
from risk_branch import solve_transition_ck_risk, bond_decomposition
from plots import OUTDIR, plot_irf, plot_default_irf, plot_bond_decomposition

# ── Run flags ────────────────────────────────────────────────────────────────
RUN_TFP                 = True   # experiment 1: TFP shock in D
RUN_SUNSPOT             = True   # experiment 2: CK sunspot + Bocola risk channel
PLOT_BOND_DECOMPOSITION = True   # extra figure for experiment 2

# ── Shock parameters ─────────────────────────────────────────────────────────
TFP_SHOCK, TFP_RHO = 0.01, 0.9    # 1% impact, AR(1) decay
# Centerpiece: Bocola-style small persistent risk shock (user-selected
# 2026-07-15).  A 10%/0.9 configuration remains available as a stress demo —
# it prices a 64% cumulative default probability; run it only with the μ
# monitor in mind (transition.py warns if the IC multiplier goes negative).
SUN_SHOCK, SUN_RHO = 0.01, 0.95  # peak priced default prob (quarterly), decay


# ═════════════════════════════════════════════════════════════════════════════
#  Diagnostics: steady-state table, transition residuals, CK table
# ═════════════════════════════════════════════════════════════════════════════

def print_ss_table(ss, cal):
    # Steady-state moments, calibrated parameters, and residual checks.
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

    # Parameters calibrated inside the SS solve (these overwrite the warm
    # starts in calibration.py; see steady_state.py)
    print(f"{'':─<65}")
    print(f"  {'Calibrated in SS solve':<{W}} {'D':>10}  {'F':>10}  Pins / target")
    print(f"{'':─<65}")
    print(row("lambda (single, IC)", f"{cal['lambda_K_D']:.5f}", f"{cal['lambda_K_F']:.5f}", "leverage + credit-spread targets"))
    print(row("omega_ent (entrants)", f"{cal['omega_ent_D']:.6f}", f"{cal['omega_ent_F']:.6f}", "solved jointly with lambda"))
    print(row("Z_ss (rescaled)",  f"{cal['Z_ss_D']:.6f}", f"{cal['Z_ss_F']:.6f}", "pins Y_ss = 1"))
    print(row("chi (GHH)",        f"{cal['chi_D']:.4f}",  f"{cal['chi_F']:.4f}",  "pins N_ss = 1"))
    print(row("beta_ss",          f"{ss['beta_D_ss']:.6f}", f"{ss['beta_F_ss']:.6f}", "deposit-market clearing"))
    print(row("xr_for (bps ann)", f"{cal['excess_return_F_D_ss']*4e4:.2f}", f"{cal['excess_return_D_F_ss']*4e4:.2f}", "foreign-bond FOC anchor"))

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


def print_transition_residuals(out, cal):
    # goods_F is Walras-redundant (never imposed) — the honest accuracy check
    # (acceptance: goods_D ≤ 1e-9, goods_F ≤ 2e-6, see CLAUDE.md).
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


def print_ck_table(out, ss, cal):
    # Cole-Kehoe / Bocola pass-through diagnostics for the sunspot run.
    bk_D = ss["ss_bank_D"];  bk_F = ss["ss_bank_F"]
    fm_D = ss["ss_firm_D"]

    dec = bond_decomposition(out, ss, cal)
    sov_spread_ann = dec["total_yield"]
    rdep_lag = np.concatenate([[cal["r_dep_D_target"]], out["rdep_D"][:-1]])
    lend_on  = 4e4 * ((out["rk_D"] - rdep_lag) - ss["rk_D_ss"])
    i_peak = int(np.argmax(sov_spread_ann))
    ip_l   = int(np.argmax(lend_on))
    br     = out["branch"]

    W2 = 30
    def ck_row(label, val, note=""):
        note_str = f"  {note}" if note else ""
        return f"  {label:<{W2}} {val:>10}{note_str}"

    print(f"\n{'':─<72}")
    print(f"  Cole-Kehoe / Bocola pass-through  (xi_0={SUN_SHOCK:.0%}, rho={SUN_RHO})")
    print(f"  {'Statistic':<{W2}} {'Risk-on':>10}  Note")
    print(f"{'':─<72}")
    print(ck_row("Q_bD[0]  (% dev from SS)",
                 f"{(out['Q_bD'][0]/ss['Q_bD_ss']-1)*100:+.2f}%",
                 "MTM repricing"))
    print(ck_row("n_D[0]  (% dev)",
                 f"{(out['n_D'][0]/bk_D['n_ss']-1)*100:+.2f}%",
                 "no default on base path"))
    print(ck_row("n_F[0]  (% dev)",
                 f"{(out['n_F'][0]/bk_F['n_ss']-1)*100:+.2f}%",
                 "contagion leg"))
    print(ck_row("Y_D trough (% dev)",
                 f"{np.min(out['Y_D']/fm_D['Y_ss']-1)*100:+.3f}%"))
    print(ck_row("I_D[0]  (% dev)",
                 f"{(out['I_D'][0]/fm_D['I_ss']-1)*100:+.3f}%"))
    print(ck_row("Sov spread peak (bps ann)",
                 f"{sov_spread_ann[i_peak]:+.0f}",
                 f"t={i_peak}: def {dec['defcomp'][i_peak]:+.0f} + risk {dec['risk'][i_peak]:+.0f} + liq {dec['liquidity'][i_peak]:+.0f}"))
    print(ck_row("Lending spread peak (bps)",
                 f"{lend_on[ip_l]:+.0f}"))
    print(ck_row("b_gov_D peak",
                 f"{np.max(out['b_gov_D']):.3f}",
                 f"SS={cal['B_gov_D_ss']:.3f}; Tax peak dev {np.max(out['Tax_D'])-ss['Tax_D_ss']:+.4f}"))
    print(ck_row("def_real ≡ 0 on base path",
                 str(np.all(out['def_real_D'] == 0))))
    print(ck_row("min mu (IC mult) D / F",
                 f"{out['mu_min_D']:+.4f}",
                 f"F {out['mu_min_F']:+.4f}; negative = always-binding IC violated"))
    print(ck_row("[branch] n_D(0)/n_ss",
                 f"{br['n_D'][0]/bk_D['n_ss']:.3f}",
                 f"Y_D(0) dev {(br['Y_D'][0]/fm_D['Y_ss']-1)*100:+.2f}%  (feared default state)"))
    print(ck_row("[branch] rescue mode",
                 br.get("rescue_mode", "full"),
                 f"haircut_scale={br.get('haircut_scale', 1.0):.3f}"))
    print(f"{'':─<72}")


# ═════════════════════════════════════════════════════════════════════════════
#  Experiments
# ═════════════════════════════════════════════════════════════════════════════

def run_tfp(ss, cal):
    print("\n" + "=" * 65)
    print(f"  TFP shock in D: {TFP_SHOCK:.0%} on impact, rho = {TFP_RHO}")
    print("=" * 65)
    T = cal["T"]
    Z_D = cal["Z_ss_D"] * np.exp(TFP_SHOCK * TFP_RHO ** np.arange(T))
    Z_F = np.full(T, cal["Z_ss_F"])

    t0 = time.perf_counter()
    out = solve_transition(ss, cal, Z_D, Z_F, verbose=False)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_transition_residuals(out, cal)

    plot_irf(out, ss, cal)
    return out


def run_sunspot(ss, cal):
    print("\n" + "=" * 65)
    print(f"  Cole-Kehoe sunspot in D: peak default prob {SUN_SHOCK:.0%} "
          f"per quarter, rho = {SUN_RHO}")
    print("=" * 65)
    T = cal["T"]
    sunspot = SUN_SHOCK * SUN_RHO ** np.arange(T)
    Z_D = np.full(T, cal["Z_ss_D"])
    Z_F = np.full(T, cal["Z_ss_F"])

    t0 = time.perf_counter()
    out = solve_transition_ck_risk(ss, cal, Z_D, Z_F,
                                   sunspot_D_path=sunspot, verbose=True)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_transition_residuals(out, cal)
    print_ck_table(out, ss, cal)

    plot_default_irf(out, ss, cal)
    if PLOT_BOND_DECOMPOSITION:
        plot_bond_decomposition(out, ss, cal)
    return out


def main():
    t0 = time.perf_counter()
    os.makedirs(OUTDIR, exist_ok=True)
    cal = get_calibration()

    print("=" * 65)
    print("  Two-country HANK-GK monetary union: steady state")
    print("=" * 65)
    ss = solve_steady_state(cal)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_ss_table(ss, cal)

    if RUN_TFP:
        run_tfp(ss, cal)
    if RUN_SUNSPOT:
        run_sunspot(ss, cal)

    print(f"\nFigures saved to {OUTDIR}")
    print(f"TOTAL  {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
