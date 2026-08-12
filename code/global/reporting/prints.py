# CONSOLE REPORTING: STEADY-STATE TABLE. The projection experiments
# (solver_recursive/) print their own IRF tables inline; this keeps output
# formatting out of the model code.


def banner(text, width=65):
    # FULL-WIDTH SECTION HEADER.
    print("\n" + "=" * width)
    print(f"  {text}")
    print("=" * width)


def _rule(width):
    # HORIZONTAL TABLE RULE.
    print(f"{'':─<{width}}")


def _row(label, v1, v2="", note="", label_w=26):
    # ONE TABLE LINE: LABEL, ONE OR TWO VALUE COLUMNS, OPTIONAL TRAILING NOTE.
    tail = f"  {note}" if note else ""
    return f"  {label:<{label_w}} {v1:>10}  {v2:>10}{tail}"


def print_ss_table(ss, cal):
    # STEADY-STATE MOMENTS, CALIBRATED PARAMETERS, AND RESIDUAL CHECKS.
    bk_D = ss["ss_bank_D"];  bk_F = ss["ss_bank_F"]
    fm_D = ss["ss_firm_D"];  fm_F = ss["ss_firm_F"]

    print()
    _rule(65)
    print(_row("Variable", "D (Greece)", "F (Germany)", "Note"))
    _rule(65)

    print(_row("Y_ss",   f"{fm_D['Y_ss']:.4f}",   f"{fm_F['Y_ss']:.4f}", "normalised (Z rescaled)"))
    print(_row("K_ss",   f"{ss['Kap_D_ss']:.4f}", f"{ss['Kap_F_ss']:.4f}"))
    print(_row("I_ss",   f"{fm_D['I_ss']:.4f}",   f"{fm_F['I_ss']:.4f}", "= delta*K"))
    print(_row("C_ss",   f"{ss['C_D_ss']:.4f}",   f"{ss['C_F_ss']:.4f}", "HA aggregate"))
    print(_row("A_ss (HH deposits)", f"{ss['A_D_ss']:.4f}", f"{ss['A_F_ss']:.4f}"))
    print(_row("w_ss",   f"{fm_D['w_ss']:.4f}",   f"{fm_F['w_ss']:.4f}"))
    print(_row("rk_ss (ann %)", f"{ss['rk_D_ss']*400:.3f}", f"{ss['rk_F_ss']*400:.3f}",
               "target ~ 1.8% ann"))

    _rule(65)
    print(_row("n_ss (net worth)", f"{bk_D['n_ss']:.4f}", f"{bk_F['n_ss']:.4f}"))
    print(_row("theta (leverage)", f"{bk_D['theta_ss']:.4f}", f"{bk_F['theta_ss']:.4f}",
               "target 4-6 (GK11)"))
    print(_row("kappa (K/n)", f"{bk_D['kappa_ss']:.4f}", f"{bk_F['kappa_ss']:.4f}"))
    print(_row("phi_bdom (dom bond/n)", f"{bk_D['phi_bdom_ss']:.4f}", f"{bk_F['phi_bdom_ss']:.4f}"))
    print(_row("phi_bfor (for bond/n)", f"{bk_D['phi_bfor_ss']:.4f}", f"{bk_F['phi_bfor_ss']:.4f}"))
    print(_row("Q*b_dom / n",
               f"{ss['Q_bD_ss'] * ss['b_D_D_ss'] / bk_D['n_ss']:.3f}",
               f"{ss['Q_bF_ss'] * ss['b_F_F_ss'] / bk_F['n_ss']:.3f}",
               "≈ 0.08 of assets (Bocola exp^bg 7.6%)"))
    print(_row("alpha (V/n)", f"{bk_D['alpha_ss']:.4f}", f"{bk_F['alpha_ss']:.4f}",
               "franchise value"))
    print(_row("mu (IC mult)", f"{bk_D['mu_ss']:.6f}", f"{bk_F['mu_ss']:.6f}"))
    print(_row("Dep_supply", f"{bk_D['Dep_supply_ss']:.4f}", f"{bk_F['Dep_supply_ss']:.4f}"))
    print(_row("rb_dom (ann %)", f"{bk_D['rb_dom_ss']*400:.3f}", f"{bk_F['rb_dom_ss']*400:.3f}",
               "= rdep + IC spread"))

    _rule(65)
    print(_row("p_ss (RER)", f"{ss['p_ss']:.6f}", "—", "1 = symmetric SS"))
    print(_row("Q_bD_ss / Q_bF_ss", f"{ss['Q_bD_ss']:.5f}", f"{ss['Q_bF_ss']:.5f}",
               "IC-consistent prices"))
    print(_row("b_D_D / b_F_D (D-bank)", f"{ss['b_D_D_ss']:.4f}", f"{ss['b_F_D_ss']:.4f}",
               "dom / for holdings"))
    print(_row("b_F_F / b_D_F (F-bank)", f"{ss['b_F_F_ss']:.4f}", f"{ss['b_D_F_ss']:.4f}"))
    print(_row("F-bank share of D-debt", f"{ss['b_D_F_ss']/cal['B_gov_D_ss']:.1%}", "—",
               "contagion leg, target 20%"))
    print(_row("B_gov / 4Y (debt/GDP)",
               f"{cal['B_gov_D_ss']/(4*fm_D['Y_ss']):.1%}",
               f"{cal['B_gov_F_ss']/(4*fm_F['Y_ss']):.1%}", "≈ 24% (Bocola exposure)"))

    _rule(65)
    print(_row("Calibrated in SS solve", "D", "F", "Pins / target"))
    _rule(65)
    print(_row("lambda (single, IC)", f"{cal['lambda_K_D']:.5f}", f"{cal['lambda_K_F']:.5f}",
               "leverage + credit-spread targets"))
    print(_row("omega_ent (entrants)", f"{cal['omega_ent_D']:.6f}", f"{cal['omega_ent_F']:.6f}",
               "solved jointly with lambda"))
    print(_row("Z_ss (rescaled)", f"{cal['Z_ss_D']:.6f}", f"{cal['Z_ss_F']:.6f}", "pins Y_ss = 1"))
    print(_row("chi (GHH)", f"{cal['chi_D']:.4f}", f"{cal['chi_F']:.4f}", "pins N_ss = 1"))
    print(_row("beta_ss", f"{ss['beta_D_ss']:.6f}", f"{ss['beta_F_ss']:.6f}",
               "deposit-market clearing"))
    print(_row("xr_for (bps ann)", f"{cal['excess_return_F_D_ss']*4e4:.2f}",
               f"{cal['excess_return_D_F_ss']*4e4:.2f}", "foreign-bond FOC anchor"))

    _rule(65)
    print(_row("Residual", "D", "F", "Threshold"))
    _rule(65)
    print(_row("IC  (n_IC/n_ACCUM - 1)",
               f"{bk_D['n_ss_IC'] / bk_D['n_ss_ACCUM'] - 1:.2e}",
               f"{bk_F['n_ss_IC'] / bk_F['n_ss_ACCUM'] - 1:.2e}", "<= 1e-9"))
    print(_row("Deposit (A - Dep)",
               f"{ss['A_D_ss'] - bk_D['Dep_supply_ss']:.2e}",
               f"{ss['A_F_ss'] - bk_F['Dep_supply_ss']:.2e}", "<= 1e-9"))
    print(_row("Walras (Y-C-I-G)",
               f"{fm_D['Y_ss'] - ss['C_D_ss'] - fm_D['I_ss'] - cal['G_D']:.2e}",
               f"{fm_F['Y_ss'] - ss['C_F_ss'] - fm_F['I_ss'] - cal['G_F']:.2e}",
               "F = diagnostic only"))
    _rule(65)
