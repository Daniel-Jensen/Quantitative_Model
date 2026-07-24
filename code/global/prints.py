# ALL CONSOLE REPORTING FOR THE PIPELINE: SS TABLE, RESIDUAL CHECKS, EXPERIMENT DIAGNOSTICS.
# main.py orchestrates the solves; every table and headline number is printed
# from here, so output formatting never mixes with the model code.
import numpy as np

from transition import market_residuals
from risk_branch import bond_decomposition


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


def lending_spread_bps(out, ss, cal, country="D"):
    # LENDING SPREAD rk_t - rdep_{t-1} IN ANNUALIZED BPS DEVIATION FROM THE SS.
    # The deposit rate paid at t was locked at t-1, hence the lag; the SS
    # reference is the SS spread rk_ss - rdep_ss, so an unshocked path reads 0.
    rdep_ss = cal[f"r_dep_{country}_target"]
    rdep = np.asarray(out[f"rdep_{country}"])
    rdep_lag = np.concatenate([[rdep_ss], rdep[:-1]])
    return 4e4 * ((np.asarray(out[f"rk_{country}"]) - rdep_lag)
                  - (ss[f"rk_{country}_ss"] - rdep_ss))


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
               "Bocola GIPS target ~ 0.93"))
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
               f"{cal['B_gov_F_ss']/(4*fm_F['Y_ss']):.1%}", "target ~ 93%"))

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


def print_transition_residuals(out, cal, ss=None):
    # MARKET-CLEARING RESIDUALS OF A SOLVED PATH (goods_F IS THE WALRAS ACCURACY CHECK).
    # Pass ss whenever TPI is active, so the CB rebate's cross-border transfer
    # enters both goods residuals exactly as it does in the imposed system.
    r = market_residuals(out, cal, ss=ss)
    print(f"  max|mu*slack D| (complementarity) = {r['cap_D']:.2e}"
          f"  (min mu {r['mu_min_D']:+.2e}, min slack {r['slack_min_D']:+.2e})")
    print(f"  max|mu*slack F| (complementarity) = {r['cap_F']:.2e}"
          f"  (min mu {r['mu_min_F']:+.2e}, min slack {r['slack_min_F']:+.2e})")
    print(f"  max|union deposit clearing|       = {r['dep_union']:.2e}")
    print(f"  max|deposit UIP|                  = {r['uip']:.2e}")
    print(f"  max|goods mkt D|                  = {r['goods_D']:.2e}")
    print(f"  max|goods mkt F| [diagnostic]     = {r['goods_F']:.2e}")
    print(f"  max|NFA deposit position|         = {r['nfa_dep_D']:.4f}"
          "  (cross-border margin, 0 pre-union)")


def print_risk_table(out, ss, cal, pi_shock, pi_rho):
    # BOCOLA PASS-THROUGH DIAGNOSTICS FOR THE SOVEREIGN-RISK RUN.
    # Works for both pricing modes: liquidity-only (plain solve_transition,
    # no branch) and the two-branch risk channel (solve_transition_risk).
    bk_D = ss["ss_bank_D"];  bk_F = ss["ss_bank_F"];  fm_D = ss["ss_firm_D"]
    dec = bond_decomposition(out, ss, cal)
    lend = lending_spread_bps(out, ss, cal)
    i_peak = int(np.argmax(dec["total_yield"]))
    br = out.get("branch")
    mode = ("liquidity channel (risk-neutral pricing)" if br is None
            else "two-branch risk channel")

    def line(label, val, note=""):
        # ONE SINGLE-VALUE DIAGNOSTIC LINE.
        return _row(label, val, "", note, label_w=30)

    print()
    _rule(72)
    print(f"  Bocola pass-through — {mode}  (pi_0={pi_shock:.0%}, rho={pi_rho})")
    print(_row("Statistic", "Value", "", "Note", label_w=30))
    _rule(72)
    print(line("Q_bD[0]  (% dev from SS)",
               f"{(out['Q_bD'][0]/ss['Q_bD_ss']-1)*100:+.2f}%", "MTM repricing"))
    print(line("n_D[0]  (% dev)",
               f"{(out['n_D'][0]/bk_D['n_ss']-1)*100:+.2f}%", "no default on base path"))
    print(line("n_F[0]  (% dev)",
               f"{(out['n_F'][0]/bk_F['n_ss']-1)*100:+.2f}%", "contagion leg"))
    print(line("Y_D trough (% dev)", f"{np.min(out['Y_D']/fm_D['Y_ss']-1)*100:+.3f}%"))
    print(line("I_D[0]  (% dev)", f"{(out['I_D'][0]/fm_D['I_ss']-1)*100:+.3f}%"))
    print(line("Sov spread peak (bps ann)", f"{dec['total_yield'][i_peak]:+.0f}",
               f"t={i_peak}: def {dec['defcomp'][i_peak]:+.0f} + "
               f"risk {dec['risk'][i_peak]:+.0f} + liq {dec['liquidity'][i_peak]:+.0f}"))
    print(line("Lending spread peak (bps)", f"{np.max(lend):+.0f}"))
    print(line("b_gov_D peak", f"{np.max(out['b_gov_D']):.3f}",
               f"SS={cal['B_gov_D_ss']:.3f}; Tax peak dev "
               f"{np.max(out['Tax_D'])-ss['Tax_D_ss']:+.4f}"))
    print(line("def_real = 0 on base path", str(bool(np.all(out['def_real_D'] == 0)))))
    print(line("min mu (IC mult) D / F", f"{out['mu_min_D']:+.4f}",
               f"F {out['mu_min_F']:+.4f}; 0 = IC slack somewhere (occ. binding)"))
    print(line("min IC slack D / F", f"{out['slack_min_D']:+.2e}",
               f"F {out['slack_min_F']:+.2e}; negative = complementarity violated"))
    print(line("rdep_D - rdep_F [0] (bps ann)",
               f"{(out['rdep_D'][0] - out['rdep_F'][0]) * 4e4:+.1f}",
               "= p-growth (real-rate parity in the union)"))
    print(line("NFA deposit position peak",
               f"{out['nfa_dep_D'][int(np.argmax(np.abs(out['nfa_dep_D'])))]:+.4f}",
               "D claim on the union interbank (D-goods)"))
    if br is not None:
        print(line("[branch] n_D(0)/n_ss", f"{br['n_D'][0]/bk_D['n_ss']:.3f}",
                   f"Y_D(0) dev {(br['Y_D'][0]/fm_D['Y_ss']-1)*100:+.2f}%  (feared default state)"))
    _rule(72)


def print_tpi_table(out, out_base, ss, cal, pi_shock, s_on_date, doubt_peak):
    # TPI DIAGNOSTICS: INTERVENTION SIZE, CB CASH FLOW, AND COMPRESSION VS THE NO-TPI RUN.
    bk_D = ss["ss_bank_D"];  fm_D = ss["ss_firm_D"]
    dec_on   = bond_decomposition(out, ss, cal)
    dec_base = bond_decomposition(out_base, ss, cal)
    lend_on   = lending_spread_bps(out, ss, cal)
    lend_base = lending_spread_bps(out_base, ss, cal)

    def line(label, val, note=""):
        # ONE SINGLE-VALUE DIAGNOSTIC LINE.
        return _row(label, val, "", note, label_w=30)

    print()
    _rule(76)
    print(f"  TPI backstop  (same pi_0={pi_shock:.0%} shock; backstop on from "
          f"t={s_on_date}, doubt peak {doubt_peak:.0%})")
    print(_row("Statistic", "TPI-on", "", "Note (vs no-TPI)", label_w=30))
    _rule(76)
    print(line("cb_buy_D peak", f"{np.max(out['cb_buy_D']):.4f}",
               f"t={int(np.argmax(out['cb_buy_D']))}; "
               f"{np.max(out['cb_buy_D'])/cal['B_gov_D_ss']*100:.2f}% of B_gov_D_ss"))
    print(line("rem_cb_D peak / trough", f"{np.max(out['rem_cb_D']):+.4f}",
               f"trough {np.min(out['rem_cb_D']):+.4f}  "
               "(+ = CB net income, - = net outlay)"))
    print(line("Q_bD[0]  (% dev from SS)",
               f"{(out['Q_bD'][0]/ss['Q_bD_ss']-1)*100:+.2f}%",
               f"no-TPI {(out_base['Q_bD'][0]/ss['Q_bD_ss']-1)*100:+.2f}%"))
    print(line("Sov spread peak (bps ann)", f"{np.max(dec_on['total_yield']):+.0f}",
               f"no-TPI {np.max(dec_base['total_yield']):+.0f}"))
    print(line("Lending spread peak (bps)", f"{np.max(lend_on):+.0f}",
               f"no-TPI {np.max(lend_base):+.0f}"))
    print(line("Y_D trough (% dev)", f"{np.min(out['Y_D']/fm_D['Y_ss']-1)*100:+.3f}%",
               f"no-TPI {np.min(out_base['Y_D']/fm_D['Y_ss']-1)*100:+.3f}%"))
    print(line("n_D[0]  (% dev)", f"{(out['n_D'][0]/bk_D['n_ss']-1)*100:+.2f}%",
               f"no-TPI {(out_base['n_D'][0]/bk_D['n_ss']-1)*100:+.2f}%"))
    print(line("min mu (IC mult) D / F", f"{out['mu_min_D']:+.4f}",
               f"F {out['mu_min_F']:+.4f}"))
    print(line("min IC slack D / F", f"{out['slack_min_D']:+.2e}",
               f"F {out['slack_min_F']:+.2e}; negative = complementarity violated"))
    _rule(76)
