# CONSOLE REPORTING: STEADY-STATE TABLE. The projection experiments
# (solver_recursive/) print their own IRF tables inline; this keeps output
# formatting out of the model code.
import numpy as np


# UNIT CONVENTION FOR EVERY NUMBER THE EXPERIMENTS PRINT (2026-08-28).
# The model is quarterly. Three different objects used to be reported in three
# different units with no label, and the Bocola comparison was read off the wrong one.
#   RATES (deposit, working-capital, credit spread, sovereign yield) -> ANNUALISED
#     basis points, bp_ann(). 4 x the quarterly rate x 1e4.
#   PROBABILITIES (p^d) -> both the QUARTERLY figure, which is what Bocola's Figure 7
#     plots and what the s-process is calibrated in, and the annualised
#     1 - (1-p)^4, ann_prob().
#   LEVEL RESPONSES (Y, C, I, hours, K, net worth) -> the % deviation from the
#     no-shock path, which is what Bocola's Figures 5 and 7 plot, AND for the FLOW
#     variables an annualised companion ann_pct() = 4 x that.
# ann_pct IS Bocola's Table 5 unit. His output losses are cumsum(g_s - g_ns)*400 where
# g is a quarterly log growth rate, so the cumulated object is the log LEVEL gap and
# the 400 is 100 (to %) x 4 (to an annual rate). His -1.05 / -1.44 / -1.53 are therefore
# level gaps of -0.26 / -0.36 / -0.38%. Reporting both columns is what makes the
# comparison exact in either unit.
# THE LIKE-FOR-LIKE IRF TARGETS, in level %: -0.295 (his closed benchmark) and -0.186
# (his SS V.C open economy, the version whose GHH + working-capital transmission this
# model shares), both at his p^d ~ 2.5-3.0%/qtr shock; -0.222 and -0.157 rescaled to
# a p^d = 1.98% shock.
BOCOLA_IRF_CLOSED = -0.2225      # level %, his closed benchmark at p^d = 1.98%/qtr
BOCOLA_IRF_OPEN = -0.157         # level %, his open economy at the same shock
BOCOLA_EPISODE_LEVEL = -0.36     # level %, 2011Q4 Italian episode (Table 5 / 4)


def bp_ann(rate_q):
    # QUARTERLY RATE -> ANNUALISED BASIS POINTS.
    return 4e4 * rate_q


def ann_pct(dev_pct):
    # QUARTERLY-FLOW LEVEL GAP IN % -> THE SAME GAP AT AN ANNUAL RATE (Bocola x400).
    return 4.0 * dev_pct


def ann_prob(p_q):
    # QUARTERLY PROBABILITY -> ANNUAL.
    return 1.0 - (1.0 - p_q) ** 4


def banner(text, width=65):
    # FULL-WIDTH SECTION HEADER.
    print("\n" + "=" * width)
    print(f"  {text}")
    print("=" * width)


def print_solve_stage(label, ok, its, worst, n_fail, n_pts):
    # ONE LINE PER TIME-ITERATION STAGE, INCLUDING FAILURE. A stage that does not
    # converge must SAY SO: the exit test needs both a settled rule and every point
    # clearing, so a run can look finished while part of the grid is frozen on its
    # cold start. Silence here is what let that go unnoticed.
    # The exit test is BOTH a settled rule AND every point clearing, so a bare "did not
    # converge" conflates two very different states: a stage sitting at 1e-14 with zero
    # frozen points that merely ran out of sweep budget, and a stage with points that
    # never solved. Only the second is a failure to act on -- say which.
    if ok:
        tag, note = "converged", ""
    elif n_fail == 0 and worst < 1e-6:
        tag = "residuals OK"
        note = "  (rule-change tol not reached in the sweep budget; no point unsolved)"
    else:
        tag = "DID NOT CONVERGE"
        note = f"  <-- {n_fail}/{n_pts} points frozen (unsolved)"
    print(f"    {label:<34s} {tag:>16s} in {its:3d} sweeps   "
          f"max|F|={worst:.2e}{note}")


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
    print(_row("C_ss",   f"{ss['C_D_ss']:.4f}",   f"{ss['C_F_ss']:.4f}", "HA aggregate, basket units"))
    print(_row("A_ss (HH deposits)", f"{ss['A_D_ss']:.4f}", f"{ss['A_F_ss']:.4f}", "basket units"))
    print(_row("w_ss",   f"{fm_D['w_ss']:.4f}",   f"{fm_F['w_ss']:.4f}"))
    print(_row("rk_ss (ann %)", f"{ss['rk_D_ss']*400:.3f}", f"{ss['rk_F_ss']*400:.3f}",
               "= r_dep + spread instrument"))
    print(_row("country mass", f"{cal['size_D']:.1f}", f"{cal['size_F']:.1f}",
               "F/D = 8: D is a small member"))
    print(_row("omega_home", f"{cal['omega_home_D']:.5f}", f"{cal['omega_home_F']:.5f}",
               "size-consistent: (1-w_F)=(1-w_D)/8"))

    _rule(65)
    print(_row("n_ss (net worth)", f"{bk_D['n_ss']:.4f}", f"{bk_F['n_ss']:.4f}"))
    print(_row("theta (leverage)", f"{bk_D['theta_ss']:.4f}", f"{bk_F['theta_ss']:.4f}",
               "instrument; 8 at the rest point"))
    print(_row("kappa (K/n)", f"{bk_D['kappa_ss']:.4f}", f"{bk_F['kappa_ss']:.4f}"))
    print(_row("phi_bdom (dom bond/n)", f"{bk_D['phi_bdom_ss']:.4f}", f"{bk_F['phi_bdom_ss']:.4f}"))
    print(_row("phi_bfor (for bond/n)", f"{bk_D['phi_bfor_ss']:.4f}", f"{bk_F['phi_bfor_ss']:.4f}"))
    own_D = ss["Q_bD_ss"] * ss["b_D_D_ss"]
    own_F = ss["Q_bF_ss"] * ss["b_F_F_ss"]
    print(_row("own sovereign / assets", f"{100 * own_D / (bk_D['theta_ss'] * bk_D['n_ss']):.3f}%",
               f"{100 * own_F / (bk_F['theta_ss'] * bk_F['n_ss']):.3f}%",
               "instrument; EBA 14.4% / 6.5% at rest"))
    print(_row("own sovereign / n", f"{own_D / bk_D['n_ss']:.3f}", f"{own_F / bk_F['n_ss']:.3f}",
               "GR sovereign / CT1: 2.39"))
    print(_row("alpha (V/n)", f"{bk_D['alpha_ss']:.4f}", f"{bk_F['alpha_ss']:.4f}",
               "franchise value"))
    print(_row("mu (IC mult)", f"{bk_D['mu_ss']:.6f}", f"{bk_F['mu_ss']:.6f}"))
    print(_row("Dep_supply (goods)", f"{bk_D['Dep_supply_ss']:.4f}", f"{bk_F['Dep_supply_ss']:.4f}"))
    print(_row("rb_dom (ann %)", f"{bk_D['rb_dom_ss']*400:.3f}", f"{bk_F['rb_dom_ss']*400:.3f}",
               "= rdep + IC spread"))

    _rule(65)
    print(_row("p_ss (RER)", f"{ss['p_ss']:.6f}", "—", "D-goods per F-good; free"))
    print(_row("P_CES_ss", f"{ss['P_CES_D_ss']:.6f}", f"{ss['P_CES_F_ss']:.6f}",
               "consumption-basket price"))
    print(_row("NX / Y (per capita)", f"{100 * ss['NX_D_ss'] / fm_D['Y_ss']:+.4f}%",
               f"{100 * ss['NX_F_ss'] / fm_F['Y_ss']:+.4f}%", "pays the net coupons"))
    print(_row("net coupons from abroad / Y", f"{100 * ss['NFI_D_ss'] / fm_D['Y_ss']:+.4f}%",
               f"{100 * ss['NFI_F_ss'] / fm_F['Y_ss']:+.4f}%", "NX + NFI = 0"))
    print(_row("Q_bD_ss / Q_bF_ss", f"{ss['Q_bD_ss']:.5f}", f"{ss['Q_bF_ss']:.5f}",
               "IC-consistent prices"))
    # SOVEREIGN HOLDINGS ARE IN THE ISSUER'S PER-CAPITA UNITS, so b_D_F_ss is the slice
    # of D's own stock held abroad and the F BANK's own book carries b_D_F_ss/sz of it.
    # Printing the raw state next to the F bank's balance sheet would overstate the
    # F bank's exposure by the mass ratio.
    sz = cal["size_F"] / cal["size_D"]
    print(_row("b_D_D / b_F_D (D-bank)", f"{ss['b_D_D_ss']:.4f}", f"{ss['b_F_D_ss']:.4f}",
               "dom / for holdings, per D capita"))
    print(_row("b_F_F / b_D_F (F-bank)", f"{ss['b_F_F_ss']:.4f}",
               f"{ss['b_D_F_ss']/sz:.4f}", "per F capita (= b_D_F/sz)"))
    print(_row("held by the other's banks", f"{ss['b_D_F_ss']/cal['B_gov_D_ss']:.1%}",
               f"{ss['b_F_D_ss']/sz/cal['B_gov_F_ss']:.2%}", "of each stock (EBA: 20%, ~0)"))
    print(_row("B_gov / 4Y (face)",
               f"{cal['B_gov_D_ss']/(4*fm_D['Y_ss']):.1%}",
               f"{cal['B_gov_F_ss']/(4*fm_F['Y_ss']):.1%}", "bank-held debt"))
    print(_row("Q*B_gov / 4Y (market)",
               f"{ss['Q_bD_ss']*cal['B_gov_D_ss']/(4*fm_D['Y_ss']):.1%}",
               f"{ss['Q_bF_ss']*cal['B_gov_F_ss']/(4*fm_F['Y_ss']):.1%}", "EBA-held GR debt/GDP: 40%"))

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
               f"{cal['excess_return_D_F_ss']*4e4:.2f}", "foreign bond over own spread"))
    print(_row("friction ref - book", f"{cal['b_F_D_ref'] - cal['b_F_D_ss']:+.2e}",
               f"{cal['b_D_F_ref'] - cal['b_D_F_ss']:+.2e}", "makes the SS books optimal"))
    print(_row("debt root (Bohn)",
               f"{1 - cal['delta_b_D'] + (cal['delta_b_D'] - ss['gs_D']['gamma_tau']) / ss['gs_D']['Q_B_ss']:.6f}",
               f"{1 - cal['delta_b_F'] + (cal['delta_b_F'] - ss['gs_F']['gamma_tau']) / ss['gs_F']['Q_B_ss']:.6f}",
               f"gamma_tau {ss['gs_D']['gamma_tau']:.5f}"))

    _rule(65)
    print(_row("Residual", "D", "F", "Threshold"))
    _rule(65)
    print(_row("IC  (n_IC/n_ACCUM - 1)",
               f"{bk_D['n_ss_IC'] / bk_D['n_ss_ACCUM'] - 1:.2e}",
               f"{bk_F['n_ss_IC'] / bk_F['n_ss_ACCUM'] - 1:.2e}", "<= 1e-10"))
    print(_row("Deposit (A*P_CES - Dep)",
               f"{ss['A_D_ss'] * ss['P_CES_D_ss'] - bk_D['Dep_supply_ss']:.2e}",
               f"{ss['A_F_ss'] * ss['P_CES_F_ss'] - bk_F['Dep_supply_ss']:.2e}", "<= 1e-9"))
    print(_row("Walras (Y-PC-I-G-NX)", f"{ss['goods_gap_D_ss']:.2e}", f"{ss['goods_gap_F_ss']:.2e}",
               "<= 1e-8, F = diagnostic"))
    print(_row("current account", f"{ss['NX_D_ss'] + ss['NFI_D_ss']:.2e}",
               f"{ss['NX_F_ss'] + ss['NFI_F_ss']:.2e}", "<= 1e-8"))
    _rule(65)


def _fmt_value(v, unit):
    # ONE MODEL NUMBER IN ITS TABLE UNIT.
    if unit == "x":
        return f"{v:.3f}"
    if unit == "bp":
        return f"{v:.1f}bp"
    if unit == "years":
        return f"{v:.2f} yrs"
    return f"{v:.3g}%" if abs(v) < 1.0 else f"{v:.2f}%"


def print_calibration_table(vals, moments, cal):
    # THE CALIBRATION TABLE AT THE STOCHASTIC REST POINT, AND THE SIX TARGETS AGAINST TOLERANCE.
    from global_projection.solver_recursive.calibrate_stochastic import PAIRS, TABLE_ROWS
    print("\n  CALIBRATION TABLE (model at the stochastic rest point; rates annualised, "
          "p^d quarterly)")
    print(f"   {'moment':<52s}{'data':>15s}{'model':>12s}   role")
    block = None
    for blk, label, data, key, unit, role in TABLE_ROWS:
        if blk != block:
            print(f"   {blk}")
            block = blk
        print(f"     {label:<50s}{data:>15s}{_fmt_value(vals[key], unit):>12s}   {role}")
    print("   targets at the rest point:")
    for inst, tgt, mom, tol, *_ in PAIRS:
        gap = moments[mom] - cal[tgt]
        print(f"     {mom:<9s} {moments[mom]:.8g} vs {cal[tgt]:.8g}  gap {gap:+.2e} "
              f"(tol {tol:.1e})  {'OK' if abs(gap) <= tol else 'MISS'}   instrument "
              f"{inst} = {cal[inst]:.8g}")
    print(f"   the solution box stops at p^d = {vals['pd_box_top']:.2f}%/qtr, so the 9.1% "
          f"exceedance is the AR(1)'s own ergodic share (z = {vals['crisis_z']:.2f}); "
          f"p at rest {vals['p_rest']:.5f}, mu_D {vals['mu_D']:.5f}, mu_F {vals['mu_F']:.5f}, "
          f"read residual {vals['read_resid']:.0e}")


def print_sovereign_spread(legs, label=""):
    # THE SOVEREIGN SPREAD, LEG BY LEG -- WHAT ANY GIVEN INSTRUMENT COULD EVER REACH.
    # Each number is the yield the bond would LOSE if that leg were removed entirely, so
    # it is the CEILING on an instrument that acts only through that leg. The reason this
    # table exists: a bank-liquidity facility can touch the liquidity leg and nothing
    # else, and on this calibration that leg is ~2% of the D-F spread, which is why the
    # retired LTRO moved the credit spread by tens of basis points and the sovereign spread
    # by single digits. Legs are removals, not a partition -- y is convex in q.
    from global_projection.solver_recursive.output_decomposition import SOVEREIGN_LEGS
    print(f"\n  SOVEREIGN SPREAD DECOMPOSITION{(' - ' + label) if label else ''}"
          f"   (annualised bp)")
    print(f"   {'leg':<32s}{'y_D':>10s}{'y_F':>10s}{'SPREAD':>10s}{'% of spread':>13s}")
    for k, nm in SOVEREIGN_LEGS:
        sh = 100 * legs[f"spread_{k}"] / legs["spread"] if legs["spread"] else float("nan")
        print(f"   {nm:<32s}{legs[f'y_D_{k}']:10.1f}{legs[f'y_F_{k}']:10.1f}"
              f"{legs[f'spread_{k}']:10.1f}{sh:12.1f}%")
    print(f"   {'ACTUAL YIELD':<32s}{legs['y_D']:10.1f}{legs['y_F']:10.1f}"
          f"{legs['spread']:10.1f}")
    print(f"   FOC closure off-node: D {100*legs['foc_closure_D']:+.3f}%, "
          f"F {100*legs['foc_closure_F']:+.3f}%  (the solve is exact only AT the nodes)")


def print_output_decomposition(dec, channels):
    # IMPACT OF EACH ACTIVE OUTPUT CHANNEL AND THE TOTAL, LEVEL % AND ANNUALISED.
    for k, lab in channels:
        print(f"    {lab:<22s} impact {dec[k][0]:+8.4f}%  "
              f"({4 * dec[k][0]:+8.4f}% annualised)")
    print(f"    {'TOTAL':<22s} impact {dec['total'][0]:+8.4f}%  "
          f"({4 * dec['total'][0]:+8.4f}% annualised)")


def print_bond_decomposition(bdec, channels):
    # BOCOLA'S SPLIT OF THE D BOND'S EXPECTED EXCESS RETURN ON IMPACT, bp/yr FROM THE NO-SHOCK PATH.
    print("\n  D-bond expected excess return E[R_b] - R (Bocola 2016 split; bp/yr, deviation from the no-shock path)")
    for k, lab in channels:
        print(f"    {lab:<40s} impact {bdec[k][0]:+9.2f} bp")
    print(f"    {'TOTAL':<40s} impact {bdec['total'][0]:+9.2f} bp")


def print_tpi_report(res, bond_channels):
    # THE TPI AGAINST NO TPI: REST POINT, IMPACT, FOOTPRINT, PRICE IMPACT, KT CONDITIONS.
    # modes "relief" and "cap" (tpi_experiment.intervention): unanticipated purchases, so no rest-point gap and no rule
    # to check
    relief = res.get("mode") in ("relief", "cap")
    if relief:
        r = res["relief"]
        if res.get("mode") == "cap":
            m_path = np.asarray(res["irf_on"]["m_cb"], dtype=float)
            later = ", ".join(f"q{t} {m_path[t]:.1f}%" for t in range(1, len(m_path)) if m_path[t] > 1e-6) or "none"
            print(f"\n  TPI: spread cap {res['cap_bp']:.0f} bp/yr over the F bond, held by unanticipated purchases held to "
                  f"maturity: {100 * res['purchase']:.1f}% of the SS D debt at impact, later {later}; capital key D "
                  f"{res['key_D']:.3f}")
        else:
            print(f"\n  TPI: relief intervention at the shock, unanticipated: {100 * res['purchase']:.1f}% of the SS D "
                  f"debt bought and held to maturity, capital key D {res['key_D']:.3f}")
        print(f"   on impact: mu_D {r['mu_off']:.4f} -> {r['mu_on']:.4f}, lending spread {r['lend_off']:.0f} -> "
              f"{r['lend_on']:.0f} bp/yr, sovereign spread {r['sov_off']:.0f} -> {r['sov_on']:.0f} bp/yr")
    else:
        cap = res["cap_bp"]
        solved = res["cap_solved"]
        print(f"\n  TPI: spread cap {cap:.0f} bp/yr over the F bond, capital key D "
              f"{res['key_D']:.3f}"
              + ("" if res["solve_ok"] else f"   [NOT ROOTED at the target: last cap "
                                             f"solved {solved}]"))
        off, on = res["rest_off"], res["rest_on"]
        print("\n  REST POINT (nothing is bought here: the gap is the ANNOUNCEMENT)")
        print(f"   {'object':<26s}{'no TPI':>12s}{'TPI':>12s}{'diff':>12s}")
        for k, lab in (("Q_bD", "Q_bD"), ("sov_bp", "sovereign spread bp/yr"),
                       ("spread_bp", "credit spread bp/yr"), ("mu", "mu_D"),
                       ("E_Om", "E[Omega_D]"), ("alpha", "alpha_D (franchise)"),
                       ("n", "n_D (bank net worth)"), ("Y", "Y_D"), ("C", "C_D"),
                       ("I", "I_D"), ("M", "Eurosystem book, % of B"),
                       ("m", "purchases, % of B")):
            print(f"   {lab:<26s}{off[k]:12.5f}{on[k]:12.5f}{on[k] - off[k]:+12.5f}")
    io, im = res["irf_off"], res["irf_on"]
    print("\n  HEADLINE RISK SHOCK, impact (q0) and trough        no TPI        TPI")
    for k, lab in (("Y", "GDP % (level)"), ("sov_bp", "sovereign spread bp/yr"),
                   ("spread", "credit spread bp/yr"), ("n", "bank net worth %"),
                   ("dQ_bD", "Q_bD % vs no-shock")):
        print(f"   {lab:<26s} q0  {io[k][0]:+10.4f} {im[k][0]:+10.4f}"
              f"    min {np.min(io[k]):+10.4f} {np.min(im[k]):+10.4f}")
    print(f"   footprint: peak purchases {np.max(im['m_cb']):.2f}% of B in one quarter, "
          f"peak book {np.max(im['M_cb']):.2f}% of B, "
          f"worst P&L {np.min(im['Pi_cb']):+.4f}% of quarterly GDP")
    print("\n  THE TPI's EFFECT ON THE D-BOND EXCESS RETURN, BOCOLA SPLIT (bp/yr, TPI minus no TPI)")
    for k, lab in bond_channels:
        print(f"    {lab:<40s} q0 {res['bond_legs'][k][0]:+9.2f} bp")
    print(f"    {'TOTAL':<40s} q0 {res['bond_legs']['total'][0]:+9.2f} bp")
    for tag, rows in (("at the rest point", res["impact_rest"]),
                      ("at the headline shock", res["impact_shock"])):
        q0 = rows[0]["Q_bD"]
        print(f"\n  PRICE IMPACT OF A ONE-OFF PURCHASE {tag} (rules fixed)")
        for r in rows:
            print(f"    m = {100 * r['m']:5.1f}% of B: Q_bD {100 * (r['Q_bD'] / q0 - 1):+7.3f}%"
                  f"  spread {r['sov_bp']:7.1f} bp  mu_D {r['mu']:.4f}  b_DD {r['b_DD']:.3f}"
                  f"  (resid {r['resid']:.0e})")
    if relief:
        return
    kt = res["kt"]
    print(f"\n  COMPLEMENTARITY on the solved rules ({kt['n']} no-default nodes): floor binds "
          f"at {kt['n_binding']}, corner b_DD = 0 at {kt['n_corner']}")
    print(f"    m >= 0: min {kt['min_m']:+.1e}   Q >= floor: min gap {kt['min_gap']:+.1e}"
          f"   max|m*gap| {kt['max_m_gap']:.1e}")
    print(f"    b_DD >= 0: min {kt['min_b']:+.1e}   FOC gap <= 0: max {kt['max_foc']:+.1e}"
          f"   max|b_DD*gap| {kt['max_b_foc']:.1e}")
