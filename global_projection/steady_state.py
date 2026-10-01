# TWO-COUNTRY STEADY-STATE SOLVER: STAGE 1 {rk_D, rk_F, p}, STAGE 2 {beta_D, beta_F}.
# The steady state must be SYMMETRIC: country asymmetries enter through shocks
# only. An asymmetric SS shifts p_ss off 1 and opens an O(1e-4) goods-market
# wedge, because p is only weakly identified by external balance at a trade
# elasticity of 0.5.
#
# THE ANATOMY (solve_steady_state):
#   _calibrate_bank_frictions   lambda, omega_ent from the leverage and spread targets
#   _sovereign_holdings         SS bond books in issuer per-capita units
#   _solve_stage1               {rk_D, rk_F, p}: bank n_IC = n_ACCUM + external balance
#   _rescale_tfp                Z_ss so that Y_ss = 1
#   _country_blocks             capital -> firm -> working-capital loan -> bank, per country
#   _solve_stage2               {beta_D, beta_F}: household deposits = bank deposit supply
#   _goods_market_check         Walras check with the stage-2 consumption
import numpy as np
from scipy.optimize import brentq, root

from global_projection.blocks.rouwenhorst import rouwenhorst
from global_projection.blocks.household import make_asset_grid, solve_steady_state_household
from global_projection.blocks.distribution import stationary_distribution, aggregate_assets, aggregate_consumption
from global_projection.blocks.firms import steady_state_firm, markup_ss
from global_projection.blocks.capital import capital_demand
from global_projection.blocks.bank import steady_state_bank, calibrate_bank_targets
from global_projection.blocks.government import govt_steady_state
from global_projection.blocks.trade import ces_price, import_demand, trade_balance, size_ratio


def _calibrate_bank_frictions(cal, verbose):
    # SINGLE lambda AND omega_ent PER COUNTRY FROM LEVERAGE + SPREAD, WRITTEN INTO cal.
    for c in ("D", "F"):
        lam, om, *_ = calibrate_bank_targets(
            cal[f"beta_inter_{c}"], cal[f"f_{c}"], cal[f"r_dep_{c}_target"],
            cal[f"leverage_target_{c}"], cal[f"credit_spread_target_{c}"],
        )
        cal[f"lambda_K_{c}"] = cal[f"lambda_bD_{c}"] = cal[f"lambda_bF_{c}"] = lam
        cal[f"omega_ent_{c}"] = om
        cal[f"rk_{c}_guess"]  = cal[f"r_dep_{c}_target"] + cal[f"credit_spread_target_{c}"]
        if verbose:
            print(f"[bank-cal {c}] lambda={lam:.6f}  omega_ent={om:.6f}  "
                  f"(target theta={cal[f'leverage_target_{c}']:.2f}, "
                  f"spread={cal[f'credit_spread_target_{c}']*4e4:.0f} bps/yr)")


def _sovereign_holdings(cal):
    # SS BOND BOOKS, CARRIED IN THE ISSUER'S PER-CAPITA UNITS.
    # b_D_F_ss is the slice of D's own per-capita stock held abroad, so the F bank's
    # per-capita book is b_D_F_ss/sz (sz = size_F/size_D). At sz = 1 this is the
    # symmetric block exactly; at sz = 8 both banks still end up with IDENTICAL
    # per-capita balance sheets.
    sz = size_ratio(cal)
    b_F_D_ss = cal["b_F_D_ss"]           # D's per-capita holding of the F sovereign
    b_D_F_ss = cal["b_D_F_ss"]           # F's holding of the D sovereign, in D units
    return dict(b_F_D_ss=b_F_D_ss, b_D_F_ss=b_D_F_ss,
                b_D_D_ss=cal["B_gov_D_ss"] - b_D_F_ss,
                b_F_F_ss=cal["B_gov_F_ss"] - b_F_D_ss / sz,
                b_D_F_ss_pc=b_D_F_ss / sz)   # the same holding per F capita


def _country_blocks(cal, rk, p, mc, Q_bD, Q_bF, b_own, b_for, c):
    # CAPITAL -> FIRM -> WORKING-CAPITAL LOAN -> BANK FOR ONE COUNTRY AT (rk, p).
    # The firm block comes first: the loan the bank holds is zeta * the wage bill, and
    # w_ss is a firm object (N_ss = 1).
    Kap = capital_demand(rk, mc, cal, country=c)
    fm = steady_state_firm(cal, Kap, country=c)
    L_wc = cal[f"zeta_wc_{c}"] * fm["w_ss"]
    bk = steady_state_bank(cal, rk, Kap, Q_bD, Q_bF, b_own, b_for, p, country=c,
                           L_wc_ss=L_wc)
    return Kap, fm, bk


def _stage1_residual(x, cal, mc, Q_b, hold):
    # CAPITAL-MARKET (n_IC = n_ACCUM) + EXTERNAL-BALANCE RESIDUALS IN {rk_D, rk_F, p}.
    rk_D, rk_F, p = x
    _, fm_D, bk_D = _country_blocks(cal, rk_D, p, mc["D"], Q_b["D"], Q_b["F"],
                                    hold["b_D_D_ss"], hold["b_F_D_ss"], "D")
    _, fm_F, bk_F = _country_blocks(cal, rk_F, p, mc["F"], Q_b["D"], Q_b["F"],
                                    hold["b_F_F_ss"], hold["b_D_F_ss_pc"], "F")
    res_cap_D = (bk_D["n_ss_IC"] - bk_D["n_ss_ACCUM"]) / bk_D["n_ss_ACCUM"]
    res_cap_F = (bk_F["n_ss_IC"] - bk_F["n_ss_ACCUM"]) / bk_F["n_ss_ACCUM"]
    P_CES_D = ces_price(p, cal, country="D")
    P_CES_F = ces_price(p, cal, country="F")
    IM_D = import_demand(p, fm_D["C_ss"], P_CES_D, cal, country="D")
    IM_F = import_demand(p, fm_F["C_ss"], P_CES_F, cal, country="F")
    NX_D, _ = trade_balance(p, IM_D, IM_F, cal)
    # cross-border coupon income, already priced into Q; both legs are in D-per-capita
    # units already, so no mass ratio
    rb_D_mkt = cal["r_dep_D_target"] + bk_D["IC_spread_dom"]
    rb_F_mkt = cal["r_dep_F_target"] + bk_F["IC_spread_dom"]
    income_in_D  = rb_F_mkt * p * bk_F["Q_bdom_IC"] * hold["b_F_D_ss"]
    income_out_D = rb_D_mkt * bk_D["Q_bdom_IC"] * hold["b_D_F_ss"]
    res_ext = (NX_D + income_in_D - income_out_D) / fm_D["Y_ss"]
    return res_cap_D, res_cap_F, res_ext


def _solve_stage1(cal, mc, Q_b, hold, verbose):
    # ROOT-FIND {rk_D, rk_F, p}; FAIL LOUDLY IF THE BANK IDENTITY DOES NOT CLOSE.
    ncalls = [0]

    def resid(x):
        # ONE STAGE-1 EVALUATION, WITH THE FAILURE SENTINEL AND THE PROGRESS LINE.
        ncalls[0] += 1
        try:
            res_cap_D, res_cap_F, res_ext = _stage1_residual(x, cal, mc, Q_b, hold)
        except (RuntimeError, ValueError, FloatingPointError, ZeroDivisionError):
            return [1e3, 1e3, 1e3]
        if verbose:
            print(f"  stage1 call {ncalls[0]:3d}: rk_D={x[0]:.5f}  rk_F={x[1]:.5f}  "
                  f"p={x[2]:.4f}  |resid|=[{res_cap_D:.3e},{res_cap_F:.3e},{res_ext:.3e}]")
        return [res_cap_D, res_cap_F, res_ext]

    if verbose:
        print("=== Stage 1: capital markets + external balance {rk_D, rk_F, p} ===")
    sol1 = root(resid, [cal["rk_D_guess"], cal["rk_F_guess"], 1.0],
                method="hybr", options={"xtol": cal["tol_mkt"], "maxfev": 2000})
    if not sol1.success and verbose:
        print(f"  Warning: stage1 hybr did not flag success "
              f"(resid={np.max(np.abs(sol1.fun)):.2e})")
    # the capital-market residuals ARE the bank n_IC/n_ACCUM identity: a large value
    # means the (f, spread, leverage) targets sit on the franchise fold's UPPER root
    assert np.max(np.abs(sol1.fun[:2])) < 1e-6, (
        f"stage-1 n_IC/n_ACCUM inconsistent (max|res_cap|={np.max(np.abs(sol1.fun[:2])):.2e}): "
        f"the (f, spread, leverage) targets are fold-blocked -- raise f until leverage is the "
        f"least root (f>=~0.14 at spread 720bp).")
    return sol1.x


def _rescale_tfp(cal, rk, mc):
    # SET Z_ss SO THAT Y_ss = 1 (EXACT: THE STAGE-1 SOLUTION IS Z-INDEPENDENT).
    for c in ("D", "F"):
        a_c = cal[f"alpha_{c}"]
        cal[f"Z_ss_{c}"] = ((rk[c] + cal[f"delta_{c}"]) / (mc[c] * a_c)) ** a_c


def _check_targets(cal, bk, rk, fm):
    # THE CALIBRATION TARGETS MUST HOLD EXACTLY AT THE SOLVED STEADY STATE.
    for c in ("D", "F"):
        rdep_c = cal[f"r_dep_{c}_target"]
        assert abs(bk[c]["theta_ss"] - cal[f"leverage_target_{c}"]) < 1e-6, \
            f"[{c}] leverage {bk[c]['theta_ss']:.6f} != target {cal[f'leverage_target_{c}']}"
        assert abs((rk[c] - rdep_c) - cal[f"credit_spread_target_{c}"]) < 1e-6, \
            f"[{c}] spread {(rk[c] - rdep_c):.6f} != target {cal[f'credit_spread_target_{c}']}"
    assert abs(fm["D"]["Y_ss"] - 1.0) < 1e-9, f"Y_ss_D={fm['D']['Y_ss']:.8f} != 1"
    assert abs(fm["F"]["Y_ss"] - 1.0) < 1e-9, f"Y_ss_F={fm['F']['Y_ss']:.8f} != 1"


def _deposit_market(cal, hh, beta, c, tol, verbose):
    # HOUSEHOLD SAVING MINUS BANK DEPOSIT SUPPLY FOR ONE COUNTRY AT A GUESSED beta.
    h = hh[c]
    try:
        c_ss, a_pol = solve_steady_state_household(
            h["a_grid"], h["Pi"], h["rdep"], h["y_e"], beta, cal[f"sigma_{c}"],
            cal[f"a_min_{c}"], tol, vN_ss=h["vN_ss"])
    except RuntimeError:
        # bracketing only needs the correct sign, so retry loose before failing
        c_ss, a_pol = solve_steady_state_household(
            h["a_grid"], h["Pi"], h["rdep"], h["y_e"], beta, cal[f"sigma_{c}"],
            cal[f"a_min_{c}"], max(1e-5, cal["tol_hh"] * 1e4), vN_ss=h["vN_ss"])
    D_ss = stationary_distribution(a_pol, h["a_grid"], h["Pi"], h["pi_e"], cal["tol_dist"])
    A_ss = aggregate_assets(D_ss, h["a_grid"])
    if verbose:
        print(f"  beta_{c}={beta:.8f}  A - Dep = {A_ss - h['Dep_supply']:.4e}")
    return A_ss - h["Dep_supply"], (c_ss, D_ss, A_ss)


def _household_problem(cal, c, fm, bk, gs, p_ss, mc):
    # EVERYTHING THE COUNTRY-c HOUSEHOLD PROBLEM NEEDS AT THE SOLVED FIRM/BANK SS.
    # The working-capital financing income stays with the BANK (inside div_ss), NOT the
    # household: as a household dividend it made the spread an intra-period transfer,
    # which under GHH turned the risk channel expansionary.
    e, Pi, pi_e = rouwenhorst(cal[f"rho_e_{c}"], cal[f"sigma_e_{c}"], n=cal[f"n_e_{c}"])
    Div_ss = (1 - mc) * fm["Y_ss"] + bk["div_ss"]
    P_CES = ces_price(p_ss, cal, c)
    return dict(e=e, Pi=Pi, pi_e=pi_e, a_grid=make_asset_grid(cal, country=c),
                rdep=cal[f"r_dep_{c}_target"], Dep_supply=bk["Dep_supply_ss"],
                vN_ss=cal[f"chi_{c}"] / (1 + 1 / cal[f"frisch_{c}"]),  # GHH v(N_ss = 1)
                P_CES=P_CES,
                y_e=fm["w_ss"] / P_CES * e + (Div_ss - gs["Tax_ss"]) / P_CES)


def _solve_stage2(cal, hh, verbose):
    # {beta_D, beta_F}: HOUSEHOLD DEPOSITS CLEAR THE BANK DEPOSIT SUPPLY.
    if verbose:
        print("\n=== Stage 2: deposit markets {beta_D, beta_F} ===")
    beta_upper_D = 1 / (1 + hh["D"]["rdep"]) - 1e-4   # keep rdep positive
    beta_D_ss = brentq(lambda b: _deposit_market(cal, hh, b, "D", cal["tol_hh"], verbose)[0],
                       0.5, beta_upper_D, xtol=1e-11)
    # union deposit market: the symmetric-SS doctrine pins beta_F = beta_D, where the
    # union clearing coincides with each national market and the cross-border position
    # is zero (an asymmetric SS would need a portfolio-split condition)
    beta_F_ss = beta_D_ss
    resid_F_chk, _ = _deposit_market(cal, hh, beta_F_ss, "F", cal["tol_hh"], verbose)
    assert abs(resid_F_chk) < 5e-6, (
        f"F deposit market off by {resid_F_chk:.2e} at beta_D — asymmetric SS? "
        "(the union stage 2 requires a symmetric SS)")
    _, sol_D = _deposit_market(cal, hh, beta_D_ss, "D", cal["tol_hh"], verbose)
    _, sol_F = _deposit_market(cal, hh, beta_F_ss, "F", cal["tol_hh"], verbose)
    return beta_D_ss, beta_F_ss, sol_D, sol_F


def _goods_market_check(cal, p_ss, C, P_CES, fm_D, verbose):
    # WALRAS CHECK ON THE D GOODS MARKET WITH THE TRUE STAGE-2 CONSUMPTION.
    IM_D_chk = import_demand(p_ss, C["D"], P_CES["D"], cal, "D")
    IM_F_chk = import_demand(p_ss, C["F"], P_CES["F"], cal, "F")
    NX_D_chk, _ = trade_balance(p_ss, IM_D_chk, IM_F_chk, cal)
    walras_D_chk = (fm_D["Y_ss"] - P_CES["D"] * C["D"] - fm_D["I_ss"]
                    - cal["G_D"] - NX_D_chk)
    if verbose:
        print(f"  SS goods-market check: walras_D = {walras_D_chk:.3e}")
    if abs(walras_D_chk) > 5e-6:
        print(f"  WARNING: SS goods market off by {walras_D_chk:.2e} — asymmetric SS "
              "calibration? (a symmetric SS is required)")


def solve_steady_state(cal, verbose=True):
    # SOLVE THE SYMMETRIC TWO-COUNTRY STEADY STATE (BANK CALIBRATION, STAGE 1, STAGE 2).
    _calibrate_bank_frictions(cal, verbose)
    mc = {c: markup_ss(cal, c) for c in ("D", "F")}
    gs = {c: govt_steady_state(cal, cal[f"r_dep_{c}_target"], c) for c in ("D", "F")}
    hold = _sovereign_holdings(cal)

    # stage 1 at the risk-free bond prices
    Q_b = {c: gs[c]["Q_B_ss"] for c in ("D", "F")}
    rk_D_ss, rk_F_ss, p_ss = _solve_stage1(cal, mc, Q_b, hold, verbose)
    rk = {"D": rk_D_ss, "F": rk_F_ss}
    _rescale_tfp(cal, rk, mc)

    # re-evaluate every SS object at the solution with the corrected Z
    books = {"D": (hold["b_D_D_ss"], hold["b_F_D_ss"]),
             "F": (hold["b_F_F_ss"], hold["b_D_F_ss_pc"])}
    Kap, fm, bk = {}, {}, {}
    for c in ("D", "F"):
        Kap[c], fm[c], bk[c] = _country_blocks(cal, rk[c], p_ss, mc[c], Q_b["D"], Q_b["F"],
                                               *books[c], c)
    # the traded bond prices are the IC-consistent ones, not the risk-free ones
    Q_bD_ss, Q_bF_ss = bk["D"]["Q_bdom_IC"], bk["F"]["Q_bdom_IC"]
    for c, Q in (("D", Q_bD_ss), ("F", Q_bF_ss)):
        gs[c] = dict(gs[c], Q_B_ss=Q,
                     Tax_ss=cal[f"G_{c}"] + cal[f"delta_b_{c}"] * cal[f"B_gov_{c}_ss"] * (1.0 - Q))
    cal["chi_D"] = fm["D"]["chi"]
    cal["chi_F"] = fm["F"]["chi"]
    cal["p_ss"]  = p_ss
    # foreign-bond FOC anchors (excess returns above the IC-required spread)
    cal["excess_return_F_D_ss"] = (bk["D"]["rb_for_ss"] - cal["r_dep_D_target"]
                                   - bk["D"]["IC_spread_for"])
    cal["excess_return_D_F_ss"] = (bk["F"]["rb_for_ss"] - cal["r_dep_F_target"]
                                   - bk["F"]["IC_spread_for"])
    _check_targets(cal, bk, rk, fm)
    if verbose:
        print(f"\nStage 1 solution: rk_D={rk_D_ss:.5f}  rk_F={rk_F_ss:.5f}  p_ss={p_ss:.4f}")
        print(f"  Z_ss (rescaled): D={cal['Z_ss_D']:.6f}  F={cal['Z_ss_F']:.6f}")
        print(f"  Kap_D={Kap['D']:.3f}  Kap_F={Kap['F']:.3f}  "
              f"n_ss_D={bk['D']['n_ss']:.4f}  n_ss_F={bk['F']['n_ss']:.4f}")

    # stage 2: household discount factors that clear the deposit markets
    hh = {c: _household_problem(cal, c, fm[c], bk[c], gs[c], p_ss, mc[c]) for c in ("D", "F")}
    beta_D_ss, beta_F_ss, (c_D_ss, D_D_ss, A_D_ss), (c_F_ss, D_F_ss, A_F_ss) = \
        _solve_stage2(cal, hh, verbose)
    C_D_ss = aggregate_consumption(D_D_ss, c_D_ss)
    C_F_ss = aggregate_consumption(D_F_ss, c_F_ss)
    if verbose:
        print(f"\nStage 2 solution: beta_D={beta_D_ss:.6f}  beta_F={beta_F_ss:.6f}")
        print(f"  A_D={A_D_ss:.4f}  Dep_supply_D={bk['D']['Dep_supply_ss']:.4f}")
        print(f"  A_F={A_F_ss:.4f}  Dep_supply_F={bk['F']['Dep_supply_ss']:.4f}")
    _goods_market_check(cal, p_ss, {"D": C_D_ss, "F": C_F_ss},
                        {c: hh[c]["P_CES"] for c in ("D", "F")}, fm["D"], verbose)

    return dict(
        beta_D_ss=beta_D_ss, beta_F_ss=beta_F_ss,
        rk_D_ss=rk_D_ss, rk_F_ss=rk_F_ss, p_ss=p_ss,
        Kap_D_ss=Kap["D"], Kap_F_ss=Kap["F"],
        ss_bank_D=bk["D"], ss_bank_F=bk["F"],
        ss_firm_D=fm["D"], ss_firm_F=fm["F"],
        gs_D=gs["D"], gs_F=gs["F"],
        e_D=hh["D"]["e"], Pi_D=hh["D"]["Pi"], e_F=hh["F"]["e"], Pi_F=hh["F"]["Pi"],
        a_grid_D=hh["D"]["a_grid"], a_grid_F=hh["F"]["a_grid"],
        c_D_ss=c_D_ss, D_D_ss=D_D_ss, c_F_ss=c_F_ss, D_F_ss=D_F_ss,
        A_D_ss=A_D_ss, C_D_ss=C_D_ss, A_F_ss=A_F_ss, C_F_ss=C_F_ss,
        Tax_D_ss=gs["D"]["Tax_ss"], Tax_F_ss=gs["F"]["Tax_ss"],
        Q_bD_ss=Q_bD_ss, Q_bF_ss=Q_bF_ss,
        b_D_D_ss=hold["b_D_D_ss"], b_F_D_ss=hold["b_F_D_ss"],
        b_F_F_ss=hold["b_F_F_ss"], b_D_F_ss=hold["b_D_F_ss"],
    )
