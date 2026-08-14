# TWO-COUNTRY STEADY-STATE SOLVER: STAGE 1 {rk_D, rk_F, p}, STAGE 2 {beta_D, beta_F}.
# The steady state must be SYMMETRIC: country asymmetries enter through shocks
# only. An asymmetric SS shifts p_ss off 1 and opens an O(1e-4) goods-market
# wedge, because p is only weakly identified by external balance at a trade
# elasticity of 0.5.
import numpy as np
from scipy.optimize import brentq, root

from blocks.rouwenhorst import rouwenhorst
from blocks.household import make_asset_grid, solve_steady_state_household
from blocks.distribution import stationary_distribution, aggregate_assets, aggregate_consumption
from blocks.firms import steady_state_firm, markup_ss
from blocks.capital import capital_demand
from blocks.bank import steady_state_bank, calibrate_bank_targets
from blocks.government import govt_steady_state
from blocks.trade import ces_price, import_demand, trade_balance


def solve_steady_state(cal, verbose=True):
    # SOLVE THE SYMMETRIC TWO-COUNTRY STEADY STATE (BANK CALIBRATION, STAGE 1, STAGE 2).

    # bank agency friction: solve lambda and omega_ent to hit the leverage and
    # spread targets, then write them back into cal
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

    e_D, Pi_D, pi_e_D = rouwenhorst(cal["rho_e_D"], cal["sigma_e_D"], n=cal["n_e_D"])
    e_F, Pi_F, pi_e_F = rouwenhorst(cal["rho_e_F"], cal["sigma_e_F"], n=cal["n_e_F"])
    a_grid_D = make_asset_grid(cal, country="D")
    a_grid_F = make_asset_grid(cal, country="F")
    mc_D = markup_ss(cal, "D")
    mc_F = markup_ss(cal, "F")

    rdep_D_tgt = cal["r_dep_D_target"]
    rdep_F_tgt = cal["r_dep_F_target"]
    gs_D = govt_steady_state(cal, rdep_D_tgt, "D")
    gs_F = govt_steady_state(cal, rdep_F_tgt, "F")
    Q_bD_ss = gs_D["Q_B_ss"]
    Q_bF_ss = gs_F["Q_B_ss"]
    b_F_D_ss = cal["b_F_D_ss"]
    b_D_F_ss = cal["b_D_F_ss"]
    b_D_D_ss = cal["B_gov_D_ss"] - b_D_F_ss
    b_F_F_ss = cal["B_gov_F_ss"] - b_F_D_ss

    ncalls = [0]

    def stage1_resid(x):
        # CAPITAL-MARKET (n_IC = n_ACCUM) + EXTERNAL-BALANCE RESIDUALS IN {rk_D, rk_F, p}.
        rk_D, rk_F, p = x
        ncalls[0] += 1
        try:
            Kap_D = capital_demand(rk_D, mc_D, cal, country="D")
            Kap_F = capital_demand(rk_F, mc_F, cal, country="F")

            bk_D = steady_state_bank(cal, rk_D, Kap_D, Q_bD_ss, Q_bF_ss,
                                     b_D_D_ss, b_F_D_ss, p, country="D")
            bk_F = steady_state_bank(cal, rk_F, Kap_F, Q_bD_ss, Q_bF_ss,
                                     b_F_F_ss, b_D_F_ss, p, country="F")

            res_cap_D = (bk_D["n_ss_IC"] - bk_D["n_ss_ACCUM"]) / bk_D["n_ss_ACCUM"]
            res_cap_F = (bk_F["n_ss_IC"] - bk_F["n_ss_ACCUM"]) / bk_F["n_ss_ACCUM"]

            fm_D = steady_state_firm(cal, Kap_D, country="D")
            fm_F = steady_state_firm(cal, Kap_F, country="F")

            P_CES_D = ces_price(p, cal, country="D")
            P_CES_F = ces_price(p, cal, country="F")
            IM_D = import_demand(p, fm_D["C_ss"], P_CES_D, cal, country="D")
            IM_F = import_demand(p, fm_F["C_ss"], P_CES_F, cal, country="F")
            NX_D, _ = trade_balance(p, IM_D, IM_F)

            # cross-border coupon income, already priced into Q
            rb_D_mkt = rdep_D_tgt + bk_D["IC_spread_dom"]
            rb_F_mkt = rdep_F_tgt + bk_F["IC_spread_dom"]
            income_in_D  = rb_F_mkt * p * bk_F["Q_bdom_IC"] * b_F_D_ss
            income_out_D = rb_D_mkt * bk_D["Q_bdom_IC"] * b_D_F_ss

            res_ext = (NX_D + income_in_D - income_out_D) / fm_D["Y_ss"]

        except (RuntimeError, ValueError, FloatingPointError, ZeroDivisionError):
            return [1e3, 1e3, 1e3]

        if verbose:
            print(f"  stage1 call {ncalls[0]:3d}: rk_D={rk_D:.5f}  rk_F={rk_F:.5f}  "
                  f"p={p:.4f}  |resid|=[{res_cap_D:.3e},{res_cap_F:.3e},{res_ext:.3e}]")

        return [res_cap_D, res_cap_F, res_ext]

    if verbose:
        print("=== Stage 1: capital markets + external balance {rk_D, rk_F, p} ===")
    sol1 = root(stage1_resid, [cal["rk_D_guess"], cal["rk_F_guess"], 1.0],
                method="hybr", options={"xtol": cal["tol_mkt"], "maxfev": 2000})
    if not sol1.success and verbose:
        print(f"  Warning: stage1 hybr did not flag success "
              f"(resid={np.max(np.abs(sol1.fun)):.2e})")
    # HARD GUARD: the capital-market residuals ARE the bank n_IC/n_ACCUM identity. A large
    # value means the (f, spread, leverage) targets sit on the franchise fold's UPPER root
    # while the dynamics rest on the least root -- an inconsistent SS (see main.py). Fail
    # loudly rather than propagate a silently-off steady state into the projection solve.
    assert np.max(np.abs(sol1.fun[:2])) < 1e-6, (
        f"stage-1 n_IC/n_ACCUM inconsistent (max|res_cap|={np.max(np.abs(sol1.fun[:2])):.2e}): "
        f"the (f, spread, leverage) targets are fold-blocked -- raise f until leverage is the "
        f"least root (f>=~0.14 at spread 720bp).")
    rk_D_ss, rk_F_ss, p_ss = sol1.x

    # rescale TFP so Y_ss = 1 (exact: the stage-1 solution is Z-independent)
    for c, rk_c, mc_c in (("D", rk_D_ss, mc_D), ("F", rk_F_ss, mc_F)):
        a_c = cal[f"alpha_{c}"]
        cal[f"Z_ss_{c}"] = ((rk_c + cal[f"delta_{c}"]) / (mc_c * a_c)) ** a_c

    # re-evaluate every SS object at the solution with the corrected Z
    Kap_D_ss = capital_demand(rk_D_ss, mc_D, cal, "D")
    Kap_F_ss = capital_demand(rk_F_ss, mc_F, cal, "F")
    bk_D_ss  = steady_state_bank(cal, rk_D_ss, Kap_D_ss, Q_bD_ss, Q_bF_ss,
                                 b_D_D_ss, b_F_D_ss, p_ss, "D")
    bk_F_ss  = steady_state_bank(cal, rk_F_ss, Kap_F_ss, Q_bD_ss, Q_bF_ss,
                                 b_F_F_ss, b_D_F_ss, p_ss, "F")
    fm_D_ss  = steady_state_firm(cal, Kap_D_ss, "D")
    fm_F_ss  = steady_state_firm(cal, Kap_F_ss, "F")

    # the traded bond prices are the IC-consistent ones, not the risk-free ones
    Q_bD_ss = bk_D_ss["Q_bdom_IC"]
    Q_bF_ss = bk_F_ss["Q_bdom_IC"]
    gs_D = dict(gs_D, Q_B_ss=Q_bD_ss,
                Tax_ss=cal["G_D"] + cal["delta_b_D"] * cal["B_gov_D_ss"] * (1.0 - Q_bD_ss))
    gs_F = dict(gs_F, Q_B_ss=Q_bF_ss,
                Tax_ss=cal["G_F"] + cal["delta_b_F"] * cal["B_gov_F_ss"] * (1.0 - Q_bF_ss))

    cal["chi_D"] = fm_D_ss["chi"]
    cal["chi_F"] = fm_F_ss["chi"]
    cal["p_ss"]  = p_ss

    # foreign-bond FOC anchors (excess returns above the IC-required spread)
    cal["excess_return_F_D_ss"] = (bk_D_ss["rb_for_ss"] - rdep_D_tgt
                                   - bk_D_ss["IC_spread_for"])
    cal["excess_return_D_F_ss"] = (bk_F_ss["rb_for_ss"] - rdep_F_tgt
                                   - bk_F_ss["IC_spread_for"])

    for c, bk, rk, rdep_c in (("D", bk_D_ss, rk_D_ss, rdep_D_tgt),
                              ("F", bk_F_ss, rk_F_ss, rdep_F_tgt)):
        assert abs(bk["theta_ss"] - cal[f"leverage_target_{c}"]) < 1e-6, \
            f"[{c}] leverage {bk['theta_ss']:.6f} != target {cal[f'leverage_target_{c}']}"
        assert abs((rk - rdep_c) - cal[f"credit_spread_target_{c}"]) < 1e-6, \
            f"[{c}] spread {(rk - rdep_c):.6f} != target {cal[f'credit_spread_target_{c}']}"
    assert abs(fm_D_ss["Y_ss"] - 1.0) < 1e-9, f"Y_ss_D={fm_D_ss['Y_ss']:.8f} != 1"
    assert abs(fm_F_ss["Y_ss"] - 1.0) < 1e-9, f"Y_ss_F={fm_F_ss['Y_ss']:.8f} != 1"

    if verbose:
        print(f"\nStage 1 solution: rk_D={rk_D_ss:.5f}  rk_F={rk_F_ss:.5f}  p_ss={p_ss:.4f}")
        print(f"  Z_ss (rescaled): D={cal['Z_ss_D']:.6f}  F={cal['Z_ss_F']:.6f}")
        print(f"  Kap_D={Kap_D_ss:.3f}  Kap_F={Kap_F_ss:.3f}  "
              f"n_ss_D={bk_D_ss['n_ss']:.4f}  n_ss_F={bk_F_ss['n_ss']:.4f}")

    # stage 2: dividends include the working-capital financing income
    # (r_wc_ss = rdep + the credit-spread target, a constant at the SS)
    wc_D_ss = (cal["zeta_wc_D"]
               * (rdep_D_tgt + cal["credit_spread_target_D"]) * fm_D_ss["w_ss"])
    wc_F_ss = (cal["zeta_wc_F"]
               * (rdep_F_tgt + cal["credit_spread_target_F"]) * fm_F_ss["w_ss"])
    Div_D_ss = (1 - mc_D) * fm_D_ss["Y_ss"] + bk_D_ss["div_ss"] + wc_D_ss
    Div_F_ss = (1 - mc_F) * fm_F_ss["Y_ss"] + bk_F_ss["div_ss"] + wc_F_ss

    P_CES_D_ss = ces_price(p_ss, cal, "D")
    P_CES_F_ss = ces_price(p_ss, cal, "F")
    vN_D_ss = cal["chi_D"] / (1 + 1 / cal["frisch_D"])   # GHH disutility at N_ss = 1
    vN_F_ss = cal["chi_F"] / (1 + 1 / cal["frisch_F"])

    if verbose:
        print("\n=== Stage 2: deposit markets {beta_D, beta_F} ===")

    def _deposit_resid(beta, country, tol):
        # DEPOSIT-MARKET RESIDUAL A - Dep_supply FOR ONE COUNTRY AT A GUESSED beta.
        if country == "D":
            a_grid, Pi, pi_e, e = a_grid_D, Pi_D, pi_e_D, e_D
            rdep_tgt, vN_ss, bk = rdep_D_tgt, vN_D_ss, bk_D_ss
            y_e = fm_D_ss["w_ss"] / P_CES_D_ss * e + (Div_D_ss - gs_D["Tax_ss"]) / P_CES_D_ss
        else:
            a_grid, Pi, pi_e, e = a_grid_F, Pi_F, pi_e_F, e_F
            rdep_tgt, vN_ss, bk = rdep_F_tgt, vN_F_ss, bk_F_ss
            y_e = fm_F_ss["w_ss"] / P_CES_F_ss * e + (Div_F_ss - gs_F["Tax_ss"]) / P_CES_F_ss

        try:
            c_ss, a_pol = solve_steady_state_household(
                a_grid, Pi, rdep_tgt, y_e, beta, cal[f"sigma_{country}"],
                cal[f"a_min_{country}"], tol, vN_ss=vN_ss)
        except RuntimeError:
            # bracketing only needs the correct sign, so retry loose before failing
            c_ss, a_pol = solve_steady_state_household(
                a_grid, Pi, rdep_tgt, y_e, beta, cal[f"sigma_{country}"],
                cal[f"a_min_{country}"], max(1e-5, cal["tol_hh"] * 1e4), vN_ss=vN_ss)

        D_ss = stationary_distribution(a_pol, a_grid, Pi, pi_e, cal["tol_dist"])
        A_ss = aggregate_assets(D_ss, a_grid)
        if verbose:
            print(f"  beta_{country}={beta:.8f}  A - Dep = {A_ss - bk['Dep_supply_ss']:.4e}")
        return A_ss - bk["Dep_supply_ss"], (c_ss, D_ss, A_ss)

    beta_upper_D = 1 / (1 + rdep_D_tgt) - 1e-4   # keep rdep positive
    beta_D_ss = brentq(lambda b: _deposit_resid(b, "D", cal["tol_hh"])[0],
                       0.5, beta_upper_D, xtol=1e-11)
    # Union deposit market: the symmetric-SS doctrine pins beta_F = beta_D. At a
    # symmetric SS the union clearing coincides with each national market and the
    # cross-border deposit position is zero; an asymmetric SS would instead need
    # the union clearing plus a portfolio-split condition.
    beta_F_ss = beta_D_ss
    resid_F_chk, _ = _deposit_resid(beta_F_ss, "F", cal["tol_hh"])
    assert abs(resid_F_chk) < 5e-6, (
        f"F deposit market off by {resid_F_chk:.2e} at beta_D — asymmetric SS? "
        "(the union stage 2 requires a symmetric SS)")

    _, (c_D_ss, D_D_ss, A_D_ss) = _deposit_resid(beta_D_ss, "D", cal["tol_hh"])
    _, (c_F_ss, D_F_ss, A_F_ss) = _deposit_resid(beta_F_ss, "F", cal["tol_hh"])
    C_D_ss = aggregate_consumption(D_D_ss, c_D_ss)
    C_F_ss = aggregate_consumption(D_F_ss, c_F_ss)

    if verbose:
        print(f"\nStage 2 solution: beta_D={beta_D_ss:.6f}  beta_F={beta_F_ss:.6f}")
        print(f"  A_D={A_D_ss:.4f}  Dep_supply_D={bk_D_ss['Dep_supply_ss']:.4f}")
        print(f"  A_F={A_F_ss:.4f}  Dep_supply_F={bk_F_ss['Dep_supply_ss']:.4f}")

    # goods-market check with the true stage-2 consumption
    IM_D_chk = import_demand(p_ss, C_D_ss, P_CES_D_ss, cal, "D")
    IM_F_chk = import_demand(p_ss, C_F_ss, P_CES_F_ss, cal, "F")
    NX_D_chk, _ = trade_balance(p_ss, IM_D_chk, IM_F_chk)
    walras_D_chk = (fm_D_ss["Y_ss"] - P_CES_D_ss * C_D_ss - fm_D_ss["I_ss"]
                    - cal["G_D"] - NX_D_chk)
    if verbose:
        print(f"  SS goods-market check: walras_D = {walras_D_chk:.3e}")
    if abs(walras_D_chk) > 5e-6:
        print(f"  WARNING: SS goods market off by {walras_D_chk:.2e} — asymmetric SS "
              "calibration? (a symmetric SS is required)")

    return dict(
        beta_D_ss=beta_D_ss, beta_F_ss=beta_F_ss,
        rk_D_ss=rk_D_ss, rk_F_ss=rk_F_ss, p_ss=p_ss,
        Kap_D_ss=Kap_D_ss, Kap_F_ss=Kap_F_ss,
        ss_bank_D=bk_D_ss, ss_bank_F=bk_F_ss,
        ss_firm_D=fm_D_ss, ss_firm_F=fm_F_ss,
        gs_D=gs_D, gs_F=gs_F,
        e_D=e_D, Pi_D=Pi_D, e_F=e_F, Pi_F=Pi_F,
        a_grid_D=a_grid_D, a_grid_F=a_grid_F,
        c_D_ss=c_D_ss, D_D_ss=D_D_ss, c_F_ss=c_F_ss, D_F_ss=D_F_ss,
        A_D_ss=A_D_ss, C_D_ss=C_D_ss, A_F_ss=A_F_ss, C_F_ss=C_F_ss,
        Tax_D_ss=gs_D["Tax_ss"], Tax_F_ss=gs_F["Tax_ss"],
        Q_bD_ss=Q_bD_ss, Q_bF_ss=Q_bF_ss,
        b_D_D_ss=b_D_D_ss, b_F_D_ss=b_F_D_ss,
        b_F_F_ss=b_F_F_ss, b_D_F_ss=b_D_F_ss,
    )
