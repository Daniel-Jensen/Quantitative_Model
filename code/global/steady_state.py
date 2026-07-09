"""Two-country steady-state solver.

Two sequential stages:
Stage 1 — joint 3D solve over {rk_D, rk_F, p}:
Stage 2 — joint 2D solve over {beta_D, beta_F}:
"""
import numpy as np
from scipy.optimize import brentq, root

from rouwenhorst import rouwenhorst
from household import make_asset_grid, solve_steady_state_household
from distribution import stationary_distribution, aggregate_assets, aggregate_consumption
from firms import steady_state_firm, markup_ss
from capital import capital_demand
from bank import steady_state_bank
from government import govt_steady_state, hm_bond_price_ss
from trade import ces_price, import_demand, trade_balance


def solve_steady_state(cal, verbose=True):
    "For now only works for symmetric countries! Solve the two-country steady state"

    # ── Income processes and asset grids ──────────────────────────────────────
    # Income grids
    e_D, Pi_D, pi_D = rouwenhorst(cal["rho_e_D"], cal["sigma_e_D"], n=cal["n_e_D"])
    e_F, Pi_F, pi_F = rouwenhorst(cal["rho_e_F"], cal["sigma_e_F"], n=cal["n_e_F"])
    
    # Asset grids 
    a_grid_D = make_asset_grid(cal, country="D")
    a_grid_F = make_asset_grid(cal, country="F")

    # Markups 
    mc_D = markup_ss(cal, "D")
    mc_F = markup_ss(cal, "F")

    # Government SS
    rdep_D_tgt = cal["r_dep_D_target"]
    rdep_F_tgt = cal["r_dep_F_target"]
    #Dictionary of the government block
    gs_D = govt_steady_state(cal, rdep_D_tgt, "D")
    gs_F = govt_steady_state(cal, rdep_F_tgt, "F")
    Q_bD_ss = gs_D["Q_B_ss"]
    Q_bF_ss = gs_F["Q_B_ss"]
    b_F_D_ss = cal["b_F_D_ss"]
    b_D_F_ss = cal["b_D_F_ss"]
    B_gov_D  = cal["B_gov_D_ss"]
    B_gov_F  = cal["B_gov_F_ss"]

    b_D_D_ss = B_gov_D - b_D_F_ss
    b_F_F_ss = B_gov_F - b_F_D_ss

    ncalls = [0]

    def stage1_resid(x):
        rk_D, rk_F, p = x
        ncalls[0] += 1
        try:
            #Get capital values for each country given the current guess of rk_D and rk_F
            Kap_D = capital_demand(rk_D, mc_D, cal, country="D")
            Kap_F = capital_demand(rk_F, mc_F, cal, country="F")

            bk_D = steady_state_bank(cal, rk_D, Kap_D, Q_bD_ss, Q_bF_ss,
                                      b_D_D_ss, b_F_D_ss, p, country="D")
            bk_F = steady_state_bank(cal, rk_F, Kap_F, Q_bD_ss, Q_bF_ss,
                                      b_F_F_ss, b_D_F_ss, p, country="F")
            
            # Capital market residuals
            res_cap_D = (bk_D["n_ss_IC"] - bk_D["n_ss_ACCUM"]) / bk_D["n_ss_ACCUM"]
            res_cap_F = (bk_F["n_ss_IC"] - bk_F["n_ss_ACCUM"]) / bk_F["n_ss_ACCUM"]

            #Dictionary of the firm block 
            fm_D = steady_state_firm(cal, Kap_D, country="D")
            fm_F = steady_state_firm(cal, Kap_F, country="F")
            

            C_D_approx = fm_D["C_ss"]
            C_F_approx = fm_F["C_ss"]
            # Price of consumption basket
            P_CES_D = ces_price(p, cal, country="D")
            P_CES_F = ces_price(p, cal, country="F")
            # Import demand for each country
            IM_D = import_demand(p, C_D_approx, P_CES_D, cal, country="D")
            IM_F = import_demand(p, C_F_approx, P_CES_F, cal, country="F")
            #Net exports 
            NX_D, _ = trade_balance(p, IM_D, IM_F)

            # Market prices and returns
            Q_bD_mkt = bk_D["Q_bdom_IC"]
            Q_bF_mkt = bk_F["Q_bdom_IC"]
            rb_D_mkt = cal["r_dep_D_target"] + bk_D["IC_spread_dom"]
            rb_F_mkt = cal["r_dep_F_target"] + bk_F["IC_spread_dom"]
            income_in_D  = rb_F_mkt * p * Q_bF_mkt * b_F_D_ss   #coupon is priced in price of the bond. 
            income_out_D = rb_D_mkt * Q_bD_mkt * b_D_F_ss       

            res_ext = (NX_D + income_in_D - income_out_D) / fm_D["Y_ss"]

        except (RuntimeError, ValueError, FloatingPointError, ZeroDivisionError):
            return [1e3, 1e3, 1e3]

        if verbose:
            print(f"  stage1 call {ncalls[0]:3d}: rk_D={rk_D:.5f}  rk_F={rk_F:.5f}  "
                  f"p={p:.4f}  |resid|=[{res_cap_D:.3e},{res_cap_F:.3e},{res_ext:.3e}]")

        return [res_cap_D, res_cap_F, res_ext]

    ### -- ACTUAL SOLVER FOR THE FIRST STAGE ---
    x0 = [cal["rk_D_guess"], cal["rk_F_guess"], 1.0]
    if verbose:
        print("=== Stage 1: capital markets + external balance {rk_D, rk_F, p} ===")
    sol1 = root(stage1_resid, x0, method="hybr",
                options={"xtol": cal["tol_mkt"], "maxfev": 2000})
    if not sol1.success and verbose:
        print(f"  Warning: stage1 hybr did not flag success (resid={np.max(np.abs(sol1.fun)):.2e})")
    rk_D_ss, rk_F_ss, p_ss = sol1.x

    # Re-evaluate at solution to get all SS objects
    # capital
    Kap_D_ss = capital_demand(rk_D_ss, mc_D, cal, "D")
    Kap_F_ss = capital_demand(rk_F_ss, mc_F, cal, "F")

    #bank dict
    bk_D_ss  = steady_state_bank(cal, rk_D_ss, Kap_D_ss, Q_bD_ss, Q_bF_ss,
                                  b_D_D_ss, b_F_D_ss, p_ss, "D")
    bk_F_ss  = steady_state_bank(cal, rk_F_ss, Kap_F_ss, Q_bD_ss, Q_bF_ss,
                               b_F_F_ss, b_D_F_ss, p_ss, "F")

    #firm dict 
    fm_D_ss  = steady_state_firm(cal, Kap_D_ss, "D")
    fm_F_ss  = steady_state_firm(cal, Kap_F_ss, "F")

    # Two price inconsistency -> use IC-consistent bond prices 
    # Q_B_risk_free = δ_b / (rdep + δ_b) vs Q_bdom_IC = δ_b / (rdep + δ_b + IC_spread)

    Q_bD_IC = bk_D_ss["Q_bdom_IC"]   
    Q_bF_IC = bk_F_ss["Q_bdom_IC"]   
    gs_D = dict(gs_D)   # mutable copy
    gs_F = dict(gs_F)
    gs_D["Q_B_ss"] = Q_bD_IC
    gs_F["Q_B_ss"] = Q_bF_IC
    gs_D["Tax_ss"] = (cal["G_D"]
                      + cal["delta_b_D"] * cal["B_gov_D_ss"] * (1.0 - Q_bD_IC))
    gs_F["Tax_ss"] = (cal["G_F"]
                      + cal["delta_b_F"] * cal["B_gov_F_ss"] * (1.0 - Q_bF_IC))
    Q_bD_ss = Q_bD_IC 
    Q_bF_ss = Q_bF_IC

    # Store GHH chi and p back into cal
    cal["chi_D"] = fm_D_ss["chi"]
    cal["chi_F"] = fm_F_ss["chi"]
    cal["p_ss"]  = p_ss

    # Compute Excessive Returns 
    ic_spread_bF_D = bk_D_ss["IC_spread_for"]   # = lambda_bF_D * mu_D / Omega_D
    ic_spread_bD_F = bk_F_ss["IC_spread_for"]   # = lambda_bD_F * mu_F / Omega_F
    cal["excess_return_F_D_ss"] = bk_D_ss["rb_for_ss"] - rdep_D_tgt - ic_spread_bF_D
    cal["excess_return_D_F_ss"] = bk_F_ss["rb_for_ss"] - rdep_F_tgt - ic_spread_bD_F

    if verbose:
        print(f"\nStage 1 solution: rk_D={rk_D_ss:.5f}  rk_F={rk_F_ss:.5f}  p_ss={p_ss:.4f}")
        print(f"  Kap_D={Kap_D_ss:.3f}  Kap_F={Kap_F_ss:.3f}")
        print(f"  n_ss_D={bk_D_ss['n_ss']:.4f}  n_ss_F={bk_F_ss['n_ss']:.4f}")
        print(f"  IC resid D={(bk_D_ss['n_ss_IC']-bk_D_ss['n_ss_ACCUM'])/bk_D_ss['n_ss_ACCUM']:.2e}")
        print(f"  IC resid F={(bk_F_ss['n_ss_IC']-bk_F_ss['n_ss_ACCUM'])/bk_F_ss['n_ss_ACCUM']:.2e}")

    # ── Stage 2: deposit markets — solve {beta_D, beta_F} ────────────────────
    # Government taxes and dividends
    Tax_D_ss = gs_D["Tax_ss"]
    Tax_F_ss = gs_F["Tax_ss"]
    Div_D_ss = (1 - mc_D) * fm_D_ss["Y_ss"] + bk_D_ss["div_ss"]
    Div_F_ss = (1 - mc_F) * fm_F_ss["Y_ss"] + bk_F_ss["div_ss"]

    P_CES_D_ss = ces_price(p_ss, cal, "D")
    P_CES_F_ss = ces_price(p_ss, cal, "F")

    # GHH SS labour disutility
    vN_D_ss = cal["chi_D"] * 1.0 ** (1 + 1 / cal["frisch_D"]) / (1 + 1 / cal["frisch_D"])
    vN_F_ss = cal["chi_F"] * 1.0 ** (1 + 1 / cal["frisch_F"]) / (1 + 1 / cal["frisch_F"])

    if verbose:
        print("\n=== Stage 2: deposit markets {beta_D, beta_F} ===")


    _tol_tight  = cal["tol_hh"]
    _tol_sign   = max(1e-5, _tol_tight * 1e4)   # fallback: only need correct sign

    def _egm_solve(a_grid, Pi, rdep_tgt, y_e, beta, country, tol):
        return solve_steady_state_household(
            a_grid, Pi, rdep_tgt, y_e, beta,
            cal[f"sigma_{country}"], cal[f"a_min_{country}"], tol,
            vN_ss=(vN_D_ss if country == "D" else vN_F_ss),
        )

    def _egm_solve_robust(a_grid, Pi, rdep_tgt, y_e, beta, country, tol):
        #EGM SOLVE WITH A FALLBACK ON LOWER TOLERENCE IF THE FIRST ATTEMPT FAILS
        try:
            return _egm_solve(a_grid, Pi, rdep_tgt, y_e, beta, country, tol)
        except RuntimeError:
            return _egm_solve(a_grid, Pi, rdep_tgt, y_e, beta, country, _tol_sign)

    def deposit_resid_D(beta_D, tol=_tol_tight):
        # CALCULATE THE DEPOSITS MARKET RESIDUAL FROM THE STEADY STATE HOUSEHOLD PROBLEM
        w_real_D = fm_D_ss["w_ss"] / P_CES_D_ss
        y_e_D = w_real_D * e_D + (Div_D_ss - Tax_D_ss) / P_CES_D_ss
        c_ss_D, a_pol_D = _egm_solve_robust(a_grid_D, Pi_D, rdep_D_tgt, y_e_D, beta_D, "D", tol)
        D_ss_D = stationary_distribution(a_pol_D, a_grid_D, Pi_D, pi_D, cal["tol_dist"])
        A_ss_D = aggregate_assets(D_ss_D, a_grid_D)
        val = A_ss_D - bk_D_ss["Dep_supply_ss"]
        if verbose:
            print(f"  beta_D={beta_D:.8f}  A_D − Dep_D = {val:.4e}")
        return val, (c_ss_D, a_pol_D, D_ss_D, A_ss_D, y_e_D)

    def deposit_resid_F(beta_F, tol=_tol_tight):
        # CALCULATE THE DEPOSITS MARKET RESIDUAL FROM THE STEADY STATE HOUSEHOLD PROBLEM
        w_real_F = fm_F_ss["w_ss"] / P_CES_F_ss
        y_e_F = w_real_F * e_F + (Div_F_ss - Tax_F_ss) / P_CES_F_ss
        c_ss_F, a_pol_F = _egm_solve_robust(a_grid_F, Pi_F, rdep_F_tgt, y_e_F, beta_F, "F", tol)
        D_ss_F = stationary_distribution(a_pol_F, a_grid_F, Pi_F, pi_F, cal["tol_dist"])
        A_ss_F = aggregate_assets(D_ss_F, a_grid_F)
        val = A_ss_F - bk_F_ss["Dep_supply_ss"]
        if verbose:
            print(f"  beta_F={beta_F:.8f}  A_F − Dep_F = {val:.4e}")
        return val, (c_ss_F, a_pol_F, D_ss_F, A_ss_F, y_e_F)
    
    #calculate the upper bound for beta to ensure that the deposit rate is positive
    beta_upper_D = 1 / (1 + rdep_D_tgt) - 1e-4
    beta_upper_F = 1 / (1 + rdep_F_tgt) - 1e-4

    #ACTUALLY SOLVE FOR THE STEADY STATE DISCOUNT FACTORS USING THE BRENTQ METHOD
    beta_D_ss = brentq(lambda b: deposit_resid_D(b)[0], 0.5, beta_upper_D, xtol=1e-11)
    beta_F_ss = brentq(lambda b: deposit_resid_F(b)[0], 0.5, beta_upper_F, xtol=1e-11)

    # Re-solve at tight tolerance to get accurate policy functions
    _, (c_D_ss, a_pol_D_ss, D_D_ss, A_D_ss, y_e_D) = deposit_resid_D(beta_D_ss, tol=cal["tol_hh"])
    _, (c_F_ss, a_pol_F_ss, D_F_ss, A_F_ss, y_e_F) = deposit_resid_F(beta_F_ss, tol=cal["tol_hh"])

    C_D_ss = aggregate_consumption(D_D_ss, c_D_ss)
    C_F_ss = aggregate_consumption(D_F_ss, c_F_ss)

    P_CES_D_ss = ces_price(p_ss, cal, "D")
    P_CES_F_ss = ces_price(p_ss, cal, "F")

    if verbose:
        print(f"\nStage 2 solution: beta_D={beta_D_ss:.6f}  beta_F={beta_F_ss:.6f}")
        print(f"  A_D={A_D_ss:.4f}  Dep_supply_D={bk_D_ss['Dep_supply_ss']:.4f}")
        print(f"  A_F={A_F_ss:.4f}  Dep_supply_F={bk_F_ss['Dep_supply_ss']:.4f}")

    # ── Post-solve goods-market check with TRUE stage-2 consumption ──────────
    IM_D_chk = import_demand(p_ss, C_D_ss, P_CES_D_ss, cal, "D")
    IM_F_chk = import_demand(p_ss, C_F_ss, P_CES_F_ss, cal, "F")
    NX_D_chk, _ = trade_balance(p_ss, IM_D_chk, IM_F_chk)
    walras_D_chk = fm_D_ss["Y_ss"] - P_CES_D_ss * C_D_ss - fm_D_ss["I_ss"] - cal["G_D"] - NX_D_chk
    if verbose:
        print(f"  SS goods-market check: walras_D = {walras_D_chk:.3e}")
    if abs(walras_D_chk) > 5e-6:
        print(f"  WARNING: SS goods market off by {walras_D_chk:.2e} — asymmetric SS "
              "calibration? (see solve_steady_state docstring)")

    return dict(
        # Solved scalars
        beta_D_ss=beta_D_ss, beta_F_ss=beta_F_ss,
        rk_D_ss=rk_D_ss, rk_F_ss=rk_F_ss, p_ss=p_ss,
        Kap_D_ss=Kap_D_ss, Kap_F_ss=Kap_F_ss,
        # Steady-state blocks
        ss_bank_D=bk_D_ss, ss_bank_F=bk_F_ss,
        ss_firm_D=fm_D_ss, ss_firm_F=fm_F_ss,
        gs_D=gs_D, gs_F=gs_F,
        # Income processes
        e_D=e_D, Pi_D=Pi_D, pi_D=pi_D,
        e_F=e_F, Pi_F=Pi_F, pi_F=pi_F,
        # Asset grids
        a_grid_D=a_grid_D, a_grid_F=a_grid_F,
        # Household SS solutions
        c_D_ss=c_D_ss, a_pol_D_ss=a_pol_D_ss, D_D_ss=D_D_ss,
        c_F_ss=c_F_ss, a_pol_F_ss=a_pol_F_ss, D_F_ss=D_F_ss,
        # Aggregates
        A_D_ss=A_D_ss, C_D_ss=C_D_ss, y_e_D=y_e_D,
        A_F_ss=A_F_ss, C_F_ss=C_F_ss, y_e_F=y_e_F,
        Div_D_ss=Div_D_ss, Tax_D_ss=Tax_D_ss,
        Div_F_ss=Div_F_ss, Tax_F_ss=Tax_F_ss,
        P_CES_D_ss=P_CES_D_ss, P_CES_F_ss=P_CES_F_ss,
        # Bond SS quantities
        Q_bD_ss=Q_bD_ss, Q_bF_ss=Q_bF_ss,
        b_D_D_ss=b_D_D_ss, b_F_D_ss=b_F_D_ss,
        b_F_F_ss=b_F_F_ss, b_D_F_ss=b_D_F_ss,
        B_gov_D=B_gov_D, B_gov_F=B_gov_F,
    )
