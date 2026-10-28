# TWO-COUNTRY STEADY-STATE SOLVER: THE CAPITAL FOC, THE EXTERNAL BALANCE IN p, THEN {beta_D, beta_F}.
# The steady state is ASYMMETRIC (2026-10-19 calibration): the banks' own-sovereign exposures
# differ, and F banks hold a fifth of the D stock against D banks' sliver of the F stock, so D pays
# net coupons abroad, runs a trade surplus and p_ss != 1. Every block below holds at any p:
# consumption is the national-income identity P_CES*C = Y - I - G + NFI, deposits clear in
# consumption-basket units, chi is set at w/P_CES and the government is built from the solved debt.
# p is weakly identified: a trade elasticity of 0.5 sits on the Marshall-Lerner boundary, so the
# trade balance moves ~0.012 of output per unit of log p. It is bracketed, never Newton-solved.
#
# THE ANATOMY (solve_steady_state):
#   _calibrate_bank_frictions   lambda, omega_ent from the leverage and spread instruments
#   rk = r + spread             the capital FOC at the SS (single lambda: every class earns it)
#   _rescale_tfp                Z_ss so that Y_ss = 1
#   _solve_external_balance     p, with the sovereign books sized to the exposure instruments
#   _country_blocks             capital -> firm -> working-capital loan -> bank, per country
#   _friction_references        cross-border friction anchors at which the SS books are optimal
#   _solve_stage2               {beta_D, beta_F}: each country's deposits clear its own bank
#   _goods_market_check         both goods markets and the current account, asserted
import numpy as np
from scipy.optimize import brentq

from global_projection.blocks.rouwenhorst import rouwenhorst
from global_projection.blocks.household import make_asset_grid, solve_steady_state_household
from global_projection.blocks.distribution import stationary_distribution, aggregate_assets, aggregate_consumption
from global_projection.blocks.firms import steady_state_firm, markup_ss
from global_projection.blocks.capital import capital_demand
from global_projection.blocks.bank import steady_state_bank, calibrate_bank_targets
from global_projection.blocks.government import govt_steady_state
from global_projection.blocks.trade import ces_price, import_demand, trade_balance, size_ratio

_OTHER = {"D": "F", "F": "D"}
# the bracket for p; the external balance is monotone on it at this calibration
P_BRACKET = (0.5, 2.0)


def _calibrate_bank_frictions(cal, verbose, share_for=None):
    # SINGLE lambda AND omega_ent PER COUNTRY FROM LEVERAGE + SPREAD, WRITTEN INTO cal.
    # share_for: each bank's foreign-bond share of assets; that leg earns its issuer's spread, so
    # omega_ent is re-solved once the books are sized (lambda does not depend on them)
    for c in ("D", "F"):
        lam, om, *_ = calibrate_bank_targets(
            cal[f"beta_inter_{c}"], cal[f"f_{c}"], cal[f"r_dep_{c}_target"],
            cal[f"leverage_target_{c}"], cal[f"credit_spread_target_{c}"],
            spread_for=cal[f"credit_spread_target_{_OTHER[c]}"],
            share_for=0.0 if share_for is None else share_for[c],
            mu_eps=float(cal.get("mu_eps", 0.0)))
        cal[f"lambda_K_{c}"] = cal[f"lambda_bD_{c}"] = cal[f"lambda_bF_{c}"] = lam
        cal[f"omega_ent_{c}"] = om
        if verbose and share_for is not None:
            print(f"[bank-cal {c}] lambda={lam:.6f}  omega_ent={om:.6f}  "
                  f"(theta={cal[f'leverage_target_{c}']:.4f}, "
                  f"spread={cal[f'credit_spread_target_{c}']*4e4:.2f} bps/yr, "
                  f"foreign-bond share {100*share_for[c]:.3f}% of assets)")


def _sovereign_holdings(cal):
    # SS BOND BOOKS, CARRIED IN THE ISSUER'S PER-CAPITA UNITS.
    # b_D_F_ss is the slice of D's own per-capita stock held abroad, so the F bank's
    # per-capita book is b_D_F_ss/sz (sz = size_F/size_D); b_F_D_ss is the D bank's
    # per-capita F book, b_F_D_ss/sz of each F capita's stock.
    sz = size_ratio(cal)
    b_F_D_ss = cal["b_F_D_ss"]           # D's per-capita holding of the F sovereign
    b_D_F_ss = cal["b_D_F_ss"]           # F's holding of the D sovereign, in D units
    return dict(b_F_D_ss=b_F_D_ss, b_D_F_ss=b_D_F_ss,
                b_D_D_ss=cal["B_gov_D_ss"] - b_D_F_ss,
                b_F_F_ss=cal["B_gov_F_ss"] - b_F_D_ss / sz,
                b_D_F_ss_pc=b_D_F_ss / sz)   # the same holding per F capita


def _books(hold):
    # EACH BANK'S (OWN, FOREIGN) SOVEREIGN BOOK IN ITS OWN PER-CAPITA UNITS.
    return {"D": (hold["b_D_D_ss"], hold["b_F_D_ss"]),
            "F": (hold["b_F_F_ss"], hold["b_D_F_ss_pc"])}


def _country_blocks(cal, rk, p, mc, b_own, b_for, c):
    # CAPITAL -> FIRM -> WORKING-CAPITAL LOAN -> BANK FOR ONE COUNTRY AT (rk, p).
    # The firm block comes first: the loan the bank holds is zeta * the wage bill, and
    # w_ss is a firm object (N_ss = 1). The foreign bond earns its issuer's spread.
    Kap = capital_demand(rk, mc, cal, country=c)
    fm = steady_state_firm(cal, Kap, country=c)
    L_wc = cal[f"zeta_wc_{c}"] * fm["w_ss"]
    bk = steady_state_bank(cal, rk, Kap, b_own, b_for, p, country=c, L_wc_ss=L_wc,
                           spread_for=cal[f"credit_spread_target_{_OTHER[c]}"])
    return Kap, fm, bk


def _rescale_tfp(cal, rk, mc):
    # SET Z_ss SO THAT Y_ss = 1 (EXACT: THE CAPITAL FOC DOES NOT DEPEND ON Z).
    for c in ("D", "F"):
        a_c = cal[f"alpha_{c}"]
        cal[f"Z_ss_{c}"] = ((rk[c] + cal[f"delta_{c}"]) / (mc[c] * a_c)) ** a_c


def _size_books(cal, blk, p):
    # THE SOVEREIGN STOCKS THAT PUT EACH BANK'S OWN BOND AT ITS EXPOSURE INSTRUMENT, AT A TRIAL p.
    # D bank: own (1-phi_D)*Q_D*B_D, foreign p*Q_F*b_FD with b_FD = sz*phi_F*B_F. F bank: own
    # (1-phi_F)*Q_F*B_F, foreign Q_D*b_DF/(sz*p) with b_DF = phi_D*B_D. own/assets = e is linear
    # in (B_D, B_F), so the two exposure targets are one 2x2 solve.
    sz = size_ratio(cal)
    eD, eF = cal["sov_exposure_D"], cal["sov_exposure_F"]
    phD, phF = cal["foreign_share_D"], cal["foreign_share_F"]
    QD, QF = blk["D"][2]["Q_bdom_IC"], blk["F"][2]["Q_bdom_IC"]
    KL = {c: blk[c][0] + blk[c][2]["L_wc_ss"] for c in ("D", "F")}
    A = np.array([[(1.0 - phD) * QD * (1.0 - eD), -eD * p * QF * sz * phF],
                  [-eF * QD * phD / (sz * p), (1.0 - phF) * QF * (1.0 - eF)]])
    B_D, B_F = np.linalg.solve(A, np.array([eD * KL["D"], eF * KL["F"]]))
    cal["B_gov_D_ss"], cal["B_gov_F_ss"] = float(B_D), float(B_F)
    cal["b_D_F_ss"] = phD * cal["B_gov_D_ss"]          # F banks' slice of the D stock
    cal["b_F_D_ss"] = sz * phF * cal["B_gov_F_ss"]     # D banks' per-capita F book


def _net_foreign_income(cal, blk, p):
    # D's NET COUPON INCOME FROM ABROAD PER D CAPITA (D GOODS): EACH BOND YIELDS r + ITS ISSUER'S SPREAD.
    hold = _sovereign_holdings(cal)
    yD = cal["r_dep_D_target"] + cal["credit_spread_target_D"]
    yF = cal["r_dep_F_target"] + cal["credit_spread_target_F"]
    return (yF * p * blk["F"][2]["Q_bdom_IC"] * hold["b_F_D_ss"]
            - yD * blk["D"][2]["Q_bdom_IC"] * hold["b_D_F_ss"])


def _absorption(cal, blk, p, nfi_D):
    # CONSUMPTION FROM THE NATIONAL-INCOME IDENTITY P_CES*C = Y - I - G + NFI, PER COUNTRY.
    nfi = {"D": nfi_D, "F": -nfi_D / (size_ratio(cal) * p)}
    return {c: (blk[c][1]["Y_ss"] - blk[c][1]["I_ss"] - cal[f"G_{c}"] + nfi[c])
            / ces_price(p, cal, c) for c in ("D", "F")}, nfi


def _external_residual(p, cal, blk):
    # D's CURRENT ACCOUNT OVER Y AT A TRIAL p, WITH THE BOOKS SIZED AT THAT p.
    _size_books(cal, blk, p)
    nfi_D = _net_foreign_income(cal, blk, p)
    C, _ = _absorption(cal, blk, p, nfi_D)
    IM_D = import_demand(p, C["D"], ces_price(p, cal, "D"), cal, country="D")
    IM_F = import_demand(p, C["F"], ces_price(p, cal, "F"), cal, country="F")
    NX_D, _ = trade_balance(p, IM_D, IM_F, cal)
    return (NX_D + nfi_D) / blk["D"][1]["Y_ss"]


def _solve_external_balance(cal, blk, verbose):
    # p FROM D's CURRENT ACCOUNT BY BRACKETED ROOT (MONOTONE BUT NEARLY FLAT IN p); BOOKS LEFT AT IT.
    if verbose:
        print("=== Stage 1: external balance in p, books sized to the exposure targets ===")
    lo, hi = P_BRACKET
    f_lo, f_hi = _external_residual(lo, cal, blk), _external_residual(hi, cal, blk)
    if f_lo * f_hi > 0.0:
        raise RuntimeError(f"external balance has no root for p in [{lo}, {hi}] "
                           f"(residual {f_lo:+.3e}, {f_hi:+.3e})")
    p = brentq(_external_residual, lo, hi, args=(cal, blk), xtol=1e-15, rtol=1e-15,
               maxiter=500)
    res = _external_residual(p, cal, blk)
    assert abs(res) < 1e-13, f"external balance off by {res:.2e} at p = {p:.12f}"
    if verbose:
        print(f"  p_ss = {p:.10f}   current account / Y = {res:.2e}")
    return p


def _foreign_shares(cal, blk, p):
    # EACH BANK'S FOREIGN-BOND SHARE OF ITS ASSETS AT THE SIZED BOOKS.
    hold = _sovereign_holdings(cal)
    QD, QF = blk["D"][2]["Q_bdom_IC"], blk["F"][2]["Q_bdom_IC"]
    v_for = {"D": p * QF * hold["b_F_D_ss"], "F": QD * hold["b_D_F_ss_pc"] / p}
    v_own = {"D": QD * hold["b_D_D_ss"], "F": QF * hold["b_F_F_ss"]}
    return {c: v_for[c] / (blk[c][0] + blk[c][2]["L_wc_ss"] + v_own[c] + v_for[c])
            for c in ("D", "F")}


def _friction_references(cal):
    # CROSS-BORDER FRICTION ANCHORS AT WHICH THE CALIBRATED SS BOOKS ARE OPTIMAL.
    # At the SS the holder's FOC reads 1 + r + s_issuer = (1 + r + s_holder) * adj, with
    # adj = 1 + psi*(b - b_ref)/B_issuer: equal spreads give b_ref = b exactly; a spread gap moves
    # the friction's reference point, never the calibrated book.
    r = cal["r_dep_D_target"]
    sD, sF = cal["credit_spread_target_D"], cal["credit_spread_target_F"]
    adj_F = (1.0 + r + sD) / (1.0 + r + sF)      # the F bank's D-bond FOC
    adj_D = (1.0 + r + sF) / (1.0 + r + sD)      # the D bank's F-bond FOC
    cal["b_D_F_ref"] = cal["b_D_F_ss"] - cal["B_gov_D_ss"] * (adj_F - 1.0) / cal["psi_bD_F"]
    cal["b_F_D_ref"] = cal["b_F_D_ss"] - cal["B_gov_F_ss"] * (adj_D - 1.0) / cal["psi_bF_D"]
    # each bank's excess return on its foreign bond above its own required spread
    cal["excess_return_F_D_ss"] = sF - sD
    cal["excess_return_D_F_ss"] = sD - sF


def _check_targets(cal, bk, rk, fm, hold):
    # THE INSTRUMENTS MUST HOLD EXACTLY AT THE SOLVED STEADY STATE.
    sz = size_ratio(cal)
    own = {"D": bk["D"]["Q_bdom_IC"] * hold["b_D_D_ss"],
           "F": bk["F"]["Q_bdom_IC"] * hold["b_F_F_ss"]}
    for c in ("D", "F"):
        rdep_c = cal[f"r_dep_{c}_target"]
        # n_IC = n_ACCUM IS the bank identity: a gap means the franchise fixed point folded
        gap = bk[c]["n_ss_IC"] / bk[c]["n_ss_ACCUM"] - 1.0
        assert abs(gap) < 1e-10, (
            f"[{c}] n_IC/n_ACCUM - 1 = {gap:.2e}: the (f, spread, leverage) targets are "
            f"fold-blocked -- raise f until leverage is the least root.")
        assert abs(bk[c]["theta_ss"] - cal[f"leverage_target_{c}"]) < 1e-9, \
            f"[{c}] leverage {bk[c]['theta_ss']:.9f} != target {cal[f'leverage_target_{c}']}"
        assert abs((rk[c] - rdep_c) - cal[f"credit_spread_target_{c}"]) < 1e-12, \
            f"[{c}] spread {(rk[c] - rdep_c):.6f} != target {cal[f'credit_spread_target_{c}']}"
        exposure = own[c] / (bk[c]["theta_ss"] * bk[c]["n_ss"])
        assert abs(exposure - cal[f"sov_exposure_{c}"]) < 1e-12, \
            f"[{c}] own-sovereign exposure {exposure:.12f} != {cal[f'sov_exposure_{c}']}"
        assert abs(fm[c]["Y_ss"] - 1.0) < 1e-9, f"Y_ss_{c}={fm[c]['Y_ss']:.8f} != 1"
    assert abs(hold["b_D_F_ss"] / cal["B_gov_D_ss"] - cal["foreign_share_D"]) < 1e-14
    assert abs(hold["b_F_D_ss"] / sz / cal["B_gov_F_ss"] - cal["foreign_share_F"]) < 1e-14


def _deposit_market(cal, hh, beta, c, tol, verbose):
    # HOUSEHOLD SAVING MINUS BANK DEPOSIT SUPPLY (BASKET UNITS) FOR ONE COUNTRY AT A GUESSED beta.
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
    # deposits are held in units of the consumption basket, the bank's supply is in goods
    return dict(e=e, Pi=Pi, pi_e=pi_e, a_grid=make_asset_grid(cal, country=c),
                rdep=cal[f"r_dep_{c}_target"], Dep_supply=bk["Dep_supply_ss"] / P_CES,
                vN_ss=cal[f"chi_{c}"] / (1 + 1 / cal[f"frisch_{c}"]),  # GHH v(N_ss = 1)
                P_CES=P_CES,
                y_e=fm["w_ss"] / P_CES * e + (Div_ss - gs["Tax_ss"]) / P_CES)


def _solve_stage2(cal, hh, verbose):
    # {beta_D, beta_F}: EACH COUNTRY'S HOUSEHOLD DEPOSITS CLEAR ITS OWN BANK'S DEPOSIT SUPPLY.
    # The SS cross-border deposit position is zero: deposit-UIP's debt-elastic premium admits no
    # other rest point at one deposit rate, so the union market clears nation by nation.
    if verbose:
        print("\n=== Stage 2: deposit markets {beta_D, beta_F} ===")
    beta, sol = {}, {}
    for c in ("D", "F"):
        upper = 1 / (1 + hh[c]["rdep"]) - 1e-4   # keep rdep positive
        beta[c] = brentq(lambda b: _deposit_market(cal, hh, b, c, cal["tol_hh"], verbose)[0],
                         0.5, upper, xtol=1e-11)
        _, sol[c] = _deposit_market(cal, hh, beta[c], c, cal["tol_hh"], verbose)
    return beta["D"], beta["F"], sol["D"], sol["F"]


def _goods_market_check(cal, p_ss, C, P_CES, fm, nfi, verbose):
    # BOTH GOODS MARKETS AND D's CURRENT ACCOUNT AT THE TRUE STAGE-2 CONSUMPTION, ASSERTED.
    IM_D = import_demand(p_ss, C["D"], P_CES["D"], cal, "D")
    IM_F = import_demand(p_ss, C["F"], P_CES["F"], cal, "F")
    NX = dict(zip(("D", "F"), trade_balance(p_ss, IM_D, IM_F, cal)))
    gap = {c: (fm[c]["Y_ss"] - P_CES[c] * C[c] - fm[c]["I_ss"] - cal[f"G_{c}"] - NX[c])
           for c in ("D", "F")}
    ca = NX["D"] + nfi["D"]
    if verbose:
        print(f"  SS goods markets: D {gap['D']:.3e}, F {gap['F']:.3e};  "
              f"current account D {ca:.3e}  (NX_D = {NX['D']:+.6f}, NFI_D = {nfi['D']:+.6f})")
    assert max(abs(gap["D"]), abs(gap["F"]), abs(ca)) < 1e-8, (
        f"SS goods markets / current account do not clear: D {gap['D']:.2e}, F {gap['F']:.2e}, "
        f"CA {ca:.2e}")
    return NX, gap


def solve_steady_state(cal, verbose=True):
    # SOLVE THE TWO-COUNTRY STEADY STATE: BANKS, THE EXTERNAL BALANCE IN p, THEN STAGE 2.
    # deposit-UIP at a constant p and a zero deposit position needs one SS deposit rate
    assert cal["r_dep_D_target"] == cal["r_dep_F_target"], "the SS deposit rates must coincide"
    mc = {c: markup_ss(cal, c) for c in ("D", "F")}
    _calibrate_bank_frictions(cal, verbose)
    # the capital FOC at the SS: with one lambda every class earns rdep + lambda*mu/Omega
    rk = {c: cal[f"r_dep_{c}_target"] + cal[f"credit_spread_target_{c}"] for c in ("D", "F")}
    _rescale_tfp(cal, rk, mc)

    # stage 1: p and the books (the p-free blocks first; their bond prices need no books)
    blk = {c: _country_blocks(cal, rk[c], 1.0, mc[c], 0.0, 0.0, c) for c in ("D", "F")}
    p_ss = _solve_external_balance(cal, blk, verbose)
    _calibrate_bank_frictions(cal, verbose, share_for=_foreign_shares(cal, blk, p_ss))
    hold = _sovereign_holdings(cal)

    # every SS object at the solution, the government from the solved debt
    Kap, fm, bk = {}, {}, {}
    for c in ("D", "F"):
        Kap[c], fm[c], bk[c] = _country_blocks(cal, rk[c], p_ss, mc[c], *_books(hold)[c], c)
    _check_targets(cal, bk, rk, fm, hold)
    Q_bD_ss, Q_bF_ss = bk["D"]["Q_bdom_IC"], bk["F"]["Q_bdom_IC"]
    gs = {"D": govt_steady_state(cal, "D", Q_bD_ss), "F": govt_steady_state(cal, "F", Q_bF_ss)}
    P_CES = {c: ces_price(p_ss, cal, c) for c in ("D", "F")}
    for c in ("D", "F"):
        # labour FOC chi*N^(1/frisch) = w/P_CES at N_ss = 1
        cal[f"chi_{c}"] = fm[c]["w_ss"] / P_CES[c]
    cal["p_ss"] = p_ss
    _friction_references(cal)
    nfi_D = _net_foreign_income(cal, blk, p_ss)
    _, nfi = _absorption(cal, blk, p_ss, nfi_D)
    if verbose:
        print(f"\nStage 1 solution: rk_D={rk['D']:.6f}  rk_F={rk['F']:.6f}  p_ss={p_ss:.6f}")
        print(f"  Z_ss (rescaled): D={cal['Z_ss_D']:.6f}  F={cal['Z_ss_F']:.6f}")
        print(f"  B_gov: D={cal['B_gov_D_ss']:.6f}  F={cal['B_gov_F_ss']:.6f}   "
              f"n_ss_D={bk['D']['n_ss']:.4f}  n_ss_F={bk['F']['n_ss']:.4f}")

    # stage 2: household discount factors that clear the deposit markets
    hh = {c: _household_problem(cal, c, fm[c], bk[c], gs[c], p_ss, mc[c]) for c in ("D", "F")}
    beta_D_ss, beta_F_ss, (c_D_ss, D_D_ss, A_D_ss), (c_F_ss, D_F_ss, A_F_ss) = \
        _solve_stage2(cal, hh, verbose)
    C_D_ss = aggregate_consumption(D_D_ss, c_D_ss)
    C_F_ss = aggregate_consumption(D_F_ss, c_F_ss)
    if verbose:
        print(f"\nStage 2 solution: beta_D={beta_D_ss:.6f}  beta_F={beta_F_ss:.6f}")
        print(f"  A_D*P_D={A_D_ss * P_CES['D']:.4f}  Dep_supply_D={bk['D']['Dep_supply_ss']:.4f}")
        print(f"  A_F*P_F={A_F_ss * P_CES['F']:.4f}  Dep_supply_F={bk['F']['Dep_supply_ss']:.4f}")
    NX, gap = _goods_market_check(cal, p_ss, {"D": C_D_ss, "F": C_F_ss}, P_CES, fm, nfi, verbose)

    return dict(
        beta_D_ss=beta_D_ss, beta_F_ss=beta_F_ss,
        rk_D_ss=rk["D"], rk_F_ss=rk["F"], p_ss=p_ss,
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
        P_CES_D_ss=P_CES["D"], P_CES_F_ss=P_CES["F"],
        NX_D_ss=NX["D"], NX_F_ss=NX["F"], NFI_D_ss=nfi["D"], NFI_F_ss=nfi["F"],
        goods_gap_D_ss=gap["D"], goods_gap_F_ss=gap["F"],
    )
