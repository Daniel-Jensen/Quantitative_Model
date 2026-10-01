# SINGLE-POINT PERIOD MAP FOR THE RECURSIVE SOLUTION (THE IMAGE OF _inner_economy).
# Evaluates the full economy at ONE grid point (d, S) in QUANTITY FORM (realized
# payoffs = quantities x current prices).
#
# STATE (10): [K_D, K_F, P_D, P_F, b_DD, b_DF, b_FD, V_dep, s, Z_D] -- see state_grid.
#   The CB backstop needs NO state: an LTRO is a change in the COMPOSITION of the bank's
#   funding at an unchanged rate, so no stock is carried and no budget identity moves.
#   b_DD/b_DF are the two banks' CARRIED holdings of the D sovereign (B_D is their
#   sum). They were one state while the split was a fixed SS share; once both banks
#   bid through their own FOCs, last period's split is part of the state.
#   V_dep is the carried CROSS-BORDER deposit position, W_D - P_D. Each household's
#   claim used to be identified with its own bank's obligation -- true only under
#   NATIONAL clearing, where it is force-fed that bank's funding need. Under the union
#   deposit market they differ, and V is what they differ by: W_D = P_D + V,
#   W_F = P_F - V/p, so the union identity holds by construction.
#
# UNKNOWNS (13): [N_D, N_F, Kp_D, Kp_F, rdep_D, rdep_F, p, Q_bD, b_DF', Q_bF, b_FD',
#   A_D, A_F]
# RESIDUALS (13): 2 bank capital-Euler (occasionally binding), 2 labour, BOTH household
#   deposit Eulers, deposit-UIP, goods-D, BOTH banks' D-bond FOCs, union deposit
#   clearing, BOTH banks' F-bond FOCs. The D-sovereign is not force-fed to anyone: the
#   two bond FOCs are demand schedules, b_DD = B' - b_DF clears the market, and Q_bD is
#   the price that does it. The banker valuations alpha are still READ OFF the recursions
#   given a FROZEN continuation; under the collocation solve they are unknowns with their
#   own identity residual, so the freeze is exact at every node.
#
# OCCASIONALLY-BINDING IC -- Bocola Prop. 1 CLOSED FORM: mu = max{1 - E[Om]*R*n /
# (lambda*divertable assets), 0}, explicit and bounded in [0,1); the capital Euler
# E[Om(R_K-R)] = lambda_K*mu is the residual pinning K'. Same KKT as transition.py's
# Fischer-Burmeister, in the recursive-stable form (user-approved 2026-07-30;
# transition.py / bank.py FB left untouched).
#
# REGIMES ARE COMPOUND AND DATA-DRIVEN. decision_rules.regime_table maps the index the
# rules are stored under to a (default d', CB-active m') pair, and _regime_weights turns
# it into one probability per (s' node, regime') cell. Every expectation below is then
# the SAME weighted sum over that list, so adding a regime is a table entry rather than
# a new hand-written branch -- which is how the honoured/reneged TPI fork used to be
# written, with its own blended continuation and a separate _E3. The SDF in each regime
# is the genuine state-contingent GHH-composite kernel against that regime's own
# continuation. A regime with zero weight is never evaluated and aliases regime 0, so
# pi = 0 (no_default) and phi = 0 nest the smaller models EXACTLY, not to tolerance.
#
# Tier-3 cuts (documented, reversible): the F-sovereign split stays at SS shares (only
# the RISKY D bond is market-cleared); B_F at SS; the wc-wedge rate component uses
# current rdep.
import numpy as np

from blocks.firms import solve_firm_path, markup_ss
from blocks.capital import solve_capital_path
from blocks.trade import ces_price, import_demand, trade_balance, size_ratio
from solver_recursive.state_grid import default_prob

from solver_recursive.state_grid import (IK_D, IK_F, IP_D, IP_F, IBDD, IBDF,
                                         IBFD, IV, IS, IZ, NSTATE)

# the TEN pointwise Newton unknowns per regime -- single-sourced from decision_rules
# so the solve order, the writeback order and this unpacking cannot drift apart.
from solver_recursive.decision_rules import SOLVE as SOLVE7, DERIVED as DERIVED4


def gh_nodes(n=7):
    # GAUSS-HERMITE NODES/WEIGHTS FOR ONE STANDARD-NORMAL INNOVATION.
    # hermegauss is the PROBABILISTS' rule (weight exp(-x^2/2)), so the nodes are
    # already in standard-deviation units and s' = mean + sigma*node is exact.
    # (Bocola's GaussHermite.m returns the PHYSICISTS' rule and then uses the same
    # sigma*node map, which silently rescales his innovation by 1/sqrt(2).)
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return x, w / w.sum()


# Smoothing scale for the BACKSTOP guards (not for the net-worth floor, whose _smax
# keeps its own calibrated 1e-3). Away from the bound the bias is eps^2/(4*gap), so
# 1e-5 leaves the SS rest point at ~1e-11 -- three orders below the acceptance -- while
# still sitting ~3 decades above hybr's ~1.5e-8 forward-difference step, which is what
# the smoothing has to hide the kink from.
_GUARD_EPS = 1e-5


def _smax(x, floor, eps=1e-3):
    # SMOOTH MAX(x, floor) -- differentiable everywhere so the FD-Jacobian stays valid.
    return floor + 0.5 * ((x - floor) + np.sqrt((x - floor) ** 2 + eps ** 2))


def _smin(x, cap, eps=1e-3):
    # SMOOTH MIN(x, cap) -- the mirror of _smax, same differentiability argument.
    return cap - 0.5 * ((cap - x) + np.sqrt((cap - x) ** 2 + eps ** 2))


def _sclip(x, lo, hi, eps=None):
    # SMOOTH clip: the guards below must not be plateaus. A hard np.clip on a value
    # that is then FITTED puts a kink inside the box, and one saturated node moves the
    # interpolant at the ergodic centre by ~26% and can drive it negative (measured).
    return _smin(_smax(x, lo, eps or _GUARD_EPS), hi, eps or _GUARD_EPS)


def _regime_weights(wq, pd, phi, reg):
    # PROBABILITY OF EACH (s' NODE x REGIME') CELL, one weight vector per regime.
    # d' = 1 with the priced default probability pd, and the CB is active with
    # probability phi CONDITIONAL on d'. Reading phi off the rows that share a given d
    # is what makes one function correct for every table: where a d has a single row it
    # takes all of that d's mass (so phi cannot leak into a table that has no CB regime
    # for that d), and where it has two the mass splits phi / 1 - phi. The cells sum to
    # wq identically in every case, so the quadrature measure is preserved by
    # construction rather than by arithmetic that has to be re-checked per table.
    out = []
    for d_n, m_n in reg:
        p_d = pd if d_n else (1.0 - pd)
        n_rows = sum(1 for dd, _ in reg if dd == d_n)
        p_m = 1.0 if n_rows == 1 else (phi if m_n else 1.0 - phi)
        out.append(wq * p_d * p_m)
    return out


def _firm_capital(N, K, Kp, Z, cal, c):
    # ELEMENTWISE FIRM + CAPITAL BLOCK AT ONE POINT (untouched blocks).
    # Kp is floored to the Jermann-feasible band exactly as _cont_capital floors the
    # continuation. Under POINTWISE solving an infeasible Kp raised and the outer root
    # finder simply scored the guess badly; under the GLOBAL collocation solve
    # (collocation.py) an exception inside one finite-difference column kills the whole
    # Jacobian, so the period map has to be evaluable everywhere. This is Bocola's own
    # device -- his residual_model.m floors investment (inve = max(., 0.01)) and net
    # worth (N_tom = max(., 0.65)) for the same reason. The guard is slack by ~5% of K
    # at ksi = 0.5 and never binds near the fixed point.
    delta, ksi = cal[f"delta_{c}"], cal[f"ksi_{c}"]
    floor_ratio = (1.0 - delta) - delta * ksi / (1.0 - ksi)
    Kp = max(Kp, 1.0001 * floor_ratio * K)
    f = solve_firm_path(np.array([N]), np.array([K]), np.array([Z]), cal, c)
    cp = solve_capital_path(np.array([Kp]), K, 1.0, f["mpk"], cal, c,
                            Kap_lag_path=np.array([K]))
    return (f["Y"][0], f["w"][0], f["mpk"][0], cp["Q"][0], cp["rk"][0],
            cp["I"][0], cp["cap_profit"][0])


def _cont_capital(N_next, Kp, Kpp, Z, cal, c):
    # NEXT-PERIOD mpk' AND Q_K' FROM CONTINUATION RULES (for rk_next = t->t+1).
    # Clip the continuation next-capital to the Jermann-feasible band so a
    # transient-iterate Kp rule that dips too low does not raise inside the
    # untouched capital block (the guard is slack near the fixed point).
    delta, ksi = cal[f"delta_{c}"], cal[f"ksi_{c}"]
    floor_ratio = (1.0 - delta) - delta * ksi / (1.0 - ksi)
    Kpp = np.maximum(Kpp, 1.0001 * floor_ratio * Kp)
    f = solve_firm_path(N_next, np.full_like(N_next, Kp),
                        np.full_like(N_next, Z), cal, c)
    cp = solve_capital_path(Kpp, Kp, 1.0, f["mpk"], cal, c,
                            Kap_lag_path=np.full_like(Kpp, Kp))
    return f["mpk"], cp["Q"], f["Y"]


def point_residuals(S, d, x, cont, cal, ss, sproc, n_gh=7, no_default=False,
                    no_cb=False):
    # TEN UNIT-CONSISTENT RESIDUALS AT ONE POINT. x follows SOLVE; cont is the
    # FROZEN continuation RuleSet (previous iterate). alpha/Q_b/mu AND the
    # household aggregates C/A are computed here and returned in `out` to be
    # stored as the next-iterate rules.
    (N_D, N_F, Kp_D, Kp_F, rdep_D, rdep_F, p,
     Q_bD, b_DF_new, Q_bF, b_FD_new, A_D, A_F) = x
    K_D, K_F, P_D, P_F, b_D_D_lag, b_D_F_lag, b_F_D_lag, V_dep, s, Z_D = S
    # THE ACTIVATION PROBABILITY, OPTIONALLY STATE-CONTINGENT.
    # A CONSTANT phi offers the facility in EVERY state, including the ergodic centre
    # where the constraint barely binds -- and measured, that is where it does most of
    # its work: at phi = 0.5 the multiplier falls 37% at the rest point but only 8.9% at
    # the headline shock. A backstop that bites hardest in normal times is a permanent
    # liquidity subsidy, and it moves the steady state, which is exactly what makes the
    # cross-phi IRF comparison awkward (every activation rests somewhere different).
    # ltro_s_thr turns it into a real backstop: the offer probability is logistic in the
    # EXOGENOUS risk factor, so it is ~0 in normal times and ~phi_ltro in a crisis. It is
    # smooth in s, adds no state and no kink, and phi is a one-period-ahead probability
    # conditioned on today's s -- the same convention as the default probability pd.
    # ltro_s_thr = None keeps the constant-phi design exactly.
    phi = float(cal.get("phi_ltro", 0.0))
    _thr = cal.get("ltro_s_thr", None)
    if _thr is not None:
        phi *= 1.0 / (1.0 + np.exp(-(s - float(_thr)) / float(cal["ltro_s_width"])))
    B_D = b_D_D_lag + b_D_F_lag
    # THE REGIME INDEX IS COMPOUND: d indexes (default d', CB-active m') through the
    # table in decision_rules, because the TPI backstop is a second discrete regime the
    # continuation has to be conditioned on. At n_regimes = 2 the table is the identity
    # on the default indicator and every line below is the pre-TPI model unchanged.
    reg, nreg = cont.reg, cont.n_regimes
    d_reg, m_reg = reg[d]
    # no_cb switches the backstop OFF ENTIRELY -- zero probability in the expectation
    # AND zero facility in this regime -- exactly as no_default does for the default fork.
    # Both halves are needed: zeroing only phi would leave the m = 1 coefficient sets
    # solving a different economy from their twins while carrying no weight, which is not
    # the no-backstop model. It is what the first ladder stages and the nesting gate use.
    if no_cb:
        phi, m_reg = 0.0, 0
    # COUNTRY MASSES. Every variable here is PER CAPITA of its own country; sz =
    # size_F/size_D is the only place the asymmetry enters. Sovereign holdings are
    # carried in the ISSUER's per-capita units (so b_D_D + b_D_F = B_D still clears
    # the D market), which means the F bank's own book holds b_D_F/sz of the D bond
    # -- the same aggregate split across sz times as many agents. sz = 1 reproduces
    # the symmetric model bit for bit.
    sz = size_ratio(cal)

    surv = 1.0 - float(d_reg) * (1.0 - cal["recovery_rate_D"])

    # CARRIED HOUSEHOLD CLAIMS FROM THE CROSS-BORDER POSITION. V = W_D - P_D, so the
    # union identity W_D + p*W_F = P_D + p*P_F holds BY CONSTRUCTION rather than being
    # imposed on fitted states that each carry their own error (see state_grid).
    W_D = P_D + V_dep
    bkD, bkF = ss["ss_bank_D"], ss["ss_bank_F"]
    # GHH composite floors: 5% of SS consumption (x_ss/C_ss = 0.50, so ~10% of x_ss)
    _X_FLOOR_D, _X_FLOOR_F = 0.05 * ss["C_D_ss"], 0.05 * ss["C_F_ss"]
    _X_EPS = 1e-3 * ss["C_D_ss"]

    # --- firms + capital (current) ----------------------------------------
    # Z_D is the 7th state (deterministic TFP); F is never shocked.
    Y_D, w_D0, mpk_D, Q_D, _, I_D, capprof_D = _firm_capital(
        N_D, K_D, Kp_D, Z_D, cal, "D")
    Y_F, w_F0, mpk_F, Q_F, _, I_F, capprof_F = _firm_capital(
        N_F, K_F, Kp_F, cal["Z_ss_F"], cal, "F")
    P_CES_D = ces_price(np.array([p]), cal, "D")[0]
    P_CES_F = ces_price(np.array([p]), cal, "F")[0]

    # --- government (per-period body, quantity form) ----------------------
    # Bohn/Bocola fiscal rule on the SURVIVING stock, LINEAR in the debt level with
    # gamma_tau solved from a target debt root (government.govt_steady_state). Both
    # earlier readings of Bocola's gamma_tau = 1 (Table 1) failed here, for opposite
    # reasons. Reading it as a UNIT level coefficient made the default event a
    # fiscal windfall of 53% of GDP at recovery 0.45 -- taxes went to -53% of GDP against
    # a steady-state tax of 0.69%, and a 55% haircut left the sovereign MORE indebted
    # (post-haircut stock 1.09*B_ss). That killed 28-41 of 41 d=1 points.
    # The ELASTICITY form that replaced it fixed the default node but does not stabilise
    # DEBT: G_D = 0 makes Tax_ss the net interest bill alone (0.00297 against B_ss =
    # 0.98), so at phi = 1 dTax/dB = 0.003 and the debt root is 1.0002 -- the risk IRF
    # walked B_D into the box wall at q7 and the "trough then flat recovery" in the
    # 2026-08-25 figure WAS that wall. Raising the ELASTICITY is not the fix either:
    # phi = 15 swings taxes 0.29x-3.17x across the +-8% B band and to 6e-6 at the default
    # node, a convexity the mu=1 Chebyshev basis cannot carry, and the aliasing lands on
    # the centre node (measured: SS spread 129bp -> 250bp, C_D response flips positive).
    # LINEAR IN THE LEVEL with gamma_tau solved from a target debt root is what works:
    # it is exactly representable in the basis and its root is the calibrated object.
    # THE ANCHOR MUST BE REGIME-DEPENDENT. government.govt_transition already re-anchors
    # its branches to the post-haircut stock, "else the haircut becomes a tax-cut
    # windfall -> default expansionary, wrong-signed risk premium" -- this map used a
    # FIXED b_gov_ss in both regimes and never mirrored it. Latent under the old
    # elasticity rule (Tax_ss*0.45 is a 0.16%-of-Y windfall), it bit hard under the
    # linear rule, which multiplies gamma_tau by the default node's large deviation
    # (B*surv - B = -0.539): taxes went to -3.44% of Y against a steady-state +0.30%, a
    # transfer to households worth 12.6x the SS tax level ON IMPACT, against a bank
    # haircut loss that reaches them only through dividends (f_D = 0.04). The feared
    # state became one the household is BETTER OFF in, so pricing more of it RAISED
    # C_D (+0.19% on impact) and left hours at exactly +0.00%. Re-anchoring puts the
    # relief leg back to 0.16% of Y. d is a discrete regime index with its own rule set,
    # so the switch adds no kink to any Newton solve, and d = 0 is unchanged.
    if d_reg:
        anchor = cal["recovery_rate_D"] * ss["gs_D"]["b_gov_ss"]
        Tax_base = cal["G_D"] + cal["delta_b_D"] * anchor * (1.0 - ss["gs_D"]["Q_B_ss"])
    else:
        anchor, Tax_base = ss["gs_D"]["b_gov_ss"], ss["gs_D"]["Tax_ss"]
    Tax_D = Tax_base + ss["gs_D"]["gamma_tau"] * (B_D * surv - anchor)
    Tax_F = ss["gs_F"]["Tax_ss"]
    coupon_D = cal["delta_b_D"] * B_D * surv

    # --- end-of-period stocks (payoffs at current prices; alpha/mu-free) ---
    # b_D_D_lag / b_D_F_lag now come STRAIGHT FROM THE STATE -- no fixed share.
    # b_F_D_lag comes from the state; the F bank holds the remainder of a FIXED F stock
    # (B_F has no debt dynamics -- Tier-3 -- so only the SPLIT is a state).
    b_F_F_lag = cal["B_gov_F_ss"] - b_F_D_lag / sz
    db_D, db_F = cal["delta_b_D"], cal["delta_b_F"]

    eps, wq = gh_nodes(n_gh)
    pd = 0.0 if no_default else float(default_prob(s))
    # deterministic AR(1) for the TFP state (perfect foresight; no innovation)
    Z_next = (1.0 - sproc["rho_z"]) * sproc["z_star"] + sproc["rho_z"] * Z_D
    s_next = ((1.0 - sproc["rho_s"]) * sproc["s_star"] + sproc["rho_s"] * s
              + sproc["sigma_s"] * eps)
    m = eps.size

    # current-period bank valuations = FROZEN previous-iterate rules at this state
    # (breaks the S'<->Q_b knot; coincide at the fixed point)
    Sm = np.atleast_2d(S)
    # Q_bD IS THE SOLVED UNKNOWN, so the current D-bond price no longer comes off the
    # frozen previous iterate -- the Q_b/alpha knot is broken for the risky bond by
    # solving it, not by freezing it. Q_bF (safe, no market to clear here) keeps the
    # frozen-iterate convention.
    Q_bD_cur = Q_bD
    Q_bF_cur = Q_bF
    alpha_D_cur = float(cont.eval("alpha_D", d, Sm)[0])
    alpha_F_cur = float(cont.eval("alpha_F", d, Sm)[0])
    # WORKING-CAPITAL LOANS ARE BANK ASSETS (Bocola SV.C). Their SIZE needs the wage,
    # which needs r_wc, which needs mu, which needs the balance sheet -- a cycle that
    # runs through the (expensive) continuation. Break it exactly as the Q_b/alpha knot
    # above is broken: take r_wc from the FROZEN previous iterate for the loan QUANTITY.
    # The labour residual still uses the contemporaneous r_wc, and the two coincide at
    # the fixed point.
    r_wc_D_cur = float(cont.eval("r_wc_D", d, Sm)[0])
    r_wc_F_cur = float(cont.eval("r_wc_F", d, Sm)[0])
    zD_, zF_ = cal["zeta_wc_D"], cal["zeta_wc_F"]
    L_wc_D = zD_ * (w_D0 / (1.0 + zD_ * r_wc_D_cur)) * N_D
    L_wc_F = zF_ * (w_F0 / (1.0 + zF_ * r_wc_F_cur)) * N_F
    new_D = (cal["G_D"] + coupon_D - Tax_D) / Q_bD_cur
    Bp_D = (1.0 - cal["delta_b_D"]) * B_D * surv + new_D

    payD = surv * (db_D + (1.0 - db_D) * Q_bD_cur)
    payF = db_F + (1.0 - db_F) * Q_bF_cur
    X_D = ((mpk_D + (1.0 - cal["delta_D"]) * Q_D) * K_D
           + payD * b_D_D_lag + p * payF * b_F_D_lag)
    X_F = ((mpk_F + (1.0 - cal["delta_F"]) * Q_F) * K_F
           + payF * b_F_F_lag + payD * (b_D_F_lag / sz) / p)
    ng_D, ng_F = X_D - P_D, X_F - P_F
    # MARKET CLEARING BY CONSTRUCTION: whatever the F bank does not take, the D bank
    # holds. b_DF_new is a Newton unknown pinned by the F bank's own D-bond FOC below,
    # and Q_bD is the unknown pinned by the D bank's. Floored so a transient iterate
    # cannot hand the D bank a negative book (which would flip the sign of lev_D).
    b_D_F_new = _smax(b_DF_new, 1e-4, 1e-5)
    b_D_D_new = _smax(Bp_D - b_D_F_new, 1e-4, 1e-5)
    # SAME FOR THE F SOVEREIGN: b_FD is the D bank's chosen holding, the F bank takes
    # the rest of the fixed F stock. Floored identically.
    b_F_D_new = _smax(b_FD_new, 1e-4, 1e-5)
    b_F_F_new = _smax(cal["B_gov_F_ss"] - b_F_D_new / sz, 1e-4, 1e-5)
    # end-of-period portfolio valued at current PRICES (Q_b), not payoffs
    assets_D = (Q_D * Kp_D + Q_bD_cur * b_D_D_new + p * Q_bF_cur * b_F_D_new
                + L_wc_D)
    assets_F = (Q_F * Kp_F + Q_bF_cur * b_F_F_new
                + Q_bD_cur * (b_D_F_new / sz) / p + L_wc_F)
    n_D = (1.0 - cal["f_D"]) * ng_D + cal["omega_ent_D"] * assets_D
    n_F = (1.0 - cal["f_F"]) * ng_F + cal["omega_ent_F"] * assets_F
    # BOCOLA FEASIBILITY FLOOR (his N_tom = max(.,0.65)): a deep post-haircut default
    # corner drives net worth negative, where mu saturates the crude cap and the pointwise
    # solve fails, poisoning the global fit. Smooth-floor net worth at a small fraction of
    # n_ss so it is INACTIVE in the ergodic region (nw_floor_frac=0 => baseline unchanged)
    # and only catches the deep default corners.
    nwf = cal.get("nw_floor_frac", 0.0)
    if nwf > 0.0:
        n_D = _smax(n_D, nwf * bkD["n_ss"])
        n_F = _smax(n_F, nwf * bkF["n_ss"])
    dep_D, dep_F = assets_D - n_D, assets_F - n_F
    # deposits fund the whole book INCLUDING the working-capital loan; the loan is
    # repaid at the lending rate, so the obligation carried forward is net of it --
    # Bocola: P' = R*(Q*K' + q*B' + L - N') - R_W*L
    Pp_D = (1.0 + rdep_D) * dep_D - (1.0 + r_wc_D_cur) * L_wc_D
    Pp_F = (1.0 + rdep_F) * dep_F - (1.0 + r_wc_F_cur) * L_wc_F
    # UNION DEPOSIT MARKET. A_D is the D household's CHOSEN real saving (Newton unknown,
    # pinned by euler_D); A_F follows from ONE union-wide clearing in D-good units --
    # total household saving must fund the whole union book. Under the old NATIONAL
    # clearing each household held exactly its own bank's dep_X, which is what made
    # consumption a bookkeeping residual of the balance sheet and left euler_F violated.
    # The gap A_D*P_CES_D - dep_D is the cross-border deposit position: the absorption
    # margin that breaks the national S = I trap.
    dep_union = dep_D + sz * p * dep_F
    save_union = A_D * P_CES_D + sz * p * A_F * P_CES_F
    nfa_dep_D = A_D * P_CES_D - dep_D
    # THE HOUSEHOLD'S CARRIED CLAIM MIRRORS THE BANK'S OBLIGATION, WC NETTING INCLUDED.
    # P_state is NET of the working-capital receivable (bank.py), and in the pre-union
    # model the household's carried wealth WAS that net object. Grossing A_D without the
    # same deduction makes W_D drift away from P_D by (1+r_wc)*L_wc EVERY period, so the
    # union identity W_D + p*W_F = P_D + p*P_F fails out of the SS -- measured as the
    # rest point walking p to 0.932 and mu to 0 with every point still clearing at 1e-14.
    # The deduction is charged to the household whose bank carries the loan (exact at
    # nfa = 0, first-order in the cross-border position otherwise).
    Wp_D = (1.0 + rdep_D) * A_D * P_CES_D - (1.0 + r_wc_D_cur) * L_wc_D
    # THE CARRIED CROSS-BORDER POSITION, the only wealth state. V' = W_D' - P_D' and both
    # legs carry the SAME working-capital deduction, so it cancels exactly and V' is just
    # this period's cross-border flow grossed up at the D deposit rate. That exactness is
    # the point: the union identity can no longer drift, which is what was poisoning
    # euler_F when W_F was derived from it.
    Vp_dep = (1.0 + rdep_D) * nfa_dep_D

    # --- per-regime continuation objects (FROZEN cont) --------------------
    def _cont(j):
        # CONTINUATION RULE VALUES + NEXT-PERIOD CAPITAL RETURN, REGIME j.
        Sn = np.empty((m, NSTATE))
        Sn[:, IK_D] = Kp_D; Sn[:, IK_F] = Kp_F
        Sn[:, IP_D] = Pp_D; Sn[:, IP_F] = Pp_F
        Sn[:, IBDD] = b_D_D_new; Sn[:, IBDF] = b_D_F_new
        Sn[:, IBFD] = b_F_D_new
        Sn[:, IV] = Vp_dep
        Sn[:, IS] = s_next
        Sn[:, IZ] = Z_next
        r = cont.eval_all(j, cont.grid.clip(Sn))
        # guard continuation outputs before they enter fractional powers / sqrt
        # (a deep default-regime iterate can push N, p, C negative -> NaN)
        for k in ("N_D", "N_F"):
            r[k] = np.maximum(r[k], 0.05)
        for k in ("C_D", "C_F"):
            r[k] = np.maximum(r[k], 1e-3)
        r["p"] = np.maximum(r["p"], 1e-2)
        mpkD_n, QKD_n, _ = _cont_capital(r["N_D"], Kp_D, r["Kp_D"],
                                         Z_next, cal, "D")
        mpkF_n, QKF_n, _ = _cont_capital(r["N_F"], Kp_F, r["Kp_F"],
                                         cal["Z_ss_F"], cal, "F")
        rkD_n = (mpkD_n + (1.0 - cal["delta_D"]) * QKD_n) / Q_D - 1.0
        rkF_n = (mpkF_n + (1.0 - cal["delta_F"]) * QKF_n) / Q_F - 1.0
        return r, rkD_n, rkF_n

    # QUADRATURE WEIGHTS OVER (s' node x regime'), and the continuation in each regime.
    # A regime carrying ZERO weight is never evaluated and is ALIASED to regime 0: that
    # is what makes pi = 0 (and, with the CB, phi = 0) nest the smaller model EXACTLY
    # rather than to solver tolerance, and it is also where the saved continuation
    # evaluation comes from -- each one costs a full grid interpolation plus two
    # capital blocks.
    wgt = _regime_weights(wq, pd, phi, reg)
    R, RKD, RKF = [], [], []
    for j in range(nreg):
        if j > 0 and not np.any(wgt[j] > 0.0):
            R.append(R[0]); RKD.append(RKD[0]); RKF.append(RKF[0])
            continue
        rj, rkDj, rkFj = _cont(j)
        R.append(rj); RKD.append(rkDj); RKF.append(rkFj)

    # BRANCH SDFs -- the GENUINE state-contingent household kernel, in EVERY regime.
    # Bocola's banker discounts with Lambda' = beta*c/c' (his exp_1/exp_2 in
    # residual_model.m), which is what makes his constraint TIGHTEN with risk: when the
    # sovereign shock hits, next-period consumption moves, c/c' falls, Omega falls and
    # mu = 1 - E[Om]*R*n/lev RISES. The old convention used a CONSTANT beta_inter on the
    # no-default branch (documented in CLAUDE.md as "Lambda^nd = beta_inter"), which
    # cannot fall -- so nothing offset the rising alpha', E[Om] climbed, and the
    # constraint went SLACK exactly when it should bind. Measured at three calibrations
    # including Bocola's own (mu_ss 0.001, alpha_ss 1.026): mu -> 0 as soon as p^d rises.
    # The kernel is the GHH composite x = C - chi*N^(1+1/nu)/(1+1/nu), the same object
    # the deposit Euler below uses, evaluated on the frozen current rules (the Q_b/alpha
    # knot-breaking convention) against each branch's continuation.
    biD, biF = cal["beta_inter_D"], cal["beta_inter_F"]
    frD, frF = cal["frisch_D"], cal["frisch_F"]
    sgD, sgF = cal["sigma_D"], cal["sigma_F"]

    def _ghh(C, N, chi, fr):
        # GHH CONSUMPTION COMPOSITE x = C - v(N), SMOOTH-floored. The old hard
        # np.maximum(., 1e-9) is a plateau with ZERO gradient, and x_ss = 0.391 is eight
        # orders above it: a corner iterate that pushed x down landed on the plateau,
        # the FD-Jacobian saw nothing, and the point stayed stuck at |F| ~ 4e8 forever
        # (5-10 of 19 points, measured). Flooring at a fraction of C_ss with the standard
        # smoothing keeps a usable derivative and is inactive anywhere near the ergodic
        # set, where x/C = 0.50.
        return _smax(C - chi * N ** (1 + 1 / fr) / (1 + 1 / fr),
                     _X_FLOOR_D, _X_EPS)

    x_cur_D = _ghh(float(cont.eval("C_D", d, Sm)[0]), N_D, cal["chi_D"], frD)
    x_cur_F = _ghh(float(cont.eval("C_F", d, Sm)[0]), N_F, cal["chi_F"], frF)
    XN_D = [_ghh(R[j]["C_D"], R[j]["N_D"], cal["chi_D"], frD) for j in range(nreg)]
    XN_F = [_ghh(R[j]["C_F"], R[j]["N_F"], cal["chi_F"], frF) for j in range(nreg)]
    Lam_D = [biD * (x_cur_D / XN_D[j]) ** sgD for j in range(nreg)]
    Lam_F = [biF * (x_cur_F / XN_F[j]) ** sgF for j in range(nreg)]
    Om_D = [Lam_D[j] * (cal["f_D"] + (1 - cal["f_D"]) * R[j]["alpha_D"])
            for j in range(nreg)]
    Om_F = [Lam_F[j] * (cal["f_F"] + (1 - cal["f_F"]) * R[j]["alpha_F"])
            for j in range(nreg)]

    def _E(v):
        # EXPECTATION OVER (s' NODE x REGIME'): v carries ONE ARRAY PER REGIME.
        return float(sum(np.dot(wgt[j], v[j]) for j in range(nreg)))

    # --- banker FOC block (Bocola Prop. 1 CLOSED-FORM multiplier) -----------
    lKD, lKF = cal["lambda_K_D"], cal["lambda_K_F"]
    lbDD, lbFF = cal["lambda_bD_D"], cal["lambda_bF_F"]
    E_Om_D = _E(Om_D)
    E_Om_F = _E(Om_F)
    # divertable assets (same leverage term as the FB slack), frozen-price valued
    lev_D = max(lKD * Q_D * Kp_D + lbDD * Q_bD_cur * b_D_D_new
                + cal["lambda_bF_D"] * p * Q_bF_cur * b_F_D_new
                + lKD * L_wc_D, 1e-6)
    lev_F = max(lKF * Q_F * Kp_F + lbFF * Q_bF_cur * b_F_F_new
                + cal["lambda_bD_F"] * Q_bD_cur * (b_D_F_new / sz) / p
                + lKF * L_wc_F, 1e-6)
    # THE LTRO BACKSTOP -- BOCOLA'S OWN, residual_model_ltro_firstperiod.m:
    #     mu_ratio = N'/(lambda*A')   ->   (N' + m)/(lambda*(A' - m))
    # With per-period probability phi (cal["phi_ltro"], a per-experiment scalar and NOT a
    # state) the central bank offers collateralised credit of size m. The facility is
    # lent at the DEPOSIT RATE, so it is a change in the COMPOSITION of the bank's
    # funding -- divertable deposits for non-divertable central-bank credit -- at an
    # unchanged rate. Every budget identity in the model is therefore untouched: with
    # r_ltro = rdep the deposit obligation P' = (1+rdep)(dep + m) - (1+r_wc)L_wc is
    # exactly (1+rdep)(assets - n) - (1+r_wc)L_wc, the household swaps one claim for
    # another at the same rate so union clearing and nfa are unchanged, and the CB lends
    # at the rate it pays so its carry is identically zero and NO remittance is needed.
    # The whole effect passes through the incentive constraint, which is why this is a
    # two-line change and why the Walras leak that dogged a purchase design cannot arise.
    #
    # It does TWO things where a bond purchase does one: the assets it funds leave the
    # divertable base AND the funding counts as equity. To first order that is
    # (1 + A/n) = 1 + leverage = 6x the constraint relief of a purchase of the same size,
    # and the numerator term is a margin no quantity of bond-buying can reach.
    #
    # SIZE MATTERS AND BOCOLA'S OWN IS A TRAP. 2.0% of quarterly GDP already unbinds the
    # constraint at the SS and 3.4% unbinds it in the crisis state, against his 40%. At
    # his size mu = 0 with huge margin in every relieved regime, so the whole m = 1
    # coefficient set sits ON the KKT kink -- exactly where a Chebyshev interpolant of a
    # C0 multiplier is least reliable. Size the envelope to RELIEVE, not to unbind.
    m_ltro_D = cal.get("ltro_D", 0.0) if m_reg else 0.0
    m_ltro_F = cal.get("ltro_F", 0.0) if m_reg else 0.0
    n_IC_D, n_IC_F = n_D + m_ltro_D, n_F + m_ltro_F
    # the pledged collateral is the bank's own sovereign, so it leaves the base at that
    # asset's lambda; under the single-lambda doctrine every lambda is equal and this is
    # Bocola's lambda*(A - m) exactly, but writing it per-asset keeps it right if they
    # ever diverge
    lev_IC_D = max(lev_D - lbDD * m_ltro_D, 1e-6)
    lev_IC_F = max(lev_F - lbFF * m_ltro_F, 1e-6)
    # closed-form mu in [0, mu_cap]: the lower max is Bocola's KKT switch (hard, as in
    # his residual_model.m -- it IS the complementarity, not a numerical guard); the
    # upper cap keeps alpha = E_Om R/(1-mu) finite when a default-regime net worth n
    # goes negative, and is SMOOTHED because unlike the KKT switch it is a guard that
    # would otherwise plant a plateau in a fitted object.
    _MU_CAP = 0.95
    mu_D = float(_smin(max(1.0 - E_Om_D * (1.0 + rdep_D) * n_IC_D / lev_IC_D, 0.0),
                       _MU_CAP, _GUARD_EPS))
    mu_F = float(_smin(max(1.0 - E_Om_F * (1.0 + rdep_F) * n_IC_F / lev_IC_F, 0.0),
                       _MU_CAP, _GUARD_EPS))
    # capital-Euler surplus E[Om(R_K - R)] - lambda_K*mu (the residual pinning K').
    # Normalise by the O(1) discount kernel E_Om so the residual is in return units
    # (dividing by lambda_K*mu_ss amplifies ~130x and wrecks the solver step).
    E_Om_rk_D = _E([Om_D[j] * (RKD[j] - rdep_D) for j in range(nreg)])
    E_Om_rk_F = _E([Om_F[j] * (RKF[j] - rdep_F) for j in range(nreg)])
    cap_eul_D = (E_Om_rk_D - lKD * mu_D) / E_Om_D
    cap_eul_F = (E_Om_rk_F - lKF * mu_F) / E_Om_F
    # alpha follows from mu, so _MU_CAP already bounds it at E_Om*R/(1-_MU_CAP) ~ 38.
    # The OLD hard clip at 8.0 bound FIRST (at mu ~ 0.76) and bit exactly where
    # nw_floor_frac puts the deep default corner: n = 0.15*n_ss gives mu ~ 0.85 and
    # alpha ~ 13, and Bocola's own solution reaches alpha/alpha_ss = 13.6. _ALPHA_CAP
    # is now a smooth backstop above that, so the two guards no longer contradict.
    # NB the alpha RECURSION alpha = E[Om]*R/(1-mu) with Om = beta*[f + (1-f)*alpha'] has
    # slope beta*(1-f)*R/(1-mu), which exceeds 1 once mu > 1 - beta*(1-f)*R ~ 0.04. Above
    # that the cap is not cosmetic -- it is what arrests a divergent fixed point. Keep it
    # tunable so the trade-off (truncating a region Bocola reaches vs letting the deep
    # default corners run away) can be measured rather than assumed.
    _ALPHA_CAP = float(cal.get("alpha_cap", 40.0))
    alpha_D_new = float(_sclip(E_Om_D * (1.0 + rdep_D) / (1.0 - mu_D), 0.05, _ALPHA_CAP))
    alpha_F_new = float(_sclip(E_Om_F * (1.0 + rdep_F) / (1.0 - mu_F), 0.05, _ALPHA_CAP))

    # bond Euler (D risky: haircut wherever the regime defaults; F safe), same mu.
    # THE D-BOND PAYOFF PER REGIME: the gross HM perpetuity payoff at that regime's own
    # continuation price, times the haircut in the regimes that default.
    hc = [cal["recovery_rate_D"] if d_n else 1.0 for d_n, _ in reg]
    payD_gross = [db_D + (1.0 - db_D) * R[j]["Q_bD"] for j in range(nreg)]
    payD = [hc[j] * payD_gross[j] for j in range(nreg)]
    E_Om_payD = _E([Om_D[j] * payD[j] for j in range(nreg)])
    # BOND-PRICE DECOMPOSITION LEGS (diagnostics only -- no residual uses them).
    # The D bank's FOC is E[Om*pay] = Q*(E[Om]*R + lambda_bD*mu), so
    #   Q = E[pay]/R                                  risk-free discounting of the
    #                                                 EXPECTED payoff (actuarial), and
    #   Q = that  x  [E[Om*pay]/(E[Om]*E[pay])]       the RISK premium (the Om-payoff
    #                                                 covariance: Om is high exactly
    #                                                 when the D bond pays little), and
    #   Q = that  x  [E[Om]R/(E[Om]R + lambda*mu)]    the LIQUIDITY premium (the
    #                                                 constraint's own wedge).
    # Splitting E[pay] further at the NO-DEFAULT payoff isolates the EXPECTED LOSS.
    # In logs the four legs are exactly additive to log Q_bD, which is what makes the
    # figure a decomposition rather than an attribution. The naming follows Bocola's
    # Table 4 ("risk premia" vs the multiplier's "liquidity" component).
    E_payD = _E(payD)                                       # physical-measure payoff
    # the same with DEFAULT SWITCHED OFF: each no-default regime keeps its own price,
    # and a defaulting regime is replaced by the plain no-default one (regime 0), so
    # the leg isolates the haircut AND the default state's own repricing together
    E_payD_nodef = _E([payD_gross[j] if not reg[j][0] else payD_gross[0]
                       for j in range(nreg)])
    payF = [db_F + (1.0 - db_F) * R[j]["Q_bF"] for j in range(nreg)]
    E_Om_payF = _E([Om_F[j] * payF[j] for j in range(nreg)])
    # THE D BANK'S F-BOND FOC -- DIAGNOSTIC ONLY, and that is the point. Under a single
    # lambda the banker's portfolio problem yields E[Om(R_j - R)] = lambda_j*mu for EVERY
    # asset j. Capital and the D bond both have their FOC in the residual system; the D
    # bank's F-bond position does NOT -- b_F_D is pinned at its SS value (see the holdings
    # block above), so nothing holds this condition to zero. Its size is the measure of
    # how far the fixed F-leg is from an optimising portfolio.
    # THE F BANK'S F-BOND FOC -> Q_bF (home leg), and THE D BANK'S -> b_FD (foreign leg,
    # carrying the cross-border portfolio adjustment cost psi_bF_D exactly as the D-bond
    # market carries psi_bD_F). Both are RESIDUALS now. Previously only the F bank's was
    # used, and only to READ OFF Q_bF at a quantity nobody chose: b_F_D was pinned at its
    # SS value forever, so the D bank held F bonds with no first-order condition behind
    # the position. That is an incomplete optimum -- under a single lambda the banker has
    # an FOC for EVERY asset -- and it is also why Q_bF could not move: with both F
    # quantities frozen there was no demand schedule to shift, hence no flight to safety
    # and no contagion through F repricing.
    E_Om_payF_D = _E([Om_D[j] * payF[j] for j in range(nreg)])
    dmd_F_home = E_Om_F * (1.0 + rdep_F) + lbFF * mu_F
    dmd_F_for = E_Om_D * (1.0 + rdep_D) + cal["lambda_bF_D"] * mu_D
    adj_D = 1.0 + cal["psi_bF_D"] * (b_F_D_new - cal["b_F_D_ss"]) / cal["B_gov_F_ss"]
    # smooth bounds (Bocola's own q bottoms at 0.639, so these never bind in the
    # ergodic region; they are numerical backstops, not economics)
    # THE D BANK'S D-BOND FOC, kept as an implicit residual (bondD below) instead of
    # being solved for the price. Same equation as before -- it is the SOLUTION METHOD
    # that changes: the price is now the market-clearing unknown rather than a value
    # read off a recursion at a quantity nobody chose.
    dmd_D = E_Om_D * (1.0 + rdep_D) + lbDD * mu_D
    # THE F BANK'S D-BOND FOC, with Bocola's cross-border portfolio adjustment cost on
    # the position (psi_bD_F, already calibrated at 0.05 and previously used only by the
    # deleted transition solver). Without it the two banks' schedules for the SAME bond
    # differ only through mu_X, and lambda_bD*Q*b is ~7% of divertable assets, so both
    # demand curves are near-flat: the split is numerically indeterminate and the Newton
    # is ill-conditioned. The cost is what gives the foreign schedule a definite slope,
    # and it is zero at the SS position so it does not distort the steady state.
    E_Om_payD_F = _E([Om_F[j] * payD[j] for j in range(nreg)])
    dmd_F = E_Om_F * (1.0 + rdep_F) + cal["lambda_bD_F"] * mu_F
    adj_F = 1.0 + cal["psi_bD_F"] * (b_D_F_new - cal["b_D_F_ss"]) / cal["B_gov_D_ss"]
    # --- working-capital wedge + dividends (untouched formulas) ------------
    r_wc_D = rdep_D + lKD * mu_D / E_Om_D
    r_wc_F = rdep_F + lKF * mu_F / E_Om_F
    zD, zF = cal["zeta_wc_D"], cal["zeta_wc_F"]
    w_D = w_D0 / (1.0 + zD * r_wc_D)
    w_F = w_F0 / (1.0 + zF * r_wc_F)
    wc_inc_D = zD * r_wc_D * w_D * N_D
    wc_inc_F = zF * r_wc_F * w_F * N_F
    mcD, mcF = markup_ss(cal, "D"), markup_ss(cal, "F")
    div_bank_D = cal["f_D"] * ng_D - cal["omega_ent_D"] * assets_D
    div_bank_F = cal["f_F"] * ng_F - cal["omega_ent_F"] * assets_F
    # WHERE THE WORKING-CAPITAL FINANCING INCOME GOES. cal["wc_rebate"] = 1.0 (default,
    # unchanged) hands it to households as dividends, which makes the spread a pure
    # INTRA-PERIOD TRANSFER: firms pay zeta*r_wc*w*N, households receive the same amount,
    # and under GHH there is no wealth effect to offset it, so hours and output RISE with
    # the spread. Measured: at f=0.12 the whole +0.499% output response is this rebate --
    # switching the wedge off entirely (zeta_wc=0) collapses it to -0.001%. wc_rebate = 0
    # instead treats the wedge as a real resource cost, keeping the wage channel while
    # removing the offsetting income. (The third option -- routing it to BANK net worth,
    # which is where it economically belongs -- changes the closed-form leverage/spread
    # calibration and is not wired here.)
    reb = float(cal.get("wc_rebate", 0.0))
    Div_D = (1 - mcD) * Y_D + capprof_D + div_bank_D + reb * wc_inc_D
    Div_F = (1 - mcF) * Y_F + capprof_F + div_bank_F + reb * wc_inc_F

    # --- household aggregate (BUDGET + DEPOSIT EULER, rep-agent GHH) ---------
    # Deposits clear BY QUANTITY: household composite deposits A = bank funding
    # dep / P_CES. Carried wealth = P/P_CES (the gross deposit claim = the bank's
    # gross obligation state). Consumption is the residual of the budget; the
    # deposit RATE is pinned by the household Euler (residual 5). beta_eff =
    # 1/(1+rdep_ss) makes the rep-agent Euler reproduce the HA aggregate at the SS;
    # hh_T anchors the SS budget so C = C_ss, A = A_ss exactly.
    frisch_D, frisch_F = cal["frisch_D"], cal["frisch_F"]
    sigD, sigF = cal["sigma_D"], cal["sigma_F"]
    inc_D = (w_D / P_CES_D) * N_D + (Div_D - Tax_D) / P_CES_D + ss.get("hh_T_D", 0.0)
    inc_F = (w_F / P_CES_F) * N_F + (Div_F - Tax_F) / P_CES_F + ss.get("hh_T_F", 0.0)
    # smooth bounds; eps scales with C_ss so the smoothing is relatively as tight as
    # the O(1) guards above (Bocola's c stays within +-12%, so these never bind)
    # eps at 1e-3*C_ss, not _GUARD_EPS: at 1e-5 the bound was effectively HARD and a
    # corner iterate parked on it with no gradient (C_F sat at exactly 2*C_ss).
    _ce = 1e-3 * ss["C_D_ss"]
    # F's carried claim is its own bank's obligation LESS the cross-border position, in
    # F-good units. This used to be the union residual (P_D + p*P_F - W_D)/p, which made
    # W_F inherit the fit error of three large states at once; C_F is a ~0.79 difference
    # of ~8-sized terms, so a 0.06% error in W_F became a 1.1% error in euler_F -- 22x
    # the D-side error, corr 0.97. V is small, so both sides now carry one large state's
    # error plus the same small one.
    W_F = P_F - V_dep / (sz * p)
    C_D = float(_sclip(W_D / P_CES_D + inc_D - A_D,
                       0.15 * ss["C_D_ss"], 3.0 * ss["C_D_ss"], _ce))
    C_F = float(_sclip(W_F / P_CES_F + inc_F - A_F,
                       0.15 * ss["C_F_ss"], 3.0 * ss["C_F_ss"], _ce))

    # deposit Euler (GHH composite x = C - vN): (x)^-sigma = beta_eff*E[(1+r')x'^-sigma]
    vN_D = cal["chi_D"] * N_D ** (1 + 1 / frisch_D) / (1 + 1 / frisch_D)
    vN_F = cal["chi_F"] * N_F ** (1 + 1 / frisch_F) / (1 + 1 / frisch_F)
    vNp_D = [cal["chi_D"] * R[j]["N_D"] ** (1 + 1 / frisch_D) / (1 + 1 / frisch_D)
             for j in range(nreg)]
    vNp_F = [cal["chi_F"] * R[j]["N_F"] ** (1 + 1 / frisch_F) / (1 + 1 / frisch_F)
             for j in range(nreg)]
    xp_D = [_smax(R[j]["C_D"] - vNp_D[j], _X_FLOOR_D, _X_EPS) for j in range(nreg)]
    xp_F = [_smax(R[j]["C_F"] - vNp_F[j], _X_FLOOR_F, _X_EPS) for j in range(nreg)]
    rp_D = [(1.0 + rdep_D) * P_CES_D / ces_price(R[j]["p"], cal, "D") - 1.0
            for j in range(nreg)]
    rp_F = [(1.0 + rdep_F) * P_CES_F / ces_price(R[j]["p"], cal, "F") - 1.0
            for j in range(nreg)]
    beff_D = 1.0 / (1.0 + cal["r_dep_D_target"])
    beff_F = 1.0 / (1.0 + cal["r_dep_F_target"])
    E_mu_D = _E([(1.0 + rp_D[j]) * xp_D[j] ** (-sigD) for j in range(nreg)])
    E_mu_F = _E([(1.0 + rp_F[j]) * xp_F[j] ** (-sigF) for j in range(nreg)])
    xC_D = float(_smax(C_D - vN_D, _X_FLOOR_D, _X_EPS))
    xC_F = float(_smax(C_F - vN_F, _X_FLOOR_F, _X_EPS))
    euler_D = xC_D ** (-sigD) / (beff_D * E_mu_D) - 1.0
    euler_F = xC_F ** (-sigF) / (beff_F * E_mu_F) - 1.0

    IM_D = import_demand(np.array([p]), np.array([C_D]), np.array([P_CES_D]), cal, "D")[0]
    IM_F = import_demand(np.array([p]), np.array([C_F]), np.array([P_CES_F]), cal, "F")[0]
    NX_D, NX_F = (lambda a, b: (a[0], b[0]))(
        *trade_balance(np.array([p]), np.array([IM_D]), np.array([IM_F]), cal))

    # E[p'] for deposit-UIP: the SAME measure over regimes as every other expectation
    # in this block (it used to average the no-default branch alone, which priced the
    # cross-border deposit leg as if default never happened)
    Ep_next = _E([R[j]["p"] for j in range(nreg)])

    # --- the SEVEN residuals (transition.py's per-period image) -------------
    # slack is read off the SAME constraint the multiplier came from, so the
    # complementarity mu*slack = 0 still holds point-wise once the facility is on
    slack_D = alpha_D_cur * n_IC_D - lev_IC_D
    slack_F = alpha_F_cur * n_IC_F - lev_IC_F
    res = np.array([
        cap_eul_D,                                               # 1 cap Euler D -> Kp_D
        cap_eul_F,                                               # 2 cap Euler F -> Kp_F
        (cal["chi_D"] * N_D ** (1 / frisch_D) - w_D / P_CES_D)    # 3 lab_D -> N_D
        / (w_D / P_CES_D),
        (cal["chi_F"] * N_F ** (1 / frisch_F) - w_F / P_CES_F)    # 4 lab_F -> N_F
        / (w_F / P_CES_F),
        euler_D,                                                 # 5 D deposit Euler
        # 6 deposit-UIP; residuals 5, 6 and 10 are the JOINT deposit block pinning
        # (rdep_D, rdep_F, A_D) under one union-wide market. Plus Bocola's SGU debt-elastic premium on the net
        # external position (proxied by the P_D-P_F wealth imbalance; +ve => the deficit
        # side F pays more, saves more, P_F reverts up). 0 at the symmetric SS (P_D=P_F).
        # 6 deposit-UIP, OR the union-rate diagnostic. cal["union_nominal_rate"] = True
        # replaces real UIP with a literal rdep_D = rdep_F. THIS IS NOT AN EQUILIBRIUM:
        # with own-good deposit legs the real-exchange-rate valuation profit is
        # unassigned and Walras leaks (watch goods_F in the accuracy report). It is the
        # cheap FALSIFICATION TEST in docs/nominal_block_scope.md S4 -- if pinning the
        # two rates together does NOT move the output response most of the way to
        # -0.15%, the terms-of-trade channel is not what the arithmetic says and the
        # nominal block should not be built. Diagnostic only; default off.
        ((rdep_D - rdep_F) / (1.0 + cal["r_dep_D_target"])
         if cal.get("union_nominal_rate", False) else
         (1.0 + rdep_D) - (1.0 + rdep_F) * Ep_next / p
         + cal.get("kappa_nfa", 0.0) * nfa_dep_D / ss["ss_firm_D"]["Y_ss"]),
        (Y_D - P_CES_D * C_D - I_D - NX_D - cal["G_D"])          # 7 goods_D -> p
        / ss["ss_firm_D"]["Y_ss"],
        # 8 D-bank D-bond FOC -> Q_bD. In return units (divide by the O(1) kernel), the
        # same normalisation the capital Euler uses.
        (E_Om_payD - dmd_D * Q_bD) / E_Om_D,
        # 9 F-bank D-bond FOC (with the adjustment cost) -> b_DF: the FOREIGN demand
        # schedule. Together with 8 and b_DD = B' - b_DF the D-sovereign market clears.
        (E_Om_payD_F - dmd_F * Q_bD * adj_F) / E_Om_F,
        euler_F,                                                 # 10 F deposit Euler
        # 11 UNION DEPOSIT CLEARING: total household saving funds the whole union book,
        # in D-good units, scaled by SS deposits so it is O(1) like every other residual.
        (save_union - dep_union) / ((1.0 + sz) * bkD["Dep_supply_ss"]),
        # 12 F-bank F-bond FOC -> Q_bF (the home leg prices the safe perpetuity)
        (E_Om_payF - dmd_F_home * Q_bF) / E_Om_F,
        # 13 D-bank F-bond FOC (with the adjustment cost) -> b_FD: the FOREIGN demand
        # schedule for the safe bond. With 12 and b_FF = B_F - b_FD the F market clears.
        (E_Om_payF_D - dmd_F_for * Q_bF * adj_D) / E_Om_D,
    ])
    out = dict(E_payD=E_payD, E_payD_nodef=E_payD_nodef, E_Om_payD=E_Om_payD,
               lam_bD_mu_D=lbDD * mu_D,
               # the F-side legs, so the D-F SPREAD can be decomposed and not just the
               # D yield: the spread is what the policy targets, and it moves only by
               # the difference between two yields that can both shift
               E_Om_F=E_Om_F, E_payF=_E(payF), E_Om_payF=E_Om_payF,
               lam_bF_mu_F=lbFF * mu_F, mu_F_=mu_F,
               mu_D=mu_D, mu_F=mu_F, n_D=n_D, n_F=n_F, Y_D=Y_D, I_D=I_D, I_F=I_F,
               alpha_D=alpha_D_new, alpha_F=alpha_F_new,
               Q_bF=Q_bF, C_D=C_D, C_F=C_F,
               b_F_D_new=b_F_D_new, b_F_F_new=b_F_F_new,
               b_D_D_new=b_D_D_new, b_D_F_new=b_D_F_new, nfa_dep_D=nfa_dep_D,
               W_D=W_D, W_F=W_F, Wp_D=Wp_D, Vp_dep=Vp_dep, A_D=A_D, A_F=A_F,
               Q_bD=Q_bD, Bp_tot=Bp_D, dep_union=dep_union, save_union=save_union,
               phi=phi, B_D=B_D, m_ltro_D=m_ltro_D, m_ltro_F=m_ltro_F,
               n_IC_D=n_IC_D, n_IC_F=n_IC_F, lev_IC_D=lev_IC_D, lev_IC_F=lev_IC_F,
               inc_D=inc_D, inc_F=inc_F, w_D=w_D, dep_D=dep_D, dep_F=dep_F,
               Pp_D=Pp_D, Pp_F=Pp_F, Bp_D=Bp_D, slack_D=slack_D, slack_F=slack_F,
               # accounting legs consumed by the output decomposition and the
               # heterogeneous-agent welfare overlay (never by the residuals)
               N_D=N_D, Kap_prod_D=K_D, Z_D=Z_D, Kp_D=Kp_D, P_CES_D=P_CES_D,
               E_Om_D=E_Om_D, r_wc_D=r_wc_D, wedge_sp_D=lKD * mu_D / E_Om_D,
               rdep_D=rdep_D, rdep_F=rdep_F, Div_D=Div_D, Tax_D=Tax_D, p=p, Y_F=Y_F,
               r_wc_F=r_wc_F, L_wc_D=L_wc_D, L_wc_F=L_wc_F,
               # DIAGNOSTICS, never residuals. euler_F is the F household's
               # intertemporal condition, which this system COMPUTES AND DROPS: with
               # A_F = dep_F/P_CES_F the F household is force-fed its own bank's funding
               # need and rdep_F comes from UIP, so nothing holds euler_F to zero. Its
               # size is the measure of how far national deposit clearing is from the
               # union deposit market the model is documented to run. C_D_terms are the
               # three legs of C = carried claim + income - new deposits, which is a
               # ~0.78 residual of ~8-sized gross flows.
               euler_F_resid=euler_F, euler_D_resid=euler_D,
               C_D_terms=(P_D / P_CES_D, inc_D, A_D))
    return res, out
