# SINGLE-POINT PERIOD MAP FOR THE RECURSIVE SOLUTION (THE IMAGE OF _inner_economy).
# Evaluates the full economy at ONE grid point (d, S) in QUANTITY FORM (realized
# payoffs = quantities x current prices).
#
# STATE (12): [K_D, K_F, P_D, P_F, b_DD, b_DF, b_FD, V_dep, s, Z_D, M_cb, O_cb] -- see
#   state_grid.
#   b_DD/b_DF are the two banks' CARRIED holdings of the D sovereign (B_D is their
#   sum). They were one state while the split was a fixed SS share; once both banks
#   bid through their own FOCs, last period's split is part of the state.
#   V_dep is the carried CROSS-BORDER deposit position, W_D - P_D. Each household's
#   claim used to be identified with its own bank's obligation -- true only under
#   NATIONAL clearing, where it is force-fed that bank's funding need. Under the union
#   deposit market they differ, and V is what they differ by: W_D = P_D + V,
#   W_F = P_F - V/p, so the union identity holds by construction.
#
# UNKNOWNS (14): [N_D, N_F, Kp_D, Kp_F, rdep_D, rdep_F, p, Q_bD, b_DF', Q_bF, b_FD',
#   A_D, A_F, m_cb]
# RESIDUALS (14): 2 bank capital-Euler (occasionally binding), 2 labour, BOTH household
#   deposit Eulers, deposit-UIP, goods-D, BOTH banks' D-bond FOCs, union deposit
#   clearing, BOTH banks' F-bond FOCs, and the TPI purchase rule. The D-sovereign is not
#   force-fed to anyone: the two bond FOCs are demand schedules, b_DD = B' - b_DF - M
#   clears the market, and Q_bD is the price that does it.
#
# THE TPI BOOK. The Eurosystem holds M of the D bond, bought at market prices and held
#   to maturity, M = (1-delta_b)*surv*M_lag + m. It pays for them with a SAFE claim
#   Z = Q_bD*M on itself, held by the D banks and paying (1+rdep_D)*Z next period (the
#   state O_cb). A purchase therefore SWAPS bonds for Z at the same value inside the
#   bank's book, so assets, the divertable base (single lambda), deposits and P' do not
#   move within the period: the whole effect is the RISK TRANSFER to next period, where
#   the bank is owed O instead of the bond's payoff. The Eurosystem keeps no equity and
#   remits its P&L Pi = Xi*M_lag - O_lag every period, a share tpi_key_D to the D
#   treasury and the rest to F. With M = O = m = 0 every new term is an exact 0.0. The banker valuations alpha are still READ OFF the recursions
#   given a FROZEN continuation; under the collocation solve they are unknowns with their
#   own identity residual, so the freeze is exact at every node.
#
# OCCASIONALLY-BINDING IC -- Bocola Prop. 1 CLOSED FORM: mu = max{1 - E[Om]*R*n /
# (lambda*divertable assets), 0}, explicit and bounded in [0,1); the capital Euler
# E[Om(R_K-R)] = lambda_K*mu is the residual pinning K'.
#
# REGIMES. The index the rules are stored under is the default indicator d' (0 = no
# default, 1 = the haircut is realised); _regime_weights turns it into one probability
# per (s' node, regime') cell, and every expectation below is the SAME weighted sum over
# that list. The SDF in each regime is the genuine state-contingent GHH-composite kernel
# against that regime's own continuation. A regime with zero weight is never evaluated and
# aliases regime 0, so pi = 0 (no_default) nests the risk-free model EXACTLY.
#
# THE ANATOMY. point_residuals runs the stages below in order on one per-point record v;
# each stage reads what earlier stages wrote and adds its own objects:
#   _read_state_and_regime   state, unknowns, the default regime
#   _production              firms + capital blocks, CES prices (current period)
#   _government              Bohn tax on the surviving stock, new issuance
#   _balance_sheets          bank payoffs, net worth, deposits, household claims
#   _continuation            quadrature over (s' node x regime') on the frozen rules
#   _discount_kernels        GHH-composite SDFs and the banker's Omega per regime
#   _incentive_constraint    Bocola closed-form mu, capital Euler, franchise value alpha
#   _bond_demands            both banks' FOCs for both sovereigns
#   _wages_and_dividends     working-capital wedge, wage, dividends
#   _households              budget, consumption, deposit Eulers, trade, E[p']
#   _residual_vector         the 13 equilibrium conditions
#   _outputs                 everything stored as rules or read by the reporting layer
#
# Tier-3 cuts (documented, reversible): B_F at SS; the wc-wedge rate component uses
# current rdep.
from types import SimpleNamespace

import numpy as np

from global_projection.blocks.firms import solve_firm_path, markup_ss
from global_projection.blocks.capital import solve_capital_path
from global_projection.blocks.trade import ces_price, import_demand, trade_balance, size_ratio
from global_projection.solver_recursive.state_grid import default_prob

from global_projection.solver_recursive.state_grid import (IK_D, IK_F, IP_D, IP_F, IBDD, IBDF,
                                         IBFD, IV, IS, IZ, IM, IO, NSTATE)



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
# upper cap on mu: keeps alpha = E_Om R/(1-mu) finite when a default-regime net worth
# goes negative; smoothed because unlike the KKT switch it is a guard, not economics
_MU_CAP = 0.95


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


def _regime_weights(wq, pd, reg):
    # PROBABILITY OF EACH (s' NODE x REGIME') CELL: d' = 1 with the priced default
    # probability pd. The cells sum to wq, so the quadrature measure is preserved.
    return [wq * (pd if d_n else (1.0 - pd)) for d_n in reg]


def _expect(wgt, vals):
    # EXPECTATION OVER (s' NODE x REGIME'): vals carries ONE ARRAY PER REGIME.
    return float(sum(np.dot(wgt[j], vals[j]) for j in range(len(wgt))))


def _ghh(C, N, chi, fr, floor, eps):
    # GHH CONSUMPTION COMPOSITE x = C - v(N), SMOOTH-floored.
    # A hard np.maximum(., 1e-9) is a plateau with ZERO gradient, and x_ss = 0.391 is
    # eight orders above it: a corner iterate that pushed x down landed on the plateau,
    # the FD-Jacobian saw nothing, and the point stayed stuck at |F| ~ 4e8 forever.
    # Flooring at a fraction of C_ss with the standard smoothing keeps a usable
    # derivative and is inactive anywhere near the ergodic set, where x/C = 0.50.
    return _smax(C - chi * N ** (1 + 1 / fr) / (1 + 1 / fr), floor, eps)


def _jermann_floor(K, cal, c):
    # LOWEST FEASIBLE NEXT CAPITAL UNDER JERMANN ADJUSTMENT COSTS, times a 1e-4 margin.
    delta, ksi = cal[f"delta_{c}"], cal[f"ksi_{c}"]
    floor_ratio = (1.0 - delta) - delta * ksi / (1.0 - ksi)
    return 1.0001 * floor_ratio * K


def _firm_capital(N, K, Kp, Z, cal, c):
    # ELEMENTWISE FIRM + CAPITAL BLOCK AT ONE POINT.
    # Kp is floored to the Jermann-feasible band exactly as _cont_capital floors the
    # continuation: under the GLOBAL collocation solve an exception inside one
    # finite-difference column kills the whole Jacobian, so the period map has to be
    # evaluable everywhere. This is Bocola's own device (his residual_model.m floors
    # investment and net worth). The guard is slack by ~5% of K and never binds near
    # the fixed point.
    Kp = max(Kp, _jermann_floor(K, cal, c))
    f = solve_firm_path(np.array([N]), np.array([K]), np.array([Z]), cal, c)
    cp = solve_capital_path(np.array([Kp]), K, 1.0, f["mpk"], cal, c,
                            Kap_lag_path=np.array([K]))
    return (f["Y"][0], f["w"][0], f["mpk"][0], cp["Q"][0], cp["rk"][0],
            cp["I"][0], cp["cap_profit"][0])


def _cont_capital(N_next, Kp, Kpp, Z, cal, c):
    # NEXT-PERIOD mpk' AND Q_K' FROM CONTINUATION RULES (for rk_next = t->t+1).
    # The continuation next-capital is clipped to the same Jermann band so a transient
    # iterate that dips too low does not raise inside the capital block.
    Kpp = np.maximum(Kpp, _jermann_floor(Kp, cal, c))
    f = solve_firm_path(N_next, np.full_like(N_next, Kp),
                        np.full_like(N_next, Z), cal, c)
    cp = solve_capital_path(Kpp, Kp, 1.0, f["mpk"], cal, c,
                            Kap_lag_path=np.full_like(Kpp, Kp))
    return f["mpk"], cp["Q"], f["Y"]


def _read_state_and_regime(S, d, x, cont, cal, ss):
    # UNPACK STATE AND UNKNOWNS, RESOLVE THE DEFAULT REGIME, SET THE PER-POINT FLOORS.
    v = SimpleNamespace()
    (v.N_D, v.N_F, v.Kp_D, v.Kp_F, v.rdep_D, v.rdep_F, v.p,
     v.Q_bD, v.b_DF_new, v.Q_bF, v.b_FD_new, v.A_D, v.A_F, v.m_cb) = x
    (v.K_D, v.K_F, v.P_D, v.P_F, v.b_D_D_lag, v.b_D_F_lag, v.b_F_D_lag,
     v.V_dep, v.s, v.Z_D, v.M_lag, v.O_lag) = S
    v.S = S
    # GROSS debt: the banks' holdings plus the Eurosystem's
    v.B_D = v.b_D_D_lag + v.b_D_F_lag + v.M_lag
    # the regime index is the default indicator (decision_rules.regime_table)
    v.reg, v.nreg = cont.reg, cont.n_regimes
    v.d = d
    v.d_reg = v.reg[d]
    # every variable is PER CAPITA of its own country; sz = size_F/size_D is the only
    # place the asymmetry enters, and sovereign holdings are carried in the ISSUER's
    # per-capita units (the F bank's own book holds b_D_F/sz of the D bond)
    v.sz = size_ratio(cal)
    v.surv = 1.0 - float(v.d_reg) * (1.0 - cal["recovery_rate_D"])
    # carried household claim: V = W_D - P_D, so the union identity holds by construction
    v.W_D = v.P_D + v.V_dep
    v.bkD, v.bkF = ss["ss_bank_D"], ss["ss_bank_F"]
    # GHH composite floors: 5% of SS consumption (x_ss/C_ss = 0.50, so ~10% of x_ss)
    v.X_FLOOR_D, v.X_FLOOR_F = 0.05 * ss["C_D_ss"], 0.05 * ss["C_F_ss"]
    v.X_EPS = 1e-3 * ss["C_D_ss"]
    v.Sm = np.atleast_2d(S)
    return v


def _production(v, cal):
    # FIRMS + CAPITAL AT THE CURRENT STATE (Z_D is the TFP state; F is never shocked).
    (v.Y_D, v.w_D0, v.mpk_D, v.Q_D, _, v.I_D, v.capprof_D) = _firm_capital(
        v.N_D, v.K_D, v.Kp_D, v.Z_D, cal, "D")
    (v.Y_F, v.w_F0, v.mpk_F, v.Q_F, _, v.I_F, v.capprof_F) = _firm_capital(
        v.N_F, v.K_F, v.Kp_F, cal["Z_ss_F"], cal, "F")
    v.P_CES_D = ces_price(np.array([v.p]), cal, "D")[0]
    v.P_CES_F = ces_price(np.array([v.p]), cal, "F")[0]


def _government(v, cal, ss):
    # BOHN/BOCOLA FISCAL RULE ON THE SURVIVING STOCK, LINEAR IN THE DEBT LEVEL.
    # gamma_tau is solved from a target debt root (government.govt_steady_state): a
    # UNIT level coefficient made the default event a 53%-of-GDP fiscal windfall, and
    # an ELASTICITY rule left the debt root at 1.0002 (the IRF walked B_D into the box
    # wall) or, raised, a convexity the mu=1 basis cannot carry. Linear in the level is
    # exactly representable and its root is the calibrated object.
    # THE ANCHOR IS REGIME-DEPENDENT: with a fixed b_gov_ss the haircut became a
    # tax-cut windfall (taxes -3.44% of Y on impact, C_D RISING with priced default),
    # so the default regime re-anchors to the post-haircut stock, exactly as
    # government.govt_transition does.
    if v.d_reg:
        anchor = cal["recovery_rate_D"] * ss["gs_D"]["b_gov_ss"]
        Tax_base = cal["G_D"] + cal["delta_b_D"] * anchor * (1.0 - ss["gs_D"]["Q_B_ss"])
    else:
        anchor, Tax_base = ss["gs_D"]["b_gov_ss"], ss["gs_D"]["Tax_ss"]
    v.Tax_D = Tax_base + ss["gs_D"]["gamma_tau"] * (v.B_D * v.surv - anchor)
    # THE EUROSYSTEM's P&L: what its bonds paid (coupon + the surviving stock at today's
    # price, haircut included -- it is pari passu) less what it owes the banks. Remitted
    # by capital key: the D share lowers D issuance, the F share lowers F taxes (B_F is
    # fixed), per F capita and in F goods.
    db_D = cal["delta_b_D"]
    v.Xi_D = v.surv * (db_D + (1.0 - db_D) * v.Q_bD)
    v.Pi_cb = v.Xi_D * v.M_lag - v.O_lag
    kD = cal["tpi_key_D"]
    v.Tax_F = ss["gs_F"]["Tax_ss"] - (1.0 - kD) * v.Pi_cb / (v.sz * v.p)
    v.coupon_D = cal["delta_b_D"] * v.B_D * v.surv
    v.new_D = (cal["G_D"] + v.coupon_D - v.Tax_D - kD * v.Pi_cb) / v.Q_bD
    v.Bp_D = (1.0 - cal["delta_b_D"]) * v.B_D * v.surv + v.new_D


def _balance_sheets(v, cont, cal):
    # BANK PAYOFFS, NET WORTH, DEPOSIT OBLIGATIONS AND HOUSEHOLD CLAIMS (alpha/mu-free).
    db_D, db_F = cal["delta_b_D"], cal["delta_b_F"]
    # b_F_D_lag is a state; the F bank holds the remainder of the FIXED F stock
    v.b_F_F_lag = cal["B_gov_F_ss"] - v.b_F_D_lag / v.sz
    # current bank valuations are the FROZEN previous-iterate rules at this state
    # (breaks the S'<->Q_b knot; they coincide at the fixed point). Q_bD is solved, so
    # it is not read off the frozen iterate.
    v.alpha_D_cur = float(cont.eval("alpha_D", v.d, v.Sm)[0])
    v.alpha_F_cur = float(cont.eval("alpha_F", v.d, v.Sm)[0])
    # WORKING-CAPITAL LOANS ARE BANK ASSETS (Bocola SV.C). Their SIZE needs r_wc, which
    # needs mu, which needs the balance sheet: take r_wc from the FROZEN previous iterate
    # for the loan QUANTITY. The labour residual uses the contemporaneous r_wc; the two
    # coincide at the fixed point.
    v.r_wc_D_cur = float(cont.eval("r_wc_D", v.d, v.Sm)[0])
    v.r_wc_F_cur = float(cont.eval("r_wc_F", v.d, v.Sm)[0])
    zD_, zF_ = cal["zeta_wc_D"], cal["zeta_wc_F"]
    v.L_wc_D = zD_ * (v.w_D0 / (1.0 + zD_ * v.r_wc_D_cur)) * v.N_D
    v.L_wc_F = zF_ * (v.w_F0 / (1.0 + zF_ * v.r_wc_F_cur)) * v.N_F

    payD_now = v.Xi_D
    payF_now = db_F + (1.0 - db_F) * v.Q_bF
    # the D bank is also owed O_lag on last period's TPI claim
    X_D = ((v.mpk_D + (1.0 - cal["delta_D"]) * v.Q_D) * v.K_D
           + payD_now * v.b_D_D_lag + v.p * payF_now * v.b_F_D_lag + v.O_lag)
    X_F = ((v.mpk_F + (1.0 - cal["delta_F"]) * v.Q_F) * v.K_F
           + payF_now * v.b_F_F_lag + payD_now * (v.b_D_F_lag / v.sz) / v.p)
    v.ng_D, v.ng_F = X_D - v.P_D, X_F - v.P_F
    # MARKET CLEARING BY CONSTRUCTION: whatever the F bank does not take, the D bank
    # holds (and the same for the F sovereign). Floored so a transient iterate cannot
    # hand a bank a negative book, which would flip the sign of lev.
    v.b_D_F_new = _smax(v.b_DF_new, 1e-4, 1e-5)
    # THE TPI BOOK, held to maturity, and the safe claim it is paid for with
    v.M_cb_new = (1.0 - db_D) * v.surv * v.M_lag + v.m_cb
    v.Z_cb = v.Q_bD * v.M_cb_new
    v.b_D_D_new = _smax(v.Bp_D - v.b_D_F_new - v.M_cb_new, 1e-4, 1e-5)
    v.b_F_D_new = _smax(v.b_FD_new, 1e-4, 1e-5)
    v.b_F_F_new = _smax(cal["B_gov_F_ss"] - v.b_F_D_new / v.sz, 1e-4, 1e-5)
    # end-of-period portfolio valued at current PRICES (Q_b), not payoffs
    v.assets_D = (v.Q_D * v.Kp_D + v.Q_bD * v.b_D_D_new + v.p * v.Q_bF * v.b_F_D_new
                  + v.L_wc_D + v.Z_cb)
    v.assets_F = (v.Q_F * v.Kp_F + v.Q_bF * v.b_F_F_new
                  + v.Q_bD * (v.b_D_F_new / v.sz) / v.p + v.L_wc_F)
    v.n_D = (1.0 - cal["f_D"]) * v.ng_D + cal["omega_ent_D"] * v.assets_D
    v.n_F = (1.0 - cal["f_F"]) * v.ng_F + cal["omega_ent_F"] * v.assets_F
    # BOCOLA FEASIBILITY FLOOR (his N_tom = max(.,0.65)): inactive in the ergodic region
    # (nw_floor_frac = 0 => baseline unchanged), catches only the deep default corners
    nwf = cal.get("nw_floor_frac", 0.0)
    if nwf > 0.0:
        v.n_D = _smax(v.n_D, nwf * v.bkD["n_ss"])
        v.n_F = _smax(v.n_F, nwf * v.bkF["n_ss"])
    v.dep_D, v.dep_F = v.assets_D - v.n_D, v.assets_F - v.n_F
    # deposits fund the whole book INCLUDING the working-capital loan, repaid at the
    # lending rate -- Bocola: P' = R*(Q*K' + q*B' + L - N') - R_W*L
    v.Pp_D = (1.0 + v.rdep_D) * v.dep_D - (1.0 + v.r_wc_D_cur) * v.L_wc_D
    v.Pp_F = (1.0 + v.rdep_F) * v.dep_F - (1.0 + v.r_wc_F_cur) * v.L_wc_F
    # UNION DEPOSIT MARKET: A_D, A_F are the households' chosen savings; one union-wide
    # clearing in D-good units; the gap A_D*P_CES_D - dep_D is the cross-border position
    v.dep_union = v.dep_D + v.sz * v.p * v.dep_F
    v.save_union = v.A_D * v.P_CES_D + v.sz * v.p * v.A_F * v.P_CES_F
    v.nfa_dep_D = v.A_D * v.P_CES_D - v.dep_D
    # the household's carried claim mirrors the bank's obligation, WC netting included:
    # without the same deduction W_D drifts from P_D by (1+r_wc)*L_wc every period
    v.Wp_D = (1.0 + v.rdep_D) * v.A_D * v.P_CES_D - (1.0 + v.r_wc_D_cur) * v.L_wc_D
    # the carried cross-border position: both legs carry the same WC deduction, so it
    # cancels and V' is this period's cross-border flow grossed up at the D deposit rate
    v.Vp_dep = (1.0 + v.rdep_D) * v.nfa_dep_D
    # what the Eurosystem will owe the D banks next period
    v.Op_cb = (1.0 + v.rdep_D) * v.Z_cb


def _next_states(v, sproc, eps):
    # THE m NEXT-PERIOD STATES, ONE PER s' QUADRATURE NODE.
    Z_next = (1.0 - sproc["rho_z"]) * sproc["z_star"] + sproc["rho_z"] * v.Z_D
    s_next = ((1.0 - sproc["rho_s"]) * sproc["s_star"] + sproc["rho_s"] * v.s
              + sproc["sigma_s"] * eps)
    Sn = np.empty((eps.size, NSTATE))
    Sn[:, IK_D] = v.Kp_D; Sn[:, IK_F] = v.Kp_F
    Sn[:, IP_D] = v.Pp_D; Sn[:, IP_F] = v.Pp_F
    Sn[:, IBDD] = v.b_D_D_new; Sn[:, IBDF] = v.b_D_F_new
    Sn[:, IBFD] = v.b_F_D_new
    Sn[:, IV] = v.Vp_dep
    Sn[:, IS] = s_next
    Sn[:, IZ] = Z_next
    Sn[:, IM] = v.M_cb_new; Sn[:, IO] = v.Op_cb
    return Sn, Z_next


def _regime_continuation(v, j, Sn, Z_next, cont, cal):
    # CONTINUATION RULE VALUES + NEXT-PERIOD CAPITAL RETURN IN REGIME j.
    r = cont.eval_all(j, cont.grid.clip(Sn))
    # guard continuation outputs before they enter fractional powers / sqrt
    # (a deep default-regime iterate can push N, p, C negative -> NaN)
    for k in ("N_D", "N_F"):
        r[k] = np.maximum(r[k], 0.05)
    for k in ("C_D", "C_F"):
        r[k] = np.maximum(r[k], 1e-3)
    r["p"] = np.maximum(r["p"], 1e-2)
    mpkD_n, QKD_n, _ = _cont_capital(r["N_D"], v.Kp_D, r["Kp_D"], Z_next, cal, "D")
    mpkF_n, QKF_n, _ = _cont_capital(r["N_F"], v.Kp_F, r["Kp_F"], cal["Z_ss_F"], cal, "F")
    rkD_n = (mpkD_n + (1.0 - cal["delta_D"]) * QKD_n) / v.Q_D - 1.0
    rkF_n = (mpkF_n + (1.0 - cal["delta_F"]) * QKF_n) / v.Q_F - 1.0
    return r, rkD_n, rkF_n


def _continuation(v, cont, cal, sproc, n_gh, no_default):
    # QUADRATURE WEIGHTS OVER (s' NODE x REGIME') AND THE CONTINUATION IN EACH REGIME.
    # A regime carrying ZERO weight is never evaluated and is ALIASED to regime 0: that is
    # what makes pi = 0 nest the risk-free model EXACTLY rather than to
    # solver tolerance, and it saves a full grid interpolation plus two capital blocks.
    eps, wq = gh_nodes(n_gh)
    pd = 0.0 if no_default else float(default_prob(v.s))
    Sn, Z_next = _next_states(v, sproc, eps)
    v.wgt = _regime_weights(wq, pd, v.reg)
    v.R, v.RKD, v.RKF = [], [], []
    for j in range(v.nreg):
        if j > 0 and not np.any(v.wgt[j] > 0.0):
            v.R.append(v.R[0]); v.RKD.append(v.RKD[0]); v.RKF.append(v.RKF[0])
            continue
        rj, rkDj, rkFj = _regime_continuation(v, j, Sn, Z_next, cont, cal)
        v.R.append(rj); v.RKD.append(rkDj); v.RKF.append(rkFj)


def _discount_kernels(v, cont, cal):
    # BRANCH SDFs AND THE BANKER'S OMEGA -- the genuine state-contingent GHH kernel.
    # Bocola's banker discounts with Lambda' = beta*c/c', which is what makes his
    # constraint TIGHTEN with risk: next-period consumption falls, Omega falls and
    # mu = 1 - E[Om]*R*n/lev RISES. A CONSTANT beta_inter cannot fall, so the constraint
    # went SLACK exactly when it should bind (measured at three calibrations). The kernel
    # is evaluated on the frozen current rules against each regime's continuation.
    # NB the F kernel floors at the D floor (_X_FLOOR_D), as it always has; the two are
    # equal under the symmetric per-capita steady state.
    nreg = v.nreg
    frD, frF = cal["frisch_D"], cal["frisch_F"]
    sgD, sgF = cal["sigma_D"], cal["sigma_F"]
    x_cur_D = _ghh(float(cont.eval("C_D", v.d, v.Sm)[0]), v.N_D, cal["chi_D"], frD,
                   v.X_FLOOR_D, v.X_EPS)
    x_cur_F = _ghh(float(cont.eval("C_F", v.d, v.Sm)[0]), v.N_F, cal["chi_F"], frF,
                   v.X_FLOOR_D, v.X_EPS)
    v.XN_D = [_ghh(v.R[j]["C_D"], v.R[j]["N_D"], cal["chi_D"], frD, v.X_FLOOR_D, v.X_EPS)
              for j in range(nreg)]
    XN_F = [_ghh(v.R[j]["C_F"], v.R[j]["N_F"], cal["chi_F"], frF, v.X_FLOOR_D, v.X_EPS)
            for j in range(nreg)]
    Lam_D = [cal["beta_inter_D"] * (x_cur_D / v.XN_D[j]) ** sgD for j in range(nreg)]
    Lam_F = [cal["beta_inter_F"] * (x_cur_F / XN_F[j]) ** sgF for j in range(nreg)]
    v.Om_D = [Lam_D[j] * (cal["f_D"] + (1 - cal["f_D"]) * v.R[j]["alpha_D"])
              for j in range(nreg)]
    v.Om_F = [Lam_F[j] * (cal["f_F"] + (1 - cal["f_F"]) * v.R[j]["alpha_F"])
              for j in range(nreg)]
    v.E_Om_D = _expect(v.wgt, v.Om_D)
    v.E_Om_F = _expect(v.wgt, v.Om_F)


def _incentive_constraint(v, cal):
    # BOCOLA PROP. 1 CLOSED-FORM MULTIPLIER, CAPITAL EULER AND FRANCHISE VALUE.
    lKD, lKF = cal["lambda_K_D"], cal["lambda_K_F"]
    lbDD, lbFF = cal["lambda_bD_D"], cal["lambda_bF_F"]
    # divertable assets (same leverage term as the FB slack), frozen-price valued
    # the TPI claim Z sits in the base at the bond's lambda, so a swap leaves lev unchanged
    lev_D = max(lKD * v.Q_D * v.Kp_D + lbDD * v.Q_bD * v.b_D_D_new
                + cal["lambda_bF_D"] * v.p * v.Q_bF * v.b_F_D_new
                + lKD * v.L_wc_D + lbDD * v.Z_cb, 1e-6)
    lev_F = max(lKF * v.Q_F * v.Kp_F + lbFF * v.Q_bF * v.b_F_F_new
                + cal["lambda_bD_F"] * v.Q_bD * (v.b_D_F_new / v.sz) / v.p
                + lKF * v.L_wc_F, 1e-6)
    v.lev_D, v.lev_F = lev_D, lev_F
    # closed-form mu in [0, _MU_CAP]: the lower max is Bocola's KKT switch (hard -- it IS
    # the complementarity); the upper cap is a smoothed guard
    v.mu_D = float(_smin(max(1.0 - v.E_Om_D * (1.0 + v.rdep_D) * v.n_D / v.lev_D, 0.0),
                         _MU_CAP, _GUARD_EPS))
    v.mu_F = float(_smin(max(1.0 - v.E_Om_F * (1.0 + v.rdep_F) * v.n_F / v.lev_F, 0.0),
                         _MU_CAP, _GUARD_EPS))
    # capital-Euler surplus E[Om(R_K - R)] - lambda_K*mu, normalised by the O(1) kernel
    # E_Om so the residual is in return units (dividing by lambda_K*mu_ss amplifies ~130x)
    E_Om_rk_D = _expect(v.wgt, [v.Om_D[j] * (v.RKD[j] - v.rdep_D) for j in range(v.nreg)])
    E_Om_rk_F = _expect(v.wgt, [v.Om_F[j] * (v.RKF[j] - v.rdep_F) for j in range(v.nreg)])
    v.cap_eul_D = (E_Om_rk_D - lKD * v.mu_D) / v.E_Om_D
    v.cap_eul_F = (E_Om_rk_F - lKF * v.mu_F) / v.E_Om_F
    # alpha = E_Om*R/(1-mu), smooth-capped. The recursion's slope beta*(1-f)*R/(1-mu)
    # exceeds 1 once mu > ~0.04, so above that the cap is what arrests a divergent fixed
    # point; tunable (alpha_cap) so the truncation can be measured, not assumed.
    alpha_cap = float(cal.get("alpha_cap", 40.0))
    v.alpha_D_new = float(_sclip(v.E_Om_D * (1.0 + v.rdep_D) / (1.0 - v.mu_D), 0.05, alpha_cap))
    v.alpha_F_new = float(_sclip(v.E_Om_F * (1.0 + v.rdep_F) / (1.0 - v.mu_F), 0.05, alpha_cap))
    # slack is read off the SAME constraint the multiplier came from, so mu*slack = 0
    v.slack_D = v.alpha_D_cur * v.n_D - v.lev_D
    v.slack_F = v.alpha_F_cur * v.n_F - v.lev_F


def _bond_demands(v, cal):
    # BOTH BANKS' FOCs FOR BOTH SOVEREIGNS, AND THE BOND-PRICE DECOMPOSITION LEGS.
    # D bond: the gross HM perpetuity payoff at each regime's continuation price, times
    # the haircut in the regimes that default. F bond: safe.
    nreg, R, wgt = v.nreg, v.R, v.wgt
    db_D, db_F = cal["delta_b_D"], cal["delta_b_F"]
    hc = [cal["recovery_rate_D"] if d_n else 1.0 for d_n in v.reg]
    payD_gross = [db_D + (1.0 - db_D) * R[j]["Q_bD"] for j in range(nreg)]
    payD = [hc[j] * payD_gross[j] for j in range(nreg)]
    payF = [db_F + (1.0 - db_F) * R[j]["Q_bF"] for j in range(nreg)]
    v.E_Om_payD = _expect(wgt, [v.Om_D[j] * payD[j] for j in range(nreg)])
    # DECOMPOSITION LEGS (diagnostics only). The D bank's FOC is
    # E[Om*pay] = Q*(E[Om]*R + lambda_bD*mu), so Q splits into actuarial discounting of
    # E[pay], the RISK premium E[Om*pay]/(E[Om]*E[pay]) and the LIQUIDITY premium
    # E[Om]R/(E[Om]R + lambda*mu); E_payD_nodef isolates the expected loss (a defaulting
    # regime is replaced by the plain no-default one). The legs are exactly additive in
    # logs -- Bocola's Table 4 split.
    v.E_payD = _expect(wgt, payD)
    v.E_payD_nodef = _expect(wgt, [payD_gross[j] if not v.reg[j] else payD_gross[0]
                                   for j in range(nreg)])
    v.E_payF = _expect(wgt, payF)
    v.E_Om_payF = _expect(wgt, [v.Om_F[j] * payF[j] for j in range(nreg)])
    # the D bank's F-bond FOC (home leg prices Q_bF; foreign leg pins b_FD with the
    # cross-border adjustment cost psi_bF_D)
    v.E_Om_payF_D = _expect(wgt, [v.Om_D[j] * payF[j] for j in range(nreg)])
    v.dmd_F_home = v.E_Om_F * (1.0 + v.rdep_F) + cal["lambda_bF_F"] * v.mu_F
    v.dmd_F_for = v.E_Om_D * (1.0 + v.rdep_D) + cal["lambda_bF_D"] * v.mu_D
    v.adj_D = 1.0 + cal["psi_bF_D"] * (v.b_F_D_new - cal["b_F_D_ss"]) / cal["B_gov_F_ss"]
    # the D bank's D-bond FOC is an implicit residual: Q_bD is the market-clearing unknown
    v.dmd_D = v.E_Om_D * (1.0 + v.rdep_D) + cal["lambda_bD_D"] * v.mu_D
    # the F bank's D-bond FOC, with Bocola's cross-border adjustment cost psi_bD_F: without
    # it both demand curves are near-flat, the split is indeterminate and the Newton is
    # ill-conditioned; the cost is zero at the SS position
    v.E_Om_payD_F = _expect(wgt, [v.Om_F[j] * payD[j] for j in range(nreg)])
    v.dmd_F = v.E_Om_F * (1.0 + v.rdep_F) + cal["lambda_bD_F"] * v.mu_F
    v.adj_F = 1.0 + cal["psi_bD_F"] * (v.b_D_F_new - cal["b_D_F_ss"]) / cal["B_gov_D_ss"]


def _wages_and_dividends(v, cal):
    # WORKING-CAPITAL WEDGE r_wc = rdep + lambda_K*mu/E[Om], NET WAGE AND DIVIDENDS.
    lKD, lKF = cal["lambda_K_D"], cal["lambda_K_F"]
    v.r_wc_D = v.rdep_D + lKD * v.mu_D / v.E_Om_D
    v.r_wc_F = v.rdep_F + lKF * v.mu_F / v.E_Om_F
    zD, zF = cal["zeta_wc_D"], cal["zeta_wc_F"]
    v.w_D = v.w_D0 / (1.0 + zD * v.r_wc_D)
    v.w_F = v.w_F0 / (1.0 + zF * v.r_wc_F)
    wc_inc_D = zD * v.r_wc_D * v.w_D * v.N_D
    wc_inc_F = zF * v.r_wc_F * v.w_F * v.N_F
    mcD, mcF = markup_ss(cal, "D"), markup_ss(cal, "F")
    div_bank_D = cal["f_D"] * v.ng_D - cal["omega_ent_D"] * v.assets_D
    div_bank_F = cal["f_F"] * v.ng_F - cal["omega_ent_F"] * v.assets_F
    # where the working-capital financing income goes: wc_rebate = 1 hands it to
    # households, which makes the spread an intra-period transfer and (no GHH wealth
    # effect) turns the risk channel expansionary; the default 0 leaves it with the bank
    reb = float(cal.get("wc_rebate", 0.0))
    v.Div_D = (1 - mcD) * v.Y_D + v.capprof_D + div_bank_D + reb * wc_inc_D
    v.Div_F = (1 - mcF) * v.Y_F + v.capprof_F + div_bank_F + reb * wc_inc_F


def _households(v, cal, ss):
    # REP-AGENT GHH HOUSEHOLDS: BUDGET, CONSUMPTION, DEPOSIT EULERS, TRADE, E[p'].
    # beta_eff = 1/(1+rdep_ss) makes the rep-agent Euler reproduce the HA aggregate at the
    # SS; hh_T anchors the SS budget so C = C_ss, A = A_ss exactly.
    nreg, R = v.nreg, v.R
    frisch_D, frisch_F = cal["frisch_D"], cal["frisch_F"]
    sigD, sigF = cal["sigma_D"], cal["sigma_F"]
    v.inc_D = (v.w_D / v.P_CES_D) * v.N_D + (v.Div_D - v.Tax_D) / v.P_CES_D + ss.get("hh_T_D", 0.0)
    v.inc_F = (v.w_F / v.P_CES_F) * v.N_F + (v.Div_F - v.Tax_F) / v.P_CES_F + ss.get("hh_T_F", 0.0)
    # F's carried claim is its own bank's obligation LESS the cross-border position, in
    # F-good units (deriving it as the union residual inherited three states' fit error)
    v.W_F = v.P_F - v.V_dep / (v.sz * v.p)
    # smooth bounds with eps at 1e-3*C_ss: at _GUARD_EPS the bound was effectively hard
    v.C_D = float(_sclip(v.W_D / v.P_CES_D + v.inc_D - v.A_D,
                         0.15 * ss["C_D_ss"], 3.0 * ss["C_D_ss"], v.X_EPS))
    v.C_F = float(_sclip(v.W_F / v.P_CES_F + v.inc_F - v.A_F,
                         0.15 * ss["C_F_ss"], 3.0 * ss["C_F_ss"], v.X_EPS))
    # deposit Euler on the GHH composite: x^-sigma = beta_eff*E[(1+r')x'^-sigma]; the
    # D continuation composite is the kernel's own XN_D (same expression, same floor)
    xp_F = [_ghh(R[j]["C_F"], R[j]["N_F"], cal["chi_F"], frisch_F, v.X_FLOOR_F, v.X_EPS)
            for j in range(nreg)]
    rp_D = [(1.0 + v.rdep_D) * v.P_CES_D / ces_price(R[j]["p"], cal, "D") - 1.0
            for j in range(nreg)]
    rp_F = [(1.0 + v.rdep_F) * v.P_CES_F / ces_price(R[j]["p"], cal, "F") - 1.0
            for j in range(nreg)]
    beff_D = 1.0 / (1.0 + cal["r_dep_D_target"])
    beff_F = 1.0 / (1.0 + cal["r_dep_F_target"])
    E_mu_D = _expect(v.wgt, [(1.0 + rp_D[j]) * v.XN_D[j] ** (-sigD) for j in range(nreg)])
    E_mu_F = _expect(v.wgt, [(1.0 + rp_F[j]) * xp_F[j] ** (-sigF) for j in range(nreg)])
    xC_D = float(_ghh(v.C_D, v.N_D, cal["chi_D"], frisch_D, v.X_FLOOR_D, v.X_EPS))
    xC_F = float(_ghh(v.C_F, v.N_F, cal["chi_F"], frisch_F, v.X_FLOOR_F, v.X_EPS))
    v.euler_D = xC_D ** (-sigD) / (beff_D * E_mu_D) - 1.0
    v.euler_F = xC_F ** (-sigF) / (beff_F * E_mu_F) - 1.0
    IM_D = import_demand(np.array([v.p]), np.array([v.C_D]), np.array([v.P_CES_D]), cal, "D")[0]
    IM_F = import_demand(np.array([v.p]), np.array([v.C_F]), np.array([v.P_CES_F]), cal, "F")[0]
    NX_D, NX_F = trade_balance(np.array([v.p]), np.array([IM_D]), np.array([IM_F]), cal)
    v.NX_D = NX_D[0]
    # THE F GOODS MARKET, never a residual (Walras) -- the check that the TPI's
    # cross-border P&L closes through the existing budgets. Its baseline is NOT zero off
    # the SS: B_F and Tax_F are fixed (a Tier-3 cut), so the F treasury is short
    # delta_F*B_F*(Q_ss - Q_bF) whenever Q_bF moves, and goods_F carries that leak.
    v.goods_F = ((v.Y_F - v.P_CES_F * v.C_F - v.I_F - NX_F[0] - cal["G_F"])
                 / ss["ss_firm_F"]["Y_ss"])
    # E[p'] for deposit-UIP, under the same measure over regimes as every expectation
    v.Ep_next = _expect(v.wgt, [R[j]["p"] for j in range(nreg)])


def _residual_vector(v, cal, ss):
    # THE 14 EQUILIBRIUM CONDITIONS, EACH IN O(1) UNITS.
    frisch_D, frisch_F = cal["frisch_D"], cal["frisch_F"]
    # 6: deposit-UIP (with Bocola's SGU debt-elastic premium on the cross-border
    # position), OR -- cal["union_nominal_rate"] -- a literal rdep_D = rdep_F. The latter
    # is NOT an equilibrium (the real-exchange-rate valuation profit is unassigned and
    # Walras leaks); it is the falsification test in Claude files/docs/nominal_block_scope.md S4.
    if cal.get("union_nominal_rate", False):
        uip = (v.rdep_D - v.rdep_F) / (1.0 + cal["r_dep_D_target"])
    else:
        uip = ((1.0 + v.rdep_D) - (1.0 + v.rdep_F) * v.Ep_next / v.p
               + cal.get("kappa_nfa", 0.0) * v.nfa_dep_D / ss["ss_firm_D"]["Y_ss"])
    return np.array([
        v.cap_eul_D,                                                   # 1 cap Euler D -> Kp_D
        v.cap_eul_F,                                                   # 2 cap Euler F -> Kp_F
        (cal["chi_D"] * v.N_D ** (1 / frisch_D) - v.w_D / v.P_CES_D)    # 3 lab_D -> N_D
        / (v.w_D / v.P_CES_D),
        (cal["chi_F"] * v.N_F ** (1 / frisch_F) - v.w_F / v.P_CES_F)    # 4 lab_F -> N_F
        / (v.w_F / v.P_CES_F),
        v.euler_D,                                                     # 5 D deposit Euler
        uip,                                                           # 6 deposit-UIP
        (v.Y_D - v.P_CES_D * v.C_D - v.I_D - v.NX_D - cal["G_D"])       # 7 goods_D -> p
        / ss["ss_firm_D"]["Y_ss"],
        (v.E_Om_payD - v.dmd_D * v.Q_bD) / v.E_Om_D,                   # 8 D-bank D-bond -> Q_bD
        (v.E_Om_payD_F - v.dmd_F * v.Q_bD * v.adj_F) / v.E_Om_F,       # 9 F-bank D-bond -> b_DF
        v.euler_F,                                                     # 10 F deposit Euler
        (v.save_union - v.dep_union) / ((1.0 + v.sz) * v.bkD["Dep_supply_ss"]),  # 11 union clearing
        (v.E_Om_payF - v.dmd_F_home * v.Q_bF) / v.E_Om_F,              # 12 F-bank F-bond -> Q_bF
        (v.E_Om_payF_D - v.dmd_F_for * v.Q_bF * v.adj_D) / v.E_Om_D,   # 13 D-bank F-bond -> b_FD
        v.m_cb / cal["B_gov_D_ss"],                                    # 14 TPI off: m = 0
    ])


def _outputs(v, cal):
    # EVERY OBJECT STORED AS A RULE OR READ BY THE EXPERIMENTS AND REPORTING LAYER.
    lKD, lbDD, lbFF = cal["lambda_K_D"], cal["lambda_bD_D"], cal["lambda_bF_F"]
    return dict(E_payD=v.E_payD, E_payD_nodef=v.E_payD_nodef, E_Om_payD=v.E_Om_payD,
                lam_bD_mu_D=lbDD * v.mu_D,
                # the F-side legs, so the D-F SPREAD can be decomposed, not just the D yield
                E_Om_F=v.E_Om_F, E_payF=v.E_payF, E_Om_payF=v.E_Om_payF,
                lam_bF_mu_F=lbFF * v.mu_F, mu_F_=v.mu_F,
                mu_D=v.mu_D, mu_F=v.mu_F, n_D=v.n_D, n_F=v.n_F, Y_D=v.Y_D, I_D=v.I_D,
                I_F=v.I_F, alpha_D=v.alpha_D_new, alpha_F=v.alpha_F_new,
                Q_bF=v.Q_bF, C_D=v.C_D, C_F=v.C_F,
                b_F_D_new=v.b_F_D_new, b_F_F_new=v.b_F_F_new,
                b_D_D_new=v.b_D_D_new, b_D_F_new=v.b_D_F_new, nfa_dep_D=v.nfa_dep_D,
                W_D=v.W_D, W_F=v.W_F, Wp_D=v.Wp_D, Vp_dep=v.Vp_dep, A_D=v.A_D, A_F=v.A_F,
                Q_bD=v.Q_bD, Bp_tot=v.Bp_D, dep_union=v.dep_union,
                save_union=v.save_union, B_D=v.B_D,
                inc_D=v.inc_D, inc_F=v.inc_F, w_D=v.w_D, dep_D=v.dep_D, dep_F=v.dep_F,
                Pp_D=v.Pp_D, Pp_F=v.Pp_F, Bp_D=v.Bp_D, slack_D=v.slack_D, slack_F=v.slack_F,
                # the TPI book: purchases, holdings, the banks' claim, next period's
                # obligation and this period's remitted P&L
                m_cb=v.m_cb, M_cb_new=v.M_cb_new, Z_cb=v.Z_cb, Op_cb=v.Op_cb, Pi_cb=v.Pi_cb,
                Tax_F=v.Tax_F, goods_F=v.goods_F,
                # accounting legs for the output decomposition and the welfare overlay
                N_D=v.N_D, Kap_prod_D=v.K_D, Z_D=v.Z_D, Kp_D=v.Kp_D, P_CES_D=v.P_CES_D,
                E_Om_D=v.E_Om_D, r_wc_D=v.r_wc_D, wedge_sp_D=lKD * v.mu_D / v.E_Om_D,
                rdep_D=v.rdep_D, rdep_F=v.rdep_F, Div_D=v.Div_D, Tax_D=v.Tax_D, p=v.p,
                Y_F=v.Y_F, r_wc_F=v.r_wc_F, L_wc_D=v.L_wc_D, L_wc_F=v.L_wc_F,
                # diagnostics, never residuals; C_D_terms are the three legs of
                # C = carried claim + income - new deposits
                euler_F_resid=v.euler_F, euler_D_resid=v.euler_D,
                C_D_terms=(v.P_D / v.P_CES_D, v.inc_D, v.A_D))


def point_residuals(S, d, x, cont, cal, ss, sproc, n_gh=7, no_default=False):
    # THE PERIOD MAP AT ONE POINT: 14 RESIDUALS AND THE OBJECTS STORED AS RULES.
    # x follows SOLVE; cont is the FROZEN continuation RuleSet (previous iterate, or the
    # current guess under collocation).
    v = _read_state_and_regime(S, d, x, cont, cal, ss)
    _production(v, cal)
    _government(v, cal, ss)
    _balance_sheets(v, cont, cal)
    _continuation(v, cont, cal, sproc, n_gh, no_default)
    _discount_kernels(v, cont, cal)
    _incentive_constraint(v, cal)
    _bond_demands(v, cal)
    _wages_and_dividends(v, cal)
    _households(v, cal, ss)
    return _residual_vector(v, cal, ss), _outputs(v, cal)
