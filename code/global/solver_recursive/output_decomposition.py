# FACTOR DECOMPOSITION OF THE OUTPUT RESPONSE ALONG A SIMULATED RISK-SHOCK PATH.
# Two pieces, both read off ALREADY-CONVERGED decision rules (no extra solving):
#   1. `simulate` -- forward-iterate the 7-state vector under a deterministic
#      s-path, so the endogenous states (K, P, B) MOVE. The IRFs in
#      recursive_experiment.py freeze them at the SS, which zeroes the capital
#      channel by construction; here capital accumulation is alive.
#   2. `decompose_output` -- split d log Y_D into the factors that produce it.
#      Production is Cobb-Douglas and capital is predetermined, so
#          d log Y = d log Z + alpha*d log K + (1-alpha)*d log N,
#      and the GHH labour FOC chi*N^(1/frisch) = w/P_CES with the Neumeyer-Perri
#      wedge w = mc(1-alpha)(Y/N)/(1+zeta*r_wc) inverts to
#          d log N = theta*[d log Z + alpha*d log K
#                           - d log(1+zeta*r_wc) - d log P_CES],
#      theta = 1/(1/frisch + alpha). r_wc = rdep + lambda_K*mu/E[Omega] splits the
#      wedge into a DEPOSIT-RATE leg and a CREDIT-SPREAD leg (the Bocola channel);
#      the split is the symmetric two-path (Shapley) one, so it is exactly additive.
# Every deviation is taken against the SAME rules simulated with NO shock, so the
# projection's own approximation drift cancels out of the decomposition.
import numpy as np
from scipy.optimize import root

from solver_recursive.point_map import point_residuals, SOLVE7
from solver_recursive.recursive_main import ss_state
from solver_recursive.state_grid import (IK_D, IK_F, IP_D, IP_F, IBDD, IBDF, IBFD, IV,
                                         IS, IZ)

# per-period objects recorded along a simulated path
REC_KEYS = ("Y_D", "C_D", "I_D", "N_D", "Kap_prod_D", "Z_D", "P_CES_D", "p",
            "Q_bD", "mu_D", "n_D", "rdep_D", "r_wc_D", "wedge_sp_D", "E_Om_D",
            "alpha_D", "w_D", "Div_D", "Tax_D", "dep_D", "Y_F", "C_F",
            # the four legs of the D-bond FOC, for decompose_bond_price
            "E_payD", "E_payD_nodef", "E_Om_payD", "lam_bD_mu_D")


def _read(rules, cal, ss, sproc, S):
    # EVALUATE THE CONVERGED RULES AT ONE STATE (BINDING BRANCH, read_at's contract).
    Sm = np.atleast_2d(S)
    x = np.array([float(rules.eval(k, 0, Sm)[0]) for k in SOLVE7])
    res, o = point_residuals(S, 0, x, rules, cal, ss, sproc, n_gh=rules.n_gh or 5,
                             no_default=False)
    o["_x"] = x
    o["_res"] = res
    o["_resid"] = float(np.max(np.abs(res)))
    return o


def _refine(rules, cal, ss, sproc, S, x0):
    # RE-SOLVE THE SEVEN UNKNOWNS AT S, WARM-STARTED AT THE RULES. OFF BY DEFAULT.
    # Tempting (it would drive the labour-FOC residual to ~0) but WRONG at the
    # barely-binding SS: measured at mu=1 the Newton step slips onto the slack
    # equilibrium in 49 of 60 periods -- mu_D walks 0.103 -> 0 and the whole path
    # turns expansionary (Y_D +2.5% at impact against -0.004% for the rule read).
    # This is exactly the slip read_at is written to avoid. Kept, off, because the
    # obvious "fix" for a large residual channel is to re-solve, and it must not be
    # re-tried blind: the residual is real approximation error, cured by a finer
    # grid (mu=2), not by re-solving the current period.
    def f(x):
        try:
            return point_residuals(S, 0, x, rules, cal, ss, sproc, n_gh=rules.n_gh or 5,
                                   no_default=False)[0]
        except (ValueError, RuntimeError, FloatingPointError):
            return np.full(len(SOLVE7), 10.0)

    sol = root(f, x0, method="hybr", tol=1e-12)
    return sol.x, float(np.max(np.abs(sol.fun)))


def s_decay_path(sproc, s_shock, T):
    # DETERMINISTIC s-PATH: ONE-OFF INNOVATION DECAYING AT rho_s (the IRF convention).
    t = np.arange(T)
    return sproc["s_star"] + sproc["rho_s"] ** t * (s_shock - sproc["s_star"])


def simulate(rules, cal, ss, sproc, s_path, endogenous_states=True, refine=False,
             S_init=None):
    # FORWARD-SIMULATE THE ECONOMY UNDER s_path, CARRYING THE ENDOGENOUS STATES.
    # The next state is the period map's own end-of-period stocks
    # [Kp_D, Kp_F, Pp_D, Pp_F, Bp_D], so K/P/B evolve exactly as the rules imply.
    # States are clipped into the Smolyak box (Chebyshev extrapolation is unsafe);
    # `off_box` reports the worst excursion so a drifting path is never silent.
    # refine=True clears the period map exactly at each visited state (see _refine);
    # `n_slack`/`n_fail` count the periods where that was rejected and the plain
    # rule read stands instead.
    # S_init defaults to the DETERMINISTIC SS, which is NOT where the solved model
    # rests: with risk priced the no-shock economy walks to a different point (measured
    # Y_D -0.13%, n_D +2.1%, mu_D 0.0072 -> 0). Both the shocked and the reference path
    # start here, so the DIFFERENCE is valid either way -- but starting both at
    # recursive_experiment.stochastic_rest_point() decomposes the response AROUND the
    # point the model inhabits rather than along its transition to it, which is what
    # Bocola's generate_irf.m does.
    T = len(s_path)
    S0 = ss_state(ss, cal, sproc) if S_init is None else np.asarray(S_init, dtype=float)
    S = S0.copy()
    rec = {k: np.empty(T) for k in REC_KEYS}
    # the whole residual vector and, separately, the LABOUR FOC -- the only
    # equation the decomposition leans on, so its own accuracy is reported
    rec["resid"] = np.empty(T)
    rec["resid_lab"] = np.empty(T)
    lo, hi = rules.grid.lo, rules.grid.hi
    off_box = 0.0
    n_slack = n_fail = 0
    for t in range(T):
        S[IS] = s_path[t]
        span = np.maximum(hi - lo, 1e-12)
        off_box = max(off_box, float(np.max(np.maximum(
            (lo - S) / span, (S - hi) / span))))
        S = rules.grid.clip(S)[0]
        o = _read(rules, cal, ss, sproc, S)
        if refine and o["_resid"] > 1e-11:
            x1, fn = _refine(rules, cal, ss, sproc, S, o["_x"])
            if fn > 1e-9:
                n_fail += 1
            else:
                res1, o1 = point_residuals(S, 0, x1, rules, cal, ss, sproc,
                                           n_gh=rules.n_gh or 5, no_default=False)
                if o1["mu_D"] <= 1e-9:      # slipped off the binding branch
                    n_slack += 1
                else:
                    o1["_x"], o1["_res"] = x1, res1
                    o1["_resid"] = float(np.max(np.abs(res1)))
                    o = o1
        for k in REC_KEYS:
            rec[k][t] = o[k]
        rec["resid"][t] = o["_resid"]
        rec["resid_lab"][t] = abs(float(o["_res"][2]))
        if endogenous_states:
            S = S0.copy()
            S[IK_D], S[IK_F] = o["_x"][2], o["_x"][3]
            S[IP_D], S[IP_F] = o["Pp_D"], o["Pp_F"]
            S[IBDD], S[IBDF] = o["b_D_D_new"], o["b_D_F_new"]
            S[IBFD] = o["b_F_D_new"]
            S[IV] = o["Vp_dep"]
            S[IS], S[IZ] = s_path[t], S0[IZ]
        else:
            S = S0.copy()
    rec["off_box"] = off_box
    rec["n_slack"] = n_slack
    rec["n_fail"] = n_fail
    rec["s"] = np.asarray(s_path, dtype=float)
    return rec


def _wedge_shapley(rd, sp, rd_r, sp_r, zeta):
    # SYMMETRIC (SHAPLEY) SPLIT OF d log(1+zeta*r_wc) INTO ITS TWO LEGS.
    # Averaging the two orderings makes the two legs sum EXACTLY to the total,
    # with no first-order approximation and no ordering choice to defend.
    def L(a, b):
        return np.log(1.0 + zeta * (a + b))
    d_rd = 0.5 * ((L(rd, sp_r) - L(rd_r, sp_r)) + (L(rd, sp) - L(rd_r, sp)))
    d_sp = 0.5 * ((L(rd_r, sp) - L(rd_r, sp_r)) + (L(rd, sp) - L(rd, sp_r)))
    return d_rd, d_sp


# display order and labels of the output-decomposition channels
# THE BOND-PRICE LEGS, in the order they compose. Colours are reused from the output
# decomposition where the object is the same (the deposit rate and the constraint), so
# the same channel is the same hue on both figures.
BOND_CHANNELS = (("deposit_rate", "Deposit rate (discounting)"),
                 ("continuation", "Continuation price (duration)"),
                 ("expected_loss", "Expected default loss"),
                 ("risk_premium", "Risk premium"),
                 ("liquidity_premium", "Liquidity premium (bank constraint)"),
                 ("residual", "FOC residual"))


def decompose_bond_price(sim, ref, cal):
    # EXACT DECOMPOSITION OF THE D-SOVEREIGN BOND-PRICE RESPONSE.
    # The D bank's first-order condition is
    #     E[Om*pay] = Q_bD * (E[Om]*R + lambda_bD*mu),      R = 1 + rdep,
    # which factors, with no approximation, into
    #     log Q = -log R                                    the risk-free discount
    #             + log E[pay | no default]                 the continuation price
    #             + log( E[pay] / E[pay | no default] )      the expected default loss
    #             + log( E[Om*pay] / (E[Om]*E[pay]) )        the RISK premium: Om is
    #                                                        high exactly where the D
    #                                                        bond pays little, and that
    #                                                        covariance is a price
    #             + log( E[Om]R / (E[Om]R + lambda*mu) )     the LIQUIDITY premium, the
    #                                                        constraint's own wedge.
    # Every leg is a log difference against the SAME quarter of the no-shock path, so
    # the projection's own approximation drift cancels and the legs sum to the total
    # identically. The vocabulary is Bocola's Table 4 (risk premium vs the multiplier's
    # liquidity component); `residual` carries whatever the FOC misses at the read
    # rather than smuggling it into a leg.
    dl = lambda a, b: np.log(np.asarray(a) / np.asarray(b))
    R_s, R_r = 1.0 + sim["rdep_D"], 1.0 + ref["rdep_D"]
    EOR_s = sim["E_Om_D"] * R_s
    EOR_r = ref["E_Om_D"] * R_r
    out = {
        "deposit_rate": -100.0 * dl(R_s, R_r),
        "continuation": 100.0 * dl(sim["E_payD_nodef"], ref["E_payD_nodef"]),
        "expected_loss": 100.0 * (dl(sim["E_payD"], sim["E_payD_nodef"])
                                  - dl(ref["E_payD"], ref["E_payD_nodef"])),
        "risk_premium": 100.0 * (dl(sim["E_Om_payD"], sim["E_Om_D"] * sim["E_payD"])
                                 - dl(ref["E_Om_payD"], ref["E_Om_D"] * ref["E_payD"])),
        "liquidity_premium": 100.0 * (dl(EOR_s, EOR_s + sim["lam_bD_mu_D"])
                                      - dl(EOR_r, EOR_r + ref["lam_bD_mu_D"])),
    }
    out["total"] = 100.0 * dl(sim["Q_bD"], ref["Q_bD"])
    out["residual"] = out["total"] - sum(out[k] for k, _ in BOND_CHANNELS
                                         if k != "residual")
    return out


CHANNELS = (("credit_spread", "Credit spread (bank leverage)"),
            ("deposit_rate", "Deposit rate"),
            ("capital", "Capital stock"),
            ("rel_price", "Relative price / terms of trade"),
            ("tfp", "TFP"),
            ("residual", "Labour-FOC residual"))


SOVEREIGN_LEGS = (("continuation", "continuation / duration"),
                  ("expected_loss", "expected default loss"),
                  ("risk_premium", "risk premium (Om-payoff cov)"),
                  ("liquidity_premium", "liquidity premium (the IC)"))


def sovereign_spread_legs(o, cal):
    # EXACT LEVEL DECOMPOSITION OF BOTH SOVEREIGN YIELDS AND OF THE D-F SPREAD.
    # decompose_bond_price splits the RESPONSE of q_D against a reference path; this
    # splits the LEVEL of each yield, which is what says how much of the spread a given
    # instrument could ever reach. Each bank's own FOC factors with no approximation into
    #     q = E[pay|nd]/R x E[pay]/E[pay|nd] x E[Om pay]/(E[Om]E[pay])
    #                     x E[Om]R/(E[Om]R + lam mu)
    #       = continuation x expected loss   x risk premium        x liquidity premium
    # and the yield attributed to a leg is what the yield would LOSE if that leg were
    # removed, y(q/leg) - y(q). Legs are NOT additive in yield (y is convex in q), so
    # they are reported as removals rather than as a partition.
    # EACH BOND IS DISCOUNTED AT ITS OWN DEPOSIT RATE. rdep_D and rdep_F coincide at the
    # symmetric steady state but diverge under the shock through deposit-UIP, and using
    # the D rate for the F bond puts that divergence into the F liquidity leg.
    RD, RF = 1.0 + o["rdep_D"], 1.0 + o["rdep_F"]
    D = {"continuation": o["E_payD_nodef"] / RD,
         "expected_loss": o["E_payD"] / o["E_payD_nodef"],
         "risk_premium": o["E_Om_payD"] / (o["E_Om_D"] * o["E_payD"]),
         "liquidity_premium": o["E_Om_D"] * RD / (o["E_Om_D"] * RD + o["lam_bD_mu_D"])}
    F = {"continuation": o["E_payF"] / RF,
         "expected_loss": 1.0,                      # the F sovereign never defaults
         "risk_premium": o["E_Om_payF"] / (o["E_Om_F"] * o["E_payF"]),
         "liquidity_premium": o["E_Om_F"] * RF / (o["E_Om_F"] * RF + o["lam_bF_mu_F"])}
    qD, qF = float(np.prod(list(D.values()))), float(np.prod(list(F.values())))
    dbD, dbF = cal["delta_b_D"], cal["delta_b_F"]

    def y(q, d):
        return 4e4 * d * (1.0 - q) / q             # HM perpetuity flow yield, ann. bp

    out = {"y_D": y(qD, dbD), "y_F": y(qF, dbF), "spread": y(qD, dbD) - y(qF, dbF),
           "q_D": qD, "q_F": qF,
           # how far the FITTED rules are from satisfying the FOC at this state: the
           # solve is exact only AT the collocation nodes and every IRF reads off-node,
           # so this is the honest error bar on every number in the table
           "foc_closure_D": qD / o["Q_bD"] - 1.0, "foc_closure_F": qF / o["Q_bF"] - 1.0}
    for k, _ in SOVEREIGN_LEGS:
        out[f"D_{k}"], out[f"F_{k}"] = D[k], F[k]
        out[f"y_D_{k}"] = y(qD / D[k], dbD) - y(qD, dbD)
        out[f"y_F_{k}"] = y(qF / F[k], dbF) - y(qF, dbF)
        out[f"spread_{k}"] = out[f"y_D_{k}"] - out[f"y_F_{k}"]
    return out


def decompose_output(sim, ref, cal):
    # EXACT FACTOR DECOMPOSITION OF THE OUTPUT DEVIATION (in % of the no-shock path).
    # The channels sum to `total` identically: the production function is used as an
    # identity, the labour FOC only to SPLIT d log N, and whatever the FOC misses at
    # the read is carried as `residual` rather than smuggled into a channel.
    a, phi, zeta = cal["alpha_D"], cal["frisch_D"], cal["zeta_wc_D"]
    theta = 1.0 / (1.0 / phi + a)

    dlZ = np.log(sim["Z_D"] / ref["Z_D"])
    dlK = np.log(sim["Kap_prod_D"] / ref["Kap_prod_D"])
    dlN = np.log(sim["N_D"] / ref["N_D"])
    dlP = np.log(sim["P_CES_D"] / ref["P_CES_D"])
    d_rd, d_sp = _wedge_shapley(sim["rdep_D"], sim["wedge_sp_D"],
                                ref["rdep_D"], ref["wedge_sp_D"], zeta)

    n_Z, n_K = theta * dlZ, theta * a * dlK
    n_rd, n_sp, n_P = -theta * d_rd, -theta * d_sp, -theta * dlP
    n_res = dlN - (n_Z + n_K + n_rd + n_sp + n_P)

    out = {
        "tfp": 100.0 * (dlZ + (1 - a) * n_Z),
        "capital": 100.0 * (a * dlK + (1 - a) * n_K),
        "credit_spread": 100.0 * (1 - a) * n_sp,
        "deposit_rate": 100.0 * (1 - a) * n_rd,
        "rel_price": 100.0 * (1 - a) * n_P,
        "residual": 100.0 * (1 - a) * n_res,
    }
    out["total"] = 100.0 * np.log(sim["Y_D"] / ref["Y_D"])
    out["hours"] = 100.0 * dlN
    return out


def active_channels(dec, tol=1e-4):
    # THE CHANNELS WORTH DRAWING (a channel flat at zero, e.g. TFP under a risk
    # shock, is dropped rather than shown as an empty band).
    return [(k, lab) for k, lab in CHANNELS
            if np.max(np.abs(dec[k])) > tol or k == "residual"]
