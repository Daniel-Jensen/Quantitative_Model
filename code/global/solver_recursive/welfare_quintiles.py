# CONSUMPTION-EQUIVALENT WELFARE OF THE OMT/TPI BACKSTOP, BY INCOME QUINTILE.
# A DISTRIBUTIONAL OVERLAY on the projection solution, not a second equilibrium:
# the recursive solver closes with the rep-agent GHH deposit Euler, so the general
# equilibrium is taken as given and the incomplete-markets block (the SAME EGM +
# Young lottery used to pin beta at the steady state) is run along the resulting
# aggregate paths. Household saving therefore does NOT feed back into clearing --
# the standard sequence-space distributional accounting step, and the one caveat
# to state when the numbers are used.
#
# Construction. For each priced activation a: simulate the risk shock and the SAME
# rules with NO shock, and feed the household block the SS aggregates plus their
# DIFFERENCE (so a = anything with no shock returns exactly the steady state, and
# the projection's approximation drift cancels). Households then face
#   real return   1+r_t = (1+rdep_{t-1}) * P_CES_{t-1}/P_CES_t
#   income        y_t(e) = (w_t N_t / P_CES_t)*e + (Div_t - Tax_t)/P_CES_t
# with the aggregate GHH disutility v(N_t) common across households (labour income
# is exogenous to the household here, exactly as in the steady-state block).
# Welfare is the date-0 value V_0(a,e) under the shock; the consumption-equivalent
# cost of the shock is lambda(a,e) from scaling the GHH composite c - v(N) in every
# period and state, and the OMT GAIN is lambda(activation) - lambda(no backstop).
# Quintiles are cut on SS TOTAL income -- labour + lump-sum + asset income -- over
# the joint (a,e) stationary distribution, with the boundary cell's mass SPLIT so
# every quintile carries exactly 20% of households.
import numpy as np

from blocks.firms import markup_ss
from blocks.household import solve_backward_transition
from blocks.distribution import forward_iterate
from blocks.trade import ces_price


def ss_household_inputs(ss, cal):
    # THE STEADY-STATE HOUSEHOLD ENVIRONMENT, REBUILT EXACTLY AS steady_state.py DID.
    # Must match that block line for line: c_D_ss / D_D_ss / beta_D_ss were solved
    # against THIS y_e, so any other anchor would make V_ss inconsistent with them.
    p_ss = ss["p_ss"]
    P_CES = float(ces_price(np.array([p_ss]), cal, "D")[0])
    mc = markup_ss(cal, "D")
    w_ss = ss["ss_firm_D"]["w_ss"]
    r_wc_ss = cal["r_dep_D_target"] + cal["credit_spread_target_D"]
    wc_ss = cal["zeta_wc_D"] * r_wc_ss * w_ss
    Div_ss = (1 - mc) * ss["ss_firm_D"]["Y_ss"] + ss["ss_bank_D"]["div_ss"] + wc_ss
    Tax_ss = ss["gs_D"]["Tax_ss"]
    vN_ss = cal["chi_D"] / (1 + 1 / cal["frisch_D"])          # GHH v(N) at N_ss = 1
    return dict(P_CES=P_CES, wage=w_ss / P_CES, lump=(Div_ss - Tax_ss) / P_CES,
                r=cal["r_dep_D_target"], vN=vN_ss)


def _u(x, sigma):
    # PERIOD UTILITY OVER THE GHH COMPOSITE x = c - v(N).
    x = np.maximum(x, 1e-11)
    return np.log(x) if abs(sigma - 1.0) < 1e-12 else x ** (1 - sigma) / (1 - sigma)


def _expected_value(V_next, a_pol, a_grid, Pi):
    # E_e'[V_next(a'(a,e), e')] UNDER THE PRODUCTIVITY TRANSITION.
    n_a, n_e = a_pol.shape
    M = np.empty((n_a, n_e, n_e))
    for ep in range(n_e):
        M[:, :, ep] = np.interp(a_pol, a_grid, V_next[:, ep])
    return np.einsum("aek,ek->ae", M, Pi)


def steady_state_value(ss, cal, env, tol=1e-10, maxiter=200_000):
    # THE STEADY-STATE VALUE FUNCTION ON (a, e), ITERATED ON THE SS POLICIES.
    a_grid, Pi, e = ss["a_grid_D"], ss["Pi_D"], ss["e_D"]
    c_ss, beta, sigma = ss["c_D_ss"], ss["beta_D_ss"], cal["sigma_D"]
    y_e = env["wage"] * e + env["lump"]
    a_pol = np.maximum((1 + env["r"]) * a_grid[:, None] + y_e[None, :] - c_ss,
                       cal["a_min_D"])
    u = _u(c_ss - env["vN"], sigma)
    V = u / (1.0 - beta)
    for _ in range(maxiter):
        V_new = u + beta * _expected_value(V, a_pol, a_grid, Pi)
        if np.max(np.abs(V_new - V)) < tol:
            return V_new, a_pol
        V = V_new
    raise RuntimeError("steady-state value iteration did not converge")


def aggregate_inputs(sim, ref, ss, cal, env):
    # SS AGGREGATES PLUS THE SIMULATED DEVIATION -> THE HOUSEHOLD'S (r, y, v(N)) PATHS.
    # Deviations, not levels: the no-shock path is then the steady state EXACTLY,
    # whatever small level offset the projection carries at its own rest point.
    T = len(sim["Y_D"])
    wage = env["wage"] + (sim["w_D"] * sim["N_D"] / sim["P_CES_D"]
                          - ref["w_D"] * ref["N_D"] / ref["P_CES_D"])
    lump = env["lump"] + ((sim["Div_D"] - sim["Tax_D"]) / sim["P_CES_D"]
                          - (ref["Div_D"] - ref["Tax_D"]) / ref["P_CES_D"])
    # the return earned at t was locked at t-1 (predetermined deposit rate), and is
    # deflated by the CES basket's own move between t-1 and t
    rd_lag = np.concatenate(([cal["r_dep_D_target"]], sim["rdep_D"][:-1]))
    rd_lag_r = np.concatenate(([cal["r_dep_D_target"]], ref["rdep_D"][:-1]))
    P_lag = np.concatenate(([env["P_CES"]], sim["P_CES_D"][:-1]))
    P_lag_r = np.concatenate(([env["P_CES"]], ref["P_CES_D"][:-1]))
    r = env["r"] + ((1 + rd_lag) * P_lag / sim["P_CES_D"]
                    - (1 + rd_lag_r) * P_lag_r / ref["P_CES_D"])
    N = sim["N_D"] / ref["N_D"]                                # N_ss = 1
    vN = cal["chi_D"] * N ** (1 + 1 / cal["frisch_D"]) / (1 + 1 / cal["frisch_D"])
    y = wage[:, None] * ss["e_D"][None, :] + lump[:, None]
    r_full = np.concatenate((r, [env["r"]]))                   # EGM needs r_{T+1}
    drift = float(max(np.max(np.abs(y[-20:] - (env["wage"] * ss["e_D"] + env["lump"]))),
                      np.max(np.abs(r[-20:] - env["r"]))))
    return dict(r=r_full, y=y, vN=vN, T=T, drift=drift)


def transition_welfare(inp, ss, cal, V_ss):
    # DATE-0 VALUE ON (a, e) ALONG THE TRANSITION, PLUS THE PATH OF DISTRIBUTIONS.
    # Backward: the same EGM the SS block used, terminal policy c_ss; then the value
    # recursion on the realised consumption path with terminal V_ss.
    a_grid, Pi, beta = ss["a_grid_D"], ss["Pi_D"], ss["beta_D_ss"]
    c_path, a_pol_path = solve_backward_transition(
        a_grid, Pi, inp["r"], inp["y"], ss["c_D_ss"], beta, cal["sigma_D"],
        cal["a_min_D"], vN_path=inp["vN"], use_fast=cal["use_numba"])
    V = V_ss
    for t in range(inp["T"] - 1, -1, -1):
        V = (_u(c_path[t] - inp["vN"][t], cal["sigma_D"])
             + beta * _expected_value(V, a_pol_path[t], a_grid, Pi))
    return V, c_path, a_pol_path


def cev(V_shock, V_base, ss, cal):
    # CONSUMPTION-EQUIVALENT DEVIATION: THE PERMANENT SCALING OF c - v(N) THAT
    # MAKES THE SHOCK PATH AS GOOD AS THE REFERENCE (negative = welfare cost).
    beta, sigma = ss["beta_D_ss"], cal["sigma_D"]
    if abs(sigma - 1.0) < 1e-12:
        return np.exp((1.0 - beta) * (V_shock - V_base)) - 1.0
    return (V_shock / V_base) ** (1.0 / (1.0 - sigma)) - 1.0


def income_quintile_weights(ss, cal, env, n_q=5):
    # QUINTILE WEIGHT MASKS ON (a, e), CUT ON SS TOTAL INCOME, EXACT 20% MASS EACH.
    # Cells are ranked by income and filled into buckets; the cell straddling a
    # boundary has its MASS SPLIT, so a point mass at the borrowing constraint
    # (which alone can exceed a fifth of the population) cannot distort the cut.
    a_grid, e, D = ss["a_grid_D"], ss["e_D"], ss["D_D_ss"]
    inc = (env["wage"] * e[None, :] + env["lump"]
           + env["r"] * a_grid[:, None])
    order = np.argsort(inc.ravel(), kind="stable")
    mass = D.ravel()[order]
    total = mass.sum()
    W = np.zeros((n_q, inc.size))
    edge = total / n_q
    q, filled = 0, 0.0
    for j, m in enumerate(mass):
        while m > 1e-15:
            if q >= n_q - 1:                       # last bucket absorbs the remainder
                W[q, order[j]] += m
                filled += m
                break
            room = max(edge * (q + 1) - filled, 0.0)
            take = min(m, room)
            W[q, order[j]] += take
            filled += take
            m -= take
            if filled >= edge * (q + 1) - 1e-15:
                q += 1
    W = W.reshape(n_q, *D.shape)
    inc_q = np.array([float(np.sum(W[k] * inc) / np.sum(W[k])) for k in range(n_q)])
    return W, inc_q, inc


def quintile_cev(V_shock, V_base, W, ss, cal):
    # PER-QUINTILE CEV: THE UTILITARIAN GROUP AGGREGATE (mass-weighted mean of dV,
    # then converted), which is the CEV of the quintile as a single welfare unit.
    dV = V_shock - V_base
    beta, sigma = ss["beta_D_ss"], cal["sigma_D"]
    out = np.empty(W.shape[0])
    for k in range(W.shape[0]):
        w = W[k] / np.sum(W[k])
        if abs(sigma - 1.0) < 1e-12:
            out[k] = np.exp((1.0 - beta) * float(np.sum(w * dV))) - 1.0
        else:
            out[k] = (float(np.sum(w * V_shock)) / float(np.sum(w * V_base))) \
                ** (1.0 / (1.0 - sigma)) - 1.0
    return 100.0 * out


def cohort_consumption(ss, c_path, a_pol_path, W, T):
    # PER-QUINTILE MEAN CONSUMPTION OVER T PERIODS, TRACKING THE DATE-0 COHORTS.
    # W[k] already IS the quintile's mass per (a, e) cell; forwarding it with the
    # same lottery operator follows that cohort as it moves through the asset grid
    # (the cohort is cut once, at date 0, and never re-cut). The cohort drifts even
    # with NO shock, so this is only meaningful against the no-shock cohort path.
    a_grid, Pi = ss["a_grid_D"], ss["Pi_D"]
    n_q = W.shape[0]
    C = np.empty((n_q, T))
    Dq = [W[k].copy() for k in range(n_q)]
    for t in range(T):
        for k in range(n_q):
            C[k, t] = float(np.sum(Dq[k] * c_path[t]) / np.sum(Dq[k]))
            Dq[k] = forward_iterate(Dq[k], a_pol_path[t], a_grid, Pi)
    return C


def ss_cohort_consumption(ss, a_pol_ss, W, T):
    # THE SAME COHORT PATH WITH NO SHOCK: SS POLICIES HELD FIXED FOR T PERIODS.
    c_path = np.broadcast_to(ss["c_D_ss"], (T,) + ss["c_D_ss"].shape)
    a_path = np.broadcast_to(a_pol_ss, (T,) + a_pol_ss.shape)
    return cohort_consumption(ss, c_path, a_path, W, T)
