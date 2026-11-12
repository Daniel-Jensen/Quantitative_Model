# THE TPI EXPERIMENT: THE ECONOMY WITH AND WITHOUT THE EUROSYSTEM'S SPREAD CAP.
import numpy as np
from scipy.optimize import root

from global_projection.solver_recursive.state_grid import IS, IM
from global_projection.solver_recursive.decision_rules import SOLVE
from global_projection.solver_recursive.point_map import point_residuals
from global_projection.solver_recursive.recursive_experiment import (
    solve_recursive, dynamic_irf, stochastic_rest_point, read_at, read_exact, s_from_pd,
    S_REFINE)
from global_projection.solver_recursive.output_decomposition import decompose_bond_price

PD_SHOCK, T_IRF = 0.0198, 25
# bond-FOC legs read by decompose_bond_price
_BOND_KEYS = ("rdep_D", "E_Om_D", "E_payD", "E_payD_nodef", "E_Om_payD", "lam_bD_mu_D", "Q_bD")
# purchases for the price-impact report, as shares of the SS debt stock
IMPACT_M = (0.05, 0.10, 0.20, 0.30)


def _sov_bp(o, cal):
    # THE D-F SOVEREIGN SPREAD, ANNUALISED bp.
    yD = cal["delta_b_D"] * (1.0 - o["Q_bD"]) / o["Q_bD"]
    yF = cal["delta_b_F"] * (1.0 - o["Q_bF"]) / o["Q_bF"]
    return 4e4 * (yD - yF)


def rest_point_row(rules, cal, ss, sproc):
    # THE REST POINT, WHERE NOTHING IS BOUGHT.
    S = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    o = read_at(rules, cal, ss, sproc, S.copy())
    return dict(Q_bD=o["Q_bD"], sov_bp=_sov_bp(o, cal), mu=o["mu_D"], E_Om=o["E_Om_D"],
                alpha=o["alpha_D"], n=o["n_D"], Y=o["Y_D"], C=o["C_D"], I=o["I_D"],
                spread_bp=4e4 * cal["lambda_K_D"] * o["mu_D"] / max(o["alpha_D"], 1e-6),
                M=100 * S[IM] / cal["B_gov_D_ss"], m=100 * o["m_cb"] / cal["B_gov_D_ss"])


def price_impact(rules, cal, ss, sproc, S, ms=IMPACT_M):
    # THE PRICE IMPACT OF A FORCED ONE-OFF PURCHASE, RULES HELD FIXED.
    c = dict(cal, tpi_on=False)
    ngh = rules.n_gh or 5
    x0 = np.array([float(rules.eval(k, 0, np.atleast_2d(S))[0]) for k in SOLVE])[:-1]
    rows = []
    for m in (0.0,) + tuple(ms):
        # with the TPI off, the last unknown is the purchase itself
        f = lambda z: point_residuals(S, 0, np.append(z, m), rules, c, ss, sproc,
                                      n_gh=ngh)[0][:-1]
        sol = root(f, x0, method="hybr", tol=1e-12)
        _, o = point_residuals(S, 0, np.append(sol.x, m), rules, c, ss, sproc, n_gh=ngh)
        rows.append(dict(m=m, Q_bD=o["Q_bD"], sov_bp=_sov_bp(o, cal), mu=o["mu_D"],
                         b_DD=o["b_D_D_new"], resid=float(np.max(np.abs(sol.fun)))))
        x0 = sol.x
    return rows


def complementarity(rules, cal, ss, sproc):
    # THE TWO COMPLEMENTARITY PAIRS AT EVERY NO-DEFAULT NODE.
    rows = []
    for i in range(rules.grid.n):
        x = np.array([rules.vals[k][0][i] for k in SOLVE])
        _, o = point_residuals(rules.grid.points[i], 0, x, rules, cal, ss, sproc,
                               n_gh=rules.n_gh or 5)
        B, Q = cal["B_gov_D_ss"], ss["Q_bD_ss"]
        rows.append((o["m_cb"] / B, (o["Q_bD"] - o["Q_floor"]) / Q,
                     o["b_D_D_new"] / B, o["foc_D"]))
    r = np.array(rows)
    return dict(n=len(r), n_binding=int((r[:, 0] > 1e-6).sum()),
                n_corner=int((r[:, 2] < 1e-6).sum()),
                min_m=float(r[:, 0].min()), min_gap=float(r[:, 1].min()),
                max_m_gap=float(np.max(np.abs(r[:, 0] * r[:, 1]))),
                min_b=float(r[:, 2].min()), max_foc=float(r[:, 3].max()),
                max_b_foc=float(np.max(np.abs(r[:, 2] * r[:, 3]))))


def run(cal, ss, sproc, rules_off, base, s_refine=S_REFINE, mu_vec=None,
        pd_shock=PD_SHOCK):
    # SOLVE THE TPI ECONOMY AND COMPARE IT WITH THE NO-TPI ONE.
    rules_on = solve_recursive(dict(cal, tpi_on=True), ss, sproc, mu_vec=mu_vec,
                               s_refine=s_refine, base=base)
    return compare(cal, ss, sproc, rules_on, rules_off, pd_shock), rules_on


def compare(cal, ss, sproc, rules_on, rules_off, pd_shock=PD_SHOCK):
    # EVERY COMPARISON OF THE TWO SOLVED ECONOMIES.
    c_on, c_off = dict(cal, tpi_on=True), dict(cal, tpi_on=False)
    out = dict(cap_bp=cal["tpi_cap_bp"], key_D=cal["tpi_key_D"],
               cap_solved=getattr(rules_on, "tpi_cap_solved", None),
               solve_ok=bool(getattr(rules_on, "solve_ok", True)))
    out["rest_off"] = rest_point_row(rules_off, c_off, ss, sproc)
    out["rest_on"] = rest_point_row(rules_on, c_on, ss, sproc)
    # both paths are solved exactly each quarter: fitted reads break the corner
    out["irf_off"] = dynamic_irf(rules_off, c_off, ss, sproc, pd_shock=pd_shock, T=T_IRF,
                                 rest_verbose=False, exact=True)
    out["irf_on"] = dynamic_irf(rules_on, c_on, ss, sproc, pd_shock=pd_shock, T=T_IRF,
                                rest_verbose=False, exact=True)
    # bond-price legs on impact, each economy shocked from its own rest point
    reads = []
    for r, c in ((rules_on, c_on), (rules_off, c_off)):
        S = stochastic_rest_point(r, c, ss, sproc, verbose=False)
        S[IS] = s_from_pd(pd_shock)
        o = read_exact(r, c, ss, sproc, S)
        reads.append({k: np.array([o[k]]) for k in _BOND_KEYS})
    out["bond_legs"] = decompose_bond_price(reads[0], reads[1], cal)
    S_rest = stochastic_rest_point(rules_on, c_on, ss, sproc, verbose=False)
    S_shock = S_rest.copy(); S_shock[IS] = s_from_pd(pd_shock)
    out["impact_rest"] = price_impact(rules_on, c_on, ss, sproc, S_rest)
    out["impact_shock"] = price_impact(rules_on, c_on, ss, sproc, S_shock)
    out["kt"] = complementarity(rules_on, c_on, ss, sproc)
    return out
