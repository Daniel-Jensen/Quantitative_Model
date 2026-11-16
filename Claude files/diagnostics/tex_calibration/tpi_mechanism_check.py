# WHY THE TPI CANNOT HOLD THE CAP (2026-10-30): AT EACH RISK LEVEL (EXACT READS, OTHER STATES AT THE REST POINT) THE D
# CONSTRAINT'S TIGHTNESS, THE BOND'S EXCESS RETURN SPLIT INTO EXPECTED DEFAULT LOSS + RISK PREMIUM + LIQUIDITY PREMIUM, AND
# THE PRICE EFFECT OF A ONE-OFF PURCHASE OF 5% OF THE D DEBT STOCK (RULES HELD FIXED), WITH THE PART THAT GOES THROUGH mu.
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
import numpy as np
from scipy.optimize import root


def main():
    import run
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals, _tpi_floor
    from global_projection.solver_recursive.output_decomposition import bond_premia
    cal = get_calibration(); cal.update(nw_floor_frac=run.NW_FLOOR, jacobian=run.JACOBIAN)
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    r = pickle.load(open(ROOT / "results/GLOBAL/data/rules_risk.pkl", "rb"))
    c = dict(cal, tpi_on=False)
    S0 = RE.stochastic_rest_point(r, c, ss, sp, verbose=False)
    B, ngh = cal["B_gov_D_ss"], r.n_gh or 5
    lam = cal["lambda_bD_D"]

    def solve_at(S, m, x0):
        # THE PERIOD MAP CLEARED AT S WITH A ONE-OFF PURCHASE m (TPI OFF: THE LAST UNKNOWN IS THE PURCHASE).
        f = lambda z: point_residuals(S, 0, np.append(z, m), r, c, ss, sp, n_gh=ngh)[0][:-1]
        sol = root(f, x0, method="hybr", tol=1e-12)
        return sol.x, point_residuals(S, 0, np.append(sol.x, m), r, c, ss, sp, n_gh=ngh)[1], float(np.max(np.abs(sol.fun)))

    m5 = 0.05 * B
    print(f"{'p^d':>6s} {'mu_D':>7s} {'slack':>7s} {'E[Om]':>6s} {'n/n_ss':>6s} | {'excess':>7s} {'=loss':>7s} {'+risk':>7s} "
          f"{'+liq':>6s} (bp/yr) | {'Q_bD':>6s} {'floor200':>8s} | buy 5% of B: {'dQ/Q%':>6s} {'d mu':>7s} {'d liq bp':>8s} "
          f"{'liq relief max %Q':>17s}")
    for pd in (None, 0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.12):
        S = S0.copy()
        if pd:
            S[IS] = RE.s_from_pd(pd)
        x0 = np.array([float(r.eval_read(k, 0, S[None])[0]) for k in SOLVE])[:-1]
        x, o, res0 = solve_at(S, 0.0, x0)
        _, o5, res5 = solve_at(S, m5, x)
        ex, liq, risk = (float(v) for v in bond_premia(o))
        # the actuarial expected default loss rate: (E[pay | no default] - E[pay]) / Q, annualised
        loss = 4e4 * (o["E_payD_nodef"] - o["E_payD"]) / o["Q_bD"]
        _, liq5, _ = (float(v) for v in bond_premia(o5))
        slack = o["E_Om_D"] * (1 + o["rdep_D"]) * o["n_D"] / max(o.get("lev_D", 0) or 1e-12, 1e-12) if "lev_D" in o else np.nan
        floor = _tpi_floor(o["Q_bF"], dict(cal, tpi_cap_bp=200.0))
        # the most a full relief of today's constraint could lift the price: lambda*mu/(E[Om]R)
        relief = 100 * lam * o["mu_D"] / (o["E_Om_D"] * (1 + o["rdep_D"]))
        print(f"{'rest' if pd is None else f'{100*pd:.0f}%':>6s} {o['mu_D']:7.4f} {slack:7.3f} {o['E_Om_D']:6.3f} "
              f"{o['n_D']/ss['ss_bank_D']['n_ss']:6.3f} | {ex + loss:7.1f} {loss:7.1f} {risk:7.1f} {liq:6.1f}         | "
              f"{o['Q_bD']:6.3f} {floor:8.3f} |              {100*(o5['Q_bD']/o['Q_bD']-1):+6.2f} {o5['mu_D']-o['mu_D']:+7.4f} "
              f"{liq5-liq:+8.1f} {relief:17.2f}   (resid {max(res0, res5):.0e})")


if __name__ == "__main__":
    main()
