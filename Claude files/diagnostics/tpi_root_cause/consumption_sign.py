# WHY D CONSUMPTION RISES IN THE NO-TPI CRISIS (2026-11-06, READ-ONLY; 12% BOX, p^d 9%): THE D HOUSEHOLD'S BUDGET
# C = W/P + inc - A AND ITS INCOME LEGS (WAGES, BANK AND FIRM DIVIDENDS, TAXES, THE WORKING-CAPITAL FLOW) ALONG THE PATH,
# AGAINST THE UNSHOCKED PATH, WITH THE D DEPOSIT RATE AND THE REAL RATE IN CONSUMPTION UNITS.
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

T = 9


def legs(rules, cal, ss, sproc, S, o):
    # THE BUDGET LEGS AT A SOLVED POINT, IN D CONSUMPTION UNITS.
    from global_projection.solver_recursive import point_map as PM
    v = PM._read_state_and_regime(S, 0, o["_x"], rules, cal, ss)
    for step in (PM._production, ):
        step(v, cal)
    PM._government(v, cal, ss); PM._balance_sheets(v, rules, cal)
    PM._continuation(v, rules, cal, sproc, rules.n_gh or 5, False); PM._discount_kernels(v, rules, cal)
    PM._incentive_constraint(v, cal); PM._bond_demands(v, cal); PM._tpi_rule(v, cal, ss)
    PM._wages_and_dividends(v, cal); PM._households(v, cal, ss)
    P = v.P_CES_D
    mcD = 1.0 - (v.Div_D - v.capprof_D - (cal["f_D"] * v.ng_D - cal["omega_ent_D"] * v.assets_D)) / v.Y_D
    bank = (cal["f_D"] * v.ng_D - cal["omega_ent_D"] * v.assets_D) / P
    wcf = (cal["zeta_wc_D"] * v.r_wc_D * v.w_D * v.N_D + v.L_wc_D) / P
    return dict(C=v.C_D, claim=v.W_D / P, A=v.A_D, wages=v.w_D * v.N_D / P, bank_div=bank,
                firm_div=(v.Div_D / P) - bank, taxes=-v.Tax_D / P, wc_flow=wcf, rdep=4e4 * v.rdep_D,
                P_CES=P, inc=v.inc_D)


if __name__ == "__main__":
    from crisis_shock import setup12
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    cal, ss, sproc, rules = setup12()
    S0 = np.load(HERE / "runs" / "crisis_rest_12pct.npy")
    paths = {}
    for lab, shock in (("ref", None), ("crisis", RE.s_from_pd(0.09))):
        S = S0.copy()
        if shock is not None:
            S[IS] = shock
        rows = []
        for t in range(T):
            o = RE.read_exact(rules, cal, ss, sproc, S)
            rows.append(legs(rules, cal, ss, sproc, S, o))
            S = rules.grid.clip(RE.advance(o, S, sproc))[0]
        paths[lab] = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    a, r = paths["crisis"], paths["ref"]
    Y = ss["ss_firm_D"]["Y_ss"]
    print("D HOUSEHOLD BUDGET, p^d 9% crisis minus unshocked path, % of quarterly D GDP (rates bp/yr), quarters 0..8")
    for k in ("C", "claim", "A", "inc", "wages", "bank_div", "firm_div", "taxes", "wc_flow"):
        print(f"  {k:9s} " + " ".join(f"{100*(a[k][t]-r[k][t])/Y:+8.3f}" for t in range(T)))
    print(f"  {'rdep':9s} " + " ".join(f"{a['rdep'][t]-r['rdep'][t]:+8.1f}" for t in range(T)))
    print(f"  {'CPI %':9s} " + " ".join(f"{100*(a['P_CES'][t]/r['P_CES'][t]-1):+8.3f}" for t in range(T)))
