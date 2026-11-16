# IMPACT WITH NEXT QUARTER SOLVED EXACTLY (2026-11-06, READ-ONLY): EACH NEXT STATE (5 GH NODES x 2 REGIMES) IS CLEARED
# EXACTLY, WITH THE FITTED RULES ONLY FROM t+2 ON. HOW MUCH OF THE IMPACT RESPONSE IS THE t+1 RULE ERROR.
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import root

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))


class ExactNext:
    # THE FITTED RULES FOR TODAY'S FROZEN OBJECTS; EXACT SOLVES AT THE NEXT STATES.
    def __init__(self, rules, cal, ss, sproc):
        self.r, self.cal, self.ss, self.sproc = rules, cal, ss, sproc
        self.reg, self.n_regimes, self.n_gh, self.grid = rules.reg, rules.n_regimes, rules.n_gh, rules.grid
        self.home = self.post = None

    def eval(self, name, d, x):
        return self.r.eval(name, d, x)

    def eval_cont(self, j, Sn):
        from regime_anatomy import solve_regime
        from global_projection.solver_recursive.decision_rules import SOLVE, STORE_RULES
        Sc = self.r.grid.clip(np.real(np.atleast_2d(Sn)))
        out = {k: np.empty(len(Sc)) for k in STORE_RULES}
        for i, S in enumerate(Sc):
            o = solve_regime(self.r, self.cal, self.ss, self.sproc, S, j)
            for n, k in enumerate(SOLVE):
                out[k][i] = o["_x"][n]
            out["Q_bD"][i] = o["Q_bD"]
            for k in ("alpha_D", "alpha_F", "C_D", "C_F", "r_wc_D", "r_wc_F"):
                out[k][i] = o[k]
        return out


def solve_with(cont, rules, cal, ss, sproc, S, x0, m=None):
    # CLEAR THE PERIOD MAP AT S WITH CONTINUATION cont; m FIXES THE PURCHASE (TPI OFF).
    from global_projection.solver_recursive.point_map import point_residuals
    ngh = rules.n_gh or 5
    if m is None:
        f = lambda x: point_residuals(S, 0, x, cont, cal, ss, sproc, n_gh=ngh)[0]
        sol = root(f, x0, method="hybr", tol=1e-12)
        x = sol.x
    else:
        f = lambda z: point_residuals(S, 0, np.append(z, m), cont, cal, ss, sproc, n_gh=ngh)[0][:-1]
        sol = root(f, x0[:-1], method="hybr", tol=1e-12)
        x = np.append(sol.x, m)
    _, o = point_residuals(S, 0, x, cont, cal, ss, sproc, n_gh=ngh)
    o["_x"] = x
    o["_resid"] = float(np.max(np.abs(sol.fun)))
    return o


def summary(rules, cal, ss, S, o):
    # THE IMPACT OBJECTS.
    from one_off_purchase import capital_price
    yD = cal["delta_b_D"] * (1 - o["Q_bD"]) / o["Q_bD"]
    yF = cal["delta_b_F"] * (1 - o["Q_bF"]) / o["Q_bF"]
    return dict(sov_bp=4e4 * (yD - yF), Q_bD=o["Q_bD"], Q_K=capital_price(rules, cal, ss, S, o), I_D=o["I_D"],
                n_D=o["n_D"], C_D=o["C_D"], N_D=o["_x"][0], Y_D=o["Y_D"], rdep_D=4e4 * o["rdep_D"], E_Om=o["E_Om_D"],
                p=o["p"])


if __name__ == "__main__":
    from post_impact import setup
    from one_off_purchase import read_forced
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    cal, ss, sproc, r_off, r_on, c_off, c_on, run = setup()
    rests = np.load(HERE / "runs" / "exact_rest_points.npy", allow_pickle=True).item()
    s_sh = RE.s_from_pd(run.RISK_SHOCK_PD)
    res = {}
    for lab, rules, c, m in (("A no TPI", r_off, c_off, None), ("B TPI", r_on, c_on, None),
                             ("C no TPI + 41% buy", r_off, c_off, 0.4071)):
        S_rest = rests["off" if rules is r_off else "on"].copy()
        S_imp = S_rest.copy(); S_imp[IS] = s_sh
        xn = ExactNext(rules, c, ss, sproc)
        row = {}
        for st, S, mm in (("rest", S_rest, None), ("impact", S_imp, m)):
            t0 = time.time()
            of = RE.read_exact(rules, c, ss, sproc, S.copy()) if mm is None else read_forced(rules, c, ss, sproc, S, mm)
            ox = solve_with(xn, rules, c, ss, sproc, S.copy(), of["_x"], mm)
            row[st] = (summary(rules, c, ss, S, of), summary(rules, c, ss, S, ox))
            print(f"  {lab} {st}: fitted-continuation resid {of['_resid']:.0e}, exact-next resid {ox['_resid']:.0e}"
                  f" ({time.time()-t0:.0f}s)", flush=True)
        res[lab] = row
    np.save(HERE / "runs" / "one_step_exact.npy", res, allow_pickle=True)
    print("\nIMPACT RESPONSE (impact vs rest, same continuation method): fitted t+1 | exact t+1")
    for k in ("sov_bp", "Q_bD", "Q_K", "I_D", "n_D", "C_D", "N_D", "Y_D", "rdep_D", "E_Om", "p"):
        line = []
        for lab, row in res.items():
            vals = []
            for i in (0, 1):
                a, b = row["rest"][i][k], row["impact"][i][k]
                vals.append(b - a if k in ("sov_bp", "rdep_D") else 100 * (b / a - 1))
            line.append(f"{lab}: {vals[0]:+8.3f} | {vals[1]:+8.3f}")
        print(f"  {k:7s} " + "    ".join(line))
    print("\nREST LEVELS: fitted t+1 | exact t+1")
    for k in ("sov_bp", "Q_bD", "n_D", "C_D", "I_D", "E_Om"):
        print(f"  {k:7s} " + "    ".join(f"{lab}: {row['rest'][0][k]:.5f} | {row['rest'][1][k]:.5f}" for lab, row in res.items()))
