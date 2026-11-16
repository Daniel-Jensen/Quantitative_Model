# THE HEADLINE SHOCK WITH AND WITHOUT A REALIZED DEFAULT, IN THE TPI AND NO-TPI ECONOMIES OF ONE RUN.
#   WT=dff3 python3 default_path.py <tag> [default quarter, default 1]
# Every path starts at the economy's own rest point and is differenced against its unshocked path. The default
# path takes the regime d = 1 (the haircut) in one quarter only, as the model's i.i.d. default regime does.
import json
import os
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / os.environ.get("WT", "dff3")))
sys.path.insert(0, str(HERE))

import numpy as np

KEYS = ("Y_D", "C_D", "I_D", "n_D", "Q_bD", "B_D", "Bp_D", "Tax_D", "Tax_F", "Pi_cb", "m_cb", "M_cb_new",
        "C_F", "Y_F", "n_F", "rdep_D", "r_wc_D", "mu_D", "p")


def read_d(rules, cal, ss, sproc, S, d):
    # CLEAR THE PERIOD MAP EXACTLY AT S IN REGIME d.
    from scipy.optimize import root
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    ngh = rules.n_gh or 5
    x0 = np.array([float(rules.eval(k, d, np.atleast_2d(S))[0]) for k in SOLVE])

    def f(x):
        try:
            return point_residuals(S, d, x, rules, cal, ss, sproc, n_gh=ngh)[0]
        except (ValueError, RuntimeError, ArithmeticError):
            return np.full(len(SOLVE), 10.0)

    sol = root(f, x0, method="hybr", tol=1e-13)
    _, o = point_residuals(S, d, sol.x, rules, cal, ss, sproc, n_gh=ngh)
    o["_x"] = sol.x
    o["_resid"] = float(np.max(np.abs(sol.fun)))
    return o


def path(rules, cal, ss, sproc, S0, T, s_shock=None, d_at=None):
    # ONE SIMULATED PATH; s_shock SETS s AT q0, d_at THE QUARTER THAT DEFAULTS.
    from global_projection.solver_recursive.recursive_experiment import advance
    from global_projection.solver_recursive.state_grid import IS
    S, rows, worst = S0.copy(), [], 0.0
    if s_shock is not None:
        S[IS] = s_shock
    for t in range(T):
        o = read_d(rules, cal, ss, sproc, S, 1 if t == d_at else 0)
        worst = max(worst, o["_resid"])
        rows.append({k: float(o[k]) for k in KEYS})
        S = advance(o, S, sproc, rules.grid)
    return {k: np.array([r[k] for r in rows]) for k in KEYS}, worst


def run(tag, d_at=1, T=12, pd=0.02):
    # BOTH ECONOMIES OF ONE SAVED RUN: SHOCK WITHOUT AND WITH DEFAULT.
    import argparse
    import tpi_national as tn
    from global_projection.solver_recursive.recursive_experiment import (
        stochastic_rest_point, s_from_pd)
    meta = json.loads((HERE / "tpi_runs" / f"{tag}.json").read_text())["args"]
    a = argparse.Namespace(**{k: meta.get(k) for k in ("market", "n_jobs", "cap", "lamZ", "phi", "keyD",
                                                       "theta", "warranted", "keyLoss")})
    cal, ss, sproc = tn.setup(a)
    rr = pickle.load(open(HERE / "tpi_runs" / f"{tag}_rules.pkl", "rb"))
    cal["tpi_cap_bp"] = float(getattr(rr["on"], "tpi_cap_solved", None) or cal["tpi_cap_bp"])
    out = {}
    for name in ("off", "on"):
        c = dict(cal, tpi_on=(name == "on"))
        S0 = stochastic_rest_point(rr[name], c, ss, sproc, verbose=False)
        ref, w0 = path(rr[name], c, ss, sproc, S0, T)
        nod, w1 = path(rr[name], c, ss, sproc, S0, T, s_from_pd(pd))
        dfl, w2 = path(rr[name], c, ss, sproc, S0, T, s_from_pd(pd), d_at)
        Y = ss["ss_firm_D"]["Y_ss"]
        dev = {}
        for lab, pth in (("nodefault", nod), ("default", dfl)):
            dev[lab] = {
                "Y": 100 * (pth["Y_D"] / ref["Y_D"] - 1), "C": 100 * (pth["C_D"] / ref["C_D"] - 1),
                "I": 100 * (pth["I_D"] / ref["I_D"] - 1), "n": 100 * (pth["n_D"] / ref["n_D"] - 1),
                "C_F": 100 * (pth["C_F"] / ref["C_F"] - 1), "n_F": 100 * (pth["n_F"] / ref["n_F"] - 1),
                # debt and fiscal flows in % of quarterly D GDP
                "B": 100 * (pth["Bp_D"] - ref["Bp_D"]) / Y,
                "Tax_D": 100 * (pth["Tax_D"] - ref["Tax_D"]) / Y,
                "Tax_F_in_D_GDP": 100 * (pth["Tax_F"] - ref["Tax_F"]) * meta_size(cal) * pth["p"] / Y,
                "Pi_cb": 100 * (pth["Pi_cb"] - ref["Pi_cb"]) / Y,
                "M": 100 * pth["M_cb_new"] / cal["B_gov_D_ss"]}
        out[name] = dict(dev=dev, resid=max(w0, w1, w2))
    return out


def meta_size(cal):
    # F'S MASS RELATIVE TO D, TO PUT F'S TAXES IN D-GDP UNITS.
    return cal["size_F"] / cal["size_D"]


def main():
    # PRINT AND SAVE THE COMPARISON.
    tag = sys.argv[1]
    d_at = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    out = run(tag, d_at)
    for name in ("off", "on"):
        print(f"\n  [{tag}] TPI {name}: worst residual {out[name]['resid']:.1e}")
        for lab in ("nodefault", "default"):
            d = out[name]["dev"][lab]
            print(f"   {lab:9s} " + "  ".join(
                f"{k} " + " ".join(f"{x:+.2f}" for x in d[k][:6]) for k in ("Y", "C", "n", "B", "Tax_D",
                                                                            "Tax_F_in_D_GDP", "Pi_cb", "M")))
    clean = {n: {"resid": out[n]["resid"], "dev": {l: {k: v.tolist() for k, v in out[n]["dev"][l].items()}
                                                   for l in out[n]["dev"]}} for n in out}
    (HERE / "tpi_runs" / f"{tag}_default_q{d_at}.json").write_text(json.dumps(clean, indent=1))


if __name__ == "__main__":
    main()
