# THE REST POINT AS A FIXED POINT OF THE EXACT DYNAMICS, AND BOTH IRFS FROM IT (2026-11-06, READ-ONLY): run.py's rest
# point iterates the fitted reads, and in the TPI economy the fitted x_cb sets the bond price.
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import root

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

# the states the rest point solves for (s and Z_D are exogenous)
FREE = (0, 1, 2, 3, 4, 5, 6, 7, 10, 11)


def exact_rest(rules, cal, ss, sproc, S_start):
    # SOLVE S = advance(read_exact(S)) FOR THE ENDOGENOUS STATES.
    from global_projection.solver_recursive import recursive_experiment as RE
    idx = np.array(FREE)
    scale = np.maximum(np.abs(S_start[idx]), 1e-2)

    def f(z):
        S = S_start.copy(); S[idx] = z * scale
        o = RE.read_exact(rules, cal, ss, sproc, S)
        return (RE.advance(o, S, sproc)[idx] - S[idx]) / scale

    sol = root(f, S_start[idx] / scale, method="hybr", tol=1e-13)
    S = S_start.copy(); S[idx] = sol.x * scale
    return S, float(np.max(np.abs(sol.fun)))


def row(o, cal):
    # THE REST-POINT OBJECTS THE REPORT QUOTES.
    yD = cal["delta_b_D"] * (1 - o["Q_bD"]) / o["Q_bD"]
    yF = cal["delta_b_F"] * (1 - o["Q_bF"]) / o["Q_bF"]
    return dict(sov_bp=4e4 * (yD - yF), Q_bD=o["Q_bD"], n_D=o["n_D"], E_Om=o["E_Om_D"], C_D=o["C_D"],
                Y_D=o["Y_D"], I_D=o["I_D"], B=o["B_D"], Bp=o["Bp_D"], rdep=4e4 * o["rdep_D"])


if __name__ == "__main__":
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import STATE_NAMES
    from post_impact import setup
    cal, ss, sproc, r_off, r_on, c_off, c_on, run = setup()

    rests = {}
    for lab, rules, c in (("off", r_off, c_off), ("on", r_on, c_on)):
        S_fit = RE.stochastic_rest_point(rules, c, ss, sproc, verbose=False)
        of = RE.read_at(rules, c, ss, sproc, S_fit.copy())
        ox = RE.read_exact(rules, c, ss, sproc, S_fit.copy())
        print(f"\n[{lab}] AT THE FITTED-READ REST POINT: fitted read resid {of['_resid']:.1e}")
        for k, v in row(of, c).items():
            print(f"   {k:7s} fitted {v: .6f}   exact {row(ox, c)[k]: .6f}")
        S_ex, err = exact_rest(rules, c, ss, sproc, S_fit)
        rests[lab] = S_ex
        print(f"[{lab}] EXACT REST POINT (|S - g(S)| = {err:.1e}); state change vs the fitted-read rest point:")
        for i, nm in enumerate(STATE_NAMES):
            if abs(S_fit[i]) > 1e-9:
                print(f"   {nm:6s} {S_fit[i]: .6e} -> {S_ex[i]: .6e}  ({100*(S_ex[i]/S_fit[i]-1):+.4f}%)")
            else:
                print(f"   {nm:6s} {S_fit[i]: .6e} -> {S_ex[i]: .6e}")
        o = RE.read_exact(rules, c, ss, sproc, S_ex.copy())
        print(f"   at the exact rest point: " + ", ".join(f"{k} {v:.6f}" for k, v in row(o, c).items()))

    print("\nREST POINTS, EXACT, TPI vs NO TPI")
    oa = row(RE.read_exact(r_off, c_off, ss, sproc, rests["off"].copy()), cal)
    ob = row(RE.read_exact(r_on, c_on, ss, sproc, rests["on"].copy()), cal)
    for k in oa:
        print(f"   {k:7s} off {oa[k]: .6f}  on {ob[k]: .6f}  diff {ob[k]-oa[k]:+.6f}"
              + (f" ({100*(ob[k]/oa[k]-1):+.3f}%)" if abs(oa[k]) > 1e-9 else ""))

    # both IRFs from the exact rest points (dynamic_irf reads the cached rest point)
    irf = {}
    for lab, rules, c in (("off", r_off, c_off), ("on", r_on, c_on)):
        rules._rest_point = rests[lab].copy()
        irf[lab] = RE.dynamic_irf(rules, c, ss, sproc, pd_shock=run.RISK_SHOCK_PD, T=25, rest_verbose=False,
                                  exact=True)
    np.save(Path(__file__).resolve().parent / "runs" / "exact_rest_irf.npy", irf, allow_pickle=True)
    print("\nIRFS FROM THE EXACT REST POINTS")
    for k in ("sov_bp", "Q_bD", "n", "E_Om", "rdep", "I", "K", "Y", "C", "N", "m_cb", "M_cb", "Pi_cb", "n_F", "C_F"):
        for lab in ("off", "on"):
            print(f"  {k:7s} {lab:3s} " + " ".join(f"{v:8.3f}" for v in irf[lab][k][[0, 1, 2, 4, 8, 12, 16, 24]]))
