# OUT-OF-SAMPLE ACCURACY OF A SOLVED RULE SET (Bocola's generate_euler_errors.m).
# The regression tests check the grid, the kernels, the SS identities and the pi=0
# nesting -- all necessary, none of them an accuracy measure on the ergodic set. This
# module answers the question a referee asks: how wrong is the solution WHERE THE MODEL
# ACTUALLY SPENDS ITS TIME. It simulates the solved rules forward, and at every visited
# state evaluates the seven equilibrium conditions at the RULE-IMPLIED policy (no
# re-solving), reporting the residual distribution in log10 units per equation.
# It also reports how often the simulated path leaves the collocation box, which is the
# invariance check the box design turns on.
import numpy as np

from solver_recursive.point_map import point_residuals, SOLVE7
from solver_recursive.state_grid import default_prob
from solver_recursive.decision_rules import DERIVED, to_fit

# RESIDUAL NAMES, in point_map order. Single-sourced so the report cannot fall out of
# step with the system the way the hard-wired 7 did when it grew to 11.
# THE FULL 19-EQUATION SYSTEM the global collocation solve roots (collocation.py):
# the 13 residuals point_map returns, plus Bocola's identity residual log(guess/implied)
# for the six objects that used to be READ OFF the recursions against a frozen
# continuation (alpha, C, r_wc per country). Those are unknowns now, so leaving them out
# of the accuracy report would have measured only part of the system.
_EQN = ("cap_D", "cap_F", "lab_D", "lab_F", "euler_D", "uip", "goods_D",
        "bondD_D", "bondD_F", "euler_F", "dep_clear", "bondF_F", "bondF_D"
        ) + tuple(f"id_{k}" for k in DERIVED)
from solver_recursive.state_grid import (IK_D, IK_F, IP_D, IP_F, IBDD,
                                          IBDF, IBFD, IV, IS, IZ, STATE_NAMES)

EQ_NAMES = ("cap_D", "cap_F", "lab_D", "lab_F", "euler_D", "uip", "goods_D")


def simulate(rules, cal, ss, sproc, S0, T=2000, burn=200, seed=0, no_default=False):
    # SIMULATE THE SOLVED RULES FORWARD, RETURNING THE VISITED STATES.
    # s follows its own AR(1) with drawn innovations; the endogenous stocks follow the
    # period map's own end-of-period values, so the path is the model's, not the box's.
    # States are NOT clipped here -- leaving the box is a result to be reported.
    rng = np.random.default_rng(seed)
    ngh = rules.n_gh or 5
    S = np.asarray(S0, dtype=float).copy()
    out, off = [], []
    for t in range(T + burn):
        Sm = np.atleast_2d(S)
        x = np.array([float(rules.eval(k, 0, Sm)[0]) for k in SOLVE7])
        try:
            _, o = point_residuals(S, 0, x, rules, cal, ss, sproc,
                                   n_gh=ngh, no_default=no_default)
        except (ValueError, RuntimeError, FloatingPointError):
            break
        if t >= burn:
            out.append(S.copy())
            off.append(rules.grid.outside(S)[0])
        Sn = S.copy()
        Sn[IK_D], Sn[IK_F] = x[2], x[3]                 # K' = Kp
        Sn[IP_D], Sn[IP_F] = o["Pp_D"], o["Pp_F"]
        Sn[IBDD], Sn[IBDF] = o["b_D_D_new"], o["b_D_F_new"]
        Sn[IBFD] = o["b_F_D_new"]
        Sn[IV] = o["Vp_dep"]
        Sn[IS] = ((1.0 - sproc["rho_s"]) * sproc["s_star"] + sproc["rho_s"] * S[IS]
                  + sproc["sigma_s"] * rng.standard_normal())
        Sn[IZ] = (1.0 - sproc["rho_z"]) * sproc["z_star"] + sproc["rho_z"] * S[IZ]
        # the SOLVER only knows the box; a path that leaves it is evaluated on the
        # clipped state, exactly as the continuation would be
        S = rules.grid.clip(Sn)[0]
    return np.array(out), np.array(off)


def euler_errors(rules, cal, ss, sproc, states, no_default=False):
    # THE RESIDUALS AT THE RULE-IMPLIED POLICY, AT EVERY VISITED STATE.
    ngh = rules.n_gh or 5
    err = np.full((states.shape[0], len(_EQN)), np.nan)
    nres = len(_EQN) - len(DERIVED)
    for i, S in enumerate(states):
        Sm = np.atleast_2d(S)
        x = np.array([float(rules.eval(k, 0, Sm)[0]) for k in SOLVE7])
        try:
            res, out = point_residuals(S, 0, x, rules, cal, ss, sproc,
                                       n_gh=ngh, no_default=no_default)
        except (ValueError, RuntimeError, ArithmeticError):
            continue
        err[i, :nres] = np.abs(res)
        for q, k in enumerate(DERIVED):
            g = float(rules.eval(k, 0, Sm)[0])
            err[i, nres + q] = abs(to_fit(k, g) - to_fit(k, out[k]))
    return err


def accuracy_report(rules, cal, ss, sproc, S0, T=1500, burn=200, seed=0,
                    no_default=False, label=""):
    # SIMULATE, MEASURE, PRINT. Returns the per-equation log10 mean/max for reuse.
    states, off = simulate(rules, cal, ss, sproc, S0, T=T, burn=burn, seed=seed,
                           no_default=no_default)
    if states.size == 0:
        print("  accuracy: simulation failed at the first state")
        return None
    err = euler_errors(rules, cal, ss, sproc, states, no_default=no_default)
    ok = np.isfinite(err).all(axis=1)
    if not ok.any():
        print(f"  accuracy: no simulated state evaluated ({states.shape[0]} visited)")
        return None
    e = np.maximum(err[ok], 1e-16)
    l10 = np.log10(e)
    head = f"  EULER-EQUATION ERRORS on the ergodic set{(' — ' + label) if label else ''}"
    print(f"\n{head}")
    print(f"   {states.shape[0]} simulated states ({ok.sum()} evaluable), "
          f"n_gh={rules.n_gh or 5}")
    print("   equation     mean log10|R|   max log10|R|   90th pct")
    for j, nm in enumerate(_EQN):
        print(f"   {nm:10s} {l10[:, j].mean():13.2f} {l10[:, j].max():14.2f}"
              f" {np.percentile(l10[:, j], 90):10.2f}")
    print(f"   {'ALL':10s} {l10.mean():13.2f} {l10.max():14.2f}"
          f" {np.percentile(l10, 90):10.2f}")
    names = STATE_NAMES      # not a local copy: this one still said W_D after the
                             # slot changed to carry V_dep
    frac = (off > 0).mean(axis=0)
    print("   box escapes (share of simulated periods outside, per state):")
    print("     " + "  ".join(f"{n}={100*f:.1f}%" for n, f in zip(names, frac)))
    print(f"     worst overshoot = {100*off.max():.2f}% of box width"
          f"   any-dimension escape = {100*(off > 0).any(axis=1).mean():.1f}% of periods")
    pd = default_prob(states[:, IS])
    print(f"   simulated p^d: mean {100*pd.mean():.3f}%  "
          f"p5 {100*np.percentile(pd, 5):.4f}%  p95 {100*np.percentile(pd, 95):.3f}%")
    return dict(l10_mean=l10.mean(axis=0), l10_max=l10.max(axis=0),
                off_frac=frac, states=states)
