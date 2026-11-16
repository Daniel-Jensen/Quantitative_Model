# HETEROGENEOUS-HOUSEHOLD DEPOSIT CLEARING ALONG AN IMPULSE PATH OF THE GLOBAL MODEL.
#   WT=dff4 python3 ha_path.py <run tag> [T] [workers]
#
# The global rules (banks, firms, government, Eurosystem, expectations) are kept. D's representative
# household is replaced by the steady state's heterogeneous households; F (8 times D) keeps its
# representative household, whose Euler anchors the union deposit rate. The outer unknowns are D's
# deposit-rate path and D's deposit path. Given both, the period map drops D's deposit Euler and the
# deposit-UIP condition and is cleared every quarter (12 unknowns); the outer residuals are deposit-UIP and
# household clearing (heterogeneous-household deposits = the deposit path) in every quarter. Pegging a
# deposit quantity alone made the map explosive (a 1e-15 error grew ~75-fold a quarter); with the price
# given too it is stable to 1e-15 over 30 quarters. For a guess of those paths the period map is cleared every quarter
# (12 unknowns, 12 equations); its prices and incomes go to the steady state's heterogeneous-household
# block (EGM backward from the SS policy, Young lottery forward from the SS distribution); Newton updates
# the deposit paths until household deposits equal them in every quarter. Each household receives the
# realized real deposit return on its own deposits, labour income in proportion to productivity, and an
# equal share of all other household resources in the period map, so the households' aggregate budget is
# the map's budget and the goods market keeps clearing.
import json
import os
import pickle
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / os.environ.get("WT", "dff4")))
sys.path.insert(0, str(HERE))

import numpy as np

_CTX = {}  # per-process model context (the pool's workers rebuild it)
FD_H = 1e-8  # finite-difference step on rates (deposits get 1e3 times it)


def context(tag):
    # RULES, CALIBRATION, SS AND HOUSEHOLD OBJECTS FOR ONE SAVED RUN.
    if _CTX.get("tag") == tag:
        return _CTX
    import argparse
    import tpi_national as tn
    meta = json.loads((HERE / "tpi_runs" / f"{tag}.json").read_text())["args"]
    a = argparse.Namespace(**{k: meta.get(k) for k in ("market", "n_jobs", "cap", "lamZ", "phi", "keyD",
                                                       "theta", "warranted", "keyLoss")})
    cal, ss, sproc = tn.setup(a)
    rr = pickle.load(open(HERE / "tpi_runs" / f"{tag}_rules.pkl", "rb"))
    cal["tpi_cap_bp"] = float(getattr(rr["on"], "tpi_cap_solved", None) or cal["tpi_cap_bp"])
    _CTX.clear()
    _CTX.update(tag=tag, cal=cal, ss=ss, sproc=sproc, rules=rr)
    return _CTX


def _idx():
    # THE GIVEN UNKNOWNS (D'S DEPOSIT RATE AND DEPOSITS) AND THE DROPPED RESIDUALS (D'S EULER, UIP).
    from global_projection.solver_recursive.decision_rules import SOLVE
    return (SOLVE.index("rdep_D"), SOLVE.index("A_D")), (4, 5)


def clear(rules, cal, ss, sproc, S, rD, aD, x0):
    # THE PERIOD MAP AT S WITH D'S DEPOSIT RATE AND HOUSEHOLD DEPOSITS GIVEN.
    from scipy.optimize import root
    from global_projection.solver_recursive.point_map import point_residuals
    (ir, ia), drop = _idx()
    keep_x = [i for i in range(len(x0)) if i not in (ir, ia)]
    keep_r = [i for i in range(len(x0)) if i not in drop]
    ngh = rules.n_gh or 5

    def full(z):
        x = np.array(x0, dtype=float)
        x[keep_x] = z
        x[ir], x[ia] = rD, aD
        return x

    def f(z):
        try:
            return point_residuals(S, 0, full(z), rules, cal, ss, sproc, n_gh=ngh)[0][keep_r]
        except (ValueError, RuntimeError, ArithmeticError):
            return np.full(len(keep_r), 10.0)

    sol = root(f, np.asarray(x0, dtype=float)[keep_x], method="hybr", tol=1e-12)
    x = full(sol.x)
    res, o = point_residuals(S, 0, x, rules, cal, ss, sproc, n_gh=ngh)
    o["_x"] = x
    o["_uip"] = float(res[drop[1]])
    o["_resid"] = float(np.max(np.abs(sol.fun)))
    return o


def simulate(econ, r_D, A_D, s_shock, start=0, base=None):
    # THE MAP PATH FOR GIVEN DEPOSIT PATHS; QUARTERS BEFORE start ARE COPIED FROM base.
    from global_projection.solver_recursive.recursive_experiment import advance
    from global_projection.solver_recursive.state_grid import IS
    from global_projection.solver_recursive.decision_rules import SOLVE
    c = _CTX
    rules, cal = c["rules"][econ], dict(c["cal"], tpi_on=(econ == "on"))
    T = len(A_D)
    if base is None or start == 0:
        S = c["S0"][econ].copy()
        if s_shock is not None:
            S[IS] = s_shock
        outs, start = [], 0
    else:
        outs = list(base["outs"][:start])
        S = base["S"][start].copy()
    Ss = list(base["S"][:start]) if (base is not None and start > 0) else []
    for t in range(start, T):
        Ss.append(S.copy())
        x0 = (base["outs"][t]["_x"] if base is not None else
              np.array([float(rules.eval(k, 0, np.atleast_2d(S))[0]) for k in SOLVE]))
        o = clear(rules, cal, c["ss"], c["sproc"], S, r_D[t], A_D[t], x0)
        outs.append(o)
        S = advance(o, S, c["sproc"], rules.grid)
    Ss.append(S.copy())
    return dict(outs=outs, S=Ss)


def households(path, country, A_guess):
    # HETEROGENEOUS-HOUSEHOLD DEPOSITS AND CONSUMPTION FOR THE PRICES AND INCOMES OF ONE MAP PATH.
    from global_projection.blocks.household import solve_backward_transition
    from global_projection.blocks.distribution import forward_paths, aggregate_assets
    c = _CTX
    cal, ss = c["cal"], c["ss"]
    outs = path["outs"]
    T = len(outs)
    P = np.array([o[f"P_CES_{country}"] for o in outs])
    R = 1.0 + np.array([o[f"rdep_{country}"] for o in outs])
    N = np.array([o["_x"][0 if country == "D" else 1] for o in outs])
    w = np.array([o[f"w_{country}"] for o in outs])
    W = np.array([o[f"W_{country}"] for o in outs])
    inc = np.array([o[f"inc_{country}"] for o in outs])
    rest = c["rest"][country]
    # the realized basket return on deposits carried into t
    rc = np.empty(T)
    rc[0] = rest["R"] * rest["P"] / P[0] - 1.0
    rc[1:] = R[:-1] * P[:-1] / P[1:] - 1.0
    # both boundary conditions at the rest point's own prices (the SS's would promise a price jump)
    D0, c_term = rest["D0"], rest["c_term"]
    a_grid = ss[f"a_grid_{country}"]
    A_lag = np.concatenate([[aggregate_assets(D0, a_grid)], A_guess[:-1]])
    labour = (w / P) * N
    # everything else in the map's household resources, shared equally
    tau = W / P + inc - labour - (1.0 + rc) * A_lag
    e = ss[f"e_{country}"]
    y = labour[:, None] * e[None, :] + tau[:, None]
    fr, chi = cal[f"frisch_{country}"], cal[f"chi_{country}"]
    vN = chi * N ** (1 + 1 / fr) / (1 + 1 / fr)
    # the deposit contracted at T-1 pays its own rate at T, where the basket price is back at rest
    r_path = np.concatenate([rc, [R[-1] * P[-1] / rest["P"] - 1.0]])
    c_path, a_pol = solve_backward_transition(a_grid, ss[f"Pi_{country}"], r_path, y,
                                              c_term, rest["beta"],
                                              cal[f"sigma_{country}"], cal[f"a_min_{country}"],
                                              vN_path=vN)
    A_path, C_path, _ = forward_paths(D0, a_pol, c_path, a_grid, ss[f"Pi_{country}"])
    return A_path, C_path


def residual(econ, u, s_shock, base=None, start=0):
    # DEPOSIT-UIP (x100) AND HOUSEHOLD CLEARING (RELATIVE) IN EVERY QUARTER, AND THE MAP PATH.
    T = len(u) // 2
    r_D, A_D = u[:T], u[T:]
    path = simulate(econ, r_D, A_D, s_shock, start, base)
    hD, _ = households(path, "D", A_D)
    uip = np.array([o["_uip"] for o in path["outs"]])
    return np.concatenate([100.0 * uip, (hD - A_D) / _CTX["ss"]["A_D_ss"]]), path


def _column(args):
    # ONE FINITE-DIFFERENCE COLUMN, RE-SIMULATING THE MAP ONLY FROM THE PERTURBED QUARTER ON.
    tag, econ, u, j, h, s_shock, base, f0 = args
    context(tag)
    _setup_rest(econ)
    T = len(u) // 2
    hj = h if j < T else 1e3 * h  # rates are ~1e-3, deposits ~8
    up = u.copy()
    up[j] += hj
    fj, _ = residual(econ, up, s_shock, base, start=j % T)
    return j, (fj - f0) / hj


def _setup_rest(econ):
    # THE REST POINT OF ONE ECONOMY AND ITS PRICES (THE PATH'S QUARTER -1).
    from global_projection.solver_recursive.recursive_experiment import (
        stochastic_rest_point, read_exact)
    c = _CTX
    c.setdefault("S0", {})
    c.setdefault("rest_all", {})
    if econ in c["S0"]:
        c["rest"] = c["rest_all"][econ]
        return
    cal = dict(c["cal"], tpi_on=(econ == "on"))
    S0 = stochastic_rest_point(c["rules"][econ], cal, c["ss"], c["sproc"], verbose=False)
    o = read_exact(c["rules"][econ], cal, c["ss"], c["sproc"], S0)
    c["S0"][econ] = S0
    c["rest_all"][econ] = {k: dict(R=1.0 + o[f"rdep_{k}"], P=o[f"P_CES_{k}"], A=o[f"A_{k}"])
                           for k in ("D", "F")}
    c["rest_all"][econ]["D"].update(_rest_households(o, "D"))
    c["rest"] = c["rest_all"][econ]


def _rest_households(o, country):
    # THE HOUSEHOLDS' STATIONARY POLICY AND DISTRIBUTION AT THE REST POINT'S PRICES AND INCOMES.
    from global_projection.blocks.household import solve_steady_state_household
    from global_projection.blocks.distribution import stationary_distribution, aggregate_assets
    c = _CTX
    cal, ss = c["cal"], c["ss"]
    P, N = o[f"P_CES_{country}"], o["_x"][0 if country == "D" else 1]
    r = o[f"rdep_{country}"]  # prices constant at rest, so the realized basket return is the rate
    labour = o[f"w_{country}"] / P * N
    tau = o[f"W_{country}"] / P + o[f"inc_{country}"] - labour - (1.0 + r) * o[f"A_{country}"]
    e = ss[f"e_{country}"]
    fr, chi = cal[f"frisch_{country}"], cal[f"chi_{country}"]
    vN = chi * N ** (1 + 1 / fr) / (1 + 1 / fr)
    a_grid, Pi = ss[f"a_grid_{country}"], ss[f"Pi_{country}"]
    pi_e = np.full(len(e), 1.0 / len(e))
    y = labour * e + tau

    def stationary(beta):
        c_pol, a_pol = solve_steady_state_household(a_grid, Pi, r, y, beta, cal[f"sigma_{country}"],
                                                    cal[f"a_min_{country}"], 1e-12, vN_ss=vN)
        D = stationary_distribution(a_pol, a_grid, Pi, pi_e, 1e-13)
        return c_pol, D, aggregate_assets(D, a_grid)

    # the SS's calibration rule applied at the rest point: beta clears household deposits there
    from scipy.optimize import brentq
    A_target = o[f"A_{country}"]
    b0 = ss[f"beta_{country}_ss"]
    beta = brentq(lambda b: stationary(b)[2] - A_target, b0 - 2e-3, min(b0 + 2e-3, 1 / (1 + r) - 1e-6),
                  xtol=1e-12)
    c_pol, D, A = stationary(beta)
    return dict(r=r, c_term=c_pol, D0=D, A_HA=A, beta=beta, beta_ss=b0,
                A_HA_at_ss_beta=stationary(b0)[2])


def solve_path(tag, econ, s_shock, T, workers, tol=1e-6, maxit=60, verbose=True, J0=None, u0=None):
    # NEWTON ON THE 2T DEPOSIT PATHS: ONE FINITE-DIFFERENCE JACOBIAN, THEN BROYDEN UPDATES.
    from concurrent.futures import ProcessPoolExecutor
    context(tag)
    _setup_rest(econ)
    c = _CTX
    # start from the representative-household path of the same shock (exact reads)
    ra = ra_path(tag, econ, s_shock, T)["outs"]
    u = np.array([float(o["rdep_D"]) for o in ra] + [float(o["A_D"]) for o in ra]) if u0 is None else u0.copy()
    f, path = residual(econ, u, s_shock)
    t0 = time.perf_counter()
    if verbose:
        print(f"    [{econ}, shock={s_shock is not None}] start: max|A_HA - A| = {np.max(np.abs(f)):.2e}",
              flush=True)
    def jacobian(u, path, f):
        # THE FINITE-DIFFERENCE JACOBIAN, COLUMNS IN PARALLEL.
        jobs = [(tag, econ, u, j, FD_H, s_shock, path, f) for j in range(2 * T)]
        J = np.empty((2 * T, 2 * T))
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for j, col in ex.map(_column, jobs, chunksize=max(1, (2 * T) // (4 * workers))):
                J[:, j] = col
        if verbose:
            print(f"    Jacobian {2*T}x{2*T} at {time.perf_counter() - t0:.0f}s", flush=True)
        return J

    # a Jacobian handed in (the reference path's) is used as the first quasi-Newton matrix, as SSJ does
    J = jacobian(u, path, f) if J0 is None else J0.copy()
    refreshed = J0 is None
    last_fresh = 0
    # Levenberg-Marquardt on the FD Jacobian (Broyden-updated on accepted steps): the HA deposit supply is
    # very rate-elastic, so J is ill-conditioned and the plain Newton step overshoots into the nonlinear region
    lam = 1e-3
    for it in range(1, maxit + 1):
        JTJ, JTf = J.T @ J, J.T @ f
        d = np.diag(JTJ).copy()
        accepted = False
        for _ in range(10):
            du = -np.linalg.solve(JTJ + lam * np.diag(d), JTf)
            fn, pn = residual(econ, u + du, s_shock, path)
            # accept on the sum of squares the LM step minimises (the max-norm stalls on spread-out residuals)
            if float(fn @ fn) < float(f @ f):
                accepted = True
                break
            lam *= 4.0
        if not accepted:
            if refreshed:
                if verbose:
                    print(f"    it {it}: no descent at lambda {lam:.1e} after a refresh; stopping", flush=True)
                break
            J, refreshed, lam, last_fresh = jacobian(u, path, f), True, 1e-3, it
            continue
        refreshed = False
        # nearly undamped steps that still cut |F| by less than 20%: Broyden has drifted
        slow = float(fn @ fn) > 0.64 * float(f @ f) and lam < 1e-6 and it - last_fresh >= 5
        J += np.outer((fn - f) - J @ du, du) / float(du @ du)
        u, f, path = u + du, fn, pn
        if slow and np.max(np.abs(f)) >= tol:
            J, last_fresh = jacobian(u, path, f), it
        lam = max(lam / 3.0, 1e-9)
        if verbose:
            print(f"    it {it}: max|F| = {np.max(np.abs(f)):.2e}  (lambda {lam:.1e}, "
                  f"{time.perf_counter() - t0:.0f}s)", flush=True)
        if np.max(np.abs(f)) < tol:
            break
    solve_path.last_J = J
    return u, path, float(np.max(np.abs(f)))


def ra_path(tag, econ, s_shock, T):
    # THE REPRESENTATIVE-HOUSEHOLD PATH (ALL 14 EQUATIONS, EXACT READS) FOR COMPARISON.
    from global_projection.solver_recursive.recursive_experiment import read_exact, advance
    from global_projection.solver_recursive.state_grid import IS
    c = _CTX
    rules, cal = c["rules"][econ], dict(c["cal"], tpi_on=(econ == "on"))
    S = c["S0"][econ].copy()
    if s_shock is not None:
        S[IS] = s_shock
    outs = []
    for _ in range(T):
        o = read_exact(rules, cal, c["ss"], c["sproc"], S)
        outs.append(o)
        S = advance(o, S, c["sproc"], rules.grid)
    return dict(outs=outs)


KEYS = ("Y_D", "C_D", "I_D", "n_D", "Q_bD", "A_D", "rdep_D", "r_wc_D", "mu_D", "m_cb", "M_cb_new",
        "Pi_cb", "C_F", "Y_F", "A_F", "rdep_F", "p")


def table(path):
    # THE RECORDED OBJECTS OF ONE PATH, AS ARRAYS.
    out = {k: np.array([float(o[k]) for o in path["outs"]]) for k in KEYS}
    db = _CTX["cal"]["delta_b_D"]
    qD, qF = out["Q_bD"], np.array([float(o["Q_bF"]) for o in path["outs"]])
    out["sov_bp"] = 4e4 * (db * (1 - qD) / qD - _CTX["cal"]["delta_b_F"] * (1 - qF) / qF)
    out["spread_bp"] = 4e4 * (out["r_wc_D"] - out["rdep_D"])
    return out


def irf(shock, ref):
    # DIFFERENCES OF A SHOCKED PATH AGAINST ITS REFERENCE PATH.
    pct = ("Y_D", "C_D", "I_D", "n_D", "Q_bD", "A_D", "C_F", "Y_F", "A_F")
    d = {k: 100 * (shock[k] / ref[k] - 1) for k in pct}
    for k in ("rdep_D", "rdep_F", "r_wc_D"):
        d[k] = 4e4 * (shock[k] - ref[k])
    for k in ("sov_bp", "spread_bp"):
        d[k] = shock[k] - ref[k]
        d[k + "_level"] = shock[k]
    d["m_cb"] = 100 * shock["m_cb"] / _CTX["cal"]["B_gov_D_ss"]
    d["M_cb"] = 100 * shock["M_cb_new"] / _CTX["cal"]["B_gov_D_ss"]
    d["Pi_cb"] = 100 * shock["Pi_cb"] / _CTX["ss"]["ss_firm_D"]["Y_ss"]
    return d


def main():
    # BOTH ECONOMIES OF ONE RUN: HA AND RA IRFS OF THE HEADLINE SHOCK.
    tag = sys.argv[1]
    T = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    econs = sys.argv[4].split(",") if len(sys.argv) > 4 else ["off", "on"]
    from global_projection.solver_recursive.recursive_experiment import s_from_pd
    context(tag)
    s1 = s_from_pd(0.02)
    res = {}
    for econ in econs:
        _setup_rest(econ)
        print(f"\n[{tag}] economy TPI {econ}: heterogeneous-household paths, T = {T}", flush=True)
        u_ref, ref, e0 = solve_path(tag, econ, None, T, workers)
        # start the shocked solve from the RA shocked path plus the HA correction of the reference path
        ra0, ra1 = (ra_path(tag, econ, s, T)["outs"] for s in (None, s1))
        vec = lambda outs: np.array([float(o["rdep_D"]) for o in outs] + [float(o["A_D"]) for o in outs])
        _, sho, e1 = solve_path(tag, econ, s1, T, workers, u0=vec(ra1) + (u_ref - vec(ra0)))
        ha = irf(table(sho), table(ref))
        ra = irf(table(ra_path(tag, econ, s1, T)), table(ra_path(tag, econ, None, T)))
        res[econ] = dict(ha=ha, ra=ra, err=(e0, e1),
                         map_resid=max(max(o["_resid"] for o in p["outs"]) for p in (ref, sho)))
        for lab, d in (("HA", ha), ("RA", ra)):
            print(f"   {lab}: Y " + " ".join(f"{x:+.3f}" for x in d["Y_D"][:8])
                  + f" | cumY25 {d['Y_D'][:25].sum():+.2f}  C0 {d['C_D'][0]:+.3f}  I0 {d['I_D'][0]:+.3f}"
                  f"  n0 {d['n_D'][0]:+.2f}  rdep0 {d['rdep_D'][0]:+.1f}bp  sov0 {d['sov_bp_level'][0]:.1f}bp"
                  f"  m0 {d['m_cb'][0]:.1f}%", flush=True)
    clean = {e: {"err": r["err"], "map_resid": r["map_resid"],
                 "ha": {k: v.tolist() for k, v in r["ha"].items()},
                 "ra": {k: v.tolist() for k, v in r["ra"].items()}} for e, r in res.items()}
    (HERE / "tpi_runs" / f"{tag}_ha_T{T}.json").write_text(json.dumps(clean, indent=1))


if __name__ == "__main__":
    main()
