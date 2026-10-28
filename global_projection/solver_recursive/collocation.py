# GLOBAL CHEBYSHEV COLLOCATION: ONE NEWTON ON THE RULES' VALUES AT THE NODES (BOCOLA'S DESIGN).
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import newton_krylov
from scipy.optimize import NoConvergence

from global_projection.solver_recursive.decision_rules import (SOLVE, DERIVED, STORE_RULES,
                                             to_fit, from_fit)
from global_projection.solver_recursive.point_map import point_residuals

# the unknowns are stacked rule by rule, then regime, then grid point
N_RES_POINT = len(SOLVE)
N_RES = N_RES_POINT + len(DERIVED)  # plus an identity residual for each derived rule
RES_NAMES = ("cap_D", "cap_F", "lab_D", "lab_F", "euler_D", "uip", "goods_D",
             "bondD_D", "bondD_F", "euler_F", "dep_clear", "bondF_F", "bondF_D", "tpi"
             ) + tuple(f"id_{k}" for k in DERIVED)
assert len(RES_NAMES) == N_RES, (len(RES_NAMES), N_RES)
_BIG = 1e3  # residual for a point that cannot be evaluated


# the rule sets a solve can take as unknowns: the normal box's, and the post-default box's
SETS = ("normal", "post")


def _sets(rules, sets=SETS):
    # THE RULE SETS WHOSE NODE VALUES ARE UNKNOWNS (THE POST-DEFAULT ONE ONLY IF ATTACHED).
    return ([rules] if "normal" in sets else []) + (
        [rules.post] if "post" in sets and rules.post is not None else [])


def pack(rules, regimes=(0, 1), sets=SETS):
    # RULE VALUES -> THE FLAT UNKNOWN VECTOR (IN THE FITTED TRANSFORM), SET BY SET.
    return np.concatenate([to_fit(k, rs.vals[k][d]) for rs in _sets(rules, sets)
                           for k in STORE_RULES for d in regimes])


def unpack(theta, rules, regimes=(0, 1), sets=SETS):
    # THE FLAT VECTOR -> RULE VALUES, ONE DICT PER SET (THE INVERSE OF pack).
    outs, i = [], 0
    for rs in _sets(rules, sets):
        n, out = rs.grid.n, {}
        for k in STORE_RULES:
            out[k] = {}
            for d in regimes:
                out[k][d] = from_fit(k, theta[i:i + n])
                i += n
        outs.append(out)
    return outs


def write_back(theta, rules, regimes=(0, 1), sets=SETS):
    # WRITE A SOLUTION INTO THE RULE SETS.
    for rs, vals in zip(_sets(rules, sets), unpack(theta, rules, regimes, sets)):
        _install_guess(rs, vals, regimes)
    return rules


def _install_guess(rules, vals, regimes):
    # FIT EVERY RULE EXACTLY TO THE CURRENT GUESS.
    for k in STORE_RULES:
        for d in regimes:
            rules.set_values(k, d, vals[k][d])


def _point_block(i, d, vals, rules, cal, ss, sproc, n_gh, no_default):
    # THE RESIDUALS AT ONE GRID POINT OF ONE RULE SET IN ONE REGIME.
    out_row = np.empty(N_RES)
    x = np.array([vals[k][d][i] for k in SOLVE])
    try:
        r, out = point_residuals(rules.grid.points[i], d, x, rules, cal, ss, sproc,
                                 n_gh=n_gh, no_default=no_default)
    except (ValueError, RuntimeError, ArithmeticError):
        # one bad point should cost the step, not the solve
        out_row[:] = _BIG
        return out_row
    out_row[:N_RES_POINT] = r
    # Bocola's identity residual: log(guess / implied)
    for q, k in enumerate(DERIVED):
        out_row[N_RES_POINT + q] = to_fit(k, vals[k][d][i]) - to_fit(k, out[k])
    return out_row


def make_residual(rules, cal, ss, sproc, regimes=(0, 1), no_default=False, n_gh=5,
                  scale=None, sets=SETS):
    # THE GLOBAL RESIDUAL F(theta); EACH CALL REFITS THE RULES TO theta.
    sw = np.ones(N_RES) if scale is None else np.asarray(scale, dtype=float)

    def F(theta):
        # UNPACK, FIT, RUN THE PERIOD MAP AT EVERY POINT OF EVERY SET, STACK.
        sets_ = _sets(rules, sets)
        allvals = unpack(theta, rules, regimes, sets)
        for rs, vals in zip(sets_, allvals):
            _install_guess(rs, vals, regimes)
        blocks = []
        for rs, vals in zip(sets_, allvals):
            res = np.empty((len(regimes), rs.grid.n, N_RES))
            for jd, d in enumerate(regimes):
                for i in range(rs.grid.n):
                    res[jd, i, :] = _point_block(i, d, vals, rs, cal, ss, sproc,
                                                 n_gh, no_default)
            blocks.append(res.ravel())
        # set by set (normal, then post-default), each (regime, point, equation)
        res = np.concatenate(blocks).reshape(-1, N_RES)
        res = np.where(np.isfinite(res), res, _BIG)
        return (res * sw).ravel()

    return F


def residual_table(theta, F):
    # THE WORST RESIDUAL OF EACH EQUATION.
    r = np.abs(F(theta)).reshape(-1, N_RES)
    return dict(zip(RES_NAMES, r.max(axis=0)))


def _fd_jacobian(F, x, f, eps):
    # FORWARD-DIFFERENCE JACOBIAN, COLUMN BY COLUMN.
    J = np.empty((f.size, x.size))
    for i in range(x.size):
        xp = x.copy()
        xp[i] += eps
        J[:, i] = (F(xp) - f) / eps
    return J


# the Jacobian's columns run on every core; workers re-import the model, so don't edit code mid-run
_WORKER = {}


def _worker_init(args):
    # BUILD THE RESIDUAL ONCE IN EACH WORKER.
    _WORKER["F"] = make_residual(*args)


def _worker_columns(x, f, eps, cols):
    # A BLOCK OF JACOBIAN COLUMNS, COMPUTED AS _fd_jacobian DOES.
    F = _WORKER["F"]
    block = np.empty((f.size, len(cols)))
    for j, i in enumerate(cols):
        xp = x.copy()
        xp[i] += eps
        block[:, j] = (F(xp) - f) / eps
    return cols, block


class _PoolJacobian:
    # THE FD JACOBIAN ON A POOL OF SPAWNED PROCESSES.

    def __init__(self, args, n):
        # START n WORKERS.
        self.n = n
        self.pool = ProcessPoolExecutor(max_workers=n, mp_context=mp.get_context("spawn"),
                                        initializer=_worker_init, initargs=(args,))

    def __call__(self, F, x, f, eps):
        # COLUMNS IN 4n CHUNKS, ASSEMBLED IN SERIAL ORDER.
        chunks = np.array_split(np.arange(x.size), 4 * self.n)
        J = np.empty((f.size, x.size))
        jobs = [self.pool.submit(_worker_columns, x, f, eps, c) for c in chunks if c.size]
        for job in jobs:
            cols, block = job.result()
            J[:, cols] = block
        return J

    def close(self):
        # STOP THE WORKERS.
        self.pool.shutdown()


def _n_workers(cal):
    # THE WORKER COUNT: 0 = EVERY CORE, 1 = SERIAL.
    n = int(cal.get("n_jobs", 0))
    return (os.cpu_count() or 1) if n <= 0 else n


def _factor(J):
    # LU FACTORS, OR None IF THE JACOBIAN IS SINGULAR.
    try:
        return lu_factor(J)
    except (np.linalg.LinAlgError, ValueError):
        return None


def _newton_direction(lu, J, f):
    # THE NEWTON STEP, FALLING BACK TO LEAST SQUARES.
    try:
        dx = lu_solve(lu, f) if lu is not None else np.linalg.lstsq(J, f, rcond=None)[0]
    except (np.linalg.LinAlgError, ValueError):
        return np.linalg.lstsq(J, f, rcond=None)[0]
    # a near-singular factor returns inf/nan instead of raising
    return dx if np.all(np.isfinite(dx)) else np.linalg.lstsq(J, f, rcond=None)[0]


def _line_search(F, x, dx, gap, cc, backtrack):
    # DAMPED STEP, HALVING THE DAMPING UNTIL THE RESIDUAL FALLS.
    step = cc
    while True:
        xn = x - step * dx
        fn = F(xn)
        gn = float(fn @ fn)
        if (not backtrack) or gn < gap or step < 1e-5:
            return xn, fn, gn, step
        step *= 0.5


# the Levenberg-Marquardt fallback tries this many regularisations, from (1e-8 s0)^2 up by 10x
LM_TRIES = 30


def _lm_step(F, x, f, gap, J):
    # REGULARISED LEAST-SQUARES STEPS ON THE COLUMN-SCALED JACOBIAN, LAMBDA RAISED UNTIL ONE DESCENDS.
    s = 1.0 / np.maximum(np.linalg.norm(J, axis=0), 1e-300)
    U, S, Vt = np.linalg.svd(J * s, full_matrices=False)
    uf = U.T @ f
    lam = (1e-8 * S[0]) ** 2
    for _ in range(LM_TRIES):
        xt = x - s * (Vt.T @ (S / (S ** 2 + lam) * uf))
        ft = F(xt)
        gt = float(ft @ ft)
        if np.isfinite(gt) and gt < gap:
            return xt, ft, gt, lam / S[0] ** 2
        lam *= 10.0
    return x, f, gap, np.inf


def parsolve(F, x0, cc=1.0, tol=1e-20, maxcount=60, eps=1e-6, verbose=True,
             blowup=100.0, backtrack=True, label="", stall_step=1e-3, jac_every=1,
             floor_tol=1e-8, jacobian=None):
    # DAMPED NEWTON WITH A FINITE-DIFFERENCE JACOBIAN (A PORT OF BOCOLA'S parsolve.m) AND A
    # LEVENBERG-MARQUARDT FALLBACK.
    x = np.asarray(x0, dtype=float).copy()
    f = F(x)
    gap = float(f @ f)
    lu = None
    for count in range(1, maxcount + 1):
        if gap <= tol:
            break
        if lu is None or (count - 1) % jac_every == 0:
            J = (_fd_jacobian if jacobian is None else jacobian)(F, x, f, eps)
            lu = _factor(J)
        dx = _newton_direction(lu, J, f)
        xn, fn, gn, step = _line_search(F, x, dx, gap, cc, backtrack)
        how, lm = f"step = {step:.3g}", False
        # the fallback (2026-10-21): the cold start converged or failed on round-off, and LM
        # converged it in every case tested; where Newton succeeds it never runs
        if not gn < gap or step < stall_step:
            xl, fl, gl, lam = _lm_step(F, x, f, gap, J)
            if gl < min(gn, gap):
                xn, fn, gn, how, lm = xl, fl, gl, f"LM lambda = {lam:.0e} s0^2", True
        improved = gn < gap
        if improved:
            x, f, gap = xn, fn, gn
        if verbose:
            print(f"    [parsolve{label} {count:2d}] sum|F|^2 = {gap:.3e}   "
                  f"max|F| = {np.max(np.abs(f)):.3e}   {how}"
                  f"{'' if improved else '   (no improvement)'}")
        if not np.isfinite(gap) or gap > blowup:
            return x, False, count, np.max(np.abs(f))
        if not improved or (step < stall_step and not lm):
            # no step helps: that is convergence only if the residual is at its arithmetic floor
            worst = float(np.max(np.abs(f)))
            return x, worst <= floor_tol, count, worst
    return x, gap <= tol, count, float(np.max(np.abs(f)))


def krylov_solve(F, x0, f_tol=1e-11, maxiter=60, verbose=True, label=""):
    # JACOBIAN-FREE NEWTON-KRYLOV (NOT USED: IT STALLS ON THIS SYSTEM).
    it = {"n": 0}

    def cb(x, fx):
        it["n"] += 1
        if verbose:
            print(f"    [krylov{label} {it['n']:2d}] max|F| = {np.max(np.abs(fx)):.3e}")
    try:
        x = newton_krylov(F, x0, f_tol=f_tol, maxiter=maxiter, callback=cb,
                          method="lgmres", line_search="armijo", verbose=False)
        f = F(x)
        return x, bool(np.max(np.abs(f)) <= f_tol * 10), it["n"], float(np.max(np.abs(f)))
    except NoConvergence as e:
        x = np.asarray(e.args[0], dtype=float)
        f = F(x)
        return x, False, it["n"], float(np.max(np.abs(f)))


# acceptance: max|F| <= 1e-9, just above the period map's arithmetic floor
TOL_MAXF = 1e-9
# above this many unknowns "auto" would use Krylov; set out of reach on purpose
BACKEND_DENSE_MAX = 10 ** 9


def solve_collocation(rules, cal, ss, sproc, regimes=(0, 1), no_default=False,
                      n_gh=5, backend="auto", tol=TOL_MAXF, maxit=60, cc=1.0,
                      verbose=True, label="", jac_every=1, sets=SETS):
    # SOLVE FOR ALL NODE VALUES AT ONCE AND WRITE THE ANSWER INTO rules.
    F = make_residual(rules, cal, ss, sproc, regimes=regimes,
                      no_default=no_default, n_gh=n_gh, sets=sets)
    x0 = pack(rules, regimes, sets)
    m = x0.size
    if backend == "auto":
        backend = "parsolve" if m <= BACKEND_DENSE_MAX else "krylov"
    if verbose:
        npts = " + ".join(str(rs.grid.n) for rs in _sets(rules, sets))
        print(f"  collocation solve{label}: {m} unknowns "
              f"({len(STORE_RULES)} rules x {len(regimes)} regime(s) x {npts} "
              f"points), backend = {backend}")
    if backend == "parsolve":
        # parsolve stops on the sum of squares: m*tol^2 matches max|F| = tol
        nj = _n_workers(cal)
        jac = (_PoolJacobian((rules, cal, ss, sproc, regimes, no_default, n_gh, None, sets), nj)
               if nj > 1 else None)
        try:
            x, ok, its, worst = parsolve(F, x0, cc=cc, tol=m * tol ** 2, maxcount=maxit,
                                         verbose=verbose, label=label, jac_every=jac_every,
                                         floor_tol=10.0 * tol, jacobian=jac)
        finally:
            if jac is not None:
                jac.close()
    else:
        x, ok, its, worst = krylov_solve(F, x0, f_tol=tol, maxiter=maxit,
                                         verbose=verbose, label=label)
        if not ok:
            # fall back to the dense Newton
            if verbose:
                print(f"    krylov stopped at max|F| = {worst:.2e}; "
                      f"falling back to the dense-Jacobian Newton (jac_every=3)")
            x, ok, its2, worst = parsolve(F, x, cc=cc, tol=m * tol ** 2,
                                          maxcount=max(6, maxit // 4), verbose=verbose,
                                          label=label + " dense", jac_every=3)
            its += its2
    write_back(x, rules, regimes, sets)
    if verbose:
        tab = residual_table(x, F)
        worst_eq = max(tab, key=tab.get)
        print(f"    -> {'converged' if ok else 'STOPPED'} in {its} iterations, "
              f"max|F| = {worst:.2e}  (worst equation: {worst_eq} {tab[worst_eq]:.1e})")
    return ok, its, worst
