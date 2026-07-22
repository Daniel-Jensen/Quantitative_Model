# STACKED-SYSTEM SOLVER UTILITIES: PARALLEL FD JACOBIAN + DAMPED NEWTON WITH REUSE.
# The point is Jacobian REUSE: jac_cache is a caller-owned dict, so the expensive
# FD Jacobian (7T+1 evaluations) is built once per system kind and shared across
# warm re-solves of the risk fixed point (scipy hybr rebuilds it every call).
# Workers are spawned (macOS), so every entry point needs its __main__ guard.
import os
import numpy as np
from scipy.linalg import lu_factor, lu_solve

_SQRT_EPS = np.sqrt(np.finfo(float).eps)   # ~1.49e-8, MINPACK's default step
_WALL = 5.0        # |dF| this large means the perturbation hit a penalty wall
_WORKER = {}       # per-worker residual context (set by _init_worker)


def _fd_steps(y0):
    # PER-COMPONENT FORWARD-DIFFERENCE STEPS: RELATIVE WITH AN O(1) FLOOR.
    # The floor keeps the step sane for deposit rates near 0 (else eps*|x| ~ FD noise).
    return _SQRT_EPS * np.maximum(np.abs(y0), 1.0)


def _fd_columns(residual, y0, F0, cols, steps):
    # FD COLUMNS FOR `cols`, BACKING OFF THE STEP IF A PERTURBATION HITS THE PENALTY WALL.
    out = np.empty((F0.size, len(cols)))
    for k, j in enumerate(cols):
        h = steps[j]
        y = y0.copy()
        y[j] += h
        Fj = residual(y)
        if np.max(np.abs(Fj - F0)) >= _WALL:          # penalty wall - back off
            y = y0.copy()
            y[j] += h / 8.0
            Fj_small = residual(y)
            if np.max(np.abs(Fj_small - F0)) < _WALL:
                Fj, h = Fj_small, h / 8.0
        out[:, k] = (Fj - F0) / h
    return out


def _init_worker(spec):
    # MULTIPROCESSING WORKER INIT: REBUILD THE RESIDUAL FROM THE PICKLABLE SPEC.
    import transition
    _WORKER["residual"] = transition.make_residual(spec)


def _worker_columns(args):
    # WORKER ENTRY: COMPUTE ITS ASSIGNED FD COLUMNS.
    y0, F0, cols, steps = args
    return cols, _fd_columns(_WORKER["residual"], y0, F0, cols, steps)


def fd_jacobian(residual, y0, F0, spec=None, n_jobs=0, verbose=False):
    # FORWARD-DIFFERENCE JACOBIAN dF/dy AT y0 (PARALLEL OVER COLUMNS, SERIAL FALLBACK).
    n = y0.size
    steps = _fd_steps(y0)
    if n_jobs == 0:
        n_jobs = os.cpu_count() or 1
    n_jobs = min(n_jobs, n)

    if spec is not None and n_jobs > 1:   # spec != None enables the parallel path
        try:
            import multiprocessing as mp
            ctx = mp.get_context("spawn")
            col_chunks = np.array_split(np.arange(n), n_jobs)
            J = np.empty((F0.size, n))
            with ctx.Pool(n_jobs, initializer=_init_worker,
                          initargs=(spec,)) as pool:
                for cols, block in pool.imap_unordered(
                        _worker_columns,
                        [(y0, F0, chunk, steps) for chunk in col_chunks]):
                    J[:, cols] = block
            return J
        except Exception as e:                        # noqa: BLE001
            print(f"  [fd_jacobian] parallel build failed ({e}); "
                  "falling back to serial.")

    J = np.empty((F0.size, n))
    J[:] = _fd_columns(residual, y0, F0, np.arange(n), steps)
    return J


def newton_solve(residual, y0, F0=None, jac_cache=None, accept_tol=1e-10,
                 build_jac=None, max_iter=80, max_rebuilds=8, verbose=False):
    # DAMPED NEWTON ON max|F| WITH A REUSABLE LU-FACTORIZED JACOBIAN (NEVER RAISES).
    # jac_cache (caller-owned) may carry a Jacobian from a nearby system; Broyden
    # updates + stall-triggered rebuilds absorb the mismatch. Failure is declared
    # only when the line search fails on a FRESHLY built Jacobian.
    if jac_cache is None:
        jac_cache = {}
    y = np.array(y0, dtype=float)
    F = residual(y) if F0 is None else np.asarray(F0)
    fnorm = np.max(np.abs(F))
    if fnorm <= accept_tol:
        return y, F, True
    if fnorm >= _WALL:
        # starting point on the penalty wall -> FD would difference wall vs wall
        # (zero Jacobian); let the caller's trust-region fallback take over
        return y, F, False

    def _factor(J):
        # STORE AND LU-FACTORIZE A JACOBIAN IN THE CACHE.
        if not np.all(np.isfinite(J)):
            raise np.linalg.LinAlgError("non-finite Jacobian")
        jac_cache["J"] = J
        jac_cache["lu"] = lu_factor(J)

    try:
        if jac_cache.get("J") is None or jac_cache["J"].shape[0] != y.size:
            _factor(build_jac(y, F))
        elif "lu" not in jac_cache:
            _factor(jac_cache["J"])
    except np.linalg.LinAlgError:
        jac_cache.pop("J", None)
        jac_cache.pop("lu", None)
        return y, F, False

    rebuilds = 0
    just_rebuilt = False
    for it in range(max_iter):
        dy = -lu_solve(jac_cache["lu"], F)
        accepted = False
        for lam in (1.0, 0.5, 0.25, 0.1):   # backtracking line search
            y_new = y + lam * dy
            F_new = residual(y_new)
            fn_new = np.max(np.abs(F_new))
            if np.isfinite(fn_new) and fn_new < fnorm:
                accepted = True
                break
        if accepted:
            just_rebuilt = False
            s = y_new - y
            dF = F_new - F
            J = jac_cache["J"]
            J += np.outer((dF - J @ s) / (s @ s), s)   # Broyden "good" update
            try:
                _factor(J)
            except np.linalg.LinAlgError:
                _factor(build_jac(y_new, F_new))
            y, F, fnorm = y_new, F_new, fn_new
            if verbose:
                print(f"    [newton] it {it + 1:2d}: max|F|={fnorm:.3e} (step {lam})")
            if fnorm <= accept_tol:
                return y, F, True
        else:
            if just_rebuilt or rebuilds >= max_rebuilds:
                return y, F, False   # line search failed on a fresh Jacobian -> true stall
            rebuilds += 1
            if verbose:
                print(f"    [newton] stall at max|F|={fnorm:.3e}; "
                      f"rebuilding Jacobian ({rebuilds}/{max_rebuilds})")
            try:
                _factor(build_jac(y, F))
            except np.linalg.LinAlgError:
                jac_cache.pop("J", None)
                jac_cache.pop("lu", None)
                return y, F, False
            just_rebuilt = True
    return y, F, fnorm <= accept_tol
