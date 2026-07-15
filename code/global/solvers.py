"""Stacked-system solver utilities for the 7T transition path.

Two pieces:
  fd_jacobian  — forward-difference Jacobian of the market-clearing residual,
                 with the independent columns farmed out to a multiprocessing
                 pool (each worker rebuilds the residual from a picklable
                 spec via transition.make_residual).
  newton_solve — damped Newton with an LU-factorized Jacobian, backtracking
                 line search, Broyden (good) rank-one updates, and
                 stall-triggered Jacobian rebuilds.

The point of this module is Jacobian REUSE: `jac_cache` is a plain dict owned
by the caller and updated in place, so the expensive FD Jacobian (7T+1
residual evaluations) is built once per system kind (base path / default
branch) and then shared across all warm re-solves of the risk fixed point.
scipy's hybr cannot do this — it rebuilds its internal FD Jacobian on every
call, which made each warm re-solve cost ~1400 evaluations for ~6 iterations
of actual work.

Multiprocessing note: workers are spawned (macOS default), which re-imports
the caller's __main__ module — every entry-point script must keep its
`if __name__ == "__main__":` guard (they all do).
"""
import os
import numpy as np
from scipy.linalg import lu_factor, lu_solve

_SQRT_EPS = np.sqrt(np.finfo(float).eps)   # ~1.49e-8, MINPACK's default step
_WALL = 5.0        # |ΔF| this large means the perturbation hit a penalty wall
_WORKER = {}       # per-worker residual context (set by _init_worker)


# ─────────────────────────────────────────────────────────────────────────────
# Finite-difference Jacobian (parallel columns)
# ─────────────────────────────────────────────────────────────────────────────

def _fd_steps(y0):
    """Per-component forward-difference steps: relative with an O(1) floor.

    All unknowns are O(1e-3)–O(10) here (N, Kap, rdep, p); the floor keeps the
    step sane for the deposit rates near 0 (MINPACK would use eps·|x| ≈ 1e-11
    there, uncomfortably close to the FD noise floor)."""
    return _SQRT_EPS * np.maximum(np.abs(y0), 1.0)


def _fd_columns(residual, y0, F0, cols, steps):
    """FD columns for index list `cols`; retries a column with a smaller step
    if the perturbation hits the solver's penalty wall (all-10.0 residual)."""
    out = np.empty((F0.size, len(cols)))
    for k, j in enumerate(cols):
        h = steps[j]
        y = y0.copy()
        y[j] += h
        Fj = residual(y)
        if np.max(np.abs(Fj - F0)) >= _WALL:          # penalty wall — back off
            y = y0.copy()
            y[j] += h / 8.0
            Fj_small = residual(y)
            if np.max(np.abs(Fj_small - F0)) < _WALL:
                Fj, h = Fj_small, h / 8.0
        out[:, k] = (Fj - F0) / h
    return out


def _init_worker(spec):
    import transition
    _WORKER["residual"] = transition.make_residual(spec)


def _worker_columns(args):
    y0, F0, cols, steps = args
    return cols, _fd_columns(_WORKER["residual"], y0, F0, cols, steps)


def fd_jacobian(residual, y0, F0, spec=None, n_jobs=0, verbose=False):
    """Forward-difference Jacobian dF/dy at y0 (F0 = residual(y0)).

    spec   : picklable residual spec (transition.make_residual) enabling the
             parallel path; None → serial in-process.
    n_jobs : worker count; 0 → os.cpu_count().  Falls back to serial on any
             multiprocessing failure.
    """
    n = y0.size
    steps = _fd_steps(y0)
    if n_jobs == 0:
        n_jobs = os.cpu_count() or 1
    n_jobs = min(n_jobs, n)

    if spec is not None and n_jobs > 1:
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


# ─────────────────────────────────────────────────────────────────────────────
# Damped Newton with Jacobian reuse
# ─────────────────────────────────────────────────────────────────────────────

def newton_solve(residual, y0, F0=None, jac_cache=None, accept_tol=1e-10,
                 build_jac=None, max_iter=80, max_rebuilds=8, verbose=False):
    """Damped Newton on max|F| with a reusable LU-factorized Jacobian.

    jac_cache : dict persisted BY THE CALLER across solves; keys "J" and "lu"
                are (re)filled here.  An empty dict means "build fresh and
                cache it".  The cached J may come from a slightly different
                system (e.g. risk inputs moved between rounds) — Broyden
                updates plus the stall-triggered rebuild absorb that.
    build_jac : callable (y, F) -> J.

    Stall policy: on a failed line search the Jacobian is rebuilt at the
    current point (strongly nonlinear stretches — e.g. the default-branch
    impact period — go stale within a few steps, and each rebuild is cheap
    and parallel).  Failure is declared only when the line search fails on a
    FRESHLY built Jacobian: at that point even exact-J Newton cannot descend
    and the caller's trust-region fallback takes over.
    Returns (y, F, converged).  Never raises on non-convergence.
    """
    if jac_cache is None:
        jac_cache = {}
    y = np.array(y0, dtype=float)
    F = residual(y) if F0 is None else np.asarray(F0)
    fnorm = np.max(np.abs(F))
    if fnorm <= accept_tol:
        return y, F, True
    if fnorm >= _WALL:
        # Starting point sits on the penalty wall: FD columns would difference
        # wall against wall (zero Jacobian).  Let the caller's trust-region
        # fallback find the feasible region instead.
        return y, F, False

    def _factor(J):
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
        for lam in (1.0, 0.5, 0.25, 0.1):
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
            # Broyden "good" update keeps the cached J tracking the true
            # Jacobian along the path; refactor is ~50 ms at 1400², trivial
            # against a 45 ms residual evaluation.  A wall-poisoned update
            # (F_new on the penalty wall never gets accepted, so dF is always
            # finite here) cannot corrupt the cache.
            J = jac_cache["J"]
            J += np.outer((dF - J @ s) / (s @ s), s)
            try:
                _factor(J)
            except np.linalg.LinAlgError:
                _factor(build_jac(y_new, F_new))
            y, F, fnorm = y_new, F_new, fn_new
            if verbose:
                print(f"    [newton] it {it + 1:2d}: max|F|={fnorm:.3e} "
                      f"(step {lam})")
            if fnorm <= accept_tol:
                return y, F, True
        else:
            if just_rebuilt or rebuilds >= max_rebuilds:
                # line search failed on a fresh Jacobian — true stall
                return y, F, False
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
