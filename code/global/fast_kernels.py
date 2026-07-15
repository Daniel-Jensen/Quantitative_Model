"""JIT-compiled hot kernels for the transition residual (numba, optional).

Two kernels dominate a residual evaluation (~88% of 45.5 ms measured):
  hh_backward  — the household EGM backward induction over the whole path
                 (household.solve_backward_transition's loop);
  dist_forward — the distribution forward simulation with Young-lottery
                 scatter plus the C/A aggregations (distribution.forward_paths'
                 loop).

Both replicate the pure-numpy reference implementations operation-for-
operation (no fastmath, np.interp semantics, side='right' searchsorted,
denom>0 lottery guard), and tests/test_fast_kernels.py asserts agreement to
float precision on random inputs.  When numba is not importable the callers
fall back to the numpy paths automatically; cal["use_numba"]=False forces the
fallback.  cache=True persists compiled artifacts to __pycache__, so spawned
Jacobian workers (solvers.fd_jacobian) load them instead of recompiling.
"""
import numpy as np

try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:                                   # pragma: no cover
    HAVE_NUMBA = False

    def njit(*args, **kwargs):                        # no-op decorator
        def wrap(f):
            return f
        return wrap


@njit(cache=True)
def hh_backward(a_grid, Pi_T, r_path, y_path, c_terminal, beta, sigma, a_min,
                vN_path):
    """EGM backward induction, terminal condition c_terminal (= SS policy).

    Pi_T is Pi.T (C-contiguous).  r_path has length T+1 (today's return plus
    the return entering each period, Fisher-adjusted upstream).  Mirrors
    household.egm_step + solve_backward_transition exactly.
    """
    T, n_e = y_path.shape
    n_a = a_grid.shape[0]
    c_path = np.empty((T, n_a, n_e))
    a_pol_path = np.empty((T, n_a, n_e))

    c_next = c_terminal
    vN_next = vN_path[T - 1]
    for t in range(T - 1, -1, -1):
        r_today = r_path[t]
        r_next = r_path[t + 1]
        vN_today = vN_path[t]

        # GHH composite tomorrow and expected marginal utility
        x_next = np.maximum(c_next - vN_next, 1e-11)
        Eu_next = (x_next ** (-sigma)) @ Pi_T                    # (n_a, n_e)

        # EGM inversion and endogenous asset grid
        x_endo = (beta * (1.0 + r_next) * Eu_next) ** (-1.0 / sigma)
        c_endo = x_endo + vN_today

        for e in range(n_e):
            y_e = y_path[t, e]
            xp = np.empty(n_a)                        # a_endo[:, e]
            fp = np.empty(n_a)                        # c_endo[:, e]
            for i in range(n_a):
                fp[i] = c_endo[i, e]
                xp[i] = (fp[i] + a_grid[i] - y_e) / (1.0 + r_today)
            c_col = np.interp(a_grid, xp, fp)
            a0 = xp[0]
            for i in range(n_a):
                a_i = a_grid[i]
                if a_i < a0:                          # borrowing-constrained
                    c_i = (1.0 + r_today) * a_i + y_e - a_min
                    a_p = a_min
                else:
                    c_i = c_col[i]
                    a_p = (1.0 + r_today) * a_i + y_e - c_i
                    if a_p < a_min:
                        a_p = a_min
                c_path[t, i, e] = c_i
                a_pol_path[t, i, e] = a_p

        c_next = c_path[t]
        vN_next = vN_today

    return c_path, a_pol_path


@njit(cache=True)
def dist_forward(D0, a_pol_path, c_path, a_grid, Pi):
    """Distribution forward simulation with Young (2010) lottery weights.

    Timing convention (matches the transition solver): C_t over the
    start-of-period distribution, A_t over the end-of-period one; D_start[t]
    is the distribution entering period t.
    """
    T = a_pol_path.shape[0]
    n_a, n_e = D0.shape
    A_path = np.empty(T)
    C_path = np.empty(T)
    D_start = np.empty((T + 1, n_a, n_e))

    a_lo_grid = a_grid[0]
    a_hi_grid = a_grid[n_a - 1]

    D = D0.copy()
    for t in range(T):
        c_sum = 0.0
        for i in range(n_a):
            for e in range(n_e):
                D_start[t, i, e] = D[i, e]
                c_sum += D[i, e] * c_path[t, i, e]
        C_path[t] = c_sum

        pre = np.zeros((n_a, n_e))
        for i in range(n_a):
            for e in range(n_e):
                ap = a_pol_path[t, i, e]
                if ap < a_lo_grid:
                    ap = a_lo_grid
                elif ap > a_hi_grid:
                    ap = a_hi_grid
                # searchsorted(a_grid, ap, side='right'), clipped to [1, n_a-1]
                lo, hi = 0, n_a
                while lo < hi:
                    mid = (lo + hi) // 2
                    if a_grid[mid] <= ap:
                        lo = mid + 1
                    else:
                        hi = mid
                idx_hi = lo
                if idx_hi < 1:
                    idx_hi = 1
                elif idx_hi > n_a - 1:
                    idx_hi = n_a - 1
                idx_lo = idx_hi - 1
                denom = a_grid[idx_hi] - a_grid[idx_lo]
                if denom > 0.0:
                    w_hi = (ap - a_grid[idx_lo]) / denom
                else:
                    w_hi = 0.0
                mass = D[i, e]
                pre[idx_lo, e] += mass * (1.0 - w_hi)
                pre[idx_hi, e] += mass * w_hi

        D = pre @ Pi

        a_sum = 0.0
        for i in range(n_a):
            for e in range(n_e):
                a_sum += D[i, e] * a_grid[i]
        A_path[t] = a_sum

    for i in range(n_a):
        for e in range(n_e):
            D_start[T, i, e] = D[i, e]

    return A_path, C_path, D_start
