# NUMBA JIT KERNELS FOR THE TWO HOT LOOPS (EGM BACKWARD + DISTRIBUTION FORWARD).
# These replicate the numpy reference ops-for-op (equivalence is tested in
# test_fast_kernels.py); callers fall back to numpy when numba is absent or
# cal["use_numba"]=False. cache=True lets spawned Jacobian workers load artifacts.
import numpy as np

try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:                                   # pragma: no cover
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        # NO-OP STAND-IN FOR numba.njit WHEN NUMBA IS UNAVAILABLE.
        def wrap(f):
            # RETURN THE FUNCTION UNCOMPILED.
            return f
        return wrap


@njit(cache=True)
def hh_backward(a_grid, Pi_T, r_path, y_path, c_terminal, beta, sigma, a_min,
                vN_path):
    # EGM BACKWARD INDUCTION OVER THE PATH (MIRRORS household.egm_step EXACTLY).
    # Pi_T = Pi.T; r_path has length T+1 (Fisher-adjusted upstream).
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

        x_next = np.maximum(c_next - vN_next, 1e-11)             # GHH composite tomorrow
        Eu_next = (x_next ** (-sigma)) @ Pi_T                     # (n_a, n_e)

        x_endo = (beta * (1.0 + r_next) * Eu_next) ** (-1.0 / sigma)   # Euler inversion
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
    # DISTRIBUTION FORWARD SIMULATION WITH YOUNG (2010) LOTTERY WEIGHTS.
    # C_t over the start-of-period dist, A_t over the end-of-period one;
    # D_start[t] is the dist entering period t.
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
