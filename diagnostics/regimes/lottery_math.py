"""
lottery_math.py — pure-numpy regime machinery on cached response matrices.

Conventions: A_def, A_cb are T x T maps from date-0-anticipated input paths
(shock_def_D, cb_buy_D) to the spread_rb deviation path. Stage A closed loop:
spread = (I - gamma*A_cb)^{-1} A_def @ eps. Stage B-lite (spec §10.2): the CB's
type s (gamma_s) is unknown until date k; agents hold beliefs pi; the CB is
silent pre-k. Time-invariance of the linearised model makes shift_k(A) the exact
response operator for news arriving at date k.

First-order certainty equivalence (no precautionary/Jensen term at this order)
lets branch solutions be mixed linearly under pi; time-invariance makes shift_k
the exact news operator *within the linearisation*. Both properties are
first-order facts of the model's linear IRFs and do not extend to a nonlinear
solve.
"""
import numpy as np

def shift_k(M, k):
    if k == 0:
        return M.copy()
    S = np.zeros_like(M)
    S[k:, k:] = M[:-k, :-k]
    return S

def closed_loop(A_def, A_cb, eps, gamma):
    T = len(eps)
    spread = np.linalg.solve(np.eye(T) - gamma * A_cb, A_def @ eps)
    return spread, gamma * spread

def peak(x, n=100):
    """Signed max over the first n periods (not abs-max) — correct for the
    positive spread response this module targets."""
    return float(np.asarray(x)[:n].max())

def gamma_for_compression(A_def, A_cb, eps, target, lo=0.0, hi=40.0, tol=1e-8):
    """Bisect for the gamma whose closed-loop peak spread is (1-target) x passive peak.
    Verifies monotonicity of peak(gamma) over a scan grid first (spec §14).

    hi=40 (2026-08-07; was 25 from 2026-07-31, 60 before that). Bisection's validity
    condition is monotonicity on the BRACKETING INTERVAL, not on an arbitrarily wide
    scan; hi is just the upper bracket, and the assert below is what tells you to move
    it.

    RAISED because the 50% target stopped bracketing on [0,25] after the country-size
    asymmetry and the rem_cb_F conduit fix (see docs/STATE.md). TPI is materially less
    effective than before both changes -- the same compression now needs ~8x the
    intervention:

        target      gamma before      gamma now
        25%              ~1.6            ~13
        50%              ~5.1            ~34

    Measured on cache_G_main_v3_psilam3p01 (passive peak 149.9bp):

        gamma      0     10     20     30     40   40.5     42    42.5
        peak bp  149.9  116.2   94.4   79.3   69.6   69.6  538.8    57.0

    peak(gamma) declines monotonically to 53.6% compression at gamma=40, turns at
    ~40.5, and there is a CLOSED-LOOP SINGULARITY at gamma ~ 42 (I - gamma*A_cb going
    near-singular; compute_tpi_irfs prints a conditioning warning there) with a
    separate branch beyond it. Do NOT set hi past 40 -- a coarse scan that samples 40
    and 50 reads the two branches as one smooth decline and hides the pole.

    NOTE FOR THE PAPER: the aggressive regime now sits ~8 gamma-units below that
    singularity, where before the fixes it was at gamma~5.1 with the pole far away.
    How much intervention the model can represent is now a live constraint, not a
    formality."""
    p0 = peak(closed_loop(A_def, A_cb, eps, 0.0)[0])
    grid = np.linspace(lo, hi, 61)
    peaks = np.array([peak(closed_loop(A_def, A_cb, eps, g)[0]) for g in grid])
    if not np.all(np.diff(peaks) < 0):
        bad = grid[1:][np.diff(peaks) >= 0]
        raise RuntimeError(f"peak spread not monotone in gamma near {bad[:5]} — "
                           f"bisection invalid (spec §14 monotonicity check)")
    f = lambda g: 1.0 - peak(closed_loop(A_def, A_cb, eps, g)[0]) / p0 - target
    assert f(lo) < 0 < f(hi), f"target {target} not bracketed on [{lo},{hi}]"
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return 0.5 * (lo + hi)

def solve_lottery(A_def, A_cb, eps, gammas, pi, k):
    """Spec §10.2 fixed point. Returns (spreads[n,T], cbs[n,T], cb_e[T]).

    Per type s:  spread^s = A_def eps + A_cb cb^e + shift_k(A_cb) (cb^s - cb^e)
    with cb^s = gamma_s * 1{t>=k} * spread^s  and  cb^e = sum_r pi_r cb^r.
    Jointly linear in the stacked branch spreads."""
    T, n = len(eps), len(gammas)
    pi = np.asarray(pi, dtype=float)
    assert abs(pi.sum() - 1.0) < 1e-12 and np.all(pi >= 0)
    assert 0 <= k < T, f"k={k} must satisfy 0 <= k < T={T} (k>=T silently degrades to open-loop)"
    Pi_k = np.diag((np.arange(T) >= k).astype(float))
    Sc = shift_k(A_cb, k)          # unanticipated-at-0 (news at k) component
    Dc = A_cb - Sc                 # anticipated-at-0 component
    b = A_def @ eps
    B = np.zeros((n * T, n * T))
    for s in range(n):
        for r in range(n):
            blk = -pi[r] * gammas[r] * (Dc @ Pi_k)
            if r == s:
                blk = blk + np.eye(T) - gammas[s] * (Sc @ Pi_k)
            B[s*T:(s+1)*T, r*T:(r+1)*T] = blk
    x = np.linalg.solve(B, np.tile(b, n))
    spreads = x.reshape(n, T)
    cbs = gammas[:, None] * (spreads @ Pi_k.T)
    cb_e = (pi[:, None] * cbs).sum(axis=0)
    return spreads, cbs, cb_e

def branch_output(M_def, M_cb, eps, cb_e, cb_s, k):
    """Any cached output's branch-s path: expected component + revision at k."""
    return M_def @ eps + M_cb @ cb_e + shift_k(M_cb, k) @ (cb_s - cb_e)
