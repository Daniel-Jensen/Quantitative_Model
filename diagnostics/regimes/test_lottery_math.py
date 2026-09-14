"""Unit tests for lottery_math on small synthetic matrices — no model solve, <1s."""
import numpy as np
import pytest
from lottery_math import shift_k, closed_loop, gamma_for_compression, solve_lottery, branch_output

T = 40
rng = np.random.default_rng(7)
A_def = 0.5 * np.tril(rng.standard_normal((T, T))) / T + 0.3 * np.eye(T)
A_cb  = -(0.4 * rng.standard_normal((T, T)) / T + 0.2 * np.eye(T))  # purchases compress spreads
eps = np.zeros(T); eps[0] = 1.0; eps = 0.8 ** np.arange(T) * eps[0]

def test_shift_k_zero_is_identity():
    assert np.array_equal(shift_k(A_cb, 0), A_cb)

def test_shift_k_structure():
    S = shift_k(A_cb, 3)
    assert np.all(S[:3, :] == 0) and np.all(S[:, :3] == 0)
    assert np.array_equal(S[3:, 3:], A_cb[:-3, :-3])

def test_closed_loop_gamma0_is_open_loop():
    sp, cb = closed_loop(A_def, A_cb, eps, 0.0)
    assert np.allclose(sp, A_def @ eps) and np.allclose(cb, 0)

def test_compression_bisection():
    g = gamma_for_compression(A_def, A_cb, eps, target=0.5)
    sp0, _ = closed_loop(A_def, A_cb, eps, 0.0)
    spg, _ = closed_loop(A_def, A_cb, eps, g)
    assert abs(1.0 - spg[:100].max() / sp0[:100].max() - 0.5) < 1e-6

def test_lottery_k0_nests_stage_a_any_pi():
    gammas = np.array([3.0, 1.0, 0.0])
    for pi in ([0.2, 0.3, 0.5], [1.0, 0.0, 0.0], [1/3, 1/3, 1/3]):
        spreads, cbs, cb_e = solve_lottery(A_def, A_cb, eps, gammas, np.array(pi), k=0)
        for s, g in enumerate(gammas):
            sp_known, _ = closed_loop(A_def, A_cb, eps, g)
            assert np.allclose(spreads[s], sp_known, atol=1e-9), f"k=0 nesting failed s={s} pi={pi}"

def test_lottery_degenerate_pi_equals_delayed_known():
    gammas = np.array([3.0, 1.0, 0.0]); k = 4
    for s in range(3):
        pi = np.zeros(3); pi[s] = 1.0
        spreads, cbs, cb_e = solve_lottery(A_def, A_cb, eps, gammas, pi, k)
        Pi_k = np.diag((np.arange(T) >= k).astype(float))
        sp_delayed = np.linalg.solve(np.eye(T) - gammas[s] * A_cb @ Pi_k, A_def @ eps)
        assert np.allclose(spreads[s], sp_delayed, atol=1e-9)

def test_lottery_pre_k_branch_identity_and_re_jump():
    gammas = np.array([3.0, 1.0, 0.0]); pi = np.array([0.25, 0.35, 0.40]); k = 5
    spreads, cbs, cb_e = solve_lottery(A_def, A_cb, eps, gammas, pi, k)
    assert np.max(np.abs(spreads[:, :k] - spreads[0, :k])) < 1e-10, "pre-k branch identity"
    assert np.max(np.abs((pi[:, None] * (cbs - cb_e)).sum(0))) < 1e-10, "E[revision]=0"
    y_jump = np.stack([shift_k(A_cb, k) @ (cbs[s] - cb_e) for s in range(3)])
    assert np.max(np.abs((pi[:, None] * y_jump).sum(0))) < 1e-10, "pi-weighted output jump != 0"

def test_defining_equation_residual_k_positive():
    gammas = np.array([3.0, 1.0, 0.0]); pi = np.array([0.25, 0.35, 0.40]); k = 5
    spreads, cbs, cb_e = solve_lottery(A_def, A_cb, eps, gammas, pi, k)
    for s in range(3):
        rhs = branch_output(A_def, A_cb, eps, cb_e, cbs[s], k)
        assert np.max(np.abs(spreads[s] - rhs)) < 1e-9, f"defining-equation residual failed s={s}"
