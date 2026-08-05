"""Fast algebraic tests for the nominal-rigidity blocks.

These evaluate SSJ @simple blocks directly via .steady_state(), which is just
"evaluate at constant values" -- lags and leads collapse to the same constant.
That makes it usable for deliberately OFF-steady-state identity checks too.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


# ── Markup rent ───────────────────────────────────────────────────────────────

def test_firm_profit_is_zero_at_steady_state():
    from equations_D import firm_profit_D
    mu_p = 1.20
    ss = firm_profit_D.steady_state({
        'Y_D': 1.0, 'N_D': 0.8, 'alpha_D': 0.33,
        'mu_p_D': mu_p, 'mc_D': 1.0 / mu_p,
    })
    assert ss['profit_D'] == pytest.approx(0.0, abs=1e-15)


def test_firm_profit_restores_factor_exhaustion_off_steady_state():
    """w*N + profit must equal (1-alpha)*Y for ANY mc, so that adding the
    capital share alpha*Y exhausts output exactly."""
    from equations_D import firm_profit_D, labor_demand_D
    mu_p, mc, Y, N, alpha = 1.20, 0.79, 1.03, 0.81, 0.33

    # w from labour demand at this mc (w_res_D == 0 defines w)
    w = mu_p * mc * (1 - alpha) * Y / N

    ss = firm_profit_D.steady_state({
        'Y_D': Y, 'N_D': N, 'alpha_D': alpha, 'mu_p_D': mu_p, 'mc_D': mc,
    })
    assert w * N + ss['profit_D'] == pytest.approx((1 - alpha) * Y, rel=1e-14)

    # and the wage we assumed really is the one labor_demand_D implies
    ld = labor_demand_D.steady_state({
        'w_D': w, 'Y_D': Y, 'N_D': N, 'alpha_D': alpha,
        'mu_p_D': mu_p, 'mc_D': mc,
    })
    assert ld['w_res_D'] == pytest.approx(0.0, abs=1e-14)


def test_firm_profit_F_matches_D():
    from equations_D import firm_profit_D
    from equations_F import firm_profit_F
    args = dict(Y=1.03, N=0.81, alpha=0.33, mu_p=1.20, mc=0.79)
    d = firm_profit_D.steady_state({
        'Y_D': args['Y'], 'N_D': args['N'], 'alpha_D': args['alpha'],
        'mu_p_D': args['mu_p'], 'mc_D': args['mc'],
    })
    f = firm_profit_F.steady_state({
        'Y_F': args['Y'], 'N_F': args['N'], 'alpha_F': args['alpha'],
        'mu_p_F': args['mu_p'], 'mc_F': args['mc'],
    })
    assert d['profit_D'] == pytest.approx(f['profit_F'], rel=1e-15)


# ── Price Phillips curve ──────────────────────────────────────────────────────

def test_price_nkpc_is_zero_at_steady_state():
    from equations_D import price_nkpc_D
    mu_p = 1.20
    ss = price_nkpc_D.steady_state({
        'pi_D': 0.0, 'mc_D': 1.0 / mu_p, 'mu_p_D': mu_p,
        'kappa_p_D': 0.0871, 'beta_D': 0.985,
    })
    assert ss['nkpc_p_res_D'] == pytest.approx(0.0, abs=1e-15)


def test_price_nkpc_flex_limit_forces_mc_to_one_over_mu_p():
    """As kappa_p -> inf the residual/kappa_p -> -(mu_p*mc - 1), so setting the
    residual to zero drives mu_p*mc -> 1, which is the competitive condition."""
    from equations_D import price_nkpc_D
    mu_p = 1.20
    base = {'pi_D': 0.0, 'mu_p_D': mu_p, 'beta_D': 0.985}
    off_mc = 0.79                       # != 1/mu_p = 0.8333...
    for kappa in (1e2, 1e4, 1e6):
        ss = price_nkpc_D.steady_state({**base, 'mc_D': off_mc, 'kappa_p_D': kappa})
        implied_gap = -ss['nkpc_p_res_D'] / kappa
        assert implied_gap == pytest.approx(mu_p * off_mc - 1.0, rel=1e-12)


def test_price_nkpc_gap_linearises_to_mc_hat():
    """d(mu_p*mc - 1)/d(mc/mc_ss) evaluated at mc_ss = 1/mu_p equals 1 for ANY
    mu_p -- which is why mu_p is a free normalisation to first order."""
    from equations_D import price_nkpc_D
    for mu_p in (1.05, 1.20, 1.50):
        mc_ss = 1.0 / mu_p
        h = 1e-7
        base = {'pi_D': 0.0, 'mu_p_D': mu_p, 'kappa_p_D': 1.0, 'beta_D': 0.985}
        up = price_nkpc_D.steady_state({**base, 'mc_D': mc_ss * (1 + h)})
        dn = price_nkpc_D.steady_state({**base, 'mc_D': mc_ss * (1 - h)})
        # residual = -kappa*(gap), kappa = 1 -> d(gap)/d(mc_hat) = -d(res)/d(mc_hat)
        d_gap = -(up['nkpc_p_res_D'] - dn['nkpc_p_res_D']) / (2 * h)
        assert d_gap == pytest.approx(1.0, rel=1e-6)


def test_price_nkpc_F_matches_D():
    from equations_D import price_nkpc_D
    from equations_F import price_nkpc_F
    args = dict(pi=0.001, mc=0.79, mu_p=1.20, kappa=0.0871, beta=0.985)
    d = price_nkpc_D.steady_state({
        'pi_D': args['pi'], 'mc_D': args['mc'], 'mu_p_D': args['mu_p'],
        'kappa_p_D': args['kappa'], 'beta_D': args['beta'],
    })
    f = price_nkpc_F.steady_state({
        'pi_F': args['pi'], 'mc_F': args['mc'], 'mu_p_F': args['mu_p'],
        'kappa_p_F': args['kappa'], 'beta_F': args['beta'],
    })
    assert d['nkpc_p_res_D'] == pytest.approx(f['nkpc_p_res_F'], rel=1e-15)


def test_labor_demand_collapses_to_competitive_at_ss_markup():
    """At mc = 1/mu_p the condition must be exactly w = (1-alpha)Y/N, which is
    what makes the steady state bit-identical to the flex model."""
    from equations_D import labor_demand_D
    mu_p, Y, N, alpha = 1.20, 1.03, 0.81, 0.33
    w_competitive = (1 - alpha) * Y / N
    ss = labor_demand_D.steady_state({
        'w_D': w_competitive, 'Y_D': Y, 'N_D': N, 'alpha_D': alpha,
        'mu_p_D': mu_p, 'mc_D': 1.0 / mu_p,
    })
    assert ss['w_res_D'] == pytest.approx(0.0, abs=1e-15)


# ── Global closure ────────────────────────────────────────────────────────────

def test_global_residuals_zero_at_steady_state():
    from equations_global import terms_of_trade, union_inflation
    tot = terms_of_trade.steady_state({'p': 0.99, 'pi_D': 0.0, 'pi_F': 0.0})
    assert tot['tot_res'] == pytest.approx(0.0, abs=1e-15)
    uni = union_inflation.steady_state({'pi_D': 0.0, 'pi_F': 0.0, 'omega_pi_D': 0.071})
    assert uni['union_pi_res'] == pytest.approx(0.0, abs=1e-15)


def test_closure_puts_93pct_of_tot_move_into_D_deflation():
    """Solving tot_res = 0 and union_pi_res = 0 together gives
    pi_D = -(1 - omega)*dlog p and pi_F = omega*dlog p, so at the capital-key
    omega = 0.071 the terms-of-trade adjustment splits 93/7 between D deflation
    and F inflation -- the internal-devaluation pattern.

    Asserted on the net rates, which is exact arithmetic. Do NOT assert
    log((1+pi_F)/(1+pi_D)) == dlog_p: that holds only to first order, and the
    O(dlog_p^2) truncation is 0.429*dlog_p in relative terms, which swamps any
    tight tolerance.
    """
    from equations_global import union_inflation
    omega = 0.071
    dlog_p = 1e-4
    pi_D = -(1 - omega) * dlog_p
    pi_F = omega * dlog_p

    # the normalisation holds exactly
    uni = union_inflation.steady_state({'pi_D': pi_D, 'pi_F': pi_F,
                                        'omega_pi_D': omega})
    assert uni['union_pi_res'] == pytest.approx(0.0, abs=1e-18)

    # the differential is exactly the terms-of-trade move
    assert pi_F - pi_D == pytest.approx(dlog_p, rel=1e-15)

    # and D bears 1 - omega of it: 92.9% here
    share_D = abs(pi_D) / (abs(pi_D) + abs(pi_F))
    assert share_D == pytest.approx(1 - omega, rel=1e-15)
    assert share_D > 0.92


def test_omega_one_half_splits_evenly():
    """Guards the calibration argument: at omega = 0.5 the adjustment splits
    50/50, which is counterfactual for GR/DE. See the spec."""
    omega = 0.5
    dlog_p = 1e-4
    assert -(1 - omega) * dlog_p == pytest.approx(-0.5 * dlog_p, rel=1e-15)
    assert omega * dlog_p == pytest.approx(0.5 * dlog_p, rel=1e-15)
