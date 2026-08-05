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
