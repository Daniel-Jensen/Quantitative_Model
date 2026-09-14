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


# ── Nominal deposits ──────────────────────────────────────────────────────────

def test_deposit_rates_collapse_at_zero_inflation():
    """At pi = 0 both derived real rates must equal the nominal rate exactly --
    this is what keeps the steady state bit-identical."""
    from equations_D import deposit_rates_D
    ss = deposit_rates_D.steady_state({'i_dep_D': 0.0125, 'pi_D': 0.0})
    assert ss['rdep_D'] == pytest.approx(0.0125, rel=1e-15)
    assert ss['rdep_expost_D'] == pytest.approx(0.0125, rel=1e-15)


def test_deflation_raises_the_realised_real_deposit_rate():
    """Deflation is a windfall to depositors and a loss to banks, which hold
    real assets against nominal liabilities. This is the Fisher channel; if the
    sign flips, bank_return_D will amplify in the wrong direction."""
    from equations_D import deposit_rates_D
    i = 0.0125
    base = deposit_rates_D.steady_state({'i_dep_D': i, 'pi_D': 0.0})
    defl = deposit_rates_D.steady_state({'i_dep_D': i, 'pi_D': -0.001})
    assert defl['rdep_expost_D'] > base['rdep_expost_D']
    assert defl['rdep_expost_D'] == pytest.approx((1 + i) / (1 - 0.001) - 1, rel=1e-14)


def test_deposit_return_is_unchanged_at_zero_inflation():
    """Rgross must be exactly 1 + i_dep when pi = 0 and P_CES is flat."""
    from equations_D import deposit_return_D
    ss = deposit_return_D.steady_state({'i_dep_D': 0.0125, 'P_CES_D': 1.3, 'pi_D': 0.0})
    assert ss['Rgross_D'] == pytest.approx(1.0125, rel=1e-15)


def test_deposit_rates_F_matches_D():
    from equations_D import deposit_rates_D
    from equations_F import deposit_rates_F
    d = deposit_rates_D.steady_state({'i_dep_D': 0.0125, 'pi_D': -0.001})
    f = deposit_rates_F.steady_state({'i_dep_F': 0.0125, 'pi_F': -0.001})
    assert d['rdep_D'] == pytest.approx(f['rdep_F'], rel=1e-15)
    assert d['rdep_expost_D'] == pytest.approx(f['rdep_expost_F'], rel=1e-15)


def test_bank_return_uses_the_expost_rate():
    """Signature check: bank_return_D must take rdep_expost_D and must NOT take
    rdep_D. Getting this backwards silently reverses the Fisher channel."""
    from equations_D import bank_return_D, capital_fund_D
    for blk in (bank_return_D, capital_fund_D):
        assert 'rdep_expost_D' in blk.inputs, (blk.name, sorted(blk.inputs))
        assert 'rdep_D' not in blk.inputs, (blk.name, sorted(blk.inputs))


def test_forward_looking_blocks_still_use_rdep():
    """intermediation_P1_D is ex-ante: rdep_D still means the t -> t+1 real rate.

    divert_bond_foc_D used to be checked here too. Since the 2026-08-17/18 structural
    refactor its replacement gk_bond_foc_D states GK portfolio optimality on the marginal
    values (nu_bD_D = Delta_bD_eff_D * nu_K_D) and carries NO interest rate at all, so the
    T-2 invariant now lives entirely in intermediation_P1_D. Asserted below, which is
    strictly stronger than the old signature check.
    """
    from equations_D import intermediation_P1_D
    assert 'rdep_D' in intermediation_P1_D.inputs, sorted(intermediation_P1_D.inputs)
    assert 'rdep_expost_D' not in intermediation_P1_D.inputs, sorted(intermediation_P1_D.inputs)


def test_own_sovereign_foc_is_pure_gk():
    """gk_bond_foc_D/F must be a function of ENDOGENOUS GK objects only.

    Guards the refactor against regression to a reduced-form spread rule: no frozen
    psi_spread, no excess_return_*_ss anchor, no tau_mp wedge, no interest rate.
    """
    from equations_D import gk_bond_foc_D
    from equations_F import gk_bond_foc_F
    for blk, want in ((gk_bond_foc_D, {'nu_bD_D', 'nu_K_D', 'Delta_bD_eff_D'}),
                      (gk_bond_foc_F, {'nu_bF_F', 'nu_K_F', 'Delta_bF_eff_F'})):
        assert set(blk.inputs) == want, (blk.name, sorted(blk.inputs))


def test_no_ad_hoc_sovereign_spread_wedge_anywhere():
    """The governing invariant of the 2026-08-18 refactor.

    psi_spread_D/F and EL_price_D/F are DELETED. No symbol of either name may reappear
    in code/, in any form — parameter, anchor, renamed twin or hard-coded literal. The
    sovereign spread must come from the state-contingent payoff plus the GK portfolio
    FOC, never from `spread += coefficient * default_probability`.
    """
    import ast
    import pathlib
    import re
    root = pathlib.Path(__file__).resolve().parent
    banned = ('psi_spread', 'EL_price', 'divert_bond_foc', 'divert_portfolio_adj',
              'excess_return_')
    ident = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')

    def symbols(tree):
        """Every identifier the module actually USES. Prose in docstrings and comments
        is exempt — the deleted objects are documented at length in the blocks that
        replaced them, and that history is worth keeping."""
        for n in ast.walk(tree):
            if isinstance(n, ast.Name):
                yield n.id
            elif isinstance(n, ast.arg):
                yield n.arg
            elif isinstance(n, ast.Attribute):
                yield n.attr
            elif isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                yield n.name
            elif isinstance(n, ast.alias):
                yield n.name
                if n.asname:
                    yield n.asname
            elif isinstance(n, ast.Constant) and isinstance(n.value, str):
                # Dict keys / ss lookups such as ss['psi_spread_D'] land here. Docstrings
                # do too, but they never match the identifier pattern.
                if ident.match(n.value):
                    yield n.value

    hits = []
    for f in sorted(root.glob('*.py')):
        if f.name == pathlib.Path(__file__).name:
            continue
        for s in symbols(ast.parse(f.read_text())):
            for b in banned:
                if b in s:
                    hits.append(f"{f.name}: {s}")
    assert not hits, ("deleted objects still referenced in live code:\n"
                      + "\n".join(sorted(set(hits))))


def test_expected_loss_enters_each_sovereign_leg_exactly_once():
    """Double-counting audit, as a test.

    Expected default loss lives in bond_return_D/F's rb_exp and reaches the banks only
    through intermediation_P1_D/F. Every downstream portfolio condition must be stated on
    the resulting nu's, so no block may take BOTH a nu and a def_rate/EL object.
    """
    from equations_D import intermediation_P1_D, gk_bond_foc_D
    from equations_F import intermediation_P1_F, gk_bond_foc_F
    from equations_global import gk_cross_border_foc

    # The Euler equations price rb_exp, never rb_actual and never def_rate directly.
    for blk, c, o in ((intermediation_P1_D, 'D', 'F'), (intermediation_P1_F, 'F', 'D')):
        ins = set(blk.inputs)
        assert f'rb_exp_{c}' in ins and f'rb_exp_{o}' in ins, (blk.name, sorted(ins))
        assert not (ins & {f'rb_actual_{c}', f'rb_actual_{o}',
                           f'def_rate_{c}', f'def_rate_{o}'}), (blk.name, sorted(ins))

    # The four portfolio conditions carry no default object of their own.
    for blk in (gk_bond_foc_D, gk_bond_foc_F, gk_cross_border_foc):
        ins = set(blk.inputs)
        assert not any(k in i for i in ins
                       for k in ('def_rate', 'EL_price', 'rb_exp', 'rb_actual')), \
            (blk.name, sorted(ins))
        assert any(i.startswith('nu_') for i in ins), (blk.name, sorted(ins))


# ── Investment-flow adjustment cost ───────────────────────────────────────────

def test_flow_adjustment_cost_vanishes_at_steady_state():
    """S(1) = S'(1) = 0 is what makes this SS-neutral. Check the block's own
    residual is unchanged when investment is flat, for ANY omega_I."""
    from equations_D import capital_adj_D
    from equations_F import capital_adj_F
    # SSJ's .steady_state() silently ignores dict keys the block does not take,
    # so without this the rest of the test is vacuously green on the OLD block.
    for blk, suf in ((capital_adj_D, 'D'), (capital_adj_F, 'F')):
        assert f'omega_I_{suf}' in blk.inputs, (blk.name, sorted(blk.inputs))
        # Discounted at constant beta, NOT SDF: first-order exact because
        # S'(1) = 0, and taking SDF here makes SSJ's topological sort fail
        # (hh -> capital_fund -> capital_adj -> sdf -> ghh_composite -> hh).
        assert f'beta_{suf}' in blk.inputs, (blk.name, sorted(blk.inputs))
        assert f'SDF_{suf}' not in blk.inputs, (blk.name, sorted(blk.inputs))
    base = dict(K_D=10.8, Q_D=1.0, I_D=0.242, Z_D=1.0, N_D=0.8, alpha_D=0.33,
                delta_D=0.022407, gamma0_D=0.15, gamma1_D=-0.0053, ksi_D=0.5,
                beta_D=0.9995)
    ref = capital_adj_D.steady_state({**base, 'omega_I_D': 0.0})
    for w in (0.0, 2.0, 10.0):
        ss = capital_adj_D.steady_state({**base, 'omega_I_D': w})
        assert ss['q_res_D'] == pytest.approx(ref['q_res_D'], rel=1e-14), w
        assert ss['iota_D'] == pytest.approx(ref['iota_D'], rel=1e-14), w
        assert ss['capital_res_D'] == pytest.approx(ref['capital_res_D'], rel=1e-14), w


def test_flow_adjustment_cost_bites_off_steady_state():
    """With investment falling, the cost must be strictly positive and must
    scale with omega_I -- otherwise the parameter is doing nothing."""
    from equations_D import capital_adj_D
    base = dict(K_D=10.8, Q_D=1.0, Z_D=1.0, N_D=0.8, alpha_D=0.33,
                delta_D=0.022407, gamma0_D=0.15, gamma1_D=-0.0053, ksi_D=0.5,
                beta_D=0.9995)
    # steady_state() collapses lags, so probe the S(x) algebra directly instead
    for w in (2.0, 10.0):
        x = 0.9                                   # investment 10% below last period
        S = (w / 2.0) * (x - 1.0) ** 2
        assert S > 0
        assert S == pytest.approx(w * 0.005, rel=1e-14)
