import copy
import numpy as np
import sequence_jacobian as sj
from sequence_jacobian import simple, combine, create_model

from equations_D import (
    hh_init_D, hh_D, make_grids_D, income_D, hh_extended_D,
    smart_steady_D, market_clearing_D, steady_auxilliary_D,
    banker_div_D, sdf_D, sdf_ss_D, sdf_banker_ss_D, government_ss_D, labor_ss_D,
    government_default_D, bond_return_D, collateral_quality_D, gk_bond_foc_D,
    ces_price_D, import_demand_D, deposit_rates_D, deposit_return_D,
    firm_profit_D, price_nkpc_D,
)
from equations_F import (
    hh_init_F, hh_F, make_grids_F, income_F, hh_extended_F,
    smart_steady_F, market_clearing_F, steady_auxilliary_F,
    banker_div_F, sdf_F, sdf_ss_F, sdf_banker_ss_F, government_ss_F, labor_ss_F,
    government_default_F, bond_return_F, collateral_quality_F, gk_bond_foc_F,
    ces_price_F, import_demand_F, deposit_rates_F, deposit_return_F,
    firm_profit_F, price_nkpc_F,
)
from equations_global import (
    trade_balance, domestic_bond_clearing,
    portfolio_level_anchors, bond_yield,
    global_goods_mkt, external_account_D,
    terms_of_trade, union_inflation,
)
# NB: import the MODULE, not the flag. `from calibration import EBA_CALIBRATION`
# binds the value at import time, so a sweep that flips the switch afterwards
# would silently keep the old portfolio targets — the same stale-binding trap as
# the regimes cache key and PSILAM_MAIN.
import calibration
from calibration import load_eba_targets, load_eba_foreign_shares


def assert_gk_well_posed(ss_in):
    """Hard guard: the Gertler-Karadi block must be economically well-posed.

    `lambda_gk` is the IC multiplier and `Omega` the banker's marginal value of
    net worth; both must be strictly positive. When they are not, the solver
    still converges, every Walras residual is machine-zero, and the IRFs look
    plausible — but the banker's continuation value is negative and nothing
    computed from the block means anything. This is the C-1 failure mode:
    silent degeneracy that passes every check the pipeline previously ran.

    From `steady_auxilliary_D/F`,
        lambda_gk = f / (D_target/(beta_inter*(1+rn)) - (1-f)*theta),
        D_target  = theta - (1-Delta_own)*phi_own - (1-Delta_cross)*phi_cross,
    so lambda_gk > 0 requires a positive denominator, i.e. approximately

        f * theta  >  (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross.

    Read: the banker's franchise value (left) must cover the non-divertable
    ("good collateral") part of the sovereign book (right). High measured
    concentration `phi_own` therefore puts a LOWER BOUND on the divertability
    `Delta_own` — which is how the EBA data partially identifies a parameter
    that has no direct empirical counterpart. See docs/eba_calibration.md.
    """
    vals, bad = [], []
    for c in ("D", "F"):
        lam, om = float(ss_in[f"lambda_gk_{c}"]), float(ss_in[f"Omega_{c}"])
        vals.append(f"{c}: lambda_gk={lam:+.6f}, Omega={om:+.6f}")
        if not (lam > 0.0 and om > 0.0):
            bad.append(c)
    if bad:
        raise ValueError(
            f"GK block is NOT well-posed — negative IC multiplier / franchise "
            f"value in {'/'.join(bad)}:\n  " + "\n  ".join(vals)
            + "\n\nThe steady state will still 'solve' and every Walras residual will be "
              "machine-zero, but the banker's continuation value is negative and all "
              "IRFs are uninterpretable. Required (approximately):\n"
              "    f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross\n"
              "Raise Delta_own (bond divertability), raise f, or lower the sovereign "
              "concentration. See assert_gk_well_posed.__doc__ and "
              "docs/eba_calibration.md 'GK feasibility'.")


def gk_feasibility_margin(theta, f, phi_own, phi_cross, Delta_own, Delta_cross):
    """f*theta - [(1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross]; must be > 0."""
    return f * theta - ((1.0 - Delta_own) * phi_own + (1.0 - Delta_cross) * phi_cross)


def min_Delta_own(theta, f, phi_own, phi_cross, Delta_cross, margin=0.0):
    """Smallest own-bond divertability consistent with a well-posed GK block."""
    return 1.0 - (f * theta - (1.0 - Delta_cross) * phi_cross - margin) / phi_own


def gk_cross_wedges(ss_in):
    """SS shadow cost of the calibrated cross-border sovereign positions, in return units.

    ``(nu_cross/nu_K - Delta_cross_eff) * (rk - rdep)``, i.e. how far the EBA-measured
    cross-border book sits from GK proportionality at the solved steady state. Enters
    ``gk_cross_border_foc`` as a CONSTANT — it carries no ``def_rate`` and is not a spread
    loading; it is the level at which the portfolio adjustment cost sits so that the
    measured position is optimal.

    At a riskless steady state (``def_rate_ss = 0``) with ``rk_D = rk_F`` and
    ``rdep_D = rdep_F``, GK optimality forces the cross divertability to equal the own
    one, and these wedges vanish. They are printed and range-checked on every solve so an
    inconsistent ``Delta_bF_D``/``Delta_bD_F`` cannot hide inside them — which is exactly
    what the pre-refactor ``excess_return_F_D_ss`` / ``excess_return_D_F_ss`` anchors did,
    since those absorbed the FULL excess return rather than its deviation from the FOC.
    """
    exc_D = float(ss_in['rk_D']) - float(ss_in['rdep_D'])
    exc_F = float(ss_in['rk_F']) - float(ss_in['rdep_F'])
    w_FD = (float(ss_in['nu_bF_D']) / float(ss_in['nu_K_D'])
            - float(ss_in['Delta_bF_eff_D'])) * exc_D
    w_DF = (float(ss_in['nu_bD_F']) / float(ss_in['nu_K_F'])
            - float(ss_in['Delta_bD_eff_F'])) * exc_F
    return w_FD, w_DF


def report_gk_steady_state(ss_in, cal, tol_own=1e-9, tol_cross=5e-4):
    """Print the GK/sovereign steady state and NUMERICALLY VERIFY the portfolio FOCs.

    The own-sovereign conditions ``nu_bD_D/nu_K_D = Delta_bD_eff_D`` and
    ``nu_bF_F/nu_K_F = Delta_bF_eff_F`` are SS TARGETS, so they must hold to solver
    tolerance; a failure means ``gk_bond_foc_D/F`` is no longer wired to ``q_b_D/q_b_F``.
    The cross-border conditions are pinned by quantities the calibration takes from EBA,
    so they hold only up to ``gk_wedge_*_ss``; that wedge is reported and range-checked
    (``tol_cross``, on the dimensionless ``nu_cross/nu_K - Delta_cross_eff``) so an
    inconsistent cross divertability cannot hide inside it. At the live calibration all
    four Delta are 0.20 and the cross residual is ~2e-13.
    """
    # Quarterly fraction -> annualised basis points. NOT 400: that is the *percentage-point*
    # scaling, and mixing the two is how a 205bp spread gets reported as 2bp.
    BP_ANN = 4.0e4
    g = lambda k: float(ss_in[k])
    rows = []
    for c, own, cross in (("D", "bD", "bF"), ("F", "bF", "bD")):
        r_own = g(f"nu_{own}_{c}") / g(f"nu_K_{c}")
        r_cr = g(f"nu_{cross}_{c}") / g(f"nu_K_{c}")
        rows.append((c, r_own, g(f"Delta_{own}_eff_{c}"), r_cr, g(f"Delta_{cross}_eff_{c}")))

    print("\n=== GK steady state: sovereign block ===")
    print(f"  {'':4} {'q_b':>9} {'yield q/q':>10} {'yield ann':>11} {'def_rate':>9} "
          f"{'recovery':>9} {'EL_load':>9} {'delta_b':>8} {'duration':>9}")
    for c in ("D", "F"):
        q, db = g(f"q_b_{c}"), g(f"delta_b_{c}")
        y = db * (1.0 / q - 1.0)
        print(f"  {c:4} {q:>9.6f} {y:>10.6f} {y*BP_ANN:>9.2f}bp {g(f'def_rate_{c}'):>9.4f} "
              f"{cal[f'recovery_rate_{c}']:>9.2f} {g(f'EL_load_{c}'):>9.6f} "
              f"{db:>8.4f} {1.0/db:>8.1f}q")
    sp = g("delta_b_D") * (1.0 / g("q_b_D") - 1.0) - g("delta_b_F") * (1.0 / g("q_b_F") - 1.0)
    print(f"  SS spread (D-F) = {sp*BP_ANN:+.4f} bp annualised "
          f"(zero is correct: def_rate_ss = 0 in both countries)")

    print(f"\n  {'':4} {'lambda_gk':>10} {'Omega':>9} {'theta':>8} {'n_inter':>9} "
          f"{'K':>9} {'q_b*b_own':>10} {'q_b*b_cross':>12}")
    for c in ("D", "F"):
        own_b = "b_D_D" if c == "D" else "b_F_F"
        cr_b = "b_F_D" if c == "D" else "b_D_F"
        own_q = "q_b_D" if c == "D" else "q_b_F"
        cr_q = "q_b_F" if c == "D" else "q_b_D"
        print(f"  {c:4} {g(f'lambda_gk_{c}'):>10.4f} {g(f'Omega_{c}'):>9.4f} "
              f"{g(f'theta_{c}'):>8.4f} {g(f'n_inter_{c}'):>9.4f} {g(f'K_{c}'):>9.4f} "
              f"{g(own_q)*g(own_b):>10.4f} {g(cr_q)*g(cr_b):>12.4f}")

    print("\n  Portfolio FOC check   nu_i/nu_K  vs  Delta_i_eff")
    print(f"  {'bank':5} {'leg':7} {'nu_i/nu_K':>11} {'Delta_eff':>11} {'residual':>12} {'status':>8}")
    ok = True
    for c, r_own, D_own, r_cr, D_cr in rows:
        for leg, r, D, tol in (("own", r_own, D_own, tol_own), ("cross", r_cr, D_cr, tol_cross)):
            res = r - D
            good = abs(res) <= tol
            ok &= good
            print(f"  {c:5} {leg:7} {r:>11.6f} {D:>11.6f} {res:>12.3e} "
                  f"{'OK' if good else 'FAIL':>8}")
    # Over-identifying check the spec asks for explicitly: nu_own/nu_cross must equal
    # Delta_own/Delta_cross. Implied by the two rows above, printed because it is the
    # form the portfolio condition is usually stated in.
    for c, own, cross in (("D", "bD", "bF"), ("F", "bF", "bD")):
        nu_r = g(f"nu_{own}_{c}") / g(f"nu_{cross}_{c}")
        D_r = g(f"Delta_{own}_eff_{c}") / g(f"Delta_{cross}_eff_{c}")
        print(f"  {c:5} {'ratio':7} {nu_r:>11.6f} {D_r:>11.6f} {nu_r - D_r:>12.3e}")
    w_FD, w_DF = gk_cross_wedges(ss_in)
    print(f"  cross-border SS wedges (constant, def_rate-free):"
          f"  F-in-D = {w_FD*BP_ANN:+.3e} bp/yr   D-in-F = {w_DF*BP_ANN:+.3e} bp/yr")
    for c in ("D", "F"):
        for leg in ("bD", "bF"):
            v = float(ss_in[f"Delta_{leg}_eff_{c}"])
            if not (0.0 <= v <= 1.0):
                raise ValueError(f"Delta_{leg}_eff_{c} = {v} outside [0,1] at the steady state.")
    if not ok:
        raise ValueError(
            "GK portfolio FOC violated at the steady state. The own-sovereign legs are "
            "SS targets and must hold to solver tolerance; a cross-border failure means "
            "the cross divertability Delta is inconsistent with the calibrated "
            "cross-border position by more than gk_cross_wedges can legitimately absorb. "
            "Do NOT paper over this with a spread parameter — fix Delta or the position.")
    print("  psi_lambda_B = "
          f"{cal['psi_lambda_B_D']:.4f}/{cal['psi_lambda_B_F']:.4f}  "
          f"zeta_writeoff = {cal['zeta_writeoff_D']:.1f}/{cal['zeta_writeoff_F']:.1f}  "
          f"writeoff_enabled = {cal['writeoff_enabled_D']:.1f}/{cal['writeoff_enabled_F']:.1f}")


def _apply_ss_anchors(ss_in, cal):
    # Fires on every solved SS (this is the common path for steady_state.py and
    # depreciation_calibration.py).
    #
    # DELETED 2026-08-18 (structural GK refactor):
    #   psi_spread_D/F           -- the free sovereign-spread loading. Deleted, not
    #                               recalibrated: it existed only to absorb the
    #                               principal/continuation loss that zeta_writeoff = 0
    #                               left out of the bond payoff.
    #   EL_price_D/F             -- separate expected-loss PRICING wedge. The loss is now
    #                               inside the payoff itself (bond_return_D/F -> rb_exp).
    #                               The diagnostic loading survives as the endogenous
    #                               EL_load_D/F, an OUTPUT of bond_return_D/F.
    #   excess_return_bD_D_ss    -- own-sovereign SS excess-return anchors, consumed only
    #   excess_return_bF_F_ss       by the deleted domestic_bond_foc_D/F.
    #   excess_return_F_D_ss     -- cross-border anchors that absorbed the whole SS excess
    #   excess_return_D_F_ss        return. Replaced by gk_cross_wedges, which absorbs
    #                               only the DEVIATION from the GK portfolio FOC.
    assert_gk_well_posed(ss_in)
    w_FD, w_DF = gk_cross_wedges(ss_in)
    anchors = {
        'phi_bD_D_ss':           float(ss_in['q_b_D']) * float(ss_in['b_D_D']) / float(ss_in['n_inter_D']),
        'phi_bF_F_ss':           float(ss_in['q_b_F']) * float(ss_in['b_F_F']) / (float(ss_in['p']) * float(ss_in['n_inter_F'])),
        'b_F_D_anchor':          float(ss_in['b_F_D']),
        'b_D_F_anchor':          float(ss_in['b_D_F']),
        'gk_wedge_F_D_ss':       w_FD,
        'gk_wedge_D_F_ss':       w_DF,
        # Omega_{t+1} at SS. intermediation_P1_D/F export the dynamic Omega_p1; the SS
        # model builds Omega in steady_auxilliary_D/F instead, and theta is constant at
        # SS, so the two coincide exactly. gk_cross_border_foc needs it in ss_final.
        'Omega_p1_D': float(ss_in['Omega_D']),
        'Omega_p1_F': float(ss_in['Omega_F']),
        'q_b_D':   float(ss_in['q_b_D']),
        'q_b_F':   float(ss_in['q_b_F']),
        'p':       float(ss_in['p']),
        'C_D_ss':  float(ss_in['C_D']),
        'C_F_ss':  float(ss_in['C_F']),
    }
    cal.update(anchors)
    for k, v in anchors.items():
        ss_in.toplevel[k] = v
    ss_in.toplevel['b_F_D_ss']  = float(ss_in['b_F_D'])
    ss_in.toplevel['b_D_F_ss']  = float(ss_in['b_D_F'])
    ss_in.toplevel['Rgross_D']  = float(1 + ss_in['rdep_D'])
    ss_in.toplevel['Rgross_F']  = float(1 + ss_in['rdep_F'])
    _fr_D = float(ss_in['frisch_D']); _fr_F = float(ss_in['frisch_F'])
    ss_in.toplevel['X_D']   = (float(ss_in['C_D'])
                                - float(ss_in['vphi_D']) * float(ss_in['N_D'])**(1+1/_fr_D) / (1+1/_fr_D))
    ss_in.toplevel['X_F']   = (float(ss_in['C_F'])
                                - float(ss_in['vphi_F']) * float(ss_in['N_F'])**(1+1/_fr_F) / (1+1/_fr_F))
    ss_in.toplevel['U_D']   = ss_in.toplevel['X_D'] / float(ss_in['C_D'])
    ss_in.toplevel['U_F']   = ss_in.toplevel['X_F'] / float(ss_in['C_F'])
    ss_in.toplevel['Phi_D'] = float(ss_in['Phi_D'])
    ss_in.toplevel['Phi_F'] = float(ss_in['Phi_F'])
    ss_in.toplevel['value_D'] = (float(ss_in['beta_inter_D'])
                                  * float(ss_in['Omega_D']) * (1 + float(ss_in['rn_D'])))
    ss_in.toplevel['value_F'] = (float(ss_in['beta_inter_F'])
                                  * float(ss_in['Omega_F']) * (1 + float(ss_in['rn_F'])))
    for k, v in {
        'tau_mp_D': 0.0, 'tau_mp_F': 0.0,
        'T_D': 0.0,  'T_F': 0.0,
        'T_ls_D':   0.0, 'T_ls_F':   0.0,
        'b_F_D_res': 0.0, 'b_D_F_res': 0.0,
        'rb_D_res': 0.0,  'rb_F_res': 0.0,
        'labor_mkt_res_D': 0.0, 'labor_mkt_res_F': 0.0,
        'w_res_D': 0.0, 'w_res_F': 0.0,
    }.items():
        ss_in.toplevel[k] = v
    return anchors


def solve_steady_state(calibration_start):
    ha = sj.create_model([
        sdf_ss_D, sdf_banker_ss_D, government_default_D, bond_return_D, gk_bond_foc_D,
        sdf_ss_F, sdf_banker_ss_F, government_default_F, bond_return_F, gk_bond_foc_F,
        # Pledgeability map. Nothing in the SS consumes Delta_*_eff_* (the SS uses
        # steady_auxilliary_D/F, not intermediation_IC_D/F), but the DYNAMIC model does,
        # so ss_final must carry them. Listed here rather than injected in
        # _apply_ss_anchors so there is exactly one copy of the algebra -- a second copy
        # is how audit_artifacts/ drifted (CLAUDE.md). SS-neutral: def_rate_ss = 0.
        collateral_quality_D, collateral_quality_F,
        hh_extended_D, smart_steady_D, market_clearing_D, steady_auxilliary_D,
        banker_div_D, government_ss_D, labor_ss_D, firm_profit_D, price_nkpc_D,
        hh_extended_F, smart_steady_F, market_clearing_F, steady_auxilliary_F,
        banker_div_F, government_ss_F, labor_ss_F, firm_profit_F, price_nkpc_F,
        ces_price_D, import_demand_D, ces_price_F, import_demand_F,
        deposit_rates_D, deposit_rates_F,
        deposit_return_D, deposit_return_F,
        bond_yield,
        trade_balance, external_account_D, global_goods_mkt,
        terms_of_trade, union_inflation,
    ], name='MU HA Model 2 Country')

    # STAGE 2/3 (2026-08-17): q_b_D/q_b_F are now SS UNKNOWNS pinned by the GK
    # portfolio FOC (rb_D_res / rb_F_res from divert_bond_foc_D/F), replacing
    # bond_price_ss_D/F which priced off the banker's SDF alone and took no view on
    # collateral quality -- the reason the old SS violated nu_bD_D/nu_K_D = Delta_bD_D
    # (0.2491 vs 0.20). Cannot be a @simple block: q_b_D -> balance sheet -> K_D ->
    # rk_D -> q_b_D is a cycle in the SS DAG, so it has to go through the solver.
    unknowns_ss = {'beta_D': 0.9850, 'beta_F': 0.9850, 'p': 0.99,
                   'q_b_D': 0.9749, 'q_b_F': 0.9663}
    targets_ss  = ['deposit_mkt_D', 'deposit_mkt_F', 'ca_res_D',
                   'rb_D_res', 'rb_F_res']

    # ── Initial SS solve ──────────────────────────────────────────────────────
    print("Solving initial steady state...")
    ss = ha.solve_steady_state(calibration_start, unknowns_ss, targets_ss, solver='broyden_custom')

    anchors = {
        'phi_bD_D_ss':           float(ss['q_b_D']) * float(ss['b_D_D']) / float(ss['n_inter_D']),
        'phi_bF_F_ss':           float(ss['q_b_F']) * float(ss['b_F_F']) / (float(ss['p']) * float(ss['n_inter_F'])),
        'b_F_D_anchor':          float(ss['b_F_D']),
        'b_D_F_anchor':          float(ss['b_D_F']),
        # No excess_return_*_ss / psi_spread / EL_price anchors -- see _apply_ss_anchors.
        'gk_wedge_F_D_ss':       gk_cross_wedges(ss)[0],
        'gk_wedge_D_F_ss':       gk_cross_wedges(ss)[1],
        'Omega_p1_D':            float(ss['Omega_D']),
        'Omega_p1_F':            float(ss['Omega_F']),
        'q_b_D':                 float(ss['q_b_D']),
        'q_b_F':                 float(ss['q_b_F']),
        'p':                     float(ss['p']),
        'C_D_ss':                float(ss['C_D']),
        'C_F_ss':                float(ss['C_F']),
    }
    calibration_start.update(anchors)
    for k, v in anchors.items():
        ss.toplevel[k] = v

    ss.toplevel['b_F_D_ss'] = float(ss['b_F_D'])
    ss.toplevel['b_D_F_ss'] = float(ss['b_D_F'])
    ss.toplevel['Rgross_D'] = float(1 + ss['rdep_D'])
    ss.toplevel['Rgross_F'] = float(1 + ss['rdep_F'])
    _fr_D = float(ss['frisch_D']); _fr_F = float(ss['frisch_F'])
    ss.toplevel['X_D'] = float(ss['C_D']) - float(ss['vphi_D']) * float(ss['N_D']) ** (1 + 1/_fr_D) / (1 + 1/_fr_D)
    ss.toplevel['X_F'] = float(ss['C_F']) - float(ss['vphi_F']) * float(ss['N_F']) ** (1 + 1/_fr_F) / (1 + 1/_fr_F)
    ss.toplevel['U_D'] = ss.toplevel['X_D'] / float(ss['C_D'])
    ss.toplevel['U_F'] = ss.toplevel['X_F'] / float(ss['C_F'])
    ss.toplevel['Phi_D'] = float(ss['Phi_D'])
    ss.toplevel['Phi_F'] = float(ss['Phi_F'])
    ss.toplevel['value_D'] = float(ss['beta_inter_D']) * float(ss['Omega_D']) * (1 + float(ss['rn_D']))
    ss.toplevel['value_F'] = float(ss['beta_inter_F']) * float(ss['Omega_F']) * (1 + float(ss['rn_F']))
    for k, v in {
        'tau_mp_D': 0.0, 'tau_mp_F': 0.0,
        'T_D': 0.0,  'T_F': 0.0,
        'T_ls_D': 0.0, 'T_ls_F': 0.0,
        'b_F_D_res': 0.0, 'b_D_F_res': 0.0,
        'rb_D_res': 0.0,  'rb_F_res': 0.0,
        'labor_mkt_res_D': 0.0, 'labor_mkt_res_F': 0.0,
        'w_res_D': 0.0, 'w_res_F': 0.0,
    }.items():
        ss.toplevel[k] = v

    # ── Portfolio share targeting ─────────────────────────────────────────────
    # EBA 2011 REBUILD (2026-07-31): bank-sovereign concentration, 31 Dec 2010,
    # read from data/eba_moments.json (single source of truth — do NOT hardcode
    # a second copy here). Regenerate with `python code/eba_calibration.py`.
    #   phi_bD_D = 2.390   GR banks' Greek book / capital  (own, doom loop)
    #   phi_bF_D = 0.018   GR banks' Bund  / capital       (cross)
    #   phi_bD_F = 0.069   DE banks' Greek book / capital  (cross, contagion)
    #   phi_bF_F = 2.758   DE banks' Bund  / capital       (own)
    # Gated by calibration.EBA_CALIBRATION — see the switch note there. At the
    # measured concentration the GK block has no feasible Delta (empty set), so
    # the default is the pre-EBA symmetric placeholder, which solves.
    if calibration.EBA_CALIBRATION:
        print("Targeting portfolio shares (EBA 2011, 31 Dec 2010)...")
        _eba = load_eba_targets()
        target_phi_bD_D = _eba['phi_bD_D_ss']
        target_phi_bF_D = _eba['phi_bF_D_ss']
        target_phi_bD_F = _eba['phi_bD_F_ss']
        target_phi_bF_F = _eba['phi_bF_F_ss']
    else:
        print("Targeting portfolio shares (pre-EBA symmetric placeholder)...")
        target_phi_bD_D = 0.25
        target_phi_bF_D = 0.15
        target_phi_bD_F = 0.15
        target_phi_bF_F = 0.25

    n_D  = float(ss['n_inter_D'])
    n_F  = float(ss['n_inter_F']) * float(ss['p'])
    q_D  = float(ss['q_b_D'])
    q_F  = float(ss['q_b_F'])

    # COUNTRY SIZE (2026-08-07). Each phi is a ratio to its HOLDER's net worth, so
    # every stock below is in its holder's own per-capita units -- which is exactly
    # the convention the model now uses. The aggregation happens in
    # domestic_bond_clearing via size_F, not here.
    size_F = float(calibration_start['size_F'])

    b_D_D_new = target_phi_bD_D * n_D / q_D     # D aggregate
    b_F_D_new = target_phi_bF_D * n_D / q_F     # D aggregate
    b_D_F_new = target_phi_bD_F * n_F / q_D     # per F capita
    b_F_F_new = target_phi_bF_F * n_F / q_F     # per F capita

    # Government stocks, each in ITS OWN country's units:
    #   D debt is a D aggregate;  F debt is per F capita.
    B_D_new   = b_D_D_new + size_F * b_D_F_new
    B_F_new   = b_F_F_new + b_F_D_new / size_F

    # Over-identifying checks. The two EBA moments that the pre-size-asymmetry
    # calibration could not satisfy jointly -- portfolio composition (the phi's,
    # matched by construction above) and market structure (the foreign shares
    # below). Both should now hold. B_supply_*_qgdp and the foreign shares are
    # measured directly and were never used by the targeting, which is why the
    # inconsistency went unnoticed for so long.
    _eba_all = load_eba_targets()
    _fs_eba  = load_eba_foreign_shares()
    fs_D_new = size_F * b_D_F_new / B_D_new
    fs_F_new = (b_F_D_new / size_F) / B_F_new
    print(f"  size_F = {size_F:.4f}  (F/D GDP; every F variable is per F capita)")
    print(f"  foreign-held share of the bank-held sovereign stock:"
          f"  D = {fs_D_new:.4f} (EBA {_fs_eba['D']:.4f})"
          f"   F = {fs_F_new:.6f} (EBA {_fs_eba['F']:.6f})")
    print(f"  bank-held stock in own-country quarterly GDP:"
          f"  B_D = {B_D_new:.4f} (EBA {_eba_all['B_supply_D_qgdp']:.4f})"
          f"   B_F = {B_F_new:.4f} (EBA {_eba_all['B_supply_F_qgdp']:.4f})")

    # omega_K is MEASURED (corporate+CRE EAD / K), not back-solved. It stays at
    # its calibration value; K is then an OUTPUT of the balance sheet,
    #   K = (theta·N - q_b·bonds) / (omega_K·Q),
    # and comparing it to the conventional K/Y_annual=2.7 target (K=10.8 in
    # quarterly-GDP units) is a genuine over-identifying check on the EBA
    # balance sheet. The 2026-07-22 build had this backwards: it ASSUMED
    # theta=4.0 and solved omega_K to force K=10.8, so the check was vacuous and
    # omega_K was a free plug absorbing the theta assumption.
    omega_K_D_new = float(calibration_start['omega_K_D'])
    omega_K_F_new = float(calibration_start['omega_K_F'])

    _fr = float(calibration_start['fund_rule_D'])
    # Built from the REALISED stocks, not `target_phi * n`. Post-units-fix the two
    # differ for the cross-border legs, and using the targets would compute this
    # over-identifying check against the pre-fix bond book — overstating K_implied_F
    # by ~0.13 and quietly making a broken balance sheet look like it validated.
    # Each bank's own balance sheet, in its own per-capita units -- no size_F here:
    # the F bank holds b_D_F per F capita and funds it per F capita.
    _bank_D = (float(calibration_start['theta_D']) * n_D
               - (q_D * b_D_D_new + q_F * b_F_D_new))
    _bank_F = (float(calibration_start['theta_F']) * n_F
               - (q_F * b_F_F_new + q_D * b_D_F_new))
    K_implied_D = ((1 - _fr) * _bank_D / omega_K_D_new
                   + _fr * (_bank_D + float(calibration_start['K_fund_D'])))
    K_implied_F = ((1 - _fr) * _bank_F / omega_K_F_new
                   + _fr * (_bank_F + float(calibration_start['K_fund_F'])))

    print(f"  D-bank: phi_bD_D = {target_phi_bD_D:.3f}  phi_bF_D = {target_phi_bF_D:.3f}"
          f"  omega_K_D = {omega_K_D_new:.4f}")
    print(f"  F-bank: phi_bD_F = {target_phi_bD_F:.3f}  phi_bF_F = {target_phi_bF_F:.3f}"
          f"  omega_K_F = {omega_K_F_new:.4f}")
    print(f"  over-identifying check (K from the measured balance sheet, target 10.8):"
          f"  K_D = {K_implied_D:.3f}   K_F = {K_implied_F:.3f}")

    calibration_start.update({
        'b_D_D': b_D_D_new, 'b_F_D': b_F_D_new,
        'b_D_F': b_D_F_new, 'b_F_F': b_F_F_new,
        'b_F_D_anchor': b_F_D_new, 'b_D_F_anchor': b_D_F_new,
        'phi_bF_D_ss': target_phi_bF_D,
        'omega_K_D': omega_K_D_new, 'omega_K_F': omega_K_F_new,
        'B_supply_D': B_D_new, 'b_gov_D': B_D_new, 'b_gov_ss_D': B_D_new,
        'B_supply_F': B_F_new, 'b_gov_F': B_F_new, 'b_gov_ss_F': B_F_new,
    })

    print("Re-solving SS with new portfolio allocation...")
    # Must carry EVERY unknown in unknowns_ss, warm-started at the first solve's
    # values. q_b_D/q_b_F joined the unknown set in the stage-2/3 refactor; omitting
    # them here left 3 unknowns against 5 targets and broyden_solver ran to its
    # 100-iteration cap without converging.
    _unknowns_warm = {'beta_D': float(ss['beta_D']), 'beta_F': float(ss['beta_F']),
                      'p':      float(ss['p']),
                      'q_b_D':  float(ss['q_b_D']), 'q_b_F': float(ss['q_b_F'])}
    ss = ha.solve_steady_state(calibration_start, _unknowns_warm, targets_ss, solver='broyden_custom')
    _apply_ss_anchors(ss, calibration_start)
    print(f"SS re-solved. beta_D={float(ss['beta_D']):.8f}  p={float(ss['p']):.6f}")
    report_gk_steady_state(ss, calibration_start)

    return {
        'ss':                ss,
        'ha':                ha,
        'calibration_start': calibration_start,
        'unknowns_ss':       unknowns_ss,
        'targets_ss':        targets_ss,
    }
