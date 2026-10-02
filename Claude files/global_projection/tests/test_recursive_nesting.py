# RECURSIVE NESTING GATES N1-N4.
# N1  the SS is a rest point of the single-point map: with SS-constant rules in
#     both regimes, no default, and the state at the SS point, every residual
#     vanishes (the recursive image of "the zero-shock transition stays at the SS").
#     Requires the rep-agent household anchors.
# N2  the no-default (pi=0) d=0 block is ROOTED grid-wide by the collocation Newton.
# N3  a TPI purchase swaps bonds for the safe claim inside the period (assets, the
#     divertable base, deposits, P', net worth unchanged) and transfers the default
#     risk to the Eurosystem next period, by exactly (qR - Xi')*dm in every regime.
# N4  every budget in the period map sums to the union goods market, TPI on or off,
#     at arbitrary (not solved) points: the TPI's flows cancel exactly.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

import numpy as np

from common import get_ss
from global_projection.solver_recursive.state_grid import (build_state_box, s_process_params,
                                                           IM, IO)
from global_projection.blocks.trade import import_demand, trade_balance, size_ratio
from global_projection.solver_recursive.decision_rules import RuleSet
from global_projection.solver_recursive.point_map import point_residuals
from global_projection.solver_recursive.recursive_main import (time_iteration, calibrate_household_anchors,
                            ss_state, ss_x)
from global_projection.solver_recursive.recursive_experiment import BOX_KW, advance
from global_projection.solver_recursive.collocation import (solve_collocation, TOL_MAXF, RES_NAMES,
                                          N_RES_POINT)

# SINGLE-SOURCED FROM THE SOLVER. This list used to be a hand-kept copy and went stale
# at every change to the residual system -- the assert below is what catches that now.
_LABELS = RES_NAMES[:N_RES_POINT]


def test_n1_ss_rest_point():
    # SS + SS-CONSTANT RULES + NO DEFAULT => EVERY RESIDUAL ~ 0.
    # The label list must cover the whole residual vector: when the system grew from 7
    # to 11 this printed only the first seven, so the bond FOCs, the F household's
    # Euler and union deposit clearing were silently unreported (the assert on
    # max|res| still covered them, but nothing showed WHICH one moved).
    cal, ss = get_ss()
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    grid = build_state_box(ss, cal)
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=2)

    res, out = point_residuals(ss_state(ss, cal, sp), 0, ss_x(ss, cal),
                               rules, cal, ss, sp, no_default=True)
    worst = np.max(np.abs(res))
    assert len(_LABELS) == res.size, (
        f"residual vector is {res.size} long, _LABELS covers {len(_LABELS)}")
    for lab, r in zip(_LABELS, res):
        print(f"    {lab:10s} {r:+.3e}")
    print(f"    mu_D={out['mu_D']:.6f} (ss {ss['ss_bank_D']['mu_ss']:.6f})  "
          f"C_D={out['C_D']:.5f} (ss {ss['C_D_ss']:.5f})  "
          f"A_D={out['A_D']:.5f} (ss {ss['A_D_ss']:.5f})")
    assert worst < 1e-6, f"SS not a rest point: max|res|={worst:.2e}"


def _tpi_setup():
    # TPI ON, NET-WORTH FLOOR OFF (its smooth bias would blur exact accounting), SS RULES.
    cal, ss = get_ss()
    cal = dict(cal, tpi_on=True, nw_floor_frac=0.0)
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    grid = build_state_box(ss, cal, mu=1, **BOX_KW)
    return cal, ss, sp, grid, RuleSet.from_ss(grid, ss, cal)


def test_n3_tpi_swap_and_risk_transfer():
    # N3: A PURCHASE IS A SWAP TODAY AND A RISK TRANSFER TOMORROW.
    # Today the bank sells dm bonds for Z = Q*dm of the safe claim, so its assets, its
    # divertable base (single lambda), its deposits, P' and its net worth are unchanged.
    # Tomorrow it is owed R*Q*dm instead of the bonds' payoff Xi'*dm, in each regime d':
    #   dn' = [(1-f) + omega_ent*kappa_D]*(R*Q - Xi')*dm,   dPi' = -(R*Q - Xi')*dm.
    # (the omega*kappa term: D's share of the Eurosystem P&L is new D debt the D bank
    # absorbs, and the entrant transfer is a share of assets.)
    # the purchase is FORCED (rule off: m = B*x exactly, Q free) -- the accounting of a
    # purchase does not depend on what decided it, and N4 covers the rule's own mapping
    cal, ss, sp, _, rules = _tpi_setup()
    cal["tpi_on"] = False
    S = ss_state(ss, cal, sp)
    S[IM] = 0.10 * cal["B_gov_D_ss"]
    S[IO] = 0.95 * ss["Q_bD_ss"] * S[IM]          # a legacy book bought below today's price
    x = ss_x(ss, cal); x[-1] = 0.02

    def read(S, d, x):
        _, o = point_residuals(S, d, x, rules, cal, ss, sp, n_gh=5)
        o["_x"] = x
        return o

    xb = x.copy(); xb[-1] += 0.01
    a, b = read(S, 0, x), read(S, 0, xb)
    dm, q = b["m_cb"] - a["m_cb"], a["Q_bD"]
    assert abs(dm - 0.01 * cal["B_gov_D_ss"]) < 1e-15 and b["Q_bD"] == q
    for k in ("dep_D", "n_D", "Pp_D", "slack_D", "Vp_dep", "mu_D"):
        assert abs(b[k] - a[k]) < 1e-12, f"the swap moved {k}: {b[k] - a[k]:.3e}"
    assert abs((b["b_D_D_new"] - a["b_D_D_new"]) + dm) < 1e-12
    assert abs((b["Z_cb"] - a["Z_cb"]) - q * dm) < 1e-12
    Sa, Sb = advance(a, S, sp), advance(b, S, sp)
    xn, db = ss_x(ss, cal), cal["delta_b_D"]
    R = 1.0 + x[4]
    for dn in (0, 1):
        oa, ob = read(Sa, dn, xn), read(Sb, dn, xn)
        Xi = (cal["recovery_rate_D"] if dn else 1.0) * (db + (1.0 - db) * oa["Q_bD"])
        gain = (R * q - Xi) * dm
        want = ((1.0 - cal["f_D"]) + cal["omega_ent_D"] * cal["tpi_key_D"]) * gain
        assert abs((ob["Pi_cb"] - oa["Pi_cb"]) + gain) < 1e-12, f"P&L mirror off in d'={dn}"
        assert abs((ob["n_D"] - oa["n_D"]) - want) < 1e-12, (
            f"risk transfer off in d'={dn}: {ob['n_D'] - oa['n_D']:.6e} vs {want:.6e}")
        print(f"    N3 d'={dn}: dn' = {ob['n_D'] - oa['n_D']:+.4e} = "
              f"[(1-f)+omega*kappa]*(RQ - Xi')*dm, P&L mirror exact")


def _union_identity(o, x, cal, ss):
    # UNION GOODS GAP MINUS THE KNOWN TERMS (F treasury, working capital vs the anchor,
    # deposit clearing), with consumption UNGUARDED -- zero up to the bond floors' bias.
    sz, p = size_ratio(cal), o["p"]
    PD, PF = o["P_CES_D"], o["P_CES_F"]
    CD = o["W_D"] / PD + o["inc_D"] - o["A_D"]
    CF = o["W_F"] / PF + o["inc_F"] - o["A_F"]
    a = np.array
    IMD = import_demand(a([p]), a([CD]), a([PD]), cal, "D")[0]
    IMF = import_demand(a([p]), a([CF]), a([PF]), cal, "F")[0]
    NXD, NXF = (z[0] for z in trade_balance(a([p]), a([IMD]), a([IMF]), cal))
    gap = ((o["Y_D"] - PD * CD - o["I_D"] - NXD - cal["G_D"])
           + sz * p * (o["Y_F"] - PF * CF - o["I_F"] - NXF - cal["G_F"]))
    treas = sz * p * cal["delta_b_F"] * cal["B_gov_F_ss"] * (o["Q_bF"] - ss["Q_bF_ss"])
    wcD = cal["zeta_wc_D"] * o["r_wc_D"] * o["w_D"] * x[0] + o["L_wc_D"] - PD * ss["hh_T_D"]
    wcF = cal["zeta_wc_F"] * o["r_wc_F"] * o["w_F"] * x[1] + o["L_wc_F"] - PF * ss["hh_T_F"]
    return gap - treas - wcD - sz * p * wcF - (o["save_union"] - o["dep_union"])


def test_n4_union_budget_identity():
    # N4: THE TPI's FLOWS CANCEL OUT OF THE UNION BUDGET, AT ANY POINT, ON OR OFF.
    # Not a solved point: random states (live M, O), random policies (m > 0), both
    # regimes. The identity is algebra, so it must hold wherever the period map is
    # evaluated; the remainder is the smooth bond floors' bias (~1e-10). An error in any
    # TPI flow -- the capital key, the sz*p conversion, O vs Z -- is O(1e-3) or larger.
    cal, ss, sp, grid, rules = _tpi_setup()
    rng = np.random.default_rng(0)
    for on, gz in ((False, True), (True, True), (True, False)):
        c = dict(cal, tpi_on=on, tpi_gz=gz)
        worst = 0.0
        for _ in range(20):
            S = grid.lo + (grid.hi - grid.lo) * rng.random(grid.d)
            x = ss_x(ss, c) * (1.0 + 0.01 * rng.standard_normal(len(ss_x(ss, c))))
            x[-1] = 0.05 * rng.random()
            for d in (0, 1):
                _, o = point_residuals(S, d, x, rules, c, ss, sp, n_gh=5)
                worst = max(worst, abs(_union_identity(o, x, c, ss)))
        assert worst < 1e-8, f"union budget leaks with tpi_on={on}, gz={gz}: {worst:.3e}"
        print(f"    N4 tpi_on={on} ({'GZ' if gz else 'FB'} form): union identity holds to "
              f"{worst:.1e} at 40 random points")


def test_n2_no_default_grid_solve():
    # N2: THE GRID-WIDE FIXED POINT AT pi = 0 MUST ACTUALLY BE FOUND.
    # This used to be a REPORTING probe rather than a gate, because damped time
    # iteration could not converge it: its binding mode is the franchise-value
    # recursion at 0.990 per sweep, so it reported "converged=False, worst point
    # residual 2e-14" -- every point clearing against a continuation still moving.
    # The global collocation Newton (solver_recursive/collocation.py) roots the whole
    # system instead, so this is a hard assert now.
    cal, ss = get_ss()
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    # bands MUST match what solve_recursive ships, or the probe measures a box defect
    # rather than the solver -- imported from BOX_KW so the two cannot drift apart.
    # (Hard-wired bands here went stale at the 9-state change: the old 0.12/0.20 P
    # bands admit corners the period map cannot solve now that the household claim W_D
    # is a separate state, and the probe reported their |F| = 9.0 as a solver failure.)
    grid = build_state_box(ss, cal, mu=1, **BOX_KW)
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=2)

    rules.n_gh = 5
    time_iteration(rules, cal, ss, sp, regimes=(0,), no_default=True, damp=0.5,
                   tol=1e-4, max_it=12, n_gh=5, verbose=False)                               # warm start only
    ok, its, worst = solve_collocation(rules, cal, ss, sp, regimes=(0,),
                                       no_default=True, n_gh=5, backend="parsolve",
                                       maxit=20, verbose=False, label=" N2")
    assert ok, f"N2: collocation solve did not converge (max|F| = {worst:.2e})"
    assert worst <= 10 * TOL_MAXF, f"N2: max|F| = {worst:.2e}"
    print(f"    N2 (grid-wide solve at pi=0): converged in {its} Newton steps, "
          f"max|F| = {worst:.2e}")


if __name__ == "__main__":
    test_n1_ss_rest_point()
    print("test_recursive_nesting N1 (SS rest point): PASSED")
    test_n3_tpi_swap_and_risk_transfer()
    print("test_recursive_nesting N3 (TPI swap and risk transfer): PASSED")
    test_n4_union_budget_identity()
    print("test_recursive_nesting N4 (union budget identity): PASSED")
    test_n2_no_default_grid_solve()
    print("test_recursive_nesting N2 (grid-wide solve): PASSED")
