# RECURSIVE SOLUTION EXPERIMENT: THE PASS-THROUGH OF SOVEREIGN RISK (BOCOLA).
# End-to-end demonstration that the global recursive solution (two-branch
# quadrature + Bocola closed-form multiplier) delivers the sign the representative
# branch could not: an elevated PRICED probability of default lowers output AND
# consumption, persistently. Orchestration only -- the economics live in the
# untouched blocks and in point_map.py; the solver is recursive_main.time_iteration.
#
# Pipeline: solve the d=0 (no-default) rules, warm-start d=1 (post-default) from
# them and solve, refine jointly at pi>0, then read (i) the IMPACT response as the
# one-quarter-ahead default probability rises and (ii) the persistence IRF as an
# s-shock decays (rho_s). Spawn-guarded (the pointwise solves are serial, but the
# guard is kept so the module is safe to import).
import numpy as np
from scipy.optimize import root

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import build_state_box, s_process_params, default_prob
from solver_recursive.state_grid import (IK_D, IK_F, IP_D, IP_F, IBDD,
                                          IBDF, IBFD, IV, IS, IZ, STATE_NAMES)
from solver_recursive.decision_rules import RuleSet, STORE_RULES, regime_table
from solver_recursive.collocation import solve_collocation
from solver_recursive.recursive_main import (time_iteration, calibrate_household_anchors,
                            ss_state, ss_x, p_block_rotation)
from reporting.prints import (print_solve_stage, bp_ann, ann_pct, ann_prob,
                              BOCOLA_IRF_CLOSED, BOCOLA_IRF_OPEN,
                              BOCOLA_EPISODE_LEVEL)
from solver_recursive.point_map import point_residuals, SOLVE7


# STATE NAMES, for the box-escape report in dynamic_irf -- taken from state_grid so a
# state added there cannot silently relabel the report (this list had been left at nine
# names, and at "W_D", through two state-vector changes).
_SNAMES = STATE_NAMES

# THE COLLOCATION BOX, SHARED BY EVERY EXPERIMENT. It used to be passed only to the
# sovereign-risk solve, leaving solve_tfp on build_state_box's bare defaults, whose
# +-25% P band has no solution in the period map. One box for both experiments, so a
# band that is feasible for one is feasible for the other.
BOX_KW = dict(k_band=0.02, p_band_D=0.04, p_band_F=0.04, b_band=0.12, w_band=0.04)

# QUADRATURE ORDER FOR THE SOLVE. Every reader takes it off rules.n_gh, so the
# reported IRFs are evaluated under the same measure the rules were solved under.
N_GH = 5

# DENSE CHEBYSHEV NODES TENSORED ONTO THE s DIMENSION (state_grid.SmolyakGrid refine).
# All the curvature in this model is the logistic p^d(s); the other nine states are
# near-linear, and raising the Smolyak level to reach s would pay for resolution in all
# of them. Measured relative RMS error on this model's curvature profile:
#   isotropic mu=1   21 pts  1.9e-1        s_refine=5   95 pts  2.5e-2
#   isotropic mu=2  221 pts  3.9e-2        s_refine=9  171 pts  1.1e-3
# 5 IS THE SHIPPED DEFAULT and 9 is Bocola's own resolution (his mu = 3 grid carries
# m(mu+1) = 9 nodes per dimension). The difference is cost, not correctness: the dense
# Jacobian is m+1 residual evaluations and m = 19*2*n, so the solve scales as n^2 --
# ~17 min per Jacobian at 95 points against ~55 min at 171, i.e. ~70 min against ~4 h
# for the refinement stage. The ladder in solve_recursive walks 5 -> 9 automatically
# when S_REFINE = 9, seeding each rung from the last. Set 0 or 1 for the plain grid.
S_REFINE = 5

# WARM-START SWEEPS BEFORE EACH NEWTON. Time iteration is globally stable but converges
# at 0.990 per sweep on the franchise-value mode; it is used here ONLY to get inside the
# Newton's basin, never to converge. Bocola does the same thing with a warm start from a
# previously solved model (model_solution_mean.m seeds the 6-state solve from the solved
# 5-state no-default one).
WARM_SWEEPS = 12
# The refined stages need a LONGER warm start than the coarse one -- see solve_recursive.
REFINE_WARM_SWEEPS = 20


# THE FACILITY HOMOTOPY. m = 0 is IDENTICALLY the no-facility model (the nesting gate
# proves it), so the walk exists only because a LARGE envelope drives mu to zero over much
# of the grid, which is a big move to ask of one Newton step from a seed where mu > 0
# everywhere. At the shipped envelope (~1-1.5% of quarterly GDP, sized in
# docs/ltro_backstop_plan.md S5 to HALVE the crisis multiplier rather than eliminate it)
# a single rung is normally enough; the ladder is insurance for the oversized variants.
# One rung at the shipped envelope: it is small (2% of quarterly GDP) and each facility
# regime is seeded from its OWN no-facility twin, so the Newton starts close. Add rungs
# for the oversized variants, where mu is driven to zero over much of the grid and that
# is a large move to ask of one step. A failed rung still leaves the final joint polish.
LTRO_LADDER = (0.34, 0.67, 1.0)


def _seed_from(rules_fine, rules_coarse):
    # EVALUATE A SOLVED COARSE RULE SET AT THE FINE GRID'S POINTS.
    # Both grids are drawn on the SAME box, so this is interpolation, not extrapolation.
    pts = rules_fine.grid.points
    for k in STORE_RULES:
        for d in range(rules_coarse.n_regimes):
            rules_fine.set_values(k, d, rules_coarse.eval(k, d, pts))
    rules_fine.n_gh = rules_coarse.n_gh
    return rules_fine


def _stage(rules, cal, ss, sproc, regimes, no_default, label, verbose,
           backend="auto", warm=WARM_SWEEPS, maxit=40, no_cb=False):
    # ONE SOLVE STAGE: a short time-iteration warm start, then the GLOBAL NEWTON.
    if warm:
        time_iteration(rules, cal, ss, sproc, regimes=regimes, no_default=no_default,
                       damp=0.5, tol=1e-4, max_it=warm, n_gh=N_GH, verbose=False,
                       no_cb=no_cb)
    return solve_collocation(rules, cal, ss, sproc, regimes=regimes,
                             no_default=no_default, n_gh=N_GH, backend=backend,
                             maxit=maxit, verbose=verbose, no_cb=no_cb,
                             label=f" {label}")


def liquidity_ceiling_report(rules, cal, ss, sproc, target=None, regime=0,
                             verbose=True):
    # HOW MUCH OF THE BOND PRICE CAN A CONSTRAINT-RELIEF INSTRUMENT DELIVER? THE CEILING.
    # The D bank's own FOC is  E[Om*payD] = Q*(E[Om]*R + lambda_bD*mu). ANY policy that
    # works by relaxing the incentive constraint -- a bond purchase shrinking the
    # divertable base, or an LTRO doing that AND adding to the constraint's numerator --
    # raises Q only by pushing mu down, and mu is floored at zero. So, HOLDING THE
    # CONTINUATION FIXED, no such instrument at any size can lift the price above
    #     Q_max = E[Om*payD] / (E[Om] * R),
    # the same claim priced with the LIQUIDITY premium removed and NOTHING else. The
    # expected loss and the risk premium are untouched by any amount of constraint relief
    # at a given continuation.
    # This is the ceiling on channel (a) in docs/ltro_backstop_plan.md S3, and it is why
    # the interesting question is channel (b): the ANNOUNCEMENT raises E[Om*payD] itself,
    # which moves the ceiling rather than approaching it. Measured at the 100 bp
    # calibration on the coarse grid, the liquidity premium is 0.2-0.7% of the price over
    # the ergodic states and 0.63% at the crisis corner -- which is why a yield peg set at
    # the rest-point price (a 22.2% gap at that corner) is not deliverable by quantity,
    # and why the instrument had to change.
    rows = []
    for i in np.argsort(rules.grid.points[:, IS]):
        S = rules.grid.points[i]
        x = np.array([rules.vals[k][regime][i] for k in SOLVE7])
        _, o = point_residuals(S, regime, x, rules, cal, ss, sproc,
                               n_gh=rules.n_gh or N_GH)
        q_max = o["E_Om_payD"] / (o["E_Om_D"] * (1.0 + o["rdep_D"]))
        rows.append((float(default_prob(S[IS])), o["Q_bD"], q_max,
                     q_max / o["Q_bD"] - 1.0,
                     np.nan if target is None else target / o["Q_bD"] - 1.0))
    if verbose:
        head = ("    p^d %/q    Q_free   Q_max(mu=0)   liquidity"
                + ("" if target is None else "   gap to target   reachable?"))
        print("\n  LIQUIDITY CEILING: the most any constraint-relief instrument can "
              "deliver at a\n  FIXED continuation (mu -> 0). Anything beyond it must "
              "come from the announcement.")
        print(head)
        for pd_, q_f, q_m, liq, gap in rows:
            line = f"    {100*pd_:7.3f} {q_f:9.5f} {q_m:13.5f} {100*liq:10.2f}%"
            if target is not None:
                line += f" {100*gap:14.2f}% {'YES' if q_m >= target else 'no':>12s}"
            print(line)
        print(f"    -> median liquidity premium "
              f"{100*np.median([r[3] for r in rows]):.2f}% of the price")
    return rows


def solve_recursive(cal, ss, sproc, mu=1, verbose=True, mu_vec=None, rotate=False,
                    s_refine=S_REFINE, backend="auto", with_cb=False,
                    base=None, base_out=None):
    # SOLVE EVERY REGIME BY GLOBAL COLLOCATION (Bocola's model_solution_mean.m).
    #
    # THE LADDER IS HIS. He warm-starts the 6-state default model from the solved
    # 5-state no-default one, walks the haircut up in seven steps re-solving at each
    # (rec = 0, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55), and calls parsolve -- a damped
    # Newton on the WHOLE coefficient vector -- at every step. Here:
    #   1. coarse grid (isotropic mu = 1), d = 0 at pi = 0 with the facility off, Newton;
    #   2. the default regime by haircut homotopy 0.85 -> 0.70 -> 0.55 -> recovery_rate_D;
    #   3. the risk-priced baseline, still with the facility off -- the model the backstop
    #      is measured against;
    #   4. the facility regimes, each seeded from its OWN no-facility twin, walked up in
    #      size (LTRO_LADDER). m = 0 is IDENTICALLY the no-facility model, so the first
    #      rung is free; the walk only exists because a large facility drives mu to zero,
    #      which is a big move for one Newton step;
    #   5. joint polish over every regime;
    #   6. if s_refine, rebuild on the s-refined grid, SEED from the coarse solution,
    #      and re-run the joint Newton there.
    # Every stage is a genuine root of the collocation system, not a damped fixed point.
    #
    # phi is a per-experiment SCALAR (cal["phi_ltro"]), not a state and not a stage:
    # each activation intensity is its own solve, exact at its own phi.
    #
    # with_cb DEFAULTS TO FALSE, AND THAT DEFAULT IS LOAD-BEARING. The facility regimes
    # DOUBLE the regime count, and the dense Jacobian is (unknowns+1) residual
    # evaluations each costing points x regimes, so switching them on costs 4x -- at
    # s_refine = 5 the risk solve goes from 3610 unknowns and ~29 min a Jacobian to 7220
    # and ~114 min. Because phi_ltro = 0 nests the no-backstop model EXACTLY, a caller
    # that gets the facility regimes by accident pays that 4x and gets an identical
    # answer, with nothing in the output to say so. That is what happened when the
    # default was True: main.py's headline risk experiment ran four regimes for hours.
    # Only the backstop experiment should opt in.
    #
    # rotate=True collocates the P block on the eigenbasis of its own transition
    # Jacobian (Bocola's V-transform). The theory says it should win -- rho(|J|) = 1.96,
    # so no axis-aligned box is one-step invariant -- but it MEASURES WORSE, because the
    # rotated box is a parallelogram in natural coordinates that reaches P_F states the
    # model cannot solve. Kept, tested and off by default.
    rot, centre = (None, None)
    box_kw = dict(BOX_KW)
    if rotate:
        rot, centre, J, evs = p_block_rotation(ss, cal, sproc, mu=mu, mu_vec=mu_vec,
                                               probe_kw=box_kw)
        if verbose:
            print(f"  P-block rotation: |lambda| = {evs[0]:.3f}, {evs[1]:.3f}"
                  f"   rho(|J|) = {np.max(np.abs(np.linalg.eigvals(np.abs(J)))):.3f}")
    nreg = 4 if with_cb else 2
    grid = build_state_box(ss, cal, mu=mu, mu_vec=mu_vec, rot=rot, centre=centre,
                           **box_kw)
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=nreg)
    rules.n_gh = N_GH
    # REGIME INDICES BY MEANING, NEVER BY POSITION. With the facility available in the
    # default state the table is (0,0),(0,1),(1,0),(1,1), so the pure-default regime is
    # index 2 and NOT the last one -- reading it as nreg-1 would silently run the haircut
    # homotopy on the facility regime instead.
    reg = regime_table(nreg)
    D_REG = next(j for j, (d, m) in enumerate(reg) if d and not m)
    CB_REGS = [j for j, (d, m) in enumerate(reg) if m]
    if verbose:
        print(f"  coarse grid: mu={mu}, {grid.n} points x {nreg} regimes, n_gh={N_GH}")

    # THE BASELINE STAGES DO NOT DEPEND ON phi. An activation sweep re-solves them once
    # per point unless the caller hands them back in, which at eight activations is about
    # an hour of identical arithmetic -- hence the `base` argument.
    if base is None:
        ok0, it0, w0 = _stage(rules, cal, ss, sproc, (0,), True, "d0", verbose,
                              no_cb=True)

        for k in STORE_RULES:                   # warm-start the default regime from d=0
            rules.set_values(k, D_REG, rules.vals[k][0].copy())
        rec_target = cal["recovery_rate_D"]
        ok1, w1 = False, np.nan
        for rec in (0.85, 0.70, 0.55, rec_target):
            cal["recovery_rate_D"] = rec
            ok1, it1, w1 = _stage(rules, cal, ss, sproc, (D_REG,), False,
                                  f"d1 rec={rec:.2f}", verbose, no_cb=True)
        cal["recovery_rate_D"] = rec_target

        # THE RISK-PRICED BASELINE, still with the facility switched off. This is the
        # model the backstop is measured against, and every facility regime is seeded
        # from it.
        okb, itb, wb = _stage(rules, cal, ss, sproc, (0, D_REG), False,
                              "joint (no facility)", verbose, backend=backend,
                              no_cb=True)
        if base_out is not None:                # hand the caller a reusable snapshot
            base_out.append(rules.copy())
    else:
        assert base.n_regimes == nreg and base.grid.n == grid.n, \
            "reused baseline must carry the same grid and regime count"
        rules = base.copy()
        ok0 = ok1 = okb = True
        w0 = w1 = wb = np.nan
        if verbose:
            print("  reusing the solved phi-independent baseline")

    if with_cb:
        # EACH FACILITY REGIME IS SEEDED FROM ITS OWN NO-FACILITY TWIN -- (0,1) from
        # (0,0) and (1,1) from (1,0) -- so the seed already carries the right default
        # state and only the constraint has to move.
        for j in CB_REGS:
            twin = next(k for k, (d, m) in enumerate(reg)
                        if d == reg[j][0] and not m)
            for k in STORE_RULES:
                rules.set_values(k, j, rules.vals[k][twin].copy())
        m0 = (cal["ltro_D"], cal["ltro_F"])
        for frac in LTRO_LADDER:
            cal["ltro_D"], cal["ltro_F"] = frac * m0[0], frac * m0[1]
            _stage(rules, cal, ss, sproc, tuple(range(nreg)), False,
                   f"ltro {100 * frac:.0f}% of envelope", verbose, backend=backend,
                   warm=4, maxit=12)
        cal["ltro_D"], cal["ltro_F"] = m0

    okj, itj, wj = _stage(rules, cal, ss, sproc, tuple(range(nreg)), False, "joint",
                          verbose, backend=backend, no_cb=not with_cb)

    if s_refine and s_refine > 1:
        # THE REFINEMENT LADDER. Same box, a dense Chebyshev factor in s, seeded from
        # the previous solution each time -- one more continuation step in Bocola's
        # style. Going straight from 3 nodes in s to 9 asks the Newton to start from a
        # seed that is badly wrong at the new interior nodes (the mu=1 quadratic reads
        # p^d = 1.82% where the truth is 0.67%), so the node count is walked up.
        ladder = [m for m in (5, 9, 17) if 1 < m < s_refine] + [s_refine]
        for m_s in ladder:
            gfine = build_state_box(ss, cal, mu=mu, mu_vec=mu_vec, rot=rot,
                                    centre=centre, refine=(IS, m_s), **box_kw)
            if verbose:
                print(f"  s-refined grid: {gfine.n} points x {nreg} regimes "
                      f"({m_s} nodes, degree {m_s - 1} in s)")
            fine = _seed_from(RuleSet(gfine, nreg), rules)
            # THE WARM START IS NOT OPTIONAL HERE. The coarse mu = 1 grid puts only ONE
            # coordinate off centre per point; the refined grid is its product with the
            # dense s factor, so it visits (P_D at +1, s at +1)-type combinations the
            # coarse interpolant has no cross term for. Measured: the raw seed sits at
            # max|F| = 4.0e-2 with 62% of points above 1e-3, and 20 time-iteration
            # sweeps (132 s at 95 points, every point clearing at 3e-14) bring it to
            # 1.9e-3 -- inside the Newton's basin.
            okj, itj, wj = _stage(fine, cal, ss, sproc, tuple(range(nreg)), False,
                                  f"joint (s={m_s})", verbose, backend=backend,
                                  warm=REFINE_WARM_SWEEPS, maxit=12,
                                  no_cb=not with_cb)
            rules = fine

    # STAMP THE VERDICT ON THE RULE SET. A caller that plots or tabulates several solves
    # has no other way to know which of them actually rooted -- the NOTE below goes to the
    # console and the figure does not see it. Without this the certainty curve would draw
    # a stopped solve with the same solid marker as a converged one.
    rules.solve_ok = bool(ok0 and ok1 and okb and okj)
    rules.solve_worst = float(np.nanmax([w0, w1, wb, wj]))
    if not rules.solve_ok:
        print(f"    NOTE: a collocation stage did not reach the acceptance floor: "
              f"d0 ok={ok0} ({w0:.1e}), d1 ok={ok1} ({w1:.1e}), "
              f"base ok={okb} ({wb:.1e}), joint ok={okj} ({wj:.1e}). "
              f"Read the IRFs as indicative.")
    return rules


def _solve_point(rules, cal, ss, sproc, S, x0):
    # DIRECT d=0 IMPACT SOLVE AT STATE S (both regimes as continuation).
    ngh = rules.n_gh or N_GH

    def f(x):
        try:
            return point_residuals(S, 0, x, rules, cal, ss, sproc,
                                   n_gh=ngh, no_default=False)[0]
        except (ValueError, RuntimeError, FloatingPointError):
            return np.full(len(SOLVE7), 10.0)
    best = None
    for g in (x0, ss_x(ss, cal)):
        sol = root(f, g, method="hybr", tol=1e-12)
        fn = np.max(np.abs(sol.fun))
        if best is None or fn < best[1]:
            best = (sol.x, fn)
    _, o = point_residuals(S, 0, best[0], rules, cal, ss, sproc,
                           n_gh=ngh, no_default=False)
    return best[0], o


def read_at(rules, cal, ss, sproc, S):
    # READ THE CONVERGED RULES AT STATE S (binding branch), returning the implied
    # allocation + the point residual there (accuracy). Evaluating the period map
    # at the rules' OWN policy values stays on the binding branch -- re-solving
    # with a root finder can slip onto the nearby slack equilibrium at the barely-
    # binding SS. This is the standard way to read a global solution's IRF.
    Sm = np.atleast_2d(S)
    x = np.array([float(rules.eval(k, 0, Sm)[0]) for k in SOLVE7])
    res, o = point_residuals(S, 0, x, rules, cal, ss, sproc,
                             n_gh=rules.n_gh or N_GH, no_default=False)
    o["_x"] = x
    o["_resid"] = float(np.max(np.abs(res)))
    return o


def read_exact(rules, cal, ss, sproc, S, x0=None):
    # THE PERIOD MAP CLEARED EXACTLY AT S, against the same (interpolated) continuation.
    # read_at returns the INTERPOLANT's own values, which is what the collocation
    # solution is and what Bocola's simul.m reads. This clears the 13 period-map
    # residuals at S instead, warm-started at the interpolant.
    #
    # THE TWO DISAGREE BY MORE THAN THE RESPONSE, AND THAT IS THE MODEL'S OWN KKT KINK.
    # mu = max{1 - E[Om]R n / (lambda*assets), 0} is C0, and this economy RESTS ON the
    # kink: mu = 0 exactly at the stochastic rest point. A Chebyshev interpolant cannot
    # represent max(.,0), so near the kink it returns mu > 0 where the truth is 0 -- the
    # fitted read then prices a credit spread that is not there and output falls;
    # cleared exactly, mu stays 0 until p^d ~ 1.5% and output RISES (the deposit rate
    # falls with no spread to offset it). Measured Y_D at p^d = 1.98%: -0.081% fitted
    # against -0.008% cleared. The gap does shrink with resolution (0.087 -> 0.073 pp
    # going from 21 to 95 points) but slowly, as a Gibbs phenomenon does.
    # THIS IS NOT PECULIAR TO THIS MODEL. The same measurement on Bocola's own solved
    # coefficients: his fitted mu policy returns a 28.4 bp liquidity premium on impact
    # where the exact multiplier gives 2.1 bp -- 13x, and it is his published number.
    # Neither read is "the truth"; the honest object is the pair, and every reader
    # prints both so the range cannot hide.
    ngh = rules.n_gh or N_GH
    if x0 is None:
        x0 = np.array([float(rules.eval(k, 0, np.atleast_2d(S))[0]) for k in SOLVE7])

    def f(x):
        try:
            return point_residuals(S, 0, x, rules, cal, ss, sproc, n_gh=ngh,
                                   no_default=False)[0]
        except (ValueError, RuntimeError, ArithmeticError):
            return np.full(len(SOLVE7), 10.0)

    sol = root(f, x0, method="hybr", tol=1e-13)
    _, o = point_residuals(S, 0, sol.x, rules, cal, ss, sproc, n_gh=ngh,
                           no_default=False)
    o["_x"] = sol.x
    o["_resid"] = float(np.max(np.abs(sol.fun)))
    return o


def _spread_bp(o, cal, c="D"):
    # LENDING-SPREAD PROXY IN ANNUALISED BASIS POINTS (lambda_K * mu / alpha).
    # c selects the country: the paper figures plot both, because F is the control --
    # the D shock reaches the F bank only through the union deposit market, so the F
    # line is how much of the D move is a union-wide repricing rather than the shock.
    return 4e4 * cal[f"lambda_K_{c}"] * o[f"mu_{c}"] / max(o[f"alpha_{c}"], 1e-6)


def report_rest_point(rules, cal, ss, sproc):
    # PRINT THE MODEL'S OWN REST POINT NEXT TO THE DETERMINISTIC STEADY STATE.
    # print_ss_table reports the DETERMINISTIC SS -- the object steady_state.py solves,
    # and the grid centre. It is exact: at pi == 0 the model sits on it for 200 quarters
    # to six decimals. But the SOLVED rules price risk, and a risk-pricing economy does
    # not rest where its risk-free counterpart does. Every IRF below is read against the
    # REST POINT, so the two have to be shown together or the steady-state table quietly
    # describes a state the model never visits.
    S0 = ss_state(ss, cal, sproc)
    Sr = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    bd = read_at(rules, cal, ss, sproc, S0.copy())
    br = read_at(rules, cal, ss, sproc, Sr.copy())
    print("\n  DETERMINISTIC SS vs THE MODEL'S OWN REST POINT (the IRF baseline)")
    print(f"   {'object':<22s} {'at the det-SS state':>20s} {'at the rest point':>18s}"
          f" {'diff':>10s}")
    for k, lab in (("mu_D", "mu_D (IC multiplier)"), ("Q_bD", "Q_bD"),
                   ("n_D", "n_D (bank net worth)"), ("Y_D", "Y_D"),
                   ("C_D", "C_D"), ("I_D", "I_D")):
        d = 100.0 * (br[k] / bd[k] - 1.0) if abs(bd[k]) > 1e-12 else float("nan")
        print(f"   {lab:<22s} {bd[k]:20.6f} {br[k]:18.6f} {d:+9.3f}%")
    print(f"   {'credit spread bp/yr':<22s} {_spread_bp(bd, cal):20.1f} "
          f"{_spread_bp(br, cal):18.1f} {_spread_bp(br, cal) - _spread_bp(bd, cal):+9.1f}")
    print(f"   {'rdep_D bp/yr':<22s} {bp_ann(bd['rdep_D']):20.1f} "
          f"{bp_ann(br['rdep_D']):18.1f} {bp_ann(br['rdep_D'] - bd['rdep_D']):+9.1f}")
    print(f"   {'r_wc_D bp/yr':<22s} {bp_ann(bd['r_wc_D']):20.1f} "
          f"{bp_ann(br['r_wc_D']):18.1f} {bp_ann(br['r_wc_D'] - bd['r_wc_D']):+9.1f}")
    dev = ", ".join(f"{_SNAMES[i]} {100 * (Sr[i] / S0[i] - 1):+.3f}%"
                    for i in (IK_D, IP_D, IBDD) if abs(S0[i]) > 1e-12)
    print(f"   states at the rest point: {dev}")
    print("   (the gap is PRICED RISK, not solver error: solved at pi == 0 the model "
          "rests on\n    the deterministic SS exactly -- see CLAUDE.md)")
    return Sr


def impact_table(rules, cal, ss, sproc):
    # RESPONSE AT THE SS-LEVEL STATE AS THE PRICED DEFAULT PROBABILITY RISES.
    # UNITS (reporting.prints): p^d both quarterly and annual; every RATE in annualised
    # basis points; level responses in % with the annualised (Bocola x400) companion for
    # output. The three rate columns are the audit's decomposition of the labour wedge:
    # r_wc = rdep + lambda*mu/E[Om], and it is the FALL in rdep that used to cancel most
    # of the rise in the credit spread before it reached any firm's wage bill.
    # BASELINE AT THE MODEL'S OWN REST POINT, not the deterministic SS. The endogenous
    # states are frozen here (only s moves), so there is no drift to difference away --
    # but the level the deviations are taken from must still be the state the economy
    # inhabits, and the two differ by more than the response being measured
    # (Y_D -0.111%, mu_D 0.0072 -> 0). See stochastic_rest_point.
    S0 = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    base = read_at(rules, cal, ss, sproc, S0.copy())
    Yb, Cb, Ib, Nb = base["Y_D"], base["C_D"], base["I_D"], base["_x"][0]
    print("\n  IMPACT of priced default risk (deviation from the rest point)")
    base_x = read_exact(rules, cal, ss, sproc, S0.copy())
    print("   pd_q%  pd_a%     Y%    Y_ann%   Y_exact%    C%    hours%     I%    "
          "rdep_bp  spread_bp  r_wc_bp    muD   muD_ex     Q_bD   resid")
    s_hi = float(rules.grid.hi[IS])
    # Y_exact / muD_ex clear the period map at the state instead of reading the
    # interpolant -- see read_exact. The pair BRACKETS the response; they agree at the
    # collocation nodes and diverge near the mu = max(.,0) kink, which is where the
    # model's own rest point sits.
    for s_val in (sproc["s_star"], -6.0, -5.2, -4.5, -3.9,
                  0.5 * (-3.9 + s_hi), s_hi):
        S = S0.copy(); S[IS] = s_val
        o = read_at(rules, cal, ss, sproc, S)
        ox = read_exact(rules, cal, ss, sproc, S, o["_x"])
        pq = float(default_prob(s_val))
        dY = 100 * (o["Y_D"] / Yb - 1)
        dYx = 100 * (ox["Y_D"] / base_x["Y_D"] - 1)
        print(f"  {100*pq:6.2f} {100*ann_prob(pq):6.2f} {dY:+8.4f} {ann_pct(dY):+9.4f} "
              f"{dYx:+9.4f} {100*(o['C_D']/Cb-1):+8.4f} {100*(o['_x'][0]/Nb-1):+8.4f} "
              f"{100*(o['I_D']/Ib-1):+7.3f} {bp_ann(o['rdep_D']):8.1f} "
              f"{_spread_bp(o, cal):10.1f} {bp_ann(o['r_wc_D']):8.1f} "
              f"{o['mu_D']:7.5f} {ox['mu_D']:8.5f}  {o['Q_bD']:.4f} {o['_resid']:.0e}")


def s_from_pd(pd):
    # THE RISK STATE s THAT PRICES A ONE-QUARTER-AHEAD DEFAULT PROBABILITY pd.
    # default_prob is the logistic of s, so this is its inverse -- it lets the driver
    # state the shock in the units the result is read in (p^d), not in logit units.
    return float(np.log(pd / (1.0 - pd)))


def persistence_irf(rules, cal, ss, sproc, pd_shock=0.0198, T=21, s_shock=None):
    # IRF AS AN s-SHOCK DECAYS (rho_s), endogenous states held at SS so the path
    # stays on-grid -- the shock-persistence channel (a lower bound; the
    # endogenous net-worth dynamics amplify it). Reads the binding-branch rules.
    # the shock is stated as a TARGET p^d (main.py's RISK_SHOCK_PD); s_shock overrides
    # it in raw logit units. The old default s_shock = -3.9 was labelled "+2 sigma" and
    # is +4.77 sigma at the calibrated sigma_s = 0.63.
    if s_shock is None:
        s_shock = s_from_pd(pd_shock)
    # endogenous states frozen at the REST POINT (see impact_table)
    S0 = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    base = read_at(rules, cal, ss, sproc, S0.copy())
    Yr, Cr, Ir, Nr = base["Y_D"], base["C_D"], base["I_D"], base["_x"][0]
    print(f"\n  PERSISTENCE IRF (one-off risk shock to p^d = "
          f"{100*default_prob(s_shock):.2f}%/qtr, rho_s = {sproc['rho_s']} decay)")
    print("   qtr  pd_q%  pd_a%     Y%    Y_ann%      C%     I%   hours%  "
          "rdep_bp  spread_bp  r_wc_bp   Q_bD")
    path = {k: [] for k in ("pd", "pd_ann", "Y", "Y_ann", "C", "I", "N", "spread",
                            "rdep", "r_wc", "Q_bD")}
    for t in range(T):
        s_t = sproc["s_star"] + sproc["rho_s"] ** t * (s_shock - sproc["s_star"])
        S = S0.copy(); S[IS] = s_t
        o = read_at(rules, cal, ss, sproc, S)
        pq = float(default_prob(s_t))
        dY = 100 * (o["Y_D"] / Yr - 1)
        path["pd"].append(100 * pq); path["pd_ann"].append(100 * ann_prob(pq))
        path["Y"].append(dY); path["Y_ann"].append(ann_pct(dY))
        path["C"].append(100 * (o["C_D"] / Cr - 1))
        path["I"].append(100 * (o["I_D"] / Ir - 1))
        path["N"].append(100 * (o["_x"][0] / Nr - 1))
        path["spread"].append(_spread_bp(o, cal))
        path["rdep"].append(bp_ann(o["rdep_D"]))
        path["r_wc"].append(bp_ann(o["r_wc_D"]))
        path["Q_bD"].append(o["Q_bD"])
        if t in (0, 1, 2, 4, 6, 8, 12, 16, 20):
            print(f"   {t:3d} {path['pd'][-1]:6.2f} {path['pd_ann'][-1]:6.2f} "
                  f"{path['Y'][-1]:+8.4f} {path['Y_ann'][-1]:+9.4f} "
                  f"{path['C'][-1]:+8.4f} {path['I'][-1]:+7.3f} {path['N'][-1]:+7.3f} "
                  f"{path['rdep'][-1]:8.1f} {path['spread'][-1]:10.1f} "
                  f"{path['r_wc'][-1]:8.1f}  {path['Q_bD'][-1]:.4f}")
    return {k: np.array(v) for k, v in path.items()}


def advance(o, S, sproc, grid=None):
    # ONE STEP OF THE MODEL'S OWN LAW OF MOTION FROM A READ `o` AT STATE S.
    # Single-sourced: dynamic_irf, stochastic_rest_point and the no-shock reference
    # path must advance the state IDENTICALLY, or the difference between a shocked and
    # an unshocked path picks up the discrepancy instead of the shock.
    Sn = S.copy()
    Sn[IK_D], Sn[IK_F] = o["_x"][2], o["_x"][3]      # K' = Kp
    Sn[IP_D], Sn[IP_F] = o["Pp_D"], o["Pp_F"]
    Sn[IBDD], Sn[IBDF] = o["b_D_D_new"], o["b_D_F_new"]
    Sn[IBFD] = o["b_F_D_new"]
    Sn[IV] = o["Vp_dep"]
    Sn[IS] = (1 - sproc["rho_s"]) * sproc["s_star"] + sproc["rho_s"] * S[IS]
    return Sn if grid is None else grid.clip(Sn)[0]


def stochastic_rest_point(rules, cal, ss, sproc, tol=1e-11, max_it=4000, verbose=True):
    # THE STATE THE SOLVED MODEL ACTUALLY RESTS AT -- Bocola's generate_irf.m step 1
    # ("e = zeros(2000,3); [state,obs,STATE] = simul(...); initial = STATE(:,end-1)").
    #
    # WHY THIS IS NEEDED AND IS NOT A BUG. The DETERMINISTIC steady state is an exact
    # rest point of the period map at pi = 0: solved at pi == 0, the model sits on it
    # for 200 quarters to six decimals with mu pinned at 0.001001 and max|F| ~ 1e-7
    # (measured). But the SOLVED rules PRICE RISK, and a risk-pricing economy does not
    # rest where its risk-free counterpart does -- Bocola's own solution has the same
    # gap (his ergodic q = 0.979 against a deterministic 1.000, debt +2.0%). Measured
    # here: Y_D -0.111%, C_D -0.130%, I_D +0.246%, n_D +2.24%, K_D -0.541%,
    # b_DD +2.87%, and mu_D falls from 0.0072 to EXACTLY 0 -- the constraint is SLACK
    # at the point the model inhabits, so the calibrated 8 bp steady-state credit
    # spread is not a property of the ergodic economy (nor is it in Bocola's: his
    # constraint binds on 1.2% of his ergodic set).
    #
    # WHAT IT COSTS TO IGNORE IT. Starting an IRF at the deterministic SS and
    # differencing against a FIXED base charges that walk to the shock. Measured at
    # this calibration: at q12 the GDP response reads +0.0651% where the true
    # (differenced) response is +0.0300% -- 54% of the reported hump was drift -- and
    # bank net worth reads +4.20% against a true +1.55%.
    cached = getattr(rules, "_rest_point", None)
    if cached is not None:
        return cached.copy()
    S = ss_state(ss, cal, sproc)
    d = np.inf
    for it in range(1, max_it + 1):
        o = read_at(rules, cal, ss, sproc, S)
        Sn = advance(o, S, sproc, rules.grid)
        d = float(np.max(np.abs(Sn - S)))
        S = Sn
        if d < tol:
            break
    if verbose:
        S0 = ss_state(ss, cal, sproc)
        dev = ", ".join(f"{_SNAMES[i]} {100 * (S[i] / S0[i] - 1):+.3f}%"
                        for i in (IK_D, IP_D, IBDD) if abs(S0[i]) > 1e-12)
        print(f"  stochastic rest point: {it} no-shock quarters, |dS| = {d:.1e}"
              f"   (vs the deterministic SS: {dev})")
        if d > 1e-8:
            print(f"    WARNING: the no-shock path had not settled ({d:.1e} > 1e-8). "
                  "The IRF below is still differenced, so it is valid, but the "
                  "starting point is not the model's rest point.")
    rules._rest_point = S.copy()
    return S


def dynamic_irf(rules, cal, ss, sproc, pd_shock=0.0198, T=25, rest_verbose=True):
    # DYNAMIC IRF: THE STATE VECTOR ITERATES FORWARD, IT IS NOT HELD AT THE SS.
    # persistence_irf varies only s and pins K, P and B at their steady-state values --
    # its own docstring calls that "a lower bound". Bocola's Table 5 output losses are
    # the LEVEL of output over six quarters, driven by capital and bank net worth
    # accumulating downward, so a frozen-state impact reading cannot be compared with
    # them. Here every endogenous state follows the period map's own law of motion
    # (K' = Kp, P' = Pp, B' = Bp) while s decays at rho_s, which is the object his
    # numbers describe.
    # BOCOLA'S generate_irf.m, BOTH HALVES. (1) Start at the STOCHASTIC rest point,
    # not the deterministic SS -- the two differ because the solved rules price risk
    # (see stochastic_rest_point). (2) Difference the shocked path against an UNSHOCKED
    # path from the SAME state, rather than against a frozen base: his
    # `gdp = mean(gdp_s) - mean(gdp_nos)`. Either alone removes most of the artifact;
    # he does both, so the reported response cannot contain the no-shock transition
    # even if the rest point is imperfectly converged.
    S0 = stochastic_rest_point(rules, cal, ss, sproc, verbose=rest_verbose)
    S = S0.copy()
    escapes = []
    base = read_at(rules, cal, ss, sproc, S0.copy())
    # the NO-SHOCK reference path, advanced with the same law of motion
    ref, Sr = [], S0.copy()
    for _ in range(T):
        o_r = read_at(rules, cal, ss, sproc, Sr)
        ref.append(o_r)
        Sr = advance(o_r, Sr, sproc, rules.grid)
    S[IS] = s_from_pd(pd_shock)
    print(f"\n  DYNAMIC IRF (states evolve; p^d shock to {100*pd_shock:.2f}%/qtr, "
          f"rho_s = {sproc['rho_s']})")
    print("   qtr  pd_q%  pd_a%    GDP%   GDP_ann%     C%      I%   hours%  "
          "rdep_bp  bank_bp  r_wc_bp  sov_bp    K%      n%")
    # TWO SPREADS, RECORDED SEPARATELY BECAUSE THEY ARE DIFFERENT OBJECTS. "spread" is the
    # BANK CREDIT spread lambda_K*mu/alpha, identically zero once mu hits the KKT switch;
    # "sov_bp" is the SOVEREIGN spread y_D - y_F out of the bond Euler, which carries no mu
    # and persists for as long as p^d is elevated. The figure plots the sovereign one.
    # THE F COUNTERPARTS (Y_F, C_F, n_F, spread_F) ARE RECORDED, NOT DERIVED: the paper
    # figures plot both countries, and the F line is the control for the D response --
    # the shock is D's alone and reaches F only through the union deposit market.
    path = {k: [] for k in ("pd", "pd_ann", "Y", "Y_ann", "C", "I", "N", "spread",
                            "rdep", "r_wc", "K", "n", "Q_bD", "Q_bF", "sov_bp",
                            "d_rdep", "d_spread", "d_r_wc", "dQ_bD", "m_ltro",
                            "mu", "E_Om", "Y_F", "C_F", "n_F", "spread_F", "I_F")}
    # THE WEDGE DECOMPOSITION, in annualised bp DEVIATIONS from the no-shock state.
    # r_wc = rdep + lambda*mu/E[Om] is the only channel from the financial block into
    # output under GHH, and its two legs move in OPPOSITE directions: the credit spread
    # rises with the constraint while the union deposit rate falls as banks delever.
    # Reporting only the spread hides the netting -- which is what the 2026-08-28 audit
    # found was costing most of the output response before country size was made
    # asymmetric.
    Sr = S0.copy()                    # state of the no-shock path, for K
    K_ref = [Sr[IK_D]]
    for t in range(T - 1):
        Sr = advance(ref[t], Sr, sproc, rules.grid)
        K_ref.append(Sr[IK_D])
    for t in range(T):
        o = read_at(rules, cal, ss, sproc, S)
        r = ref[t]                    # the SAME quarter of the unshocked path
        pq = float(default_prob(S[IS]))
        dY = 100 * (o["Y_D"] / r["Y_D"] - 1)
        path["pd"].append(100 * pq); path["pd_ann"].append(100 * ann_prob(pq))
        path["Y"].append(dY); path["Y_ann"].append(ann_pct(dY))
        path["C"].append(100 * (o["C_D"] / r["C_D"] - 1))
        path["I"].append(100 * (o["I_D"] / r["I_D"] - 1))
        path["I_F"].append(100 * (o["I_F"] / r["I_F"] - 1))
        path["N"].append(100 * (o["_x"][0] / r["_x"][0] - 1))
        path["spread"].append(_spread_bp(o, cal))
        path["Y_F"].append(100 * (o["Y_F"] / r["Y_F"] - 1))
        path["C_F"].append(100 * (o["C_F"] / r["C_F"] - 1))
        path["n_F"].append(100 * (o["n_F"] / r["n_F"] - 1))
        path["spread_F"].append(_spread_bp(o, cal, "F"))
        path["rdep"].append(bp_ann(o["rdep_D"]))
        path["r_wc"].append(bp_ann(o["r_wc_D"]))
        path["d_rdep"].append(bp_ann(o["rdep_D"] - r["rdep_D"]))
        path["d_spread"].append(_spread_bp(o, cal) - _spread_bp(r, cal))
        path["d_r_wc"].append(bp_ann(o["r_wc_D"] - r["r_wc_D"]))
        path["K"].append(100 * (S[IK_D] / K_ref[t] - 1))
        path["n"].append(100 * (o["n_D"] / r["n_D"] - 1))
        path["Q_bD"].append(o["Q_bD"])
        path["Q_bF"].append(o["Q_bF"])
        path["dQ_bD"].append(100 * (o["Q_bD"] / r["Q_bD"] - 1))
        # THE BACKSTOP'S FOOTPRINT AND ITS TWO OPPOSING CHANNELS. m_ltro is the facility
        # actually drawn -- ZERO along the never-fired path, which is the headline read.
        # mu and E_Om are recorded because they move in OPPOSITE directions: the facility
        # relieves the constraint directly, but by lowering the franchise value it lowers
        # E[Om], which RAISES mu. Reporting only the first would hide the offset that
        # decides the sign (docs/ltro_backstop_plan.md S3).
        path["m_ltro"].append(100 * o["m_ltro_D"] / ss["ss_firm_D"]["Y_ss"])
        path["mu"].append(o["mu_D"])
        path["E_Om"].append(o["E_Om_D"])
        _yD = cal["delta_b_D"] * (1.0 - o["Q_bD"]) / o["Q_bD"]   # HM perpetuity flow yield
        _yF = cal["delta_b_F"] * (1.0 - o["Q_bF"]) / o["Q_bF"]
        path["sov_bp"].append(bp_ann(_yD - _yF))
        if t in (0, 1, 2, 4, 6, 8, 12, 16, 20, 24):
            print(f"   {t:3d} {path['pd'][-1]:6.2f} {path['pd_ann'][-1]:6.2f} "
                  f"{path['Y'][-1]:+8.4f} {path['Y_ann'][-1]:+9.4f} "
                  f"{path['C'][-1]:+8.4f} {path['I'][-1]:+7.3f} {path['N'][-1]:+7.3f} "
                  f"{path['rdep'][-1]:8.1f} {path['spread'][-1]:8.1f} "
                  f"{path['r_wc'][-1]:8.1f} {path['sov_bp'][-1]:7.1f} "
                  f"{path['K'][-1]:+7.3f} {path['n'][-1]:+7.2f}")
        Sn = advance(o, S, sproc)
        # BOX ESCAPES ARE REPORTED, NOT SWALLOWED. Clipping keeps the read on-grid, but a
        # state pinned to a band turns a divergent law of motion into a flat IRF that
        # looks like convergence: the pre-2026-08-25 figure's "trough at q7 then recovery"
        # was B_D pinned at +8.00% from q7 to q24 (debt root 0.9929; see phi_lamb).
        S = rules.grid.clip(Sn)[0]
        esc = np.abs(Sn - S) > 1e-12
        if esc.any():
            escapes.append((t, [(_SNAMES[i], 100 * (Sn[i] / S0[i] - 1),
                                 100 * (S[i] / S0[i] - 1)) for i in np.flatnonzero(esc)]))
    # THE BENCHMARK, IN BOTH UNITS. The line this replaces compared a LEVEL IRF against
    # Bocola's Table 5, which is a cumulated quarterly GROWTH gap x400 -- four times too
    # demanding, and a different experiment (an 8-quarter estimated shock sequence, not
    # one shock). His single-shock IRFs, rescaled to this p^d, are the right targets.
    tr = min(path["Y"])
    print(f"   trough GDP = {tr:+.4f}% level = {ann_pct(tr):+.4f}% annualised "
          f"(Bocola Table 5 units)")
    # THE SAME IMPACT, CLEARED EXACTLY. The interpolant and the exactly-cleared period
    # map bracket the response, and near the mu = max(.,0) kink -- where this model's
    # rest point sits -- the bracket is wide. Printing one number alone would be a
    # false precision. See read_exact.
    Sx = S0.copy(); Sx[IS] = s_from_pd(pd_shock)
    bx = read_exact(rules, cal, ss, sproc, S0.copy())
    ox = read_exact(rules, cal, ss, sproc, Sx)
    trx = 100 * (ox["Y_D"] / bx["Y_D"] - 1)
    print(f"   impact GDP brackets [{min(path['Y'][0], trx):+.4f}%, "
          f"{max(path['Y'][0], trx):+.4f}%]: {path['Y'][0]:+.4f}% reading the fitted "
          f"rules, {trx:+.4f}% clearing the period map exactly at that state "
          f"(mu {bx['mu_D']:.5f} -> {ox['mu_D']:.5f})")
    # THE IDENTITY THAT PRODUCES THAT TROUGH, on impact. Under GHH,
    #   dlogY = -(1-alpha)/(1/nu+alpha) * [dlog(1+zeta*r_wc) + dlog P_CES],
    # so the output response IS the working-capital wedge response, and the wedge is
    # the credit spread NET of the deposit rate.
    print(f"   impact wedge: credit spread {path['d_spread'][0]:+.1f} bp/yr, "
          f"deposit rate {path['d_rdep'][0]:+.1f} bp/yr, "
          f"NET r_wc {path['d_r_wc'][0]:+.1f} bp/yr "
          f"({100*path['d_rdep'][0]/max(abs(path['d_spread'][0]), 1e-12):+.0f}% of the "
          f"spread is cancelled by the funding rate)")
    print(f"   like-for-like targets at this shock: {BOCOLA_IRF_OPEN:+.3f}% "
          f"(his open economy -- GHH + working capital, our own structure), "
          f"{BOCOLA_IRF_CLOSED:+.3f}% (his closed benchmark), "
          f"{BOCOLA_EPISODE_LEVEL:+.3f}% (2011Q4 episode)")
    if escapes:
        first = escapes[0]
        print(f"   WARNING: {len(escapes)}/{T} quarters left the collocation box "
              f"(first at q{first[0]}). The path beyond it is the BOX WALL, not the model.")
        for t_e, items in escapes[:3]:
            for nm, want, got in items:
                print(f"     q{t_e:<3d} {nm:<4s} law of motion {want:+7.2f}%  ->  clipped "
                      f"{got:+7.2f}%")
        if len(escapes) > 3:
            print(f"     ... and {len(escapes) - 3} more quarters")
    else:
        print(f"   box: no escapes in {T} quarters "
              f"(final |dev| from SS: " +
              ", ".join(f"{_SNAMES[i]} {100 * (S[i] / S0[i] - 1):+.2f}%"
                        for i in (0, 2, 3, 4)) + ")")
    return {k: np.array(v) for k, v in path.items()}


def _tfp_read(rules, cal, ss, sproc, S):
    # READ THE NO-DEFAULT RULES AT STATE S (TFP experiment: no sovereign risk).
    Sm = np.atleast_2d(S)
    x = np.array([float(rules.eval(k, 0, Sm)[0]) for k in SOLVE7])
    res, o = point_residuals(S, 0, x, rules, cal, ss, sproc,
                             n_gh=rules.n_gh or N_GH, no_default=True)
    o["_x"] = x
    return o


def solve_tfp(cal, ss, sproc, mu=1):
    # SOLVE THE NO-DEFAULT (d=0) RULES FOR THE TFP EXPERIMENT, SAME GLOBAL NEWTON.
    # No s-refinement here: with pi = 0 the risk dimension carries no curvature, and
    # the TFP state Z_D enters the period map linearly through the production function.
    grid = build_state_box(ss, cal, mu=mu, **BOX_KW)
    rules = RuleSet.from_ss(grid, ss, cal)
    rules.n_gh = N_GH
    _stage(rules, cal, ss, sproc, (0,), True, "TFP d0", True)
    return rules


def tfp_irf(rules, cal, ss, sproc, dz=0.01, T=21):
    # TFP IRF read off the no-default rules along the Z_D-decay path (rho_z from
    # sproc), endogenous states held at SS so the read stays on-grid -- the exact
    # image of persistence_irf, with the TFP state Z_D in place of the risk state s.
    S0 = ss_state(ss, cal, sproc)
    Z_ss = S0[IZ]
    base = _tfp_read(rules, cal, ss, sproc, S0.copy())
    Yb, Cb, Ib, Nb = base["Y_D"], base["C_D"], base["I_D"], base["_x"][0]
    YbF, CbF, IbF = base["Y_F"], base["C_F"], base["I_F"]
    nbF, nbD = base["n_F"], base["n_D"]
    print(f"\n  TFP IRF (one-off {dz:.0%} shock, rho_z={sproc['rho_z']} decay)")
    print("   qtr    Z%     Y_D%     C_D%     I_D%    hours%")
    # same five paper series as dynamic_irf, so both figures read off one panel spec
    path = {k: [] for k in ("Z", "Y", "C", "I", "N", "Y_F", "C_F", "I_F", "n", "n_F",
                            "spread", "spread_F", "Q_bD", "Q_bF")}
    for t in range(T):
        z_t = dz * sproc["rho_z"] ** t
        S = S0.copy(); S[IZ] = Z_ss * np.exp(z_t)
        o = _tfp_read(rules, cal, ss, sproc, S)
        path["Z"].append(100 * z_t)
        path["Y"].append(100 * (o["Y_D"] / Yb - 1))
        path["C"].append(100 * (o["C_D"] / Cb - 1))
        path["I"].append(100 * (o["I_D"] / Ib - 1))
        path["I_F"].append(100 * (o["I_F"] / IbF - 1))
        path["N"].append(100 * (o["_x"][0] / Nb - 1))
        path["Y_F"].append(100 * (o["Y_F"] / YbF - 1))
        path["C_F"].append(100 * (o["C_F"] / CbF - 1))
        path["n"].append(100 * (o["n_D"] / nbD - 1))
        path["n_F"].append(100 * (o["n_F"] / nbF - 1))
        path["spread"].append(_spread_bp(o, cal))
        path["spread_F"].append(_spread_bp(o, cal, "F"))
        path["Q_bD"].append(o["Q_bD"])
        path["Q_bF"].append(o["Q_bF"])
        if t in (0, 1, 2, 4, 6, 8, 12, 16, 20):
            print(f"   {t:3d}  {path['Z'][-1]:5.2f}  {path['Y'][-1]:+7.3f}  "
                  f"{path['C'][-1]:+7.3f}  {path['I'][-1]:+7.2f}  {path['N'][-1]:+6.2f}")
    return {k: np.array(v) for k, v in path.items()}


def tfp_main():
    cal = get_calibration()
    cal["nw_floor_frac"] = 0.15      # match main.py
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print("=== TFP shock — recursive projection (Z_D as the 7th state) ===")
    rules = solve_tfp(cal, ss, sproc)
    tfp_irf(rules, cal, ss, sproc)


def main():
    cal = get_calibration()
    cal["nw_floor_frac"] = 0.15      # match main.py
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print("=== recursive global solution: pass-through of sovereign risk ===")
    rules = solve_recursive(cal, ss, sproc)
    try:
        import pickle
        with open("/private/tmp/claude-501/-Users-Huawei-Quantitative-Model/"
                  "239042af-4c74-4ebe-a83f-92681158d4c3/scratchpad/mu2_rules.pkl",
                  "wb") as fh:
            pickle.dump((rules, cal, ss, sproc), fh)
    except Exception:
        pass
    impact_table(rules, cal, ss, sproc)
    persistence_irf(rules, cal, ss, sproc)


if __name__ == "__main__":
    main()
