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

from calibration.global_projection import get_calibration
from global_projection.steady_state import solve_steady_state
from global_projection.solver_recursive.state_grid import build_state_box, s_process_params, default_prob
from global_projection.solver_recursive.state_grid import (IK_D, IK_F, IP_D, IP_F, IBDD,
                                          IBDF, IBFD, IV, IS, IZ, STATE_NAMES)
from global_projection.solver_recursive.decision_rules import RuleSet, SOLVE7, STORE_RULES
from global_projection.solver_recursive.collocation import solve_collocation
from global_projection.solver_recursive.recursive_main import (time_iteration, calibrate_household_anchors,
                            ss_state, p_block_rotation)
from global_projection.reporting.prints import (bp_ann, ann_pct, ann_prob,
                              BOCOLA_IRF_CLOSED, BOCOLA_IRF_OPEN,
                              BOCOLA_EPISODE_LEVEL)
from global_projection.solver_recursive.point_map import point_residuals


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
           backend="auto", warm=WARM_SWEEPS, maxit=40):
    # ONE SOLVE STAGE: a short time-iteration warm start, then the GLOBAL NEWTON.
    if warm:
        time_iteration(rules, cal, ss, sproc, regimes=regimes, no_default=no_default,
                       damp=0.5, tol=1e-4, max_it=warm, n_gh=N_GH, verbose=False)
    return solve_collocation(rules, cal, ss, sproc, regimes=regimes,
                             no_default=no_default, n_gh=N_GH, backend=backend,
                             maxit=maxit, verbose=verbose, label=f" {label}")


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
    # This is the ceiling on channel (a) in Claude files/docs/ltro_backstop_plan.md S3, and it is why
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


def _build_rules(cal, ss, sproc, mu, mu_vec, rotate, nreg, verbose):
    # THE COARSE COLLOCATION BOX (OPTIONALLY P-ROTATED) AND AN SS COLD-START RULE SET.
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
    grid = build_state_box(ss, cal, mu=mu, mu_vec=mu_vec, rot=rot, centre=centre,
                           **box_kw)
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=nreg)
    rules.n_gh = N_GH
    if verbose:
        print(f"  coarse grid: mu={mu}, {grid.n} points x {nreg} regimes, n_gh={N_GH}")
    return rules, dict(mu=mu, mu_vec=mu_vec, rot=rot, centre=centre, box_kw=box_kw)


def _solve_baseline(rules, cal, ss, sproc, D_REG, backend, verbose):
    # THE BASELINE: d = 0, THE HAIRCUT HOMOTOPY, THEN THE JOINT RISK SOLVE.
    # Bocola's ladder: he warm-starts the default model from the solved no-default one
    # and walks the haircut up re-solving at each step.
    ok0, it0, w0 = _stage(rules, cal, ss, sproc, (0,), True, "d0", verbose)
    for k in STORE_RULES:                   # warm-start the default regime from d=0
        rules.set_values(k, D_REG, rules.vals[k][0].copy())
    rec_target = cal["recovery_rate_D"]
    ok1, w1 = False, np.nan
    for rec in (0.85, 0.70, 0.55, rec_target):
        cal["recovery_rate_D"] = rec
        ok1, it1, w1 = _stage(rules, cal, ss, sproc, (D_REG,), False,
                              f"d1 rec={rec:.2f}", verbose)
    cal["recovery_rate_D"] = rec_target
    okb, itb, wb = _stage(rules, cal, ss, sproc, (0, D_REG), False,
                          "joint (no facility)", verbose, backend=backend)
    return (ok0, ok1, okb), (w0, w1, wb)


def _refine_s(rules, cal, ss, sproc, box, s_refine, backend, verbose):
    # THE REFINEMENT LADDER: A DENSE CHEBYSHEV FACTOR IN s, SEEDED FROM THE LAST SOLVE.
    # Same box each time -- one more continuation step in Bocola's style. Going straight
    # from 3 nodes in s to 9 asks the Newton to start from a seed that is badly wrong at
    # the new interior nodes (the mu=1 quadratic reads p^d = 1.82% where the truth is
    # 0.67%), so the node count is walked up.
    # THE WARM START IS NOT OPTIONAL HERE: the refined grid visits (P_D at +1, s at +1)
    # combinations the coarse interpolant has no cross term for. Measured: the raw seed
    # sits at max|F| = 4.0e-2, and 20 time-iteration sweeps bring it to 1.9e-3 -- inside
    # the Newton's basin.
    nreg = rules.n_regimes
    okj = wj = None
    ladder = [m for m in (5, 9, 17) if 1 < m < s_refine] + [s_refine]
    for m_s in ladder:
        gfine = build_state_box(ss, cal, mu=box["mu"], mu_vec=box["mu_vec"],
                                rot=box["rot"], centre=box["centre"], refine=(IS, m_s),
                                **box["box_kw"])
        if verbose:
            print(f"  s-refined grid: {gfine.n} points x {nreg} regimes "
                  f"({m_s} nodes, degree {m_s - 1} in s)")
        fine = _seed_from(RuleSet(gfine, nreg), rules)
        okj, itj, wj = _stage(fine, cal, ss, sproc, tuple(range(nreg)), False,
                              f"joint (s={m_s})", verbose, backend=backend,
                              warm=REFINE_WARM_SWEEPS, maxit=12)
        rules = fine
    return rules, okj, wj


def _stamp_verdict(rules, oks, worsts):
    # RECORD WHETHER EVERY STAGE ROOTED ON THE RULE SET ITSELF.
    # A caller that plots or tabulates several solves has no other way to know which of
    # them actually rooted -- without this the certainty curve would draw a stopped solve
    # with the same solid marker as a converged one.
    ok0, ok1, okb, okj = oks
    w0, w1, wb, wj = worsts
    rules.solve_ok = bool(ok0 and ok1 and okb and okj)
    rules.solve_worst = float(np.nanmax([w0, w1, wb, wj]))
    if not rules.solve_ok:
        print(f"    NOTE: a collocation stage did not reach the acceptance floor: "
              f"d0 ok={ok0} ({w0:.1e}), d1 ok={ok1} ({w1:.1e}), "
              f"base ok={okb} ({wb:.1e}), joint ok={okj} ({wj:.1e}). "
              f"Read the IRFs as indicative.")


def solve_recursive(cal, ss, sproc, mu=1, verbose=True, mu_vec=None, rotate=False,
                    s_refine=S_REFINE, backend="auto", base=None, base_out=None):
    # SOLVE EVERY REGIME BY GLOBAL COLLOCATION (Bocola's model_solution_mean.m).
    # THE LADDER IS HIS -- every stage is a genuine root of the collocation system:
    #   1. coarse grid, d = 0 at pi = 0                               (_solve_baseline)
    #   2. the default regime by haircut homotopy 0.85 -> recovery_rate_D
    #   3. the risk-priced joint solve
    #   4. joint polish over every regime
    #   5. rebuild on the s-refined grid, seeded from the coarse one   (_refine_s)
    # A sweep over a policy parameter hands the solved baseline back in through `base`
    # (and collects it through `base_out`).
    nreg = 2
    rules, box = _build_rules(cal, ss, sproc, mu, mu_vec, rotate, nreg, verbose)
    D_REG = 1                                   # the regime index is the default indicator

    if base is None:
        (ok0, ok1, okb), (w0, w1, wb) = _solve_baseline(rules, cal, ss, sproc, D_REG,
                                                        backend, verbose)
        if base_out is not None:                # hand the caller a reusable snapshot
            base_out.append(rules.copy())
    else:
        assert base.n_regimes == nreg and base.grid.n == rules.grid.n, \
            "reused baseline must carry the same grid and regime count"
        rules = base.copy()
        ok0 = ok1 = okb = True
        w0 = w1 = wb = np.nan
        if verbose:
            print("  reusing the solved baseline")

    okj, itj, wj = _stage(rules, cal, ss, sproc, tuple(range(nreg)), False, "joint",
                          verbose, backend=backend)

    if s_refine and s_refine > 1:
        rules, okj, wj = _refine_s(rules, cal, ss, sproc, box, s_refine, backend, verbose)

    _stamp_verdict(rules, (ok0, ok1, okb, okj), (w0, w1, wb, wj))
    return rules


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
    # IRF AS AN s-SHOCK DECAYS (rho_s), endogenous states held at the rest point so the
    # path stays on-grid -- the shock-persistence channel (a lower bound; the endogenous
    # net-worth dynamics amplify it). The shock is stated as a TARGET p^d; s_shock
    # overrides it in raw logit units.
    if s_shock is None:
        s_shock = s_from_pd(pd_shock)
    S0 = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    base = read_at(rules, cal, ss, sproc, S0.copy())
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
        _append(path, _persistence_row(o, base, s_t, cal))
        if t in (0, 1, 2, 4, 6, 8, 12, 16, 20):
            print(f"   {t:3d} {path['pd'][-1]:6.2f} {path['pd_ann'][-1]:6.2f} "
                  f"{path['Y'][-1]:+8.4f} {path['Y_ann'][-1]:+9.4f} "
                  f"{path['C'][-1]:+8.4f} {path['I'][-1]:+7.3f} {path['N'][-1]:+7.3f} "
                  f"{path['rdep'][-1]:8.1f} {path['spread'][-1]:10.1f} "
                  f"{path['r_wc'][-1]:8.1f}  {path['Q_bD'][-1]:.4f}")
    return {k: np.array(v) for k, v in path.items()}


def _append(path, row):
    # APPEND ONE QUARTER'S VALUES TO EVERY SERIES OF AN IRF PATH.
    for k, v in row.items():
        path[k].append(v)


def _persistence_row(o, base, s_t, cal):
    # ONE QUARTER OF THE PERSISTENCE IRF, AS DEVIATIONS FROM THE FROZEN REST-POINT READ.
    pq = float(default_prob(s_t))
    dY = 100 * (o["Y_D"] / base["Y_D"] - 1)
    return {"pd": 100 * pq, "pd_ann": 100 * ann_prob(pq), "Y": dY, "Y_ann": ann_pct(dY),
            "C": 100 * (o["C_D"] / base["C_D"] - 1),
            "I": 100 * (o["I_D"] / base["I_D"] - 1),
            "N": 100 * (o["_x"][0] / base["_x"][0] - 1),
            "spread": _spread_bp(o, cal), "rdep": bp_ann(o["rdep_D"]),
            "r_wc": bp_ann(o["r_wc_D"]), "Q_bD": o["Q_bD"]}


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


def _no_shock_path(rules, cal, ss, sproc, S0, T):
    # THE UNSHOCKED REFERENCE PATH FROM S0: THE READS AND THE CAPITAL STATE PER QUARTER.
    ref, Sr = [], S0.copy()
    for _ in range(T):
        o_r = read_at(rules, cal, ss, sproc, Sr)
        ref.append(o_r)
        Sr = advance(o_r, Sr, sproc, rules.grid)
    Sr = S0.copy()
    K_ref = [Sr[IK_D]]
    for t in range(T - 1):
        Sr = advance(ref[t], Sr, sproc, rules.grid)
        K_ref.append(Sr[IK_D])
    return ref, K_ref


def _dynamic_row(o, r, S, K_ref_t, cal, ss):
    # ONE QUARTER OF THE DYNAMIC IRF: SHOCKED READ o AGAINST THE SAME QUARTER r UNSHOCKED.
    # TWO SPREADS, DIFFERENT OBJECTS: "spread" is the BANK CREDIT spread lambda_K*mu/alpha
    # (zero once mu hits the KKT switch); "sov_bp" is the SOVEREIGN spread y_D - y_F out
    # of the bond Euler, which persists while p^d is elevated. The d_* legs are the wedge
    # decomposition in annualised bp deviations: r_wc = rdep + lambda*mu/E[Om] is the only
    # financial channel into output under GHH and its two legs move in OPPOSITE directions.
    # mu and E_Om are recorded because a backstop moves them in opposite directions.
    pq = float(default_prob(S[IS]))
    dY = 100 * (o["Y_D"] / r["Y_D"] - 1)
    _yD = cal["delta_b_D"] * (1.0 - o["Q_bD"]) / o["Q_bD"]   # HM perpetuity flow yield
    _yF = cal["delta_b_F"] * (1.0 - o["Q_bF"]) / o["Q_bF"]
    return {"pd": 100 * pq, "pd_ann": 100 * ann_prob(pq), "Y": dY, "Y_ann": ann_pct(dY),
            "C": 100 * (o["C_D"] / r["C_D"] - 1),
            "I": 100 * (o["I_D"] / r["I_D"] - 1),
            "I_F": 100 * (o["I_F"] / r["I_F"] - 1),
            "N": 100 * (o["_x"][0] / r["_x"][0] - 1),
            "spread": _spread_bp(o, cal),
            "Y_F": 100 * (o["Y_F"] / r["Y_F"] - 1),
            "C_F": 100 * (o["C_F"] / r["C_F"] - 1),
            "n_F": 100 * (o["n_F"] / r["n_F"] - 1),
            "spread_F": _spread_bp(o, cal, "F"),
            "rdep": bp_ann(o["rdep_D"]), "r_wc": bp_ann(o["r_wc_D"]),
            "d_rdep": bp_ann(o["rdep_D"] - r["rdep_D"]),
            "d_spread": _spread_bp(o, cal) - _spread_bp(r, cal),
            "d_r_wc": bp_ann(o["r_wc_D"] - r["r_wc_D"]),
            "K": 100 * (S[IK_D] / K_ref_t - 1),
            "n": 100 * (o["n_D"] / r["n_D"] - 1),
            "Q_bD": o["Q_bD"], "Q_bF": o["Q_bF"],
            "dQ_bD": 100 * (o["Q_bD"] / r["Q_bD"] - 1),
            "mu": o["mu_D"], "E_Om": o["E_Om_D"],
            "sov_bp": bp_ann(_yD - _yF)}


def _print_dynamic_summary(path, rules, cal, ss, sproc, S0, S, pd_shock, T, escapes):
    # TROUGH, THE FITTED-vs-EXACT IMPACT BRACKET, THE WEDGE NETTING, BENCHMARKS, BOX.
    # The benchmark is in both units: Bocola's Table 5 is a cumulated quarterly GROWTH
    # gap x400 over an 8-quarter shock sequence; his single-shock IRFs, rescaled to this
    # p^d, are the like-for-like targets.
    tr = min(path["Y"])
    print(f"   trough GDP = {tr:+.4f}% level = {ann_pct(tr):+.4f}% annualised "
          f"(Bocola Table 5 units)")
    # the same impact cleared exactly: near the mu = max(.,0) kink, where this model's
    # rest point sits, the fitted and exact reads bracket the response (see read_exact)
    Sx = S0.copy(); Sx[IS] = s_from_pd(pd_shock)
    bx = read_exact(rules, cal, ss, sproc, S0.copy())
    ox = read_exact(rules, cal, ss, sproc, Sx)
    trx = 100 * (ox["Y_D"] / bx["Y_D"] - 1)
    print(f"   impact GDP brackets [{min(path['Y'][0], trx):+.4f}%, "
          f"{max(path['Y'][0], trx):+.4f}%]: {path['Y'][0]:+.4f}% reading the fitted "
          f"rules, {trx:+.4f}% clearing the period map exactly at that state "
          f"(mu {bx['mu_D']:.5f} -> {ox['mu_D']:.5f})")
    # under GHH dlogY = -(1-alpha)/(1/nu+alpha) * [dlog(1+zeta*r_wc) + dlog P_CES], so
    # the output response IS the wedge response, and the wedge is the credit spread NET
    # of the deposit rate
    print(f"   impact wedge: credit spread {path['d_spread'][0]:+.1f} bp/yr, "
          f"deposit rate {path['d_rdep'][0]:+.1f} bp/yr, "
          f"NET r_wc {path['d_r_wc'][0]:+.1f} bp/yr "
          f"({100*path['d_rdep'][0]/max(abs(path['d_spread'][0]), 1e-12):+.0f}% of the "
          f"spread is cancelled by the funding rate)")
    print(f"   like-for-like targets at this shock: {BOCOLA_IRF_OPEN:+.3f}% "
          f"(his open economy -- GHH + working capital, our own structure), "
          f"{BOCOLA_IRF_CLOSED:+.3f}% (his closed benchmark), "
          f"{BOCOLA_EPISODE_LEVEL:+.3f}% (2011Q4 episode)")
    _print_box_escapes(escapes, S0, S, T)


def _print_box_escapes(escapes, S0, S, T):
    # BOX ESCAPES ARE REPORTED, NOT SWALLOWED: a state pinned to a band turns a divergent
    # law of motion into a flat IRF that looks like convergence.
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


def dynamic_irf(rules, cal, ss, sproc, pd_shock=0.0198, T=25, rest_verbose=True):
    # DYNAMIC IRF: THE STATE VECTOR ITERATES FORWARD, IT IS NOT HELD AT THE SS.
    # Every endogenous state follows the period map's own law of motion while s decays
    # at rho_s -- the object Bocola's Table 5 describes. BOCOLA'S generate_irf.m, BOTH
    # HALVES: (1) start at the STOCHASTIC rest point, not the deterministic SS, and
    # (2) difference the shocked path against an UNSHOCKED path from the SAME state (his
    # `gdp = mean(gdp_s) - mean(gdp_nos)`), so the response cannot contain the no-shock
    # transition even if the rest point is imperfectly converged.
    S0 = stochastic_rest_point(rules, cal, ss, sproc, verbose=rest_verbose)
    S = S0.copy()
    escapes = []
    ref, K_ref = _no_shock_path(rules, cal, ss, sproc, S0, T)
    S[IS] = s_from_pd(pd_shock)
    print(f"\n  DYNAMIC IRF (states evolve; p^d shock to {100*pd_shock:.2f}%/qtr, "
          f"rho_s = {sproc['rho_s']})")
    print("   qtr  pd_q%  pd_a%    GDP%   GDP_ann%     C%      I%   hours%  "
          "rdep_bp  bank_bp  r_wc_bp  sov_bp    K%      n%")
    path = {k: [] for k in ("pd", "pd_ann", "Y", "Y_ann", "C", "I", "N", "spread",
                            "rdep", "r_wc", "K", "n", "Q_bD", "Q_bF", "sov_bp",
                            "d_rdep", "d_spread", "d_r_wc", "dQ_bD",
                            "mu", "E_Om", "Y_F", "C_F", "n_F", "spread_F", "I_F")}
    for t in range(T):
        o = read_at(rules, cal, ss, sproc, S)
        _append(path, _dynamic_row(o, ref[t], S, K_ref[t], cal, ss))
        if t in (0, 1, 2, 4, 6, 8, 12, 16, 20, 24):
            print(f"   {t:3d} {path['pd'][-1]:6.2f} {path['pd_ann'][-1]:6.2f} "
                  f"{path['Y'][-1]:+8.4f} {path['Y_ann'][-1]:+9.4f} "
                  f"{path['C'][-1]:+8.4f} {path['I'][-1]:+7.3f} {path['N'][-1]:+7.3f} "
                  f"{path['rdep'][-1]:8.1f} {path['spread'][-1]:8.1f} "
                  f"{path['r_wc'][-1]:8.1f} {path['sov_bp'][-1]:7.1f} "
                  f"{path['K'][-1]:+7.3f} {path['n'][-1]:+7.2f}")
        Sn = advance(o, S, sproc)
        S = rules.grid.clip(Sn)[0]
        esc = np.abs(Sn - S) > 1e-12
        if esc.any():
            escapes.append((t, [(_SNAMES[i], 100 * (Sn[i] / S0[i] - 1),
                                 100 * (S[i] / S0[i] - 1)) for i in np.flatnonzero(esc)]))
    _print_dynamic_summary(path, rules, cal, ss, sproc, S0, S, pd_shock, T, escapes)
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


def _tfp_row(o, base, z_t, cal):
    # ONE QUARTER OF THE TFP IRF, AS DEVIATIONS FROM THE SS READ.
    return {"Z": 100 * z_t,
            "Y": 100 * (o["Y_D"] / base["Y_D"] - 1),
            "C": 100 * (o["C_D"] / base["C_D"] - 1),
            "I": 100 * (o["I_D"] / base["I_D"] - 1),
            "I_F": 100 * (o["I_F"] / base["I_F"] - 1),
            "N": 100 * (o["_x"][0] / base["_x"][0] - 1),
            "Y_F": 100 * (o["Y_F"] / base["Y_F"] - 1),
            "C_F": 100 * (o["C_F"] / base["C_F"] - 1),
            "n": 100 * (o["n_D"] / base["n_D"] - 1),
            "n_F": 100 * (o["n_F"] / base["n_F"] - 1),
            "spread": _spread_bp(o, cal), "spread_F": _spread_bp(o, cal, "F"),
            "Q_bD": o["Q_bD"], "Q_bF": o["Q_bF"]}


def tfp_irf(rules, cal, ss, sproc, dz=0.01, T=21):
    # TFP IRF read off the no-default rules along the Z_D-decay path (rho_z from
    # sproc), endogenous states held at SS so the read stays on-grid -- the exact
    # image of persistence_irf, with the TFP state Z_D in place of the risk state s.
    S0 = ss_state(ss, cal, sproc)
    Z_ss = S0[IZ]
    base = _tfp_read(rules, cal, ss, sproc, S0.copy())
    print(f"\n  TFP IRF (one-off {dz:.0%} shock, rho_z={sproc['rho_z']} decay)")
    print("   qtr    Z%     Y_D%     C_D%     I_D%    hours%")
    # same paper series as dynamic_irf, so both figures read off one panel spec
    path = {k: [] for k in ("Z", "Y", "C", "I", "N", "Y_F", "C_F", "I_F", "n", "n_F",
                            "spread", "spread_F", "Q_bD", "Q_bF")}
    for t in range(T):
        z_t = dz * sproc["rho_z"] ** t
        S = S0.copy(); S[IZ] = Z_ss * np.exp(z_t)
        o = _tfp_read(rules, cal, ss, sproc, S)
        _append(path, _tfp_row(o, base, z_t, cal))
        if t in (0, 1, 2, 4, 6, 8, 12, 16, 20):
            print(f"   {t:3d}  {path['Z'][-1]:5.2f}  {path['Y'][-1]:+7.3f}  "
                  f"{path['C'][-1]:+7.3f}  {path['I'][-1]:+7.2f}  {path['N'][-1]:+6.2f}")
    return {k: np.array(v) for k, v in path.items()}


def main():
    cal = get_calibration()
    cal["nw_floor_frac"] = 0.15      # match main.py
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print("=== recursive global solution: pass-through of sovereign risk ===")
    rules = solve_recursive(cal, ss, sproc)
    impact_table(rules, cal, ss, sproc)
    persistence_irf(rules, cal, ss, sproc)


if __name__ == "__main__":
    main()
