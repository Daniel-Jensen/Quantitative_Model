# THE SOVEREIGN-RISK EXPERIMENT: THE SOLVE LADDER, THE READS AND THE IMPULSE RESPONSES.
import json
from pathlib import Path

import numpy as np
from scipy.optimize import root

from calibration.global_projection import get_calibration
from global_projection.steady_state import solve_steady_state
from global_projection.solver_recursive.state_grid import (build_state_box, build_aligned_box,
                                                           align_weights, _box_centre, ALIGN_BLOCK, BOX_STOCKS,
                                                           _STOCK_INDEX,
                                                           s_process_params, default_prob)
from global_projection.solver_recursive.state_grid import (IK_D, IK_F, IP_D, IP_F, IBDD,
                                          IBDF, IBFD, IV, IS, IZ, IM, IO, STATE_NAMES)
from global_projection.solver_recursive.decision_rules import RuleSet, SOLVE7, STORE_RULES, to_fit, from_fit
from global_projection.solver_recursive.collocation import solve_collocation, SETS, _point_block
from global_projection.solver_recursive.recursive_main import (time_iteration, calibrate_household_anchors,
                            ss_state, p_block_rotation)
from global_projection.reporting.prints import (bp_ann, ann_pct, ann_prob,
                              BOCOLA_IRF_CLOSED, BOCOLA_IRF_OPEN,
                              BOCOLA_EPISODE_LEVEL)
from global_projection.solver_recursive.point_map import point_residuals, _tpi_floor


# state names for the box-escape report
_SNAMES = STATE_NAMES

# the collocation box, shared by every experiment. Its bands are in balance-sheet units, so at
# leverage 8 (2026-10-20) the leverage-5 widths made the edge nodes swing net worth 2-3x more and
# the solve failed there; three bands are cut to their net-worth-equivalent widths: deposit
# obligations 4% -> 2.26%, the cross-border deposit position 4% -> 2%, the TPI book 100% -> 32% of
# the D bank's book (the book then fits ~26% of the D stock; check the TPI's footprint against it)
BOX_KW = dict(k_band=0.02, p_band_D=0.0226, p_band_F=0.0226, b_band=0.12, w_band=0.02,
              m_band=0.32)

# the two-box design (2026-10-21): a normal box aligned with the states the economy visits and a
# coarse post-default box with its own rules, blended into one continuation (box_specs.json, from
# Claude files/diagnostics/state_box/measure_box.py). SWITCHED OFF: on the coarse grid its
# post-default stage plateaued at 8e-3 and drove the normal default-regime rules off (STATE.md),
# so the model solves on the V5 box above, reading post-default states at its edge
TWO_BOXES = False
BOX_SPECS = json.loads(Path(__file__).with_name("box_specs.json").read_text())


def normal_box(ss, cal, mu=1, mu_vec=None, refine=None, at=None):
    # THE NORMAL-TIMES BOX: THE V5 BOX, OR THE MEASURED ALIGNED ONE WITH TWO_BOXES.
    if not TWO_BOXES:
        return build_state_box(ss, cal, mu=mu, mu_vec=mu_vec, refine=refine, at=at, **BOX_KW)
    assert at is None, "an instrument walk runs on the V5 box (TWO_BOXES off)"
    return build_aligned_box(ss, cal, BOX_SPECS["normal"], mu=mu, mu_vec=mu_vec, refine=refine)


def post_box(ss, cal):
    # THE POST-DEFAULT BOX, COARSE: ITS RULES REACH NORMAL TIMES ONLY THROUGH THE DEFAULT BRANCH.
    return build_aligned_box(ss, cal, BOX_SPECS["post"], mu=1)


def make_blend(ss, cal):
    # THE BLEND BETWEEN THE TWO RULE SETS AT THIS STEADY STATE (RuleSet.blend_weight).
    c, b = _box_centre(ss, cal), BOX_SPECS["blend"]
    return dict(block=list(ALIGN_BLOCK), centre=np.array([c["K_D"], c["P_D"], c["b_DD"], 0.0]),
                weights=align_weights(ss, cal), direction=np.asarray(b["direction"], dtype=float),
                post_until=float(b["post_until"]), normal_from=float(b["normal_from"]))


def _attach_post(rules, ss, cal):
    # THE POST-DEFAULT RULE SET, SEEDED WITH THE NORMAL RULES READ THROUGH THEIR BOX EDGE (WHAT THE
    # CONTINUATION USED BEFORE IT EXISTED), ATTACHED WITH ITS BLEND.
    post = RuleSet(post_box(ss, cal), rules.n_regimes)
    post.n_gh = rules.n_gh
    for d in range(rules.n_regimes):
        seed = rules.eval_cont(d, post.grid.points)
        for k in STORE_RULES:
            post.set_values(k, d, seed[k])
    rules.attach_post(post, make_blend(ss, cal))


# quadrature order; every reader reuses rules.n_gh
N_GH = 5

# dense Chebyshev nodes in s: 5 is the production setting (converged for every reported object on
# the 2026-08 tensor grid), 9 is Bocola's
S_REFINE = 5
# how the s nodes are added: "axis" (on the s axis alone, 27 points at 5 nodes), "tensor" (times
# the whole sparse grid, 115 points) or "cross" (2026-10-23): the axis grid on a box whose bank stocks
# are re-centred on the rest point (_recentre_box), then the s x state cross terms (75 points at 9
# nodes). At leverage 8 the tensor grid's s x bank-state nodes have no collocation solution near the
# coarse one (2026-10-21); at p^d 0.02% the axis grid on the SS box mispriced the D bond by 0.036 at
# rest, the re-centred one by 0.011 and the cross grid by 0.0009 (Claude files/docs/STATE.md)
S_REFINE_MODE = "cross"
# Newton steps per refined stage in the cross mode: the first axis pass on the SS box stops at ~2e-5
# and is only used to locate the rest point
CROSS_MAXIT = 60

# time-iteration sweeps that put each Newton inside its basin
WARM_SWEEPS = 12

# the TPI cap is walked down from 700 bp in 50 bp rungs, halving a rung that fails
TPI_CAP_START = 700.0
TPI_CAP_STEP = 50.0
TPI_MIN_STEP = 6.25
# the TPI smoothing is eased in on the refined grid, ending at cal["tpi_eps"]
TPI_EPS_LADDER = (1e-2, 1e-3)


def _seed_from(rules_fine, rules_coarse):
    # EVALUATE A SOLVED COARSE RULE SET AT A FINER GRID'S POINTS.
    pts = rules_fine.grid.points
    for k in STORE_RULES:
        for d in range(rules_coarse.n_regimes):
            rules_fine.set_values(k, d, rules_coarse.eval(k, d, pts))
    rules_fine.n_gh = rules_coarse.n_gh
    return rules_fine


def _stage(rules, cal, ss, sproc, regimes, no_default, label, verbose,
           backend="auto", warm=WARM_SWEEPS, maxit=40, sets=SETS):
    # ONE SOLVE STAGE: A SHORT WARM START, THEN THE GLOBAL NEWTON (sets: WHICH RULE SETS SOLVE).
    if warm:
        time_iteration(rules, cal, ss, sproc, regimes=regimes, no_default=no_default,
                       damp=0.5, tol=1e-4, max_it=warm, n_gh=N_GH, verbose=False)
    return solve_collocation(rules, cal, ss, sproc, regimes=regimes,
                             no_default=no_default, n_gh=N_GH, backend=backend,
                             maxit=maxit, verbose=verbose, label=f" {label}", sets=sets)


def liquidity_ceiling_report(rules, cal, ss, sproc, target=None, regime=0,
                             verbose=True):
    # THE LIQUIDITY CEILING: THE BOND PRICE WITH mu = 0, CONTINUATION HELD FIXED.
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


def _build_rules(cal, ss, sproc, mu, mu_vec, rotate, nreg, verbose, at=None):
    # THE COARSE GRID (OPTIONALLY P-ROTATED) AND A STEADY-STATE COLD START.
    rot, centre = (None, None)
    box_kw = dict(BOX_KW)
    if rotate:
        rot, centre, J, evs = p_block_rotation(ss, cal, sproc, mu=mu, mu_vec=mu_vec,
                                               probe_kw=box_kw)
        if verbose:
            print(f"  P-block rotation: |lambda| = {evs[0]:.3f}, {evs[1]:.3f}"
                  f"   rho(|J|) = {np.max(np.abs(np.linalg.eigvals(np.abs(J)))):.3f}")
    grid = (build_state_box(ss, cal, mu=mu, mu_vec=mu_vec, rot=rot, centre=centre, at=at, **box_kw)
            if rotate else normal_box(ss, cal, mu=mu, mu_vec=mu_vec, at=at))
    rules = RuleSet.from_ss(grid, ss, cal, n_regimes=nreg)
    rules.n_gh = N_GH
    if verbose:
        print(f"  coarse grid: mu={mu}, {grid.n} points x {nreg} regimes, n_gh={N_GH}")
    # the box spec travels with the rules, for the refinement and for a reused baseline (an
    # instrument walk holds its first economy's bank-stock centre, at)
    rules.box = dict(mu=mu, mu_vec=mu_vec, rot=rot, centre=centre, box_kw=box_kw,
                     aligned=not rotate, at=at)
    return rules, rules.box


WC_LADDER = (0.0, 0.2, 0.5, 1.0)
# the default regime's recovery, walked down from a mild haircut to the target in 5 pp rungs: at
# leverage 8 the D bank turns insolvent in default past ~0.6 and the funded floor binds at most
# nodes, and from 15 pp rungs the Newton started outside its basin (2026-10-20)
HAIRCUT_LADDER = (0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50)
# the cross-border p'/p valuation, walked in on the first rung: from the cold start the Newton
# stalls on the valued economy (max|F| ~4e-3) but finds the unvalued one, and from that it
# takes ~3 steps per rung (measured on the old calibration, 2026-10-19)
FX_LADDER = (0.0, 0.5, 1.0)
# the KKT smoothing ladder (2026-10-20): every solve starts at MU_EPS_START and walks
# cal["mu_eps"] down to its calibrated value. At leverage 8 no cold start converges at the
# calibrated smoothing (the joint stage lands on a near-fold, not a root); the walk follows
# the regular branch. The log-step starts at MU_EPS_STEP, grows x1.5 on success, halves on a
# failure, and the walk stops below MU_EPS_MIN_STEP
MU_EPS_START = 1e-2
MU_EPS_STEP = 1.4
MU_EPS_MIN_STEP = 1.02
# the instrument walk (2026-10-22): the SS cold start solves near these instruments (the targets'
# own values), not at the rest-point calibration's: above a D SS-spread instrument of ~125 bp its
# d = 0 rungs stall (at 200 bp the joint stage stops at max|F| 1.6e-4), and every 10 bp moves the
# SS ~0.8% in capital. solve_recursive solves this economy cold, then walks the instruments to the
# requested ones on its box, one joint solve per rung (_walk_instruments)
WALK_FROM = dict(leverage_target_D=8.0, leverage_target_F=8.0,
                 credit_spread_target_D=0.0025, credit_spread_target_F=0.0025,
                 sov_exposure_D=0.144, sov_exposure_F=0.065)
# the largest move per rung (10 bp/yr of spread converges in 8-23 Newton steps)
WALK_STEP = dict(leverage_target_D=0.1, leverage_target_F=0.1,
                 credit_spread_target_D=10 / 4e4, credit_spread_target_F=10 / 4e4,
                 sov_exposure_D=0.005, sov_exposure_F=0.005)


def _restore(rules, good):
    # PUT A SAVED RULE SET'S VALUES BACK INTO rules (BOTH BOXES' SETS).
    for rs, gs in ((rules, good), (rules.post, good.post)):
        if rs is None:
            continue
        for k in STORE_RULES:
            for d in range(rs.n_regimes):
                rs.set_values(k, d, gs.vals[k][d])


def _walk_mu_eps(rules, cal, ss, sproc, regimes, no_default, target, label, verbose,
                 backend="auto"):
    # WALK THE KKT SMOOTHING DOWN TO target, THEN RE-ANCHOR THE HOUSEHOLDS AND POLISH THERE.
    assert target > 0.0, "the smoothing walk needs cal['mu_eps'] > 0 (a log walk never reaches 0)"
    eps, lstep = cal["mu_eps"], np.log(MU_EPS_STEP)
    while eps > target * (1.0 + 1e-12):
        nxt, good = max(target, eps / np.exp(lstep)), rules.copy()
        cal["mu_eps"] = nxt
        ok, _, _ = _stage(rules, cal, ss, sproc, regimes, no_default,
                          f"{label} mu_eps={nxt:.2e}", verbose, backend=backend, warm=0, maxit=20)
        if ok:
            eps, lstep = nxt, min(1.5 * lstep, np.log(3.0))
            continue
        _restore(rules, good)
        cal["mu_eps"], lstep = eps, 0.5 * lstep
        if lstep < np.log(MU_EPS_MIN_STEP):
            break
    # the anchors absorb the SS household budget, which moves with the smoothing (via r_wc)
    calibrate_household_anchors(cal, ss, sproc)
    ok, _, w = _stage(rules, cal, ss, sproc, regimes, no_default, f"{label} mu_eps={eps:.2e}",
                      verbose, backend=backend, warm=0, maxit=20)
    return bool(ok and eps <= target * (1.0 + 1e-12)), w


# the working smoothings a genuine (zero) target is solved at, before the final hard polish
# (2026-10-26: the IC is a genuine max, the book pairs and the funded floor are hard too)
HARD_WORK = dict(mu_eps=1e-4, book_eps=1e-4, floor_eps=1e-3)


def _soften(cal):
    # SET EVERY ZERO (GENUINE) SMOOTHING TARGET TO ITS WORKING VALUE; RETURNS THE TARGETS.
    targets = {k: float(cal.setdefault(k, v)) for k, v in HARD_WORK.items()}
    for k, v in HARD_WORK.items():
        if targets[k] <= 0.0:
            cal[k] = v
    return targets


def _hard_polish(rules, cal, ss, sproc, targets, regimes, no_default, label, verbose,
                 backend="auto"):
    # THE FINAL POLISH AT THE GENUINE CONDITIONS: ALL AT ONCE, ELSE HARDENED ONE BY ONE IN
    # targets' ORDER; RETURNS (ok, max|F|) AND LEAVES cal AT WHAT THE RULES SOLVE.
    todo = [k for k, v in targets.items() if cal[k] != v]
    if not todo:
        return True, np.nan
    good, work = rules.copy(), {k: cal[k] for k in todo}
    for k in todo:
        cal[k] = targets[k]
    calibrate_household_anchors(cal, ss, sproc)
    ok, _, w = _stage(rules, cal, ss, sproc, regimes, no_default, f"{label} hard polish",
                      verbose, backend=backend, warm=0, maxit=20)
    if ok:
        return True, w
    # one at a time from the smoothed solution: the first that fails stays at its working value
    _restore(rules, good)
    cal.update(work)
    for k in todo:
        good = rules.copy()
        cal[k] = targets[k]
        calibrate_household_anchors(cal, ss, sproc)
        ok, _, w = _stage(rules, cal, ss, sproc, regimes, no_default, f"{label} hard {k}",
                          verbose, backend=backend, warm=0, maxit=20)
        if not ok:
            print(f"    NOTE: the hard polish stalls on {k} (max|F| {w:.1e}); "
                  f"the rules solve {k} = {work[k]:.1e}")
            _restore(rules, good)
            cal[k] = work[k]
            calibrate_household_anchors(cal, ss, sproc)
            return False, w
    return True, w


def _economy(cal, over):
    # cal WITH THE INSTRUMENTS IN over, ITS SS, ITS s/Z PROCESS (z_star IS THE SS's Z) AND ANCHORS.
    cal_k = dict(cal, **over)
    ss_k = solve_steady_state(cal_k, verbose=False)
    sp_k = s_process_params(cal_k)
    calibrate_household_anchors(cal_k, ss_k, sp_k)
    return cal_k, ss_k, sp_k


def _walk_start(cal):
    # THE ECONOMY THE COLD START SOLVES FOR cal (cal AT THE WALK_FROM INSTRUMENTS); NONE WHEN cal
    # IS ALREADY THERE.
    if all(abs(cal[k] - v) <= 1e-12 for k, v in WALK_FROM.items()):
        return None
    return _economy(cal, WALK_FROM)


def _walk_instruments(rules, cal0, cal, ss, sproc, verbose, backend="auto"):
    # WALK THE INSTRUMENTS FROM cal0'S TO cal'S, ONE JOINT SOLVE PER RUNG, HALVING A RUNG THAT
    # FAILS; THE LAST RUNG IS (cal, ss, sproc) ITSELF. RETURNS (rules, ok, max|F|).
    # each rung's box holds the bank stocks where rules.box["at"] puts them and centres s and Z_D
    # on the rung's own economy (its SS moves Z_D), so the rules are re-seeded on its nodes
    box = rules.box
    gap = {k: cal[k] - cal0[k] for k in WALK_FROM}
    n = max(1, int(np.ceil(max(abs(g) / WALK_STEP[k] for k, g in gap.items()) - 1e-9)))
    t, h, ok, w = 0.0, 1.0 / n, False, np.nan
    while t < 1.0 - 1e-12:
        nxt = min(1.0, t + h)
        if nxt >= 1.0 - 1e-12:
            cal_k, ss_k, sp_k = cal, ss, sproc
        else:
            cal_k, ss_k, sp_k = _economy(cal, {k: cal0[k] + nxt * g for k, g in gap.items()})
        grid = normal_box(ss_k, cal_k, mu=box["mu"], mu_vec=box["mu_vec"], at=box["at"])
        trial = _seed_from(RuleSet(grid, rules.n_regimes), rules)
        trial.box = box
        ok, _, w = _stage(trial, cal_k, ss_k, sp_k, (0, 1), False, f"walk {nxt:.3f}", verbose,
                          backend=backend, warm=0)
        if ok:
            rules, t, h = trial, nxt, min(1.0 / n, 2.0 * h)
            continue
        h *= 0.5
        if h < 1.0 / (16 * n):
            break
    return rules, bool(ok and t >= 1.0 - 1e-12), w


def _stage_wc(rules, cal, ss, sproc, regimes, no_default, label, verbose):
    # THE FIRST d = 0 SOLVE, WALKING THE FX VALUATION, THEN THE WORKING-CAPITAL FLOW, INTO PLACE.
    # From the steady-state seed the Newton cannot find the corrected economy directly; from
    # the solution with the flow held at its SS value it takes 3-5 steps per rung.
    first = True
    for k, w in enumerate(WC_LADDER):
        cal["wc_flow_D"] = cal["wc_flow_F"] = w
        calibrate_household_anchors(cal, ss, sproc)  # the anchor absorbs the SS flow at w < 1
        for fx in (FX_LADDER if k == 0 else (1.0,)):
            cal["fx_val"] = fx
            tag = f"{label} wc={w:.1f}" + (f" fx={fx:.1f}" if k == 0 else "")
            ok, it, wst = _stage(rules, cal, ss, sproc, regimes, no_default, tag, verbose,
                                 **({} if first else {"warm": 0}))
            first = False
    for key in ("fx_val", "wc_flow_D", "wc_flow_F"):
        cal.pop(key, None)
    calibrate_household_anchors(cal, ss, sproc)
    return ok, it, wst


def solve_no_default(rules, cal, ss, sproc, label, verbose):
    # THE NO-DEFAULT BLOCK AT pi = 0 (THE TFP RULES, NESTING GATE N2): THE FIRST STAGE AT
    # MU_EPS_START, THEN THE SMOOTHING WALKED DOWN; RETURNS (ok, max|F|).
    hard = _soften(cal)
    target = cal["mu_eps"]
    cal["mu_eps"] = max(target, MU_EPS_START)
    _stage_wc(rules, cal, ss, sproc, (0,), True, label, verbose)
    ok, w = _walk_mu_eps(rules, cal, ss, sproc, (0,), True, target, label, verbose)
    okh, wh = _hard_polish(rules, cal, ss, sproc, hard, (0,), True, label, verbose)
    return ok and okh, (w if np.isnan(wh) else wh)


# a normal-times node on the self-fulfilling crisis branch (2026-10-28, appendix B: the K_F+ node at n_D 11% of
# n_ss and alpha_D 8, reading its own continuation with weight 0.97; it stalled every stage): the signature is D
# equity below CRISIS_N of its SS level or alpha_D above CRISIS_ALPHA, at a node whose priced default probability is at
# most CRISIS_PD_MAX: above it low equity is the genuine distressed solution, and the centre-seeded root there is the
# spurious one (2026-10-29: at p^d 12% it gave n_D 138% of n_ss and cascaded through the s x state terms)
CRISIS_N, CRISIS_ALPHA, CRISIS_PD_MAX = 0.5, 3.0, 0.01
# node-by-node polish sweeps after a switch, so the nodes agree before the global Newton, each node's move
# capped at CALM_JUMP in fit space (2026-10-29: uncapped sweeps sent converged high-risk nodes to another root)
CALM_SWEEPS = 4
CALM_JUMP = 0.05


def _node_values(rules, i, regimes):
    # NODE i's VALUES OF EVERY RULE IN FIT SPACE, RULE BY RULE THEN REGIME.
    return np.array([to_fit(k, rules.vals[k][d][i]) for k in STORE_RULES for d in regimes], dtype=float)


def _node_residual(rules, cal, ss, sproc, i, regimes, no_default, z):
    # NODE i's RESIDUALS WITH ITS VALUES SET TO z AND REFITTED, THE OTHER NODES HELD.
    j = 0
    for k in STORE_RULES:
        for d in regimes:
            v = rules.vals[k][d].copy()
            v[i] = from_fit(k, z[j])
            rules.set_values(k, d, v)
            j += 1
    vals = {k: {d: rules.vals[k][d] for d in regimes} for k in STORE_RULES}
    return np.concatenate([_point_block(i, d, vals, rules, cal, ss, sproc, N_GH, no_default) for d in regimes])


def _node_solve(rules, cal, ss, sproc, i, regimes, no_default, z0):
    # NODE i SOLVED ALONE FROM z0; THE RULES KEEP THE ROOT. RETURNS ITS max|residual|.
    # a non-finite residual counts as far from a root, so the caller keeps the node's values
    f = lambda z: np.nan_to_num(_node_residual(rules, cal, ss, sproc, i, regimes, no_default, z),
                                nan=1e3, posinf=1e3, neginf=-1e3)
    try:
        sol = root(f, z0, method="hybr", tol=1e-13)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError):
        return np.inf
    return float(np.max(np.abs(f(sol.x))))


def _node_nD(rules, cal, ss, sproc, i, no_default):
    # THE D BANK'S NET WORTH AT NODE i IN NORMAL TIMES, AT THE RULES' CURRENT VALUES.
    x = np.array([rules.vals[k][0][i] for k in SOLVE7])
    return point_residuals(rules.grid.points[i], 0, x, rules, cal, ss, sproc, n_gh=N_GH,
                           no_default=no_default)[1]["n_D"]


def _crisis_nodes(rules, cal, ss, sproc, no_default):
    # THE NODES WHOSE NORMAL-TIMES SOLUTION CARRIES THE CRISIS SIGNATURE.
    nss, out = ss["ss_bank_D"]["n_ss"], []
    for i in range(rules.grid.n):
        if default_prob(rules.grid.points[i, IS]) > CRISIS_PD_MAX:
            continue
        x = np.array([rules.vals[k][0][i] for k in SOLVE7])
        o = point_residuals(rules.grid.points[i], 0, x, rules, cal, ss, sproc, n_gh=N_GH,
                            no_default=no_default)[1]
        if o["n_D"] < CRISIS_N * nss or o["alpha_D"] > CRISIS_ALPHA:
            out.append(i)
    return out


def _calm_branch(rules, cal, ss, sproc, regimes, no_default, label, verbose, backend="auto"):
    # EVERY CRISIS NODE SWITCHED TO ITS CALM ROOT (SOLVED ALONE FROM THE CENTRE NODE'S VALUES, CAPITAL SCALED TO
    # ITS OWN STATE), THE NODES POLISHED IN SWEEPS, THE STAGE RE-SOLVED. RETURNS (ok, max|F|), OR None IF NONE.
    bad = _crisis_nodes(rules, cal, ss, sproc, no_default)
    if not bad:
        return None
    g = rules.grid
    c = int(np.argmin(np.abs(g.to_unit(g.points)).sum(axis=1)))
    kp = {"Kp_D": IK_D, "Kp_F": IK_F}
    switched = []
    for i in bad:
        z_now = _node_values(rules, i, regimes)
        z0 = _node_values(rules, c, regimes)
        j = 0
        for k in STORE_RULES:
            for _ in regimes:
                if k in kp:
                    z0[j] += np.log(g.points[i, kp[k]] / g.points[c, kp[k]])
                j += 1
        n_now = _node_nD(rules, cal, ss, sproc, i, no_default)
        ok = _node_solve(rules, cal, ss, sproc, i, regimes, no_default, z0) <= 1e-10
        n_new = _node_nD(rules, cal, ss, sproc, i, no_default)
        # a flagged node without a calmer converged root is a genuine distressed state (high p^d)
        if ok and n_new > n_now * (1.0 + 1e-3):
            switched.append(i)
        else:
            _node_residual(rules, cal, ss, sproc, i, regimes, no_default, z_now)
    if not switched:
        return None
    for _ in range(CALM_SWEEPS):
        for i in range(g.n):
            z_now = _node_values(rules, i, regimes)
            f_now = float(np.max(np.abs(_node_residual(rules, cal, ss, sproc, i, regimes, no_default, z_now))))
            f_new = _node_solve(rules, cal, ss, sproc, i, regimes, no_default, z_now)
            if (not np.isfinite(f_new) or f_new >= f_now
                    or np.max(np.abs(_node_values(rules, i, regimes) - z_now)) > CALM_JUMP):
                _node_residual(rules, cal, ss, sproc, i, regimes, no_default, z_now)
    if verbose:
        print(f"    {label}: crisis branch at nodes {switched} (flagged {bad}); calm roots selected "
              f"(left: {_crisis_nodes(rules, cal, ss, sproc, no_default) or 'none'})")
    ok, _, w = _stage(rules, cal, ss, sproc, regimes, no_default, f"{label} (calm)", verbose,
                      backend=backend, warm=0)
    return ok, w


def _solve_baseline(rules, cal, ss, sproc, D_REG, backend, verbose):
    # THE NO-TPI BASELINE: d = 0, THEN THE HAIRCUT HOMOTOPY, THEN THE JOINT SOLVE, ALL AT
    # MU_EPS_START, THEN THE SMOOTHING WALKED DOWN TO ITS CALIBRATED VALUE.
    target = cal["mu_eps"]
    cal["mu_eps"] = max(target, MU_EPS_START)
    ok0, it0, w0 = _stage_wc(rules, cal, ss, sproc, (0,), True, "d0", verbose)
    calm = _calm_branch(rules, cal, ss, sproc, (0,), True, "d0", verbose, backend)
    if calm is not None:
        ok0, w0 = calm
    for k in STORE_RULES:  # start the default regime from d = 0
        rules.set_values(k, D_REG, rules.vals[k][0].copy())
    rec_target = cal["recovery_rate_D"]
    ok1, w1 = False, np.nan
    for rec in HAIRCUT_LADDER + (rec_target,):
        cal["recovery_rate_D"] = rec
        ok1, it1, w1 = _stage(rules, cal, ss, sproc, (D_REG,), False,
                              f"d1 rec={rec:.2f}", verbose)
    cal["recovery_rate_D"] = rec_target
    # straight from the default-regime solve, no time-iteration sweeps: at leverage 8 they can
    # carry the start onto a second configuration at the capital edges, one that folds before
    # the calibrated smoothing (2026-10-20)
    okj, itj, wj = _stage(rules, cal, ss, sproc, (0, D_REG), False,
                          "joint (no TPI)", verbose, backend=backend, warm=0)
    calm = _calm_branch(rules, cal, ss, sproc, (0, D_REG), False, "joint (no TPI)", verbose, backend)
    if calm is not None:
        okj, wj = calm
    if TWO_BOXES:
        # the post-default rules: seeded from the edge reads the joint stage just used, solved
        # with the normal rules held, then both boxes together
        _attach_post(rules, ss, cal)
        _stage(rules, cal, ss, sproc, (0, D_REG), False, "post-default box", verbose,
               backend=backend, warm=0, sets=("post",))
        okj, _, wj = _stage(rules, cal, ss, sproc, (0, D_REG), False, "joint (both boxes)",
                            verbose, backend=backend, warm=0)
    okw, ww = _walk_mu_eps(rules, cal, ss, sproc, (0, D_REG), False, target,
                           "joint (no TPI)", verbose, backend=backend)
    return (ok0, ok1, okj and okw), (w0, w1, ww)


def _tpi_form(rules, cal, ss, gz):
    # SWITCH THE TPI VARIABLE BETWEEN THE FB (SOLVE) AND GZ (READ) FORMS AT THE NODES (BOTH SETS).
    for rs in (rules, rules.post):
        if rs is None:
            continue
        x = rs.vals["x_cb"][0]
        if gz:
            gap = (rs.vals["Q_bD"][0] - _tpi_floor(rs.vals["Q_bF"][0], cal)) / ss["Q_bD_ss"]
            x = np.where(x > gap, x, -np.maximum(gap, 0.0))
        else:
            x = np.maximum(x, 0.0)
        rs.set_values("x_cb", 0, x)
    cal["tpi_gz"] = bool(gz)


def _tpi_polish(rules, cal, ss, sproc, label, verbose, backend):
    # CONVERT AN FB SOLUTION TO THE GZ FORM AND RE-ROOT IT (1-2 STEPS).
    _tpi_form(rules, cal, ss, gz=True)
    return _stage(rules, cal, ss, sproc, tuple(range(rules.n_regimes)), False,
                  f"{label} (GZ form)", verbose, backend=backend, warm=0, maxit=6)


def _solve_tpi(rules, cal, ss, sproc, backend, verbose):
    # WALK THE TPI CAP DOWN TO ITS TARGET ON A SOLVED BASELINE.
    target = cal["tpi_cap_bp"]
    regimes = tuple(range(rules.n_regimes))
    cap, last, step = max(TPI_CAP_START, target), None, TPI_CAP_STEP
    _tpi_form(rules, cal, ss, gz=False)  # the walk runs in the FB form
    good, ok, w = rules.copy(), False, np.nan
    while True:
        cal["tpi_cap_bp"] = cap
        ok, _, w = _stage(rules, cal, ss, sproc, regimes, False, f"TPI cap={cap:.1f}bp",
                          verbose, backend=backend, warm=0)
        if ok:
            good, last = rules.copy(), cap
            if cap <= target:
                break
            cap = max(cap - step, target)
            continue
        step *= 0.5
        if last is None or step < TPI_MIN_STEP:
            break
        rules = good.copy()
        cap = max(last - step, target)
    if not ok:
        rules = good
    if last is not None:  # hand back in the form the reads use
        cal["tpi_cap_bp"] = last
        ok_gz, _, w = _tpi_polish(rules, cal, ss, sproc, f"TPI cap={last:.1f}bp", verbose,
                                  backend)
        ok = ok and ok_gz
    else:
        _tpi_form(rules, cal, ss, gz=True)
    cal["tpi_cap_bp"] = target
    rules.tpi_cap_solved = last
    return rules, bool(ok and last is not None and last <= target), w


def _refine_ladder(s_refine, tpi, cross=True):
    # THE (NODES IN s, GRID MODE) RUNGS UP TO s_refine; cross=False STOPS THE CROSS MODE AT ITS AXIS RUNGS.
    ms = [m for m in (5, 9, 17) if 1 < m < s_refine] + [s_refine]
    if S_REFINE_MODE == "cross":
        return [(m, "axis") for m in ms] + ([(s_refine, "cross")] if cross else [])
    # with the TPI on, start from 3 nodes in s, which seeds the 5-node solve well (on the axis grid
    # 3 nodes is the coarse grid itself)
    if tpi and s_refine > 3 and S_REFINE_MODE == "tensor":
        ms = [3] + ms
    return [(m, S_REFINE_MODE) for m in ms]


def _refine_s(rules, cal, ss, sproc, box, s_refine, backend, verbose, ladder=None):
    # REFINE THE s DIMENSION STEP BY STEP, EACH STEP SEEDED FROM THE LAST.
    nreg = rules.n_regimes
    okj = wj = None
    tpi = bool(cal["tpi_on"])
    maxit = CROSS_MAXIT if S_REFINE_MODE == "cross" else 20
    for m_s, mode in (_refine_ladder(s_refine, tpi) if ladder is None else ladder):
        ref = (IS, m_s, mode)
        gfine = (normal_box(ss, cal, mu=box["mu"], mu_vec=box["mu_vec"], refine=ref,
                            at=box.get("at"))
                 if box.get("aligned") else
                 build_state_box(ss, cal, mu=box["mu"], mu_vec=box["mu_vec"], rot=box["rot"],
                                 centre=box["centre"], refine=ref, at=box.get("at"),
                                 **box["box_kw"]))
        if verbose:
            print(f"  s-refined grid: {gfine.n} points x {nreg} regimes "
                  f"({m_s} nodes {dict(axis='on the s axis', cross='in s, with s x state terms').get(mode, 'in s')}, "
                  f"degree {m_s - 1} in s)")
        fine = _seed_from(RuleSet(gfine, nreg), rules)
        fine.box = rules.box
        if rules.post is not None:  # the post-default set stays coarse (2026-10-21)
            fine.attach_post(rules.post.copy(), rules.blend)
        # straight from the seed: time-iteration sweeps damaged it at leverage 8 (2026-10-20)
        if not tpi:
            okj, itj, wj = _stage(fine, cal, ss, sproc, tuple(range(nreg)), False,
                                  f"joint (s={m_s} {mode})", verbose, backend=backend,
                                  warm=0, maxit=maxit)
        else:  # root in the FB form, easing the smoothing in, then hand back in GZ
            _tpi_form(fine, cal, ss, gz=False)
            eps, warm = cal["tpi_eps"], 0
            for e in [e for e in TPI_EPS_LADDER if e > eps] + [eps]:
                cal["tpi_eps"] = e
                okj, itj, wj = _stage(fine, cal, ss, sproc, tuple(range(nreg)), False,
                                      f"joint (s={m_s} {mode}, eps={e:.0e})", verbose,
                                      backend=backend, warm=warm, maxit=maxit)
                warm = 0
            cal["tpi_eps"] = eps
            ok_gz, _, wj = _tpi_polish(fine, cal, ss, sproc, f"joint (s={m_s})", verbose,
                                       backend)
            okj = okj and ok_gz
        rules = fine
    return rules, okj, wj


def _recentre_box(rules, cal, ss, sproc, s_refine, backend, verbose):
    # THE CROSS MODE'S FIRST PASS: THE AXIS REFINEMENT ON THE SS BOX, ITS REST POINT, AND THE COARSE GRID
    # SOLVED ON A BOX WHOSE BANK STOCKS ARE CENTRED THERE. RETURNS (coarse rules, ok, max|F|).
    # the mu = 1 rules are additive, accurate only on the axes through the box centre, and the long bond
    # sums their error over its duration; the first pass need not converge, it only locates the rest point
    first, _, w1 = _refine_s(rules, cal, ss, sproc, rules.box, s_refine, backend, verbose,
                             ladder=_refine_ladder(s_refine, False, cross=False))
    S0 = stochastic_rest_point(first, cal, ss, sproc, verbose=verbose)
    at = {k: float(S0[_STOCK_INDEX[k]]) for k in BOX_STOCKS}
    box = dict(rules.box, at=at)
    coarse = _seed_from(RuleSet(normal_box(ss, cal, mu=box["mu"], mu_vec=box["mu_vec"], at=at),
                                rules.n_regimes), first)
    coarse.box = box
    if verbose:
        print(f"  box re-centred on that rest point (first pass max|F| {w1:.1e}): " +
              ", ".join(f"{k} {100 * (v / _box_centre(ss, cal)[k] - 1):+.2f}%" for k, v in at.items()))
    ok, _, w = _stage(coarse, cal, ss, sproc, tuple(range(rules.n_regimes)), False,
                      "joint (re-centred box)", verbose, backend=backend, warm=0, maxit=40)
    return coarse, ok, w


def _stamp_verdict(rules, oks, worsts):
    # RECORD ON THE RULES WHETHER EVERY STAGE CONVERGED.
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
    # SOLVE EVERY REGIME BY GLOBAL COLLOCATION (BOCOLA'S LADDER).
    nreg = 2
    D_REG = 1  # the regime index is the default indicator
    tpi = bool(cal["tpi_on"])
    # every stage runs at the working smoothings; the genuine conditions are polished in at the end
    hard = _soften(cal)

    if base is None:
        # the cold start solves the WALK_FROM economy on its own box; the walk then brings the
        # instruments to cal's on that box
        start = _walk_start(cal)
        cal0, ss0, sp0 = (cal, ss, sproc) if start is None else start
        rules, box = _build_rules(cal0, ss0, sp0, mu, mu_vec, rotate, nreg, verbose)
        cal["tpi_on"] = cal0["tpi_on"] = False  # the baseline is always the no-TPI economy
        (ok0, ok1, okb), (w0, w1, wb) = _solve_baseline(rules, cal0, ss0, sp0, D_REG,
                                                        backend, verbose)
        if start is not None:
            # the walk holds the base economy's bank-stock centre; s and Z_D follow each rung
            c0 = _box_centre(ss0, cal0)
            rules.box = box = dict(box, at={k: c0[k] for k in BOX_STOCKS})
            rules, okw, ww = _walk_instruments(rules, cal0, cal, ss, sproc, verbose, backend)
            okb, wb = okb and okw, float(np.nanmax([wb, ww]))
        cal["tpi_on"] = tpi
        if base_out is not None:  # a reusable snapshot for the caller
            base_out.append(rules.copy())
    else:
        rules = base.copy()
        box = rules.box
        if box is None:  # a baseline pickled before the box travelled with the rules
            box = _build_rules(cal, ss, sproc, mu, mu_vec, rotate, nreg, False)[1]
            rules.box = box
        assert rules.n_regimes == nreg, "reused baseline must carry the same regime count"
        ok0 = ok1 = okb = True
        w0 = w1 = wb = np.nan
        if verbose:
            print("  reusing the solved baseline")

    cap_solved, target = None, cal["tpi_cap_bp"]
    if tpi:  # with the TPI on, the cap walk replaces the joint polish
        rules, okj, wj = _solve_tpi(rules, cal, ss, sproc, backend, verbose)
        cap_solved = rules.tpi_cap_solved
    else:
        # a polish of the walked-down baseline: no time-iteration sweeps, which can leave the
        # regular branch the smoothing walk followed
        okj, itj, wj = _stage(rules, cal, ss, sproc, tuple(range(nreg)), False, "joint",
                              verbose, backend=backend, warm=0)
        calm = _calm_branch(rules, cal, ss, sproc, tuple(range(nreg)), False, "joint", verbose, backend)
        if calm is not None:
            okj, wj = calm

    if s_refine and s_refine > 1 and S_REFINE_MODE == "cross" and base is None and not tpi:
        # the cross mode re-centres the box first; the TPI reuses the re-centred coarse baseline
        rules, ok_rc, w_rc = _recentre_box(rules, cal, ss, sproc, s_refine, backend, verbose)
        box = rules.box
        okb, wb = okb and ok_rc, float(np.nanmax([wb, w_rc]))
        if base_out is not None:
            base_out[:] = [rules.copy()]

    if s_refine and s_refine > 1:
        # a cap walk that stopped short refines at the cap it reached
        if tpi and cap_solved is not None:
            cal["tpi_cap_bp"] = cap_solved
        rules, okj2, wj = _refine_s(rules, cal, ss, sproc, box, s_refine, backend, verbose)
        okj = okj2 and (okj or not tpi)  # a cap walk short of its target counts as a failure
        cal["tpi_cap_bp"] = target
    rules.tpi_cap_solved = cap_solved
    okh, wh = _hard_polish(rules, cal, ss, sproc, hard, tuple(range(nreg)), False, "joint",
                           verbose, backend=backend)
    okj, wj = okj and okh, (wj if np.isnan(wh) else wh)
    # the smoothings the returned rules solve (all zero when the hard polish converged)
    rules.smoothing = {k: cal[k] for k in hard}

    _stamp_verdict(rules, (ok0, ok1, okb, okj), (w0, w1, wb, wj))
    return rules


def read_at(rules, cal, ss, sproc, S):
    # READ THE FITTED RULES AT STATE S.
    Sm = np.atleast_2d(S)
    x = np.array([float(rules.eval_read(k, 0, Sm)[0]) for k in SOLVE7])
    res, o = point_residuals(S, 0, x, rules, cal, ss, sproc,
                             n_gh=rules.n_gh or N_GH, no_default=False)
    o["_x"] = x
    o["_resid"] = float(np.max(np.abs(res)))
    return o


def read_exact(rules, cal, ss, sproc, S, x0=None):
    # CLEAR THE PERIOD MAP EXACTLY AT S; NEAR THE mu KINK THIS DIFFERS FROM THE FITTED READ.
    ngh = rules.n_gh or N_GH
    if x0 is None:
        x0 = np.array([float(rules.eval_read(k, 0, np.atleast_2d(S))[0]) for k in SOLVE7])

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
    # THE BANK LENDING SPREAD lambda_K*mu/alpha, ANNUALISED bp.
    return 4e4 * cal[f"lambda_K_{c}"] * o[f"mu_{c}"] / max(o[f"alpha_{c}"], 1e-6)


def report_rest_point(rules, cal, ss, sproc):
    # PRINT THE REST POINT NEXT TO THE DETERMINISTIC STEADY STATE.
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
    # THE IMPACT RESPONSE AS THE DEFAULT PROBABILITY RISES, FROM THE REST POINT.
    S0 = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    base = read_at(rules, cal, ss, sproc, S0.copy())
    Yb, Cb, Ib, Nb = base["Y_D"], base["C_D"], base["I_D"], base["_x"][0]
    print("\n  IMPACT of priced default risk (deviation from the rest point)")
    base_x = read_exact(rules, cal, ss, sproc, S0.copy())
    print("   pd_q%  pd_a%     Y%    Y_ann%   Y_exact%    C%    hours%     I%    "
          "rdep_bp  spread_bp  r_wc_bp    muD   muD_ex     Q_bD   resid")
    s_hi = float(rules.grid.hi[IS])
    # Y_exact and muD_ex clear the period map exactly (see read_exact)
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
    # THE RISK STATE s FOR A GIVEN DEFAULT PROBABILITY (THE INVERSE LOGISTIC).
    return float(np.log(pd / (1.0 - pd)))


def persistence_irf(rules, cal, ss, sproc, pd_shock=0.0198, T=21, s_shock=None):
    # THE IRF AS s DECAYS, OTHER STATES HELD AT THE REST POINT.
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
    # APPEND ONE QUARTER TO EVERY SERIES.
    for k, v in row.items():
        path[k].append(v)


def _persistence_row(o, base, s_t, cal):
    # ONE QUARTER OF THE PERSISTENCE IRF.
    pq = float(default_prob(s_t))
    dY = 100 * (o["Y_D"] / base["Y_D"] - 1)
    return {"pd": 100 * pq, "pd_ann": 100 * ann_prob(pq), "Y": dY, "Y_ann": ann_pct(dY),
            "C": 100 * (o["C_D"] / base["C_D"] - 1),
            "I": 100 * (o["I_D"] / base["I_D"] - 1),
            "N": 100 * (o["_x"][0] / base["_x"][0] - 1),
            "spread": _spread_bp(o, cal), "rdep": bp_ann(o["rdep_D"]),
            "r_wc": bp_ann(o["r_wc_D"]), "Q_bD": o["Q_bD"]}


def advance(o, S, sproc, grid=None):
    # ONE STEP OF THE LAW OF MOTION; EVERY IRF AND THE REST POINT USE THIS.
    Sn = S.copy()
    Sn[IK_D], Sn[IK_F] = o["_x"][2], o["_x"][3]
    Sn[IP_D], Sn[IP_F] = o["Pp_D"], o["Pp_F"]
    Sn[IBDD], Sn[IBDF] = o["b_D_D_new"], o["b_D_F_new"]
    Sn[IBFD] = o["b_F_D_new"]
    Sn[IV] = o["Vp_dep"]
    Sn[IM], Sn[IO] = o["M_cb_new"], o["Op_cb"]
    Sn[IS] = (1 - sproc["rho_s"]) * sproc["s_star"] + sproc["rho_s"] * S[IS]
    return Sn if grid is None else grid.clip(Sn)[0]


def stochastic_rest_point(rules, cal, ss, sproc, tol=1e-11, max_it=4000, verbose=True):
    # WHERE THE SOLVED MODEL RESTS WITHOUT SHOCKS (BOCOLA'S generate_irf.m, STEP 1).
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


def _no_shock_path(rules, cal, ss, sproc, S0, T, read=read_at):
    # THE UNSHOCKED REFERENCE PATH FROM S0.
    ref, Sr = [], S0.copy()
    for _ in range(T):
        o_r = read(rules, cal, ss, sproc, Sr)
        ref.append(o_r)
        Sr = advance(o_r, Sr, sproc, rules.grid)
    Sr = S0.copy()
    K_ref = [Sr[IK_D]]
    for t in range(T - 1):
        Sr = advance(ref[t], Sr, sproc, rules.grid)
        K_ref.append(Sr[IK_D])
    return ref, K_ref


def _dynamic_row(o, r, S, K_ref_t, cal, ss):
    # ONE QUARTER OF THE DYNAMIC IRF, AGAINST THE SAME QUARTER UNSHOCKED.
    pq = float(default_prob(S[IS]))
    dY = 100 * (o["Y_D"] / r["Y_D"] - 1)
    _yD = cal["delta_b_D"] * (1.0 - o["Q_bD"]) / o["Q_bD"]
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
            # the loan spread r_wc - rdep = lambda*mu/E[Omega], the calibration table's lending spread
            "loan_D": bp_ann(o["r_wc_D"] - o["rdep_D"]),
            "loan_F": bp_ann(o["r_wc_F"] - o["rdep_F"]),
            "d_loan_D": bp_ann((o["r_wc_D"] - o["rdep_D"]) - (r["r_wc_D"] - r["rdep_D"])),
            "d_loan_F": bp_ann((o["r_wc_F"] - o["rdep_F"]) - (r["r_wc_F"] - r["rdep_F"])),
            "K": 100 * (S[IK_D] / K_ref_t - 1),
            "n": 100 * (o["n_D"] / r["n_D"] - 1),
            "Q_bD": o["Q_bD"], "Q_bF": o["Q_bF"],
            "dQ_bD": 100 * (o["Q_bD"] / r["Q_bD"] - 1),
            "mu": o["mu_D"], "E_Om": o["E_Om_D"],
            "sov_bp": bp_ann(_yD - _yF),
            # the TPI footprint: purchases and book in % of SS debt, P&L in % of quarterly GDP
            "m_cb": 100 * o["m_cb"] / cal["B_gov_D_ss"],
            "M_cb": 100 * S[IM] / cal["B_gov_D_ss"],
            "Pi_cb": 100 * o["Pi_cb"] / ss["ss_firm_D"]["Y_ss"]}


def _print_dynamic_summary(path, rules, cal, ss, sproc, S0, S, pd_shock, T, escapes):
    # THE TROUGH, THE FITTED-vs-EXACT IMPACT, THE WEDGE AND THE BENCHMARKS.
    tr = min(path["Y"])
    print(f"   trough GDP = {tr:+.4f}% level = {ann_pct(tr):+.4f}% annualised "
          f"(Bocola Table 5 units)")
    # the same impact, cleared exactly
    Sx = S0.copy(); Sx[IS] = s_from_pd(pd_shock)
    bx = read_exact(rules, cal, ss, sproc, S0.copy())
    ox = read_exact(rules, cal, ss, sproc, Sx)
    trx = 100 * (ox["Y_D"] / bx["Y_D"] - 1)
    print(f"   impact GDP brackets [{min(path['Y'][0], trx):+.4f}%, "
          f"{max(path['Y'][0], trx):+.4f}%]: {path['Y'][0]:+.4f}% reading the fitted "
          f"rules, {trx:+.4f}% clearing the period map exactly at that state "
          f"(mu {bx['mu_D']:.5f} -> {ox['mu_D']:.5f})")
    # output moves with the working-capital wedge: the credit spread net of the deposit rate
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
    # REPORT ANY QUARTERS THAT LEFT THE COLLOCATION BOX.
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


def dynamic_irf(rules, cal, ss, sproc, pd_shock=0.0198, T=25, rest_verbose=True,
                exact=False):
    # THE DYNAMIC IRF FROM THE REST POINT, NET OF THE UNSHOCKED PATH.
    # exact=True solves each quarter exactly: fitted reads break the TPI's corner
    read = read_exact if exact else read_at
    S0 = stochastic_rest_point(rules, cal, ss, sproc, verbose=rest_verbose)
    S = S0.copy()
    escapes = []
    ref, K_ref = _no_shock_path(rules, cal, ss, sproc, S0, T, read)
    S[IS] = s_from_pd(pd_shock)
    print(f"\n  DYNAMIC IRF (states evolve; p^d shock to {100*pd_shock:.2f}%/qtr, "
          f"rho_s = {sproc['rho_s']})")
    print("   qtr  pd_q%  pd_a%    GDP%   GDP_ann%     C%      I%   hours%  "
          "rdep_bp  bank_bp  r_wc_bp  sov_bp    K%      n%")
    path = {k: [] for k in ("pd", "pd_ann", "Y", "Y_ann", "C", "I", "N", "spread",
                            "rdep", "r_wc", "K", "n", "Q_bD", "Q_bF", "sov_bp",
                            "d_rdep", "d_spread", "d_r_wc", "dQ_bD",
                            "mu", "E_Om", "Y_F", "C_F", "n_F", "spread_F", "I_F",
                            "m_cb", "M_cb", "Pi_cb", "loan_D", "loan_F", "d_loan_D", "d_loan_F",
                            "at_wall")}
    walled = False
    for t in range(T):
        o = read(rules, cal, ss, sproc, S)
        _append(path, _dynamic_row(o, ref[t], S, K_ref[t], cal, ss))
        # 1 when this quarter's state was put back on the box edge, not reached by the law of motion
        path["at_wall"].append(float(walled))
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
        walled = bool(esc.any())
        if esc.any():
            # the TPI book is zero at rest, so its escapes are reported in levels
            escapes.append((t, [(_SNAMES[i], *((100 * (Sn[i] / S0[i] - 1),
                                                100 * (S[i] / S0[i] - 1))
                                               if abs(S0[i]) > 1e-12 else (Sn[i], S[i])))
                                for i in np.flatnonzero(esc)]))
    _print_dynamic_summary(path, rules, cal, ss, sproc, S0, S, pd_shock, T, escapes)
    return {k: np.array(v) for k, v in path.items()}


def _tfp_read(rules, cal, ss, sproc, S):
    # READ THE NO-DEFAULT RULES AT S (THE TFP EXPERIMENT).
    Sm = np.atleast_2d(S)
    x = np.array([float(rules.eval(k, 0, Sm)[0]) for k in SOLVE7])
    res, o = point_residuals(S, 0, x, rules, cal, ss, sproc,
                             n_gh=rules.n_gh or N_GH, no_default=True)
    o["_x"] = x
    return o


def solve_tfp(cal, ss, sproc, mu=1):
    # SOLVE THE NO-DEFAULT RULES FOR THE TFP EXPERIMENT.
    grid = normal_box(ss, cal, mu=mu)
    rules = RuleSet.from_ss(grid, ss, cal)
    rules.n_gh = N_GH
    solve_no_default(rules, cal, ss, sproc, "TFP d0", True)
    return rules


def _tfp_row(o, base, z_t, cal):
    # ONE QUARTER OF THE TFP IRF.
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
    # THE TFP IRF ALONG THE Z_D DECAY PATH, OTHER STATES AT THE SS.
    S0 = ss_state(ss, cal, sproc)
    Z_ss = S0[IZ]
    base = _tfp_read(rules, cal, ss, sproc, S0.copy())
    print(f"\n  TFP IRF (one-off {dz:.0%} shock, rho_z={sproc['rho_z']} decay)")
    print("   qtr    Z%     Y_D%     C_D%     I_D%    hours%")
    # the same series as dynamic_irf, so both figures share one panel layout
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
    import run
    cal["nw_floor_frac"] = run.NW_FLOOR  # the pipeline's floor
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print("=== recursive global solution: pass-through of sovereign risk ===")
    rules = solve_recursive(cal, ss, sproc)
    impact_table(rules, cal, ss, sproc)
    persistence_irf(rules, cal, ss, sproc)


if __name__ == "__main__":
    main()
