# THE SOVEREIGN-RISK EXPERIMENT: THE SOLVE LADDER, THE READS AND THE IMPULSE RESPONSES.
import numpy as np
from scipy.optimize import root

from calibration.global_projection import get_calibration
from global_projection.steady_state import solve_steady_state
from global_projection.solver_recursive.state_grid import build_state_box, s_process_params, default_prob
from global_projection.solver_recursive.state_grid import (IK_D, IK_F, IP_D, IP_F, IBDD,
                                          IBDF, IBFD, IV, IS, IZ, IM, IO, STATE_NAMES)
from global_projection.solver_recursive.decision_rules import RuleSet, SOLVE7, STORE_RULES
from global_projection.solver_recursive.collocation import solve_collocation
from global_projection.solver_recursive.recursive_main import (time_iteration, calibrate_household_anchors,
                            ss_state, p_block_rotation)
from global_projection.reporting.prints import (bp_ann, ann_pct, ann_prob,
                              BOCOLA_IRF_CLOSED, BOCOLA_IRF_OPEN,
                              BOCOLA_EPISODE_LEVEL)
from global_projection.solver_recursive.point_map import point_residuals, _tpi_floor


# state names for the box-escape report
_SNAMES = STATE_NAMES

# the collocation box, shared by every experiment
BOX_KW = dict(k_band=0.02, p_band_D=0.04, p_band_F=0.04, b_band=0.12, w_band=0.04)

# quadrature order; every reader reuses rules.n_gh
N_GH = 5

# dense Chebyshev nodes in s: 5 is converged for every reported object, 9 is Bocola's
S_REFINE = 5

# time-iteration sweeps that put each Newton inside its basin
WARM_SWEEPS = 12
# the refined grid needs a longer warm start
REFINE_WARM_SWEEPS = 20

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
           backend="auto", warm=WARM_SWEEPS, maxit=40):
    # ONE SOLVE STAGE: A SHORT WARM START, THEN THE GLOBAL NEWTON.
    if warm:
        time_iteration(rules, cal, ss, sproc, regimes=regimes, no_default=no_default,
                       damp=0.5, tol=1e-4, max_it=warm, n_gh=N_GH, verbose=False)
    return solve_collocation(rules, cal, ss, sproc, regimes=regimes,
                             no_default=no_default, n_gh=N_GH, backend=backend,
                             maxit=maxit, verbose=verbose, label=f" {label}")


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


def _build_rules(cal, ss, sproc, mu, mu_vec, rotate, nreg, verbose):
    # THE COARSE GRID (OPTIONALLY P-ROTATED) AND A STEADY-STATE COLD START.
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
    # THE NO-TPI BASELINE: d = 0, THEN THE HAIRCUT HOMOTOPY, THEN THE JOINT SOLVE.
    ok0, it0, w0 = _stage(rules, cal, ss, sproc, (0,), True, "d0", verbose)
    for k in STORE_RULES:  # start the default regime from d = 0
        rules.set_values(k, D_REG, rules.vals[k][0].copy())
    rec_target = cal["recovery_rate_D"]
    ok1, w1 = False, np.nan
    for rec in (0.85, 0.70, 0.55, rec_target):
        cal["recovery_rate_D"] = rec
        ok1, it1, w1 = _stage(rules, cal, ss, sproc, (D_REG,), False,
                              f"d1 rec={rec:.2f}", verbose)
    cal["recovery_rate_D"] = rec_target
    okb, itb, wb = _stage(rules, cal, ss, sproc, (0, D_REG), False,
                          "joint (no TPI)", verbose, backend=backend)
    return (ok0, ok1, okb), (w0, w1, wb)


def _tpi_form(rules, cal, ss, gz):
    # SWITCH THE TPI VARIABLE BETWEEN THE FB (SOLVE) AND GZ (READ) FORMS AT THE NODES.
    x = rules.vals["x_cb"][0]
    if gz:
        gap = (rules.vals["Q_bD"][0] - _tpi_floor(rules.vals["Q_bF"][0], cal)) / ss["Q_bD_ss"]
        x = np.where(x > gap, x, -np.maximum(gap, 0.0))
    else:
        x = np.maximum(x, 0.0)
    rules.set_values("x_cb", 0, x)
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


def _refine_s(rules, cal, ss, sproc, box, s_refine, backend, verbose):
    # REFINE THE s DIMENSION STEP BY STEP, EACH STEP SEEDED FROM THE LAST.
    nreg = rules.n_regimes
    okj = wj = None
    tpi = bool(cal["tpi_on"])
    ladder = [m for m in (5, 9, 17) if 1 < m < s_refine] + [s_refine]
    # with the TPI on, start from 3 nodes in s, which seeds the 5-node solve well
    if tpi and s_refine > 3:
        ladder = [3] + ladder
    for m_s in ladder:
        gfine = build_state_box(ss, cal, mu=box["mu"], mu_vec=box["mu_vec"],
                                rot=box["rot"], centre=box["centre"], refine=(IS, m_s),
                                **box["box_kw"])
        if verbose:
            print(f"  s-refined grid: {gfine.n} points x {nreg} regimes "
                  f"({m_s} nodes, degree {m_s - 1} in s)")
        fine = _seed_from(RuleSet(gfine, nreg), rules)
        if not tpi:
            okj, itj, wj = _stage(fine, cal, ss, sproc, tuple(range(nreg)), False,
                                  f"joint (s={m_s})", verbose, backend=backend,
                                  warm=REFINE_WARM_SWEEPS, maxit=12)
        else:  # root in the FB form, easing the smoothing in, then hand back in GZ
            _tpi_form(fine, cal, ss, gz=False)
            eps, warm = cal["tpi_eps"], REFINE_WARM_SWEEPS
            for e in [e for e in TPI_EPS_LADDER if e > eps] + [eps]:
                cal["tpi_eps"] = e
                okj, itj, wj = _stage(fine, cal, ss, sproc, tuple(range(nreg)), False,
                                      f"joint (s={m_s}, eps={e:.0e})", verbose,
                                      backend=backend, warm=warm, maxit=20)
                warm = 0
            cal["tpi_eps"] = eps
            ok_gz, _, wj = _tpi_polish(fine, cal, ss, sproc, f"joint (s={m_s})", verbose,
                                       backend)
            okj = okj and ok_gz
        rules = fine
    return rules, okj, wj


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
    rules, box = _build_rules(cal, ss, sproc, mu, mu_vec, rotate, nreg, verbose)
    D_REG = 1  # the regime index is the default indicator
    tpi = bool(cal["tpi_on"])

    if base is None:
        cal["tpi_on"] = False  # the baseline is always the no-TPI economy
        (ok0, ok1, okb), (w0, w1, wb) = _solve_baseline(rules, cal, ss, sproc, D_REG,
                                                        backend, verbose)
        cal["tpi_on"] = tpi
        if base_out is not None:  # a reusable snapshot for the caller
            base_out.append(rules.copy())
    else:
        assert base.n_regimes == nreg and base.grid.n == rules.grid.n, \
            "reused baseline must carry the same grid and regime count"
        rules = base.copy()
        ok0 = ok1 = okb = True
        w0 = w1 = wb = np.nan
        if verbose:
            print("  reusing the solved baseline")

    cap_solved, target = None, cal["tpi_cap_bp"]
    if tpi:  # with the TPI on, the cap walk replaces the joint polish
        rules, okj, wj = _solve_tpi(rules, cal, ss, sproc, backend, verbose)
        cap_solved = rules.tpi_cap_solved
    else:
        okj, itj, wj = _stage(rules, cal, ss, sproc, tuple(range(nreg)), False, "joint",
                              verbose, backend=backend)

    if s_refine and s_refine > 1:
        # a cap walk that stopped short refines at the cap it reached
        if tpi and cap_solved is not None:
            cal["tpi_cap_bp"] = cap_solved
        rules, okj2, wj = _refine_s(rules, cal, ss, sproc, box, s_refine, backend, verbose)
        okj = okj2 and (okj or not tpi)  # a cap walk short of its target counts as a failure
        cal["tpi_cap_bp"] = target
    rules.tpi_cap_solved = cap_solved

    _stamp_verdict(rules, (ok0, ok1, okb, okj), (w0, w1, wb, wj))
    return rules


def read_at(rules, cal, ss, sproc, S):
    # READ THE FITTED RULES AT STATE S.
    Sm = np.atleast_2d(S)
    x = np.array([float(rules.eval(k, 0, Sm)[0]) for k in SOLVE7])
    res, o = point_residuals(S, 0, x, rules, cal, ss, sproc,
                             n_gh=rules.n_gh or N_GH, no_default=False)
    o["_x"] = x
    o["_resid"] = float(np.max(np.abs(res)))
    return o


def read_exact(rules, cal, ss, sproc, S, x0=None):
    # CLEAR THE PERIOD MAP EXACTLY AT S; NEAR THE mu KINK THIS DIFFERS FROM THE FITTED READ.
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
                            "m_cb", "M_cb", "Pi_cb")}
    for t in range(T):
        o = read(rules, cal, ss, sproc, S)
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
    grid = build_state_box(ss, cal, mu=mu, **BOX_KW)
    rules = RuleSet.from_ss(grid, ss, cal)
    rules.n_gh = N_GH
    _stage(rules, cal, ss, sproc, (0,), True, "TFP d0", True)
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
    cal["nw_floor_frac"] = 0.15  # as in run.py
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print("=== recursive global solution: pass-through of sovereign risk ===")
    rules = solve_recursive(cal, ss, sproc)
    impact_table(rules, cal, ss, sproc)
    persistence_irf(rules, cal, ss, sproc)


if __name__ == "__main__":
    main()
