# M1 (2026-10-26): THE COARSE RISK SOLVE AT APPENDIX B, FROM THE d0 CONTINUATION PATH'S ENDPOINT: THE DEFAULT
# REGIME (HAIRCUT LADDER), THE JOINT SOLVE AT eps 1e-2, THE SMOOTHING WALKED DOWN WITH THE SS RE-SOLVED AT
# EVERY RUNG, THEN THE GENUINE CONDITIONS (IC MAX, BOOK PAIRS, FLOOR); READS AT WHATEVER WAS REACHED.
#   python3 -u m1_coarse.py <d0 endpoint checkpoint>   (e.g. diag_path_wedge.pkl)
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np
from diag_homotopy import economy
from diag_path import relative_seed
from measure import omega_legs, row, PDS
from diag_path import FLOOR

EPS_RUNGS = (5e-3, 2.5e-3, 1e-3, 5e-4, 2.5e-4, 1e-4)


def save(tag, **kw):
    # A CHECKPOINT IN runs/.
    with open(HERE / "runs" / f"m1{os.environ.get('M1_TAG', '')}_{tag}.pkl", "wb") as fh:
        pickle.dump(kw, fh)


def stage(RE, rules, cal, ss, sp, label, regimes=(0, 1), maxit=30):
    # ONE JOINT NEWTON STAGE, REPORTED.
    ok, it, w = RE._stage(rules, cal, ss, sp, regimes, False, label, True, warm=0, maxit=maxit)
    print(f"[{label}] ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
    return ok, w


def calm_guard(RE, rules, cal, ss, sp, label, maxit=60):
    # AFTER A JOINT STAGE: ANY NORMAL-TIMES NODE ON THE CRISIS BRANCH IS RE-SELECTED, POLISHED AND RE-SOLVED ONCE
    # (2026-10-28 PLAN, STEP 2); RETURNS (ok, max|F|) OF THE RE-SOLVE, OR None WHEN NO NODE WAS ON IT.
    import branch_select as B
    bad = B.crisis_nodes(rules, cal, ss, sp, False)
    print(f"[{label}] crisis nodes {bad}", flush=True)
    if not bad:
        return None
    B.calm_scan(rules, cal, ss, sp, (0, 1), False, only=[b[0] for b in bad])
    B.polish_sweeps(rules, cal, ss, sp, (0, 1), False)
    return stage(RE, rules, cal, ss, sp, f"{label} calm re-solve", maxit=maxit)


def reads(RE, rules, cal, ss, sp, tag):
    # REST POINT, ITS BOX POSITION, AND THE IMPACT ROWS WITH THE E[OMEGA] LEGS, FITTED AND EXACT.
    from global_projection.solver_recursive.state_grid import STATE_NAMES, IS
    S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)
    u = rules.grid.to_unit(np.atleast_2d(S0))[0]
    print(f"\n[{tag}] rest point box position u: " +
          ", ".join(f"{k} {v:+.2f}" for k, v in zip(STATE_NAMES, u) if abs(v) > 0.05), flush=True)
    for name, read in (("fitted", RE.read_at), ("exact", RE.read_exact)):
        b = omega_legs(rules, cal, ss, sp, S0, read)
        o = b["o"]
        lev = (o["dep_D"] + o["n_D"]) / o["n_D"]
        print(f"[{tag}] {name} rest: mu_D {b['mu_D']:.5f}  loan spread {4e4 * (o['r_wc_D'] - o['rdep_D']):.1f}bp  "
              f"leverage {lev:.3f}  sov spread "
              f"{4e4 * (cal['delta_b_D'] * (1 / o['Q_bD'] - 1) - cal['delta_b_F'] * (1 / o['Q_bF'] - 1)):.1f}bp  "
              f"E[Om] {b['E_Om']:.4f}  goods_F {o['goods_F']:+.1e}", flush=True)
        for pd in PDS:
            S = S0.copy(); S[IS] = RE.s_from_pd(pd)
            g = omega_legs(rules, cal, ss, sp, S, read)
            r = row(g["o"], o)
            print(f"   p^d {100 * pd:.2f}%: dQ_bD {r['Q_bD']:+.2f}%  dn_D {r['n_D']:+.2f}%  dn_F {r['n_F']:+.2f}%  "
                  f"dY_D {r['Y_D']:+.4f}%  dC_D {r['C_D']:+.4f}%  dspread {r['spread_bp']:+.1f}bp  mu_D {g['mu_D']:.5f}  "
                  f"E[Om] {g['E_Om']:.4f} = nd {g['Om_nd']:.4f}/d {g['Om_d']:.3f} (alpha' nd {g['alpha_nd']:.3f}, "
                  f"d {g['alpha_d']:.3f})", flush=True)


def main():
    from calibration.global_projection import get_calibration
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive import recursive_experiment as RE
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    rules, cal, ss, vals = sv["rules"], sv["cal"], sv["ss"], dict(sv["vals"])
    # the collocation Jacobian (M1_JAC: fd or complex, 2026-10-28)
    jac = os.environ.get("M1_JAC", "fd")
    cal.update(nw_floor_frac=FLOOR, jacobian=jac)
    sp = s_process_params(cal)
    base = get_calibration()
    base.update(nw_floor_frac=FLOOR, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3, jacobian=jac)
    print(f"[load] {sys.argv[1]}: lambda {cal['lambda_K_D']:.5f}  mu_ss {ss['ss_bank_D']['mu_ss']:.5f}  "
          f"spread instrument {4e4 * cal['credit_spread_target_D']:.1f}bp", flush=True)

    # d = 0 re-solved first (the floor became default-only on 2026-10-28; M1_D0=0 skips it)
    if os.environ.get("M1_D0", "1") == "1":
        ok, it, w = RE._stage(rules, cal, ss, sp, (0,), True, "d0", False, warm=0, maxit=40)
        print(f"[d0] ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
    # the default regime from d = 0, then the haircut walked to the calibrated recovery
    for k in rules.vals:
        rules.set_values(k, 1, rules.vals[k][0].copy())
    rec = cal["recovery_rate_D"]
    for r in RE.HAIRCUT_LADDER + (rec,):
        cal["recovery_rate_D"] = r
        ok, w = stage(RE, rules, cal, ss, sp, f"d1 rec={r:.2f}", regimes=(1,))
    cal["recovery_rate_D"] = rec
    ok, w = stage(RE, rules, cal, ss, sp, "joint eps 1e-2", maxit=60)
    g = calm_guard(RE, rules, cal, ss, sp, "joint eps 1e-2")
    if g is not None:
        ok, w = g
    save("joint_eps1e-2", rules=rules, cal=cal, ss=ss, vals=vals)
    if os.environ.get("M1_RECENTRE", "0") == "1":
        rules, ok = recentre_until(RE, rules, cal, ss, sp, vals, ok)
    if not ok:
        reads(RE, rules, cal, ss, sp, "joint eps 1e-2 (unconverged)")
        return
    finish(RE, base, rules, cal, ss, sp, vals)


def recentre_until(RE, rules, cal, ss, sp, vals, ok, tries=3, maxit=60):
    # WHILE THE JOINT STAGE FAILS: RE-CENTRE THE BOX'S BANK STOCKS ON THE CURRENT REST POINT AND RE-SOLVE.
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive.state_grid import _box_centre
    for k in range(tries):
        if ok:
            break
        S0 = RE.stochastic_rest_point(rules, cal, ss, sp, verbose=False)
        at = {key: float(S0[RE._STOCK_INDEX[key]]) for key in RE.BOX_STOCKS}
        c0 = _box_centre(ss, cal)
        print(f"[recentre {k + 1}] " + ", ".join(f"{key} {100 * (v / c0[key] - 1):+.2f}%" for key, v in at.items()),
              flush=True)
        new = RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1, at=at), 2), rules)
        new.box = dict(mu=1, mu_vec=None, rot=None, centre=None, box_kw=dict(RE.BOX_KW), aligned=True, at=at)
        rules = new
        ok, w = stage(RE, rules, cal, ss, sp, f"joint re-centred {k + 1}", maxit=maxit)
        save(f"joint_recentred{k + 1}", rules=rules, cal=cal, ss=ss, vals=vals)
    return rules, ok


def finish(RE, base, rules, cal, ss, sp, vals):
    # FROM A CONVERGED JOINT SOLVE AT eps 1e-2: THE SMOOTHING WALK, THE GENUINE CONDITIONS, THE READS.
    best = ("joint eps 1e-2", rules, cal, ss, sp)

    # the smoothing walked down, the SS re-solved consistently at every rung
    for e in EPS_RUNGS:
        v = dict(vals, mu_eps=e)
        cal_k, ss_k, sp_k = economy(base, v)
        trial = relative_seed(rules, cal, ss, cal_k, ss_k)
        ok, w = stage(RE, trial, cal_k, ss_k, sp_k, f"joint eps {e:.1e}")
        g = calm_guard(RE, trial, cal_k, ss_k, sp_k, f"joint eps {e:.1e}", maxit=30)
        if g is not None:
            ok, w = g
        if not ok:
            break
        rules, cal, ss, sp, vals = trial, cal_k, ss_k, sp_k, v
        best = (f"joint eps {e:.1e}", rules, cal, ss, sp)
        save(f"eps{e:.0e}", rules=rules, cal=cal, ss=ss, vals=vals)

    # the genuine conditions on the hard-max SS: all at once, else one at a time
    hard = dict(mu_eps=0.0, book_eps=0.0, floor_eps=0.0)
    if True:  # attempted from wherever the walk stopped
        for keys in (tuple(hard), ("mu_eps",), ("mu_eps", "book_eps"), ("mu_eps", "book_eps", "floor_eps")):
            v = dict(vals, **{k: hard[k] for k in keys})
            cal_k, ss_k, sp_k = economy(base, v)
            trial = relative_seed(rules, cal, ss, cal_k, ss_k)
            ok, w = stage(RE, trial, cal_k, ss_k, sp_k, f"hard {'+'.join(keys)}")
            g = calm_guard(RE, trial, cal_k, ss_k, sp_k, f"hard {'+'.join(keys)}", maxit=30)
            if g is not None:
                ok, w = g
            if ok:
                rules, cal, ss, sp, vals = trial, cal_k, ss_k, sp_k, v
                best = (f"hard {'+'.join(keys)}", rules, cal, ss, sp)
                save("hard_" + "_".join(keys), rules=rules, cal=cal, ss=ss, vals=vals)
                if len(keys) == 3:
                    break
            elif keys == tuple(hard):
                continue
            else:
                break
    tag, rules, cal, ss, sp = best
    print(f"\n[best] {tag}", flush=True)
    reads(RE, rules, cal, ss, sp, tag)


if __name__ == "__main__":
    main()
