# B, THE BINDING CALIBRATION (2026-10-31): THE D SS LENDING-SPREAD INSTRUMENT (credit_spread_target_D) RAISED IN RUNGS
# FROM APPENDIX B's 5 bp/yr, EACH RUNG A WARM FINE HARD SOLVE (m2_warm.solve_at) WITH THE OTHER INSTRUMENTS AND f_D HELD.
# AT EACH RUNG, EXACT READS OF WHAT THE PAPER NEEDS: (1) THE D IC BINDING AT THE REST POINT, (2) THE ACCEPTANCE SIGNS OF
# THE p^d 2% IMPACT (Q_bD, n_D, n_F, Y_D, C_D FALL, THE LENDING SPREAD RISES), (3) REAL EFFECTS OF A PURCHASE THERE (FIX A:
# RESERVES NOT DIVERTABLE); THEN A BISECTION FOR THE SMALLEST INSTRUMENT AT WHICH ALL THREE HOLD. THE FIRST RUNG (5 bp)
# IS THE APPENDIX-B SOLUTION RE-SOLVED UNDER FIX A.
#   python3 -u b_spread_sweep.py [bp,bp,...]
import json, os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np

RUNGS = (5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0)
BISECT_TOL = 5.0
SHOCK_PD = 0.02


def measure(RE, M):
    # THE REST POINT, THE p^d 2% IMPACT AND THE PURCHASE EFFECTS AT THE CURRENT SOLUTION (EXACT READS).
    from m2_ob import slack, threshold
    from tpi_purchase_tests import buy
    from global_projection.solver_recursive.calibrate_stochastic import rest_point_moments
    from global_projection.solver_recursive.state_grid import IS, STATE_NAMES
    r, cal, ss, sp = (M.CUR[k] for k in ("rules", "cal", "ss", "sp"))
    S0 = RE.stochastic_rest_point(r, cal, ss, sp, verbose=False)
    S2 = S0.copy(); S2[IS] = RE.s_from_pd(SHOCK_PD)
    o0, o2 = RE.read_exact(r, cal, ss, sp, S0), RE.read_exact(r, cal, ss, sp, S2)
    spread = lambda o: 4e4 * (o["r_wc_D"] - o["rdep_D"])
    sov = lambda o: 4e4 * (cal["delta_b_D"] * (1 / o["Q_bD"] - 1) - cal["delta_b_F"] * (1 / o["Q_bF"] - 1))
    m = rest_point_moments(o0)
    u = r.grid.to_unit(S0[None])[0]
    p_star, where = threshold(RE, M)
    out = dict(mu_rest=float(o0["mu_D"]), spread_rest=spread(o0), slack_rest=float(slack(RE, r, cal, ss, sp, S0)),
               E_Om_rest=float(o0["E_Om_D"]), lev_D=float(m["lev_D"]), exp_D=float(m["exp_D"]), sov_rest=sov(o0),
               outside=[k for k, ui in zip(STATE_NAMES, u) if abs(ui) > 1.0], threshold=float(p_star), where=where,
               mu_2=float(o2["mu_D"]), spread_rise_2=spread(o2) - spread(o0),
               **{f"d{k}_2": 100 * (o2[k] / o0[k] - 1) for k in ("Q_bD", "n_D", "n_F", "Y_D", "C_D", "I_D")})
    for tag, S in (("rest", S0), ("2", S2)):
        for share in (0.05, 0.15):
            a, b, ok = buy(r, cal, ss, sp, S, share)
            key = f"buy{100 * share:.0f}_{tag}"
            out.update({f"{key}_price": 100 * (b["Q_bD"] / a["Q_bD"] - 1), f"{key}_spread": spread(b) - spread(a),
                        f"{key}_Y": 100 * (b["Y_D"] / a["Y_D"] - 1), f"{key}_I": 100 * (b["I_D"] / a["I_D"] - 1),
                        f"{key}_ok": bool(ok)})
    out["binding_rest"] = out["mu_rest"] > 1e-6
    out["signs_2"] = (all(out[f"d{k}_2"] < 0 for k in ("Q_bD", "n_D", "n_F", "Y_D", "C_D"))
                      and out["spread_rise_2"] > 0)
    # a purchase must relieve the constraint measurably: the lending spread down by at least 1 bp, output up
    out["real_2"] = out["buy5_2_spread"] <= -1.0 and out["buy5_2_Y"] > 0
    out["paper"] = out["binding_rest"] and out["signs_2"] and out["real_2"]
    return out


def row(bp, o):
    # ONE PRINTED LINE PER RUNG.
    return (f"[B] {bp:6.1f} bp: rest mu_D {o['mu_rest']:.4f} spread {o['spread_rest']:6.1f} bp slack {o['slack_rest']:.3f} "
            f"E[Om] {o['E_Om_rest']:.3f} lev {o['lev_D']:.2f} exp {100 * o['exp_D']:.1f}% sov {o['sov_rest']:.0f} bp | binds "
            f"from p^d {100 * o['threshold']:.2f}% ({o['where']}) | p^d 2%: spread {o['spread_rise_2']:+.0f} bp, Q_bD "
            f"{o['dQ_bD_2']:+.2f}%, n_D {o['dn_D_2']:+.2f}%, n_F {o['dn_F_2']:+.3f}%, Y_D {o['dY_D_2']:+.3f}%, C_D "
            f"{o['dC_D_2']:+.3f}%, I_D {o['dI_D_2']:+.2f}% | buy 5% at 2%: price {o['buy5_2_price']:+.2f}%, spread "
            f"{o['buy5_2_spread']:+.0f} bp, Y {o['buy5_2_Y']:+.3f}%; buy 15%: spread {o['buy15_2_spread']:+.0f} bp, Y "
            f"{o['buy15_2_Y']:+.3f}% | binding {o['binding_rest']} signs {o['signs_2']} real {o['real_2']} -> PAPER "
            f"{o['paper']}; outside the box: {o['outside'] or 'none'}")


def main():
    from calibration.global_projection import get_calibration
    import global_projection.solver_recursive.calibrate_stochastic as CS
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import s_process_params, IS
    from diag_path import FLOOR
    import m2_warm as M
    rungs = tuple(float(v) for v in sys.argv[1].split(",")) if len(sys.argv) > 1 else RUNGS
    M.REFINE = (IS, 9, "cross")
    sv = pickle.load(open(HERE / "runs" / "appendixB_final.pkl", "rb"))
    base = get_calibration()
    base.update(nw_floor_frac=FLOOR, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3, jacobian="complex")
    M.CUR.update(rules=sv["rules"], cal=sv["cal"], ss=sv["ss"], sp=s_process_params(sv["cal"]),
                 vals=dict(sv["vals"], **M.HARD), base=base)
    print(f"[B] start: s box [{sv['rules'].grid.lo[IS]:.4f}, {sv['rules'].grid.hi[IS]:.4f}]", flush=True)
    x0 = np.array([M.CUR["cal"][p[0]] for p in CS.PAIRS])
    i = [p[0] for p in CS.PAIRS].index("credit_spread_target_D")
    hist, sols = [], {}
    log = HERE / "runs" / "b_spread_sweep.json"

    def solve(bp, start):
        # THE WARM SOLVE AT bp FROM THE SOLVED RUNG start (None: THE CURRENT SOLUTION); MEASURED AND SAVED ON SUCCESS.
        if start is not None:
            M.CUR.update(sols[start])
        x = x0.copy(); x[i] = bp / 4e4
        ok, w = M.solve_at(RE, x)
        if not ok:
            print(f"[B] {bp:.1f} bp: solve failed ({w:.1e})", flush=True)
            return None
        sols[bp] = {k: M.CUR[k] for k in ("rules", "cal", "ss", "sp", "vals")}
        o = dict(measure(RE, M), bp=bp, solve=float(w))
        hist.append(o)
        print(row(bp, o), flush=True)
        if bp == rungs[0]:
            print(f"[B] s box now [{M.CUR['rules'].grid.lo[IS]:.4f}, {M.CUR['rules'].grid.hi[IS]:.4f}]", flush=True)
        with open(HERE / "runs" / f"b_sweep_{bp:.0f}bp.pkl", "wb") as fh:
            pickle.dump({k: M.CUR[k] for k in ("rules", "cal", "ss", "vals")}, fh)
        log.write_text(json.dumps(hist, indent=1, default=lambda z: z.item()))
        return o

    # 1. the ladder, a failed rung retried once from halfway
    last = None
    for bp in rungs:
        o = solve(bp, last)
        if o is None and last is not None:
            mid = 0.5 * (last + bp)
            if solve(mid, last) is None:
                break
            last = mid
            o = solve(bp, mid)
        if o is None:
            break
        last = bp
    # 2. the smallest instrument at which the paper goes through, by bisection between solved rungs
    done = sorted((o["bp"], o["paper"]) for o in hist)
    passing = [bp for bp, ok in done if ok]
    if not passing:
        print("[B] no solved rung satisfies all three", flush=True)
        return
    hi = passing[0]
    below = [bp for bp, ok in done if not ok and bp < hi]
    lo = below[-1] if below else None
    while lo is not None and hi - lo > BISECT_TOL:
        mid = 0.5 * (lo + hi)
        near = min(sols, key=lambda b: abs(b - mid))
        o = solve(mid, near)
        if o is None:
            break
        lo, hi = (lo, mid) if o["paper"] else (mid, hi)
    print(f"[B] the paper goes through from an SS spread instrument of {hi:.1f} bp/yr "
          f"(fails at {lo if lo is None else f'{lo:.1f}'} bp)", flush=True)


if __name__ == "__main__":
    main()
