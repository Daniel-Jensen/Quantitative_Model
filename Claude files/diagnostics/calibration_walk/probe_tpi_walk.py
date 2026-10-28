# THE COARSE TPI CAP WALK WITH A CHOSEN SMOOTHING AND TPI-BOOK BAND, RECORDING EVERY SOLVED CAP (2026-10-25). From the
# saved re-centred coarse no-TPI baseline, on the box with BOX_KW m_band = M_BAND (the baseline is only re-seeded, not
# re-solved, when the band changes: the TPI economy alone uses the wider box), the cap is walked from 700 bp in 50 bp
# rungs (a failed rung halves the step, down to 6.25 bp) with cal["tpi_eps"] = EPS, in the FB form (_solve_tpi's walk).
# At each solved cap: the rest-point sovereign spread and, for shocks to p^d 1% and 2%, the impact purchase (% of B) and
# the quarters the exactly cleared 25-quarter path spends at the box wall. Saves runs/<TAG>_cap<c>.pkl.
# usage: python3 -u probe_tpi_walk.py TAG M_BAND EPS      e.g. t3 0.32 1e-2
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)


def reads(rules, c, ss, sp, B, RE, IS, np):
    # REST SPREAD, AND THE IMPACT PURCHASE AND WALL QUARTERS FOR TWO SHOCKS.
    c = dict(c)  # read in the form c says (cal["tpi_gz"]), which _tpi_form keeps in step with the rules
    yld = lambda Q, db: 4e4 * db * (1 - Q) / Q
    rules._rest_point = None
    S0 = RE.stochastic_rest_point(rules, c, ss, sp, verbose=False)
    o0 = RE.read_exact(rules, c, ss, sp, S0.copy())
    out = [f"rest sov {yld(o0['Q_bD'], c['delta_b_D']) - yld(o0['Q_bF'], c['delta_b_F']):6.1f} bp"]
    for pd in (0.01, 0.02):
        S = S0.copy(); S[IS] = RE.s_from_pd(pd)
        walls, m0 = 0, None
        for t in range(25):
            o = RE.read_exact(rules, c, ss, sp, S.copy())
            m0 = o["m_cb"] if m0 is None else m0
            Sn = RE.advance(o, S, sp); S = rules.grid.clip(Sn)[0]
            walls += bool(np.any(np.abs(Sn - S) > 1e-12))
        out.append(f"p^d {100 * pd:.0f}%: purchase {100 * m0 / B:5.1f}% of B, wall {walls:2d}/25")
    return "  ".join(out)


def main():
    # WALK THE CAP AND READ EVERY SOLVED RUNG.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS, IM
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.decision_rules import RuleSet
    import global_projection.solver_recursive.recursive_experiment as RE
    tag, m_band, eps = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    t0 = time.perf_counter()
    base = pickle.load(open(os.path.join(ROOT, "results", "GLOBAL", "data", "rules_base_coarse.pkl"), "rb"))
    cal = get_calibration(); cal["nw_floor_frac"] = 0.15
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    B = cal["B_gov_D_ss"]
    RE.BOX_KW["m_band"] = m_band
    rules = RE._seed_from(RuleSet(RE.normal_box(ss, cal, at=base.box["at"]), 2), base)
    rules.box = dict(base.box, box_kw=dict(RE.BOX_KW))
    M_top = rules.grid.from_unit(np.eye(rules.grid.d)[IM][None])[0][IM]
    print(f"[config] {tag}: m_band {m_band} (book top {100 * M_top / B:.1f}% of B), tpi_eps {eps}", flush=True)
    c = dict(cal, tpi_on=True, tpi_eps=eps)
    RE._tpi_form(rules, c, ss, gz=False)
    cap, step, last, good = 700.0, 50.0, None, rules.copy()
    while True:
        c["tpi_cap_bp"] = cap
        ok, it, w = RE._stage(rules, c, ss, sp, (0, 1), False, f"cap={cap:.1f}bp", True, warm=0)
        print(f"[rung] cap {cap:6.1f} ok={ok} worst={w:.2e} its={it} ({time.perf_counter() - t0:.0f}s)", flush=True)
        if ok:
            good, last = rules.copy(), cap
            print(f"[read] cap {cap:6.1f}: {reads(rules.copy(), c, ss, sp, B, RE, IS, np)}", flush=True)
            pickle.dump(dict(rules=good, cap=cap, eps=eps, m_band=m_band),
                        open(os.path.join(HERE, "runs", f"{tag}_cap{cap:.0f}.pkl"), "wb"))
            if cap <= 200.0:
                break
            cap = max(cap - step, 200.0)
            continue
        step *= 0.5
        if last is None or step < 6.25:
            break
        rules = good.copy()
        cap = max(last - step, 200.0)
    print(f"[END] {tag}: last cap solved {last} ({time.perf_counter() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
