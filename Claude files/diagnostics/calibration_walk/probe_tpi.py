# THE TPI ON THE COARSE GRID WITH A CHOSEN TPI-BOOK BAND AND CAP-WALK START (2026-10-25). From the saved re-centred
# coarse no-TPI baseline (results/GLOBAL/data/rules_base_coarse.pkl): re-seed it on the box with BOX_KW m_band = M_BAND
# and re-solve it (its rest point is compared with the saved baseline's: the M = 0 slice should not move), then walk the
# cap from CAP_START down to 200 bp (recursive_experiment._solve_tpi). A failed rung prints the twelve worst nodes with
# their unit coordinates (the first rung is solved on its own for that). Saves runs/<TAG>_base.pkl and runs/<TAG>_tpi.pkl, then reads the 2% shock: the impact purchase
# and the quarters the exact path spends at the box wall.
# usage: python3 -u probe_tpi.py TAG M_BAND CAP_START     e.g. t2 1.1 700
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)


def worst_nodes(rules, cal, ss, sp, PM, SOLVE, STATE_NAMES, np):
    # THE NODES WITH THE LARGEST PERIOD-MAP RESIDUAL, BOTH REGIMES.
    rows = []
    for d in range(rules.n_regimes):
        for i, S in enumerate(rules.grid.points):
            x = np.array([rules.vals[k][d][i] for k in SOLVE])
            res, _ = PM.point_residuals(S, d, x, rules, cal, ss, sp, n_gh=rules.n_gh, no_default=False)
            rows.append((float(np.max(np.abs(res))), d, i, int(np.argmax(np.abs(res)))))
    u = rules.grid.points_unit
    for r, d, i, e in sorted(rows, reverse=True)[:12]:
        nz = ", ".join(f"{STATE_NAMES[j]} {u[i, j]:+.2f}" for j in np.flatnonzero(np.abs(u[i]) > 1e-12)) or "centre"
        print(f"    |res| {r:.2e}  d={d} node {i:2d} eq {e:2d}  [{nz}]", flush=True)


def main():
    # RE-SOLVE THE BASELINE ON THE BOX, WALK THE CAP, READ THE SHOCK.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS, IM, STATE_NAMES
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.decision_rules import RuleSet, SOLVE
    import global_projection.solver_recursive.point_map as PM
    import global_projection.solver_recursive.recursive_experiment as RE
    tag, m_band, cap0 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    t0 = time.perf_counter()
    base = pickle.load(open(os.path.join(ROOT, "results", "GLOBAL", "data", "rules_base_coarse.pkl"), "rb"))
    cal = get_calibration(); cal["nw_floor_frac"] = 0.15
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    B = cal["B_gov_D_ss"]
    RE.BOX_KW["m_band"] = m_band
    RE.TPI_CAP_START = cap0
    at = base.box["at"]
    S_saved = RE.stochastic_rest_point(base, dict(cal, tpi_on=False), ss, sp, verbose=False)
    rules = RE._seed_from(RuleSet(RE.normal_box(ss, cal, at=at), 2), base)
    rules.box = dict(base.box, box_kw=dict(RE.BOX_KW))
    M_top = rules.grid.from_unit(np.eye(rules.grid.d)[IM][None])[0][IM]
    print(f"[config] {tag}: m_band {m_band} (book top {100 * M_top / B:.1f}% of B), cap walk from {cap0} bp", flush=True)
    c_off = dict(cal, tpi_on=False)
    ok, _, w = RE._stage(rules, c_off, ss, sp, (0, 1), False, f"no-TPI base m_band={m_band}", True, warm=0, maxit=40)
    S_new = RE.stochastic_rest_point(rules, c_off, ss, sp, verbose=False)
    rel = np.max(np.abs(S_new - S_saved)[:8] / np.abs(S_saved[:8]))
    print(f"[base] ok={ok} worst={w:.2e}; rest point vs the saved baseline: max rel. gap {rel:.2e} "
          f"({time.perf_counter() - t0:.0f}s)", flush=True)
    pickle.dump(dict(rules=rules.copy(), ok=ok, worst=w), open(os.path.join(HERE, "runs", f"{tag}_base.pkl"), "wb"))
    # the first rung on its own, so a failure can be located on the failing iterate (_solve_tpi hands back the seed)
    first = rules.copy()
    c1 = dict(cal, tpi_on=True, tpi_cap_bp=cap0)
    RE._tpi_form(first, c1, ss, gz=False)
    ok1, _, w1 = RE._stage(first, c1, ss, sp, (0, 1), False, f"first rung cap={cap0:.0f}bp", True, warm=0)
    print(f"[rung 1] ok={ok1} worst={w1:.2e} ({time.perf_counter() - t0:.0f}s)", flush=True)
    if not ok1:
        worst_nodes(first, c1, ss, sp, PM, SOLVE, STATE_NAMES, np)
        print(f"[END] {tag} ({time.perf_counter() - t0:.0f}s)", flush=True)
        return
    c_on = dict(cal, tpi_on=True, tpi_cap_bp=200.0)
    tpi, ok, w = RE._solve_tpi(rules.copy(), c_on, ss, sp, "auto", True)
    print(f"[tpi] ok={ok} worst={w:.2e} last cap solved {tpi.tpi_cap_solved} ({time.perf_counter() - t0:.0f}s)", flush=True)
    if not ok:
        worst_nodes(tpi, dict(c_on, tpi_cap_bp=tpi.tpi_cap_solved or cap0), ss, sp, PM, SOLVE, STATE_NAMES, np)
    pickle.dump(dict(rules=tpi.copy(), ok=ok, worst=w, cap=tpi.tpi_cap_solved),
                open(os.path.join(HERE, "runs", f"{tag}_tpi.pkl"), "wb"))
    if tpi.tpi_cap_solved is not None:
        c_on["tpi_cap_bp"] = tpi.tpi_cap_solved
        S0 = RE.stochastic_rest_point(tpi, c_on, ss, sp, verbose=False)
        for pd in (0.01, 0.02):
            S = S0.copy(); S[IS] = RE.s_from_pd(pd)
            walls = 0
            for t in range(25):
                o = RE.read_exact(tpi, c_on, ss, sp, S.copy())
                if t == 0:
                    m0 = o["m_cb"]
                Sn = RE.advance(o, S, sp); S = tpi.grid.clip(Sn)[0]
                walls += bool(np.any(np.abs(Sn - S) > 1e-12))
            print(f"[read] cap {c_on['tpi_cap_bp']:.0f} bp, shock p^d {100 * pd:.1f}%: impact purchase {100 * m0 / B:.1f}% of B, "
                  f"{walls}/25 quarters at the box wall", flush=True)
    print(f"[END] {tag} ({time.perf_counter() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
