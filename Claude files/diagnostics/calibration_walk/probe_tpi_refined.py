# THE TPI SOLVED DIRECTLY ON THE 5-NODE AXIS GRID AT SEVERAL CAPS, EACH SEEDED FROM THE LAST (2026-10-25): THE COARSE
# CAP WALK STALLS AT ITS FIRST RUNG (~5e-6 at 700 bp, either smoothing), WHILE THE FULL RUN'S REFINED GRIDS SOLVED THE
# 200 BP CAP STRAIGHT FROM THE NO-TPI SEED. Each cap: _refine_s's TPI path (FB form, smoothing 1e-2 -> 1e-3 -> tpi_eps,
# then the GZ polish), then probe_tpi_walk.reads. Box: the run's re-centred box with BOX_KW m_band = M_BAND.
# usage: python3 -u probe_tpi_refined.py TAG M_BAND CAP[,CAP...]     e.g. t5 0.32 600,400,300,200
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path[:0] = [ROOT, HERE]


def main():
    # SOLVE EACH CAP AND READ IT.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    import global_projection.solver_recursive.recursive_experiment as RE
    from probe_tpi_walk import reads
    tag, m_band, caps = sys.argv[1], float(sys.argv[2]), [float(c) for c in sys.argv[3].split(",")]
    t0 = time.perf_counter()
    base = pickle.load(open(os.path.join(ROOT, "results", "GLOBAL", "data", "rules_base_coarse.pkl"), "rb"))
    cal = get_calibration(); cal["nw_floor_frac"] = 0.15
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    B = cal["B_gov_D_ss"]
    RE.BOX_KW["m_band"] = m_band
    rules = base
    print(f"[config] {tag}: m_band {m_band}, caps {caps}, axis grid at 5 s nodes", flush=True)
    for cap in caps:
        c = dict(cal, tpi_on=True, tpi_cap_bp=cap)
        fine, ok, w = RE._refine_s(rules, c, ss, sp, rules.box, 5, "auto", True, ladder=[(5, "axis")])
        print(f"[cap {cap:.0f}] ok={ok} worst={w:.2e} ({time.perf_counter() - t0:.0f}s)", flush=True)
        print(f"[read] cap {cap:.0f}: {reads(fine.copy(), c, ss, sp, B, RE, IS, np)}", flush=True)
        pickle.dump(dict(rules=fine.copy(), cap=cap, ok=ok, worst=w, m_band=m_band),
                    open(os.path.join(HERE, "runs", f"{tag}_cap{cap:.0f}.pkl"), "wb"))
        rules = fine
    print(f"[END] {tag} ({time.perf_counter() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
