# THE BOND-ACCURACY EXPERIMENTS (2026-10-23): FROM A SOLVED COARSE PICKLE, REFINE s ON THE "axis" OR "cross" GRID, AND
# OPTIONALLY RE-CENTRE THE BANK STOCKS OF THE BOX ON THE REST POINT AND RE-SOLVE, n TIMES. Each solve is saved as
# runs/<TAG>_it<k>.pkl for diag_fundamental.py / diag_price.py / rest_metrics.py.
# usage: python3 -u probe_fix.py TAG SEED.pkl MODE M_S RECENTRE [key=value,...]   e.g. x1 pd002.pkl axis 9 2 pd_rest=0.0002
#   a SEED written by probe_fix itself (it carries its box centre) is refined on that box straight at M_S
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))


def main():
    # SOLVE, RE-CENTRE, SAVE.
    import numpy as np
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, IS, BOX_STOCKS, _STOCK_INDEX
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.decision_rules import RuleSet
    import global_projection.solver_recursive.recursive_experiment as RE
    tag, seed, mode, m_s, n_re = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    cal_over = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in sys.argv[6].split(",")} if len(sys.argv) > 6 else {}
    t0 = time.perf_counter()
    d_ = pickle.load(open(os.path.join(HERE, "runs", seed), "rb"))
    cal = get_calibration()
    cal["nw_floor_frac"] = d_["nw_floor_frac"]
    cal.update(cal_over)
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    assert not d_.get("sp_over") or abs(d_["sp_over"]["s_star"] - sp["s_star"]) < 1e-12, "seed solved another s process"
    calibrate_household_anchors(cal, ss, sp)
    # a seed that carries its own box centre (a probe_fix solve) is refined there directly, at m_s alone
    coarse, at = d_["rules"], d_.get("at")
    resumed = at is not None
    print(f"[config] {tag}: seed {seed}, mode {mode}, m_s {m_s}, re-centre x{n_re}, cal {cal_over}", flush=True)
    for k in range(n_re + 1):
        if at is not None and not (resumed and k == 0):  # the coarse grid on the moved box, seeded from the last fine rules
            g = RE.normal_box(ss, cal, at=at)
            c = RE._seed_from(RuleSet(g, 2), fine)
            c.box = dict(fine.box, at=at)
            ok, _, w = RE._stage(c, cal, ss, sp, (0, 1), False, f"coarse it{k}", True, warm=0, maxit=40)
            print(f"[coarse it{k}] ok={ok} worst={w:.2e} ({time.perf_counter() - t0:.0f}s)", flush=True)
            coarse = c
        rules = coarse
        for m in ([] if resumed and k == 0 else [q for q in (5, 9, 17) if q < m_s]) + [m_s]:
            g = RE.normal_box(ss, cal, refine=(IS, m, mode), at=at)
            fine = RE._seed_from(RuleSet(g, 2), rules)
            fine.box = dict(coarse.box or {}, at=at)
            ok, it, w = RE._stage(fine, cal, ss, sp, (0, 1), False, f"{mode} s={m} it{k}", True, warm=0, maxit=60)
            print(f"[fine it{k} m={m}] {g.n} points ok={ok} worst={w:.2e} its={it} ({time.perf_counter() - t0:.0f}s)",
                  flush=True)
            rules = fine
        S0 = RE.stochastic_rest_point(fine, cal, ss, sp, verbose=True)
        u0 = fine.grid.to_unit(S0)[0]
        print(f"[rest it{k}] u = {np.array2string(u0, precision=3)}", flush=True)
        with open(os.path.join(HERE, "runs", f"{tag}_it{k}.pkl"), "wb") as fh:
            pickle.dump(dict(rules=fine.copy(), nw_floor_frac=cal["nw_floor_frac"], mu_eps=cal["mu_eps"],
                             cal_over=cal_over, sp_over={}, ok=bool(ok), worst=float(w), at=at), fh)
        at = {s: float(S0[_STOCK_INDEX[s]]) for s in BOX_STOCKS}
    print(f"[END] {tag} ({time.perf_counter() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
