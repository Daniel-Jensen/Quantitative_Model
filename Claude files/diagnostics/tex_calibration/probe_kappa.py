# PROBE (2026-10-27): IS THE d0 PRECISION FLOOR THE CROSS-BORDER DEPOSIT POSITION'S NEAR-UNIT ROOT? THE SAME STALLED
# ECONOMY (CHECKPOINT) RE-SOLVED AT kappa_nfa k, AND THE NEXT G RUNG AT k.   python3 -u probe_kappa.py <ckpt> <k> [G frac]
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
from diag_homotopy import OLD, economy
from diag_path import relative_seed, FLOOR


def main():
    from calibration.global_projection import get_calibration
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive import recursive_experiment as RE
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    k = float(sys.argv[2])
    rules, cal, ss, vals = sv["rules"], dict(sv["cal"]), sv["ss"], dict(sv["vals"])
    cal["kappa_nfa"] = k
    sp = s_process_params(cal)
    ok, it, w = RE._stage(rules, cal, ss, sp, (0,), True, f"same economy kappa {k}", False, warm=0, maxit=30)
    print(f"[same economy, kappa_nfa {k}] ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
    if len(sys.argv) > 3:
        base = get_calibration(); base.update(nw_floor_frac=FLOOR, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3,
                                               kappa_nfa=k)
        g = float(sys.argv[3])
        v = dict(vals, **{c: OLD[c] + g * (float(base[c]) - OLD[c]) for c in ("G_D", "G_F")})
        cal_k, ss_k, sp_k = economy(base, v)
        trial = relative_seed(rules, cal, ss, cal_k, ss_k)
        ok, it, w = RE._stage(trial, cal_k, ss_k, sp_k, (0,), True, f"G {g} kappa {k}", False, warm=0, maxit=30)
        print(f"[G {g}, kappa_nfa {k}] ok {ok}  its {it}  max|F| {w:.2e}", flush=True)


if __name__ == "__main__":
    main()
