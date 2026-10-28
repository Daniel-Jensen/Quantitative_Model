# BITWISE HARNESS (2026-10-28) FOR THE COMPLEX-STEP JACOBIAN EDIT: THE REAL RESIDUAL PATH MUST NOT MOVE.
# RECORDS F(theta) AND EVERY point_map OUTPUT AT THE STALLED k_band-4% JOINT CHECKPOINT, AT PERTURBED theta,
# ACROSS THE KINK SETTINGS (SMOOTHED / HARD IC, BOOKS, FLOOR), NO-DEFAULT AND JOINT, PLUS 20 FD COLUMNS.
#   python3 -u csj_bitwise.py record|compare
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
import numpy as np

OUT = HERE / "runs" / "csj_bitwise_baseline.pkl"


def record():
    # EVERY NUMBER THE REAL PATH PRODUCES ON THE CASES BELOW.
    from global_projection.solver_recursive import collocation as C
    from global_projection.solver_recursive.point_map import point_residuals
    from global_projection.solver_recursive.state_grid import s_process_params
    sv = pickle.load(open(HERE / "runs" / "m1_k4_joint_eps1e-2.pkl", "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    sp = s_process_params(cal)
    rng = np.random.default_rng(7)
    out = {}
    for kinks in ((1e-2, 1e-4, 1e-3), (0.0, 0.0, 0.0)):
        c = dict(cal, mu_eps=kinks[0], book_eps=kinks[1], floor_eps=kinks[2])
        for regimes, nd in (((0, 1), False), ((0,), True)):
            F = C.make_residual(rules, c, ss, sp, regimes, nd, rules.n_gh, None, ("normal",))
            xr = C.pack(rules, regimes, ("normal",))
            for k, scale in enumerate((0.0, 1e-3, 3e-2)):
                x = xr + scale * rng.standard_normal(xr.size)
                out[(kinks, regimes, k)] = F(x)
            f = F(xr)
            cols = rng.choice(xr.size, 20, replace=False)
            out[(kinks, regimes, "fd")] = np.column_stack(
                [(F(np.where(np.arange(xr.size) == i, xr + 1e-6, xr)) - f) / 1e-6 for i in cols])
            C.write_back(xr, rules, regimes, ("normal",))
            for i in (0, 7, 24):
                for d in regimes:
                    xx = np.array([rules.vals[k][d][i] for k in C.SOLVE])
                    r, o = point_residuals(rules.grid.points[i], d, xx, rules, c, ss, sp, rules.n_gh, nd)
                    out[(kinks, regimes, "pt", i, d)] = (r, {k: np.asarray(v) for k, v in o.items()
                                                            if not isinstance(v, tuple)})
    return out


def same(a, b):
    # BITWISE EQUALITY, NaN-AWARE, OVER NESTED RESULTS.
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, tuple):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    a, b = np.asarray(a), np.asarray(b)
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


if __name__ == "__main__":
    if sys.argv[1] == "record":
        pickle.dump(record(), open(OUT, "wb"))
        print(f"recorded {OUT.name}")
    else:
        base, now = pickle.load(open(OUT, "rb")), record()
        bad = [k for k in base if not same(base[k], now[k])]
        print(f"{len(base)} cases, {len(bad)} differ" + ("" if not bad else f": {bad[:6]}"))
        sys.exit(1 if bad else 0)
