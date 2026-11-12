# THE GATE FOR THE COMPLEX-STEP JACOBIAN (2026-10-28): AT A STALLED JOINT CHECKPOINT, HOW WRONG WAS THE
# FORWARD-DIFFERENCE JACOBIAN (COLUMN ERRORS, WHERE), AND DOES parsolve CONVERGE WITH THE EXACT ONE?
#   python3 -u probe_csj.py <checkpoint> [maxit]
import os, pickle, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
import numpy as np


def main():
    from global_projection.solver_recursive import collocation as C
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES
    tag = sys.argv[1]
    maxit = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    sv = pickle.load(open(HERE / "runs" / tag, "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    sp = s_process_params(cal)
    args = (rules, cal, ss, sp, (0, 1), False, rules.n_gh, None, ("normal",))
    F = C.make_residual(*args)
    x = C.pack(rules, (0, 1), ("normal",))
    f = F(x)
    m = x.size
    print(f"[{tag}] start: max|F| {np.max(np.abs(f)):.3e}  sum F^2 {f @ f:.3e}  (accept: sum <= {m * 1e-18:.1e})",
          flush=True)
    n = os.cpu_count()
    pools = {k: C._PoolJacobian(args, n, k) for k in ("fd", "complex")}
    try:
        J = {}
        for k, pool in pools.items():
            t = time.perf_counter()
            J[k] = pool(F, x, f, 1e-6)
            print(f"  {k:7s} Jacobian: {time.perf_counter() - t:6.1f}s  cond {np.linalg.cond(J[k]):.2e}", flush=True)
        # where the forward differences were wrong: per column, relative to the exact column
        err = np.abs(J["fd"] - J["complex"]).max(axis=0) / np.maximum(np.abs(J["complex"]).max(axis=0), 1e-8)
        nrule, npt = len(C.STORE_RULES), rules.grid.n
        print(f"  FD column error vs exact: median {np.median(err):.1e}, 90% {np.quantile(err, 0.9):.1e}, "
              f"max {err.max():.1e}; columns off by > 1%: {(err > 1e-2).sum()} of {m}", flush=True)
        for c in np.argsort(err)[::-1][:8]:
            k, rem = divmod(c, 2 * npt)
            d, i = divmod(rem, npt)
            u = rules.grid.to_unit(rules.grid.points[i:i + 1])[0]
            edge = ", ".join(f"{s} {v:+.0f}" for s, v in zip(STATE_NAMES, u) if abs(v) > 0.99)
            print(f"    col {c:4d}: {C.STORE_RULES[k]:7s} regime {d} node {i:2d} ({edge or 'interior'})  "
                  f"error {err[c]:.1e}", flush=True)
        # the Newton direction each Jacobian gives, and how far it gets in one damped step
        for k in J:
            dx = np.linalg.lstsq(J[k], f, rcond=None)[0]
            best = min((float(F(x - s * dx) @ F(x - s * dx)), s) for s in (1.0, 0.5, 0.25, 0.125, 1 / 16, 1 / 64))
            print(f"  one step with {k:7s}: best sum F^2 {best[0]:.3e} at step {best[1]:.3g}", flush=True)
        # the full solver with the exact Jacobian
        t = time.perf_counter()
        xs, ok, its, worst = C.parsolve(F, x, tol=m * 1e-18, maxcount=maxit, label=" complex",
                                        floor_tol=1e-8, jacobian=pools["complex"])
        print(f"[{tag}] parsolve with the complex-step Jacobian: ok {ok}  its {its}  max|F| {worst:.2e}  "
              f"({time.perf_counter() - t:.0f}s)", flush=True)
        C.write_back(xs, rules, (0, 1), ("normal",))
        with open(HERE / "runs" / f"csj_{tag}", "wb") as fh:
            pickle.dump(dict(sv, rules=rules), fh)
    finally:
        for pool in pools.values():
            pool.close()


if __name__ == "__main__":
    main()
