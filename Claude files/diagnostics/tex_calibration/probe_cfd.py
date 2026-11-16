# T1.5 PROBE (2026-10-26): IS THE d0 PRECISION FLOOR NEAR APPENDIX B (6e-8 AT ITS HOME BIAS) THE FORWARD-
# DIFFERENCE JACOBIAN? A FEW NEWTON STEPS FROM THE SAVED ENDPOINT WITH CENTRAL DIFFERENCES (step h) vs FORWARD.
#   python3 -u probe_cfd.py <checkpoint> [h]
import os, pickle, sys
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np

_W = {}


def _init(args):
    # BUILD THE RESIDUAL ONCE PER WORKER.
    from global_projection.solver_recursive import collocation as C
    _W["F"] = C.make_residual(*args)


def _central(x, h, cols):
    # CENTRAL-DIFFERENCE COLUMNS.
    F = _W["F"]
    out = np.empty((F(x).size, len(cols)))
    for j, i in enumerate(cols):
        xp, xm = x.copy(), x.copy()
        xp[i] += h; xm[i] -= h
        out[:, j] = (F(xp) - F(xm)) / (2 * h)
    return cols, out


def main():
    from global_projection.solver_recursive import collocation as C
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    h = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-5
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    from global_projection.solver_recursive.state_grid import s_process_params
    sp = s_process_params(cal)
    args = (rules, cal, ss, sp, (0,), True, rules.n_gh, None, ("normal",))
    F = C.make_residual(*args)
    x = C.pack(rules, (0,), ("normal",))
    f = F(x)
    print(f"start: max|F| {np.max(np.abs(f)):.3e}  sum F^2 {f @ f:.3e}", flush=True)
    n = os.cpu_count()
    with ProcessPoolExecutor(max_workers=n, mp_context=mp.get_context("spawn"), initializer=_init,
                             initargs=(args,)) as pool:
        for k in range(6):
            chunks = [c for c in np.array_split(np.arange(x.size), 4 * n) if c.size]
            J = np.empty((f.size, x.size))
            for cols, blk in pool.map(_central, [x] * len(chunks), [h] * len(chunks), chunks):
                J[:, cols] = blk
            dx = np.linalg.lstsq(J, f, rcond=None)[0]
            step, best = 1.0, None
            while step > 1e-4:
                fn = F(x - step * dx)
                if fn @ fn < f @ f:
                    best = (x - step * dx, fn)
                    break
                step *= 0.5
            if best is None:
                print(f"  central step {k + 1}: no descent (cond {np.linalg.cond(J):.2e})", flush=True)
                break
            x, f = best
            print(f"  central step {k + 1}: max|F| {np.max(np.abs(f)):.3e}  sum F^2 {f @ f:.3e}  "
                  f"step {step:.3g}  cond {np.linalg.cond(J):.2e}", flush=True)
            if np.max(np.abs(f)) < 1e-10:
                break


if __name__ == "__main__":
    main()
