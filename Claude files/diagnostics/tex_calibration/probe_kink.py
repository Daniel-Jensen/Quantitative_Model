# WHY THE EXACT NEWTON STEP DOES NOT DESCEND (2026-10-28): THE LINE PROFILE ALONG THE EXACT NEWTON DIRECTION,
# THE LINEARISATION ERROR AT EACH STEP LENGTH (O(t) IF THE MAP IS SMOOTH THERE, FLAT IF A KINK IS CROSSED), AND
# THE (REGIME, NODE, EQUATION) CELLS WHERE THE NONLINEARITY SITS.
#   python3 -u probe_kink.py <checkpoint> [d0]   (d0: the no-default regime alone, as the d0 stage solves it)
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
import numpy as np


def main():
    from global_projection.solver_recursive import collocation as C
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES
    tag = sys.argv[1]
    sv = pickle.load(open(HERE / "runs" / tag, "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    sp = s_process_params(cal)
    regs, nd = ((0,), True) if sys.argv[2:] == ["d0"] else ((0, 1), False)
    args = (rules, cal, ss, sp, regs, nd, rules.n_gh, None, ("normal",))
    F = C.make_residual(*args)
    x = C.pack(rules, regs, ("normal",))
    f = F(x)
    npt, nres = rules.grid.n, C.N_RES
    pool = C._PoolJacobian(args, os.cpu_count(), "complex")
    try:
        J = pool(F, x, f, 1e-6)
    finally:
        pool.close()
    U, S, Vt = np.linalg.svd(J)
    print(f"[{tag}] sum F^2 {f @ f:.3e}; singular values: max {S[0]:.2e}, then smallest "
          + ", ".join(f"{s:.1e}" for s in S[-5:]), flush=True)
    # the residual's weight on the weakest directions: Newton divides these by the singular values
    uf = U.T @ f
    print("  |u_k . f| on the 5 weakest: " + ", ".join(f"{abs(v):.1e}" for v in uf[-5:]), flush=True)
    names = [(k, d, i) for k in C.STORE_RULES for d in regs for i in range(npt)]

    def where(v, n=6):
        # THE LARGEST ENTRIES OF AN UNKNOWN-SPACE VECTOR, AS (RULE, REGIME, NODE AT THE BOX EDGE).
        out = []
        for c in np.argsort(np.abs(v))[::-1][:n]:
            k, d, i = names[c]
            u = rules.grid.to_unit(rules.grid.points[i:i + 1])[0]
            edge = ",".join(f"{s}{v_:+.0f}" for s, v_ in zip(STATE_NAMES, u) if abs(v_) > 0.99) or "in"
            out.append(f"{k}/d{d}/n{i}({edge}) {v[c]:+.2e}")
        return "; ".join(out)
    print("  weakest direction v_min: " + where(Vt[-1]), flush=True)
    dx = Vt.T @ (uf / S)
    print(f"  Newton step |dx|_inf {np.max(np.abs(dx)):.2e}: " + where(dx), flush=True)
    lin = J @ dx
    print(f"  {'t':>9s} {'sum F^2':>10s} {'lin err/t':>10s}  worst cell of the nonlinearity", flush=True)
    for k in range(0, 31, 2):
        t = 2.0 ** -k
        ft = F(x - t * dx)
        nl = ft - (f - t * lin)
        c = int(np.argmax(np.abs(nl)))
        g, e = divmod(c, nres)
        d, i = divmod(g, npt)
        print(f"  {t:9.2e} {ft @ ft:10.3e} {np.max(np.abs(nl)) / t:10.2e}  d{d} node {i} {C.RES_NAMES[e]}", flush=True)


if __name__ == "__main__":
    main()
