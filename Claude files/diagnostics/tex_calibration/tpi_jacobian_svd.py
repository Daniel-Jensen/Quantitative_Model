# IS THE FINE TPI SYSTEM NEAR A FOLD? (2026-11-03.) THE EXACT (COMPLEX-STEP) COLLOCATION JACOBIAN AT A SOLVED CAP RUNG
# (A tpi_gs_walk PICKLE, IN THE FORM AND SMOOTHINGS THE WALK SOLVES), ITS SINGULAR VALUES, AND FOR THE SMALLEST ONES THE
# UNKNOWNS (RULE, REGIME, NODE) AND RESIDUALS (EQUATION, REGIME, NODE) THEIR SINGULAR VECTORS LOAD ON.
#   python3 -u tpi_jacobian_svd.py <rung pickle> [<rung pickle> ...]
import os, pickle, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np

K_SMALL = 4
TOP = 8


def node_label(g, i):
    # NODE i's p^d AND ITS OFF-CENTRE STATES.
    from global_projection.solver_recursive.state_grid import STATE_NAMES, IS
    u = g.to_unit(g.points[i:i + 1])[0]
    off = {k: round(float(v), 2) for k, v in zip(STATE_NAMES, u) if abs(v) > 0.5 and k != "s"}
    return f"node {i:2d} (p^d {100 / (1 + np.exp(-g.points[i][IS])):.2f}%{', ' + str(off) if off else ''})"


def analyse(path):
    # THE JACOBIAN AT ONE RUNG AND ITS SMALLEST SINGULAR DIRECTIONS.
    from tpi_floor_sweep import economy
    from global_projection.solver_recursive import recursive_experiment as RE, collocation as C
    from global_projection.solver_recursive.decision_rules import STORE_RULES
    sv = pickle.load(open(HERE / "runs" / path, "rb"))
    r, cap = sv["rules"], float(sv["cap"])
    cal, ss, sp = economy(sv["floor"])
    c = dict(cal, tpi_on=True, tpi_cap_bp=cap)
    RE._soften(c)
    c.update(tpi_eps=float(sv["tpi_eps"]), **sv.get("smoothing", {}))
    # the rung pickles hold the rules in the form the walk solved them in; FB unless the walk ran with --gz
    gz = bool(np.any(r.vals["x_cb"][0] < 0.0))
    c["tpi_gz"] = gz
    args = (r, c, ss, sp, (0, 1), False, RE.N_GH, None, ("normal",))
    F = C.make_residual(*args)
    x = C.pack(r, (0, 1), ("normal",))
    f = F(x)
    t0 = time.perf_counter()
    jac = C._PoolJacobian(args, C._n_workers(c), "complex")
    J = jac(F, x, f, 1e-6)
    jac.close()
    s = np.linalg.svd(J, compute_uv=True)
    U, S, Vt = s
    g, n, nres = r.grid, r.grid.n, C.N_RES
    print(f"\n=== {path}: cap {cap:.1f} bp, tpi_eps {c['tpi_eps']:g}, form {'GZ' if gz else 'FB'}; max|F| {np.abs(f).max():.1e}; "
          f"Jacobian {J.shape} in {time.perf_counter() - t0:.0f}s", flush=True)
    print(f"  singular values: largest {S[0]:.2e}; smallest {', '.join(f'{v:.2e}' for v in S[-K_SMALL:][::-1])}; "
          f"cond {S[0] / S[-1]:.2e}", flush=True)
    for k in range(1, K_SMALL + 1):
        v, u = Vt[-k], U[:, -k]
        print(f"  direction {k} (sigma {S[-k]:.2e}):", flush=True)
        for i in np.argsort(-np.abs(v))[:TOP]:
            rule, d, node = STORE_RULES[i // (2 * n)], (i // n) % 2, i % n
            print(f"     unknown {rule:8s} regime {d} {node_label(g, node)}: {v[i]:+.3f}", flush=True)
        for j in np.argsort(-np.abs(u))[:4]:
            d, node, e = j // (n * nres), (j // nres) % n, j % nres
            print(f"     residual {C.RES_NAMES[e]:9s} regime {d} {node_label(g, node)}: {u[j]:+.3f}", flush=True)
    return S


if __name__ == "__main__":
    for p in sys.argv[1:]:
        analyse(p)
