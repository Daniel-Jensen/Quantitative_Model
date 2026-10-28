# STEP 1 OF THE 2026-10-28 PLAN: PICK THE CALM ANSWER AT EVERY GRID POINT, THEN RE-SOLVE. A GRID POINT CAN SIT ON A
# SELF-FULFILLING CRISIS BRANCH (K_F+ ON 2026-10-28: n_D 11% OF n_ss, alpha_D 8 IN NORMAL TIMES) THAT THE GLOBAL
# NEWTON CANNOT LEAVE. EACH POINT IS SOLVED ALONE, THE OTHERS HELD, FROM THE CENTRE POINT'S VALUES; THE ROOT IS
# ADOPTED WHEN IT CONVERGES AND GIVES THE D BANK MORE EQUITY THAN THE CURRENT VALUES.
#   python3 -u branch_select.py <checkpoint> d0|joint <out tag>
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np
from scipy.optimize import root

NODE_TOL = 1e-10
# the crisis signature in normal times on the coarse box: D equity below half its SS level, or alpha_D above 3
CRISIS_N, CRISIS_ALPHA = 0.5, 3.0


def _write(rules, i, regs, z):
    # NODE i'S VALUES OF EVERY RULE (FIT SPACE, RULE BY RULE THEN REGIME) INTO THE RULES, REFITTED.
    from global_projection.solver_recursive.decision_rules import STORE_RULES, from_fit
    j = 0
    for k in STORE_RULES:
        for d in regs:
            v = rules.vals[k][d].copy()
            v[i] = from_fit(k, z[j])
            rules.set_values(k, d, v)
            j += 1


def _read(rules, i, regs):
    # NODE i'S VALUES IN FIT SPACE, IN _write's ORDER.
    from global_projection.solver_recursive.decision_rules import STORE_RULES, to_fit
    return np.array([to_fit(k, rules.vals[k][d][i]) for k in STORE_RULES for d in regs], dtype=float)


def node_residual(rules, i, regs, z, cal, ss, sp, nd):
    # THE 20 RESIDUALS PER REGIME AT NODE i WITH ITS VALUES SET TO z, THE OTHER NODES HELD.
    from global_projection.solver_recursive import collocation as C
    from global_projection.solver_recursive.recursive_experiment import N_GH
    _write(rules, i, regs, z)
    vals = {k: {d: rules.vals[k][d] for d in regs} for k in C.STORE_RULES}
    return np.concatenate([C._point_block(i, d, vals, rules, cal, ss, sp, N_GH, nd) for d in regs])


def node_read(rules, i, d, cal, ss, sp, nd):
    # THE PERIOD MAP'S OUTPUTS AT NODE i IN REGIME d.
    from global_projection.solver_recursive.collocation import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    from global_projection.solver_recursive.recursive_experiment import N_GH
    x = np.array([rules.vals[k][d][i] for k in SOLVE])
    return point_residuals(rules.grid.points[i], d, x, rules, cal, ss, sp, N_GH, nd)[1]


def crisis_nodes(rules, cal, ss, sp, nd):
    # THE NODES WHOSE NORMAL-TIMES SOLUTION SHOWS THE CRISIS SIGNATURE.
    nss = ss["ss_bank_D"]["n_ss"]
    out = []
    for i in range(rules.grid.n):
        o = node_read(rules, i, 0, cal, ss, sp, nd)
        if o["n_D"] < CRISIS_N * nss or o["alpha_D"] > CRISIS_ALPHA:
            out.append((i, round(float(o["n_D"] / nss), 3), round(float(o["alpha_D"]), 2)))
    return out


def calm_scan(rules, cal, ss, sp, regs, nd, verbose=True, only=None):
    # EACH NODE (or only the listed ones) SOLVED ALONE FROM THE CENTRE'S VALUES; THE ROOT ADOPTED WHEN IT CONVERGES
    # WITH MORE D EQUITY.
    from global_projection.solver_recursive.decision_rules import STORE_RULES
    from global_projection.solver_recursive.state_grid import IK_D, IK_F
    g, nss = rules.grid, ss["ss_bank_D"]["n_ss"]
    c = int(np.argmin(np.abs(g.to_unit(g.points)).sum(axis=1)))
    switched = []
    for i in (range(g.n) if only is None else only):
        z_now = _read(rules, i, regs)
        z0 = _read(rules, c, regs)
        j = 0
        for k in STORE_RULES:
            for _ in regs:
                # capital scaled to the node's own state, in the log fit space
                if k == "Kp_D":
                    z0[j] += np.log(g.points[i, IK_D] / g.points[c, IK_D])
                if k == "Kp_F":
                    z0[j] += np.log(g.points[i, IK_F] / g.points[c, IK_F])
                j += 1
        n_now = node_read(rules, i, 0, cal, ss, sp, nd)["n_D"]
        f_now = np.max(np.abs(node_residual(rules, i, regs, z_now, cal, ss, sp, nd)))
        sol = root(lambda z: node_residual(rules, i, regs, z, cal, ss, sp, nd), z0, method="hybr", tol=1e-13)
        f_new = np.max(np.abs(node_residual(rules, i, regs, sol.x, cal, ss, sp, nd)))
        n_new = node_read(rules, i, 0, cal, ss, sp, nd)["n_D"]
        if f_new <= NODE_TOL and n_new > n_now * (1.0 + 1e-3):
            switched.append(i)
            if verbose:
                print(f"  [calm] node {i:2d}: switched, n_D/n_ss {n_now / nss:.3f} -> {n_new / nss:.3f} "
                      f"(node residual {f_now:.1e} -> {f_new:.1e})", flush=True)
        else:
            _write(rules, i, regs, z_now)
    return switched


def global_max(rules, cal, ss, sp, regs, nd):
    # THE FULL COLLOCATION RESIDUAL AT THE CURRENT RULES.
    from global_projection.solver_recursive import collocation as C
    from global_projection.solver_recursive.recursive_experiment import N_GH
    F = C.make_residual(rules, cal, ss, sp, regs, nd, N_GH, None, ("normal",))
    return float(np.max(np.abs(F(C.pack(rules, regs, ("normal",))))))


def polish_sweeps(rules, cal, ss, sp, regs, nd, sweeps=4):
    # NONLINEAR GAUSS-SEIDEL: EVERY NODE RE-SOLVED ALONE FROM ITS OWN VALUES, SWEEP AFTER SWEEP, SO THE NODES AGREE
    # BEFORE THE GLOBAL NEWTON (A NODE'S CONTINUATION READS ITS NEIGHBOURS ONLY WEAKLY).
    for s in range(sweeps):
        for i in range(rules.grid.n):
            z_now = _read(rules, i, regs)
            f_now = np.max(np.abs(node_residual(rules, i, regs, z_now, cal, ss, sp, nd)))
            sol = root(lambda z: node_residual(rules, i, regs, z, cal, ss, sp, nd), z_now, method="hybr", tol=1e-13)
            if np.max(np.abs(node_residual(rules, i, regs, sol.x, cal, ss, sp, nd))) >= f_now:
                _write(rules, i, regs, z_now)
        w = global_max(rules, cal, ss, sp, regs, nd)
        print(f"  [sweep {s + 1}] global max|F| {w:.2e}; crisis nodes {crisis_nodes(rules, cal, ss, sp, nd)}", flush=True)
        if w < 1e-8:
            break


def main():
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import s_process_params
    from diag_path import FLOOR
    tag, mode, out = sys.argv[1], sys.argv[2], sys.argv[3]
    sv = pickle.load(open(HERE / "runs" / tag, "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    cal.update(nw_floor_frac=FLOOR, jacobian="complex")
    sp = s_process_params(cal)
    regs, nd = ((0,), True) if mode == "d0" else ((0, 1), False)
    bad = crisis_nodes(rules, cal, ss, sp, nd)
    print(f"[{tag}] {mode}: global max|F| {global_max(rules, cal, ss, sp, regs, nd):.2e}; crisis nodes {bad}", flush=True)
    sw = calm_scan(rules, cal, ss, sp, regs, nd, only=[b[0] for b in bad])
    print(f"[{tag}] switched {sw}; global max|F| {global_max(rules, cal, ss, sp, regs, nd):.2e}", flush=True)
    polish_sweeps(rules, cal, ss, sp, regs, nd)
    with open(HERE / "runs" / f"{out}_scan.pkl", "wb") as fh:
        pickle.dump(dict(sv, rules=rules, cal=cal), fh)
    ok, it, w = RE._stage(rules, cal, ss, sp, regs, nd, f"{mode} calm", True, warm=0, maxit=40 if nd else 60)
    print(f"[{mode} calm] ok {ok}  its {it}  max|F| {w:.2e}; crisis nodes {crisis_nodes(rules, cal, ss, sp, nd)}",
          flush=True)
    with open(HERE / "runs" / f"{out}.pkl", "wb") as fh:
        pickle.dump(dict(sv, rules=rules, cal=cal), fh)


if __name__ == "__main__":
    main()
