# IS THE JOINT STALL A BRANCH CHOICE? (2026-10-28) AT THE K_F EDGE NODES THE STALLED SOLVE SITS ON A FUNDED-FLOOR
# DOOM LOOP (Q_bD 0.75-0.79, recap_D > n_ss, D DEBT x3, alpha_D AT ITS CAP). RE-SEED THOSE NODES FROM THE CENTRE
# NODE (CAPITAL SCALED TO THE NODE'S OWN STATE) AND SOLVE THE JOINT STAGE WITH THE CHOSEN JACOBIAN.
#   python3 -u probe_edge_reseed.py <checkpoint> [fd|complex] [maxit]
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np


def main():
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import s_process_params, IK_D, IK_F
    from global_projection.solver_recursive.collocation import STORE_RULES
    from m1_coarse import reads
    tag, mode = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else "complex")
    maxit = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    sv = pickle.load(open(HERE / "runs" / tag, "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    cal["jacobian"] = mode
    sp = s_process_params(cal)
    g = rules.grid
    u = g.to_unit(g.points)
    centre = int(np.argmin(np.abs(u).sum(axis=1)))
    edge = [i for i in range(g.n) if abs(u[i, IK_F]) > 0.99]
    print(f"[reseed] centre node {centre}, K_F edge nodes {edge}, Jacobian {mode}", flush=True)
    for d in (0, 1):
        for k in STORE_RULES:
            v = rules.vals[k][d].copy()
            for i in edge:
                v[i] = v[centre]
                if k == "Kp_D":
                    v[i] *= g.points[i, IK_D] / g.points[centre, IK_D]
                if k == "Kp_F":
                    v[i] *= g.points[i, IK_F] / g.points[centre, IK_F]
            rules.set_values(k, d, v)
    ok, its, w = RE._stage(rules, cal, ss, sp, (0, 1), False, f"joint reseeded ({mode})", True, warm=0,
                           maxit=maxit, sets=("normal",))
    print(f"[reseed] ok {ok}  its {its}  max|F| {w:.2e}", flush=True)
    with open(HERE / "runs" / f"reseed_{mode}_{tag}", "wb") as fh:
        pickle.dump(dict(sv, rules=rules, cal=cal), fh)
    for d in (0,):
        print("  edge nodes after: " + "; ".join(
            f"n{i}: Q_bD {rules.vals['Q_bD'][d][i]:.3f} alpha_D {rules.vals['alpha_D'][d][i]:.2f}" for i in edge))
    reads(RE, rules, cal, ss, sp, f"reseeded {mode}")


if __name__ == "__main__":
    main()
