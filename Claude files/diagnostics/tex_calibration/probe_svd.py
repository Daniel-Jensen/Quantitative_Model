# PROBE (2026-10-26): ONE SMALL HOME-BIAS RUNG FROM THE SAVED PATH SOLUTION, SOLVED VERBOSELY; AT THE STALL,
# THE SVD OF THE COLUMN-SCALED COLLOCATION JACOBIAN: THE NEAR-NULL DIRECTIONS (WHICH RULE AT WHICH NODE),
# THE RESIDUALS THEY MAP TO, AND HOW MUCH OF THE STALLED RESIDUAL LIES ALONG THEM.
#   python3 -u probe_svd.py <home-bias fraction>
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np
from diag_homotopy import OLD, economy
from diag_path import relative_seed


def main():
    frac = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0625
    from calibration.global_projection import get_calibration
    from global_projection.solver_recursive.state_grid import STATE_NAMES
    from global_projection.solver_recursive.decision_rules import STORE_RULES
    from global_projection.solver_recursive import collocation as C
    sv = pickle.load(open(HERE / "runs" / "diag_path_foreign_share.pkl", "rb"))
    base = get_calibration(); base.update(nw_floor_frac=0.15, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3)
    hb = OLD["omega_home_D"] + frac * (float(base["omega_home_D"]) - OLD["omega_home_D"])
    cal, ss, sp = economy(base, dict(sv["vals"], omega_home_D=hb))
    rules = relative_seed(sv["rules"], sv["cal"], sv["ss"], cal, ss)
    ok, it, w = C.solve_collocation(rules, cal, ss, sp, regimes=(0,), no_default=True, n_gh=rules.n_gh,
                                    maxit=30, verbose=True, label=f" home bias {frac}", sets=("normal",))
    F = C.make_residual(rules, cal, ss, sp, regimes=(0,), no_default=True, n_gh=rules.n_gh,
                        sets=("normal",))
    x = C.pack(rules, (0,), ("normal",))
    f = F(x)
    pool = C._PoolJacobian((rules, cal, ss, sp, (0,), True, rules.n_gh, None, ("normal",)), os.cpu_count())
    try:
        J = pool(F, x, f, 1e-6)
    finally:
        pool.close()
    n = rules.grid.n
    u = rules.grid.points_unit
    lab = lambda i: ",".join(f"{a}{'+' if u[i, j] > 0 else '-'}" for j, a in enumerate(STATE_NAMES)
                             if abs(u[i, j]) > 1e-9) or "centre"
    s = 1.0 / np.maximum(np.linalg.norm(J, axis=0), 1e-300)
    U, S, Vt = np.linalg.svd(J * s, full_matrices=False)
    print(f"\nscaled Jacobian {J.shape}: sigma_max {S[0]:.3e}, smallest relative "
          + ", ".join(f"{v:.2e}" for v in S[-6:] / S[0]))
    tot = float(f @ f)
    for r in (1, 2, 3):
        v, uu = Vt[-r], U[:, -r]
        share = float((uu @ f) ** 2) / tot
        print(f"\n  direction sigma/sigma_max {S[-r] / S[0]:.2e}: carries {100 * share:.1f}% of |F|^2")
        print("    unknowns:", "; ".join(f"{STORE_RULES[i // n]}@{lab(i % n)} {v[i]:+.2f}"
                                       for i in np.argsort(-np.abs(v))[:6]))
        print("    residuals:", "; ".join(f"{C.RES_NAMES[q % C.N_RES]}@{lab(q // C.N_RES)} {uu[q]:+.2f}"
                                        for q in np.argsort(-np.abs(uu))[:6]))


if __name__ == "__main__":
    main()
