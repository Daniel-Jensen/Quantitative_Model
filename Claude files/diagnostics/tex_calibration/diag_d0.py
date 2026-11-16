# D0 DIAGNOSTIC (2026-10-26): THE NO-DEFAULT STAGE AT pi = 0 ON THE APPENDIX-B CALIBRATION (THE N2 SOLVE),
# RUNG BY RUNG, WITH EACH NODE'S DISTANCE TO THE IC KINK AND WHERE THE RESIDUAL SITS.
#   python3 -u "Claude files/diagnostics/tex_calibration/diag_d0.py" > runs/diag_d0.log
import os
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

import numpy as np


def kink_table(rules, cal, ss, sproc, label):
    # z = 1 - E[Om] R n / (lambda A) FOR BOTH BANKS AT EVERY NODE (d = 0, pi = 0), WITH mu.
    import global_projection.solver_recursive.point_map as pm
    from global_projection.solver_recursive.decision_rules import SOLVE
    zs = []
    for i, S in enumerate(rules.grid.points):
        x = np.array([rules.vals[k][0][i] for k in SOLVE])
        v = pm._read_state_and_regime(S, 0, x, rules, cal, ss)
        pm._production(v, cal); pm._government(v, cal, ss); pm._balance_sheets(v, rules, cal)
        pm._continuation(v, rules, cal, sproc, rules.n_gh, True)
        pm._discount_kernels(v, rules, cal); pm._incentive_constraint(v, cal)
        zD = 1.0 - v.E_Om_D * (1.0 + v.rdep_D) * v.n_D / v.lev_D
        zF = 1.0 - v.E_Om_F * (1.0 + v.rdep_F) * v.n_F / v.lev_F
        zs.append((zD, zF, v.mu_D, v.mu_F))
    zs = np.array(zs)
    u = rules.grid.points_unit
    print(f"\n[{label}] IC distance z at the {len(zs)} nodes (z > 0 binds; mu = max(z, 0) when hard)")
    for i in np.argsort(np.abs(zs[:, 0]))[:8]:
        ax = [f"{a}{'+' if u[i, j] > 0 else '-'}" for j, a in enumerate(STATE) if abs(u[i, j]) > 1e-9]
        print(f"   node {i:2d} {','.join(ax) or 'centre':18s} z_D {zs[i, 0]:+.5f}  mu_D {zs[i, 2]:.5f}  "
              f"z_F {zs[i, 1]:+.5f}  mu_F {zs[i, 3]:.5f}")
    print(f"   z_D range [{zs[:, 0].min():+.4f}, {zs[:, 0].max():+.4f}];  nodes slack (z_D<0): "
          f"{int(np.sum(zs[:, 0] < 0))}; within 1e-3 of the kink: {int(np.sum(np.abs(zs[:, 0]) < 1e-3))}")
    return zs


def worst_residuals(rules, cal, ss, sproc, label, k=8):
    # THE WORST (NODE, EQUATION) RESIDUALS OF THE d = 0 SYSTEM AT pi = 0.
    from global_projection.solver_recursive import collocation as C
    F = C.make_residual(rules, cal, ss, sproc, regimes=(0,), no_default=True, n_gh=rules.n_gh,
                        sets=("normal",))
    r = F(C.pack(rules, (0,), ("normal",))).reshape(rules.grid.n, C.N_RES)
    u = rules.grid.points_unit
    print(f"[{label}] max|F| {np.max(np.abs(r)):.2e}; worst (node, equation):")
    for flat in np.argsort(-np.abs(r).ravel())[:k]:
        i, q = divmod(flat, C.N_RES)
        ax = [f"{a}{'+' if u[i, j] > 0 else '-'}" for j, a in enumerate(STATE) if abs(u[i, j]) > 1e-9]
        print(f"   node {i:2d} {','.join(ax) or 'centre':18s} {C.RES_NAMES[q]:12s} {r[i, q]:+.2e}")


def main():
    global STATE
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive import recursive_experiment as RE
    STATE = STATE_NAMES

    cal = get_calibration()
    cal["nw_floor_frac"] = 0.15
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    rules = RuleSet.from_ss(RE.normal_box(ss, cal, mu=1), ss, cal, n_regimes=2)
    rules.n_gh = RE.N_GH
    RE._soften(cal)
    cal["mu_eps"] = RE.MU_EPS_START
    print(f"working smoothings: mu_eps {cal['mu_eps']}  book_eps {cal['book_eps']}  floor_eps {cal['floor_eps']}")
    ok, it, w = RE._stage_wc(rules, cal, ss, sp, (0,), True, "d0", True)
    print(f"\n_stage_wc at mu_eps 1e-2: ok {ok}  its {it}  max|F| {w:.2e}")
    worst_residuals(rules, cal, ss, sp, "eps 1e-2")
    kink_table(rules, cal, ss, sp, "eps 1e-2")
    with open(HERE / "runs" / "diag_d0_eps1e-2.pkl", "wb") as fh:
        pickle.dump(dict(rules=rules, cal=cal, ss=ss), fh)

    # straight rungs from the 1e-2 solution, each on its own copy
    base = rules.copy()
    for eps in (7e-3, 3e-3, 1e-3, 3e-4, 0.0):
        r = base.copy()
        cal["mu_eps"] = eps
        calibrate_household_anchors(cal, ss, sp)
        ok, it, w = RE._stage(r, cal, ss, sp, (0,), True, f"d0 eps={eps:.0e}", True, warm=0, maxit=20)
        print(f"\n== straight from 1e-2 to mu_eps {eps:.0e}: ok {ok}  its {it}  max|F| {w:.2e}")
        worst_residuals(r, cal, ss, sp, f"eps {eps:.0e}")
        kink_table(r, cal, ss, sp, f"eps {eps:.0e}")


if __name__ == "__main__":
    main()
