# THE TPI WALL (2026-10-30): THE COARSE CAP WALK SOLVES DOWN TO 1150 bp AND STALLS AT 1100-1125 bp (id_C_D ~1.5e-5). SOLVE THE
# COARSE TPI AT 1150 bp, THEN READ THE 1100 bp SYSTEM AT THAT SOLUTION: WHERE IT FAILS, AND WHERE THE TPI BOOK GOES.
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
import numpy as np


def main():
    import run
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES, IM, IO, IS
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive import collocation as C
    from global_projection.solver_recursive.point_map import point_residuals
    cal = get_calibration(); cal.update(nw_floor_frac=run.NW_FLOOR, jacobian=run.JACOBIAN)
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    fine = pickle.load(open(run.WARM_START, "rb"))
    coarse = RE.solve_recursive(dict(cal), ss, sp, base=RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1), 2), fine),
                                s_refine=0)
    RE.TPI_CAP_START, RE.TPI_CAP_STEP = 1150.0, 100.0
    c = dict(cal, tpi_on=True, tpi_cap_bp=1150.0)
    on = RE.solve_recursive(c, ss, sp, base=coarse, s_refine=0)
    print(f"[1150] ok {on.solve_ok} worst {on.solve_worst:.1e} cap {on.tpi_cap_solved}", flush=True)
    with open(HERE / "runs" / "tpi_1150.pkl", "wb") as fh:
        pickle.dump(on, fh)
    # the 1100 system at the 1150 solution (FB form, working smoothings, as the walk solves it)
    r = on.copy()
    c2 = dict(c, tpi_cap_bp=1100.0, mu_eps=1e-4, book_eps=1e-4, floor_eps=1e-3)
    RE._tpi_form(r, c2, ss, gz=False)
    F = C.make_residual(r, c2, ss, sp, (0, 1), False, RE.N_GH, None, ("normal",))
    f = np.abs(F(C.pack(r, (0, 1), ("normal",)))).reshape(2, r.grid.n, C.N_RES)
    g = r.grid
    print(f"[1100 at the 1150 solution] max|F| {f.max():.2e}", flush=True)
    for i in np.argsort(-f[0].max(axis=1))[:6]:
        x = np.array([r.vals[k][0][i] for k in C.SOLVE]); _, o = point_residuals(g.points[i], 0, x, r, c2, ss, sp, RE.N_GH)
        u = g.to_unit(g.points[i:i + 1])[0]
        Mn = np.array(g.points[i]); Mn[IM] = o["M_cb_new"]; Mn[IO] = o["Op_cb"]
        un = g.to_unit(Mn[None])[0]
        print(f"  node {i:2d} {({k: round(float(v), 2) for k, v in zip(STATE_NAMES, u) if abs(v) > 0.5})}: worst "
              f"{C.RES_NAMES[int(np.argmax(f[0, i]))]} {f[0, i].max():.1e}; m_cb {o['m_cb']:.3f} Q_bD {o['Q_bD']:.3f} floor "
              f"{o['Q_floor']:.3f}; next TPI book u {un[IM]:+.2f}, O_cb u {un[IO]:+.2f}", flush=True)
    print(f"  box: M_cb [{g.lo[IM]:.3f}, {g.hi[IM]:.3f}]  O_cb [{g.lo[IO]:.3f}, {g.hi[IO]:.3f}]", flush=True)


if __name__ == "__main__":
    main()
