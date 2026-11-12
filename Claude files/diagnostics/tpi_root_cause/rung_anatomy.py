# WHAT THE TPI WALK'S SAVED RUNGS LOOK LIKE AT THE BUYING NODES: PURCHASES, THE D BOOK, AND WHERE NEXT PERIOD'S STATE LANDS
# IN THE BOX. READ-ONLY ON results/GLOBAL_tpi_walk; rebuilds run.py's step-1 economy (PD_BOX_TOP 0.062, floor 0.5).
import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from calibration.global_projection import get_calibration
from global_projection.steady_state import solve_steady_state
from global_projection.solver_recursive.state_grid import s_process_params, default_prob, IS, IM, IO, IBDD
from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
from global_projection.solver_recursive.recursive_experiment import _soften, _sov_bp
from global_projection.solver_recursive.decision_rules import SOLVE
from global_projection.solver_recursive import point_map as PM


def economy():
    # RUN.PY'S STEP 1, THEN THE TPI WALK'S SETTINGS
    cal = get_calibration()
    cal.update(nw_floor_frac=0.5, jacobian="complex", pd_box_top=0.062)
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    cal = dict(cal, tpi_on=True)
    _soften(cal)
    cal["tpi_gz"] = False
    return cal, ss, sproc


def anatomy(path, cal, ss, sproc, show=True):
    sv = pickle.load(open(path, "rb"))
    rules, cap = sv["rules"], sv["cap"]
    cal["tpi_cap_bp"] = cap
    g, B = rules.grid, cal["B_gov_D_ss"]
    M_half = 0.32 * (cal["B_gov_D_ss"] - cal["b_D_F_ss"])
    eps, wq = PM.gh_nodes(rules.n_gh or 5)
    rows, F = [], []
    for i in range(g.n):
        for d in (0, 1):
            x = np.array([rules.vals[k][d][i] for k in SOLVE])
            r, o = PM.point_residuals(g.points[i], d, x, rules, cal, ss, sproc, n_gh=rules.n_gh or 5)
            F.append(np.max(np.abs(r)))
            if d:
                continue
            v = PM._read_state_and_regime(g.points[i], 0, x, rules, cal, ss)
            rows.append(dict(i=i, pd=100 * default_prob(g.points[i][IS]), sov=_sov_bp(o, cal), m=o["m_cb"] / B,
                             M=o["M_cb_new"] / B, Mrel=o["M_cb_new"] / M_half, bDD=o["b_D_D_new"] / B,
                             bDF=o["b_D_F_new"] / B, foc=o["foc_D"], mu=o["mu_D"], gap=o["Q_bD"] - o["Q_floor"],
                             Mlag=g.points[i][IM] / M_half, res=np.max(np.abs(r)), worst=int(np.argmax(np.abs(r)))))
    buy = [w for w in rows if w["m"] > 1e-4]
    if show:
        print(f"\n=== rung {cap:.2f} bp: {len(buy)} nodes buying; max node residual {max(F):.1e} ===")
        print("  node  p^d%   sov_bp  buy%B  M'%B  M'/M_half  b_DD'%B  b_DF'%B   foc_D     mu_D   M_lag/M_half  max|r| (eq)")
        for w in sorted(rows, key=lambda w: -w["res"])[:6] + buy:
            print(f"  {w['i']:4d} {w['pd']:5.2f} {w['sov']:8.1f} {100 * w['m']:6.2f} {100 * w['M']:6.2f} {w['Mrel']:9.2f}"
                  f" {100 * w['bDD']:8.2f} {100 * w['bDF']:7.2f} {w['foc']:9.2e} {w['mu']:7.4f} {w['Mlag']:10.2f}"
                  f"     {w['res']:.1e} ({w['worst'] + 1})")
    return rules, rows


if __name__ == "__main__":
    cal, ss, sproc = economy()
    print(f"B_gov_D_ss {cal['B_gov_D_ss']:.4f}, b_DD_ss {cal['B_gov_D_ss'] - cal['b_D_F_ss']:.4f}, "
          f"M box half-width {0.32 * (cal['B_gov_D_ss'] - cal['b_D_F_ss']):.4f} "
          f"({100 * 0.32 * (cal['B_gov_D_ss'] - cal['b_D_F_ss']) / cal['B_gov_D_ss']:.1f}% of B)")
    for p in sorted((ROOT / "results" / "GLOBAL_tpi_walk").glob("tpi_rung_*.pkl"), reverse=True):
        anatomy(p, cal, ss, sproc)
