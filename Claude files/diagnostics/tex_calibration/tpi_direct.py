# THE TPI SOLVED DIRECTLY AT ITS TARGET CAP (2026-11-03). THE CAP WALK IS SLOW ONLY WHERE MANY NODES SWITCH TO BUYING
# WITHIN A FEW bp (THE 23-NODE p^d-12% SLICE OF THE FINE CROSS GRID, ~1150-950 bp); THE JACOBIAN AT THE SOLVED RUNGS IS
# WELL POSED (cond ~4e6 AT 1075 AND AT 1002.5 bp; NO FOLD). AT THE TARGET EVERY NODE IS CLEARLY BUYING OR CLEARLY NOT, SO
# HERE: FROM A SOLVED RUNG, EVERY NODE WHOSE PRICE IS BELOW THE TARGET FLOOR IS SET TO ITS CORNER (THE D BANKS SELL THEIR
# BOOK) AND SOLVED ALONE, THEN THE GLOBAL NEWTON RUNS AT THE TARGET.
#   python3 -u tpi_direct.py <rung pickle> [--cap 200] [--maxit 30] [--tag NAME]
import json, os, pickle, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np


def main():
    from tpi_floor_sweep import economy, book
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals, _tpi_floor
    args = sys.argv[2:]
    opt = lambda k, d: type(d)(args[args.index(k) + 1]) if k in args else d
    cap, maxit = opt("--cap", 200.0), opt("--maxit", 30)
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    tag = opt("--tag", f"direct{cap:.0f}_floor{100 * sv['floor']:.0f}")
    r = sv["rules"].copy()
    cal, ss, sp = economy(sv["floor"])
    c = dict(cal, tpi_on=True, tpi_cap_bp=float(sv["cap"]))
    hard = RE._soften(c)
    c["tpi_eps"] = float(sv["tpi_eps"])
    RE._tpi_form(r, c, ss, gz=False)
    c["tpi_cap_bp"] = cap
    g, B = r.grid, cal["B_gov_D_ss"]
    t0 = time.perf_counter()
    # 1. every node below the target floor to its corner, solved alone (the other nodes held)
    moved = []
    for i in range(g.n):
        x = np.array([r.vals[k][0][i] for k in SOLVE])
        _, o = point_residuals(g.points[i], 0, x, r, c, ss, sp, RE.N_GH)
        if o["Q_bD"] >= o["Q_floor"]:
            continue
        z_now = RE._node_values(r, i, (0, 1))
        for k, val in (("x_cb", (o["b_D_D_new"] + o["m_cb"]) / B), ("Q_bD", o["Q_floor"])):
            v = r.vals[k][0].copy(); v[i] = val; r.set_values(k, 0, v)
        w = RE._node_solve(r, c, ss, sp, i, (0, 1), False, RE._node_values(r, i, (0, 1)))
        x = np.array([r.vals[k][0][i] for k in SOLVE])
        _, o1 = point_residuals(g.points[i], 0, x, r, c, ss, sp, RE.N_GH)
        if not np.isfinite(w) or w > 1e-8:
            RE._node_residual(r, c, ss, sp, i, (0, 1), False, z_now)
        moved.append(dict(node=i, alone=float(w), m=float(o1["m_cb"] / B), b_DD=float(o1["b_D_D_new"] / B)))
        print(f"  [direct] node {i:2d}: alone max|res| {w:.1e}, bought {100 * o1['m_cb'] / B:.1f}% of D debt, D banks "
              f"hold {100 * o1['b_D_D_new'] / B:.1f}%", flush=True)
    print(f"[direct] {len(moved)} nodes set to the target floor ({time.perf_counter() - t0:.0f}s)", flush=True)
    # 2. the global Newton at the target
    ok, _, w = RE._stage(r, c, ss, sp, (0, 1), False, f"TPI cap={cap:.1f}bp (direct)", True, warm=0, maxit=maxit)
    rec = dict(cap=cap, floor=sv["floor"], start=sv["cap"], nodes=moved, newton_ok=bool(ok), newton_worst=float(w),
               seconds=time.perf_counter() - t0)
    if ok:
        gz, cg = r.copy(), dict(c)
        gz.tpi_cap_solved = cap
        RE._tpi_form(gz, cg, ss, gz=True)
        rec["book"] = book(RE, gz, cg, ss, sp)
        with open(HERE / "runs" / f"tpi_gs_{tag}_{cap:.0f}.pkl", "wb") as fh:
            pickle.dump(dict(rules=r, floor=sv["floor"], cap=cap, tpi_eps=c["tpi_eps"], tag=tag,
                             smoothing={k: c[k] for k in hard}), fh)
    print(f"[direct] cap {cap:.1f} bp: Newton ok {ok} ({w:.1e}); footprint {rec.get('book')} ({rec['seconds']:.0f}s)",
          flush=True)
    (HERE / "runs" / f"tpi_{tag}.json").write_text(json.dumps(rec, indent=1, default=lambda z: z.item()))


if __name__ == "__main__":
    main()
