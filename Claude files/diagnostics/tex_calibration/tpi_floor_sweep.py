# WHICH RECAPITALISATION FLOOR LETS THE TPI SOLVE? (2026-11-02, FIX A IN.) FOR EACH nw_floor_frac (THE FUNDED FLOOR UNDER
# D AND F BANK NET WORTH IN DEFAULT, A SHARE OF n_ss; run.NW_FLOOR 0.5): THE COARSE NO-TPI BASELINE, SEEDED FROM THE NEAREST
# SOLVED FLOOR (THE FIRST FROM THE FINE WARM START); ITS ECONOMICS BY EXACT READS; THEN THE COARSE TPI CAP WALK FROM ABOVE
# THE BOX-TOP SPREAD DOWN TO run.TPI_CAP. RECORDS THE LOWEST CAP SOLVED AND, IF THE WALK STALLS, THE WORST NODES OF THE
# NEXT RUNG AT THE LAST SOLVED ONE (STATE, EQUATION, PURCHASE, PRICE AGAINST THE FLOOR, D BOOK, NEXT TPI BOOK IN THE BOX).
# --fine: THE SAME ON THE FINE WARM GRID run.py SOLVES ON (THE COARSE 3-NODE s AXIS FOLDS ONE RUNG AFTER THE CAP BINDS).
#   python3 -u tpi_floor_sweep.py [floors, e.g. 0.5,0.3,0.15,0.7,0.9] [--fine]
import json, os, pickle, sys, time, traceback
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np

FLOORS = (0.5, 0.3, 0.15, 0.7, 0.9)
CAP_START, CAP_STEP, MIN_STEP = 1500.0, 100.0, 12.5
# Newton iterations per cap rung: a rung that has not converged by then is a stall (the 2026-10-30 walk burned 40)
RUNG_MAXIT = 15
PDS = (0.02, 0.05, 0.12)


def economy(floor):
    # run.py's CALIBRATION AT THIS FLOOR, ITS SS, s PROCESS AND ANCHORS.
    import run
    from calibration.global_projection import get_calibration
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    cal = get_calibration(); cal.update(nw_floor_frac=floor, jacobian=run.JACOBIAN)
    ss = solve_steady_state(cal, verbose=False); sp = s_process_params(cal); calibrate_household_anchors(cal, ss, sp)
    return cal, ss, sp


def reads(RE, r, cal, ss, sp):
    # THE NO-TPI ECONOMICS: THE REST POINT AND THE IMPACT AT EACH p^d (EXACT READS, COARSE GRID).
    from global_projection.solver_recursive.state_grid import IS
    S0 = RE.stochastic_rest_point(r, cal, ss, sp, verbose=False)
    o0 = RE.read_exact(r, cal, ss, sp, S0)
    sprd = lambda o: 4e4 * (o["r_wc_D"] - o["rdep_D"])
    sov = lambda o: 4e4 * (cal["delta_b_D"] * (1 / o["Q_bD"] - 1) - cal["delta_b_F"] * (1 / o["Q_bF"] - 1))
    out = dict(rest=dict(mu_D=float(o0["mu_D"]), Q_bD=float(o0["Q_bD"]), sov=sov(o0), E_Om=float(o0["E_Om_D"])))
    for pd in PDS:
        S = S0.copy(); S[IS] = RE.s_from_pd(pd)
        o = RE.read_exact(r, cal, ss, sp, S)
        out[f"{100 * pd:.0f}%"] = dict(mu_D=float(o["mu_D"]), Q_bD=float(o["Q_bD"]), sov=sov(o),
                                       spread=sprd(o) - sprd(o0), Q_floor_200=float(_floor200(o, cal)),
                                       **{f"d{k}": 100 * (o[k] / o0[k] - 1) for k in ("Q_bD", "n_D", "Y_D", "C_D")})
    return out


def _floor200(o, cal):
    # THE D BOND PRICE THE 200 bp CAP DEFENDS AT THIS STATE.
    from global_projection.solver_recursive.point_map import _tpi_floor
    return _tpi_floor(o["Q_bF"], dict(cal, tpi_cap_bp=200.0))


def stall(RE, on, cal, ss, sp, cap_next):
    # THE NEXT RUNG's SYSTEM AT THE LAST SOLVED RULES (FB FORM, WORKING SMOOTHINGS, AS THE WALK SOLVES IT): WORST NODES.
    from global_projection.solver_recursive import collocation as C
    from global_projection.solver_recursive.point_map import point_residuals
    from global_projection.solver_recursive.state_grid import STATE_NAMES, IM, IS
    r = on.copy()
    c2 = dict(cal, tpi_on=True, tpi_cap_bp=cap_next, **RE.HARD_WORK)
    RE._tpi_form(r, c2, ss, gz=False)
    F = C.make_residual(r, c2, ss, sp, (0, 1), False, RE.N_GH, None, ("normal",))
    f = np.abs(F(C.pack(r, (0, 1), ("normal",)))).reshape(2, r.grid.n, C.N_RES)
    g, rows = r.grid, []
    for d, i in sorted(((d, i) for d in (0, 1) for i in range(g.n)), key=lambda t: -f[t[0], t[1]].max())[:6]:
        x = np.array([r.vals[k][d][i] for k in C.SOLVE])
        _, o = point_residuals(g.points[i], d, x, r, c2, ss, sp, RE.N_GH)
        u = g.to_unit(g.points[i:i + 1])[0]
        nxt = np.array(g.points[i]); nxt[IM] = o["M_cb_new"]
        rows.append(dict(regime=d, node=int(i), max_F=float(f[d, i].max()), eq=C.RES_NAMES[int(np.argmax(f[d, i]))],
                         state={k: round(float(v), 2) for k, v in zip(STATE_NAMES, u) if abs(v) > 0.5},
                         pd=float(1 / (1 + np.exp(-g.points[i][IS]))), m_cb_share=float(o["m_cb"] / cal["B_gov_D_ss"]),
                         Q_bD=float(o["Q_bD"]), Q_floor=float(o["Q_floor"]),
                         b_DD_share=float(o["b_D_D_new"] / cal["B_gov_D_ss"]),
                         next_M_u=float(g.to_unit(nxt[None])[0][IM])))
    return dict(cap_next=cap_next, max_F=float(f.max()), worst=rows)


def book(RE, on, cal, ss, sp):
    # THE TPI's FOOTPRINT ON THE SOLVED GRID: NODES BUYING, NODES AT THE b_DD CORNER, THE LARGEST PURCHASE AND BOOK.
    from global_projection.solver_recursive.decision_rules import SOLVE
    from global_projection.solver_recursive.point_map import point_residuals
    from global_projection.solver_recursive.state_grid import IM
    g, B, rows = on.grid, cal["B_gov_D_ss"], []
    c = dict(cal, tpi_on=True, tpi_cap_bp=on.tpi_cap_solved)
    for i in range(g.n):
        x = np.array([on.vals[k][0][i] for k in SOLVE])
        _, o = point_residuals(g.points[i], 0, x, on, c, ss, sp, RE.N_GH)
        nxt = np.array(g.points[i]); nxt[IM] = o["M_cb_new"]
        rows.append((o["m_cb"] / B, o["b_D_D_new"] / B, o["M_cb_new"] / B, g.to_unit(nxt[None])[0][IM]))
    a = np.array(rows)
    return dict(buying=int((a[:, 0] > 1e-4).sum()), corner=int((a[:, 1] < 1e-3).sum()), max_m=float(a[:, 0].max()),
                max_book=float(a[:, 2].max()), max_next_M_u=float(np.abs(a[:, 3]).max()), nodes=g.n)


def main():
    import run
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.decision_rules import RuleSet
    args = sys.argv[1:]
    fine_grid = "--fine" in args
    floors = tuple(float(v) for v in args[0].split(",")) if args and not args[0].startswith("--") else FLOORS
    tag = "_fine" if fine_grid else ""
    RE.TPI_CAP_START, RE.TPI_CAP_STEP, RE.TPI_MIN_STEP = CAP_START, CAP_STEP, MIN_STEP
    stage0 = RE._stage

    def stage(*a, **kw):
        # A CAP RUNG GETS RUNG_MAXIT NEWTON ITERATIONS; EVERY OTHER STAGE ITS DEFAULT.
        if str(a[6]).startswith("TPI cap="):
            kw["maxit"] = min(kw.get("maxit", 40), RUNG_MAXIT)
        return stage0(*a, **kw)
    RE._stage = stage
    fine = pickle.load(open(run.WARM_START, "rb"))
    bases, hist = {}, []
    log = HERE / "runs" / f"tpi_floor_sweep{tag}.json"

    def baseline(floor, cal, ss, sp):
        # THE COARSE NO-TPI BASELINE AT floor, FROM THE NEAREST SOLVED FLOOR (VIA THE MIDPOINT IF THE JUMP FAILS).
        if not bases:
            seed = fine if fine_grid else RE._seed_from(RuleSet(RE.normal_box(ss, cal, mu=1), 2), fine)
            b = RE.solve_recursive(dict(cal), ss, sp, base=seed, s_refine=0)
            return b if b.solve_ok else None
        near = min(bases, key=lambda f: abs(f - floor))
        b = RE.solve_recursive(dict(cal), ss, sp, base=bases[near], s_refine=0)
        if b.solve_ok:
            return b
        mid = 0.5 * (near + floor)
        print(f"[floor {floor}] direct baseline failed ({b.solve_worst:.1e}); via {mid:.3f}", flush=True)
        cm, sm, pm_ = economy(mid)
        bm = RE.solve_recursive(dict(cm), sm, pm_, base=bases[near], s_refine=0)
        if not bm.solve_ok:
            return None
        b = RE.solve_recursive(dict(cal), ss, sp, base=bm, s_refine=0)
        return b if b.solve_ok else None

    for floor in floors:
        t0, rec = time.perf_counter(), dict(floor=floor)
        try:
            cal, ss, sp = economy(floor)
            base = baseline(floor, cal, ss, sp)
            if base is None:
                rec["baseline"] = "failed"
                print(f"[floor {floor}] the no-TPI baseline did not converge; skipped", flush=True)
                continue
            bases[floor] = base
            rec.update(baseline_worst=float(base.solve_worst), no_tpi=reads(RE, base, cal, ss, sp))
            e = rec["no_tpi"]
            print(f"[floor {floor}] baseline {base.solve_worst:.1e} ({time.perf_counter() - t0:.0f}s); rest mu_D "
                  f"{e['rest']['mu_D']:.4f} sov {e['rest']['sov']:.0f} bp E[Om] {e['rest']['E_Om']:.3f} | "
                  + " | ".join(f"p^d {k}: Q_bD {v['Q_bD']:.3f} (200 bp floor {v['Q_floor_200']:.3f}) sov {v['sov']:.0f} bp "
                               f"mu {v['mu_D']:.3f} spread {v['spread']:+.0f} bp Y {v['dY_D']:+.3f}% C {v['dC_D']:+.3f}%"
                               for k, v in e.items() if k != "rest"), flush=True)
            on = RE.solve_recursive(dict(cal, tpi_on=True, tpi_cap_bp=float(run.TPI_CAP)), ss, sp, base=base, s_refine=0)
            cap = on.tpi_cap_solved
            rec.update(cap_solved=cap, tpi_ok=bool(on.solve_ok), tpi_worst=float(on.solve_worst))
            if cap is not None:
                rec["book"] = book(RE, on, cal, ss, sp)
            print(f"[floor {floor}] TPI: cap reached {cap} bp (target {run.TPI_CAP}), ok {on.solve_ok}, worst "
                  f"{on.solve_worst:.1e}; footprint {rec.get('book')} ({time.perf_counter() - t0:.0f}s)", flush=True)
            if cap is not None and cap > run.TPI_CAP:
                rec["stall"] = stall(RE, on, cal, ss, sp, max(cap - MIN_STEP, float(run.TPI_CAP)))
                for w in rec["stall"]["worst"]:
                    print(f"   stall at {rec['stall']['cap_next']:.1f} bp: {w}", flush=True)
            with open(HERE / "runs" / f"tpi_floor_{100 * floor:.0f}{tag}.pkl", "wb") as fh:
                pickle.dump(dict(base=base, on=on, floor=floor), fh)
        except Exception:
            rec["error"] = traceback.format_exc()
            print(f"[floor {floor}] ERROR\n{rec['error']}", flush=True)
        finally:
            hist.append(rec)
            log.write_text(json.dumps(hist, indent=1, default=lambda z: z.item() if hasattr(z, "item") else str(z)))
    print("[sweep] " + "  ".join(f"floor {h['floor']}: cap {h.get('cap_solved')}" for h in hist), flush=True)


if __name__ == "__main__":
    main()
