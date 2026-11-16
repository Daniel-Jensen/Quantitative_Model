# GATE-5 MEASUREMENT (2026-10-28): DOES A HIGHER D BANKER EXIT RATE f_D BIND THE D BANK AT THE REST POINT? THE
# FRANCHISE VALUE'S FORWARD WEIGHT SCALES ~1/f, AND THE REST POINT IS 12% SLACK AT f 0.08 FOR ANY SS SPREAD (THE
# BIND PHASE OF m2_warm.py). WARM SOLVES OF THE HARD SYSTEM AT f_D IN THE LIST; THE SLACK RATIO AND THE 2% SHOCK READ.
#   python3 -u probe_f.py <hard checkpoint> f1 f2 ...
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))
import numpy as np


def read(RE, CUR):
    # THE SLACK RATIO, mu AND THE SPREAD AT THE REST POINT AND AT p^d 2% (EXACT READS), AND THE IMPACT ROW.
    from measure import omega_legs, row
    from global_projection.solver_recursive.state_grid import IS
    r, cal, ss, sp = CUR["rules"], CUR["cal"], CUR["ss"], CUR["sp"]
    S0 = RE.stochastic_rest_point(r, cal, ss, sp, verbose=False)
    out = []
    for pd in (None, 0.02):
        S = S0.copy()
        if pd:
            S[IS] = RE.s_from_pd(pd)
        b = omega_legs(r, cal, ss, sp, S, RE.read_exact)
        out.append((b, b["o"]))
    (b0, o0), (b1, o1) = out
    sr = lambda b, o: b["E_Om"] * (1 + o["rdep_D"]) * b["n_D"] / b["lev_D"]
    rw = row(o1, o0)
    return (f"rest: slack {sr(b0, o0):.4f} mu_D {b0['mu_D']:.4f} spread {4e4 * (o0['r_wc_D'] - o0['rdep_D']):.1f}bp "
            f"leverage {(o0['dep_D'] + o0['n_D']) / o0['n_D']:.2f} E[Om] {b0['E_Om']:.3f} | p^d 2%: slack {sr(b1, o1):.4f} "
            f"mu_D {b1['mu_D']:.4f}; dQ_bD {rw['Q_bD']:+.2f}% dn_D {rw['n_D']:+.2f}% dY_D {rw['Y_D']:+.3f}% "
            f"dC_D {rw['C_D']:+.3f}% dspread {rw['spread_bp']:+.1f}bp")


def main():
    from calibration.global_projection import get_calibration
    import global_projection.solver_recursive.calibrate_stochastic as CS
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import s_process_params
    from diag_path import FLOOR
    import m2_warm as M
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    base = get_calibration()
    base.update(nw_floor_frac=FLOOR, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3, jacobian="complex")
    M.CUR.update(rules=sv["rules"], cal=sv["cal"], ss=sv["ss"], sp=s_process_params(sv["cal"]),
                 vals=dict(sv["vals"], **M.HARD), base=base)
    x = np.array([M.CUR["cal"][p[0]] for p in CS.PAIRS])
    print(f"[f_D {M.CUR['cal']['f_D']:.3f}] {read(RE, M.CUR)}", flush=True)
    f_now = float(M.CUR["cal"]["f_D"])
    for f in map(float, sys.argv[2:]):
        while True:
            step = f
            M.CUR["vals"]["f_D"] = step
            ok, w = M.solve_at(RE, x)
            if ok:
                f_now = step
                break
            # halve the move in f on a failed solve
            M.CUR["vals"]["f_D"] = f_now
            f = 0.5 * (f_now + f)
            print(f"  [f] solve failed; trying f_D {f:.4f}", flush=True)
            if abs(f - f_now) < 1e-3:
                return
        print(f"[f_D {f_now:.3f}] solved {w:.1e}; {read(RE, M.CUR)}", flush=True)
        with open(HERE / "runs" / f"probe_f_{f_now:.3f}.pkl", "wb") as fh:
            pickle.dump({k: M.CUR[k] for k in ("rules", "cal", "ss", "vals")}, fh)


if __name__ == "__main__":
    main()
