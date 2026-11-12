# M3 (2026-10-28): THE FINE CHECK OF THE CLEAN COARSE APPENDIX-B SOLUTION: THE s DIMENSION REFINED TO BOCOLA'S 9
# NODES (CROSS MODE), RE-SOLVED THROUGH solve_recursive(base=...) (CALM GUARD, HARD POLISH), THEN THE SAME READS AS M1,
# SO THE COARSE FITTED-vs-EXACT GAP AT p^d 2% (C_D, Y_D) IS SETTLED ON THE FINE GRID.
#   python3 -u m3_fine.py <hard coarse checkpoint> [s nodes]
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(HERE))


def main():
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import s_process_params
    from m1_coarse import reads
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    nref = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    cal["jacobian"] = "complex"
    sp = s_process_params(cal)
    fine = RE.solve_recursive(cal, ss, sp, base=rules, s_refine=nref, verbose=True)
    print(f"[m3] solve ok {fine.solve_ok}  worst {fine.solve_worst:.2e}  grid {fine.grid.n} points", flush=True)
    with open(HERE / "runs" / f"m3_fine_s{nref}.pkl", "wb") as fh:
        pickle.dump(dict(sv, rules=fine, cal=cal), fh)
    reads(RE, fine, cal, ss, sp, f"fine s{nref}")


if __name__ == "__main__":
    main()
