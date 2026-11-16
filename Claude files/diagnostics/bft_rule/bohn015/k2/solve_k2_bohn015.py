# THE BFT PROPORTIONAL RULE AT k = 2, ANNOUNCED, AT BOHN 0.15 (2026-11-13). SAME PROCEDURE AS run_bft.py AT 0.25:
# START FROM THE SAVED, CONVERGED NO-POLICY RULES AT BOHN 0.15 (tpi_root_cause/runs/rules_12pct_bohn0.150.pkl), WALK
# k 0 -> 1 -> 2, HARD POLISH. THE RULE'S REFERENCE EXCESS RETURN IS ITS OWN NO-POLICY REST POINT'S, AS AT 0.25.
# NOTHING IN THE PRIVATE ECONOMY, THE RESERVE TREATMENT, THE RULE OR THE SOLVER IS CHANGED; THE CALIBRATION FILE IS
# NOT EDITED (BOHN IS SET IN THIS RUN'S CALIBRATION DICTIONARY ONLY).
#   python3 -u "Claude files/diagnostics/bft_rule/bohn015/k2/solve_k2_bohn015.py"   (from the worktree root)
import copy, json, pickle, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE.parents[1])); sys.path.insert(0, str(HERE.parent))
from global_projection.solver_recursive.state_grid import IM
from global_projection.solver_recursive import recursive_experiment as RE
from run_bft import rest_row, node_checks
from bohn_compare import Econ, RUNS

OUT = ROOT.parent / "QUANTITATIVE_MODEL" / "results" / "GLOBAL_2026-11-12_pr38_bft" / "bohn015" / "k2"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    t0 = time.perf_counter()

    E = Econ(0.15, RUNS / "rules_12pct_bohn0.150.pkl")
    r0 = E.rules
    B = E.cal["B_gov_D_ss"]
    # the continuation rules were pickled without the box description; the grid is identical to PR #38's (checked
    # below), so its box travels with them (only a refinement in s, not used here, reads it)
    r38 = pickle.load(open(ROOT.parent / "QUANTITATIVE_MODEL" / "results" / "GLOBAL_2026-11-12_pr38" / "data" /
                           "rules_risk.pkl", "rb"))
    assert np.allclose(r0.grid.points, r38.grid.points) and np.allclose(r0.grid.hi, r38.grid.hi)
    r0.box = r38.box
    print(f"Bohn {E.cal['bohn_annual']}: no-policy rules {r0.grid.n} points, solve_ok {getattr(r0, 'solve_ok', None)}, worst "
          f"{getattr(r0, 'solve_worst', float('nan')):.1e}, smoothing {getattr(r0, 'smoothing', None)}; book band {100 * r0.grid.hi[IM] / B:.1f}% "
          f"of the debt", flush=True)
    r25 = pickle.load(open(RUNS / "rules_12pct_bohn0.250.pkl", "rb"))
    print(f"  grid identical to the 0.25 solution: {np.allclose(r0.grid.points, r25.grid.points)}; "
          f"box hi identical: {np.allclose(r0.grid.hi, r25.grid.hi)}", flush=True)
    chk = node_checks(r0, E.cal, E.ss, E.sp)
    print(f"  saved solution reproduces at the nodes: max|F| {chk['max_resid']:.1e}", flush=True)
    rest_off, S_off = rest_row(r0, E.cal, E.ss, E.sp)
    er_ref = rest_off["er_bp"] / 4e4
    print(f"  no-policy rest point: ER {rest_off['er_bp']:.2f} bp/yr (the rule's reference; 15.38 at Bohn 0.25), "
          f"sovereign spread {rest_off['sov_bp']:.1f} bp, mu_D {rest_off['mu_D']:.2e}", flush=True)

    cal_k = dict(E.cal, tpi_on=True, tpi_rule="bft", tpi_gz=False, tpi_bft_er_ref=er_ref, tpi_bft_k=2.0,
                 tpi_bft_from=0.0)
    ss_k = copy.deepcopy(E.ss)
    r = RE.solve_recursive(cal_k, ss_k, E.sp, s_refine=0, base=r0)
    ok = bool(r.solve_ok) and r.tpi_bft_solved == 2.0
    with open(OUT / "rules_bft_k2_bohn015.pkl", "wb") as fh:
        pickle.dump(r, fh)
    mins = (time.perf_counter() - t0) / 60
    print(f"\nk = 2 at Bohn 0.15: ok {ok}, last k solved {r.tpi_bft_solved}, max|F| {r.solve_worst:.1e}, {mins:.0f} min",
          flush=True)
    info = dict(ok=ok, solved=r.tpi_bft_solved, worst=float(r.solve_worst), minutes=mins, er_ref_bp=rest_off["er_bp"],
                smoothing=getattr(r, "smoothing", None), rest_off=rest_off)
    if ok:
        nodes = node_checks(r, cal_k, ss_k, E.sp)
        info["nodes"] = {kk: v for kk, v in nodes.items() if kk != "rows"}
        print(f"  nodes: max|F| {nodes['max_resid']:.1e}; buying at {nodes['n_buy']}/{nodes['n']}; largest next book "
              f"{100 * nodes['max_M_new']:.1f}% of the debt against a band of {100 * nodes['band_M']:.1f}%; next book "
              f"outside the box at {nodes['n_out_M']} nodes; other states outside at {nodes['n_out_other']}; min D-bank "
              f"book {nodes['min_b_DD']:.3f}", flush=True)
    json.dump(info, open(OUT / "solve_info.json", "w"), indent=1,
              default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    print("done", flush=True)


if __name__ == "__main__":
    main()
