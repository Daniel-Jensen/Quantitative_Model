# THE ANNOUNCED BFT RULE AT k = 2 ON THE CANONICAL CALIBRATION (BOHN 0.15; HEADLINE SHOCK 2%), 2026-11-13.
# STARTS FROM THE NO-POLICY RULES THE GLOBAL PIPELINE JUST SOLVED (results/GLOBAL/data of this worktree), WALKS k
# 0 -> 1 -> 2 AND HARD-POLISHES, AS run_bft.py DID AT 0.25. THE RULE'S REFERENCE EXCESS RETURN IS THE NO-POLICY REST
# POINT'S. NOTHING IN THE PRIVATE ECONOMY, THE RESERVE TREATMENT, THE RULE OR THE SOLVER IS CHANGED.
#   python3 -u "Claude files/diagnostics/bft_rule/bohn015/canonical/solve_k2.py"   (from the worktree root)
import copy, json, pickle, shutil, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE.parents[1]))
from results_io import load
from global_projection.solver_recursive.state_grid import s_process_params, IM
from global_projection.solver_recursive import recursive_experiment as RE
from run_bft import rest_row, node_checks

DATA = ROOT / "results" / "GLOBAL" / "data"
OUT = ROOT.parent / "QUANTITATIVE_MODEL" / "results" / "GLOBAL_2026-11-13_bohn015"


def main():
    t0 = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    # the canonical no-policy run, copied whole beside the announced rule (the pipeline's folder is overwritten by
    # the next run)
    if not (OUT / "no_policy").exists():
        shutil.copytree(DATA.parent, OUT / "no_policy")
    cal = load(DATA / "calibration")
    ss = load(DATA / "steady_state")
    cal["tpi_on"] = False
    sp = s_process_params(cal)
    r0 = pickle.load(open(DATA / "rules_risk.pkl", "rb"))
    B = cal["B_gov_D_ss"]
    print(f"Bohn {cal['bohn_annual']}: no-policy rules {r0.grid.n} points, solve_ok {r0.solve_ok}, worst "
          f"{r0.solve_worst:.1e}; book band {100 * r0.grid.hi[IM] / B:.1f}% of the debt", flush=True)
    chk = node_checks(r0, cal, ss, sp)
    print(f"  no-policy rules at the nodes: max|F| {chk['max_resid']:.1e}", flush=True)
    rest_off, _ = rest_row(r0, cal, ss, sp)
    er_ref = rest_off["er_bp"] / 4e4
    print(f"  no-policy rest point: ER {rest_off['er_bp']:.3f} bp/yr (the rule's reference), sovereign spread "
          f"{rest_off['sov_bp']:.2f} bp, mu_D {rest_off['mu_D']:.2e}", flush=True)
    cal_k = dict(cal, tpi_on=True, tpi_rule="bft", tpi_gz=False, tpi_bft_er_ref=er_ref, tpi_bft_k=2.0,
                 tpi_bft_from=0.0)
    ss_k = copy.deepcopy(ss)
    r = RE.solve_recursive(cal_k, ss_k, sp, s_refine=0, base=r0)
    ok = bool(r.solve_ok) and r.tpi_bft_solved == 2.0
    with open(OUT / "rules_bft_k2.pkl", "wb") as fh:
        pickle.dump(r, fh)
    mins = (time.perf_counter() - t0) / 60
    print(f"\nk = 2: ok {ok}, last k solved {r.tpi_bft_solved}, max|F| {r.solve_worst:.1e}, {mins:.0f} min", flush=True)
    info = dict(ok=ok, solved=r.tpi_bft_solved, worst=float(r.solve_worst), minutes=mins, er_ref_bp=rest_off["er_bp"],
                smoothing=getattr(r, "smoothing", None), rest_off=rest_off, nodes_off={k: v for k, v in chk.items()
                                                                                       if k != "rows"})
    if ok:
        nodes = node_checks(r, cal_k, ss_k, sp)
        info["nodes"] = {kk: v for kk, v in nodes.items() if kk != "rows"}
        info["nodes"]["min_M_new_pct"] = float(100 * nodes["rows"]["M_new"].min())
        info["nodes"]["min_m_pct"] = float(100 * nodes["rows"]["m"].min())
        print(f"  nodes: max|F| {nodes['max_resid']:.1e}; buying at {nodes['n_buy']}/{nodes['n']}; book from "
              f"{info['nodes']['min_M_new_pct']:.3f}% to {100 * nodes['max_M_new']:.1f}% of the debt (band "
              f"{100 * nodes['band_M']:.1f}%); next book outside the box at {nodes['n_out_M']} nodes; other states "
              f"outside at {nodes['n_out_other']}", flush=True)
    json.dump(info, open(OUT / "solve_k2_info.json", "w"), indent=1,
              default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    print("done", flush=True)


if __name__ == "__main__":
    main()
