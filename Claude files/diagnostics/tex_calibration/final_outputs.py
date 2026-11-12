# THE RISK-SHOCK OUTPUTS OF run.py's GLOBAL PIPELINE, ON A SOLVED APPENDIX-B RULE SET (2026-10-29): THE REST POINT, THE
# IMPACT TABLE, THE PERSISTENCE AND DYNAMIC IRFs, THE CALIBRATION TABLE, THE DECOMPOSITIONS, THE EULER-ERROR ACCURACY,
# AND THE FIGURES. THE COLD-START SOLVE AND THE TFP AND TPI SOLVES OF run.py ARE LEFT OUT; results/GLOBAL IS UNTOUCHED.
#   python3 -u final_outputs.py <checkpoint> <out dir name>
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)


def main():
    import run
    import results_io
    from results_io import save
    from global_projection.solver_recursive.state_grid import s_process_params, STATE_NAMES
    from global_projection.solver_recursive.recursive_experiment import (report_rest_point, impact_table,
                                                                         persistence_irf, dynamic_irf, s_from_pd)
    from global_projection.solver_recursive.output_decomposition import (simulate, s_decay_path, decompose_output,
                                                                         active_channels, decompose_bond_price,
                                                                         exact_bond_reads, BOND_CHANNELS)
    from global_projection.solver_recursive.accuracy import accuracy_report
    from global_projection.solver_recursive.calibrate_stochastic import (rest_point_read, rest_point_moments,
                                                                         calibration_values)
    from global_projection.reporting.prints import (banner, print_output_decomposition, print_bond_decomposition,
                                                    print_calibration_table)
    from global_projection.reporting import plots
    from global_projection.reporting.export import export_calibration_tex
    sv = pickle.load(open(HERE / "runs" / sys.argv[1], "rb"))
    rules, cal, ss = sv["rules"], sv["cal"], sv["ss"]
    sproc = s_process_params(cal)
    data, figures, _ = results_io.model_dirs(HERE / sys.argv[2], fresh=True)
    pd_shock, T = run.RISK_SHOCK_PD, run.DECOMP_T
    banner(f"Appendix-B risk shock on {sys.argv[1]}: grid {rules.grid.n} points, f_D {cal['f_D']}")
    S_rest = report_rest_point(rules, cal, ss, sproc)
    impact_table(rules, cal, ss, sproc)
    persistence = persistence_irf(rules, cal, ss, sproc, pd_shock=pd_shock)
    risk = dynamic_irf(rules, cal, ss, sproc, pd_shock=pd_shock, T=25)
    banner("Calibration table")
    _, o_rest = rest_point_read(rules, cal, ss, sproc)
    moments = rest_point_moments(o_rest)
    calib = calibration_values(o_rest, cal, sproc, risk)
    print_calibration_table(calib, moments, cal)
    banner("Decompositions")
    s_path = s_decay_path(sproc, s_from_pd(pd_shock), T)
    sim = simulate(rules, cal, ss, sproc, s_path, S_init=S_rest)
    ref = simulate(rules, cal, ss, sproc, [sproc["s_star"]] * T, S_init=S_rest)
    output_dec = decompose_output(sim, ref, cal)
    print_output_decomposition(output_dec, active_channels(output_dec))
    bond_dec = decompose_bond_price(exact_bond_reads(rules, cal, ss, sproc, sim),
                                    exact_bond_reads(rules, cal, ss, sproc, ref), cal)
    print_bond_decomposition(bond_dec, BOND_CHANNELS)
    banner("Accuracy")
    accuracy = accuracy_report(rules, cal, ss, sproc, S_rest, T=run.ACCURACY_T, label="appendix-B risk rules")
    for name, obj in (("calibration", cal), ("steady_state", ss), ("rest_point", dict(zip(STATE_NAMES, S_rest))),
                      ("persistence_irf", persistence), ("risk_irf", risk),
                      ("calibration_table", dict(values=calib, moments=moments)),
                      ("output_decomposition", output_dec), ("bond_decomposition", bond_dec), ("accuracy", accuracy)):
        save(data / name, obj)
    plots.OUTDIR = str(figures)
    shock = f"p-d shock to {100 * pd_shock:.2f}%/qtr"
    paths = [plots.plot_risk_irf(risk, note=f"{shock}, states evolving"),
             plots.plot_output_decomposition([(shock, output_dec)], active_channels(output_dec),
                                             note="d log Y_D by the production function and the GHH labour FOC"),
             plots.plot_bond_decomposition(bond_dec, list(BOND_CHANNELS),
                                           note="the D bank first-order condition, split leg by leg"),
             export_calibration_tex(calib, str(figures / "calibration_table.tex"))]
    for p in paths:
        print(f"  figure -> {p}")


if __name__ == "__main__":
    main()
