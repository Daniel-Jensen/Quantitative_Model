# DRIVER: OUTPUT-RESPONSE DECOMPOSITION + OMT/TPI WELFARE BY INCOME QUINTILE.
# One solve per priced activation probability, then BOTH figures are read off the
# same converged rules and the same simulated paths -- the decomposition and the
# welfare overlay are two views of one experiment, not two experiments.
#   output_decomposition.png  -- d log Y_D split into the credit spread, the
#       deposit rate, the capital stock, the relative price and TFP, with and
#       without the backstop, plus the backstop's effect channel by channel.
#   omt_welfare_quintiles.png -- consumption-equivalent welfare cost of the risk
#       shock and the backstop's improvement, by SS income quintile.
# Serial pointwise solves (no spawn hazard); __main__ guard kept regardless.
import time

import numpy as np

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import s_process_params, default_prob
from solver_recursive.recursive_experiment import s_from_pd
from solver_recursive.recursive_main import calibrate_household_anchors
from solver_recursive.recursive_experiment import solve_recursive
from solver_recursive.output_decomposition import (simulate, s_decay_path,
                                                   decompose_output, active_channels)
from solver_recursive import welfare_quintiles as wq

ACTIVATIONS = (0.0, 0.5, 1.0)   # per-period TPI activation probabilities phi
# The shock is a TARGET one-quarter-ahead default probability (main.py's
# RISK_SHOCK_PD), not a hard-wired s. The old constant -3.9 was labelled
# "+2 sigma": at the calibrated sigma_s = 0.63 it is +4.77 sigma, and a genuine
# +2 sigma shock is p^d = 0.35%, not 2%.
S_SHOCK = s_from_pd(0.0198)     # = main.py's RISK_SHOCK_PD
T_IRF = 21                      # horizon drawn in the decomposition figure
T_WELFARE = 400                 # welfare horizon: beta^T must kill the terminal term


def _label(a):
    # SCENARIO LABEL SHARED BY THE TABLES AND BOTH FIGURES.
    return "no backstop" if a == 0.0 else f"backstop priced at {round(a * 100):.0f}%"


def run(cal, ss, sproc, mu=1, activations=ACTIVATIONS, verbose=True):
    # SOLVE EACH ACTIVATION, SIMULATE, DECOMPOSE, PRICE WELFARE, WRITE THE FIGURES.
    from reporting.plots import (plot_output_decomposition, plot_welfare_quintiles,
                                 OUTDIR)
    env = wq.ss_household_inputs(ss, cal)
    V_ss, a_pol_ss = wq.steady_state_value(ss, cal, env)
    W, inc_q, _ = wq.income_quintile_weights(ss, cal, env)
    # the date-0 cohorts drift through the asset grid even with no shock, so the
    # incidence panel is read against the cohorts' OWN no-shock path
    Cq_ref = wq.ss_cohort_consumption(ss, a_pol_ss, W, T_IRF)

    s_shock = s_decay_path(sproc, S_SHOCK, T_WELFARE)
    s_flat = np.full(T_WELFARE, sproc["s_star"])
    t0 = time.perf_counter()
    cases, cost, cons0 = [], [], None
    for a in activations:
        cal["phi_ltro"] = a
        # coarse grid: this decomposition compares SCENARIOS against each other, so
        # the common grid error differences out; s_refine=5 is available if a level
        # rather than a difference is wanted
        rules = solve_recursive(cal, ss, sproc, mu=mu, verbose=False, s_refine=0,
                                with_cb=(a > 0.0))
        cal["phi_ltro"] = a                             # solve_recursive restores cal
        sim = simulate(rules, cal, ss, sproc, s_shock)
        ref = simulate(rules, cal, ss, sproc, s_flat)
        dec = decompose_output({k: v[:T_IRF] for k, v in sim.items() if np.ndim(v)},
                               {k: v[:T_IRF] for k, v in ref.items() if np.ndim(v)},
                               cal)
        cases.append((_label(a), dec))

        inp = wq.aggregate_inputs(sim, ref, ss, cal, env)
        V0, c_path, a_pol_path = wq.transition_welfare(inp, ss, cal, V_ss)
        cost.append(wq.quintile_cev(V0, V_ss, W, ss, cal))
        if cons0 is None:
            Cq = wq.cohort_consumption(ss, c_path, a_pol_path, W, T_IRF)
            cons0 = 100.0 * (Cq / Cq_ref - 1.0)
        if verbose:
            print(f"  [{_label(a):24s}] solved ({time.perf_counter() - t0:.0f}s)  "
                  f"Y_D[0]={dec['total'][0]:+.3f}%  worst resid="
                  f"{np.max(sim['resid']):.1e}  worst labour FOC="
                  f"{np.max(sim['resid_lab']):.1e}  off-box={sim['off_box']:+.3f}  "
                  f"slack/fail={sim['n_slack']}/{sim['n_fail']}  "
                  f"terminal drift={inp['drift']:.1e}", flush=True)

    cost = np.array(cost)
    gain = cost - cost[0]
    chans = active_channels(cases[0][1])
    note = (f"recursive Chebyshev-Smolyak projection, mu={mu} "
            f"({rules.grid.n} grid points x 2 regimes), three-branch quadrature; "
            f"risk shock (priced default {100 * default_prob(S_SHOCK):.1f}% "
            f"per quarter) decaying at rho_s={sproc['rho_s']}")
    plot_output_decomposition(cases, chans, note=note)
    plot_welfare_quintiles([lab for lab, _ in cases], cost, gain, cons0, inc_q,
                           note=note)

    print("\n  OUTPUT DECOMPOSITION — contribution to Y_D (% dev. from the "
          "no-shock path)")
    for lab, dec in cases:
        print(f"   {lab}")
        print("    qtr  " + "".join(f"{n[:13]:>14s}" for _, n in chans)
              + f"{'TOTAL':>14s}")
        for t in (0, 1, 2, 4, 8, 12, 20):
            if t < T_IRF:
                print(f"    {t:3d}  "
                      + "".join(f"{dec[k][t]:+14.4f}" for k, _ in chans)
                      + f"{dec['total'][t]:+14.4f}")

    print("\n  OMT/TPI WELFARE BY SS INCOME QUINTILE (consumption-equivalent, %)")
    print("   scenario                     " + "".join(f"{f'Q{q + 1}':>10s}"
                                                       for q in range(cost.shape[1])))
    for i, (lab, _) in enumerate(cases):
        print(f"   cost of shock {lab:15s}"
              + "".join(f"{cost[i, q]:+10.4f}" for q in range(cost.shape[1])))
    for i, (lab, _) in enumerate(cases[1:], start=1):
        print(f"   GAIN vs no backstop {lab[:9]:9s}"
              + "".join(f"{gain[i, q]:+10.4f}" for q in range(cost.shape[1])))
    print("   mean SS income               "
          + "".join(f"{inc_q[q]:10.4f}" for q in range(cost.shape[1])))
    print(f"\n  figures -> {OUTDIR}/output_decomposition.png, "
          f"{OUTDIR}/omt_welfare_quintiles.png", flush=True)


def main():
    import sys
    mu = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    cal = get_calibration()
    cal["nw_floor_frac"] = 0.15          # match main.py
    # NOTE: this used to silently override f = 0.12 and the credit spread to 0.018
    # (720bp/yr) whenever mu >= 2 -- the "720bp detour" main.py's own header calls a
    # regression, and 90x the calibrated 8bp. An experiment must not re-calibrate the
    # model behind the caller's back; if that variant is wanted it belongs in
    # calibration.py where it can be seen.
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print(f"=== output decomposition + OMT welfare by income quintile (mu={mu}, "
          f"pd peak {100 * default_prob(S_SHOCK):.2f}%) ===", flush=True)
    run(cal, ss, sproc, mu=mu)


if __name__ == "__main__":
    main()
