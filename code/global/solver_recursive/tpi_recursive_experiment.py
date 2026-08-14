# OMT/TPI ACTIVATION COMPARISON IN THE RECURSIVE PROJECTION SOLVER (THREE-BRANCH).
# The SAME sovereign-risk (priced default) shock under 0% / 50% / 100% backstop
# activation, solved with Chebyshev decision rules and genuine three-branch
# quadrature (no-default / default-honoured / default-reneged) -- NO representative
# branch anywhere. Activation is the priced probability cal["tpi_activation"] that
# the CB honours the backstop and redeems the D-bond near par; a=0 nests the pure
# two-branch risk solve. Each activation level is a separate rule solve (a is priced,
# not a state), then the impact + persistence IRF is read off the converged rules.
# Serial pointwise solves (no spawn hazard); __main__ guard kept regardless.
import numpy as np

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import s_process_params, default_prob
from solver_recursive.recursive_main import ss_state, calibrate_household_anchors
from solver_recursive.recursive_experiment import solve_recursive, read_at, _spread_bp

DEFAULT_ACTIVATIONS = (0.0, 0.5, 1.0)      # priced OMT/TPI activation probabilities
_TPI_COLORS = ("#d62728", "#ff7f0e", "#2ca02c", "#1f77b4", "#9467bd", "#8c564b")
S_SHOCK, T_IRF = -3.9, 21     # +2 sigma risk shock (pd ~ 2%), decaying over T_IRF


def _activation_specs(activations):
    # MAP ACTIVATION PROBABILITIES TO (value, label, color) TRIPLES FOR SOLVE + PLOT.
    specs = []
    for i, a in enumerate(activations):
        lab = "0% — known not active" if a == 0.0 else f"{round(a * 100):.0f}% chance"
        specs.append((a, lab, _TPI_COLORS[i % len(_TPI_COLORS)]))
    return specs


def irf_series(rules, cal, ss, sproc, a):
    # PERSISTENCE IRF READ OFF THE CONVERGED RULES (s-shock decays; states at SS).
    cal["tpi_activation"] = a
    S0 = ss_state(ss, cal, sproc)
    base = read_at(rules, cal, ss, sproc, S0.copy())
    Yb, Cb, Ib, Qb = base["Y_D"], base["C_D"], base["I_D"], base["Q_bD"]
    sp_b = _spread_bp(base, cal)
    keys = ("pd", "Y_D", "C_D", "I_D", "Q_bD", "spread")
    P = {k: np.empty(T_IRF) for k in keys}
    for t in range(T_IRF):
        s_t = sproc["s_star"] + sproc["rho_s"] ** t * (S_SHOCK - sproc["s_star"])
        S = S0.copy(); S[5] = s_t
        o = read_at(rules, cal, ss, sproc, S)
        P["pd"][t] = 100 * default_prob(s_t)
        P["Y_D"][t] = 100 * (o["Y_D"] / Yb - 1)
        P["C_D"][t] = 100 * (o["C_D"] / Cb - 1)
        P["I_D"][t] = 100 * (o["I_D"] / Ib - 1)
        P["Q_bD"][t] = 100 * (o["Q_bD"] / Qb - 1)
        P["spread"][t] = _spread_bp(o, cal) - sp_b
    return P


def run(cal, ss, sproc, mu=1, activations=None):
    # SOLVE THE ACTIVATION SCENARIOS AND WRITE THE OVERLAY FIGURE (shared SS).
    from reporting.plots import plot_activation_irf, OUTDIR
    import os, time
    os.makedirs(OUTDIR, exist_ok=True)
    specs = _activation_specs(DEFAULT_ACTIVATIONS if activations is None else activations)
    print("=== OMT/TPI activation comparison — recursive projection (three-branch) ===",
          flush=True)
    t0 = time.perf_counter()
    scenarios = []
    for a, label, color in specs:
        cal["tpi_activation"] = a
        rules = solve_recursive(cal, ss, sproc, mu=mu, verbose=False)
        P = irf_series(rules, cal, ss, sproc, a)
        scenarios.append((label, P, color))
        print(f"  [{label}] solved ({time.perf_counter()-t0:.0f}s)  "
              f"impact Y_D={P['Y_D'][0]:+.3f}%  Q_bD={P['Q_bD'][0]:+.2f}%  "
              f"spread={P['spread'][0]:+.0f}bp", flush=True)

    plot_activation_irf(scenarios)

    print("\n  IMPACT (t=0) of the SAME default shock by activation probability:",
          flush=True)
    print("   scenario                 Y_D%     C_D%    I_D%    Q_bD%   spread_bp", flush=True)
    for label, P, _ in scenarios:
        print(f"   {label:24s} {P['Y_D'][0]:+6.3f}  {P['C_D'][0]:+6.3f}  "
              f"{P['I_D'][0]:+6.3f}  {P['Q_bD'][0]:+6.2f}  {P['spread'][0]:+7.0f}",
              flush=True)
    print(f"\n  figure -> {OUTDIR}/tpi_activation_recursive.png", flush=True)


def main():
    cal = get_calibration()
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    run(cal, ss, sproc)


if __name__ == "__main__":
    main()
