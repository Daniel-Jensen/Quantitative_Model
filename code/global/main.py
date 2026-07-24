# ENTRY POINT: STEADY STATE -> TFP IRF -> BOCOLA RISK PASS-THROUGH -> TPI BACKSTOP.
# Every table and headline number is printed from prints.py.
import os
import time

import numpy as np

from calibration import get_calibration
from steady_state import solve_steady_state
from transition import solve_transition
from risk_branch import solve_transition_risk
from plots import OUTDIR, plot_irf, plot_default_irf, plot_bond_decomposition
from prints import (banner, print_ss_table, print_transition_residuals,
                    print_risk_table, print_tpi_table)

RUN_TFP  = True   # experiment 1: TFP shock in D
RUN_RISK = True   # experiment 2: exogenous sovereign-risk shock (Bocola)
RUN_TPI  = False   # experiment 3: TPI backstop under the same risk shock
# headline = LIQUIDITY CHANNEL ONLY (risk-neutral pricing of pi_t): correct
# signs, tested in tests/test_signs_bocola.py. The two-branch risk channel
# currently produces an expansionary artifact (capital is the branch safe
# haven; see docs) and is being rebuilt against the Bocola (2016) global
# replication in code/bocola2016/ -- keep it opt-in until that lands.
RUN_RISK_TWO_BRANCH = False

TFP_SHOCK, TFP_RHO = 0.01, 0.9    # 1% impact, AR(1) decay

# Exogenous priced-default-probability path pi_t (Bocola's s-shock, eqs. 11-12):
# impact level x AR(1) decay. Risk is priced, never realized.
PI_SHOCK, PI_RHO = 0.02, 0.95

# TPI: TPI_S_ON_DATE is when the mechanical price-floor purchases activate
# (0 = the risk shock's own impact date). TPI_DOUBT/_RHO are the peak
# probability the backstop is doubted and its decay, mirroring PI_SHOCK/PI_RHO;
# the priced-survival path handed to the solver is 1 - doubt.
TPI_S_ON_DATE = 0
TPI_DOUBT, TPI_DOUBT_RHO = 0.03, 0.9


def _risk_paths(cal):
    # SHOCK-FREE TFP PATHS PLUS THE EXOGENOUS PRICED-DEFAULT PATH pi_t.
    T = cal["T"]
    return (np.full(T, cal["Z_ss_D"]), np.full(T, cal["Z_ss_F"]),
            PI_SHOCK * PI_RHO ** np.arange(T))


def run_tfp(ss, cal):
    # EXPERIMENT 1: TFP SHOCK IN D.
    banner(f"TFP shock in D: {TFP_SHOCK:.0%} on impact, rho = {TFP_RHO}")
    T = cal["T"]
    Z_D = cal["Z_ss_D"] * np.exp(TFP_SHOCK * TFP_RHO ** np.arange(T))
    Z_F = np.full(T, cal["Z_ss_F"])

    t0 = time.perf_counter()
    out = solve_transition(ss, cal, Z_D, Z_F, verbose=False)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_transition_residuals(out, cal)

    plot_irf(out, ss, cal)
    return out


def run_risk(ss, cal):
    # EXPERIMENT 2: EXOGENOUS SOVEREIGN-RISK SHOCK, LIQUIDITY CHANNEL HEADLINE.
    # Risk-neutral pricing of the exogenous pi_t path (MTM losses tighten the
    # IC today) -- the subset of Bocola's decomposition with verified signs.
    banner(f"Sovereign-risk shock in D (liquidity channel): peak priced "
           f"default prob {PI_SHOCK:.0%} per quarter, rho = {PI_RHO}")
    Z_D, Z_F, pi_D = _risk_paths(cal)

    t0 = time.perf_counter()
    out = solve_transition(ss, cal, Z_D, Z_F, def_price_D=pi_D, verbose=False)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_transition_residuals(out, cal)
    print_risk_table(out, ss, cal, PI_SHOCK, PI_RHO)

    plot_default_irf(out, ss, cal)
    plot_bond_decomposition(out, ss, cal)

    if RUN_RISK_TWO_BRANCH:
        banner("Two-branch risk channel [UNDER REPAIR: expansionary artifact]")
        t0 = time.perf_counter()
        out_rb = solve_transition_risk(ss, cal, Z_D, Z_F, pi_D_path=pi_D,
                                       verbose=True)
        print(f"  solved in {time.perf_counter() - t0:.0f}s")
        print_transition_residuals(out_rb, cal)
        print_risk_table(out_rb, ss, cal, PI_SHOCK, PI_RHO)
        plot_default_irf(out_rb, ss, cal, filename="default_irf_twobranch.png")
        plot_bond_decomposition(out_rb, ss, cal,
                                filename="bond_decomposition_twobranch.png")
    return out


def run_tpi(ss, cal):
    # EXPERIMENT 3: THE SAME RISK SHOCK, NOW WITH THE TPI BACKSTOP ACTIVE.
    # Financed by an explicit CB-budget-consistent rebate; tests/test_tpi.py's
    # budget-closure test is the regression guard against an accounting leak.
    banner(f"TPI backstop under the sovereign-risk shock in D: peak priced "
           f"default prob {PI_SHOCK:.0%} per quarter, rho = {PI_RHO}; backstop "
           f"on from t={TPI_S_ON_DATE}, peak doubt {TPI_DOUBT:.0%}")
    T = cal["T"]
    Z_D, Z_F, pi_D = _risk_paths(cal)
    pi_tpi_D = 1.0 - TPI_DOUBT * TPI_DOUBT_RHO ** np.arange(T)
    s_tpi_D = np.zeros(T)
    s_tpi_D[TPI_S_ON_DATE:] = 1.0

    t0 = time.perf_counter()
    out_base = solve_transition_risk(ss, cal, Z_D, Z_F, pi_D_path=pi_D, verbose=False)
    out = solve_transition_risk(ss, cal, Z_D, Z_F, pi_D_path=pi_D,
                                pi_tpi_D_path=pi_tpi_D, s_tpi_D_path=s_tpi_D,
                                verbose=True, y0=out_base["y_vec"])
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_transition_residuals(out, cal, ss=ss)
    print_tpi_table(out, out_base, ss, cal, PI_SHOCK, TPI_S_ON_DATE, TPI_DOUBT)

    plot_default_irf(out, ss, cal, filename="tpi_irf.png")
    plot_bond_decomposition(out, ss, cal, filename="tpi_bond_decomposition.png")
    return out


def main():
    # RUN THE FULL PIPELINE AND WRITE EVERY FIGURE TO output/.
    t0 = time.perf_counter()
    os.makedirs(OUTDIR, exist_ok=True)
    cal = get_calibration()

    banner("Two-country HANK-GK monetary union: steady state")
    ss = solve_steady_state(cal)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_ss_table(ss, cal)

    if RUN_TFP:
        run_tfp(ss, cal)
    if RUN_RISK:
        run_risk(ss, cal)
    if RUN_TPI:
        run_tpi(ss, cal)

    print(f"\nFigures saved to {OUTDIR}")
    print(f"TOTAL  {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
