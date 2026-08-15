# ENTRY POINT: CHEBYSHEV-SMOLYAK PROJECTION SOLVER FOR THE WHOLE MODEL.
# The two-country model is solved GLOBALLY as recursive decision rules on a Smolyak
# grid (solver_recursive/) -- there is NO perfect-foresight / representative-branch
# machinery. main() solves the steady state once, then runs the projection
# experiments, each a full time-iteration solve: a TFP shock (Z_D as the 7th state) on
# the isotropic MU grid, the Bocola sovereign-risk pass-through on the ANISOTROPIC
# RISK_MU_VEC grid (resolve s + net worth), and optionally (RUN_TPI) the OMT/TPI
# activation overlay over the selectable TPI_ACTIVATIONS. Heavy by design (~1h; more
# with RUN_TPI).
import time

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import s_process_params
from solver_recursive.recursive_main import calibrate_household_anchors
from solver_recursive.recursive_experiment import (
    solve_tfp, tfp_irf, solve_recursive, impact_table, persistence_irf)
from solver_recursive import tpi_recursive_experiment
from reporting.prints import banner, print_ss_table

# PRODUCTION SOLVER at the Bocola-faithful calibration (leverage 4, 200bp, mu_ss~0.02;
# see config/calibration.py). The SOVEREIGN-RISK + OMT/TPI channels use an ANISOTROPIC
# Smolyak grid (RISK_MU_VEC): resolve the risk state s and both banks' net worth P, keep
# the near-fixed stocks at their centre -- Bocola's low-dimensional approach, concentrate
# resolution where the mechanism lives. This (i) reprices the long bond with its true
# persistence (Q_bD ~ -4.8% vs the flat -1.3% of the isotropic mu=1 grid) and (ii) lets
# the bond loss crash net worth so the occasionally-binding mu spikes ENDOGENOUSLY (mu up,
# spread up -- correctly signed), all at the realistic 200bp with NO raised spread. The
# isotropic 720bp/mu=2 detour was a regression (removed; see git history + the Bocola-2016
# code comparison). NW_FLOOR is Bocola's net-worth floor (his N_tom=max(.,0.65)): it keeps
# the deep default-regime corners feasible so the d=1 fit does not poison the global basis.
# TFP (a Z shock, no risk dimension) stays isotropic at MU. KNOWN LIMIT: output pass-through
# is ~0 here -- the working-capital financing income to households lifts consumption about
# as much as investment falls; the ROBUST pass-through is to spreads/bond/investment.
MU = 1                                  # TFP grid (isotropic; no risk dimension to resolve)
RISK_MU_VEC = [0, 0, 1, 1, 0, 3, 0]     # risk/TPI [K_D,K_F,P_D,P_F,B_D,s,Z_D]: resolve s + net worth
NW_FLOOR = 0.15                         # Bocola net-worth floor (fraction of n_ss), default corners
RUN_TPI = False                         # append the OMT/TPI activation overlay
TPI_ACTIVATIONS = (0.00, 0.50, 1.00)    # priced activation probabilities to compare


def main():
    t0 = time.perf_counter()
    cal = get_calibration()          # standard Bocola-faithful SS: leverage 4, 200bp, mu_ss~0.02
    cal["nw_floor_frac"] = NW_FLOOR  # Bocola net-worth floor for the anisotropic default corners

    banner("Two-country HANK-GK monetary union: steady state")
    ss = solve_steady_state(cal)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_ss_table(ss, cal)

    banner("TFP shock — recursive projection (Z_D as the 7th state)")
    rules_tfp = solve_tfp(cal, ss, sproc, mu=MU)
    tfp_irf(rules_tfp, cal, ss, sproc)

    banner("Sovereign-risk pass-through — recursive projection (anisotropic s + net worth)")
    rules_risk = solve_recursive(cal, ss, sproc, mu_vec=RISK_MU_VEC)
    impact_table(rules_risk, cal, ss, sproc)
    persistence_irf(rules_risk, cal, ss, sproc)

    if RUN_TPI:
        pct = "/".join(f"{round(a * 100):.0f}%" for a in TPI_ACTIVATIONS)
        banner(f"OMT/TPI activation — recursive projection ({pct})")
        tpi_recursive_experiment.run(cal, ss, sproc, mu_vec=RISK_MU_VEC, activations=TPI_ACTIVATIONS)

    print(f"\nTOTAL  {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
