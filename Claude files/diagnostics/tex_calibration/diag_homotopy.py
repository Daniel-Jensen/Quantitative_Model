# T1.4 PROBE (2026-10-26): WALK THE d0 STAGE (pi = 0, eps 1e-2, SS CONSISTENT WITH IT, 15% FLOOR) FROM THE
# PRE-2026-10-26 CALIBRATION TOWARD APPENDIX B: EACH PARAMETER ALONE (FULL MOVE, ELSE HALF), FROM THE
# OLD SOLUTION. Which appendix values the coarse d0 solve cannot absorb on its own.
#   python3 -u "Claude files/diagnostics/tex_calibration/diag_homotopy.py" > runs/diag_homotopy.log
import os
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
sys.path.insert(0, str(HERE))

import numpy as np

# the calibration before the port (the bank instruments as on 2026-10-25; Bohn root 0.93 at q 0.9106 is
# 4*gamma = 0.275)
OLD = dict(sigma_D=1.0, sigma_F=1.0, frisch_D=2.0, frisch_F=2.0, epsilon_D=6.0, epsilon_F=6.0,
           size_F=8.0, omega_home_D=0.85, epsilon_trade=0.5, foreign_share_D=0.20,
           G_D=0.0, G_F=0.0, delta_b_F=0.056, bohn_annual=0.275,
           beta_inter_D=0.997, beta_inter_F=0.997,
           credit_spread_target_D=0.0025, credit_spread_target_F=0.0025)
GROUPS = (("fiscal", ("G_D", "G_F", "bohn_annual", "delta_b_F")),
          ("trade+size", ("size_F", "omega_home_D", "epsilon_trade", "foreign_share_D")),
          ("markup", ("epsilon_D", "epsilon_F")),
          ("households", ("sigma_D", "sigma_F", "frisch_D", "frisch_F")),
          ("banker beta", ("beta_inter_D", "beta_inter_F")),
          ("wedge", ("credit_spread_target_D", "credit_spread_target_F")))


SINGLES = (("G", ("G_D", "G_F")), ("bohn", ("bohn_annual",)), ("delta_b_F", ("delta_b_F",)),
           ("size", ("size_F",)), ("home bias", ("omega_home_D",)), ("eta", ("epsilon_trade",)),
           ("foreign share", ("foreign_share_D",)), ("markup", ("epsilon_D", "epsilon_F")),
           ("sigma", ("sigma_D", "sigma_F")), ("frisch", ("frisch_D", "frisch_F")),
           ("banker beta", ("beta_inter_D", "beta_inter_F")),
           ("wedge", ("credit_spread_target_D", "credit_spread_target_F")))


def economy(cal0, vals):
    # THE CALIBRATION WITH vals, ITS SS (CONSISTENT WITH THE 1e-2 SMOOTHING) AND ANCHORS.
    from global_projection.steady_state import solve_steady_state
    from global_projection.solver_recursive.state_grid import s_process_params
    from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
    cal = dict(cal0, **vals)
    # trade.ces_price divides by 1 - eta: the Cobb-Douglas point eta = 1 itself is a singularity
    if abs(cal["epsilon_trade"] - 1.0) < 1e-9:
        cal["epsilon_trade"] = 1.0 + 1e-6
    cal["omega_home_F"] = 1.0 - (1.0 - cal["omega_home_D"]) * cal["size_D"] / cal["size_F"]
    cal["omega_home"] = cal["omega_home_D"]
    ss = solve_steady_state(cal, verbose=False)
    sp = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sp)
    return cal, ss, sp


def main():
    # ONE PARAMETER (GROUP) AT A TIME FROM THE OLD CALIBRATION, EACH COLD-STARTED: FULL MOVE, ELSE HALF.
    from calibration.global_projection import get_calibration
    from global_projection.solver_recursive.decision_rules import RuleSet
    from global_projection.solver_recursive import recursive_experiment as RE
    base = get_calibration()
    base.update(nw_floor_frac=0.15, mu_eps=1e-2, book_eps=1e-4, floor_eps=1e-3)
    new = {k: float(base[k]) for k in OLD}
    cal, ss, sp = economy(base, OLD)
    rules0 = RuleSet.from_ss(RE.normal_box(ss, cal, mu=1), ss, cal, n_regimes=2)
    rules0.n_gh = RE.N_GH
    ok, it, w = RE._stage_wc(rules0, cal, ss, sp, (0,), True, "d0 old calibration", False)
    print(f"[start] old calibration: ok {ok}  its {it}  max|F| {w:.2e}", flush=True)
    if not ok:
        return
    for name, keys in SINGLES:
        for frac in (1.0, 0.5):
            vals = dict(OLD, **{k: OLD[k] + frac * (new[k] - OLD[k]) for k in keys})
            cal_k, ss_k, sp_k = economy(base, vals)
            # a cold start from this economy's own SS, as the production ladder does
            trial = RuleSet.from_ss(RE.normal_box(ss_k, cal_k, mu=1), ss_k, cal_k, n_regimes=2)
            trial.n_gh = RE.N_GH
            ok, it, w = RE._stage_wc(trial, cal_k, ss_k, sp_k, (0,), True, f"{name} {frac}", False)
            print(f"[{name}] alone, fraction {frac}: ok {ok}  its {it}  max|F| {w:.2e}  "
                  f"(mu_ss {ss_k['ss_bank_D']['mu_ss']:.5f}, C_D {ss_k['C_D_ss']:.3f})", flush=True)
            if ok:
                break


if __name__ == "__main__":
    main()
