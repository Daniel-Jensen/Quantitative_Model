# CALIBRATE THE BANK FRICTION AGAINST THE *STOCHASTIC* REST POINT, NOT THE DETERMINISTIC SS.
# blocks.bank.calibrate_bank_targets solves lambda and omega_ent analytically from a
# leverage and a credit-spread target, but it does so at the DETERMINISTIC steady state --
# where the D sovereign is priced at Q_bD_ss = 0.9459 and no default is priced at all.
# The solved model does not rest there. Its ergodic mean p^d is 0.405% against 0.100% at
# s*, so the bank permanently holds a bond marked ~0.906, its divertable base is smaller,
# and the SAME lambda delivers a materially tighter constraint: measured 32.8 bp at the
# joint-stage rest point against the 8.0 bp the calibration was told to hit (275 bp before
# the union-market and bond-market fixes).
#
# So the deterministic target is an INSTRUMENT, not the calibration. This module inverts
# the map numerically: secant-iterate on cal["credit_spread_target_*"] until the spread
# READ OFF THE CONVERGED RULES at the steady-state point equals the wanted one. Each
# evaluation is a full d=0 -> d=1 homotopy -> joint solve, so this is minutes per step,
# not seconds -- it is a calibration run, not something main.py does every time. The
# result is a pair of numbers to paste into calibration.py.
#
# MEASURED 2026-08-26: credit_spread_target is a DEAD instrument here. Driving it from
# 8.04 bp to 0.04 bp moved the stochastic rest point only 32.8 -> 24.6 bp and left lambda
# at 0.199980 throughout, because lambda = alpha/theta_target and alpha -> Omega*(1+rdep)
# as mu -> 0, i.e. lambda is pinned by the LEVERAGE target. The spread target enters only
# through mu inside the franchise fold, which is second order at small mu -> a flat map,
# and the secant stalls. The live instruments are leverage_target (which sets lambda
# directly) and f (which sets the fold, hence alpha).
#
# Leverage is REPORTED at every step but not targeted: with one instrument only one
# moment can be hit, and the spread is the one the risk channel runs through. If the
# stochastic leverage drifts materially from its target the second instrument
# (leverage_target, which moves omega_ent) has to come in too -- the report says so.
import time

import numpy as np

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import s_process_params
from solver_recursive.recursive_main import calibrate_household_anchors, ss_state
from solver_recursive.recursive_experiment import solve_recursive, read_at, _spread_bp

TARGET_BP = 8.0          # wanted ANNUALISED spread at the stochastic rest point
TOL_BP = 0.3             # accept within this many bp
MAX_STEPS = 6


def rest_point(spread_target=None, nw_floor=0.15, verbose=False, **over):
    # SOLVE AT A CANDIDATE CALIBRATION, RETURN THE STOCHASTIC REST POINT.
    # `over` sets any calibration key for BOTH countries (leverage_target, f, ...):
    # credit_spread_target turned out to be the WRONG instrument -- see the module note.
    cal = get_calibration()
    cal["nw_floor_frac"] = nw_floor
    if spread_target is not None:
        cal["credit_spread_target_D"] = cal["credit_spread_target_F"] = spread_target
    for k, v in over.items():
        cal[f"{k}_D"] = cal[f"{k}_F"] = v
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    # COARSE GRID ON PURPOSE (s_refine=0). This module sweeps an instrument across
    # several full solves to read a steady-state moment; the s-refinement changes the
    # ergodic rest point only in the third digit and costs ~70 min per solve.
    rules = solve_recursive(cal, ss, sproc, mu_vec=None, verbose=verbose, s_refine=0)
    o = read_at(rules, cal, ss, sproc, ss_state(ss, cal, sproc).copy())
    # leverage at the rest point, on ACTUAL net worth (bank.py's convention)
    lev = (o["dep_D"] + o["n_D"]) / max(o["n_D"], 1e-9)
    return dict(spread_bp=_spread_bp(o, cal), mu=o["mu_D"], lev=lev,
                Q_bD=o["Q_bD"], n_D=o["n_D"], resid=o["_resid"],
                lam=cal["lambda_K_D"], om=cal["omega_ent_D"])


def calibrate(target_bp=TARGET_BP, tol_bp=TOL_BP, max_steps=MAX_STEPS):
    # SECANT ON THE DETERMINISTIC TARGET UNTIL THE STOCHASTIC REST POINT HITS target_bp.
    cal0 = get_calibration()
    x0 = cal0["credit_spread_target_D"]
    # second point from the ~proportional first guess: the map is monotone and close to
    # linear through the origin, so target * (wanted / measured) lands near the root
    print(f"  target {target_bp:.1f} bp at the STOCHASTIC rest point "
          f"(tol {tol_bp:.1f} bp, <= {max_steps} solves)\n")
    print("   step  det.target(bp/yr)   stoch.spread(bp)    mu       leverage   Q_bD   resid")
    hist = []
    t0 = time.perf_counter()
    r0 = rest_point(x0)
    hist.append((x0, r0))
    print(f"   {0:4d} {4e4*x0:16.3f} {r0['spread_bp']:18.2f} {r0['mu']:9.5f}"
          f" {r0['lev']:10.3f} {r0['Q_bD']:8.4f} {r0['resid']:.0e}")
    if abs(r0["spread_bp"] - target_bp) <= tol_bp:
        return _finish(hist, target_bp, t0)
    x1 = x0 * target_bp / max(r0["spread_bp"], 1e-9)
    for k in range(1, max_steps):
        r1 = rest_point(x1)
        hist.append((x1, r1))
        print(f"   {k:4d} {4e4*x1:16.3f} {r1['spread_bp']:18.2f} {r1['mu']:9.5f}"
              f" {r1['lev']:10.3f} {r1['Q_bD']:8.4f} {r1['resid']:.0e}")
        if abs(r1["spread_bp"] - target_bp) <= tol_bp:
            break
        f0, f1 = r0["spread_bp"] - target_bp, r1["spread_bp"] - target_bp
        if abs(f1 - f0) < 1e-9:
            print("   secant stalled (flat map)"); break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x2 = float(np.clip(x2, 1e-6, 5e-3))        # keep the target physically sane
        x0, r0, x1 = x1, r1, x2
    return _finish(hist, target_bp, t0)


def _finish(hist, target_bp, t0):
    # REPORT THE CALIBRATION AND WHAT IT IMPLIES FOR calibration.py.
    x, r = hist[-1]
    print(f"\n  CONVERGED in {len(hist)} solves, {time.perf_counter()-t0:.0f}s")
    print(f"  credit_spread_target_D/F = {x:.8f}   ({4e4*x:.3f} bp/yr deterministic)")
    print(f"    -> stochastic rest point {r['spread_bp']:.2f} bp, mu = {r['mu']:.6f}")
    print(f"    -> implied lambda = {r['lam']:.6f}, omega_ent = {r['om']:.6f}")
    print(f"    -> leverage at the rest point {r['lev']:.3f}")
    if abs(r["lev"] - get_calibration()["leverage_target_D"]) > 0.25:
        print("    NOTE: stochastic leverage is off its target by more than 0.25 -- one")
        print("          instrument cannot hit both moments; bring leverage_target in too.")
    return x, r


if __name__ == "__main__":
    calibrate()
