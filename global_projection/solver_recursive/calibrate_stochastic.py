# CALIBRATE THE DETERMINISTIC-SS INSTRUMENTS SO THAT THE STOCHASTIC REST POINT HITS THE TARGETS.
# The solved model does not rest at the deterministic SS: priced default risk marks the D bond
# down, moves net worth and the books, and the rest point (Bocola's generate_irf.m step 1) is
# where every reported number is read. The calibration table's targets are rest-point moments
# (2026-10-19 decision): leverage 8 in both banks, the loan spread r_wc - r_dep = lambda*mu/E[Omega]
# at 5 bp/yr in both (appendix B, 2026-10-26; 100 bp before), and own-sovereign bonds at 14.4% (D)
# and 6.5% (F) of own-bank assets.
# Six targets, six instruments matched one to one (leverage_target, credit_spread_target,
# sov_exposure per country), solved by Broyden's method in tolerance-scaled units; every
# evaluation is a full solve on the pipeline's own grid (S_REFINE nodes in s), so the pipeline
# re-solving at the converged instruments reproduces the rest point the loop hit.
# The loop checkpoints after every evaluation and resumes from the checkpoint; its result is
# the six numbers to paste into calibration/global_projection.py.
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import norm

from calibration.global_projection import get_calibration
from global_projection.steady_state import solve_steady_state
from global_projection.solver_recursive.state_grid import s_process_params, S_COVER_SD
from global_projection.solver_recursive.recursive_main import calibrate_household_anchors
from global_projection.solver_recursive.recursive_experiment import (
    solve_recursive, stochastic_rest_point, read_at)

# (instrument, target, moment, tolerance, lower bound, upper bound), matched one to one
PAIRS = (("leverage_target_D", "srp_leverage_D", "lev_D", 0.005, 6.0, 11.0),
         ("leverage_target_F", "srp_leverage_F", "lev_F", 0.005, 6.0, 11.0),
         # wide: the D bank earns the rest point's sovereign risk premium and is close to slack
         # there at an SS spread of 100 bp (2026-10-20: 2 bp at rest), so the instrument must rise;
         # from 1 bp/yr, so the F instrument can sit near the 5 bp target (2026-10-26)
         ("credit_spread_target_D", "srp_spread_D", "spread_D", 0.1 / 4e4, 1 / 4e4, 0.0150),
         ("credit_spread_target_F", "srp_spread_F", "spread_F", 0.1 / 4e4, 1 / 4e4, 0.0150),
         ("sov_exposure_D", "srp_exposure_D", "exp_D", 1e-4, 0.08, 0.25),
         ("sov_exposure_F", "srp_exposure_F", "exp_F", 1e-4, 0.03, 0.12))
TOL = np.array([p[3] for p in PAIRS])
MAX_EVALS = 12
# the yield the book is revalued at in the repricing-loss row: 9.1% p.a., the 2010-average
# Greek 10y (FRED IRLTLT01GRM156N: 9.09%), the yield that reproduces both table entries
REPRICE_YIELD = 0.091
# the 2010-2011 average quarterly default probability of the exceedance row
CRISIS_PD = 0.091


def rest_point_read(rules, cal, ss, sproc):
    # THE STOCHASTIC REST POINT AND THE FITTED READ THERE.
    S = stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    return S, read_at(rules, cal, ss, sproc, S.copy())


def rest_point_moments(o):
    # THE SIX TARGETED MOMENTS FROM A READ AT THE REST POINT (assets = deposits + net worth).
    a_D, a_F = o["dep_D"] + o["n_D"], o["dep_F"] + o["n_F"]
    return dict(lev_D=a_D / o["n_D"], lev_F=a_F / o["n_F"],
                spread_D=o["r_wc_D"] - o["rdep_D"], spread_F=o["r_wc_F"] - o["rdep_F"],
                exp_D=o["Q_bD"] * o["b_D_D_new"] / a_D, exp_F=o["Q_bF"] * o["b_F_F_new"] / a_F)


def _instruments(cal):
    # THE INSTRUMENT VECTOR IN PAIRS ORDER.
    return np.array([cal[p[0]] for p in PAIRS])


def _targets(cal):
    # THE TARGET VECTOR IN PAIRS ORDER.
    return np.array([cal[p[1]] for p in PAIRS])


def evaluate(x, settings, verbose=True):
    # ONE FULL SOLVE AT INSTRUMENTS x: THE REST-POINT MOMENTS AND THE SOLVE'S VERDICT.
    cal = get_calibration()
    cal["nw_floor_frac"] = settings["NW_FLOOR"]
    for p, v in zip(PAIRS, x):
        cal[p[0]] = float(v)
    ss = solve_steady_state(cal, verbose=False)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    rules = solve_recursive(cal, ss, sproc, mu_vec=settings["RISK_MU_VEC"],
                            rotate=settings["ROTATE_P"], s_refine=settings["S_REFINE"],
                            verbose=verbose)
    _, o = rest_point_read(rules, cal, ss, sproc)
    m = rest_point_moments(o)
    return (np.array([m[p[2]] for p in PAIRS]),
            dict(solve_ok=bool(rules.solve_ok), solve_worst=float(rules.solve_worst),
                 read_resid=float(o["_resid"]), p_ss=float(ss["p_ss"]), p_rest=float(o["p"])))


def _load(path):
    # THE CHECKPOINT, OR A FRESH ONE.
    if path.exists():
        return json.loads(path.read_text())
    return dict(history=[], J=None)


def calibrate(settings, path, max_evals=MAX_EVALS):
    # BROYDEN ON THE SIX INSTRUMENTS IN TOLERANCE-SCALED UNITS, ONE FULL SOLVE PER EVALUATION.
    cal0 = get_calibration()
    t = _targets(cal0)
    lo, hi = np.array([p[4] for p in PAIRS]), np.array([p[5] for p in PAIRS])
    state = _load(path)
    hist = state["history"]
    # J maps scaled instrument steps to scaled moment changes; own moments move ~1:1 with
    # their own instruments, so the identity is the starting guess
    J = np.eye(len(PAIRS)) if state["J"] is None else np.array(state["J"])
    x = _instruments(cal0) if not hist else np.array(hist[-1]["x_next"])
    for k in range(len(hist), max_evals):
        t0 = time.perf_counter()
        print(f"\n=== calibration evaluation {k}: "
              + "  ".join(f"{p[0]}={v:.8g}" for p, v in zip(PAIRS, x)), flush=True)
        m, info = evaluate(x, settings)
        F = (m - t) / TOL
        print(f"  rest point: " + "  ".join(f"{p[2]}={v:.8g}" for p, v in zip(PAIRS, m))
              + f"\n  scaled gaps (|.| <= 1 converges): "
              + "  ".join(f"{g:+.2f}" for g in F)
              + f"   [{time.perf_counter() - t0:.0f}s, solve ok {info['solve_ok']}]", flush=True)
        good = [h for h in hist if not h.get("rejected")]
        if not info["solve_ok"]:
            # an unconverged solve is not a moment: halve the step from the last accepted point
            if not good:
                raise RuntimeError("the solve at the starting instruments did not converge")
            x_prev = np.array(good[-1]["x"])
            hist.append(dict(x=x.tolist(), m=m.tolist(), F=F.tolist(), info=info,
                             x_next=(x_prev + 0.5 * (x - x_prev)).tolist(), rejected=True))
            x = np.array(hist[-1]["x_next"])
            path.write_text(json.dumps(dict(history=hist, J=J.tolist()), indent=1))
            continue
        if good:
            # Broyden's good update on the last accepted step
            du = (x - np.array(good[-1]["x"])) / TOL
            dF = F - np.array(good[-1]["F"])
            if du @ du > 0:
                J = J + np.outer(dF - J @ du, du) / (du @ du)
        done = bool(np.all(np.abs(F) <= 1.0))
        step = np.zeros_like(x) if done else -np.linalg.solve(J, F) * TOL
        x_next = np.clip(x + step, lo, hi)
        hist.append(dict(x=x.tolist(), m=m.tolist(), F=F.tolist(), info=info,
                         x_next=x_next.tolist(), converged=done))
        path.write_text(json.dumps(dict(history=hist, J=J.tolist()), indent=1))
        if done:
            print("\n  CONVERGED: paste into calibration/global_projection.py")
            for p, v in zip(PAIRS, x):
                print(f"    {p[0]} = {v!r}")
            return x, m
        x = x_next
    print("\n  STOPPED at the evaluation budget without converging")
    return x, None


def _duration_years(y_q, delta_b):
    # MACAULAY DURATION OF THE HM PERPETUITY AT A QUARTERLY FLOW YIELD, IN YEARS.
    return (1.0 + y_q) / (y_q + delta_b) / 4.0


def calibration_values(o, cal, sproc, irf):
    # EVERY MODEL NUMBER OF THE CALIBRATION TABLE, BY NAME (RATES ANNUALISED, p^d QUARTERLY).
    # o is the fitted read at the rest point; irf is the dynamic IRF at p^d = 2%/qtr.
    m = rest_point_moments(o)
    dbD, dbF = cal["delta_b_D"], cal["delta_b_F"]
    # HM perpetuity flow yields, quarterly
    yD = dbD * (1.0 - o["Q_bD"]) / o["Q_bD"]
    yF = dbF * (1.0 - o["Q_bF"]) / o["Q_bF"]
    y_det = cal["r_dep_D_target"] + cal["credit_spread_target_D"]
    sd_s = sproc["sigma_s"] / np.sqrt(1.0 - sproc["rho_s"] ** 2)
    z = (np.log(CRISIS_PD / (1.0 - CRISIS_PD)) - sproc["s_star"]) / sd_s
    od = np.asarray(irf["d_loan_D"]) - np.asarray(irf["d_loan_F"])
    sz = cal["size_F"] / cal["size_D"]
    return dict(
        lev_D=m["lev_D"], lev_F=m["lev_F"],
        loan_D=4e4 * m["spread_D"], loan_F=4e4 * m["spread_F"],
        rdep_par=400 * cal["r_dep_D_target"], rdep_rest=400 * o["rdep_D"],
        pd_rest=100.0 / (1.0 + np.exp(-sproc["s_star"])),
        # the AR(1) for s is Gaussian, so its ergodic share above a threshold is exact
        crisis_share=100 * norm.sf(z), crisis_z=float(z),
        pd_box_top=max(100.0 / (1.0 + np.exp(-(sproc["s_star"] + S_COVER_SD * sd_s))),
                       100.0 * cal.get("pd_box_top", 0.0)),
        sov_spread=4e4 * (yD - yF), y_D=4e4 * yD, y_F=4e4 * yF,
        recovery=100 * cal["recovery_rate_D"],
        exp_D=100 * m["exp_D"], exp_F=100 * m["exp_F"],
        sov_ct1=o["Q_bD"] * o["b_D_D_new"] / o["n_D"],
        foreign_D=100 * o["b_D_F_new"] / o["Bp_D"],
        foreign_F=100 * o["b_F_D_new"] / (o["b_F_D_new"] + sz * o["b_F_F_new"]),
        dur_rest=_duration_years(yD, dbD), dur_det=_duration_years(y_det, dbD),
        dur_45=_duration_years(0.045 / 4.0, dbD),
        debt_mkt=100 * o["Q_bD"] * o["Bp_D"] / (4.0 * o["Y_D"]),
        debt_face=100 * o["Bp_D"] / (4.0 * o["Y_D"]),
        reprice=100 * (o["Q_bD"] - dbD / (dbD + REPRICE_YIELD / 4.0)) * o["b_D_D_new"] / o["n_D"],
        deposits=100 * o["A_D"] * o["P_CES_D"] / (4.0 * o["Y_D"]),
        overdraft_q0=float(od[0]), overdraft_peak=float(np.max(od)),
        p_rest=float(o["p"]), Q_bD_rest=float(o["Q_bD"]), mu_D=float(o["mu_D"]),
        mu_F=float(o["mu_F"]), read_resid=float(o["_resid"]))


# (block, moment, data, value key, unit, role): the user's table, one model number per row
TABLE_ROWS = (
    ("Banks", "Leverage, GR fin. accounts, 2002-07", "7.98", "lev_D", "x", "Target"),
    ("Banks", "Leverage, GR credit book / CT1, 2010", "5.5", "lev_D", "x", "Narrower measure"),
    ("Banks", "Leverage, F bank", "--", "lev_F", "x", "Target"),
    ("Banks", "Intermediation wedge r_wc - r_dep, D", "100bp", "loan_D", "bp", "Target"),
    ("Banks", "Intermediation wedge r_wc - r_dep, F", "100bp", "loan_F", "bp", "Target"),
    ("Banks", "Lending spread, 2003-07, GR", "290bp", "loan_D", "bp", "Upper bound"),
    ("Banks", "Lending spread, 2003-07, DE", "180bp", "loan_F", "bp", "Upper bound"),
    ("Banks", "Euribor-Eurepo, 2002-2012", "8bp", "loan_D", "bp", "Not matched"),
    ("Banks", "Real deposit rate, 1999-2007 (parameter)", "1.11%", "rdep_par", "%",
     "Set (Bocola)"),
    ("Banks", "Real deposit rate at the rest point, D", "1.11%", "rdep_rest", "%", "Set (Bocola)"),
    ("Sovereign risk", "Default probability, GR 1y CDS, 2008Q3", "0.092%", "pd_rest", "%",
     "Set directly"),
    ("Sovereign risk", "Share of quarters with p^d above the 2010-11 avg", "9.1%",
     "crisis_share", "%", "Validation"),
    ("Sovereign risk", "Bond-implied hazard, 1999-2007", "0.23%", "pd_rest", "%", "Upper bound"),
    ("Sovereign risk", "Bond-implied hazard, basis-adjusted", "0.02-0.05%", "pd_rest", "%",
     "Lower bound"),
    ("Sovereign risk", "GR-DE sovereign spread y_D - y_F", "25.3bp", "sov_spread", "bp",
     "Validation"),
    ("Sovereign risk", "PSI recovery, face value; NPV", "46.5%; 35-41%", "recovery", "%",
     "Bracket"),
    ("Exposure", "GR sovereign / GR bank assets", "14.4%", "exp_D", "%", "Target"),
    ("Exposure", "DE sovereign / DE bank assets", "6.5%", "exp_F", "%", "Target"),
    ("Exposure", "GR sovereign / GR Core Tier 1", "2.39", "sov_ct1", "x", "Limitation"),
    ("Exposure", "Share of the D stock held by F banks", "12.7%; 39.6%", "foreign_D", "%",
     "Bracket"),
    ("Exposure", "Share of the F stock held by D banks", "~0", "foreign_F", "%", "Set"),
    ("Exposure", "Duration, GR book, at the rest-point yield", "4.4 years", "dur_rest", "years",
     "Validation"),
    ("Exposure", "Duration at the deterministic-SS yield", "4.4 years", "dur_det", "years",
     "Validation"),
    ("Exposure", "Duration at 4.5%", "4.4 years", "dur_45", "years", "Validation"),
    ("Exposure", "Bank-held GR debt / GDP, market value", "40%", "debt_mkt", "%", "Not matched"),
    ("Exposure", "Bank-held GR debt / GDP, face value", "40%", "debt_face", "%", "Not matched"),
    ("Exposure", "Repricing loss / GR CT1, book revalued at 9.1%", "40%", "reprice", "%",
     "Validation"),
    ("Households", "Liquid deposits / GDP, Greece", "84%", "deposits", "%", "Validation"),
    ("Transmission", "GR-DE overdraft-spread rise, impact of p^d = 2%", "108bp", "overdraft_q0",
     "bp", "Upper bound"),
    ("Transmission", "GR-DE overdraft-spread rise, peak", "108bp", "overdraft_peak", "bp",
     "Upper bound"),
)


def main():
    # RUN THE CALIBRATION WITH THE PIPELINE'S OWN GRID SETTINGS (run.py), CHECKPOINTED.
    import run
    settings = dict(NW_FLOOR=run.NW_FLOOR, RISK_MU_VEC=run.RISK_MU_VEC, ROTATE_P=run.ROTATE_P,
                    S_REFINE=run.S_REFINE)
    out = Path(run.RESULTS) / "CALIBRATION"
    out.mkdir(parents=True, exist_ok=True)
    print(f"calibrating at the stochastic rest point, settings {settings}", flush=True)
    calibrate(settings, out / "broyden_state.json")


if __name__ == "__main__":
    main()
