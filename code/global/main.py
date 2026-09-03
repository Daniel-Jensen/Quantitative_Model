# ENTRY POINT: GLOBAL CHEBYSHEV COLLOCATION SOLVER FOR THE WHOLE MODEL.
# The two-country model is solved GLOBALLY as recursive decision rules on a Smolyak grid
# (solver_recursive/) -- there is NO perfect-foresight machinery. main() solves the steady
# state once, then runs each experiment as a full collocation solve: a TFP shock, the
# Bocola sovereign-risk pass-through, and (RUN_LTRO) the LTRO-backstop activation sweep.
#
# STATE (10): [K_D, K_F, P_D, P_F, b_DD, b_DF, b_FD, V_dep, s, Z_D]
# UNKNOWNS: every stored rule is a collocation unknown -- 19 per point per regime (the 13
#   market-clearing/Euler rules [N_D, N_F, Kp_D, Kp_F, rdep_D, rdep_F, p, Q_bD, b_DF,
#   Q_bF, b_FD, A_D, A_F] plus the 6 that used to be READ OFF a frozen continuation
#   [alpha_D/F, C_D/F, r_wc_D/F], which now carry Bocola's identity residual
#   log(guess/implied)). 798 unknowns on the coarse grid, 3610 s-refined.
# The D sovereign is MARKET-CLEARED (both banks' bond FOCs are residuals, b_DD = B' - b_DF,
# Q_bD is the price that clears) and deposits clear on ONE UNION market (both households'
# Eulers are residuals, V_dep carries the cross-border position).
# The two countries have DIFFERENT MASSES (size_F/size_D = 8) but an identical PER-CAPITA
# steady state; see calibration.py and CLAUDE.md.
#
# THE SOLVER IS solver_recursive/collocation.py -- one damped Newton on the whole
# coefficient vector (Bocola's residual_model.m + parsolve.m). Time iteration survives
# only as the warm start that puts the Newton inside its basin; it is not convergent
# here (0.990 per sweep on the franchise-value mode).
#
# RUNTIME, MEASURED. The cost rule is one dense finite-difference Jacobian per Newton
# step = (unknowns + 1) residual evaluations, each costing points x regimes x ~2.5 ms,
# so a solve scales as (points x regimes)^2:
#     stage                          unknowns   one Jacobian   solve
#     risk, coarse (21 pts, 2 reg)        798        ~1.5 min   ~7 min
#     risk, s-refined (95 pts, 2 reg)    3610         ~29 min   ~86 min
#     backstop, coarse (21 pts, 4 reg)   1596         ~5.6 min   ~30-55 min
#     backstop, s-refined (95 pts, 4 reg) 7220        ~114 min   ~8-10 h
# ~100 min end to end at S_REFINE = 5 WITHOUT the overlay. The eight-point activation
# sweep adds ~3 h on the coarse grid (LTRO_S_REFINE = 0, the default) because the
# phi-independent baseline is solved once and reused; it would be ~70 h refined, which
# is why the sweep is coarse and the refined grid is reserved for two or three points.
#
# READ THE SPREADS CORRECTLY -- they are two different objects:
#   BANK CREDIT spread lambda_K*mu/alpha is the wedge the leverage constraint puts on
#     capital. It is the STEADY-STATE calibration target and, through the working-capital
#     rate r_wc = rdep + lambda*mu/E[Om], it is the ONLY channel into output under GHH.
#     It is identically zero once mu hits the KKT switch (measured q4 at the headline
#     shock). NB it is a NET wedge: dynamic_irf prints the deposit-rate leg beside it,
#     because the two move in opposite directions and the netting is what sets the
#     output response.
#   SOVEREIGN spread y_D - y_F comes out of the bond Euler and carries no mu, so it
#     persists for as long as p^d is elevated. This is the one the FIGURE plots.
# Measurement record: CLAUDE.md, and docs/recursive_9state_findings.md (partly
# superseded -- see the note at its head).
import time

from config.calibration import get_calibration
from config.steady_state import solve_steady_state
from solver_recursive.state_grid import s_process_params
from solver_recursive.recursive_main import calibrate_household_anchors, ss_state
from solver_recursive.accuracy import accuracy_report
from solver_recursive.recursive_experiment import (
    solve_tfp, tfp_irf, solve_recursive, impact_table, persistence_irf, dynamic_irf,
    s_from_pd)
from solver_recursive.output_decomposition import (simulate, s_decay_path,
                                                   decompose_output, active_channels,
                                                   decompose_bond_price, BOND_CHANNELS)
from solver_recursive.recursive_experiment import (stochastic_rest_point,
                                                   report_rest_point)
from solver_recursive import ltro_experiment
from reporting.prints import banner, print_ss_table
from reporting.plots import (plot_risk_irf, plot_tfp_irf, plot_output_decomposition,
                             plot_bond_decomposition)

# CONFIG. NW_FLOOR is Bocola's net-worth floor (his N_tom = max(., 0.65)): it keeps the
# deep default-regime corners feasible so the d=1 fit does not poison the global basis.
MU = 1                                  # TFP grid (isotropic; no risk dimension to resolve)
# BASE SMOLYAK LEVEL PER STATE for the COARSE rung of the ladder. None = isotropic mu=1,
# 21 points. Refining the risk dimension through mu_vec was the old plan and is the wrong
# tool: raising ONE dimension raises the GLOBAL Smolyak budget, so [1,...,2,1] costs 165
# points, and every cheaper vector freezes CAPITAL at a single node. S_REFINE below does
# it properly with a tensor factor instead. Left here for anisotropy that is NOT about s.
RISK_MU_VEC = None                      # isotropic mu=1 base grid for the coarse stage
# DENSE CHEBYSHEV NODES TENSORED ONTO THE s DIMENSION (state_grid.SmolyakGrid refine).
# All the curvature in this model is the logistic p^d(s); the other nine states are
# near-linear. The coarse ladder solves first and SEEDS this grid, so the refinement is a
# continuation step, not a cold start. Measured relative RMS error on this curvature
# profile: isotropic mu=1 21pts 1.9e-1; isotropic mu=2 221pts 3.9e-2; s_refine=5 95pts
# 2.5e-2; s_refine=9 171pts 1.1e-3.
# 5 = degree 4 in s (95 points, ~85 min); 9 = degree 8, Bocola's own resolution
# (171 points, ~4 h, walked up through 5 by the ladder). 0 or 1 = no refinement.
S_REFINE = 5
# THE BACKSTOP OVERLAY RUNS ON THE COARSE GRID BY DEFAULT, and that is a deliberate
# trade of resolution for the SHAPE of the activation curve. One refined (95-point)
# four-regime solve is ~8-10 h, so the eight-point sweep below would be ~70 h; the same
# sweep coarse is ~3 h. Every point shares one grid and one baseline, so the CURVE --
# which is what the experiment is for -- is far better measured than any single level.
# Set LTRO_S_REFINE = 5 with two or three activations for publication levels instead.
LTRO_S_REFINE = 0                       # 0 = coarse (21 pts); 5 = refined (95 pts)
NW_FLOOR = 0.15                         # Bocola net-worth floor (fraction of n_ss), default corners
ROTATE_P = False                        # eigenbasis box: right in theory, measures worse (see solve_recursive)
# THE LTRO BACKSTOP OVERLAY. One four-regime collocation solve per activation -- phi is
# a per-experiment scalar, not a state -- reusing the phi-independent baseline across all
# of them. The headline read is the NEVER-FIRED path: the realisation on which the
# facility is announced and never drawn, which is where the OMT fact lives.
# phi = 0 .. 0.7 in 10pp steps. The range stops at 0.7 because beyond it the facility
# drives mu to zero AT THE REST POINT and the KKT kink takes the identification with it:
# at phi = 1 the solve does not reach the acceptance floor and the fitted-vs-exact output
# bracket straddles zero. 0.7 is the last point that is a number rather than a limit.
RUN_LTRO = True                         # append the LTRO-backstop activation overlay
LTRO_ACTIVATIONS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
ACCURACY_T = 1200                       # simulated periods for the Euler-error report
DECOMP_T = 25                           # quarters in the output-decomposition figure
TFP_SHOCK = 0.01                        # one-off TFP shock: +1% to Z_D, decaying at rho_z = 0.9
RISK_SHOCK_PD = 0.0198                  # one-off risk shock: p^d jumps 0.10% -> 1.98%, decays at rho_s


def main():
    t0 = time.perf_counter()
    # Bocola's leverage 5 / exposure 7.6% / recovery 0.45; the credit spread is
    # 100 bp/yr, NOT his 8 -- at 8 the rest point sits on the KKT kink and the level
    # of the output response is not identified. See config/calibration.py.
    cal = get_calibration()
    cal["nw_floor_frac"] = NW_FLOOR  # Bocola net-worth floor for the deep default corners

    banner("Two-country HANK-GK monetary union: steady state")
    ss = solve_steady_state(cal)
    sproc = s_process_params(cal)
    calibrate_household_anchors(cal, ss, sproc)
    print(f"  solved in {time.perf_counter() - t0:.0f}s")
    print_ss_table(ss, cal)

    # NB the TFP experiment needs NO rest-point rebasing, unlike the risk one. It is
    # solved at pi == 0, and at pi == 0 the deterministic SS *is* the model's rest point
    # -- verified: 200 no-shock quarters, every state within 0.0001%, mu pinned at
    # 0.001001, max|F| ~ 1e-7. The SS/rest-point gap is priced risk and nothing else.
    banner("TFP shock — global collocation (Z_D as the TFP state)")
    rules_tfp = solve_tfp(cal, ss, sproc, mu=MU)
    tfp_path = tfp_irf(rules_tfp, cal, ss, sproc, dz=TFP_SHOCK)
    print(f"  figure -> {plot_tfp_irf(tfp_path, note=f'{TFP_SHOCK:.0%} shock to Z_D')}")

    banner("Sovereign-risk pass-through — global collocation (10-state, Newton)")
    rules_risk = solve_recursive(cal, ss, sproc, mu_vec=RISK_MU_VEC, rotate=ROTATE_P,
                                 s_refine=S_REFINE)
    # WHERE THE MODEL ACTUALLY RESTS. print_ss_table above reports the DETERMINISTIC
    # steady state; the solved rules price risk and rest somewhere else, and every IRF
    # below is read against THAT point. Print both or the SS table silently describes a
    # state the model never visits.
    S_rest = report_rest_point(rules_risk, cal, ss, sproc)
    impact_table(rules_risk, cal, ss, sproc)
    # states HELD at the rest point -- the pure shock-persistence channel
    persistence_irf(rules_risk, cal, ss, sproc, pd_shock=RISK_SHOCK_PD)
    # states EVOLVE -- capital and bank net worth accumulate. NB Bocola's Table 5 is
    # NOT this object: it is a cumulated quarterly GROWTH gap x400 over an 8-quarter
    # estimated shock sequence, so its level equivalent is -0.26/-0.36/-0.38%. The
    # single-shock targets dynamic_irf prints are the like-for-like ones.
    risk_path = dynamic_irf(rules_risk, cal, ss, sproc, pd_shock=RISK_SHOCK_PD, T=25)
    print(f"  figure -> {plot_risk_irf(risk_path, note=f'p-d shock to {100*RISK_SHOCK_PD:.2f}%/qtr, states evolving')}")

    # OUTPUT DECOMPOSITION, read off the SAME converged rules -- no extra solve.
    # Under GHH the production function and the labour FOC are identities, so d log Y
    # splits EXACTLY into TFP, capital, the relative price and the two legs of the
    # working-capital wedge (credit spread and deposit rate), with whatever the FOC
    # misses at the read carried as an explicit `residual` rather than absorbed into a
    # channel. This is the figure that shows the netting the console table reports:
    # the spread leg and the deposit-rate leg point in opposite directions, and their
    # sum is the output response.
    banner("Decompositions — which channels produce the response")
    # BOTH PATHS START AT THE MODEL'S OWN REST POINT, not the deterministic SS: with
    # risk priced the two differ, and starting at the SS decomposes the response along
    # the transition between them rather than around the point the economy inhabits.
    s_path = s_decay_path(sproc, s_from_pd(RISK_SHOCK_PD), DECOMP_T)
    sim = simulate(rules_risk, cal, ss, sproc, s_path, S_init=S_rest)
    ref = simulate(rules_risk, cal, ss, sproc,
                   [sproc["s_star"]] * DECOMP_T, S_init=S_rest)   # same rules, no shock
    dec = decompose_output(sim, ref, cal)
    for k, lab in active_channels(dec):
        print(f"    {lab:<22s} impact {dec[k][0]:+8.4f}%  "
              f"({4 * dec[k][0]:+8.4f}% annualised)")
    print(f"    {'TOTAL':<22s} impact {dec['total'][0]:+8.4f}%  "
          f"({4 * dec['total'][0]:+8.4f}% annualised)")
    fig = plot_output_decomposition(
        [(f"p-d shock to {100*RISK_SHOCK_PD:.2f}%/qtr", dec)], active_channels(dec),
        note="d log Y_D split by the production function and the GHH labour FOC; "
             "the two wedge legs are a symmetric (Shapley) split, so the channels "
             "sum to the total identically")
    print(f"  figure -> {fig}")

    # WHY THE SOVEREIGN REPRICES, from the D bank's own FOC. Bocola's Table 4 splits
    # the EXCESS RETURN into a risk premium and a liquidity premium; this splits the
    # PRICE, and adds the two legs his table takes as given -- the discount rate and
    # the continuation price -- so the legs sum to the observed repricing exactly.
    bdec = decompose_bond_price(sim, ref, cal)
    print("\n  Q_bD decomposition (%, deviation from the no-shock path)")
    for k, lab in BOND_CHANNELS:
        print(f"    {lab:<36s} impact {bdec[k][0]:+8.4f}%")
    print(f"    {'TOTAL':<36s} impact {bdec['total'][0]:+8.4f}%")
    print(f"  figure -> {plot_bond_decomposition(bdec, list(BOND_CHANNELS), note='the D bank first-order condition, split leg by leg')}")

    banner("Solution accuracy — Euler errors on the ergodic set")
    # simulate from the REST POINT, not the deterministic SS: the 200-period burn-in
    # would get there anyway, but starting on the ergodic set means the whole sample
    # measures accuracy where the model lives rather than partly along its transition.
    accuracy_report(rules_risk, cal, ss, sproc, S_rest,
                    T=ACCURACY_T, label="sovereign-risk rules")

    if RUN_LTRO:
        pct = "-".join(f"{round(a * 100):.0f}" for a in LTRO_ACTIVATIONS[:1]
                       + LTRO_ACTIVATIONS[-1:])
        banner(f"LTRO backstop by activation probability ({pct}%, "
               f"{len(LTRO_ACTIVATIONS)} points)")
        ltro_experiment.run(cal, ss, sproc, mu_vec=RISK_MU_VEC,
                            activations=LTRO_ACTIVATIONS, pd_shock=RISK_SHOCK_PD,
                            s_refine=LTRO_S_REFINE, accuracy=False)

    print(f"\nTOTAL  {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
