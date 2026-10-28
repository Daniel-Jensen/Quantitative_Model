# Calibration worklog, branch `calibration` (global model), 2026-10-19 onwards

Working draft of the PROGRESS.md entry; it moves there when the branch is committed.
(Recreated 2026-10-22 after a reboot wiped the scratchpad copy.)

## 2026-10-20: The 2010 calibration table, an asymmetric steady state with p free, rest-point targets

- **The calibration table** (user, 2026-10-19), global model only: leverage 8 (GR financial
  accounts 2002-07: 7.98) in BOTH banks, the intermediation wedge r_wc - r_dep at 100 bp/yr in
  both, own-sovereign bonds at 14.4% (D) and 6.5% (F) of own-bank assets (EBA 2010), F banks
  hold 20% of the D stock, D banks 0.5% of the F stock (~0 in the data). The six targets
  (leverage, wedge, exposure, per country) are hit at the STOCHASTIC REST POINT by six
  deterministic-SS instruments (`calibrate_stochastic.py`, Broyden on full fine-grid solves,
  tolerances 0.005 / 0.1 bp / 0.01 pp). Kept: r_dep 1.2%, p^d 0.1%/qtr at rest, recovery 45%,
  delta_b 0.056 (4.1-year duration at the 2.2% yield), sigma_s 0.4455, rho_s 0.95.
- **The steady state is no longer symmetric, and p is free.** D pays net coupons abroad (F holds
  a fifth of its debt, D almost none of F's), so the current account needs a D trade surplus,
  and at a trade elasticity of 0.5 the trade balance moves only ~0.012 of output per unit of log
  p: p_ss = 1.151. Four steady-state blocks assumed p = 1 and are generalised: the stage-1 solver
  (it asserted only the capital residuals and never converged in p: now the capital FOC in
  closed form and a bracketed root in p), consumption in the external balance (now the
  national-income identity P_CES*C = Y - I - G + NFI, not Y - I - G), chi (w/P_CES, not w) and the
  deposit clearing (basket units, beta_F solved from F's own market). The government block is
  built after the debt stocks are solved (it kept the old stock as the Bohn anchor). The debt
  stocks are solved from the exposure targets (a 2x2 linear system), not hard-coded.
- **Unequal bank spreads** (the rest-point targets make the instruments differ by bank): the
  foreign bond is valued at its issuer's price and earns its issuer's spread in the bank block,
  omega_ent absorbs that leg, and the cross-border friction's reference point is set so the
  calibrated books are optimal at the SS (`_friction_references`; equal spreads nest exactly).
- **Ported from local commit 0dd0732 (user-approved):** the cross-border bond FOCs value the
  foreign payoff at p'/p (the paper's eqs. bondD-foreign / bondF-foreign; main omitted it) and
  the least-squares fallback for a non-finite Newton step.
- **Fixed:** the Bohn gamma_tau is solved at the SS bond price, not the risk-free one (the
  realised debt root was 0.9294, not 0.93); F's banker kernel used D's consumption floor.
- **Solver (each change user-approved after a diagnosis, 2026-10-20/21).**
  1. *Box-edge corner.* At the capital edges investment almost exactly replaces depreciation,
     so next-period capital lands on the box edge, where the hard clip made the residual
     non-differentiable. `state_grid.CLIP_EPS = 0.05`: outside states are read through a C1
     saturation reaching at most 5% of a half-width past the edge; inside, bit-identical.
  2. *Cross-border books in levels.* The D bank's F book (0.5% of F debt) sits at its
     no-short-sale corner in default; a log unknown's Jacobian column vanishes there. `b_DF`/`b_FD`
     are fitted in levels (`decision_rules.LOG_RULES`).
  3. *KKT smoothing* mu = (z + sqrt(z^2 + eps^2))/2, `cal["mu_eps"] = 1e-2`, with the SS bank block
     solving the same smoothed constraint (`bank.ic_slack_factor`), so pi = 0 still nests the SS.
  4. *Smoothing ladder* (`MU_EPS_START`, `_walk_mu_eps`, `solve_no_default`); joint stage and polish
     without time-iteration sweeps.
  5. *Levenberg-Marquardt fallback* in `collocation.parsolve`, used only when the damped Newton
     step fails or stalls (bit-identical otherwise).
- **Two boxes, built and switched off** (`recursive_experiment.TWO_BOXES = False`, user decision
  2026-10-21): aligned normal box + coarse post-default rule set with a blended continuation. On the
  coarse grid the aligned box's cold start does not converge and the post-default stage plateaus at
  7.8e-3; production runs on the V5 box, post-default states read at its edge.

## 2026-10-21/22: The s-refinement on the axis

At leverage 8 the tensor refinement (115 points: 5 nodes in s times the whole sparse grid) has no
collocation solution near the coarse one. Measured: the current solver crawls (Newton steps cut to
1/16-1/32); a Levenberg-Marquardt variant (best of the lambda ladder whenever Newton cuts its step)
takes 37-60x first steps but stalls at max|F| ~3e-3 on both the 3- and 5-node grids; the scaled
Jacobian's smallest singular value is 5e-6 against 4.2 (forward, backward and central differences
agree); each stalled node has a unique root with the continuation frozen (30 starts), so the fold is
global; a wider KKT smoothing (2e-2) breaks even the coarse grid; a lower box top (p^d 3% instead of
4.8%/qtr) crawls the same way. The residual sits where high risk meets the bank-state edges: D-bank
net worth on the 15% floor in the default regime, mu_D switching between 0 and 0.4-0.75 from node to
node in the normal regime, deposit rates of -300 to -580 bp/yr.

User decision: `state_grid._axis_refine` (`refine=(dim, m, "axis")`): the sparse grid plus the higher
Smolyak levels of s on its axis alone, a downward-closed index set, so the basis stays square and
nested; degree 4 in s at 27 points (31 at 9 nodes), no new s x state cross terms.
`recursive_experiment.S_REFINE_MODE = "axis"`; the tensor mode stays available. Solves in ~6 full
Newton steps. Test: `test_state_grid.test_axis_refine`.

## 2026-10-21/22: The rest point, the instrument walk

- At the target-valued instruments the D bank is slack at the stochastic rest point (mu_D ~1e-4,
  D loan spread ~2 bp): it earns the sovereign risk premium (net worth 5-9% above the SS), which
  relaxes its constraint. On the refined grid the rest point lies on the V5 capital edge.
- The deterministic SS does not exist above a D SS-spread instrument of ~240 bp at f_D = 0.08
  (`steady_state._check_targets`: fold-blocked). The SS moves ~0.8% of capital per 10 bp, and the SS
  cold start fails at 200 bp (coarse joint stage stops at 1.6e-4).
- Continuation path (user decision): a coarse walk 100 -> 150 bp in 10 bp rungs converges (8-23
  Newton steps per rung) with the box's bank stocks held at the 100 bp SS. Rest-point D loan spread
  (fitted / exactly cleared): 100 bp 2.6/2.0, 110 7.5/4.2, 120 29.7/9.3, 130 77.5/31.3,
  140 101.7/163.7, 150 74.9/59.9 -- steep and non-monotone. At 140 bp: F spread 97.3 bp, D leverage
  7.69, D exposure 14.04%.
- The stochastic rest point is not unique: from the SS, a neighbouring rung's rest point and random
  interior starts, the no-shock simulation settles at different states (at 150 bp: D spread at rest
  75, 159, 53, 4, 0.1 bp; D leverage 7.7-11.5); the slack ones sit on the box walls. The slowest
  eigenvalue of the linearised no-shock map is 0.987-0.994 (cross-border deposit position and the
  D-bank balance sheet). A 5x stronger debt-elastic premium (kappa_nfa 0.05) does not fix it.
  Cause: the KKT switch across the bank states; the quadratic rules read "slack" where the exactly
  cleared map binds (22 vs 115 bp at the capital edge).
- Remedies tested and rejected (user asked for them): the D-bank block at Smolyak level 2 (59 points)
  stalls the same way (normal regime alone: LM to max|F| ~4e-3, Newton steps 1/128); an
  always-binding constraint gives negative multipliers at 8 of 25 coarse nodes (mu_D down to -0.34),
  so the D bank is genuinely slack over a third of the box (switch removed).
- User decision: binding at rest, with the rest point defined as the one reached from the
  deterministic SS (the code's existing convention); the D spread instrument found by a 1-D search
  near 135-140 bp, leverage and exposure instruments nudged; non-uniqueness and the fitted/exact gap
  documented.
- Implemented: `recursive_experiment.WALK_FROM`/`WALK_STEP`, `_economy` (each economy's SS, its own
  s/Z process -- `z_star = Z_ss_D` moves with the instruments -- and household anchors),
  `_walk_start`, `_walk_instruments` (re-seeds each rung, halves a failed rung); `solve_recursive`
  solves the WALK_FROM economy cold and walks to the requested instruments; `build_state_box(at=...)`
  moves only the bank-stock centres (bands keep the current SS widths; s and Z_D follow the economy;
  `test_box_at`); the box spec travels with the rules (`RuleSet.box`) for the refinement and the
  TPI's reused baseline.
- IRF figures shade the quarters from the first box-edge read on (`dynamic_irf` "at_wall",
  `plots._shade_wall`).

## 2026-10-23: The 0.02% priced default probability, and the bond price's sensitivity to the s rule

- Scripts (`Claude files/diagnostics/calibration_walk/`): `diag_price.py` (the D-bond FOC term by term),
  `diag_fundamental.py` (the bond's value at the model's own discount rates: the perpetuity recursion solved on 81
  points in s with req(s), the Om-weighted s' weights and the default leg read off the solved rules; then the same
  recursion collocated in s alone at 3/5/9/17 nodes), `rest_metrics.py` (the floor sweep's rest-point row for any
  pickle), `probe_solve.py --refine`.
- The s box is centred on p^d 0.1% whatever s_star is (`build_state_box` calls `s_process_params(cal)`, which hard-codes
  0.1%): at 0.02% the rest point sits at u_s = -0.41, between nodes.
- Coarse grid (3 s nodes), 0.02%: the price rule is non-monotone in s (0.977 at the bottom node, 0.996 at u = -0.75,
  0.902 at the centre); at rest it reads 0.980 (yield 45 bp/yr, below the 115 bp deposit rate) against 0.909 at the
  model's own discount rates; sovereign spread at rest -160/-175 bp (exact/fitted). The 1-D recursion at 3 nodes alone
  reproduces the overshoot (0.946); 5 nodes 0.911, 9 nodes 0.9085. The coarse grid misprices the bond at rest at
  every calibration saved: 0.1% 0.847 vs 0.877, 0.05% 0.910 vs 0.885, 0.02% 0.980 vs 0.909.
- Axis-refined grid (5 s nodes, 27 points), 0.02%: joint stage STOPPED at max|F| 2.1e-5 (id_C_D). The rule is monotone
  but uniformly ~0.03 below the recursion's fixed point (0.896 vs 0.930 at rest, yield 260 vs 170 bp); its one-step
  FOC error off the nodes is ~0.002, and the recursion's modulus (1 - delta)/req ~0.94 amplifies that ~17x (the bond's
  duration in quarters). At rest: D wedge 3.5/5.0 bp (exact/fitted), mu_D 0.0007/0.0010 (slack), D net worth +1.8%,
  sovereign spread +37/+40 bp, F wedge 71 bp, leverage 7.85/7.98, exposures 14.31/6.50%; rest point from other
  starts differs by up to 2.07 half-widths. (Coarse grid at 0.02%, for comparison: mu_D 0.0103, binding, with the
  bond at 0.980.)

## 2026-10-23/24: Bond accuracy -- the re-centred cross grid, wired into run.py

- Diagnosis (`diag_price.py` on the 5-node axis solve): on the axes through the box centre the bond FOC holds to
  <8e-4 between nodes; at the rest point's bank state (off-axis: K_D at u 0.42) it misses by 0.7-3.7e-3, all one
  sign. The mu = 1 index set has no cross terms (by construction; `test_cross_refine` shows the axis grid misses a
  s x state product by > 1e-2). The recursion's modulus (1 - delta)/req is ~0.94, so a one-step error e bounds the
  price error by ~17e; measured: one-step ~0.002, price gap ~0.03.
- Experiments (`probe_fix.py`, all at pd_rest 0.0002, s box fixed at p^d 0.002-4.8%):
  X1 axis, 9 s-nodes: it0 (SS box) STOPPED 2.1e-5, bond gap at rest 0.036; it1 (box bank stocks re-centred on it0's
  rest point) CONVERGED 2.3e-10 in 5 steps, gap 0.011; it2 (re-centred again) STOPPED 2.3e-5 and its rest point lay
  off-centre again (K_D u 0.71).
  X2 cross (axis + the level-2 x level-2 pairs of s with each other state, 75 points), 9 s-nodes, from X1 it1 on its
  box: CONVERGED 2.0e-9 in 11 steps (56 min). Bond at rest 0.8907 vs 0.8916 at the model's own rates (2 bp).
- X2 at rest (fitted/exact): mu_D 0.021/0.025 (binding), D wedge 111/128 bp, F 136, leverage 7.89/8.03, exposures
  14.38/6.50%, sovereign spread 55/54 bp, n_D +1.7% vs SS; rest point non-unique (2.06 half-widths across starts).
  2% shock: sovereign spread +592/+594 bp on impact, +159 at q8, +33 at q24; Q_bD -19%; n_D -22%; the D loan spread
  FALLS -108/-126 bp on impact (mu_D 0.0007): fails the acceptance rule "lending spread up". Bond FOC error along
  the path up to 1.6e-3.
- Code: `state_grid` (cal["pd_rest"], S_BOX_PD, the "cross" refinement + `test_cross_refine`),
  `recursive_experiment` (S_REFINE_MODE "cross", CROSS_MAXIT 60, `_refine_ladder`, `_recentre_box`; the TPI reuses
  the re-centred coarse baseline and refines axis 5 -> 9 -> cross 9), calibration pd_rest 0.0002, run.py S_REFINE 9.
- Tests: state_grid, ss_identities, no_unbound_names pass; test_collocation's real solve and nesting N2 fail at
  pd_rest 0.001 too (7.7e-4 and 6.7e-6, so not caused by the pd_rest change); N1, N3-N5 pass.

## 2026-10-24/25: The full run at pd_rest 0.02% (`python3 run.py global`, log `Claude files/diagnostics/calibration_walk/runs/fullrun_pd002.log`)

- Compute 23,114 s (6 h 25 min; ~20.5 h wall clock, the Mac slept); no-TPI results saved at 3 h 38 min. The risk solve
  reproduces X1 + X2 step for step (first 101 Newton steps identical; cross solve 11 steps, 2.03e-9; rest sov 54.75 bp).
- No-TPI IRF (fitted): sovereign spread 647 bp at q0 (~+592 over rest), 342 q4, 214 q8, 88 q24; n_D -21.9%; GDP -0.017%
  fitted / +0.043% exact on impact (mu_D 0.0247 -> 0.0007), trough -0.21% at q24; credit spread -106 bp on impact.
  Acceptance FAILS on impact: Y_D up (exact) and the lending spread down.
- TPI (cap 200 bp): cap walk failed at its first rung (700 bp, 6.8e-6 coarse) -> refined grids solved cap 200 directly
  from the no-TPI seed (axis 5/9 and cross 9 converged; GZ polish 1.1e-8) -> flagged NOT ROOTED, solve_ok False. Rest:
  sovereign spread 55 -> 7.7 bp. Shock: 25/25 quarters at the box wall (M_cb/O_cb and b_DD clipped at q0), GDP -2.63%.

## 2026-10-25/26: The TPI footprint against the box, and the cap walk

Scripts in `Claude files/diagnostics/calibration_walk/`: `tpi_wall.py`, `tpi_reads.py` (reads on the saved run),
`probe_tpi.py` (T1, T2), `probe_tpi_walk.py` (T3, T4), `probe_tpi_refined.py` (T5); logs `runs/t1..t5.log`.
All at pd_rest 0.0002, cap target 200 bp, on the run's re-centred box; "B" is the SS D debt stock.

- Saved TPI solve (75-point cross grid, flagged): the box for the Eurosystem book M_cb is +-25.6% of B on its axis;
  b_DD box 71.6-90.8% of B. Rest: Q_bD 0.8976, purchase 0, b_DD' 80.1% of B. Impact of the 2% shock: Q_bD 0.8457,
  purchase 80.9% of B, b_DD' 0 (the D banks sell their whole book); the clip puts M' at 26.9% and b_DD' at 54.0% of B.
- Impact purchase by shock size, same saved rules (exact reads):

  | p^d shock | no-TPI sov spread (bp) | TPI sov spread (bp) | purchase (% of B) | within the M_cb box |
  |---|---|---|---|---|
  | 0.25% | 167.0 | 160.3 | 0.0 | yes |
  | 0.50% | 252.8 | 200.0 | 8.1 | yes |
  | 0.75% | 328.8 | 200.0 | 22.0 | yes |
  | 1.00% | 399.6 | 200.0 | 29.2 | no |
  | 1.50% | 529.8 | 200.0 | 80.7 | no |
  | 2.00% | 647.4 | 200.0 | 80.9 | no |

- Coarse cap walk, first rung (700 bp), from the no-TPI coarse baseline:

  | test | M_cb band (m_band; box top) | TPI smoothing | no-TPI baseline on that box | rung 700 bp |
  |---|---|---|---|---|
  | T1 | 0.32; 25.6% of B | 1e-4 | 1.98e-10, rest point moves 3.6e-12 | STOPPED 6.2e-6 (40 its) |
  | T2 | 1.1; 88.0% of B | 1e-4 | STOPPED 6.8e-4, rest point moves 12.7% | STOPPED 4.7e-4 |
  | T3 | 0.32 | 1e-2 | not re-solved | STOPPED 4.5e-6 (27 its) |
  | T4 | 1.1 (TPI only; baseline re-seeded) | 1e-2 | not re-solved | STOPPED 6.6e-4 (40 its) |

  T1's twelve worst period-map residuals (3.6e-6 down to 3.5e-7) sit at the centre node and at the +-1 nodes of
  K_D, P_D, P_F, b_FD, V_dep, s and M_cb (both signs), all in the no-default regime.
- T5: the TPI solved directly on the 5-node axis grid (27 points, m_band 0.32), caps in sequence, each seeded from the
  previous (smoothing 1e-2 -> 1e-3 -> 1e-4, then the GZ polish); every cap converged:

  | cap (bp) | final max abs F | rest sov spread (bp) | p^d 1%: purchase, wall quarters | p^d 2%: purchase, wall quarters |
  |---|---|---|---|---|
  | 600 | 9.7e-11 | -74.8 | 0.0%, 2/25 | 14.2%, 2/25 |
  | 400 | 9.6e-10 | -35.4 | 0.0%, 1/25 | 19.1%, 3/25 |
  | 300 | 1.5e-9 | -17.0 | 0.0%, 0/25 | 23.2%, 3/25 |
  | 200 | 2.3e-9 | 5.5 | 13.6%, 0/25 | 29.4%, 3/25 |

  This is the 5-node axis grid, not the production cross grid: at 2% and cap 200 it buys 29.4% of B, the cross grid 80.9%.

Verified: on the current box the production TPI shock path is clipped in every quarter (purchase 80.9% of B against an
M_cb box top of 25.6% at q0, the state clipped at all 25 quarters), so that path is not the model's law of motion.
Widening the M_cb band to 1.1 breaks the coarse no-TPI baseline (T2: 6.8e-4, rest point moved 12.7%).

## Recorded, not fixed (user decision)

welfare_quintiles adds wc_ss to SS dividends, which the SS household block does not.
