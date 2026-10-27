# Global projection: modelling choices and known limitations

Relocated from CLAUDE.md on 2026-10-26. Every entry is a deliberate choice with its reason; do not
"fix" one without reading it and `Claude files/docs/SPEC.md`. Text marked `[corrected 2026-10-26]`
or `[STATUS 2026-10-26]` was checked against the live code and the viva paper during the move;
everything else is verbatim and keeps its original date. Numbers quoted here are historical unless
they name a calibration that still matches `calibration/global_projection.py`.

## Background (the project paragraph formerly at the top of CLAUDE.md)

Two-country heterogeneous-agent model of a monetary union with Gertler-Karadi
financial intermediaries and sovereign default risk, calibrated to the
2010–2012 Greek sovereign debt crisis. The default mechanism follows
**Bocola (2016, JPE) "The Pass-Through of Sovereign Risk"**: an EXOGENOUS
rise in the *priced* probability of default π_t (his s-shock, eqs. 11–12 —
an input path, never a function of debt) makes bond prices fall, banks take
mark-to-market losses, the single-λ occasionally-binding incentive
constraint tightens, lending spreads rise and output falls — with no default
ever realized. Only D is default-risky; F bonds are safe. Application: ECB
asset purchases (TPI). Primary output is a research paper (Overleaf:
https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c).
(The 2026-07-16 Bocola-faithful rewrite replaced the earlier Cole-Kehoe
crisis-zone wrapper, the always-binding IC, and the patched default branch
*within the global-projection pipeline*; see git history on branch
`bocola-rewrite`.)

**The repository carries TWO solvers of this model.** The sequence-space (SSJ)
pipeline in `linear_ssj/` is the linearised sticky-price / nominal-deposit model
and is what `experiments/` (E1–E4), `diagnostics/regimes/` and the paper's
current figures run on. The global-projection pipeline in `global_projection/` is the
nonlinear occasionally-binding solution. They share no model code and no interpreter;
neither supersedes the other. Both calibrations live side by side in the shared
`calibration/` package (`calibration/ssj.py`, `calibration/global_projection.py`),
`run.py` is the single entry point, and `compare.py` sets the two IRFs side by side
(layout and rationale: `Claude files/REFACTOR_ARCHITECTURE.md`, 2026-09-30). Say which one you mean before quoting a number:
the two are calibrated differently and their impulse magnitudes are not
comparable.

## Key modelling choices — do not "fix" without checking Claude files/docs/SPEC.md

- **Single λ (Bocola 2016 eq. 3):** all three asset classes carry the same
  divertability. Diverging them re-opens the portfolio-substitution margin
  that made sovereign risk *expansionary* pre-rework.
- **Priced vs realised default:** risk is priced through the expectation over next period's
  default regime d′ ∈ {0,1}, weighted by the exogenous π_t, and the reported impulse response
  follows the no-default history, so risk is priced but never realised (Bocola's pass-through
  design). The default regime has its own solved rules. [corrected 2026-10-26]: the former
  `def_price` / `def_real` inputs belonged to the deleted perfect-foresight solver.
- **Endogenous debt in clearing:** the government's end-of-period stock is
  computed from its budget and the Bohn rule inside every period-map evaluation and
  absorbed by banks (and, under the TPI, the Eurosystem) [corrected 2026-10-26]: was "forward-integrated inside every residual evaluation", the transition
  solver's wording; the leak figure below was measured then.
  Clearing against a fixed `B_gov_ss` instead re-opens a Walras leak of
  ~0.5% of GDP per 5% debt deviation.
- **ASYMMETRIC COUNTRY SIZE, SYMMETRIC PER-CAPITA STEADY STATE (2026-08-28):**
  `size_F/size_D = 8`. Every variable is PER CAPITA of its own country and the
  per-capita SS is UNCHANGED (p_ss = 1, identical n_ss, leverage 5, μ_ss = 0.001 at the time [STATUS 2026-10-26: ≈0.012 at the code's 100 bp wedge; 0.001 dated from the 8 bp calibration],
  identical deposit supply); the mass ratio enters ONLY where D and F quantities
  are aggregated — goods market, union deposit clearing, both sovereign markets,
  the union wealth identity `W_F = P_F − V/(sz·p)`. Sovereign holdings are carried
  in the ISSUER's per-capita units, so `b_DD + b_DF = B_D` still clears the D
  market and the F bank's own book holds `b_DF/sz`. Home bias MUST scale with size
  or trade cannot balance: `(1−ω_F) = (1−ω_D)·size_D/size_F`, so D imports 15% of
  its basket and F imports 1.875% of its (`omega_home_F` is DERIVED in
  calibration.py). WHY: with a symmetric union D is half the union, so D's own
  sovereign shock moved the union real deposit rate 45 bp/yr and cancelled 78% of
  the credit-spread rise before it reached any firm's wage bill — the 2026-08-28
  audit's finding. Bocola's §V.C open economy has no such feedback: his
  `R = 1/β + 0.01·(B_for/gdp)` is a WORLD rate. `size_F = size_D` nests the old
  symmetric model exactly.
- **Symmetric steady state (in per-capita ratios):** other country asymmetries
  enter through shocks only.
  An asymmetric SS (e.g. δ_b_D ≠ δ_b_F) shifts p_ss off 1 and opens an
  O(1e−4) SS goods-market wedge (p is weakly identified by external balance
  at trade elasticity 0.5; see steady_state.py docstring).
  [STATUS 2026-10-26, re-verify: the viva paper's calibration is asymmetric per capita
  (delta_b 0.056/0.039, G 0.204/0.191, own-sovereign exposure 14.4%/6.5%, size 11.3), so it
  cannot be solved under this doctrine without extending the SS solve; the derivation of
  omega_home_F also assumes trade balances at p = 1.]
- **Occasionally-binding IC (Bocola):** μ is computed in closed form at every point,
  μ = max{1 − E[Ω]·(1+r)·n/(λ·assets), 0} (smoothed), so the constraint binds or is slack state
  by state, with slack = αn − λ·assets. Fischer-Burmeister pairs are used only for the D bank's
  b_DD ≥ 0 corner and the TPI price floor. Portfolio shares divide by ACTUAL net worth.
  [corrected 2026-10-26]: the former text described the deleted stacked system (an FB pair
  between μ and slack scaled by μ_ss ≈ 0.02 and n_ss, and "branch initial conditions").
- **Ω-kernel weights (Bocola Prop. 1):** f = exit/payout share, so
  Ω = β·[f + (1−f)·α′] — weight 1−f ≈ 0.95 on the franchise value α′
  (Bocola's survival ψ). beta_inter ≈ β_hh ≈ 0.99 proxies the household SDF;
  [STATUS 2026-10-26: the code has f = 0.08, so 1−f = 0.92, and beta_inter = 0.997 = 1/(1+r̄)]
  values ≪ 1/(1+rdep) drive α_ss below 1 and mute the franchise channel
  (the pre-rewrite code had the weights swapped AND beta_inter = 0.96).
- **Risk channel = genuine multi-branch quadrature (solver_recursive/), NOT a
  representative branch.** The default fork enters `point_map.py`'s banker FOCs
  as a real probability-weighted integral: Gauss-Hermite over the s-innovation ×
  the default realization d′∈{0,1} weighted by π_t (EXOGENOUS input path), where
  the default state is the SAME fitted decision rules evaluated at a reachable
  next-period point — never a frozen stand-in economy. The premium is endogenous
  (Ω^d > Ω^nd on the low default payoffs). `pi ≡ 0` nests the risk-neutral model
  exactly (test_recursive_nesting). The earlier perfect-foresight
  representative-branch pricing got the sign wrong (expansionary); the entire PF
  stack (`solver_pf/`: transition + solvers + risk_branch) was deleted 2026-08-11
  and the Chebyshev-Smolyak projection solver is now the ONLY machinery (TFP is the
  deterministic state Z_D, 10th of the 12 states [corrected 2026-10-26]: was "7th").
- **TPI = a spread-cap bond-purchase backstop (2026-10-02).** Design, derivations and solution
  method: `Claude files/docs/tpi_backstop_plan.md` (the CLAUDE.md summary was moved there
  2026-10-26). Results of the 6 Oct 2026 run: `STATE.md` Part II.
- **Predetermined deposit rate:** the rate paid at t was locked at t−1
  throughout (bank funding legs, household EGM returns, μ timing).
- **Predetermined capital (Bocola eq. 6):** the stock producing at t was
  bought at t−1 (`Kap_prod[t] = Kap[t−1]`; that name survives in
  `blocks/firms.solve_firm_path`, while the collocation solver carries the stock as a state); mpk is the marginal product of
  the bank-held vintage, so impact output moves through hours alone. The
  old contemporaneous timing let the sovereign-risk investment boom raise
  Y_0 directly — reverting it re-opens the comovement problem.
- **Union deposit market (deposit-UIP):** deposits are own-good claims at
  national rates; a frictionless union interbank replaces the two national
  clearings with ONE union-wide clearing (D-good units) plus real-rate
  parity (1+rdep_D) = (1+rdep_F)·p′/p — the flexible-price image of one
  nominal union rate + national inflation differentials (BKK/Baxter-Crucini
  single-traded-bond margin). UIP makes the interbank pass-through
  zero-profit → no Walras leak; the cross-border deposit position
  (`out["nfa_dep_D"]`) is the absorption margin that broke the national
  S=I trap (the M1 comovement mechanism). A literal rdep_D=rdep_F with
  own-good legs is WRONG (unassigned RER valuation profit → Walras leak).
  Stage-2 SS imposes β_F = β_D (symmetric-SS doctrine).
  **KEPT AS THE BASELINE 2026-10-07 (user decision).** The union-real-rate shortcut (drop
  E[p′]/p, keep the premium) is that literal equality and is shelved (`git stash`
  "SHELVED 2026-10-07"); so is the nominal block. Parity's D-crisis cushion (D's deposit rate
  falls ~48 bp at the headline shock, absorbing ~60% of the lending-spread rise) is a
  flexible-price limitation to STATE, not fix.
  [STATUS 2026-10-26, re-verify: still the code's baseline, but the viva paper (section 2,
  eq. dep-national; equilibrium condition (v)) switched to national deposit clearing on
  2026-10-14. The code has not followed, so no solved number describes the paper's deposit
  market.]
- **Hatchondo-Martinez perpetuity:** stock decays at rate 1−δ_b; duration
  ≈ 1/δ_b quarters (0.036 ⇒ ~7y; [STATUS 2026-10-26: the code has 0.056, about 4.5 years,
  and the paper 0.056/0.039]). Long duration is what makes priced risk
  generate large MTM losses — an interlude with δ_b=0.25/recovery=0.80 cut
  the repricing ~6x and made the risk channel expansionary (study §8). Duration works
  through its interaction with the persistence of s, not alone (see the SSJ entry in
  `ssj_pipeline.md`, 2026-08-28).
- **Default branch = ONE pure-haircut feared event (Bocola):** the branch
  solves a single deterministic event — a full write-down to recovery
  `recovery_rate_D` = 0.45 (Greek PSI; Bocola's D = 0.55) on the whole
  claim, with the default-state recession arising endogenously through bank
  balance sheets. There are no scarring add-ons: the old Arellano output
  cost, GK ξ_K capital-quality loss and HFSF recap FLAGS were all 0 at the
  Bocola-pure baseline and were deleted in the 2026-07-21 cleanup (see git
  history if a variant needs them back). [corrected 2026-10-26]: the former
  `_RECAP_LADDER` warm start and the raise-on-infeasibility no longer exist. In the collocation
  solver the default regime's rules are reached by the haircut homotopy in the solve ladder
  (recovery 0.85, 0.70, 0.55, then 0.45). The 2026-10-14 theoretical re-audit
  (`Claude files/diagnostics/model_audit/THEORETICAL_REAUDIT.md`) raises the default regime's
  domain as an open issue.
  Bohn taxes respond to the SURVIVING stock (taxing the pre-haircut stock at
  t=0 was a ~31%-of-GDP artifact).
- **Working capital (Neumeyer-Perri):** ζ_wc=1 × wage bill pre-financed at
  r_wc = rdep(−1) + λμ/Ω̃; the wedge is the only channel from spreads into
  impact output (without it Y_D moved −0.2% even at Q_bD −30%, n_D −20%).
  The LOAN is a bank asset inside the divertable base at the same λ, and the
  financing income accrues to BANK net worth through the deposit obligation
  `P' = R(QK'+qB'+L−N') − R_W·L` — Bocola's `residual_model_open.m` exactly.
  (`cal["wc_rebate"] = 1.0` instead hands it to households as a dividend, which
  makes the spread a pure intra-period transfer and, with no GHH wealth effect to
  offset it, turns the risk channel expansionary. Default is 0.)
- **Walras redundancy:** goods_F and the current account are *dropped* from
  the residual system and monitored as diagnostics.
- **Policy rules present:** the Bohn tax and the TPI backstop above (off unless
  `cal["tpi_on"]`). No macroprudential policy, by design.

## Known limitations (documented, next thesis phases)

- (Historical, perfect-foresight era) Comovement problem RESOLVED (2026-07-18): predetermined capital + the
  union deposit market restored the impact contraction at the headline
  shock (Y_D[0] and I_D[0] both negative at π = 1%·0.95^t). Remaining
  next-phase dial: NK/union nominal block, SHELVED 2026-10-07; real
  interest parity plays the role of the single policy rate.
- Risk channel (recursive) approximations: the banker kernel Λ is the
  rep-agent household GHH SDF in each regime. Households' deposit Euler DOES
  weight the default branch by π (`point_map._households`, `v.wgt`); the former
  note that it never does was stale (corrected 2026-10-07).
- **THE DETERMINISTIC SS IS NOT THE STOCHASTIC REST POINT** (2026-08-29). Verified:
  solved at pi == 0 the model sits on the deterministic SS for 200 quarters to six
  decimals with mu = 0.001001 and max|F| ~ 1e-7, so the SS and the solver are exact.
  But the risk-pricing rules rest elsewhere -- measured Y_D -0.13%, C_D -0.13%,
  I_D +0.23%, n_D +2.1%, K_D -0.58%, b_DD +2.6%, and **mu_D falls to EXACTLY 0**.
  Grid-independent (coarse vs s-refined agree to ~10% of the gap) and a UNIQUE GLOBAL
  ATTRACTOR (eight perturbed starts converge to the same state to 1e-12). Bocola has
  the same gap (his ergodic q = 0.979 against a deterministic 1.000, debt +2.0%) and
  locates it the same way: `recursive_experiment.stochastic_rest_point` is his
  `generate_irf.m` step 1, a zero-shock simulation to convergence. EVERY IRF starts
  there AND is differenced against an unshocked path (his
  `gdp = mean(gdp_s) - mean(gdp_nos)`). Before that fix the IRF charged the walk
  between the two rest points to the shock: 54-63% of the reported post-impact hump
  and of the +4.2% bank-net-worth overshoot was drift, and the "capital grinding
  down" was almost entirely drift.
- **THE ECONOMY RESTS ON THE KKT KINK, AND THAT IS THE BINDING ACCURACY LIMIT.**
  mu = max{.,0} is C0 and mu = 0 exactly at the rest point, so a Chebyshev interpolant
  returns mu > 0 in a neighbourhood where the truth is 0. Reading the fitted rules and
  clearing the period map exactly at the same state therefore disagree by MORE than the
  response: Y_D at p^d = 1.98% is -0.081% fitted against -0.008% cleared. The gap
  shrinks with resolution (0.087 -> 0.073 pp from 21 to 95 points) but slowly, as a
  Gibbs phenomenon does. THE SAME PATHOLOGY IS IN BOCOLA'S OWN SOLUTION: his fitted mu
  policy returns a 28.4 bp liquidity premium on impact where the exact multiplier gives
  2.1 bp, and 28.4 is his published number. `impact_table` and `dynamic_irf` print BOTH
  reads (`read_exact`); the pair is the honest object and the level of the output
  response is NOT identified at the current resolution.
  **CURED 2026-08-29** by moving `credit_spread_target` 8 -> **100 bp/yr**, which is
  where the constraint starts binding at the rest point (mu_rest 0 -> 0.0098) and the
  identification gap collapses 4x (0.087 -> 0.021 pp). f does NOT do this:
  `calibrate_bank_targets` forces alpha_ss = lambda*theta at the SS whatever f is, so
  the binding margin is ~theta*mu_ss and mu_ss is set by the SPREAD (measured: 8 bp ->
  mu_rest 0; 25 -> 0; 100 -> 0.0098; 250 -> 0.0286 with no further identification
  gain). 100 bp is Gertler-Kiyotaki's own target and inside the 100-300 bp periphery
  lending spreads of 2011-12; Bocola's 8 bp is his ESTIMATE, but it belongs to his
  closed model where mu contributes 2 bp to output and the wedge channel does not
  exist -- his SS V.C transmission and his closed-model mu^bg cannot both be imported.
  [STATUS 2026-10-26, OPEN: the cure holds for the code calibration only. The viva paper sets
  the wedge to 5 bp/yr at leverage 8, which gives a deterministic mu_bar = 0.0010, the order
  of the 8 bp case whose rest point sat at mu = 0. The numerical identification problem may
  therefore return under the paper calibration; whether it does is unknown until that
  calibration is solved and the fitted and exactly cleared responses are compared.]
- THE OUTPUT CHANNEL IS THE WORKING-CAPITAL WEDGE ALONE. GHH removes Bocola's
  closed-economy channel (the labour-supply wealth effect: in his benchmark
  `dlog l = −1.25·dlog c` exactly, and the leverage multiplier contributes 2 bp),
  so output moves only through `r_wc = rdep + λμ/E[Ω]`. BOTH legs matter, which
  is why country size is now asymmetric — see the key-choices list.
- THE BENCHMARK. Bocola's Table 5 (−1.05/−1.44/−1.53) is a cumulated quarterly
  GROWTH gap ×400 over an 8-quarter estimated shock sequence — its output LEVEL
  equivalent is −0.26/−0.36/−0.38%. The like-for-like single-shock IRF targets,
  rescaled to p^d = 1.98%/qtr, are **−0.157% (his §V.C open economy, whose GHH +
  working-capital transmission this model shares)** and −0.222% (his closed
  benchmark). `reporting/prints.py` carries these as constants and `dynamic_irf`
  prints them next to the trough.

