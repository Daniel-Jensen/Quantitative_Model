# Global projection guide (`global_projection/`)

Moved verbatim from CLAUDE.md on 2026-10-22 (section text unchanged); CLAUDE.md keeps the rules most tasks need.

## Project background

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

## Model code B — global projection (`global_projection/`)

**Package layout.** The modules are grouped into subpackages; imports are
absolute from the REPO root (`from global_projection.blocks.bank import …`). The
CHEBYSHEV-SMOLYAK PROJECTION pipeline is `compute_global()` / `plot_global()` in the root
`run.py` (the package's own `main.py` was folded into it 2026-09-30) — there is NO
perfect-foresight / representative-branch machinery (`solver_pf/` was deleted 2026-08-11,
and its last leftovers `bank_backward`, `bank_forward`, `govt_transition` on 2026-09-30;
git history preserves them).
- `steady_state.py` — the SS solve; the calibration itself is the shared
  `calibration/global_projection.py` (was `config/calibration.py`)
- `blocks/` — economic blocks (solver-agnostic): `bank.py`, `government.py`,
  `household.py`, `distribution.py`, `rouwenhorst.py`, `fast_kernels.py`,
  `firms.py`, `capital.py`, `trade.py`
- `solver_recursive/` — the ONLY solver: recursive global solution by GLOBAL
  CHEBYSHEV COLLOCATION (Bocola's own design): `state_grid.py`,
  `decision_rules.py`, `point_map.py`, `collocation.py` (the Newton),
  `recursive_main.py` (time iteration, now only a warm start),
  `recursive_experiment.py` (risk + TFP), `tpi_experiment.py` (the TPI backstop)
- `reporting/` — `prints.py` (SS table), `plots.py` (figures)
- `tests/` — regression suite

The model is solved GLOBALLY as recursive decision rules on a Smolyak sparse
grid (Chebyshev interpolation), over the 12-state vector
`[K_D, K_F, P_D, P_F, b_DD, b_DF, b_FD, V_dep, s, Z_D, M_cb, O_cb]` — two capital
stocks, two banks' gross deposit obligations, the three carried sovereign
holdings, the cross-border deposit position, the sovereign-risk factor s, the
TFP state Z_D (deterministic AR(1); the TFP experiment reads the IRF along a
Z-decay path), and the TPI book: the Eurosystem's D-bond holdings and its gross
obligation to the D banks (both zero at the SS, bands symmetric round 0, so the
TPI-off model nests the 10-state one to solver tolerance). At each grid point
FOURTEEN unknowns are solved (the thirteen market-clearing/Euler unknowns plus the
TPI variable `x_cb`; the per-period image of the old stacked system),
with Bocola's closed-form occasionally-binding μ. Expectations are genuine
multi-branch Gauss-Hermite quadrature over the s-innovation × the default regime
d′ ∈ {0,1} — see `decision_rules.regime_table`.

**Driver: GLOBAL COLLOCATION NEWTON (`solver_recursive/collocation.py`), 2026-08-28.**
The policy VALUES at the collocation points are the unknowns and there is no inner
root find — Bocola's `residual_model.m` + `parsolve.m` exactly. Every stored rule is
an unknown (20 per point per regime: the 14 market-clearing/Euler/TPI unknowns plus the
six objects that used to be READ OFF a frozen continuation — alpha, C, r_wc per
country — which now carry Bocola's identity residual `log(guess/implied)`). The whole
coefficient vector goes to one damped Newton with a finite-difference Jacobian
(`parsolve`, dense) or Newton-Krylov (`krylov`, Jacobian-free) on the refined grid.
The FD Jacobian's columns run on every core (`cal["n_jobs"]`: 0 = all, 1 = serial;
`collocation._PoolJacobian`, spawned workers, bit-identical to the serial one). The workers
re-import the model at every solve: **do not edit model code while a run is in progress.**
`cal["jacobian"] = "complex"` (`run.JACOBIAN`, 2026-10-28) replaces the forward differences with the complex
step Im F(θ + ih e_i)/h, h = 1e-20: exact to round-off, same pool, ~1.25× the FD time. The residual path must
stay complex-safe: no `float()` on a traced value, no float-typed arrays for traced values
(`tests/test_complex_jacobian.py` turns ComplexWarning into an error).
**Calm-branch selection (2026-10-28).** A normal-times node can sit on a self-fulfilling crisis branch. Its
continuation reads itself when a slow state such as K_F sits at the box edge (weight 0.97 at K_F+). At appendix B this
stalled every stage: n_D at 11% of n_ss, alpha_D 8, a near-singular Jacobian. `recursive_experiment._calm_branch`,
called after the d0 and joint stages, works as follows:
- it flags nodes with n_D < 0.5 n_ss or alpha_D > 3 in normal times;
- it solves each flagged node alone from the centre node's values, which gives the calm root;
- it polishes all nodes in four Gauss-Seidel sweeps;
- it re-solves the stage.

The selection is the calm root everywhere (author decision, 2026-10-28).
Solve ladder, also his: coarse μ=1 grid → d=0 at π=0 → d=1 by haircut homotopy
(0.85/0.70/0.55/0.45) → joint → SEED the s-refined grid and re-solve there.
Time iteration (`recursive_main.time_iteration`) survives ONLY as the warm start that
puts the Newton inside its basin — it is not a convergent solver here: its binding
mode is the franchise-value recursion at 0.990 per sweep, so runs reported
`max|F| = 1e-14` and "rule-change tol not reached" simultaneously.

**Grid: μ=1 Smolyak × a DENSE Chebyshev factor in s** (`SmolyakGrid(refine=(dim, m))`).
Raising the Smolyak level instead raises the GLOBAL budget; the tensor factor buys
degree m−1 in the one dimension that carries curvature (the logistic p^d(s)) and full
interaction with the sparse basis. Measured relative RMS error on this model's
curvature profile: μ=1 21pts **1.9e-1**, μ=2 221pts **3.9e-2**, m=5 95pts **2.5e-2**,
m=9 171pts **1.1e-3** (measured on the 10-state grid; with the two TPI states the s=5
refined grid is 115 points). `S_REFINE = 5` ships (95 points, ~70 min, before the TPI
states; ~1.5x that now); `S_REFINE = 9` is
Bocola's own resolution and the ladder walks 5 → 9, at ~4 h, because the dense Jacobian
is m+1 = 20·2·n+1 residual evaluations and the solve scales as n².
**CONVERGENCE CHECKED 2026-08-29** at the 100 bp calibration: going 5 → 9 moves the
impact output response from −0.1105% to −0.1087% (fitted) and −0.1278% to −0.1234%
(exact) — 1.6% and 3.4%, both well inside the 13%-wide identification bracket — and
every other reported number in the third digit (credit spread +90.6 → +90.8 bp/yr,
Q_bD −9.263 → −9.304%, Euler ALL −4.37 → −4.38). Re-run 2026-10-01 at 5 through `run.py` after the bit-for-bit refactor: −0.1105% / −0.1278%, +90.6 bp/yr and −9.263% reproduce exactly; Euler ALL reads −4.34, which the refactor cannot have moved (old and new code print identical accuracy tables), so the drift predates it. **The solution is converged at 5 for
every reported object**; 9 is the confirmation, not the working setting. NB this was
NOT true before the recalibration, when refining was fighting a Gibbs phenomenon at
the KKT kink. Hot kernels (household EGM backward, distribution
forward) are numba-JITed with an exact pure-numpy fallback (`cal["use_numba"]`).

| File | Contents |
|------|----------|
| `calibration/global_projection.py` (shared dir) | All parameters. Single λ per bank (Bocola IC); Bocola/Greece anchors documented inline. Credit spread 100 bp/yr (NOT his 8 — see the kink note above); leverage 5, exposure 7.6%, recovery 0.45 are his. f = exit/payout share; Ω = β·[f + (1−f)α′] (Bocola's ψ = 1−f survival weight on the franchise value). |
| `global_projection/steady_state.py` | Two-stage SS solve: {rk_D, rk_F, p} on capital markets + current account, then {β_D, β_F} on deposit markets. Symmetric SS required (see docstring). |
| `bank.py` | GK/Bocola bank block at the steady state: `calibrate_bank_targets` (single λ and ω_ent from leverage + spread), `steady_state_bank` (bond prices, net worth, portfolio, dividends). The dynamic bank problem lives in `solver_recursive/point_map.py`; only D is risky, F bonds are safe. |
| `government.py` | HM perpetuity bonds, Bohn rule: `govt_steady_state` (Bohn coefficient solved from a target debt root). The per-period fiscal block is in `point_map._government`. Default risk is exogenous (no crisis zones). |
| `solver_recursive/point_map.py` | The per-point period map (image of the old stacked system): 13 residuals at one grid point given the frozen continuation rules. Bocola closed-form μ; quadrature over the s-innovation × the default regime; Z_D read from the state. |
| `solver_recursive/state_grid.py` | Smolyak sparse grid + Chebyshev basis, with `refine=(dim, m)` for a dense tensor factor on one dimension; `build_state_box`, `default_prob`, `s_process_params`. |
| `solver_recursive/collocation.py` | THE SOLVER. `make_residual` (the global F(theta), image of `residual_model.m`), `parsolve` (port of his damped FD Newton), `krylov_solve`, `solve_collocation`. |
| `solver_recursive/recursive_main.py`, `recursive_experiment.py`, `tpi_experiment.py` | Time iteration (warm start only) + SS anchors; the risk + TFP experiments and the solve ladder (with the TPI cap homotopy); the TPI against no TPI (rest point, IRF with footprint, bond-price legs, price impact, KT check). |
| `fast_kernels.py` | numba kernels for EGM backward + distribution forward; exact numpy fallback when numba is absent (`cal["use_numba"]`). |
| `household.py`, `distribution.py` | EGM with GHH utility; stationary distribution and forward iteration. |
| `trade.py` | CES basket and bilateral flows with PER-COUNTRY home bias and the country-mass ratio (`size_ratio`); `omega_home_F` is derived from `omega_home_D` and the sizes so trade balances at p = 1. |
| `firms.py`, `capital.py` | Flexible-price production with the Neumeyer-Perri working-capital wedge (w ÷ (1+ζ·r_wc), the spread→output channel; ζ=0 nests exactly — Bocola §V.C's own open-economy fix), Jermann adjustment costs, CES/Armington trade. |
| `prints.py` | Console reporting: `banner`, `print_ss_table`, and THE UNIT CONVENTION (`bp_ann`, `ann_pct`, `ann_prob`, and the `BOCOLA_IRF_*` benchmarks). Rates are annualised bp; p^d is printed quarterly AND annual; flow responses in level % with a ×4 annualised companion — Bocola's Table 5 unit. |
| `plots.py` | Every global figure; written to `plots.OUTDIR` (default `results/GLOBAL/figures/`, set by `run.py`). |
| `export.py` | Writes `results/GLOBAL/data/comparison_irfs.json` for `compare.py` (the SSJ twin is `linear_ssj/reporting/export.py`). |
| `run.py` (repo root) | THE PIPELINE: `compute_global` (SS → TFP → risk pass-through → decompositions → accuracy → TPI → save) and `plot_global` (figures from the saved data). `RUN_TPI = True` adds one solve (the coarse cap homotopy, then the refined grid); `TPI_CAP` sets the cap (bp/yr, or `"rest"` = defend the no-TPI rest-point spread, ~56 bp). The 12-state grid costs ~1.5x the old 10-state one. |
| `tests/` | Regression suite (see below). |

## Running and testing

### B. Global projection pipeline (`global_projection/`)

```bash
# from the repo root
python3 run.py global                                 # full projection pipeline (SS+TFP+risk+TPI)
python3 run.py global --quick                         # coarse grid (preview, NOT converged)
python3 -m global_projection.solver_recursive.recursive_experiment   # risk pass-through only
python3 "Claude files/global_projection/tests/test_ss_identities.py"      # SS theory identities (fast)
python3 "Claude files/global_projection/tests/test_fast_kernels.py"       # numba/numpy kernel equivalence (fast)
python3 "Claude files/global_projection/tests/test_state_grid.py"         # Smolyak grid exactness (fast)
python3 "Claude files/global_projection/tests/test_no_unbound_names.py"   # static unbound-name gate (fast)
python3 "Claude files/global_projection/tests/test_collocation.py"        # THE SOLVER: packing, the six identity
                                                           #   residuals, the refined grid, and a real
                                                           #   d=0 solve to max|F| ~ 1e-9 (~90 s)
python3 "Claude files/global_projection/tests/test_recursive_nesting.py"  # SS rest point (N1), the pi=0 grid-wide
                                                           #   solve (N2), the TPI swap + risk transfer
                                                           #   (N3) and the union budget identity (N4)
```

**Comment convention** (enforced across `global_projection/`): every module and every
function carries exactly ONE leading ALL-CAPS comment saying what it is; any
further explanation is lowercase `#` comments attached to the specific hard
line. No docstrings, no bold markers, no prose blocks inside function bodies.
Console output lives in `prints.py`, never inside the model blocks.

**Acceptance thresholds** (all enforced in tests):
- Global collocation: Bocola's own test, `sum(F^2) <= m*(1e-9)^2` — the sum a
  uniform `max|F| = 1e-9` (`collocation.TOL_MAXF`) would give — at EVERY stage, over
  the 20 equations × points × regimes. This replaces the old two-part
  time-iteration test (settled rule AND every point clearing), which could pass on
  residuals while the rules were still moving. 1e-9 rather than machine zero because
  the period map's arithmetic floor is ~1e-10: the capital and bond Eulers difference
  O(1) expectations down to O(1e-4), and no Newton step improves on that. It is still
  four orders below any economic signal (the headline shock moves μ by 7.5e-3).
- goods_D (imposed) ≤ 1e−9; goods_F (Walras-redundant diagnostic) ≤ 2e−6 —
  including when the debt stock moves.  (The Newton solver typically lands
  goods_D near 1e−13; acceptance is `tol_transition` = 1e−10 normalized.)
- Zero-shock transition stays at SS to ≤ 1e−5.
- Risk-only shock (exogenous π): Q_bD↓, n_D↓, n_F↓, Y_D[0]↓, C_D[0]↓,
  lending spread↑, b_gov↑, Tax↑ (a positive Y or n response to sovereign
  risk = bug).
- Complementarity on every solved path: μ ≥ 0, slack = αn − λ·assets ≥ 0,
  μ·slack ≈ 0 (`out["mu_D/F"]`, `out["slack_D/F"]` from point_map.py).
  Known open item: risk-on n_D[0] can sit above
  risk-off (M1 deposit-rate channel; test warning, not assert) and
  post-impact Y_D runs mildly positive — both die with the union deposit
  market (§8 of the retired sunspot_transition_study.md; git history).

## Key modelling choices — do not "fix" without checking Claude files/docs/SPEC.md

- **Single λ (Bocola 2016 eq. 3):** all three asset classes carry the same
  divertability. Diverging them re-opens the portfolio-substitution margin
  that made sovereign risk *expansionary* pre-rework.
- **Priced vs realized default:** `def_price` enters bond pricing and
  expected-return FOCs; `def_real` enters realized returns and government
  flows. The baseline experiment prices risk but never realizes it
  (Bocola's pass-through design); a realized-default variant just passes
  `def_real ≠ 0`.
- **Endogenous debt in clearing:** the government's end-of-period stock is
  forward-integrated inside every residual evaluation and absorbed by banks.
  Clearing against a fixed `B_gov_ss` instead re-opens a Walras leak of
  ~0.5% of GDP per 5% debt deviation.
- **ASYMMETRIC COUNTRY SIZE, SYMMETRIC PER-CAPITA STEADY STATE (2026-08-28):**
  `size_F/size_D = 8`. Every variable is PER CAPITA of its own country and the
  per-capita SS is UNCHANGED (p_ss = 1, identical n_ss, leverage 5, μ_ss = 0.001,
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
- **Occasionally-binding IC (Bocola):** the leverage constraint enters the
  stacked system as the Fischer-Burmeister complementarity between μ (from
  the capital FOC, valid in both regimes) and slack = αn − λ·assets, scaled
  by μ_ss and n_ss (FB's zero set is scaling-invariant). At the SS the
  constraint binds (μ_ss ≈ 0.02, slack = 0), where FB is smooth. Portfolio
  shares and branch initial conditions divide by ACTUAL net worth, not n_IC.
- **Ω-kernel weights (Bocola Prop. 1):** f = exit/payout share, so
  Ω = β·[f + (1−f)·α′] — weight 1−f ≈ 0.95 on the franchise value α′
  (Bocola's survival ψ). beta_inter ≈ β_hh ≈ 0.99 proxies the household SDF;
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
  and the Chebyshev-Smolyak projection solver is now the ONLY machinery (TFP is a
  deterministic 7th state Z_D).
- **TPI = A SPREAD-CAP BOND-PURCHASE BACKSTOP (2026-10-02, branch `OMT-fix`; replaces the
  LTRO).** Design, derivations and results: `Claude files/docs/tpi_backstop_plan.md`. In the
  no-default regime the Eurosystem buys the D bond whenever `y_D - y_F` (HM flow yields)
  would exceed `cal["tpi_cap_bp"]`, i.e. `Q_bD >= delta_D/(delta_D + y_F + cap)`, and pays
  with a SAFE claim `Z = Q_bD*M` on itself held by the D banks (paying `1 + rdep_D`). Held to
  maturity, nothing bought in default, P&L `Pi = Xi*M_lag - O_lag` remitted by capital key
  (`tpi_key_D = 0.071`: D issuance; the rest lowers F taxes). `tpi_on = False` (default)
  nests the no-TPI model.
  **A PURCHASE IS A SWAP TODAY AND A RISK TRANSFER TOMORROW.** Assets, the divertable base
  (single lambda), deposits, `P'` and `mu` do not move within the period; next period the
  bank is owed `R*Q*dm` instead of `Xi'*dm` (N3 asserts both to 1e-12). This is the yield
  peg's failure fixed on purpose: the peg could only remove the liquidity premium.
  **STATES 12, UNKNOWNS 14.** `M_cb` (the book) and `O_cb` (the obligation; it carries last
  period's price), both zero at the SS with bands symmetric round 0, so the M = O = 0 slice
  IS the old 10-state grid and the TPI-off model nests it to solver tolerance (rules
  5.5e-10). The box is a SHEAR, `z_bDD = b_DD + M`, `z_O = O - rho*M` (`state_grid`): on the
  natural axes a forced purchase moved the price non-monotonically. The 14th unknown `x_cb`
  is Garcia-Zangwill: `m = B*max(x,0)`, `Q_bD = floor + Q_ss*max(-x,0)` in the no-default
  regime, so the stored rule is smooth across the floor (a stored `m` bought 19.7% of the
  stock where the floor does not bind, by interpolation). Residual 8 is the D-bank FOC as a
  Fischer-Burmeister KT pair with `b_DD >= 0` (the corner where the Eurosystem holds the
  bank's whole book); residual 14 ties the stored `Q_bD` rule to that price. The
  continuation reads `Q_bD'` off the stored rule: `max(-x',0)` inside the expectations
  stalls the finite-difference Newton.
  **ROOT IN ONE FORM, READ IN THE OTHER (`cal["tpi_gz"]`).** The GZ Newton also stalls
  walking the cap down, so the solver roots the system in the Fischer-Burmeister form
  (`tpi_gz = False`: x is the purchase share, residual 14 is `FB(x, (Q - floor)/Q_ss)`) and
  hands back in the GZ form (`_tpi_form` converts at the nodes, exactly -- both forms
  describe the same allocation there -- then `_tpi_polish` re-roots in 1-2 steps).
  **SOLVE:** the no-TPI baseline, then a cap homotopy (`recursive_experiment._solve_tpi`,
  700 bp -> target in 50 bp rungs, a failed rung retried at half the step), then the
  s-refined grid 3 -> 5 nodes (FB, then GZ), each refined solve walking the smoothing in
  (`TPI_EPS_LADDER` 1e-2 -> 1e-3 -> `tpi_eps` = 1e-4): the 5-node grid puts nodes on the
  floor's boundary (it binds at 43 of 115) and at the final smoothing from the start the
  Newton stalled.
  **RESULTS (full grid, cap 200 bp, `Claude files/docs/tpi_backstop_plan.md` §9.2):** the
  announcement lowers the rest-point lending spread 89 -> 53 bp, output +0.18%, nothing bought;
  on the headline shock the Eurosystem buys 43% of the stock on impact, impact output -0.080%
  (no TPI -0.128%), bank net worth -1.5% (-4.2%), but the downturn lasts longer (trough -0.147%
  at q6): the IRF conditions on no default, and there the bonds' excess return goes to the
  Eurosystem and by the capital key 92.9% on to F.
  **PRICE IMPACT PER EURO IS SMALL:** at fixed rules 30% of the stock lifts `Q_bD` 0.61% at
  the headline shock (0.24% at rest), against +3.8% for the solved TPI economy, most of it
  through the continuation; a cap far below the market spread is met only at the corner (21 of
  the 23 nodes at p^d 4.8%/qtr), while at the headline shock the purchase is interior (the D
  bank keeps 39% of the stock).
  **WALRAS:** `goods_F` is reported; it exposed a PRE-EXISTING leak (the constant household
  anchor `hh_T` stands in for the working-capital flow), see STATE Part II. N4 checks that
  the TPI's own flows cancel out of the union budget exactly (2e-10 at random points).
  **PREDECESSORS, RETIRED:** the LTRO (deleted 2026-10-01, `ltro_backstop_plan.md`) and the
  one-sided yield peg; `liquidity_ceiling_report` is the peg's diagnostic, kept.
- **Predetermined deposit rate:** the rate paid at t was locked at t−1
  throughout (bank funding legs, household EGM returns, μ timing).
- **Predetermined capital (Bocola eq. 6):** the stock producing at t was
  bought at t−1 (`Kap_prod[t] = Kap[t−1]`); mpk is the marginal product of
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
- **Hatchondo-Martinez perpetuity:** stock decays at rate 1−δ_b; duration
  ≈ 1/δ_b quarters (0.036 ⇒ ~7y). Long duration is what makes priced risk
  generate large MTM losses — an interlude with δ_b=0.25/recovery=0.80 cut
  the repricing ~6x and made the risk channel expansionary (study §8).
- **Default branch = ONE pure-haircut feared event (Bocola):** the branch
  solves a single deterministic event — a full write-down to recovery
  `recovery_rate_D` = 0.45 (Greek PSI; Bocola's D = 0.55) on the whole
  claim, with the default-state recession arising endogenously through bank
  balance sheets. There are no scarring add-ons: the old Arellano output
  cost, GK ξ_K capital-quality loss and HFSF recap FLAGS were all 0 at the
  Bocola-pure baseline and were deleted in the 2026-07-21 cleanup (see git
  history if a variant needs them back). The recap machinery survives only
  as `_RECAP_LADDER`, a warm-start continuation used when the direct branch
  solve stalls. If the event is infeasible after that, the branch RAISES.
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

### B. Global projection — known limitations

## Known limitations (documented, next thesis phases)

- Comovement problem RESOLVED (2026-07-18): predetermined capital + the
  union deposit market restored the impact contraction at the headline
  shock (Y_D[0] and I_D[0] both negative at π = 1%·0.95^t). Remaining
  next-phase dial: NK/union nominal block (needed for the TPI
  application); real interest parity currently plays the role of the
  single policy rate.
- Risk channel (recursive) approximations: Λ^nd ≡ beta_inter, rep-agent
  income-SDF proxy for Λ^d, household-side π-blindness (the deposit Euler never
  weights the default branch — faithful to Bocola, where household deposits are
  riskless too).
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
