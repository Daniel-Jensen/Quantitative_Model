# Global projection solver: reference

Relocated from CLAUDE.md on 2026-10-26. Load this file when working on `global_projection/`.
It describes the solver; dated results live in `STATE.md` Part II with their calibration and
date. Text marked `[corrected 2026-10-26]` was checked against the live code during the move;
everything else is carried over verbatim and keeps its original date.

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
  `recursive_experiment.py` (risk + TFP), `tpi_experiment.py` (the TPI backstop).
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
Convergence in the s-dimension was checked on 2026-08-29 and reproduced on 2026-10-01, both under the code calibration of that time (leverage 5, wedge 100 bp/yr); the numbers are in `STATE.md` Part II. Any new calibration needs its own check, and a rest point on the KKT kink would void the conclusion that 5 nodes suffice.
Hot kernels (household EGM backward, distribution
forward) are numba-JITed with an exact pure-numpy fallback (`cal["use_numba"]`).

| File | Contents |
|------|----------|
| `calibration/global_projection.py` (shared dir) | All parameters. Single λ per bank (Bocola IC); Bocola/Greece anchors documented inline. Credit spread 100 bp/yr (NOT his 8 — see the kink note above); leverage 5, exposure 7.6%, recovery 0.45 are his. [STATUS 2026-10-26, re-verify: these are the code's values; the viva paper now uses leverage 8, a 5 bp/yr wedge and own-sovereign exposures 14.4% (D) and 6.5% (F).] f = exit/payout share; Ω = β·[f + (1−f)α′] (Bocola's ψ = 1−f survival weight on the franchise value). |
| `global_projection/steady_state.py` | Two-stage SS solve: {rk_D, rk_F, p} on capital markets + current account, then {β_D, β_F} on deposit markets. Symmetric SS required (see docstring). [STATUS 2026-10-26: the paper's calibration is asymmetric per capita (delta_b 0.056/0.039, G 0.204/0.191, exposures 14.4%/6.5%), so this solve must be extended before the paper calibration can be solved.] |
| `bank.py` | GK/Bocola bank block at the steady state: `calibrate_bank_targets` (single λ and ω_ent from leverage + spread), `steady_state_bank` (bond prices, net worth, portfolio, dividends). The dynamic bank problem lives in `solver_recursive/point_map.py`; only D is risky, F bonds are safe. |
| `government.py` | HM perpetuity bonds, Bohn rule: `govt_steady_state` (Bohn coefficient solved from a target debt root). The per-period fiscal block is in `point_map._government`. Default risk is exogenous (no crisis zones). |
| `solver_recursive/point_map.py` | The per-point period map (image of the old stacked system): the 14 market-clearing, Euler and TPI residuals at one grid point [corrected 2026-10-26]: was "13 residuals"; the collocation system adds the six identity equations, 20 per point per regime) given the frozen continuation rules. Bocola closed-form μ; quadrature over the s-innovation × the default regime; Z_D read from the state. |
| `solver_recursive/state_grid.py` | Smolyak sparse grid + Chebyshev basis, with `refine=(dim, m)` for a dense tensor factor on one dimension; `build_state_box`, `default_prob`, `s_process_params`. |
| `solver_recursive/collocation.py` | THE SOLVER. `make_residual` (the global F(theta), image of `residual_model.m`), `parsolve` (port of his damped FD Newton), `krylov_solve`, `solve_collocation`. |
| `solver_recursive/recursive_main.py`, `recursive_experiment.py`, `tpi_experiment.py` | Time iteration (warm start only) + SS anchors; the risk + TFP experiments and the solve ladder (with the TPI cap homotopy); the TPI against no TPI (rest point, IRF with footprint, bond-price legs, price impact, KT check). |
| `fast_kernels.py` | numba kernels for EGM backward + distribution forward; exact numpy fallback when numba is absent (`cal["use_numba"]`). |
| `household.py`, `distribution.py` | EGM with GHH utility; stationary distribution and forward iteration. |
| `trade.py` | CES basket and bilateral flows with PER-COUNTRY home bias and the country-mass ratio (`size_ratio`); `omega_home_F` is derived from `omega_home_D` and the sizes so trade balances at p = 1. |
| `firms.py`, `capital.py` | Flexible-price production with the Neumeyer-Perri working-capital wedge (w ÷ (1+ζ·r_wc), the spread→output channel; ζ=0 nests exactly — Bocola §V.C's own open-economy fix), Jermann adjustment costs, CES/Armington trade. |
| `prints.py` | Console reporting: `banner`, `print_ss_table`, and THE UNIT CONVENTION (`bp_ann`, `ann_pct`, `ann_prob`, and the `BOCOLA_IRF_*` benchmarks). Rates are annualised bp; p^d is printed quarterly AND annual; flow responses in level % with a ×4 annualised companion — Bocola's Table 5 unit. |
| `plots.py` | Every global figure; written to `plots.OUTDIR` (default `results/GLOBAL/figures/`, set by `run.py`). |
| `reporting/export.py` [corrected 2026-10-26] | Writes `results/GLOBAL/data/comparison_irfs.json` for `compare.py` (the SSJ twin is `linear_ssj/reporting/export.py`). |
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

**Acceptance thresholds** (checked against the test files 2026-10-26; only the first item and the notes marked as asserted are enforced):
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
  [corrected 2026-10-26] These tolerances, `tol_transition` and the zero-shock transition
  test below belonged to the deleted perfect-foresight transition solver. In the collocation
  solver the D goods market is one of the 20 equations held to `TOL_MAXF`;
  `test_recursive_nesting` asserts that the SS is a rest point to 1e-6 (N1), that a purchase is a
  swap within the period to 1e-12 (N3) and that the union budget identity holds to 1e-8 (N4).
  goods_F is reported, about 1e-5 of GDP off the SS after the 2026-10-06 household fix.
- (Historical, perfect-foresight era) Zero-shock transition stays at SS to ≤ 1e−5.
- Risk-only shock (exogenous π), a DIAGNOSTIC EXPECTATION, not an assertion: the D bond price,
  both banks' net worth and D output and consumption on impact are expected to fall, and the
  lending spread, D debt and taxes to rise. A response of the opposite sign must be investigated
  before it is reported: it can be an accounting error, a numerical artefact near the KKT kink, or
  a genuine general-equilibrium offset (deposit-rate parity, the backstop's transfer). In this
  project's history it has usually been an error (representative-branch pricing,
  `wc_rebate = 1`, the delta_b = 0.25 interlude). No test enforces these signs.
  [corrected 2026-10-26] (was "a positive Y or n response to sovereign risk = bug")
- Complementarity on every solved path: μ ≥ 0, slack = αn − λ·assets ≥ 0,
  μ·slack ≈ 0 (`out["mu_D/F"]`, `out["slack_D/F"]` from point_map.py). μ is computed in closed
  form, so this holds by construction at the nodes; between nodes the fitted μ can violate it
  (see the KKT-kink limitation in `global_modelling_choices.md`).
- (Historical open item, perfect-foresight era, not re-checked under collocation) risk-on
  n_D[0] can sit above risk-off (M1 deposit-rate channel; test warning, not assert) and
  post-impact Y_D runs mildly positive — both die with the union deposit
  market (§8 of the retired sunspot_transition_study.md; git history).
