# Minimal fixes for the global solver, mapped onto the repository's code

All findings come from reading the code on branch `calibration` (working tree, uncommitted) on
2026-10-22. No solver was run and no repository file was edited. "Source" links are repo-relative
file paths with line numbers as they stand in the working tree. The only thing executed was a
two-line scipy check of `lu_factor`/`lu_solve` on a rectangular matrix (scipy 1.15.3), noted under Q1.
Line-of-code (LOC) counts, risks and expected effects marked **estimate** are my own judgement,
not measurements.

Common facts used below:
- The coarse grid is μ=1 Smolyak in 12 states: 1 + 2·12 = 25 points. The production "axis"
  refinement adds 2 s-axis points, giving 27 at 5 nodes in s. The tensor refinement has 115 points —
  [test_state_grid.py:139-164](Claude%20files/global_projection/tests/test_state_grid.py), [state_grid.py:82-98](global_projection/solver_recursive/state_grid.py).
- Each point has 20 stored rules (14 `SOLVE` + 6 `DERIVED`: alpha, C, r_wc per country) × 2 regimes.
  So the refined solve has 20·2·27 = 1080 unknowns — [decision_rules.py:5-12](global_projection/solver_recursive/decision_rules.py), [collocation.py:15-21](global_projection/solver_recursive/collocation.py).
- The solve ladder runs in this order: `_walk_start` (cold start at the `WALK_FROM` instruments, 25 bp
  spread) → `_solve_baseline` (d0 with the WC/FX ladders → `HAIRCUT_LADDER` → joint →
  `_walk_mu_eps`) → `_walk_instruments` → joint polish (or `_solve_tpi`) → `_refine_s` —
  [recursive_experiment.py:466-523](global_projection/solver_recursive/recursive_experiment.py).

---

## 1. Least-squares projection (more residual points than coefficients), reusing parsolve's LM

### Takeaway
This is the most self-contained solver change. One new residual builder (~40 lines), plus a
fix to `parsolve` for rectangular Jacobians and a new stopping/acceptance rule (~20 lines), plus
plumbing (~10–20 lines): **estimate 70–100 LOC**, with the default unchanged. One trap: today a
rectangular Jacobian would crash `parsolve`. Also, every downstream "converged" test assumes
max|F| ≤ 1e-9, which an over-determined system cannot reach.

### Cited Findings
- The unknowns are the rule VALUES at the nodes. Each `F(theta)` call unpacks theta and refits
  every rule exactly on the square basis (`_install_guess` → `RuleSet.set_values` → `grid.fit`, an LU
  solve). It then evaluates `_point_block` at every node `i` of `rs.grid.points`, reading the
  unknowns from `vals[k][d][i]`. — [collocation.py:62-111](global_projection/solver_recursive/collocation.py), [decision_rules.py:73-80](global_projection/solver_recursive/decision_rules.py), [state_grid.py:186-188](global_projection/solver_recursive/state_grid.py)
- `_point_block` builds the 14 period-map residuals plus 6 identity residuals
  `to_fit(guess) - to_fit(implied)` for the DERIVED rules. — [collocation.py:69-84](global_projection/solver_recursive/collocation.py)
- A weighted ridge fit already exists (`SmolyakGrid.fit_weighted`, normal equations with a
  degree-weighted ridge). Only time iteration uses it, through `cal["fit_mask"]`/`fit_ridge`. — [state_grid.py:190-202](global_projection/solver_recursive/state_grid.py), [recursive_main.py:175-183, 200-211](global_projection/solver_recursive/recursive_main.py)
- `make_residual` already accepts a per-equation `scale` vector (equation weights, unused by
  callers). — [collocation.py:87-90](global_projection/solver_recursive/collocation.py)
- `parsolve`: `_factor` calls `lu_factor(J)`, and `_newton_direction` calls `lu_solve`, falling back to
  `lstsq` only on `LinAlgError`/`ValueError`. Measured on scipy 1.15.3: `lu_factor` ACCEPTS a 6×3
  matrix, and `lu_solve` then raises `_flapack.error`, a plain `Exception` subclass. That error is
  **not caught**, so a rectangular J would crash the Newton step. — [collocation.py:180-195](global_projection/solver_recursive/collocation.py) (scipy check run in this session)
- `_lm_step` uses a thin SVD (`full_matrices=False`) on the column-scaled J, which works unchanged for M > N. — [collocation.py:214-227](global_projection/solver_recursive/collocation.py)
- Stopping: `tol = m * tol**2` with `m = x0.size`, and a stall is "converged" only if
  `max|F| <= floor_tol` (= 10·TOL_MAXF). — [collocation.py:230-267, 311-318](global_projection/solver_recursive/collocation.py)
- The FD Jacobian (serial and `_PoolJacobian`) is shape-agnostic (`J = np.empty((f.size, x.size))`). — [collocation.py:120-171](global_projection/solver_recursive/collocation.py)
- Downstream acceptance: `_stamp_verdict` sets `rules.solve_ok`, and `calibrate_stochastic.calibrate`
  rejects any evaluation with `solve_ok=False`. — [recursive_experiment.py:453-463](global_projection/solver_recursive/recursive_experiment.py), [calibrate_stochastic.py:116-126](global_projection/solver_recursive/calibrate_stochastic.py)
- A ready-made denser node set exists: the tensor refinement `refine=(IS, 5)` (115 points) contains
  the axis-grid nodes, since its sparse factor includes the centre. — [state_grid.py:101-113, 135-148](global_projection/solver_recursive/state_grid.py), [test_state_grid.py:161-163](Claude%20files/global_projection/tests/test_state_grid.py)
- History: at leverage 8 the tensor (115-point) COLLOCATION "has no collocation solution near the
  coarse one" (smallest singular value of the scaled Jacobian 5e-6 against 4.2). The axis grid was
  adopted because of this. — [calibration_worklog_2026-10.md:56-71](Claude%20files/docs/calibration_worklog_2026-10.md)

### Where it goes (estimate)
- `collocation.py`: a new `make_residual_ls(rules, ..., eval_pts)`. Same unpack/install as `F`, then
  for each regime `rs.eval_all(d, eval_pts)` gives x and the DERIVED guesses at the M eval points,
  then `point_residuals(eval_pts[i], d, x_i, rs, ...)`. Most of `_point_block` can be reused by
  passing (S, x, guesses) instead of an index i. ~40 LOC.
- `_factor`: `return None if J.shape[0] != J.shape[1] else lu_factor(J)` (1 LOC). This makes the step
  Gauss-Newton via `lstsq`, with LM as the fallback.
- `parsolve`: a least-squares stop, e.g. a relative drop in `gap` < 1e-10 or ‖Jᵀf‖ small, behind a
  flag. ~10–15 LOC. `solve_collocation`: an `eval_grid` argument and an acceptance criterion
  (e.g. max|F| at the eval points ≤ a looser tol, or Euler-error based). ~10 LOC.
- `recursive_experiment._stage` / `_refine_s`: build the eval set with `normal_box(..., refine=(IS, 5, "tensor"))`
  and pass its `.points`. ~10 LOC.
- Cost: FD columns stay at 1080 (the unknowns), but each column evaluates 115 instead of 27 points: ~4.3× per Jacobian.

### Tests that would check it
- `test_collocation.py`: `test_pack_roundtrip`, `test_lm_fallback_rescues_a_failed_line_search`,
  `test_identity_residuals_are_the_readoffs`, `test_solve_reaches_the_floor` (d=0 solve to ~1e-9),
  `test_refined_grid_is_square_and_exact`. None exercises M > N. One new test is needed: an LS fit
  with eval = basis nodes must reproduce collocation exactly (nesting). — [test_collocation.py:30-102](Claude%20files/global_projection/tests/test_collocation.py)
- N1–N5 do not touch the solver's projection, so they are unaffected. N2 calls `solve_no_default`. — [test_recursive_nesting.py:183-209](Claude%20files/global_projection/tests/test_recursive_nesting.py)

### Inferences / risks
- Bit-identity: if LS is opt-in (eval set = basis nodes, or flag off), the default path is unchanged.
  With eval set = nodes, LS equals collocation in exact arithmetic but not bit for bit (`lstsq` vs LU).
- The residual mixes O(1) equilibrium residuals and log identity residuals. In LS their relative weight
  matters (it did not in square collocation, whose zero set is scale-invariant). Use `scale`.
- LS does not create a solution where none exists. It does replace "no root" by "best fit", which
  removes the hard failure mode. The KKT-switch regions would carry the residual mass, so residuals
  at those nodes must be reported, not hidden.
- TPI: `_tpi_form` converts x_cb AT THE NODES (`rs.vals`), and `_tpi_polish` re-roots. With LS, the
  GZ/FB conversion still acts on node values, so it works. But the FB/GZ "exact at nodes" argument
  becomes "exact at the basis nodes, approximate at the eval points". — [recursive_experiment.py:347-366](global_projection/solver_recursive/recursive_experiment.py)
- calibrate_stochastic's Broyden needs `solve_ok`. A new acceptance definition must be agreed first
  (per the user's "consult every modelling question" rule).

### Gaps
- Whether the LS fit at leverage 8 lands near the coarse solution or folds as the tensor collocation did. That cannot be read from code.

---

## 2. Dropping the TPI dimensions (M_cb, O_cb) from the no-TPI baseline

### Takeaway
This is almost free. The Smolyak builder already supports per-state levels, and a level of 0
gives a dimension a single centre node and degree 0. Setting `mu_vec = (1,)*10 + (0, 0)` makes the
grid 10-dimensional inside the unchanged 12-vector. That is **1 line of config + ~10–15 LOC** to
re-seed a full-dimension grid when the TPI solve starts from the baseline. The risk is low
**(estimate)**, but the no-TPI numbers change: the M/O edge nodes stop entering the fit, so it is
not bit-identical.

### Cited Findings
- `_multi_indices` enumerates `range(1, min(budget, mu_vec[j]) + 2)`. With `mu_vec[j] = 0` that is
  level 1 only, i.e. the single node 0 with degree 0. — [state_grid.py:40-54, 21-37](global_projection/solver_recursive/state_grid.py)
- The basis multiplies `T_{deg}` per dimension, so degree 0 means the rule is constant in M and O. — [state_grid.py:173-180](global_projection/solver_recursive/state_grid.py)
- The M and O bands are symmetric around 0 (`-M_half..+M_half`, `-O_half..+O_half`), so the centre node is M = O = 0. — [state_grid.py:265-277](global_projection/solver_recursive/state_grid.py)
- `mu_vec` is already threaded through everything: `run.RISK_MU_VEC` (now `None`) → `solve_recursive(mu_vec=…)`
  → `_build_rules` → `normal_box` → `build_state_box(mu=int(max(mu_vec)), mu_vec=…)`. It is stored in
  `rules.box["mu_vec"]` and reused by `_walk_instruments` and `_refine_s`. The `refine=…"axis"` path takes `mu_vec`. — [run.py:45, 107-108](run.py), [recursive_experiment.py:151-171, 256-282, 405-422](global_projection/solver_recursive/recursive_experiment.py), [state_grid.py:82-98, 290-292](global_projection/solver_recursive/state_grid.py)
- `calibrate_stochastic.main` forwards `run.RISK_MU_VEC`. — [calibrate_stochastic.py:239-247](global_projection/solver_recursive/calibrate_stochastic.py)
- With the TPI off, `res_tpi = x_cb` forces x_cb = 0, so `m_cb = B*x_cb = 0`. The book then decays,
  `M_cb_new = (1-δ)·surv·M_lag + m_cb`. From M = O = 0 the economy never leaves the slice. — [point_map.py:123-136, 231-233, 397-402](global_projection/solver_recursive/point_map.py)
- The TPI starts from the COARSE no-TPI baseline. `run.py` passes `base[0]` to `tpi_experiment.run`, which calls
  `solve_recursive(..., base=base)`. That path does `rules = base.copy()` and goes straight to `_solve_tpi` on
  the same grid (cap walk from 700 bp). — [run.py:106-108, 160-166](run.py), [tpi_experiment.py:74-78](global_projection/solver_recursive/tpi_experiment.py), [recursive_experiment.py:491-505, 369-402](global_projection/solver_recursive/recursive_experiment.py)
- The existing re-seed helper `_seed_from(rules_fine, rules_coarse)` evaluates the solved rules at another grid's points. — [recursive_experiment.py:100-107](global_projection/solver_recursive/recursive_experiment.py)
- `solve_tfp` builds `normal_box(ss, cal, mu=mu)` without `mu_vec`, so the TFP rules also carry the TPI dims. — [recursive_experiment.py:862-868](global_projection/solver_recursive/recursive_experiment.py)
- Doc claim: the TPI-off 12-state model nests the 10-state one "to solver tolerance (risk rules 5.5e-10, IRFs 4.5e-9)". This was measured 2026-10-02 on the earlier calibration. — [STATE.md:2482-2484](Claude%20files/docs/STATE.md)

### Where it goes (estimate)
- `run.py:45`: `RISK_MU_VEC = (1,)*10 + (0, 0)` (or set it in `solve_recursive` when `not tpi`).
- `recursive_experiment.solve_recursive`, `else` (base reuse) branch: if `tpi` and
  `box["mu_vec"]` has zeros at `IM, IO`, build `normal_box(ss, cal, mu=1, mu_vec=full, at=box["at"])`,
  then `rules = _seed_from(RuleSet(grid, nreg), base)` and `rules.box = dict(box, mu_vec=full)`. ~10 LOC.
  The seed is exact on the slice and flat in M/O. The first cap rung (700 bp, nothing bought) then
  solves for the M/O dependence, as the original 10→12 extension did.
- Optional: the same `mu_vec` in `solve_tfp` (1 LOC).

### Tests that would check it
- `test_state_grid.py` has no mu_vec-with-zeros test (it tests counts, axis refine, box `at`). Add one:
  12-d with zeros at IM/IO gives 21 points, and `axis` refine gives 23. — [test_state_grid.py:15-200](Claude%20files/global_projection/tests/test_state_grid.py)
- N3/N4 build their own `normal_box(ss, cal, mu=1)` with TPI on, so they are unaffected. — [test_recursive_nesting.py:59-66](Claude%20files/global_projection/tests/test_recursive_nesting.py)
- `test_collocation.test_pack_roundtrip` asserts `theta.size == len(STORE_RULES)*regimes*grid.n`, which holds for any grid.

### Inferences / risks
- Unknowns drop from 1080 to 20·2·23 = 920 (refined) and the Jacobian cost by ~27% **(estimate, ∝ n²)**.
- The corner nodes at M = ±M_half are removed. M = −M_half is unphysical: a negative Eurosystem book.
- Risk to the TPI experiment: the cap walk now starts from a seed that is flat in M/O rather than
  solved there. That is probably benign at 700 bp, but untested. `_solve_tpi` halves failed rungs.
- Pickled baselines (`rules_base_coarse.pkl`) made before the change carry `mu_vec=None`. The
  `box is None` fallback in `solve_recursive` rebuilds with the caller's `mu_vec`.
- `p_block_rotation`'s comment records that anisotropic `mu_vec` was used before (for K_D/K_F single nodes). — [recursive_main.py:133-138](global_projection/solver_recursive/recursive_main.py)

### Gaps
- Whether the leverage-8 folds sit at the M/O nodes at all. The worklog locates them at "high risk
  × bank-state edges", not at the M/O nodes — [calibration_worklog_2026-10.md:62-65](Claude%20files/docs/calibration_worklog_2026-10.md). So this fix may improve cost more than robustness.

---

## 3. Replacing the stored kinked rule alpha by a smoother object

### Takeaway
In the code, alpha is exactly `max(E[Ω]·R, lev/n)` (smoothed). The kink is the max of two
objects. Storing the two arguments (E[Ω] and the leverage ratio `lev/n`) and rebuilding alpha′ with
the same smooth max is a contained change: **estimate ~40–60 LOC across `decision_rules.py`,
`point_map.py` and tests, +2 rules (22 per point, +10% unknowns)**. Two caveats: E[Ω] is only
partly smoother (the quadrature averages over s′ alone), and leverage has a derivative kink of its own.

### Cited Findings
- The current multiplier is `mu = smax(1 - E_Om·(1+rdep)·n/lev, 0, mu_eps)` (capped at 0.95). Alpha is
  `alpha_new = E_Om·(1+rdep)/(1-mu)`, clipped to [0.05, alpha_cap=40]. — [point_map.py:330-361](global_projection/solver_recursive/point_map.py)
- **Algebra from the code:** when binding with the hard max, 1−μ = E_Om·R·n/lev, so
  alpha = lev/n. When slack (μ = 0), alpha = E_Om·R. So alpha = max(E_Om·R, lev/n).
- The continuation uses only alpha′: `Om_D[j] = Lam_D[j]·(f + (1-f)·R[j]["alpha_D"])`, where `R[j]` is
  `cont.eval_cont(j, Sn)` at the next states. — [point_map.py:307-327, 275-288](global_projection/solver_recursive/point_map.py)
- Today's alpha (`alpha_D_cur`, from the frozen rules) is used only in the slack diagnostic
  `slack_D = alpha_D_cur·n_D − lev_D`. — [point_map.py:210-215, 359-361](global_projection/solver_recursive/point_map.py)
- alpha is a DERIVED rule with an identity residual and is fitted in logs (`LOG_RULES`). — [decision_rules.py:9-21](global_projection/solver_recursive/decision_rules.py), [collocation.py:81-83](global_projection/solver_recursive/collocation.py)
- Readers use `o["alpha_D"]` (from `_outputs`, i.e. `alpha_D_new`) in `_spread_bp`, so the readers are unaffected. — [recursive_experiment.py:558-560](global_projection/solver_recursive/recursive_experiment.py), [point_map.py:527](global_projection/solver_recursive/point_map.py)
- At the SS, `calibrate_bank_targets` forces alpha_ss = λθ, and the SS is solved consistent with the
  smoothed switch through `ic_slack_factor`. — [bank.py:42-93](global_projection/blocks/bank.py), [guide_global.md:365-369](Claude%20files/docs/guide_global.md)
- The quadrature is over the s′ innovation × regime only. The next bank states `Sn` are deterministic
  given today's choices (`_next_states` fills every quadrature row with the same K′, P′, b′). — [point_map.py:258-272, 291-304](global_projection/solver_recursive/point_map.py)
- The cold start for every rule is constant at its SS value (`RuleSet.from_ss`). New rules need SS values there. — [decision_rules.py:139-161](global_projection/solver_recursive/decision_rules.py)

### Where it goes (estimate)
- `decision_rules.DERIVED`: replace `alpha_D/F` by `EOm_D/F` and `phi_D/F` (= lev/n). Add both to `LOG_RULES`
  and their SS values to `from_ss` (`bk["Omega_ss"]`; `phi_ss` = lev/n at the SS). ~8 LOC.
- `point_map._outputs`: export `EOm_D=v.E_Om_D`, `phi_D=v.lev_D/v.n_D`, keeping `alpha_D` for readers (~4 LOC).
- `point_map._discount_kernels`: `alpha′ = R′·EOm′/(1 − smax(1 − EOm′·R′/phi′, 0, mu_eps))`, with R′ from the stored
  `rdep_D` rule (already in `R[j]`), the same caps, and phi′ floored. ~10 LOC.
- `point_map._balance_sheets`: rebuild `alpha_D_cur` the same way for the slack diagnostic (~4 LOC).
- `collocation.RES_NAMES` is derived from `DERIVED`, so it updates automatically. `recursive_main._sweep` is generic over `DERIVED`.

### Tests that would check it
- `test_collocation.test_identity_residuals_are_the_readoffs` (the DERIVED identity), `test_decision_rules.py`
  (`test_constants_are_exact`, transforms), N1 (`test_n1_ss_rest_point`: the SS must stay a rest point
  with SS-constant rules, which checks the new SS values), N2 (the grid solve), and `test_post_default` P2/P5 (generic over `STORE_RULES`).

### Inferences / risks
- The gain is partial **(estimate)**. E[Ω](S) averages alpha′ over 5 GH × 2 regime cells at ONE
  next-period bank state, so in the bank-state directions it is a sum of ten shifted kinks, not a smooth
  function. lev/n (a quantity) is continuous but generally has a derivative kink where the regime switches.
- `phi′` as a rule behaves badly where n → 0 (the default regime, where n_D sits on the 15% floor). It
  needs a floor, and a log fit of lev/n there may be steep.
- Pickled rule sets (`rules_*.pkl`) become incompatible (missing keys).
- Bit-identity is lost for every result. The SS nesting must be re-verified (N1).
- Alternative suggested by the caller, E[Ω·payoff] (the bond-FOC numerators `E_Om_payD`, …): in the
  current code these expectations are computed from the continuation each call and are never stored,
  so storing them would be a parameterised-expectations redesign (more LOC; **estimate 100+**).

### Gaps
- The smoothness of E[Ω] and lev/n across the switch is not measured anywhere in the docs.

---

## 4. Time iteration (accelerated) as main solver or selector

### Takeaway
`recursive_main.time_iteration` is a damped (0.5) fixed point. Each sweep solves the 14 unknowns
pointwise with `scipy.root(hybr)` against a frozen continuation. Its binding mode contracts at
0.990 per sweep. Anderson acceleration wraps the existing sweep (**estimate ~40–60 LOC, one new
function**). The worklog records that pointwise roots are UNIQUE at the stalled nodes. That makes
a time-iteration selector economically meaningful, but it may also simply diverge where the
collocation Newton folds.

### Cited Findings
- `time_iteration(rules, ..., damp=0.5, tol=1e-7, max_it=60, n_gh=7)` freezes `cont = rules.copy()` and calls
  `_sweep`, then `_damped_update` (a damped move plus refit; optional `fit_mask`/`fit_ridge`). It exits only on
  rule change < tol AND worst point < 1e-6. — [recursive_main.py:167-211](global_projection/solver_recursive/recursive_main.py)
- `_sweep`: each point is solved by `solve_point` (hybr from the current values, with a fallback from the SS
  guess). Failed points (fn > 1e-3) keep their old values and get fit weight 0. DERIVED are read off
  `out`. — [recursive_main.py:68-126](global_projection/solver_recursive/recursive_main.py)
- "Time iteration ... survives ONLY as the warm start ... its binding mode is the franchise-value recursion
  at 0.990 per sweep, so runs reported `max|F| = 1e-14` and 'rule-change tol not reached'
  simultaneously." — [guide_global.md:73-76](Claude%20files/docs/guide_global.md); N2's comment repeats this — [test_recursive_nesting.py:184-191](Claude%20files/global_projection/tests/test_recursive_nesting.py)
- In production it runs as a 12-sweep warm start (`WARM_SWEEPS = 12`) only in the first `_stage` call. Every
  later stage passes `warm=0` because "at leverage 8 they can carry the start onto a second configuration at
  the capital edges, one that folds before the calibrated smoothing", and the refined grid's "time-iteration
  sweeps damaged it at leverage 8". — [recursive_experiment.py:89-90, 110-118, 329-333, 431-435, 508-511](global_projection/solver_recursive/recursive_experiment.py)
- "each stalled node has a unique root with the continuation frozen (30 starts), so the fold is global." — [calibration_worklog_2026-10.md:61-62](Claude%20files/docs/calibration_worklog_2026-10.md)
- The pack/unpack helpers for a flat vector of node values already exist. — [collocation.py:35-59](global_projection/solver_recursive/collocation.py)

### Where it goes (estimate)
- `recursive_main.py`: `anderson_time_iteration(rules, ..., m=5, beta=1.0)`. Use g(θ) = `pack(new)` − θ from
  `_sweep` (in the fitted transform via `to_fit`) and a least-squares Anderson mix on the last m residual
  differences, with a safeguard that falls back to the damped step if ‖g‖ rises. ~40–60 LOC.
- Use as a selector: run it from the Newton solution. If it stays (‖g‖ ≈ 0), the Newton root is the TI-stable
  branch; if it drifts, the two disagree. ~10 LOC in `recursive_experiment._stage`.
- Cost **(estimate)**: one sweep ≈ 54 point roots × (FD Jacobian of 14 + iterations) ≈ 1–3k point
  evaluations. One collocation FD Jacobian ≈ 1081 F calls × 54 points ≈ 58k point evaluations. Plain TI at
  0.99 needs ~ln(1e-9)/ln(0.99) ≈ 2,000 sweeps, so Anderson must get it to ~100–300 sweeps to compete.

### Tests that would check it
- None test `time_iteration` convergence directly. N2 tests the collocation route. A new test: Anderson TI on
  the d=0, π=0 block must reach the N2 solution (`solve_no_default`) to ~1e-8.

### Inferences / risks
- Low risk to existing results while it is opt-in: `_stage` keeps its current call pattern.
- The economic meaning of "selector": the TI limit is the limit of finite-horizon economies, which is the
  standard equilibrium-selection argument. With documented non-uniqueness this is arguably the
  principled choice, but that is a modelling decision for the user.
- The point solves use `n_gh=7` by default, while the Newton uses `N_GH = 5`. A selector must use the same order.
- The TPI GZ form `max(-x,0)` inside the per-point solve: N/A pointwise (x_cb is solved at the point), but the
  hybr solve at TPI corner nodes may fail and be masked (retained values).

### Gaps
- No measurement of Anderson on this model exists in the repo. Its convergence in the presence of masked (failed) points is unknown.

---

## 5. Centring the box on the stochastic rest point (fixed point box ↔ rest point)

### Takeaway
The hook already exists. `build_state_box(at=…)` moves the centres of the seven bank and deposit
stocks while keeping the current SS band widths, and the instrument walk already re-seeds rules onto a
moved box. A recentring loop is **estimate ~25–40 LOC in `recursive_experiment.py`**, reusing
`_seed_from`/`_stage`/`stochastic_rest_point`. The open problem is non-uniqueness: the rest point
depends on where the no-shock simulation starts.

### Cited Findings
- `build_state_box(..., at=None)`: for each key in `at` it shifts `lo` and `hi` of `_STOCK_INDEX[k]` by `v - c[k]`.
  `BOX_STOCKS = (K_D, K_F, P_D, P_F, b_DD, b_DF, b_FD)`. V_dep, s, Z_D, M, O are NOT movable. — [state_grid.py:278-283, 342-344](global_projection/solver_recursive/state_grid.py)
- `test_box_at` checks that the own centre leaves the box unchanged and that shifts move only the bank stocks. — [test_state_grid.py:167-184](Claude%20files/global_projection/tests/test_state_grid.py)
- `rules.box["at"]` travels with the rules and is honoured by `_walk_instruments` and `_refine_s`. — [recursive_experiment.py:169-171, 271, 417-422, 482-486](global_projection/solver_recursive/recursive_experiment.py)
- `stochastic_rest_point` iterates `read_at` + `advance` from `ss_state` (the deterministic SS) up to 4000
  quarters, clips to the box (`advance(..., rules.grid)`) and caches `rules._rest_point`. — [recursive_experiment.py:666-704](global_projection/solver_recursive/recursive_experiment.py)
- "On the refined grid the rest point lies on the V5 capital edge." "The stochastic rest point is not unique ... the
  slack ones sit on the box walls. The slowest eigenvalue of the linearised no-shock map is 0.987-0.994." — [calibration_worklog_2026-10.md:75-92](Claude%20files/docs/calibration_worklog_2026-10.md)
- The user decision recorded there: "the rest point defined as the one reached from the deterministic SS (the code's existing convention)". — [calibration_worklog_2026-10.md:97-100](Claude%20files/docs/calibration_worklog_2026-10.md)
- Box bands (`BOX_KW`): k 2%, P 2.26%, b 12%, w 2%, m 32%. — [recursive_experiment.py:33-34](global_projection/solver_recursive/recursive_experiment.py)

### Where it goes (estimate)
- New `_recentre(rules, cal, ss, sproc, verbose, max_rounds=4, tol=0.1)` in `recursive_experiment.py`. Each round:
  S_r = `stochastic_rest_point`; `at = {k: S_r[_STOCK_INDEX[k]]}`; grid = `normal_box(ss, cal, mu, mu_vec, at=at)`;
  trial = `_seed_from(RuleSet(grid), rules)`; `_stage(..., warm=0)`. Stop when the shift is below tol × half-width. ~25 LOC.
- Call it in `solve_recursive` after the walk / joint polish and before `_refine_s`, so the refinement inherits `box["at"]`. ~5 LOC.
- Optional: add `V_dep` (index `IV`) to `BOX_STOCKS`/`_STOCK_INDEX`, because V is ±V_half around 0 at the SS (2 LOC). This changes `test_box_at`'s set.
- For non-uniqueness: start each round's `stochastic_rest_point` from the previous S_r, not from `ss_state`.
  That needs an `S0` argument (~3 LOC) and departs from the recorded convention, so the user should decide.

### Tests that would check it
- `test_box_at`, the `test_aligned_box_identity_axes` V5 rebuild, and `test_calibrate_stochastic` (uses a synthetic evaluate, so not affected).
  None test a recentring loop. Add a test that, with at = the SS centre, the loop is a no-op.

### Inferences / risks
- Bit-identity: changes every no-TPI and TPI result. The TPI reuses `base[0]` (the coarse baseline, which
  carries `box["at"]`), so the TPI also inherits the recentred box. That is consistent but needs checking.
- The IRF needs room below the rest point (risk shocks lower capital and net worth). Centring on the rest point
  helps exactly the documented "rest point on the capital edge" problem **(estimate)**.
- Each round costs one coarse joint solve plus a refined solve if done on the fine grid. Run the rounds coarse.
- Possible cycling between rest points with different boxes (non-unique attractors). Cap the rounds and report.

### Gaps
- Whether the rest point stays put once the box moves (the box itself shapes the fitted rules) is unmeasured.

---

## 6. Treating the default regime differently (two boxes, coarser rules, bounded reads)

### Takeaway
A full two-box implementation is in the code, switched off (`TWO_BOXES = False`). It has a
post-default `RuleSet` on an aligned coarse box, a C2 blend in net-worth units, unknown packing over
both sets, and tests P1–P5. It failed on convergence: the aligned normal box's cold start, and a
post-default plateau at 7.8e-3 to 8e-3. The minimal variant is to keep the V5 normal box and make
the post-default box a plain wider `build_state_box` (not aligned): **estimate ~10–20 LOC**.
Alternatively, coarsen the post set with `mu_vec` zeros. The risk is moderate.

### Cited Findings
- `TWO_BOXES = False` with the reason: "on the coarse grid its post-default stage plateaued at 8e-3 and drove the
  normal default-regime rules off (STATE.md)". `BOX_SPECS` is loaded from `box_specs.json`. — [recursive_experiment.py:36-42](global_projection/solver_recursive/recursive_experiment.py)
- `normal_box`/`post_box`/`make_blend`/`_attach_post`: the post set is seeded from `rules.eval_cont` at its nodes
  (the edge reads). — [recursive_experiment.py:45-75](global_projection/solver_recursive/recursive_experiment.py)
- `_solve_baseline` with TWO_BOXES: attach → solve `sets=("post",)` → joint both. `_refine_s` keeps the post set coarse. — [recursive_experiment.py:334-341, 429-430](global_projection/solver_recursive/recursive_experiment.py)
- `RuleSet.post/home/blend`, `attach_post`, `blend_weight` (quintic C2 step along a recovery line), `eval_cont` (each set
  read through its own `grid.clip`, then blended) and `eval_read`. — [decision_rules.py:54-124](global_projection/solver_recursive/decision_rules.py)
- `collocation.SETS = ("normal", "post")` and `_sets` pack both sets when present. — [collocation.py:25-59](global_projection/solver_recursive/collocation.py)
- Worklog: "On the coarse grid the aligned box's cold start does not converge and the post-default stage plateaus at 7.8e-3;
  production runs on the V5 box, post-default states read at its edge." — [calibration_worklog_2026-10.md:49-52](Claude%20files/docs/calibration_worklog_2026-10.md)
- The out-of-box read is the C1 saturation `clip` with `CLIP_EPS = 0.05`: states outside are read at most 5% of a
  half-width past the edge, whatever their distance. — [state_grid.py:208-220, 249-252](global_projection/solver_recursive/state_grid.py)
- How the default regime feeds calm times: from any node, the d′=1 cell is weighted by `pd = default_prob(s)` and
  evaluates the d=1 rules at the same `Sn`. The d=1 rules' own continuation is read at their next states (post-default
  balance sheets). — [point_map.py:291-304, 48-50](global_projection/solver_recursive/point_map.py)
- The default regime's funded floor (`nw_floor_frac` = `run.NW_FLOOR` = 0.15) binds at most default nodes at leverage 8. — [recursive_experiment.py:175-178](global_projection/solver_recursive/recursive_experiment.py), [run.py:43](run.py)
- The tests `test_post_default.py` P1 (blend weight), P2 (an identical copy changes nothing), P3 (two-set packing),
  P4 (copy re-links) and P5 (the seed reproduces the edge read) all use `_attach_post` / `BOX_SPECS`. — [test_post_default.py:20-100](Claude%20files/global_projection/tests/test_post_default.py)
- The measuring tool `measure_box.py` writes the box spec from a solved rule set. — [diagnostics/state_box/measure_box.py:1-40](Claude%20files/diagnostics/state_box/measure_box.py)
- The "3–8 half-widths outside the box" magnitude comes from the caller's brief. I did not find the number in the docs I read.

### Where it goes (estimate)
- Variant A (wide V5 post box): `post_box` → `build_state_box(ss, cal, mu=1, **wide_kw)` with P_D/b_DD/K bands scaled by
  ~5–8×, keeping the blend from `BOX_SPECS["blend"]` (or a blend keyed on n_D distance). Set `TWO_BOXES = True` but keep
  `normal_box` on the V5 box: the `assert at is None` and aligned-box branch in `normal_box` must be decoupled from the
  flag. ~10–20 LOC.
- Variant B (coarser post rules): the same, with `mu_vec` zeros for s, Z_D, M, O, so the post set has very few nodes. ~5 more LOC.
- Variant C (bounded reads only): raise `CLIP_EPS` or make it per-dimension. That is 1 LOC, but it only changes how
  far past the edge reads extrapolate. It is not a fix of the content.

### Inferences / risks
- The documented failure was the ALIGNED normal box's cold start, not the post set alone. Variant A keeps the converging V5 normal box **(estimate)**.
- Variant A on a 5–8× wider box: the period map at wide post-default nodes hits the guards (`_smax` floors, `_MU_CAP`,
  `_sclip` alpha) more often, so plateaus like 7.8e-3 could recur.
- `_walk_instruments` asserts the V5 box (`assert at is None` in `normal_box` under TWO_BOXES). The walk and the two
  boxes must be made compatible. — [recursive_experiment.py:45-50](global_projection/solver_recursive/recursive_experiment.py)
- The TPI path: `_tpi_form`, `_restore` and `RuleSet.copy` already handle `rules.post`. — [recursive_experiment.py:205-212, 347-359](global_projection/solver_recursive/recursive_experiment.py), [decision_rules.py:126-137](global_projection/solver_recursive/decision_rules.py)
- Bit-identity with TWO_BOXES off is preserved (P2 asserts the blend is neutral for an identical copy).

### Gaps
- No measurement of a non-aligned wide post box exists in the docs I read.

---

## 7. Model-side small changes

### Takeaway
Three of the four are tiny edits. The floor is a run setting. p^d at rest is one hard-coded
constant. f as an instrument is one `PAIRS` tuple, though it needs a seventh target or must replace
one. The constant intermediation cost touches both the period map and the SS bank block
(**~10–15 LOC**), and it works against the "binding at rest" decision. All are modelling choices
that, under the user's standing rule, need their approval with measured numbers.

### Cited Findings
**(a) Recapitalisation floor.**
- `_funded_floor`: `floor_D = nw_floor_frac·n_ss`. The D top-up is funded by new D debt (`Bp_D += recap/Q_bD`, 8 fixed
  passes); the F top-up by an F lump-sum tax. It is smoothed by `_smax(n, floor)` with eps = 1e-3 (absolute). — [point_map.py:178-203, 28-30](global_projection/solver_recursive/point_map.py)
- The setting `run.NW_FLOOR = 0.15` is set in `compute_global`, `recursive_experiment.main` and `calibrate_stochastic.evaluate`
  (via settings). — [run.py:43, 92](run.py), [recursive_experiment.py:910](global_projection/solver_recursive/recursive_experiment.py), [calibrate_stochastic.py:72](global_projection/solver_recursive/calibrate_stochastic.py)
- Test: N5 `test_n5_funded_floor_closes_the_budget`. — [test_recursive_nesting.py:157-181](Claude%20files/global_projection/tests/test_recursive_nesting.py)

**(b) Constant intermediation cost in the loan spread.**
- The dynamic loan rate is `r_wc = rdep + λ_K·μ/E[Ω]` in `_wages_and_dividends`. It is stored as a DERIVED rule
  (`r_wc_D/F`), and `Pp_D = (1+rdep)·dep − (1+r_wc_cur)·L_wc` hands the financing income to the bank. — [point_map.py:414-423, 243-245](global_projection/solver_recursive/point_map.py)
- SS: `IC_spread_wc = λ_K·μ_ss/Ω_ss` enters `n_ss_ACCUM`. At the SS this equals rk − rdep (the `credit_spread_target`). — [bank.py:157-170, 190-209, 26-39](global_projection/blocks/bank.py)
- `RuleSet.from_ss` sets `r_wc = r_dep + credit_spread_target`. The calibration moment is `spread_D = r_wc − rdep`. — [decision_rules.py:152-154](global_projection/solver_recursive/decision_rules.py), [calibrate_stochastic.py:51-56](global_projection/solver_recursive/calibrate_stochastic.py)

**(c) Banker exit share as an instrument.**
- `f_D = f_F = 0.08` ("vs Bocola's 0.0354"). — [calibration/global_projection.py:31-32](calibration/global_projection.py)
- `_check_targets` asserts n_IC = n_ACCUM and says "fold-blocked -- raise f until leverage is the least root". — [steady_state.py:184-196](global_projection/steady_state.py)
- "The deterministic SS does not exist above a D SS-spread instrument of ~240 bp at f_D = 0.08." — [calibration_worklog_2026-10.md:78-80](Claude%20files/docs/calibration_worklog_2026-10.md)
- "f does NOT do this [move μ at rest]: `calibrate_bank_targets` forces alpha_ss = lambda*theta at the SS whatever f is." — [guide_global.md:367-369](Claude%20files/docs/guide_global.md)
- `PAIRS` matches six instruments to six targets, and Broyden assumes a square system. — [calibrate_stochastic.py:27-37, 95-146](global_projection/solver_recursive/calibrate_stochastic.py)
- The instrument walk only walks the keys of `WALK_FROM`/`WALK_STEP`. — [recursive_experiment.py:196-202, 256-263](global_projection/solver_recursive/recursive_experiment.py)

**(d) Lower priced default probability at rest.**
- `s_process_params` hard-codes `s_star = log(0.001/(1-0.001))` and `rho_s = 0.95`; `sigma_s` comes from cal (0.4455). — [state_grid.py:385-392](global_projection/solver_recursive/state_grid.py)
- The deterministic SS prices NO default: `Q_bdom_ss = δ/(rdep + δ + IC_spread)`. The SS files do not reference `s_star`. — [bank.py:157-170](global_projection/blocks/bank.py) (grep of `steady_state.py`, `blocks/*.py` found no `s_star`/`default_prob`)
- The s box is `s_star ± S_COVER_SD·σ_s/√(1−ρ²)` = ±3.92 (top p^d 4.82%/qtr at p^d_rest 0.1%). — [state_grid.py:247-248, 358-366](global_projection/solver_recursive/state_grid.py), [global_solution_theory_summary.tex:627-630](Claude%20files/docs/global_solution_theory_summary.tex)
- `test_state_box_and_s_process` asserts `default_prob(s_star) == 0.001`. — [test_state_grid.py:41-60](Claude%20files/global_projection/tests/test_state_grid.py)
- The calibration table row "Default probability, GR 1y CDS, 2008Q3 — 0.092% — Set directly". — [calibrate_stochastic.py:205-206](global_projection/solver_recursive/calibrate_stochastic.py)

### Where it goes (estimate)
- (a) 0 LOC (`run.NW_FLOOR`). Possibly also widen the `_smax` eps relative to n_ss (1 LOC).
- (b) `cal["wc_cost_D/F"]` (new key), added in `_wages_and_dividends` (2 LOC), in `_ss_bond_prices.IC_spread_wc` (1 LOC) and in
  `RuleSet.from_ss` (2 LOC). Check the SS household anchors (`calibrate_household_anchors` re-absorbs the change).
  The identity `rk − rdep = IC_spread_wc` at the SS no longer holds, so check `test_ss_identities`. ~10–15 LOC.
- (c) Add `("f_D", <target>, <moment>, tol, lo, hi)` to `PAIRS` with a new target, OR swap it in for one instrument. Add `f_D`
  to `WALK_FROM`/`WALK_STEP` if walked. ~5 LOC.
- (d) `s_star = log(p/(1-p))` with `p = cal.get("pd_rest", 0.001)` (1–2 LOC). Update the test assert (1 LOC).

### Inferences / risks
- (a) A higher floor raises default-state n_D, which shrinks the default-state recession and the calm-time risk premium. But
  it is funded by NEW D debt the banks absorb, so `b_DD′` after default goes UP. That may push the post-default b_DD read
  further outside the box even as P_D/n_D come closer **(estimate)**. It changes the economics of the default event, which
  "Bocola-pure" defines as a pure haircut.
- (b) With a constant cost carrying most of the 100 bp wedge, μ at rest is small. The documented "binding at rest" user decision
  and the 2026-08-29 cure of the kink-at-rest accuracy problem (spread 8 → 100 bp so μ_rest > 0) both depend on μ_rest being
  well away from 0. A cost term could re-open "THE ECONOMY RESTS ON THE KKT KINK" — [guide_global.md:354-374](Claude%20files/docs/guide_global.md).
  The upside is that the SS spread instrument could stay near 100 bp, below the ~125–200 bp cold-start failure **(estimate)**.
- (c) Raising f moves the ~240 bp SS existence frontier (per `_check_targets`' message), but by the guide's algebra it does
  not change μ_ss at a given spread. It helps existence and does not help binding.
- (d) Lowering p^d at rest to 0.05% lowers the box top with it. Computed from the code's formula: s_star = −7.60, s_hi = −3.68,
  so the top p^d is about 2.47%/qtr, against the `RISK_SHOCK_PD = 0.02` headline at s = −3.89. The headline then sits ~94% of the
  way to the edge (now 77%). `s_hi` or `S_COVER_SD` would need raising. — [run.py:51](run.py) (arithmetic mine)

### Gaps
- No repo measurement of the rest-point premium against `nw_floor_frac` or p^d_rest at leverage 8.

---

## 8. Calibration targets on ergodic (simulated) means instead of the no-shock rest point

### Takeaway
The pieces exist: `accuracy.simulate` (a seeded stochastic simulation conditioned on no default)
plus `read_at` plus `rest_point_moments`. Replacing `rest_point_read` in `calibrate_stochastic.evaluate`
is **estimate ~20–30 LOC**. The catches are Monte Carlo noise against Broyden's tight tolerances (0.1 bp
on spreads), and that the simulation clips to the box and visits the poorly-approximated box top.

### Cited Findings
- `calibrate_stochastic.evaluate` → `rest_point_read` → `stochastic_rest_point` + `read_at`, then
  `rest_point_moments(o)`: leverage = (dep+n)/n, spread = r_wc − rdep, exposure = Q·b/assets. — [calibrate_stochastic.py:45-56, 69-85](global_projection/solver_recursive/calibrate_stochastic.py)
- The tolerances are 0.005 leverage, 0.1/4e4 spread and 1e-4 exposure (`PAIRS`, `TOL`). `MAX_EVALS = 12`. — [calibrate_stochastic.py:27-37](global_projection/solver_recursive/calibrate_stochastic.py)
- `accuracy.simulate(rules, cal, ss, sproc, S0, T=2000, burn=200, seed=0, no_default=False)` draws s innovations with a seeded
  RNG, always evaluates regime 0 (no default realised), clips to the box, and returns the visited states and their
  out-of-box shares. — [accuracy.py:30-68](global_projection/solver_recursive/accuracy.py)
- `run.py` already calls `accuracy_report(..., S_rest, T=ACCURACY_T=1200)` on the ergodic set. — [run.py:48, 131-134](run.py)
- The worklog's rest-point non-uniqueness: "the no-shock simulation settles at different states (at 150 bp: D spread at rest 75,
  159, 53, 4, 0.1 bp)". — [calibration_worklog_2026-10.md:86-92](Claude%20files/docs/calibration_worklog_2026-10.md)
- `test_calibrate_stochastic.py` (C1–C3) replaces `CS.evaluate` with a synthetic map, so it is insensitive to how moments are computed. — [test_calibrate_stochastic.py:19-50](Claude%20files/global_projection/tests/test_calibrate_stochastic.py)

### Where it goes (estimate)
- `calibrate_stochastic.py`: `ergodic_moments(rules, cal, ss, sproc, T, burn, seed)`. Use `states, _ = accuracy.simulate(...)`, then
  `[rest_point_moments(read_at(rules, cal, ss, sproc, S)) for S in states]`, then means (or ratios of means for leverage and
  exposure). ~20 LOC. A `settings["MOMENTS"] = "rest" | "ergodic"` switch in `evaluate` (~5 LOC). The table helpers
  (`calibration_values`) stay rest-point based unless changed.
- Fixed seed, i.e. common random numbers across Broyden evaluations, so the moment map is deterministic in x.

### Inferences / risks
- MC error **(estimate)**: with ρ_s = 0.95 the effective sample is ~T·(1−ρ)/(1+ρ) ≈ T/39. With T = 2000 that is ~50 effective draws.
  The spread at rest swings by tens of bp across states, so the standard error could exceed the 0.1 bp tolerance by an order of
  magnitude. Common random numbers keep Broyden stable, but the calibrated target then holds for one seed. Loosen the tolerances or raise T.
- Semantics: "ergodic" here means conditional on no default ever realised (regime 0 always). That matches Bocola's pass-through
  design, but must be stated.
- Clipping at the box edge biases the means wherever the simulation leaves the box (accuracy reports the share).
- It partly sidesteps the non-unique no-shock rest point, since shocks mix across basins. But if basins are absorbing, the means still depend on S0.
- The IRFs still start from the no-shock rest point (`dynamic_irf`), so targets and IRF base would then differ. That is a reporting change.

### Gaps
- Cost per evaluation: T reads of `point_residuals` are cheap next to a full solve, but no timing is recorded in the repo.
