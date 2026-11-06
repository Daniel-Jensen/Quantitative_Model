# The TPI solve against Bocola (2016): method check, the 12% box, and the experiment at Bocola's scale (2026-11-03)

Global solver only (`global_projection/`), branch `calibration`. Fix A (reserves not divertable) in; recapitalisation floor
50%; cap target 200 bp/yr (run.py's TPI_CAP).

## 1. Do we compute what Bocola computes?

Bocola's replication package is in the repository (`replication/bocola2016`).

| | Bocola (2016) | Our global solver |
|---|---|---|
| Approximation | Smolyak sparse grid, Chebyshev (Grey Gordon's code), isotropic level 3 in 6 states, ~389 points | Smolyak/Chebyshev; level 2 in 12 states plus a 9-node risk axis crossed with the other edges ("cross", 75 points) |
| Unknowns | 4 rules × 2 regimes × ~389 points ≈ 3,100 | 20 rules × 2 regimes × 75 points = 3,000 |
| Expectations | Gauss–Hermite; next-period states clipped to the box | same |
| Solver | Newton (`parsolve`) on all collocation residuals | same, with an exact complex-step Jacobian and an LM fallback |
| Default | two regimes priced with p^d; haircut walked 0 → 0.55 | same haircut ladder |
| Bank constraint | genuine max of the closed-form multiplier | same (hard polish) |
| Net-worth floor | `max(N, 0.65)` in every state, unfunded | funded, default regime only, 50% of n_ss |
| Risk range | s ∈ [−11.41, −2.71]: p^d 0.001%–6.2% a quarter | 0.003%–12% |
| Device at extreme points | "tapering": haircut / 1.8 at 71 collocation points, walked back to the full haircut | none |
| Headline shock | p^d → 2.5% a quarter (`generate_irf.m`) | 2% |
| Purchases | none (his LTRO is a perfect-foresight experiment) | TPI cap rule (complementarity) |

- **Same method.** The approach is the same: projection on a Smolyak grid, Newton on the collocation residuals, continuation in
  the haircut. The system size is the same too.
- **Two differences matter here.**
  - Our risk box is twice his in p^d. The extra range, 6–12%, is where the TPI solve fails.
  - Bocola himself softens the default loss at extreme collocation points.
- **Not copied:** a literal copy of his grid (level 3 in our 12 states, ~2,650 points, ~106,000 unknowns) is not feasible.

## 2. Why the TPI walk stalls on the 12% box (measured)

**The economics have a solution.**
- **Each binding point solves alone.** Solved alone, every binding grid point has a clean root at every cap down to 200 bp.
  - At 200 bp the central bank buys the D banks' whole book: about 100% of SS debt, with the D banks at their corner.
  - The market clears, and purchases can buy everything.
- **No fold.** The exact Jacobian at solved caps is well posed: condition number 3.6e6 at 1075 bp and 3.8e6 at 1002.5 bp.
  - The smallest singular value does not shrink as the cap tightens.
  - The weak directions are the franchise values α_D, α_F at the capital-edge points, not the TPI unknowns.

**The obstacle is numerical: a band of simultaneous switches.**
- **Where the 23 points come from.** The cross grid crosses its top risk level with every other state's edges, so the top
  level carries 23 points with nearly identical no-TPI spreads.
  - On the 12% box their switch caps lie within ~1150–950 bp.
  - Some lie within 0.1 bp of each other (1001.0, 999.7, 999.6, 998.7 bp at the 1002.5 bp solution).
- **Two failure modes:**
  - A cap step that lands on a switch leaves both arguments of the purchase complementarity at zero at the solution, and
    Newton plateaus at ~1e-5.
  - A step across several switches starts too far from the solution.
- **Remedies tried and rejected,** each measured on the fine grid at 50%:

| remedy | result |
|---|---|
| TPI smoothing 1e-2 during the walk | spurious purchases at 25+ points; stalls at 1025 bp (2.1e-3) |
| TPI smoothing 1e-3 | stalls on the same-cap re-solve (3.7e-5) |
| Gauss–Seidel point sweeps (coarse grid) | oscillate (3.7e-5 → 8.7e-5) |
| Garcia–Zangwill form in the walk | ~1e-5 plateau at 999.5 bp |
| direct solve at 200 bp from corner guesses | starts at 0.57, steps of 1/64 |
| pre-solve of switching points | inconsistent start (4.5e-3) |
| active-set warm start | 5.7e-4 after 10 iterations at 997.7 bp |
| small steps, path tolerance 3e-5, 25 iterations | each switch costs 1–3 h; some steps still fail |

- **Walk reached:** cap 1250 → 1002.5 bp. Rung pickles: `runs/tpi_gs_fine50_sharp_*.pkl`.
- **Below the band** only two single-point switches remain before 200 bp (the 9% and 3.9% points, at ~809 and ~446 bp).
- **The coarse grid cannot answer the question.** With 3 risk points its collocation system folds just after the cap first
  binds (~1140 bp, the same at every floor tried).

## 3. The experiment at Bocola's scale

**Design** (`bocola_size.py`):
- the risk box cut to Bocola's top, `pd_box_top` 0.062: p^d 0.0032%–6.20%;
- the same fine cross grid: 9 risk nodes, as in his axis resolution, and ~3,000 unknowns, as in his system;
- his 2.5% shock, alongside our 2%;
- the 50% floor and fix A;
- the no-TPI baseline re-solved on that box, seeded from the 12% rules: 2.4e-10, 50 min.

### 3.1 The no-TPI economy on Bocola's box (exact reads)

| | 12% box | Bocola's 6.2% box |
|---|---|---|
| Rest point: μ_D, sovereign spread, E[Ω] | 0, 62 bp, 1.117 | 0, 69 bp, 1.119 |
| D constraint binds from | p^d 3.93% | **never inside the box** |
| 2% shock: sovereign spread, lending spread | 264 bp, 0 | 239 bp, 0 |
| 2% shock: Q_bD, n_D, n_F | −7.7%, −11.0%, +0.03% | −6.6%, −6.2%, −0.16% |
| 2% shock: Y_D, C_D | −0.031%, +0.12% | +0.011%, +0.055% |
| 2.5% shock (Bocola's): sovereign, lending spread | — | 275 bp, 0 |
| 2.5% shock: Q_bD, n_D, n_F, Y_D, C_D | — | −7.8%, −7.5%, −0.19%, +0.013%, +0.067% |
| 5% shock: sovereign, lending spread, n_D | 496 bp, +173 bp, −21% | 439 bp, 0, −12.8% |

- **No pass-through at Bocola's scale.** On his box the D constraint never binds, so the lending spread does not move and
  output rises slightly at every shock size. F bank net worth now falls (correct sign); D consumption still rises.
- **Why the box matters.** The 12% box's binding from 3.9% depends on having the deep-distress states in the grid. With the
  top at 6.2%, next-period risk beyond the top is clipped to the edge, and banks no longer anticipate those states.
  - This is the same box dependence found before: the threshold moved 1.3% → 3.9% across boxes.
  - The occasionally binding constraint is not robust to the box.

### 3.2 Where a 200 bp cap binds on Bocola's box

| p^d level | grid points | no-TPI sovereign spread | at a 200 bp cap |
|---|---|---|---|
| 6.2% (top) | 23 | 444–618 bp | buying |
| 4.7% | 1 | 473 bp | buying |
| 2.1% | 1 | 284 bp | buying |
| 0.6% | 1 | 148 bp | not buying |
| 0.15% | 23 | 20–160 bp | not buying |
| ≤ 0.03% | 26 | −27 to 136 bp | not buying |

- **A wide band.** The top band's switch caps spread over 174 bp, against a few bp on the 12% box. Their spacing is 5–19 bp
  (618.0, 599.0, 589.8, 582.1, 576.5, 572.7, …), so the gap-aware walk can take them one at a time.
- **Nothing near its switch at the target.** At 200 bp the nearest non-buying point is at 160 bp.

### 3.3 The TPI walk on Bocola's box (running)

`tpi_gs_walk.py tpi_start_bocola.pkl --gaps --margin 1 --target 200 --step 10 --grow 1.25 --maxit 12 --rung-tol 3e-5`

| cap | result |
|---|---|
| 640 bp | solved, 9e-9 (no point buying) |
| 627.5 bp | solved, 2.5e-10 |
| 611.9 bp | **first switch crossed, solved to 1.4e-13** (one point buying 2.2% of D debt) |
| 592.3 bp | running |

**Status (2026-11-04).**
- **Steps solved.** A run of the walk integrated in run.py solved 644 → 602.6 bp (5 steps; the first switch, at 618 bp,
  crossed in 8 iterations).
- **Where it struggled.** The 587.7 bp step, which crossed two switches at once, stopped at 2.3e-4.
- **That run was stopped** so that the full computation is generated by `python3 run.py` itself
  (`Claude files/docs/RUN_GLOBAL.md`). Its checkpoints were deleted, so a click starts from scratch.

## 4. What this means

- **Not an economic failure.** The TPI walk is solvable on Bocola's box, one switch at a time. Its cost is numerical:
  ~25 switch events at 2–15 Newton iterations each, 6–8 min per full iteration on the 75-point grid. On the 12% box a
  23-point cluster of switch caps within a few bp makes it far harder.
- **A different economy on Bocola's range.** The D constraint does not bind at central states, so the Bocola
  pass-through to lending is absent at his scale in this model. The constraint binds only in stressed edge states,
  where the TPI's purchases relieve it.
- **The structural way to a fast TPI solve is a smooth purchase rule** (SPEC.md's closed-loop rule with an aggressiveness
  parameter, as in the SSJ model) instead of the hard cap. It has no switch events. It is a modelling choice, not yet
  made.
