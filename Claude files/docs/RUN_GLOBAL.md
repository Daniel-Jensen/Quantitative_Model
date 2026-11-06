# Running the global model and the TPI with run.py (2026-11-04)

## One click

Open `run.py` and press Run, or type `python3 run.py` in a terminal at the repository root. With the settings at the top
of `run.py` (`MODEL = "global"`, `QUICK = False`), that one command computes every global result, the TPI included,
and then draws every figure. Nothing has to be run by hand beforehand.

The terminal first prints the list of steps below, then announces each step as it starts. When each step ends it prints
`step k/7 done in … min, … min since the start`. Everything printed is also written to `results/GLOBAL/run.log`.

| Step | What it does | Typical time on this Mac (8 cores) | Saved to |
|---|---|---|---|
| 1/7 | calibration and the deterministic steady state | ~1 min | — |
| 2/7 | the TFP-shock solve | ~10 min (seconds when cached) | `results/GLOBAL_cache/` |
| 3/7 | the sovereign-risk solve without the TPI (75-point fine grid) | ~50 min (seconds when cached) | `results/GLOBAL_cache/` |
| 4/7 | rest point, impulse responses, calibration table, decompositions, Euler-error accuracy; every readout clears the period map exactly at each visited state (2026-11-06) | ~2 min | — |
| 5/7 | saves every no-TPI result | seconds | `results/GLOBAL/data/` |
| 6/7 | the TPI. `TPI_MODE` "relief" (now): the relief purchase at the shock against no TPI, on the no-TPI rules. "cap": the cap walk down to `TPI_CAP`, then the comparison | ~10 min ("cap": many hours, see below) | `results/GLOBAL/data/tpi*` ("cap" also `results/GLOBAL_tpi_walk/`) |
| 7/7 | every figure, drawn from the saved data | ~1 min | `results/GLOBAL/figures/` |

## Stopping and resuming

- **If the run stops** (the Mac sleeps, the window closes, power runs out), click Run again.
  - Steps 2–3 reuse the cached solves (`results/GLOBAL_cache/`, keyed on the calibration and the solve settings).
  - The TPI walk resumes from the lowest cap it had reached (`results/GLOBAL_tpi_walk/`, keyed on the economy).
  - Steps 1, 4 and 5 are always recomputed (~20 min).
- **To start everything from scratch,** delete `results/GLOBAL_cache/` and `results/GLOBAL_tpi_walk/`.
- **A changed calibration or setting** gets a new key, so it is solved fresh automatically. Old cache files are simply
  not used.
- **Keep the Mac awake.** run.py stops idle sleep while it runs (`caffeinate`), but closing the lid still sleeps a Mac
  without an external display. On power, set System Settings → Battery → Options → "Prevent automatic sleeping on power
  adapter when the display is off".

## The settings that matter (top of run.py)

| Setting | Now | Meaning |
|---|---|---|
| `RISK_SHOCK_PD` | 0.06 | the headline risk shock, p^d per quarter on impact (2026-11-06, a crisis-sized shock: the D banks' constraint binds on impact from ~5% on the 12% box; at the old 0.02 it is slack along the whole path and there is no recession) |
| `PD_BOX_TOP` | 0.12 | top of the risk grid, p^d per quarter. 0.12 = the appendix-B box, `WARM_START`'s own and the only one on which the D constraint binds on impact. 0.062 = Bocola (2016)'s range (it never binds there); 0.02 = the illustrative 200 bp cap of 2026-11-05/06 (`Claude files/reports/tpi_root_cause_2026-11-05.md` §6–8). |
| `TPI_M_BAND` | 0.32 | half-width of the Eurosystem-book box, a share of the D banks' SS book: 0.32 = `WARM_START`'s own band, ±27.9% of the D debt, which holds the 6% shock's relief purchase (11.8%). 1.0 was the 200 bp cap's band on the 2% box; wider than 1.0 the D bank's book at the M+ edge turns negative. |
| `RUN_TPI` | True | run the TPI experiment and compare it with no TPI |
| `TPI_MODE` | "relief" | "relief": at the shock, unanticipated, the Eurosystem buys just enough D bonds to switch the D banks' constraint off and holds them to maturity; the no-TPI rules continue (`tpi_experiment.intervention`, no TPI solve). "cap": the spread cap `TPI_CAP` solved as an anticipated rule (`TPI_WALK`); a 200 bp cap has no solution at crisis shocks, where the spread left after relief is expected default loss. |
| `TPI_CAP` | 200.0 | the cap on the D-over-F sovereign spread, bp/yr |
| `TPI_WALK` | direct, eps_start 1e-2 (walk: step 6, grow 1.2, margin 3, path_tol 1e-4, maxit 15, reuse 3, rescue 3) | direct: the cap solved in one go (nodes above it solved alone, then the smoothing walked 1e-2 → 1e-4); otherwise the cap walk below |
| `TPI_CHECKPOINTS` | `results/GLOBAL_tpi_walk` | where each solved TPI step is kept; None = no resume |
| `SOLVE_CACHE` | `results/GLOBAL_cache` | where the solved TFP and risk rules are kept; None = always re-solve |
| `WARM_START` | `data/global_warm_start.pkl` | the solved appendix-B rules every solve starts from (seeded onto `PD_BOX_TOP` when the boxes differ) |
| `NW_FLOOR` | 0.5 | the funded recapitalisation floor in default, a share of steady-state bank net worth |

The model's parameters live in `calibration/global_projection.py`, not in run.py.

## Why the TPI takes hours

**Since 2026-11-05 run.py solves the cap directly (`TPI_WALK` direct), not by the walk described here.** The
nodes above the cap are solved alone, then the global Newton runs with the complementarities' smoothing walked from
1e-2 down to 1e-4, with no switch crossings. Set `direct=False` to get the walk back.

**The size of one Newton iteration.**
- **The system.** The global model is 3,000 unknowns: 20 decisions, 2 regimes and 75 grid points. They must hold
  together, because tomorrow at any point depends on the decisions at every point.
- **The Jacobian.** One Newton iteration needs the response of every equation to every unknown, which means re-evaluating
  the model 3,000 times.
- **The cost.** About 6–8 minutes per iteration on 8 cores.

**Why the TPI needs so many iterations.**
- **The cap is an on/off rule at each grid point.** Either the central bank does not buy (the price is above the floor),
  or it buys exactly enough to hold the floor.
- **The corners stall Newton.** Each grid point that switches from "not buying" to "buying" is a corner, and Newton's
  straight-line steps stall at corners.
- **One jump fails.** Imposing the 200 bp cap at once switches too many grid points together.
- **So the cap is walked down a few bp at a time.** On Bocola's box about 25 grid points switch between ~620 and 200 bp,
  and each switch costs 2–15 iterations. That is 100–300 iterations, or 10–40 hours.

**How the walk is kept as short as possible** (`recursive_experiment._solve_tpi`, settings `TPI_WALK`):
- **The start.** The walk begins just above the highest no-TPI grid spread; above it the TPI does nothing.
- **Steps between switches.** Before each step it reads the cap at which every not-yet-buying grid point switches, and
  never stops within `margin` bp of one. A step that lands on a switch plateaus at about 1e-5.
- **A fast first pass.** Each step first runs Newton with the Jacobian rebuilt only every `reuse` iterations (about 2 min
  per iteration instead of 8).
- **A looser tolerance on the way.** A step on the way down only seeds the next one, so it is accepted at `path_tol`.
  The final cap needs the full accuracy.
- **Giving up early.** If the fast pass ends far from the path tolerance, the step is abandoned at once and the step size
  halved. The full Newton runs only at the final cap, or when the pass ended within `rescue × path_tol`.
- **Step size.** After a success the step grows by `grow`; after a failure it halves.
- **Checkpoints.** Every accepted step is saved, so nothing is lost when the run stops.

**The 12% box (`PD_BOX_TOP = 0.12`) is much harder.** Its top risk level has 23 grid points whose switch caps lie within
a few bp of each other, some 0.1 bp apart, between ~1150 and 950 bp. The walk then needs hours per bp in that band.
Details and every remedy tried: `Claude files/reports/tpi_bocola_scale_2026-11-03.md`.

## What the results mean on Bocola's box

- **No pass-through to lending.** On the 6.2% box the D bank's incentive constraint does not bind at central states. So
  the 2% risk shock moves the D bond's expected excess return (+237 bp) but not the lending spread, and output barely
  moves.
- **The constraint binds only in stressed edge states,** with high bank debt or capital at 6.2% risk. Those are where
  the TPI's purchases bite.
- **The 12% box behaves differently.** There the constraint binds from p^d 3.9%, so the pass-through exists.

## Other commands

```bash
python3 run.py global --plots-only   # redraw the figures from results/GLOBAL/data, no solving
python3 run.py global --quick        # the coarse grid: a fast preview, not converged
```
