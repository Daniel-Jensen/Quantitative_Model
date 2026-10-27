# CLAUDE.md

## Research question
Two-country heterogeneous-agent model of a monetary union (D = Greece, F = Germany) with
Gertler-Karadi banks and exogenous sovereign risk, following Bocola (2016, JPE). A rise in the
priced default probability π_t, driven by a latent factor s_t and never a function of debt,
marks down long D bonds; bank net worth falls, the occasionally binding incentive constraint
tightens (multiplier μ), lending spreads rise and output falls, with no default realised. Only
D bonds are risky. Policy application: a Eurosystem spread-cap bond-purchase backstop (TPI).
Deliverable: the viva paper, Overleaf project 698b4f88aeef1d0e1d08cc0c; the local clone
`698b4f88aeef1d0e1d08cc0c/VIVA/` can lag Overleaf, so pull before editing it.

## Two solvers: always say which one a number comes from
| | Global projection | Sequence space (SSJ) |
|---|---|---|
| Code | `global_projection/` | `linear_ssj/`, `experiments/` (E1-E4), `diagnostics/regimes/` |
| Economics | nonlinear, occasionally binding constraint; Chebyshev collocation, Smolyak grid with a dense factor in s | first order; sticky prices, nominal deposits |
| Calibration | `calibration/global_projection.py` | `calibration/ssj.py` |
| Interpreter | `python3` | `/opt/anaconda3/envs/ssj/bin/python` only (base env has a broken liblapack) |
| Role | the model the viva paper describes | linear policy experiments |

They share no model code and are calibrated differently; impulse magnitudes are not comparable.

## Source of truth
1. **The paper calibration is in the code but not yet solved** (`"paper"` in
   `calibration/global_projection.py`; every existing global result is on `"bocola"`; STATE
   section 1 has both, and the open kink risk under the 5 bp wedge). Until the paper calibration
   is solved, paper model numbers stay [TBC] and no `"bocola"` number goes into the paper.
2. **Code over prose.** Parameter values come from the calibration files; mechanisms from the
   period map (`global_projection/solver_recursive/point_map.py`) or the SSJ equation files.
   A number in this file, STATE or PROGRESS is historical unless it names the calibration and
   date that produced it and both still match the code.
3. **Data numbers come with provenance.** Before quoting or changing a data-based calibration
   value, check `Claude files/docs/DATA.md`: several have no script, and some do not reproduce from the stored
   files.
4. **Results are dated.** The reference global solution is the 6 Oct 2026 run on the code
   calibration; global impulse-response levels from before it are void.
5. `experiments_results.md` and `paper_draft_results.md` are generated; never hand-edit them.

## Repository map
- `TODO.md` (root): the one list of open work, for both authors, in four parts: Code, Paper,
  VIVA, Post-VIVA. At the start of every session, delete any ticked item (`- [x]`) and record
  it in PROGRESS with the date and PR or commit.
- The documents written for Claude are in `Claude files/docs/`: `STATE.md` (current results),
  `PROGRESS.md` (changelog), `SPEC.md` (research question, model, policy experiment, framing),
  `GLOSSARY.md` (symbols, code names, acronyms, branches), `DATA.md` (data sources, scripts and
  provenance of every data number). This file stays at the root because Claude Code loads it
  from there.
- `run.py`: the single entry point; holds the global pipeline's settings (grid, shocks, TPI).
  `compare.py`: the two solvers' IRFs side by side.
- `global_projection/`: `blocks/` (economic blocks), `solver_recursive/` (collocation solver,
  experiments, TPI), `reporting/`, `steady_state.py`.
- `linear_ssj/`: `model/` (equation files; `full_model.build_block_list()` is the single model
  definition), `solve/`, `reporting/`, `main.py`.
- `Empirics/`: data work; `Empirics/outputs/` is read by code.
- `results/<SSJ|GLOBAL|COMPARISON>/`: git-ignored outputs (`POLICY/` and `WELFARE/` belong to
  shelved work).
- `Claude files/`: every AI-written doc, note, report and test, at the path it would otherwise
  have (model code stays in its packages even when Claude wrote it). Quote the space in shell
  commands.

## Editing safeguards
- **No model run without an explicit go-ahead in this conversation**: no solve, trial,
  `--quick` preview or Monte Carlo. A brief or plan that says "re-solve" is not authorisation.
- Never edit model code while a run is in progress (the Jacobian workers re-import it).
- A compute run replaces that model's `results/.../data` and `figures`; `--quick` overwrites a
  full run.
- New AI-written files go under `Claude files/`. Exceptions: `CLAUDE.md`, `TODO.md`,
  `.claude/`, `.githooks/`, and outputs code reads (`data/eba_moments.json`, `Empirics/outputs/`,
  `experiments/results|paper|figures/`, `diagnostics/regimes/`, `results/`).
- The modelling choices in `Claude files/docs/reference/` are deliberate. Read the entry before
  "fixing" one: several reversals are known to make sovereign risk expansionary.
- `global_projection/` style: one leading ALL-CAPS comment per module and function, lowercase
  `#` comments on specific lines, no docstrings; console output only in `reporting/prints.py`.
- SSJ: every Jacobian goes through `full_model.solve_jacobian_padded()`; new blocks only in
  `build_block_list()`; rebuild the regime cache before running `experiments/`.
- A commit that stages any `.py` updates STATE and PROGRESS in the same commit (hooks enforce
  it; once per clone run `git config core.hooksPath .githooks`). Update `TODO.md` when open
  work changes (every item has an owner of the open task, [joint] unless clearly one author's;
  a finished item is deleted there and recorded in PROGRESS with the date and PR or commit,
  without attribution) and `DATA.md` when a data number or its script changes.
- There is no handoff file. Session state goes into STATE (results), TODO (open work) and
  PROGRESS (what was done). Delete `docs/handoff/HANDOFF.md` if the handoff plugin recreates it.
- STATE holds current truth only: replace superseded sections and move the old text, verbatim,
  to `Claude files/docs/shelved/`.

## Economic reasoning standards
- Name the spread: the intermediation wedge (r^k - r at the steady state), the lending
  (working-capital) spread r^wc - r and the sovereign spread are different objects.
- The deterministic steady state is not the stochastic rest point. Global IRFs start at the
  rest point and are differenced against an unshocked path.
- μ = max{., 0} is kinked. Where μ is near zero at rest the fitted rules misstate it: report
  the exactly cleared response beside the fitted one, and quote backstop and welfare numbers only
  on exact clearing.
- Diagnostic expectation on a pure risk shock: the D bond price, both banks' net worth and D
  output and consumption should fall; the lending spread, D debt and taxes should rise. A
  response of the opposite sign must be investigated before it is reported, and classified as
  an accounting error, a numerical artefact (often at the kink) or a genuine
  general-equilibrium offset. In this project it has usually been an error.
- Mark-to-market losses come from duration interacting with the persistence of s, not from
  duration alone.
- No balance-sheet friction is a function of the default probability, and no spread equation
  takes the form `spread += parameter * π`.
- Never write "x% fundamental, y% non-fundamental"; describe the channels (expected-loss
  pricing, balance-sheet amplification).
- Percentage deviations divide by the variable's own steady-state level. Rates in bp/yr; p^d
  quarterly and annual; flows in level %.

## Running (only with a go-ahead)
```bash
python3 run.py global [--quick]     # global pipeline, then figures
python3 run.py ssj                  # SSJ, under the ssj interpreter
python3 run.py both --plots-only    # redraw figures from saved data, no solving
```
Everything runs from the repo root; imports are package-qualified.

## Where the detail lives (read on demand)
- `Claude files/docs/reference/`: `global_solver.md` (solver, grid, tests, tolerances),
  `global_modelling_choices.md` (deliberate choices, reasons, known limitations),
  `ssj_pipeline.md` (SSJ architecture, pricing rules, experiments, iteration checklist).
- `Claude files/docs/tpi_backstop_plan.md`: the backstop's design and solution method.
- `Claude files/docs/STATE.md` (Part I global, Part II SSJ): current results and calibration
  status only. History to 27 Oct 2026: `Claude files/docs/shelved/STATE_archive_2026-10-27.md`.
- `Claude files/docs/PROGRESS.md`: changelog. Branch roles: `Claude files/docs/GLOSSARY.md` section 9.
- `Claude files/docs/SPEC.md`: research goals and the paper's framing.
- `Claude files/REFACTOR_ARCHITECTURE.md`: full file map, environment and run guide.
- `Claude files/docs/shelved/CLAUDE_pre_refactor_2026-10-26.md`: the old file, verbatim.
