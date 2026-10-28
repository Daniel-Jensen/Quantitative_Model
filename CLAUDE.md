# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
It keeps the rules most tasks need; the detail is in the files under "Where the rest is".

## Project

Two-country heterogeneous-agent model of a monetary union with Gertler-Karadi
financial intermediaries and sovereign default risk, calibrated to the
2010–2012 Greek sovereign debt crisis. The default mechanism follows
**Bocola (2016, JPE) "The Pass-Through of Sovereign Risk"**: an EXOGENOUS
rise in the *priced* probability of default makes bond prices fall, banks take
mark-to-market losses, the occasionally-binding incentive constraint tightens,
lending spreads rise and output falls, with no default ever realized. Only D is
default-risky; F bonds are safe. Application: ECB asset purchases (TPI). Primary
output is a research paper (Overleaf: https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c).

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

## Claude files/ — where AI-generated supplementary files live (2026-10-01)

Documentation, notes, reports, logs, audit write-ups and tests written by Claude live under
`Claude files/`, at the path they would otherwise have had: `Claude files/docs/STATE.md`,
`Claude files/global_projection/tests/test_ss_identities.py`, `Claude files/diagnostics/cb_audit/`.
**Every new one goes there too** — never into the model packages, `docs/` or the repo root.
Exceptions, because tools or code read them at fixed paths: `CLAUDE.md`, `.claude/`,
`.githooks/`, and research outputs that code consumes (`data/eba_moments.json`,
`Empirics/outputs/`, `experiments/results|paper|figures/`, `diagnostics/regimes/`, `results/`).
Model code stays in its packages even when Claude wrote it. The path contains a space: quote
it in shell commands (`python3 "Claude files/global_projection/tests/test_state_grid.py"`).
The layout, the full file map and the run guide are in `Claude files/REFACTOR_ARCHITECTURE.md`.

## Environment

Two solvers now live in this repository, at non-overlapping paths, and they
need DIFFERENT interpreters. Pick the one that matches the code you are editing.

**Everything runs from the repo root**, which must be on `sys.path`: imports are
package-qualified (`from calibration.ssj import …`, `from linear_ssj.model.full_model
import …`, `from global_projection.blocks.bank import …`). The unified entry point
picks the right interpreter for you:

```bash
python3 run.py                    # MODEL and QUICK from run.py's configuration block
                                  #   ("global", full solve): what an editor "Run" click does
python3 run.py ssj                # SSJ pipeline, in a subprocess under SSJ_PYTHON
python3 run.py global [--quick]   # global pipeline, in this python3 process, step by step
python3 run.py both   [--quick]   # both, then the comparison
python3 run.py both --plots-only  # redraw every figure from the saved data, no solving
python3 run.py compare            # comparison only, from the saved data
```

Every run computes first and saves everything to `results/<SSJ|GLOBAL>/data/` (plain
`.json` + `.npz` via `results_io.py`; the global rules as `.pkl`), then draws every
figure FROM THE SAVED DATA into `results/<SSJ|GLOBAL>/figures/`, and logs the console to
`results/<SSJ|GLOBAL>/run.log`; the comparison goes to `results/COMPARISON/`. A compute
run REPLACES that model's `data/` and `figures/`, so a `--quick` preview overwrites a full
run. `results/` is git-ignored. `run.py` holds the global pipeline's settings (grid
resolution, shocks, which experiments run); the model parameters stay in `calibration/`.

**A. Sequence-space (SSJ) pipeline — `linear_ssj/`.** Always use
`/opt/anaconda3/envs/ssj/bin/python`. The base Anaconda environment has a broken
`liblapack` symlink that causes silent numerical failures.

```bash
conda activate ssj
/opt/anaconda3/envs/ssj/bin/python -m linear_ssj.main
```

**B. Global projection pipeline — `global_projection/`.** Plain `python3`
(numpy/scipy/matplotlib). It does NOT use `sequence_jacobian` and must not be
run under the `ssj` environment.

```bash
python3 run.py global              # the pipeline lives in run.py: compute_global / plot_global
```

## Rules for every task

- **Global runs:** the collocation Jacobian's workers re-import the model at every
  solve, so never edit model code while a run is in progress.
- **Global comment convention** (enforced across `global_projection/`): every module and
  every function carries exactly ONE leading ALL-CAPS comment saying what it is; any further
  explanation is lowercase `#` comments attached to the specific hard line. No docstrings, no
  bold markers, no prose blocks inside function bodies. Console output lives in `prints.py`,
  never inside the model blocks.
- **SSJ model definition:** `full_model.build_block_list()` is the single model definition
  (`full_model.py`, `linear_ssj/solve/tpi.py` and `diagnostics/regimes/regime_model.py` all
  call it). Add a block once, there.
- **SSJ Jacobians:** route every call through `full_model.solve_jacobian_padded()`, never
  `Block.solve_jacobian`. This must stay empty:
  `grep -rn "\.solve_jacobian(" --include="*.py" code experiments diagnostics | grep -v solve_jacobian_padded`
- **One calibration:** every script imports `get_calibration()`; the retired
  `audit_artifacts/` harness carried its own copy and silently tested a different model.
- **SSJ experiments:** rebuild the regime cache first
  (`/opt/anaconda3/envs/ssj/bin/python diagnostics/regimes/regime_model.py --force`);
  `calibration_override` patches the module attribute, so resolve `get_calibration` at use
  time; percentages divide by their own SS level (`common.pct_of_ss`).
- **Modelling choices are deliberate:** do not "fix" one without checking
  `Claude files/docs/SPEC.md` and the solver's guide. Global: single λ, priced vs realized
  default, endogenous debt in clearing, asymmetric country size, symmetric per-capita SS,
  occasionally-binding IC, Ω-kernel weights, multi-branch risk quadrature, the TPI backstop,
  predetermined deposit rate and capital, the union deposit market, the Hatchondo-Martinez
  perpetuity, the pure-haircut default branch, working capital, Walras redundancy.
  SSJ: Y = F(K_t), the nominal predetermined deposit rate, structural sovereign pricing with
  no wedges, Walras redundancy, the p-conversion in F-bank returns, nominal rigidities.
- **Paper wording:** never write "x% fundamental / y% non-fundamental"; the channels operate
  jointly, so describe mechanisms.

## Testing — fast checks first

SSJ (`linear_ssj/`):
```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest "Claude files/linear_ssj/tests/test_nkpc_blocks.py" -v  # 17 tests, ~1 s
/opt/anaconda3/envs/ssj/bin/python -m pytest "Claude files" -v  # 51 pass; 9 in test_eba_calibration need the raw EBA files
/opt/anaconda3/envs/ssj/bin/python -m linear_ssj.main           # the full regression test (or: python3 run.py ssj)
```
Check `goods_mkt_D/F` ≤ 1e−6, `ca_res_D` ≤ 1e−13, `deposit_mkt_D/F` ≤ 1e−13, the four
sticky-price targets exactly 0 at SS, every GK FOC leg `OK`, and that on the default shock
`n_inter_D[0]` and `Y_D[0]` both fall. Full thresholds: `guide_ssj.md`.

Global (`global_projection/`):
```bash
python3 "Claude files/global_projection/tests/test_ss_identities.py"      # fast
python3 "Claude files/global_projection/tests/test_fast_kernels.py"       # fast
python3 "Claude files/global_projection/tests/test_state_grid.py"         # fast
python3 "Claude files/global_projection/tests/test_no_unbound_names.py"   # fast
python3 "Claude files/global_projection/tests/test_collocation.py"        # the solver, ~90 s
python3 "Claude files/global_projection/tests/test_recursive_nesting.py"  # nesting gates N1-N4
```
Acceptance: `sum(F^2) <= m*(1e-9)^2` at every collocation stage; a risk-only shock must
lower Q_bD, n_D, n_F, Y_D[0] and C_D[0] and raise the lending spread (a positive Y or n
response is a bug). Full thresholds: `guide_global.md`.

## Committing

After any calibration or structural change, update `Claude files/docs/STATE.md`,
`PROGRESS.md` (one changelog entry per code commit) and `HANDOFF.md`, and commit them with
the `.py` files. Two hooks enforce this when a commit stages `*.py`:
`.claude/hooks/require-docs-before-commit.sh` and `.githooks/pre-commit` (enable once per
clone: `git config core.hooksPath .githooks`); bypass a false positive with
`git commit --no-verify`. `main` is the merge target; never merge `bank-cal` (port its
calibration values only) and do not reuse `audit` or `AB-audit`. Full branch list and
history: `HANDOFF.md`.

## Where the rest is

| File | Contains |
|------|----------|
| `Claude files/docs/guide_global.md` | Global projection: layout and file map, the collocation solver, grid and convergence record, acceptance thresholds, every modelling choice with its rationale (the TPI backstop included), known limitations |
| `Claude files/docs/guide_ssj.md` | SSJ: architecture, `solve_jacobian_padded`, policy experiments, modelling choices, nominal rigidities, open-issues table, the typical iteration, the full docs reference |
| `Claude files/docs/HANDOFF.md` | Quick start, session priorities, branch convention and history |
| `Claude files/docs/STATE.md` | Current calibration tables, residuals, open issues (both solvers) |
| `Claude files/docs/PROGRESS.md` | Dated changelog |
| `Claude files/docs/SPEC.md` | Research goals, modelling choices, calibration targets, paper framing |
| `Claude files/docs/cb_mechanism.md` | The canonical SSJ central-bank mechanism |
| `Claude files/REFACTOR_ARCHITECTURE.md` | Repository layout, full file map, run guide |
