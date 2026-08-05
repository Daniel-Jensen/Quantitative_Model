# Nominal Rigidities Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Rotemberg price Phillips curves and nominal deposit contracts to the two-country monetary-union HANK model, keeping the steady state bit-identical, and make the sticky model the paper's new baseline.

**Architecture:** Four new unknowns (`mc_D, pi_D, mc_F, pi_F`) and four new targets take the solver system from 23×23 to 27×27. Inflation is closed without any policy rate: the monetary-union terms-of-trade identity `p/p(−1) = (1+π_F)/(1+π_D)` pins the inflation differential off the existing `p` unknown, and a union-inflation normalisation pins the level. Markups are neutralised by subsidy so every new residual is exactly zero at the current steady state. Deposits then become nominal, which makes banks net nominal debtors and adds a Fisher debt-deflation channel.

**Tech Stack:** Python 3.12, `sequence-jacobian` (SSJ), numpy, pytest 9.0.3. Always use `/opt/anaconda3/envs/ssj/bin/python` — the base Anaconda environment has a broken `liblapack` symlink that causes silent numerical failures.

**Spec:** `docs/superpowers/specs/2026-08-05-nominal-rigidities-design.md`

**Branch:** `add-nkpc` (already created, already carries the spec commits)

---

## Background an engineer needs before starting

**This model has no test suite.** The regression test is running the whole pipeline
(`/opt/anaconda3/envs/ssj/bin/python code/main.py`, ~8-12 min) and reading printed
residuals. This plan adds fast unit tests for the new algebra so you are not
waiting 10 minutes to discover a sign error.

**How to unit-test an SSJ `@simple` block.** A `@simple`-decorated function becomes
a `SimpleBlock`. Call `blk.steady_state({...})` with a dict of every input; it
returns a dict containing inputs *and* outputs. Lags and leads (`x(-1)`, `x(+1)`)
evaluate to the same constant, so `steady_state` is really "evaluate this block at
constant values" — which is exactly what you want for algebraic identity tests,
including deliberately *off*-steady-state ones.

```python
ss = my_block.steady_state({'x': 2.0, 'y': 0.5})
assert ss['residual'] == 0.0
```

**SSJ gotcha:** SSJ discovers a block's outputs by running a regex over the
function's *source code* looking for the last `return` line. So every `@simple`
function must have a literal `return a, b` statement, and blocks cannot be defined
inside `python -c` strings (`inspect.getsource` fails). Define them in real files.

**Acceptance thresholds — use the MEASURED baseline below, not CLAUDE.md's list.**

CLAUDE.md quotes `goods_mkt_D ≤ 1e−14`. That refers to the *dynamic* solve, where
`goods_mkt_D` is an explicit entry in `targets_tp` and so is driven to Newton
tolerance — it is never printed on its own. What `main.py` *does* print under
`SS goods residuals:` are the **steady-state** residuals, which sit around 4e−7.
Do not mistake one for the other and report a false failure.

These are the actual values from a verified run of the pre-change model on
`add-nkpc` at commit `f3711bd` (Task 1, 2026-08-05). Every later task compares
against these:

| Printed line | Baseline value | Rule |
|---|---|---|
| `SS goods residuals: goods_mkt_D` | `-4.2493506589857954e-07` | must not degrade by an order of magnitude |
| `SS goods residuals: goods_mkt_F` | `-4.1914559989475464e-07` | same |
| `SS goods residuals: ca_res_D` | `6.852157730108388e-17` | same |
| Block residual table verdict | `All residuals < 1e-8  ✓` | must still print `✓` |
| `IC_D: θ − θ_tgt` | `1.776357e-15` | OK |
| `irfs_Z_D['b_gov_D'][499]` | `-0.001701` | stability |
| `irfs_def_D['b_gov_D'][499]` | `0.000014` | stability |
| `ρ_b (partial-eq.)` | `0.8451` | target < 0.95 |
| `n_inter_D[0]` on default shock | `-3.3804%` of SS | **must stay negative** |
| `Y_D[0]` on default shock | `-0.0149%` of SS | **must stay negative** |
| TPI `max|ca_res_D|` across γ | ≤ `7.55e-08` | ≤ 1e−7 |
| TPI `max|goods_mkt_F|` across γ | ≤ `2.44e-09` | ≤ 1e−7 |

**Sign checks:** on the default shock, `n_inter_D[0]` and `Y_D[0]` must both be
negative. Positive means a timing bug.

**Doc hooks are enforced — all three docs, every code commit.**
`.claude/hooks/require-docs-before-commit.sh` (verified) and its git-native twin
`.githooks/pre-commit` (active: `core.hooksPath` is set to `.githooks`) deny any
commit that stages `code/**` or any `*.py` unless **all three** of
`docs/STATE.md`, `docs/PROGRESS.md` and `docs/HANDOFF.md` are staged in the same
commit. Staging only `PROGRESS.md` is not enough.

So before each code commit, add a line to each of the three:

- `docs/PROGRESS.md` — a changelog bullet for this specific commit.
- `docs/STATE.md` — the current state after this commit (one line is fine for
  intermediate steps; Task 16 writes the real tables).
- `docs/HANDOFF.md` — where the work now stands, so an interrupted session can
  resume.

The `git add` lines in this plan already list all three. Do **not** use
`--no-verify`; the gate is deliberate project policy.

**Long-running commands must run in the background.** `code/main.py` takes about
8–12 minutes (several T=500 Jacobian solves at ~3 min each) and
`diagnostics/regimes/regime_model.py --force` takes longer. Both exceed the Bash
tool's 10-minute maximum timeout. Run them with `run_in_background: true` and poll
the output rather than blocking, or they will be killed mid-solve and you will
mistake a timeout for a model failure.

---

## One refinement to the spec, adopted here

The spec's substitution table replaces `rdep_i` with `rdep_exante_i` in
`intermediation_P1`, `divert_bond_foc` and `divert_portfolio_adj`. This plan
achieves the same model with a strictly smaller diff by **keeping the name
`rdep_i` for the ex-ante real rate**, which is already exactly what those three
blocks mean by it.

- `i_dep_i` — **new nominal unknown**, replaces `rdep_i` in the unknowns list.
- `rdep_i` — **now a derived block output**: the ex-ante real rate for t→t+1.
- `rdep_expost_i` — **new derived output**: the realised real rate at t on
  deposits placed at t−1.

Net effect: `intermediation_P1_{D,F}`, `divert_bond_foc_{D,F}`,
`divert_portfolio_adj`, `smart_steady_{D,F}` and `steady_auxilliary_{D,F}` are
**untouched**, and `steady_state.py`'s eight `ss['rdep_D']` / `ss['rdep_F']` reads
keep working unchanged. Only `deposit_return`, `bank_return` and `capital_fund`
change, plus one new block per country.

---

## File structure

| File | Change | Responsibility |
|---|---|---|
| `code/full_model.py` | Modify | Gains `build_block_list()` — the single definition of the model's block list. Updates `unknowns_tp` / `targets_tp` to 27×27. |
| `code/tpi.py` | Modify | Drops its duplicate block list; calls `build_block_list()` with TPI overrides. |
| `diagnostics/regimes/regime_model.py` | Modify | Same. |
| `code/equations_D.py` | Modify | Adds `price_nkpc_D`, `firm_profit_D`, `deposit_rates_D`. Modifies `labor_demand_D`, `income_D`, `deposit_return_D`, `bank_return_D`, `capital_fund_D`. |
| `code/equations_F.py` | Modify | Symmetric analogues. |
| `code/equations_global.py` | Modify | Adds `terms_of_trade`, `union_inflation`. |
| `code/calibration.py` | Modify | Adds `mu_p`, `kappa_p`, `omega_pi_D`, `pi`; retargets `mc`; renames `rdep` → `i_dep`. |
| `code/steady_state.py` | Modify | Adds new blocks to the SS block list. |
| `code/test_nkpc_blocks.py` | **Create** | Fast unit tests for all new block algebra. |
| `code/dump_irfs.py` | **Create** | Saves baseline IRFs to `.npz` so the equivalence gate is numerical, not eyeballed. |

---

## Phase 0 — Refactor (behaviour-preserving)

### Task 1: Extract a single `build_block_list()`

`full_model.py:69`, `tpi.py:145` and `diagnostics/regimes/regime_model.py:160`
each hardcode the `create_model` block list. Six blocks are about to be added to
all three. CLAUDE.md records that a drifting duplicate model is exactly what
invalidated the retired `audit_artifacts/` harness. Do this first, as a **pure
no-op**, so the next task validates one change rather than two.

The TPI list differs from the baseline list in four blocks: `budget_residual_D`,
`budget_residual_F`, `external_account_D` and `domestic_bond_clearing` are swapped
for `_tpi` variants defined in `tpi.py`. The factory therefore takes an
`overrides` mapping.

**Files:**
- Modify: `code/full_model.py:36-91`
- Modify: `code/tpi.py:145-164`
- Modify: `diagnostics/regimes/regime_model.py:160-183`

- [ ] **Step 1: Capture the baseline output**

```bash
cd /Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL
/opt/anaconda3/envs/ssj/bin/python code/main.py 2>&1 | tee /tmp/nkpc_baseline_main.log
```

Expected: completes, prints `Done — all figures saved to:`. Takes ~8-12 min.
Confirm the log contains a `b_gov_D` stability block and both sign lines.

- [ ] **Step 2: Add `build_block_list()` to `full_model.py`**

Insert immediately after the import block (after line 33, before `def build_and_solve`):

```python
def build_block_list(financial_solved_D, financial_solved_F, *,
                     hh_D=None, hh_F=None, overrides=None):
    """The single definition of the model's block list.

    Every consumer (full_model, tpi, diagnostics/regimes) calls this. A second
    copy of the list is how the retired audit_artifacts/ harness drifted into
    silently testing a different model — see CLAUDE.md.

    financial_solved_D/F : the runtime-constructed GK solved blocks.
    hh_D/hh_F            : optionally REPLACE the household blocks with versions
                           carrying extra hetoutputs (experiments/e4_distribution
                           adds per-decile consumption).
    overrides            : {name: block} used by the TPI layer to swap in its
                           _tpi variants without keeping a second list.
    """
    o = overrides or {}
    hh_D = hh_extended_D if hh_D is None else hh_D
    hh_F = hh_extended_F if hh_F is None else hh_F

    def pick(name, default):
        return o.get(name, default)

    return [
        # Country D
        deposit_return_D, tax_rule_D, hh_D, ghh_composite_D,
        sdf_D, sdf_banker_D, government_default_D, financial_solved_D,
        bond_return_D, bank_return_D, capital_fund_D, cap_adj_cost_inter_D, macro_pru_tax_D,
        intermediation_P2_D, intermediation_P3_D, k_balance_sheet_D,
        capital_adj_D, capital_producer_profit_D,
        pick('budget_residual_D', budget_residual_D),
        labor_D, labor_market_D, labor_demand_D, banker_div_res_D,
        market_clearing_D, welfare_agg_D,
        # Country F
        deposit_return_F, tax_rule_F, hh_F, ghh_composite_F,
        sdf_F, sdf_banker_F, government_default_F, financial_solved_F,
        bond_return_F, bank_return_F, capital_fund_F, cap_adj_cost_inter_F, macro_pru_tax_F,
        intermediation_P2_F, intermediation_P3_F, k_balance_sheet_F,
        capital_adj_F, capital_producer_profit_F,
        pick('budget_residual_F', budget_residual_F),
        labor_F, labor_market_F, labor_demand_F, banker_div_res_F,
        market_clearing_F, welfare_agg_F,
        # Global
        ces_price_D, import_demand_D, ces_price_F, import_demand_F,
        trade_balance,
        pick('external_account_D', external_account_D),
        pick('domestic_bond_clearing', domestic_bond_clearing),
        bond_yield, portfolio_level_anchors, divert_portfolio_adj,
        divert_bond_foc_D, divert_bond_foc_F, global_goods_mkt,
    ]
```

Move the two `hh_extended` imports from the bottom of the file (lines 176-178) up
into the main import block so `build_block_list` can see them, keeping the
`# noqa: F401` re-export comment for `tpi.py`.

- [ ] **Step 3: Use the factory in `build_and_solve`**

Replace the `ha_full = sj.create_model([...])` call (`full_model.py:69-91`) with:

```python
    ha_full = sj.create_model(
        build_block_list(financial_solved_D, financial_solved_F),
        name="Full 2-Country MU HANK — GHH Preferences, Flex Price & Wage, No CB",
    )
```

- [ ] **Step 4: Add a TPI overrides helper to `tpi.py`**

Add near the top of `tpi.py`, after the `_tpi` block definitions:

```python
def tpi_overrides():
    """The four blocks the TPI layer swaps into the shared block list."""
    return {
        'budget_residual_D':     budget_residual_D_tpi,
        'budget_residual_F':     budget_residual_F_tpi,
        'external_account_D':    external_account_D_tpi,
        'domestic_bond_clearing': domestic_bond_clearing_tpi,
    }
```

Then replace `tpi.py`'s `sj.create_model([...])` (lines 145-164) with:

```python
    from full_model import build_block_list
    ha_full_tpi = sj.create_model(
        build_block_list(financial_solved_D, financial_solved_F,
                         overrides=tpi_overrides()),
        name="Full 2-Country MU HANK — TPI Extension",
    )
```

- [ ] **Step 5: Point `regime_model.py` at the factory**

Replace its `return sj.create_model([...])` (lines 160-183) with:

```python
    from full_model import build_block_list
    return sj.create_model(
        build_block_list(financial_solved_D, financial_solved_F,
                         hh_D=hh_D, hh_F=hh_F,
                         overrides=t.tpi_overrides()),
        name="Full 2-Country MU HANK — TPI Extension (regimes cache, main)",
    )
```

Delete the now-dead `hh_D = t.hh_extended_D if hh_D is None else hh_D` lines
immediately above it — `build_block_list` handles the `None` default itself.

- [ ] **Step 6: Verify bit-identical**

```bash
/opt/anaconda3/envs/ssj/bin/python code/main.py 2>&1 | tee /tmp/nkpc_refactor_main.log
diff <(grep -v 'Output directory' /tmp/nkpc_baseline_main.log) \
     <(grep -v 'Output directory' /tmp/nkpc_refactor_main.log)
```

Expected: **no output from `diff`**. Any difference means the refactor changed the
model — stop and find it. Do not proceed with a non-empty diff.

- [ ] **Step 7: Commit**

Add to `docs/PROGRESS.md` under a new dated heading:
`- Extracted build_block_list() in full_model.py; tpi.py and regime_model.py now share it. Verified bit-identical main.py output.`

```bash
git add code/full_model.py code/tpi.py diagnostics/regimes/regime_model.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "refactor: single build_block_list() shared by full_model, tpi, regimes

No-op. Verified main.py output is byte-identical before and after."
```

---

## SSJ library defect found in Task 9 — read before touching any Jacobian solve

**SSJ 1.0.0 cannot solve this system with stock `Block.solve_jacobian`.**

`CombinedBlock._jacobian` (`blocks/combined_block.py:104-119`) seeds
`total_Js = JacobianDict.identity(inputs)` from the *shock* list, visits a block
only `if (inputs & block.inputs) and (outputs & block.outputs)`, and returns
`total_Js[original_outputs & total_Js.outputs, :]`. A target reachable from no
shock is therefore silently dropped from H_Z. `Block.solve_jacobian`
(`blocks/block.py:260`) then calls `np.linalg.solve(H_U, H_Z.pack(T))` with
mismatched shapes:

```
ValueError: solve: Input operand 1 has a mismatch in its core dimension 0,
with gufunc signature (m,m),(m,n)->(m,n) (size 11500 is different from 13500)
```

`11500 = 23*500`, `13500 = 27*500`. All four new targets are pure functions of
the solver's own unknowns — `nkpc_p_res_D/F` (pi, mc), `tot_res` (p, pi_D, pi_F),
`union_pi_res` (pi_D, pi_F) — and contain no `Z_*` or `shock_def_*` symbol.

**Fix: `full_model.solve_jacobian_padded()`**, added in Task 9. It restores the
missing rows as zeros, which is **exact, not an approximation**: `dH/dZ` at fixed
unknowns is identically zero when the shock never appears in the equation. It
mirrors `Block.solve_jacobian` line-for-line and prints the padded row names on
every solve, so it cannot go silent.

**A 25x25 rewrite does NOT avoid this — do not attempt it.** Solving `tot_res`
and `union_pi_res` analytically for `pi_D`, `pi_F` (which has an exact closed
form, `pi_D = (1-g)/(g + omega/(1-omega))` with `g = p/p(-1)`) would make them
block outputs and drop the system to 25x25. But `nkpc_p_res_D` would still depend
only on `pi_D` — now a function of `p`, still an unknown — and `mc_D`, another
unknown. The gate at `combined_block.py:115` tests against the *shock* set, so
the two NKPC targets would still be dropped, giving 23 H_Z rows against a 25x25
H_U. Same defect, smaller numbers, at the cost of rewriting committed work.

**Every Jacobian call site must use the padded helper.** Task 9b converts them.

## Phase 1 — Price rigidity (deposits still real)

### Task 2: Markup rent — `firm_profit_{D,F}`

With the markup in labour demand only, `w*N = mu_p*mc*(1-alpha)*Y` while the
capital return is untouched, so off steady state factor payments no longer exhaust
output and the residual has nowhere to go. That is a Walras leak of the W-1 / W-2
class. This block routes it.

The rent is distributed to households **in proportion to productivity `e`**
(Auclert–Rognlie–Straub), *not* lump-sum. Markups are countercyclical, so a
lump-sum rebate would hand households rising income exactly when output falls.
Distributing on `e` makes household labour-plus-profit income `(1-alpha)*Y*e` —
identical to the flexible model — so the wedge affects the firm's hiring decision
only, and because the share depends on the household's *type* rather than its
hours, the marginal wage stays `w` and `labor_market_{D,F}` is unchanged.

**Files:**
- Create: `code/test_nkpc_blocks.py`
- Modify: `code/equations_D.py` (append after `labor_demand_D`, line 295)
- Modify: `code/equations_F.py` (append after `labor_demand_F`, line 249)

- [ ] **Step 1: Write the failing tests**

Create `code/test_nkpc_blocks.py`:

```python
"""Fast algebraic tests for the nominal-rigidity blocks.

These evaluate SSJ @simple blocks directly via .steady_state(), which is just
"evaluate at constant values" -- lags and leads collapse to the same constant.
That makes it usable for deliberately OFF-steady-state identity checks too.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


# ── Markup rent ───────────────────────────────────────────────────────────────

def test_firm_profit_is_zero_at_steady_state():
    from equations_D import firm_profit_D
    mu_p = 1.20
    ss = firm_profit_D.steady_state({
        'Y_D': 1.0, 'N_D': 0.8, 'alpha_D': 0.33,
        'mu_p_D': mu_p, 'mc_D': 1.0 / mu_p,
    })
    assert ss['profit_D'] == pytest.approx(0.0, abs=1e-15)


def test_firm_profit_restores_factor_exhaustion_off_steady_state():
    """w*N + profit must equal (1-alpha)*Y for ANY mc, so that adding the
    capital share alpha*Y exhausts output exactly."""
    from equations_D import firm_profit_D, labor_demand_D
    mu_p, mc, Y, N, alpha = 1.20, 0.79, 1.03, 0.81, 0.33

    # w from labour demand at this mc (w_res_D == 0 defines w)
    w = mu_p * mc * (1 - alpha) * Y / N

    ss = firm_profit_D.steady_state({
        'Y_D': Y, 'N_D': N, 'alpha_D': alpha, 'mu_p_D': mu_p, 'mc_D': mc,
    })
    assert w * N + ss['profit_D'] == pytest.approx((1 - alpha) * Y, rel=1e-14)

    # and the wage we assumed really is the one labor_demand_D implies
    ld = labor_demand_D.steady_state({
        'w_D': w, 'Y_D': Y, 'N_D': N, 'alpha_D': alpha,
        'mu_p_D': mu_p, 'mc_D': mc,
    })
    assert ld['w_res_D'] == pytest.approx(0.0, abs=1e-14)


def test_firm_profit_F_matches_D():
    from equations_D import firm_profit_D
    from equations_F import firm_profit_F
    args = dict(Y=1.03, N=0.81, alpha=0.33, mu_p=1.20, mc=0.79)
    d = firm_profit_D.steady_state({
        'Y_D': args['Y'], 'N_D': args['N'], 'alpha_D': args['alpha'],
        'mu_p_D': args['mu_p'], 'mc_D': args['mc'],
    })
    f = firm_profit_F.steady_state({
        'Y_F': args['Y'], 'N_F': args['N'], 'alpha_F': args['alpha'],
        'mu_p_F': args['mu_p'], 'mc_F': args['mc'],
    })
    assert d['profit_D'] == pytest.approx(f['profit_F'], rel=1e-15)
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v
```

Expected: FAIL — `ImportError: cannot import name 'firm_profit_D'`.

- [ ] **Step 3: Implement `firm_profit_D`**

Append to `code/equations_D.py` after `labor_demand_D` (line 295):

```python
@simple
def firm_profit_D(Y_D, alpha_D, mu_p_D, mc_D):
    # Markup rent. With sticky prices labour demand pays mu_p*mc*(1-alpha)*Y
    # while capital still earns alpha*Y (capital_adj_D is unchanged), so off SS
    # factor payments do not exhaust output. profit_D is that residual; leaving
    # it unrouted is a Walras leak of the W-1/W-2 class.
    #
    # Distributed to households in proportion to productivity e (Auclert-Rognlie-
    # Straub), NOT lump-sum: markups are countercyclical, so a lump-sum rebate
    # would raise household income exactly when output falls. On the e rule,
    # w*N*e + profit*e = (1-alpha)*Y*e -- identical to the flex model -- so the
    # wedge affects the firm's hiring decision only, and because the share
    # depends on type rather than hours the marginal wage is still w_D and
    # labor_market_D is unchanged.
    #
    # Zero at SS, where mu_p*mc = 1.
    profit_D = (1.0 - mu_p_D * mc_D) * (1.0 - alpha_D) * Y_D
    return profit_D
```

Append the symmetric block to `code/equations_F.py` after `labor_demand_F` (line 249):

```python
@simple
def firm_profit_F(Y_F, alpha_F, mu_p_F, mc_F):
    # See firm_profit_D.
    profit_F = (1.0 - mu_p_F * mc_F) * (1.0 - alpha_F) * Y_F
    return profit_F
```

`N` is deliberately *not* in the signature — it does not appear in the expression,
and an unused input would make SSJ record a spurious DAG edge. The tests above
pass `N_D` / `N_F` in their dicts anyway; SSJ's `steady_state()` ignores extra
keys (verified), so they need no change.

- [ ] **Step 4: Run the tests to verify they pass**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v
```

Expected: `test_firm_profit_is_zero_at_steady_state` and
`test_firm_profit_F_matches_D` PASS. The factor-exhaustion test still FAILS on
`labor_demand_D` not accepting `mu_p_D` / `mc_D` — that is Task 4.

- [ ] **Step 5: Commit**

```bash
git add code/equations_D.py code/equations_F.py code/test_nkpc_blocks.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: firm_profit_D/F routes the markup rent proportional to e"
```

---

### Task 3: Price Phillips curves — `price_nkpc_{D,F}`

Rotemberg form in producer-price inflation. The gap is written as a **ratio**, so
it is unit-free and linearises to exactly `mc_hat`; published Calvo slopes are
therefore directly usable for `kappa_p` with no steady-state rescaling.

**Files:**
- Modify: `code/equations_D.py` (append after `firm_profit_D`)
- Modify: `code/equations_F.py` (append after `firm_profit_F`)
- Modify: `code/test_nkpc_blocks.py`

- [ ] **Step 1: Write the failing tests**

Append to `code/test_nkpc_blocks.py`:

```python
# ── Price Phillips curve ──────────────────────────────────────────────────────

def test_price_nkpc_is_zero_at_steady_state():
    from equations_D import price_nkpc_D
    mu_p = 1.20
    ss = price_nkpc_D.steady_state({
        'pi_D': 0.0, 'mc_D': 1.0 / mu_p, 'mu_p_D': mu_p,
        'kappa_p_D': 0.0871, 'beta_D': 0.985,
    })
    assert ss['nkpc_p_res_D'] == pytest.approx(0.0, abs=1e-15)


def test_price_nkpc_flex_limit_forces_mc_to_one_over_mu_p():
    """As kappa_p -> inf the residual/kappa_p -> -(mu_p*mc - 1), so setting the
    residual to zero drives mu_p*mc -> 1, which is the competitive condition."""
    from equations_D import price_nkpc_D
    mu_p = 1.20
    base = {'pi_D': 0.0, 'mu_p_D': mu_p, 'beta_D': 0.985}
    off_mc = 0.79                       # != 1/mu_p = 0.8333...
    for kappa in (1e2, 1e4, 1e6):
        ss = price_nkpc_D.steady_state({**base, 'mc_D': off_mc, 'kappa_p_D': kappa})
        implied_gap = -ss['nkpc_p_res_D'] / kappa
        assert implied_gap == pytest.approx(mu_p * off_mc - 1.0, rel=1e-12)


def test_price_nkpc_gap_linearises_to_mc_hat():
    """d(mu_p*mc - 1)/d(mc/mc_ss) evaluated at mc_ss = 1/mu_p equals 1 for ANY
    mu_p -- which is why mu_p is a free normalisation to first order."""
    from equations_D import price_nkpc_D
    for mu_p in (1.05, 1.20, 1.50):
        mc_ss = 1.0 / mu_p
        h = 1e-7
        base = {'pi_D': 0.0, 'mu_p_D': mu_p, 'kappa_p_D': 1.0, 'beta_D': 0.985}
        up = price_nkpc_D.steady_state({**base, 'mc_D': mc_ss * (1 + h)})
        dn = price_nkpc_D.steady_state({**base, 'mc_D': mc_ss * (1 - h)})
        # residual = -kappa*(gap), kappa = 1 -> d(gap)/d(mc_hat) = -d(res)/d(mc_hat)
        d_gap = -(up['nkpc_p_res_D'] - dn['nkpc_p_res_D']) / (2 * h)
        assert d_gap == pytest.approx(1.0, rel=1e-6)


def test_price_nkpc_F_matches_D():
    from equations_D import price_nkpc_D
    from equations_F import price_nkpc_F
    args = dict(pi=0.001, mc=0.79, mu_p=1.20, kappa=0.0871, beta=0.985)
    d = price_nkpc_D.steady_state({
        'pi_D': args['pi'], 'mc_D': args['mc'], 'mu_p_D': args['mu_p'],
        'kappa_p_D': args['kappa'], 'beta_D': args['beta'],
    })
    f = price_nkpc_F.steady_state({
        'pi_F': args['pi'], 'mc_F': args['mc'], 'mu_p_F': args['mu_p'],
        'kappa_p_F': args['kappa'], 'beta_F': args['beta'],
    })
    assert d['nkpc_p_res_D'] == pytest.approx(f['nkpc_p_res_F'], rel=1e-15)
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -k price_nkpc -v
```

Expected: FAIL — `ImportError: cannot import name 'price_nkpc_D'`.

- [ ] **Step 3: Implement the blocks**

Append to `code/equations_D.py`:

```python
@simple
def price_nkpc_D(pi_D, mc_D, mu_p_D, kappa_p_D, beta_D):
    # Rotemberg NK Phillips curve in D producer-price inflation.
    #
    # The gap is a RATIO (mu_p*mc - 1), so it is unit-free and linearises to
    # exactly mc_hat for any mu_p -- published Calvo slopes are directly usable
    # for kappa_p with no SS rescaling, and mu_p is a free normalisation to
    # first order under the subsidy neutralisation.
    #
    # Subsidy-neutralised: mc_ss = 1/mu_p, so the gap and pi are both exactly
    # zero at the current SS and the SS is bit-identical to the flex model.
    # kappa_p -> inf recovers flexible prices (mu_p*mc = 1).
    #
    # Discounted at constant beta rather than SDF_D: since pi_ss = 0 the SDF
    # deviation multiplies a zero, so the two are identical to first order and
    # the model is solved by linearised solve_jacobian.
    nkpc_p_res_D = pi_D - beta_D * pi_D(+1) - kappa_p_D * (mu_p_D * mc_D - 1.0)
    return nkpc_p_res_D
```

Append to `code/equations_F.py`:

```python
@simple
def price_nkpc_F(pi_F, mc_F, mu_p_F, kappa_p_F, beta_F):
    # See price_nkpc_D.
    nkpc_p_res_F = pi_F - beta_F * pi_F(+1) - kappa_p_F * (mu_p_F * mc_F - 1.0)
    return nkpc_p_res_F
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -k price_nkpc -v
```

Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add code/equations_D.py code/equations_F.py code/test_nkpc_blocks.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: price_nkpc_D/F Rotemberg Phillips curves in PPI inflation"
```

---

### Task 4: Markup wedge in labour demand

**Files:**
- Modify: `code/equations_D.py:292-295`
- Modify: `code/equations_F.py:246-249`

- [ ] **Step 1: Confirm the already-written test fails**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py::test_firm_profit_restores_factor_exhaustion_off_steady_state -v
```

Expected: FAIL — `labor_demand_D` does not accept `mu_p_D` / `mc_D`.

- [ ] **Step 2: Add the flex-limit test**

Append to `code/test_nkpc_blocks.py`:

```python
def test_labor_demand_collapses_to_competitive_at_ss_markup():
    """At mc = 1/mu_p the condition must be exactly w = (1-alpha)Y/N, which is
    what makes the steady state bit-identical to the flex model."""
    from equations_D import labor_demand_D
    mu_p, Y, N, alpha = 1.20, 1.03, 0.81, 0.33
    w_competitive = (1 - alpha) * Y / N
    ss = labor_demand_D.steady_state({
        'w_D': w_competitive, 'Y_D': Y, 'N_D': N, 'alpha_D': alpha,
        'mu_p_D': mu_p, 'mc_D': 1.0 / mu_p,
    })
    assert ss['w_res_D'] == pytest.approx(0.0, abs=1e-15)
```

- [ ] **Step 3: Modify `labor_demand_D`**

Replace `code/equations_D.py:292-295` in full:

```python
@simple
def labor_demand_D(w_D, Y_D, N_D, alpha_D, mu_p_D, mc_D):
    # Firm FOC with a price markup: w = mu_p*mc*(1-alpha)*Y/N.
    # The mu_p factor IS the production subsidy tau_s = 1 - 1/mu_p: at the SS
    # markup mc = 1/mu_p this collapses to the competitive w = (1-alpha)Y/N
    # identically, so the steady state is unchanged. Off SS the wedge shifts
    # labour demand, which is what makes N -- and hence output -- respond to
    # demand rather than being pinned by Z, K and P_CES alone.
    # The rent (1 - mu_p*mc)(1-alpha)Y is routed by firm_profit_D.
    w_res_D = w_D - mu_p_D * mc_D * (1 - alpha_D) * Y_D / N_D
    return w_res_D
```

Replace `code/equations_F.py:246-249` in full:

```python
@simple
def labor_demand_F(w_F, Y_F, N_F, alpha_F, mu_p_F, mc_F):
    # See labor_demand_D. Pins the wage in ha_full (drop labor_mkt_res_F there).
    w_res_F = w_F - mu_p_F * mc_F * (1 - alpha_F) * Y_F / N_F
    return w_res_F
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v
```

Expected: all tests pass, including
`test_firm_profit_restores_factor_exhaustion_off_steady_state`.

- [ ] **Step 5: Commit**

```bash
git add code/equations_D.py code/equations_F.py code/test_nkpc_blocks.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: markup wedge mu_p*mc in labor_demand_D/F"
```

---

### Task 5: Global closure — `terms_of_trade` and `union_inflation`

In a monetary union the nominal exchange rate is fixed at 1, so the terms of trade
**is** the accumulated inflation differential. That identity pins the differential
off `p`, which is already an unknown. The union-inflation normalisation pins the
level — the `phi_pi -> inf` limit of a Taylor rule, stated as an abstraction
rather than a modelled rule.

**Files:**
- Modify: `code/equations_global.py` (append after `bond_yield`, line 44)
- Modify: `code/test_nkpc_blocks.py`

- [ ] **Step 1: Write the failing tests**

Append to `code/test_nkpc_blocks.py`:

```python
# ── Global closure ────────────────────────────────────────────────────────────

def test_global_residuals_zero_at_steady_state():
    from equations_global import terms_of_trade, union_inflation
    tot = terms_of_trade.steady_state({'p': 0.99, 'pi_D': 0.0, 'pi_F': 0.0})
    assert tot['tot_res'] == pytest.approx(0.0, abs=1e-15)
    uni = union_inflation.steady_state({'pi_D': 0.0, 'pi_F': 0.0, 'omega_pi_D': 0.071})
    assert uni['union_pi_res'] == pytest.approx(0.0, abs=1e-15)


def test_closure_puts_93pct_of_tot_move_into_D_deflation():
    """Solving tot_res = 0 and union_pi_res = 0 together gives
    pi_D = -(1 - omega)*dlog p and pi_F = omega*dlog p. At the capital-key
    omega = 0.071 that is a 93/7 split -- the internal-devaluation pattern.
    Verified here by residual evaluation, not by re-deriving the algebra."""
    import math
    from equations_global import terms_of_trade, union_inflation
    omega = 0.071
    dlog_p = 1e-4                      # small so the log-linear form is accurate
    pi_D = -(1 - omega) * dlog_p
    pi_F = omega * dlog_p

    uni = union_inflation.steady_state({'pi_D': pi_D, 'pi_F': pi_F,
                                        'omega_pi_D': omega})
    assert uni['union_pi_res'] == pytest.approx(0.0, abs=1e-18)

    # tot_res compares p/p(-1) against (1+pi_F)/(1+pi_D); steady_state() sets
    # p(-1) = p, so feed the implied gross growth rate directly instead.
    implied = (1 + pi_F) / (1 + pi_D)
    assert math.log(implied) == pytest.approx(dlog_p, rel=1e-6)


def test_omega_one_half_splits_evenly():
    """Guards the calibration argument: at omega = 0.5 the adjustment splits
    50/50, which is counterfactual for GR/DE. See the spec."""
    omega = 0.5
    dlog_p = 1e-4
    assert -(1 - omega) * dlog_p == pytest.approx(-0.5 * dlog_p, rel=1e-15)
    assert omega * dlog_p == pytest.approx(0.5 * dlog_p, rel=1e-15)
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -k "global_residuals or closure or omega_one" -v
```

Expected: FAIL — `ImportError: cannot import name 'terms_of_trade'`.

- [ ] **Step 3: Implement the blocks**

Append to `code/equations_global.py`:

```python
@simple
def terms_of_trade(p, pi_D, pi_F):
    # p = P_F/P_D in euro producer prices. In a monetary union the nominal
    # exchange rate is fixed at 1, so terms-of-trade movement IS the inflation
    # differential. This pins pi_D - pi_F off an unknown that already exists.
    # Zero at SS: p/p(-1) = 1 and pi_D = pi_F = 0.
    tot_res = p / p(-1) - (1.0 + pi_F) / (1.0 + pi_D)
    return tot_res


@simple
def union_inflation(pi_D, pi_F, omega_pi_D):
    # The ECB stabilises union-wide producer-price inflation -- the phi_pi -> inf
    # limit of a Taylor rule, stated as an abstraction and NOT a modelled rule.
    # Financial contracts carry no policy rate, so no Fisher relation is needed
    # to close the nominal side.
    #
    # With terms_of_trade this gives pi_D = -(1 - omega_pi_D)*dlog p. At the
    # capital-key omega_pi_D = 0.071, 93% of any terms-of-trade adjustment is D
    # producer-price deflation and 7% is F inflation -- the 2010-12 internal-
    # devaluation pattern. Do NOT use model GDP weights: the model normalises
    # Y_D_ss ~ Y_F_ss ~ 1, so those would give ~0.5 and split it evenly.
    pi_U         = omega_pi_D * pi_D + (1.0 - omega_pi_D) * pi_F
    union_pi_res = pi_U
    return pi_U, union_pi_res
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add code/equations_global.py code/test_nkpc_blocks.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: terms_of_trade and union_inflation close the nominal side"
```

---

### Task 6: Route the rent into household income

**Files:**
- Modify: `code/equations_D.py:64-69`
- Modify: `code/equations_F.py:60-66`

- [ ] **Step 1: Modify `income_D`**

Replace `code/equations_D.py:64-69` in full:

```python
def income_D(e_grid_D, w_D, N_D, div_D, div_fund_D, profit_D, tau_D, lamb_D, P_CES_D, T_ls_D):
    # div_fund_D: rebate from the passive capital fund (zero when omega_K_D=1).
    # profit_D: markup rent, distributed in proportion to productivity e (see
    # firm_profit_D). w_D*N_D*e + profit_D*e = (1-alpha)*Y_D*e exactly, so
    # household income is identical to the flex model and the markup wedge acts
    # only on the firm's hiring decision. Zero at SS.
    y_pre_D  = (w_D * N_D * e_grid_D + profit_D * e_grid_D + div_D + div_fund_D) / P_CES_D
    z_D      = lamb_D * (y_pre_D ** (1 - tau_D)) - T_ls_D
    t_paid_D = y_pre_D - z_D
    return z_D, t_paid_D
```

Replace `code/equations_F.py:60-66` in full:

```python
def income_F(e_grid_F, w_F, N_F, div_F, div_fund_F, profit_F, tau_F, lamb_F, P_CES_F, T_ls_F):
    # See income_D. profit_F is the markup rent, distributed on e; zero at SS.
    y_pre_F  = (w_F * N_F * e_grid_F + profit_F * e_grid_F + div_F + div_fund_F) / P_CES_F
    z_F      = lamb_F * (y_pre_F ** (1 - tau_F)) - T_ls_F
    t_paid_F = y_pre_F - z_F
    return z_F, t_paid_F
```

- [ ] **Step 2: Verify the hetinput picks up the new argument**

```bash
/opt/anaconda3/envs/ssj/bin/python -c "
import sys; sys.path.insert(0, 'code')
from equations_D import hh_extended_D
assert 'profit_D' in hh_extended_D.inputs, hh_extended_D.inputs
print('profit_D wired into hh_extended_D')
"
```

Expected: `profit_D wired into hh_extended_D`.

- [ ] **Step 3: Commit**

```bash
git add code/equations_D.py code/equations_F.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: markup rent enters household income proportional to e"
```

---

### Task 7: Calibration parameters

`mc_D` and `mc_F` already exist in `calibration.py:284` at 1.0, but nothing reads
them — they are dead entries. They now become live and must be retargeted to
`1/mu_p`.

**Files:**
- Modify: `code/calibration.py:281-284`

- [ ] **Step 1: Replace the wage-markup / SS-real-variables block**

Replace `code/calibration.py:281-284` in full:

```python
        # ── Wage Markups ──────────────────────────────────────────────────────
        # Unchanged: wages are flexible. mu_w = 1 is the SS-neutralising device
        # in labor_ss_D/F; there is no wage Phillips curve.
        'mu_w_D':       1.0,     'mu_w_F':       1.0,

        # ── Price Rigidity (Rotemberg) ────────────────────────────────────────
        # mu_p: gross price markup, epsilon_p = 6. FREE TO FIRST ORDER under the
        #   subsidy neutralisation -- the gap (mu_p*mc - 1) linearises to mc_hat
        #   for any mu_p -- so this needs no defending unless live markups are
        #   ever adopted.
        # mc: SS real marginal cost = 1/mu_p. The production subsidy
        #   tau_s = 1 - 1/mu_p makes labour demand collapse to the competitive
        #   w = (1-alpha)Y/N at this value, so the SS is bit-identical to flex.
        # kappa_p: Calvo theta_p = 0.75 at beta = 0.985, slope
        #   (1-theta)(1-beta*theta)/theta = 0.0871. Euro-area IPN median price
        #   duration ~4 quarters (Alvarez et al. 2006; Dhyne et al. 2006).
        #   Agrees with Bi-Foerster-Traum's implied 0.0846 to within 3%.
        # pi: SS producer-price inflation, exactly zero.
        'mu_p_D':       1.20,    'mu_p_F':       1.20,
        'mc_D':    1.0 / 1.20,   'mc_F':    1.0 / 1.20,
        'kappa_p_D':    0.0871,  'kappa_p_F':    0.0871,
        'pi_D':         0.0,     'pi_F':         0.0,

        # omega_pi_D: weight on D in the union producer-price aggregate that the
        # ECB is assumed to stabilise. = 1 - kappa_cb_F, the renormalised
        # two-country capital key (BuBa 26.1 / BoG 2.0 of the euro-area key).
        # DO NOT use model GDP weights: the model normalises Y_D_ss ~ Y_F_ss ~ 1,
        # so they would give ~0.5 and split the terms-of-trade adjustment evenly
        # between Greek deflation and German inflation -- counterfactual for
        # 2010-12. Load-bearing twice over once deposits are nominal, since it
        # scales pi_D and hence the Fisher revaluation on bank balance sheets.
        'omega_pi_D':   0.071,
```

- [ ] **Step 2: Verify the calibration loads and the SS markup is consistent**

```bash
/opt/anaconda3/envs/ssj/bin/python -c "
import sys; sys.path.insert(0, 'code')
from calibration import get_calibration
c = get_calibration()
for k in ('mu_p_D','mu_p_F','mc_D','mc_F','kappa_p_D','kappa_p_F','pi_D','pi_F','omega_pi_D'):
    print(f'  {k} = {c[k]}')
assert abs(c['mu_p_D']*c['mc_D'] - 1.0) < 1e-15, 'subsidy neutralisation broken'
assert abs(c['mu_p_F']*c['mc_F'] - 1.0) < 1e-15, 'subsidy neutralisation broken'
print('mu_p*mc == 1 in both countries')
"
```

Expected: all nine printed, then `mu_p*mc == 1 in both countries`.

- [ ] **Step 3: Commit**

```bash
git add code/calibration.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: mu_p, kappa_p, omega_pi_D calibration; mc retargeted to 1/mu_p"
```

---

### Task 8: Seed the new blocks into the steady state

The SS solve must carry `mc`, `pi`, `profit` and the two global residuals so
`ss_final` hands them to `solve_jacobian`. All are exactly zero (or exactly
`1/mu_p`) at the current SS, so **the solved steady state must not move**.

**Files:**
- Modify: `code/steady_state.py:1-40` (imports), `code/steady_state.py:151-162` (block list)

- [ ] **Step 1: Add the new blocks to the SS imports**

In `code/steady_state.py`, add to the `from equations_D import ...` list:
`firm_profit_D, price_nkpc_D`. Add to the `from equations_F import ...` list:
`firm_profit_F, price_nkpc_F`. Add to the `from equations_global import ...` list:
`terms_of_trade, union_inflation`.

- [ ] **Step 2: Add them to the SS `create_model` list**

In `solve_steady_state` (`code/steady_state.py:151`), extend the list. `labor_ss_D`
already sits on the line with `banker_div_D` and `government_ss_D`; add the new
blocks alongside:

```python
        hh_extended_D, smart_steady_D, market_clearing_D, steady_auxilliary_D,
        banker_div_D, government_ss_D, labor_ss_D, firm_profit_D, price_nkpc_D,
        hh_extended_F, smart_steady_F, market_clearing_F, steady_auxilliary_F,
        banker_div_F, government_ss_F, labor_ss_F, firm_profit_F, price_nkpc_F,
        ces_price_D, import_demand_D, ces_price_F, import_demand_F,
        deposit_return_D, deposit_return_F,
        bond_yield,
        trade_balance, external_account_D, global_goods_mkt,
        terms_of_trade, union_inflation,
```

`labor_demand_D/F` is deliberately **not** in the SS list (the SS uses
`labor_ss_D/F`), so the markup change there cannot touch the SS solve.

- [ ] **Step 3: Verify the SS is unchanged and the new residuals are zero**

```bash
/opt/anaconda3/envs/ssj/bin/python -c "
import sys; sys.path.insert(0, 'code')
from calibration import get_calibration
from steady_state import solve_steady_state
r = solve_steady_state(get_calibration())
ss = r['ss_final'] if 'ss_final' in r else r['ss']
for k in ('profit_D','profit_F','nkpc_p_res_D','nkpc_p_res_F','tot_res','union_pi_res'):
    v = float(ss[k]); print(f'  {k} = {v:.3e}'); assert abs(v) < 1e-12, k
print('all new SS residuals are zero')
print(f\"  K_D = {float(ss['K_D']):.10f}\")
print(f\"  rk_D = {float(ss['rk_D']):.10f}\")
print(f\"  w_D = {float(ss['w_D']):.10f}\")
" 2>&1 | tail -20
```

Expected: every new residual `< 1e-12`, then `all new SS residuals are zero`.
Record the printed `K_D`, `rk_D`, `w_D` — compare against
`/tmp/nkpc_baseline_main.log` if it reports them, or against `docs/STATE.md`'s
calibration table. **They must not have moved.**

- [ ] **Step 4: Commit**

```bash
git add code/steady_state.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: seed mc, pi, profit and the global residuals into the SS solve"
```

---

### Task 9: Wire the 27×27 system and pass the equivalence gate

This is the gate. At `kappa_p = 1e4` prices are effectively flexible, so the
27×27 system must reproduce the pre-change baseline. If it does not, the wiring
is wrong and nothing downstream is worth debugging.

**Files:**
- Create: `code/dump_irfs.py`
- Modify: `code/full_model.py` (imports, `build_block_list`, `unknowns_tp`, `targets_tp`)

- [ ] **Step 1: Create the IRF dump script**

Create `code/dump_irfs.py`:

```python
"""Run the pipeline through build_and_solve and save IRFs for comparison.

Usage:
    /opt/anaconda3/envs/ssj/bin/python code/dump_irfs.py OUT.npz
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from calibration import get_calibration
from steady_state import solve_steady_state
from ic_delta_calibration import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model import build_and_solve

KEYS = ('Y_D', 'C_D', 'I_D', 'n_inter_D', 'K_D', 'b_gov_D', 'w_D', 'N_D',
        'p', 'q_b_D', 'spread_rb', 'Y_F', 'C_F', 'I_F', 'n_inter_F')


def main(out_path):
    r = calibrate_depreciation(calibrate_ic_delta(
        solve_steady_state(get_calibration())))
    m = build_and_solve(r)
    payload = {}
    for tag in ('irfs_def_D', 'irfs_Z_D'):
        for k in KEYS:
            if k in m[tag]:
                payload[f'{tag}__{k}'] = np.asarray(m[tag][k])
    for k in ('Y_D', 'C_D', 'n_inter_D', 'K_D'):
        payload[f'ss__{k}'] = np.asarray(float(m['ss_final'][k]))
    np.savez(out_path, **payload)
    print(f'wrote {out_path} with {len(payload)} arrays')


if __name__ == '__main__':
    main(sys.argv[1])
```

- [ ] **Step 2: Dump the pre-change baseline from `main`**

Tasks 2–8 are all committed by now, so there are no uncommitted tracked changes
and `git checkout main` is clean. `code/dump_irfs.py` is still **untracked**, so
it survives the checkout — and it only imports `build_and_solve`, which exists on
`main`. No stash is needed.

```bash
git status --short          # must show only "?? code/dump_irfs.py"
git checkout main
/opt/anaconda3/envs/ssj/bin/python code/dump_irfs.py /tmp/nkpc_irfs_baseline.npz
git checkout add-nkpc
```

Expected: `wrote /tmp/nkpc_irfs_baseline.npz with ...`. If `git status` shows any
tracked modification, commit it before checking out — do not stash.

- [ ] **Step 3: Add the new blocks to `build_block_list`**

In `code/full_model.py`, extend the `from equations_D import (...)` list with
`price_nkpc_D, firm_profit_D`, the `from equations_F import (...)` list with
`price_nkpc_F, firm_profit_F`, and the `from equations_global import (...)` list
with `terms_of_trade, union_inflation`.

Then in `build_block_list`, add `firm_profit_D, price_nkpc_D,` to the country-D
group (immediately after `labor_D, labor_market_D, labor_demand_D,`), add
`firm_profit_F, price_nkpc_F,` to the country-F group in the same position, and
append `terms_of_trade, union_inflation,` to the global group.

- [ ] **Step 4: Update the solver system to 27×27**

Replace `unknowns_tp` and `targets_tp` in `code/full_model.py:94-107`:

```python
    # ── 27×27 system ──────────────────────────────────────────────────────────
    # +4 vs the flex model: mc and pi per country. mc is pinned by the price
    # NKPC, pi jointly by the terms-of-trade identity and the union-inflation
    # normalisation. No targets are renamed or removed -- labor_mkt_res_D/F is
    # unchanged because wages stay flexible.
    unknowns_tp = [
        'K_D', 'n_inter_D', 'div_D', 'I_D', 'Q_D', 'b_gov_D', 'N_D', 'b_F_D', 'w_D', 'rdep_D',
        'mc_D', 'pi_D',
        'K_F', 'n_inter_F', 'div_F', 'I_F', 'Q_F', 'b_gov_F', 'N_F', 'b_D_F', 'w_F', 'rdep_F',
        'mc_F', 'pi_F',
        'p', 'q_b_D', 'q_b_F',
    ]
    targets_tp = [
        'deposit_mkt_D', 'K_res_D', 'n_inter_val_D', 'div_res_D',
        'capital_res_D', 'q_res_D', 'b_gov_res_D', 'b_F_D_res',
        'labor_mkt_res_D', 'w_res_D', 'nkpc_p_res_D',
        'deposit_mkt_F', 'K_res_F', 'n_inter_val_F', 'div_res_F',
        'capital_res_F', 'q_res_F', 'b_gov_res_F', 'b_D_F_res',
        'labor_mkt_res_F', 'w_res_F', 'nkpc_p_res_F',
        'goods_mkt_D', 'rb_D_res', 'rb_F_res', 'tot_res', 'union_pi_res',
    ]
```

Also update the model name on the `sj.create_model` call in `build_and_solve` to
`"Full 2-Country MU HANK — GHH Preferences, Sticky Price, Flex Wage, No CB"`.

- [ ] **Step 5: Verify the system is square**

`unknowns_tp` is local to `build_and_solve`, so count from the source text
directly rather than importing:

```bash
/opt/anaconda3/envs/ssj/bin/python - <<'PY'
import ast, pathlib
tree = ast.parse(pathlib.Path('code/full_model.py').read_text())
found = {}
for node in ast.walk(tree):
    if isinstance(node, ast.Assign) and isinstance(node.value, ast.List):
        for tgt in node.targets:
            if isinstance(tgt, ast.Name) and tgt.id in ('unknowns_tp', 'targets_tp'):
                found[tgt.id] = [e.value for e in node.value.elts]
u, t = found['unknowns_tp'], found['targets_tp']
print('unknowns:', len(u), 'targets:', len(t))
assert len(u) == len(t) == 27, (len(u), len(t))
assert len(set(u)) == 27 and len(set(t)) == 27, 'duplicate entry'
for name in ('mc_D', 'pi_D', 'mc_F', 'pi_F'):
    assert name in u, name
for name in ('nkpc_p_res_D', 'nkpc_p_res_F', 'tot_res', 'union_pi_res',
             'labor_mkt_res_D', 'labor_mkt_res_F'):
    assert name in t, name
print('27x27 confirmed, no duplicates, all new names present')
PY
```

Expected: `unknowns: 27 targets: 27` then
`27x27 confirmed, no duplicates, all new names present`.

- [ ] **Step 6: Run the equivalence gate at `kappa_p = 1e4`**

```bash
KAPPA=1e4 /opt/anaconda3/envs/ssj/bin/python - <<'PY'
import sys, os; sys.path.insert(0, 'code')
import calibration as cal
_orig = cal.get_calibration
def patched():
    c = _orig()
    c['kappa_p_D'] = c['kappa_p_F'] = float(os.environ['KAPPA'])
    return c
cal.get_calibration = patched
import dump_irfs
dump_irfs.main('/tmp/nkpc_irfs_flexlimit.npz')
PY
```

Expected: completes and writes the npz. If the solve fails to converge, step
`KAPPA` down (3e3, 1e3, 3e2) and record the largest value that converges — a very
stiff Phillips curve can be ill-conditioned. Note the value you used.

- [ ] **Step 7: Compare against the baseline**

```bash
/opt/anaconda3/envs/ssj/bin/python - <<'PY'
import numpy as np
a = np.load('/tmp/nkpc_irfs_baseline.npz')
b = np.load('/tmp/nkpc_irfs_flexlimit.npz')
worst = 0.0
for k in sorted(set(a.files) & set(b.files)):
    d = float(np.max(np.abs(a[k] - b[k])))
    scale = max(float(np.max(np.abs(a[k]))), 1e-12)
    rel = d / scale
    worst = max(worst, rel)
    flag = 'FAIL' if rel > 1e-3 else 'ok'
    print(f'{flag:4s} {k:28s} max|abs diff| = {d:.3e}  rel = {rel:.3e}')
print(f'\nworst relative deviation = {worst:.3e}')
PY
```

Expected: every line `ok`, worst relative deviation below 1e-3. This is a
*limit* comparison, not an exact one — a finite `kappa_p` cannot reproduce the flex
model to machine precision. If the worst deviation exceeds 1e-3, re-run step 6
with a larger `KAPPA` and check whether the deviation shrinks proportionally. **If
it does not shrink with `kappa_p`, the wiring is wrong — stop and debug.**

- [ ] **Step 8: Run the full pipeline once at the flex limit to check residuals**

```bash
KAPPA=1e4 /opt/anaconda3/envs/ssj/bin/python - <<'PY' 2>&1 | tee /tmp/nkpc_flexlimit_main.log
import sys, os; sys.path.insert(0, 'code')
import calibration as cal
_orig = cal.get_calibration
def patched():
    c = _orig()
    c['kappa_p_D'] = c['kappa_p_F'] = float(os.environ['KAPPA'])
    return c
cal.get_calibration = patched
import main
main.main()
PY
grep -E "goods_mkt|ca_res|deposit_mkt|n_inter_D\[0\]|Y_D\[0\]|b_gov_D\[499\]" /tmp/nkpc_flexlimit_main.log
```

Expected: `goods_mkt_D` ≤ 1e−14, `goods_mkt_F` and `ca_res_D` ≤ 1e−7,
`deposit_mkt_D/F` ≤ 1e−13, `n_inter_D[0]` and `Y_D[0]` both negative.

- [ ] **Step 9: Commit**

```bash
git add code/full_model.py code/dump_irfs.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: 27x27 sticky-price system; passes the kappa_p -> inf equivalence gate"
```

---

### Task 9b: Convert every Jacobian call site to the padded solver

Task 9 converted `full_model.py` and `tpi.py`. **Seven call sites still use stock
`solve_jacobian` and will die with the core-dimension mismatch the moment they
see the sticky-price system.** They are not broken yet only because they have not
been re-run.

Two of them block later tasks outright:
- `diagnostics/regimes/regime_model.py:177` — Task 15's cache rebuild, and hence
  **all of E1–E4**, runs off this.
- `experiments/e4_distribution.py:255` — E4's quintile incidence.

The rest are diagnostics off the plan's critical path but must not be left as
landmines: `diagnostics/solve_configs.py:176`,
`diagnostics/psilam_moment_sweep.py:76`,
`diagnostics/psilam_breakdown_sweep.py:83`,
`diagnostics/substitution_v2/solve_v2.py:106`,
`diagnostics/substitution_v2/exp_psilam0.py:64`.

- [ ] **Step 1: Convert each call site**

In each file, replace `<model>.solve_jacobian(ss, unknowns=..., targets=...,
inputs=..., T=...)` with:

```python
from full_model import solve_jacobian_padded
G = solve_jacobian_padded(<model>, ss, <unknowns>, <targets>, <inputs>, T)
```

Preserve each site's own variable names and any extra keyword arguments it
passes. Note `solve_jacobian_padded` takes `unknowns`, `targets`, `inputs`, `T`
positionally after `model` and `ss`.

- [ ] **Step 2: Verify none remain**

```bash
grep -rn "\.solve_jacobian(" --include="*.py" code experiments diagnostics \
  | grep -v solve_jacobian_padded
```

Expected: **no output**.

- [ ] **Step 3: Smoke-test the two on the critical path**

```bash
/opt/anaconda3/envs/ssj/bin/python -c "
import sys
sys.path.insert(0,'code'); sys.path.insert(0,'diagnostics/regimes'); sys.path.insert(0,'experiments')
import regime_model, e4_distribution
print('regime_model and e4_distribution import cleanly')
"
```

- [ ] **Step 4: Commit** (all three docs staged, as always)

```bash
git add code experiments diagnostics docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "fix: route every Jacobian call site through solve_jacobian_padded

SSJ 1.0.0 drops H_Z rows for targets reachable from no shock. Seven call sites
would have hit the core-dimension mismatch on first contact with the 27x27
system; regime_model.py blocks the E1-E4 cache rebuild."
```

### Task 10: Dial `kappa_p` to 0.0871 and record the price-stickiness result

This is the clean measure of what price stickiness alone does, with deposits still
real. It is a reportable result, not just a checkpoint — the benchmark literature
does not publish a flexible-price counterfactual.

**Files:** none modified — `kappa_p` is already 0.0871 in `calibration.py`.

- [ ] **Step 1: Run the full pipeline at the calibrated slope**

```bash
/opt/anaconda3/envs/ssj/bin/python code/main.py 2>&1 | tee /tmp/nkpc_sticky_main.log
```

Expected: completes.

- [ ] **Step 2: Check every acceptance threshold and sign**

```bash
grep -E "goods_mkt|ca_res|deposit_mkt|IC-δ|IC-delta" /tmp/nkpc_sticky_main.log
grep -E "n_inter_D\[0\]|Y_D\[0\]|b_gov_D\[499\]|rho_b" /tmp/nkpc_sticky_main.log
```

Expected: thresholds as in Task 9 step 8; `n_inter_D[0] < 0`; `Y_D[0] < 0`;
`b_gov_D[499]` near zero; the IC-δ consistency check unchanged.
If `assert_gk_well_posed` raises, the GK block is no longer well-posed — stop and
report, do not paper over it.

- [ ] **Step 3: Record the headline impact numbers**

```bash
/opt/anaconda3/envs/ssj/bin/python code/dump_irfs.py /tmp/nkpc_irfs_sticky.npz
/opt/anaconda3/envs/ssj/bin/python - <<'PY'
import numpy as np
b = np.load('/tmp/nkpc_irfs_baseline.npz')
s = np.load('/tmp/nkpc_irfs_sticky.npz')
print(f"{'variable':12s} {'flex (% SS)':>14s} {'sticky (% SS)':>15s}")
for k in ('Y_D', 'C_D', 'I_D', 'n_inter_D'):
    ssv = float(b[f'ss__{k}']) if f'ss__{k}' in b.files else 1.0
    flex   = b[f'irfs_def_D__{k}'][0] / ssv * 100
    sticky = s[f'irfs_def_D__{k}'][0] / ssv * 100
    print(f'{k:12s} {flex:+14.4f} {sticky:+15.4f}')
PY
```

Expected: a four-row table. **Record it verbatim** — it goes into `docs/STATE.md`
in Task 16 and is the price-stickiness-only result. `C_D[0]` is `+0.2164%` in the
flex baseline; report whatever it is now without editorialising.

- [ ] **Step 4: Sweep `kappa_p` for the robustness table**

```bash
for K in 0.03 0.0871 0.2; do
  echo "=== kappa_p = $K ==="
  KAPPA=$K /opt/anaconda3/envs/ssj/bin/python - <<'PY' 2>&1 | grep -E "n_inter_D\[0\]|Y_D\[0\]|b_gov_D\[499\]"
import sys, os; sys.path.insert(0, 'code')
import calibration as cal
_orig = cal.get_calibration
def patched():
    c = _orig()
    c['kappa_p_D'] = c['kappa_p_F'] = float(os.environ['KAPPA'])
    return c
cal.get_calibration = patched
from steady_state import solve_steady_state
from ic_delta_calibration import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model import build_and_solve
build_and_solve(calibrate_depreciation(calibrate_ic_delta(
    solve_steady_state(cal.get_calibration()))))
PY
done
```

Expected: three blocks of output. Record which values keep `b_gov_D[499]` near
zero and both signs negative — that is the stable region for the robustness table.

- [ ] **Step 5: Commit the recorded numbers**

Add the two tables to `docs/PROGRESS.md` under the current date.

```bash
git add docs/PROGRESS.md
git commit -m "docs: record price-stickiness-only impact numbers and kappa_p sweep"
```

---

## Phase 2 — Nominal deposits

### Task 11: `deposit_rates_{D,F}` and nominal `deposit_return_{D,F}`

`rdep_i` becomes a **derived** ex-ante real rate and `i_dep_i` becomes the nominal
unknown. This is the refinement noted at the top of the plan: because
`intermediation_P1`, `divert_bond_foc` and `divert_portfolio_adj` already mean
"ex-ante real" by `rdep_i`, they need no changes at all.

**Files:**
- Modify: `code/equations_D.py:75-84`, `code/equations_F.py:69-76`
- Modify: `code/test_nkpc_blocks.py`

- [ ] **Step 1: Write the failing tests**

Append to `code/test_nkpc_blocks.py`:

```python
# ── Nominal deposits ──────────────────────────────────────────────────────────

def test_deposit_rates_collapse_at_zero_inflation():
    """At pi = 0 both derived real rates must equal the nominal rate exactly --
    this is what keeps the steady state bit-identical."""
    from equations_D import deposit_rates_D
    ss = deposit_rates_D.steady_state({'i_dep_D': 0.0125, 'pi_D': 0.0})
    assert ss['rdep_D'] == pytest.approx(0.0125, rel=1e-15)
    assert ss['rdep_expost_D'] == pytest.approx(0.0125, rel=1e-15)


def test_deflation_raises_the_realised_real_deposit_rate():
    """Deflation is a windfall to depositors and a loss to banks, which hold
    real assets against nominal liabilities. This is the Fisher channel; if the
    sign flips, bank_return_D will amplify in the wrong direction."""
    from equations_D import deposit_rates_D
    i = 0.0125
    base = deposit_rates_D.steady_state({'i_dep_D': i, 'pi_D': 0.0})
    defl = deposit_rates_D.steady_state({'i_dep_D': i, 'pi_D': -0.001})
    assert defl['rdep_expost_D'] > base['rdep_expost_D']
    assert defl['rdep_expost_D'] == pytest.approx((1 + i) / (1 - 0.001) - 1, rel=1e-14)


def test_deposit_return_is_unchanged_at_zero_inflation():
    """Rgross must be exactly 1 + i_dep when pi = 0 and P_CES is flat."""
    from equations_D import deposit_return_D
    ss = deposit_return_D.steady_state({'i_dep_D': 0.0125, 'P_CES_D': 1.3, 'pi_D': 0.0})
    assert ss['Rgross_D'] == pytest.approx(1.0125, rel=1e-15)


def test_deposit_rates_F_matches_D():
    from equations_D import deposit_rates_D
    from equations_F import deposit_rates_F
    d = deposit_rates_D.steady_state({'i_dep_D': 0.0125, 'pi_D': -0.001})
    f = deposit_rates_F.steady_state({'i_dep_F': 0.0125, 'pi_F': -0.001})
    assert d['rdep_D'] == pytest.approx(f['rdep_F'], rel=1e-15)
    assert d['rdep_expost_D'] == pytest.approx(f['rdep_expost_F'], rel=1e-15)
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -k "deposit" -v
```

Expected: FAIL — `ImportError: cannot import name 'deposit_rates_D'`.

- [ ] **Step 3: Implement the blocks**

Replace `code/equations_D.py:75-84` in full:

```python
@simple
def deposit_rates_D(i_dep_D, pi_D):
    # Deposits are NOMINAL euro contracts. i_dep_D is the nominal rate and is the
    # unknown that clears deposit_mkt_D -- there is no policy rate pinning it, so
    # no absorber or cross-border claim is needed and external_account_D is
    # untouched.
    #
    # rdep_D keeps its existing meaning: the EX-ANTE real rate for the t -> t+1
    # holding period, locked at t. That is exactly what intermediation_P1_D,
    # divert_bond_foc_D and divert_portfolio_adj already mean by rdep_D, so those
    # blocks need no changes.
    #
    # rdep_expost_D is the REALISED real rate at t on deposits placed at t-1. It
    # contains the inflation surprise: a deflation raises the real value of the
    # bank's nominal liabilities. Banks hold real assets against nominal
    # liabilities, so they are net nominal debtors and this deepens the net-worth
    # loss -- the Fisher-Bernanke channel.
    #
    # At SS pi_D = 0 and both equal i_dep_D, so the SS is bit-identical.
    rdep_D        = (1 + i_dep_D) / (1 + pi_D(+1)) - 1
    rdep_expost_D = (1 + i_dep_D(-1)) / (1 + pi_D) - 1
    return rdep_D, rdep_expost_D


@simple
def deposit_return_D(i_dep_D, P_CES_D, pi_D):
    # Bundle-real gross deposit return on a NOMINAL contract.
    # P_c_D = P_D * P_CES_D is the nominal CPI, so
    #   P_c_D(-1)/P_c_D = (P_CES_D(-1)/P_CES_D) / (1 + pi_D).
    #
    # T-2 is NOT reopened: the rate is still locked at t-1 (i_dep_D(-1)); only
    # the deflator is period-t, which this block already did via P_CES. T-2 was
    # about paying a period-t UNKNOWN rate on the t-1 deposit stock.
    #
    # At SS P_CES_D(-1)/P_CES_D = 1 and pi_D = 0, so Rgross_D = 1 + i_dep_D.
    Rgross_D = (1 + i_dep_D(-1)) * P_CES_D(-1) / P_CES_D / (1 + pi_D)
    return Rgross_D
```

Replace `code/equations_F.py:69-76` in full:

```python
@simple
def deposit_rates_F(i_dep_F, pi_F):
    # See deposit_rates_D.
    rdep_F        = (1 + i_dep_F) / (1 + pi_F(+1)) - 1
    rdep_expost_F = (1 + i_dep_F(-1)) / (1 + pi_F) - 1
    return rdep_F, rdep_expost_F


@simple
def deposit_return_F(i_dep_F, P_CES_F, pi_F):
    # See deposit_return_D. Nominal contract; T-2 timing preserved.
    Rgross_F = (1 + i_dep_F(-1)) * P_CES_F(-1) / P_CES_F / (1 + pi_F)
    return Rgross_F
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add code/equations_D.py code/equations_F.py code/test_nkpc_blocks.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: nominal deposit contracts; ex-ante and ex-post real rates"
```

---

### Task 12: Ex-post funding cost in `bank_return` and `capital_fund`

These two blocks are where the bank pays for deposits it took on at t−1, so they
must use the *realised* real rate. Everything else that touches `rdep_i` is
forward-looking and stays as it is.

**Files:**
- Modify: `code/equations_D.py:325-348`
- Modify: `code/equations_F.py:276-302`

- [ ] **Step 1: Write the failing test**

Append to `code/test_nkpc_blocks.py`:

```python
def test_bank_return_uses_the_expost_rate():
    """Signature check: bank_return_D must take rdep_expost_D and must NOT take
    rdep_D. Getting this backwards silently reverses the Fisher channel."""
    from equations_D import bank_return_D, capital_fund_D
    for blk in (bank_return_D, capital_fund_D):
        assert 'rdep_expost_D' in blk.inputs, (blk.name, blk.inputs)
        assert 'rdep_D' not in blk.inputs, (blk.name, blk.inputs)


def test_forward_looking_blocks_still_use_rdep():
    """intermediation_P1_D and divert_bond_foc_D are ex-ante and must be
    untouched -- rdep_D still means the t -> t+1 real rate."""
    from equations_D import intermediation_P1_D, divert_bond_foc_D
    for blk in (intermediation_P1_D, divert_bond_foc_D):
        assert 'rdep_D' in blk.inputs, (blk.name, blk.inputs)
        assert 'rdep_expost_D' not in blk.inputs, (blk.name, blk.inputs)
```

- [ ] **Step 2: Run to verify it fails**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -k "expost_rate or forward_looking" -v
```

Expected: `test_bank_return_uses_the_expost_rate` FAILS;
`test_forward_looking_blocks_still_use_rdep` PASSES already.

- [ ] **Step 3: Modify `bank_return_D`**

In `code/equations_D.py:325-336`, change the signature `rdep_D` → `rdep_expost_D`
and replace every `rdep_D(-1)` with `rdep_expost_D`:

```python
@simple
def bank_return_D(theta_D, rk_D, rdep_expost_D, b_D_D, b_F_D, n_inter_D,
                  rb_actual_D, rb_actual_F, q_b_D, q_b_F):
    phi_bD_lag_D = q_b_D(-1) * b_D_D(-1) / n_inter_D(-1)
    phi_bF_lag_D = q_b_F(-1) * b_F_D(-1) / n_inter_D(-1)
    kappa_lag_D  = theta_D(-1) - phi_bD_lag_D - phi_bF_lag_D
    # T-2 fix: funding cost on the t-1 balance sheet is the rate locked at t-1.
    # Under nominal deposits that realised real cost is rdep_expost_D, which
    # already carries the (-1) timing internally and contains the inflation
    # surprise -- the Fisher revaluation on the bank's nominal liabilities.
    rn_D = (kappa_lag_D  * (rk_D        - rdep_expost_D)
            + phi_bD_lag_D * (rb_actual_D - rdep_expost_D)
            + phi_bF_lag_D * (rb_actual_F - rdep_expost_D)
            + rdep_expost_D)
    return rn_D
```

- [ ] **Step 4: Modify `capital_fund_D`**

In `code/equations_D.py:339-348`:

```python
@simple
def capital_fund_D(rk_D, rdep_expost_D, Q_D, K_D, omega_K_D, fund_rule_D, K_fund_D):
    # Passive capital fund funded by deposits; rebates its spread on the lagged
    # capital value to households. Same predetermined-rate timing as
    # bank_return_D (T-2); rdep_expost_D is the realised real funding cost under
    # nominal deposits. Zero when the fund is empty (omega_K_D=1, K_fund_D=0).
    # fund_rule_D: 0 = fund holds (1-omega_K)·K, 1 = fund holds a constant K_fund.
    K_fnd_lag_D = ((1.0 - fund_rule_D) * (1.0 - omega_K_D) * K_D(-1)
                   + fund_rule_D * K_fund_D)
    div_fund_D = (rk_D - rdep_expost_D) * Q_D(-1) * K_fnd_lag_D
    return div_fund_D
```

- [ ] **Step 5: Apply the identical changes to F**

In `code/equations_F.py:276-292`, change `bank_return_F`'s signature `rdep_F` →
`rdep_expost_F` and replace all four `rdep_F(-1)` with `rdep_expost_F`. In
`code/equations_F.py:295-302`, change `capital_fund_F`'s signature `rdep_F` →
`rdep_expost_F` and replace `rdep_F(-1)` with `rdep_expost_F`. Do **not** touch
the W-2 `p(-1)/p` conversions in `bank_return_F`.

- [ ] **Step 6: Run the tests**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v
/opt/anaconda3/envs/ssj/bin/python -c "
import sys; sys.path.insert(0, 'code')
from equations_F import bank_return_F, capital_fund_F
for b in (bank_return_F, capital_fund_F):
    assert 'rdep_expost_F' in b.inputs and 'rdep_F' not in b.inputs, (b.name, b.inputs)
print('F blocks wired to rdep_expost_F')
"
```

Expected: all tests pass, then `F blocks wired to rdep_expost_F`.

- [ ] **Step 7: Commit**

```bash
git add code/equations_D.py code/equations_F.py code/test_nkpc_blocks.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: bank_return and capital_fund pay the ex-post real deposit rate"
```

---

### Task 13: Wire nominal deposits and verify the Fisher sign

**Files:**
- Modify: `code/calibration.py:70`
- Modify: `code/steady_state.py` (imports + SS block list)
- Modify: `code/full_model.py` (imports, `build_block_list`, `unknowns_tp`)

- [ ] **Step 1: Rename the calibration parameter**

Replace `code/calibration.py:70`:

```python
        # Nominal deposit rate. Deposits are nominal euro contracts; the derived
        # real rates rdep_D/F (ex-ante) and rdep_expost_D/F (realised) come from
        # deposit_rates_D/F. At SS pi = 0, so rdep = i_dep and the SS is
        # unchanged from the real-deposit calibration.
        'i_dep_D':      0.000,   'i_dep_F':      0.000,
```

- [ ] **Step 2: Add `deposit_rates` to both block lists**

In `code/steady_state.py`, add `deposit_rates_D` / `deposit_rates_F` to the
imports and place them in the SS `create_model` list immediately **before**
`deposit_return_D, deposit_return_F` (they produce `rdep_D` / `rdep_F`, which
`smart_steady_D/F` and `steady_auxilliary_D/F` consume):

```python
        deposit_rates_D, deposit_rates_F,
        deposit_return_D, deposit_return_F,
```

In `code/full_model.py`, add `deposit_rates_D` / `deposit_rates_F` to the imports
and put each immediately before its `deposit_return_*` entry in
`build_block_list`.

- [ ] **Step 3: Swap the unknown**

In `code/full_model.py`'s `unknowns_tp`, replace `'rdep_D'` with `'i_dep_D'` and
`'rdep_F'` with `'i_dep_F'`. `targets_tp` is unchanged — still 27×27.

- [ ] **Step 4: Verify the steady state has not moved**

```bash
/opt/anaconda3/envs/ssj/bin/python -c "
import sys; sys.path.insert(0, 'code')
from calibration import get_calibration
from steady_state import solve_steady_state
r = solve_steady_state(get_calibration())
ss = r['ss_final'] if 'ss_final' in r else r['ss']
i, ra, rp = float(ss['i_dep_D']), float(ss['rdep_D']), float(ss['rdep_expost_D'])
print(f'  i_dep_D={i:.12f}  rdep_D={ra:.12f}  rdep_expost_D={rp:.12f}')
assert abs(ra - i) < 1e-14 and abs(rp - i) < 1e-14, 'rates do not collapse at SS'
print(f\"  K_D = {float(ss['K_D']):.10f}\")
print(f\"  rk_D = {float(ss['rk_D']):.10f}\")
print(f\"  w_D = {float(ss['w_D']):.10f}\")
print('SS rates collapse correctly')
" 2>&1 | tail -10
```

Expected: all three rates equal, and `K_D` / `rk_D` / `w_D` identical to the
values recorded in Task 8 step 3.

- [ ] **Step 5: Run the full pipeline**

```bash
/opt/anaconda3/envs/ssj/bin/python code/main.py 2>&1 | tee /tmp/nkpc_nominal_main.log
grep -E "goods_mkt|ca_res|deposit_mkt|n_inter_D\[0\]|Y_D\[0\]|b_gov_D\[499\]" /tmp/nkpc_nominal_main.log
```

Expected: all acceptance thresholds hold; both signs negative.

- [ ] **Step 6: Verify the Fisher sign — the gate for this phase**

```bash
/opt/anaconda3/envs/ssj/bin/python code/dump_irfs.py /tmp/nkpc_irfs_nominal.npz
/opt/anaconda3/envs/ssj/bin/python - <<'PY'
import numpy as np
s = np.load('/tmp/nkpc_irfs_sticky.npz')     # Task 10: sticky prices, real deposits
n = np.load('/tmp/nkpc_irfs_nominal.npz')    # this task: + nominal deposits
ss = float(s['ss__n_inter_D'])
a = s['irfs_def_D__n_inter_D'][0] / ss * 100
b = n['irfs_def_D__n_inter_D'][0] / ss * 100
print(f'n_inter_D[0]  real deposits: {a:+.4f}% of SS')
print(f'n_inter_D[0]  nominal      : {b:+.4f}% of SS')
assert b < a, 'FISHER SIGN WRONG: nominal deposits must deepen the net-worth loss'
print('Fisher channel sign OK')
for k in ('Y_D', 'C_D', 'I_D'):
    v = float(s[f'ss__{k}']) if f'ss__{k}' in s.files else 1.0
    print(f'{k:6s} real {s[f"irfs_def_D__{k}"][0]/v*100:+.4f}%   '
          f'nominal {n[f"irfs_def_D__{k}"][0]/v*100:+.4f}%')
PY
```

Expected: `Fisher channel sign OK`. **If the assertion fires, the ex-post /
ex-ante substitution in Task 12 is backwards — do not proceed.** Record the
printed `Y_D` / `C_D` / `I_D` comparison; it is the Fisher-channel result.

- [ ] **Step 7: Commit**

```bash
git add code/calibration.py code/steady_state.py code/full_model.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "feat: nominal deposits wired; Fisher channel deepens the net-worth loss"
```

---

## Phase 3 — Recalibration, regeneration, documentation

### Task 14: Re-tune `psi_lambda_B` to the 150bp target

Spread transmission now runs through both a sticky terms of trade and a Fisher
revaluation, so the amplification dial has to be re-disciplined. The moment is
**peak annualised D−F spread ≈ 150bp on a 1pp default-probability shock**.

CLAUDE.md puts the documented breakdown region around 4–5 at `n_inter = 3.0`.
Hitting the moment is not sufficient — stability must be re-verified at whatever
value it lands on.

**Files:** Modify `code/calibration.py` (`psi_lambda_B_D` / `psi_lambda_B_F`)

- [ ] **Step 1: Find the current value**

```bash
grep -n "psi_lambda_B" code/calibration.py
```

Record it.

- [ ] **Step 2: Measure the spread at three candidate values**

```bash
for PSI in 1.5 3.0 5.0; do
  echo "=== psi_lambda_B = $PSI ==="
  PSI=$PSI /opt/anaconda3/envs/ssj/bin/python - <<'PY' 2>&1 | tail -6
import sys, os; sys.path.insert(0, 'code')
import numpy as np
import calibration as cal
_orig = cal.get_calibration
def patched():
    c = _orig()
    c['psi_lambda_B_D'] = c['psi_lambda_B_F'] = float(os.environ['PSI'])
    return c
cal.get_calibration = patched
from steady_state import solve_steady_state
from ic_delta_calibration import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model import build_and_solve
m = build_and_solve(calibrate_depreciation(calibrate_ic_delta(
    solve_steady_state(cal.get_calibration()))))
sp = m['irfs_def_D']['spread_rb']
print(f"peak spread = {np.max(np.abs(sp)) * 400 * 100:.1f} bp annualised")
print(f"b_gov_D[499] = {m['irfs_def_D']['b_gov_D'][499]:.3e}")
PY
done
```

Expected: three `peak spread` readings. The `* 400 * 100` converts a quarterly
rate deviation to annualised basis points; cross-check the first reading against
the flex baseline's documented 150.3bp to confirm the scaling before trusting it.

- [ ] **Step 3: Bisect to 150bp**

Peak spread is monotone increasing in `psi_lambda_B`. Take the bracketing pair
from step 2 and bisect, re-running step 2's script with the midpoint, until the
peak spread is within 1bp of 150. Record every (psi, spread) pair evaluated.

- [ ] **Step 4: Verify stability at the tuned value**

Run the full pipeline at the tuned `psi_lambda_B`:

```bash
/opt/anaconda3/envs/ssj/bin/python code/main.py 2>&1 | tee /tmp/nkpc_tuned_main.log
grep -E "goods_mkt|ca_res|deposit_mkt|n_inter_D\[0\]|Y_D\[0\]|b_gov_D\[499\]|rho_b" /tmp/nkpc_tuned_main.log
```

Expected: all thresholds hold, both signs negative, `b_gov_D[499]` near zero.
**If the tuned value sits in the 4–5 breakdown region, stop and report** — hitting
the moment inside a breakdown region is not a valid calibration.

- [ ] **Step 5: Sweep `omega_pi_D` as the containing parameter**

```bash
for W in 0.071 0.2 0.5; do
  echo "=== omega_pi_D = $W ==="
  OMEGA=$W /opt/anaconda3/envs/ssj/bin/python - <<'PY' 2>&1 | grep -E "n_inter_D\[0\]|Y_D\[0\]|b_gov_D\[499\]"
import sys, os; sys.path.insert(0, 'code')
import calibration as cal
_orig = cal.get_calibration
def patched():
    c = _orig(); c['omega_pi_D'] = float(os.environ['OMEGA']); return c
cal.get_calibration = patched
from steady_state import solve_steady_state
from ic_delta_calibration import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model import build_and_solve
build_and_solve(calibrate_depreciation(calibrate_ic_delta(
    solve_steady_state(cal.get_calibration()))))
PY
done
```

Expected: three blocks. Record which values stay stable — `omega_pi_D` scales
`pi_D` and hence the Fisher revaluation, so this is the containing parameter if
the Fisher channel turns out to dominate.

- [ ] **Step 6: Commit the tuned calibration**

```bash
git add code/calibration.py docs/STATE.md docs/PROGRESS.md docs/HANDOFF.md
git commit -m "calib: re-tune psi_lambda_B to the 150bp spread moment under sticky prices"
```

---

### Task 15: Rebuild the regime cache and regenerate all results

**Files:** none modified — this regenerates `docs/experiments_results.md` and figures.

- [ ] **Step 1: Rebuild the cached Jacobian response matrices**

```bash
/opt/anaconda3/envs/ssj/bin/python diagnostics/regimes/regime_model.py --force 2>&1 | tail -20
```

Expected: completes and writes new `cache_G_main_v*.npz` files. This takes a long
time (multiple Jacobian solves). The `experiments/` package reads these, so it
**must** run before `run_all.py` or the experiments will silently report the old
flex-price model.

- [ ] **Step 2: Regenerate the standard results set**

```bash
/opt/anaconda3/envs/ssj/bin/python experiments/run_all.py 2>&1 | tail -40
```

Expected: completes and rewrites `docs/experiments_results.md`. E2 self-verifies
and asserts its dY decomposition closes at 1e−7 — since `market_clearing_D` is
untouched by this work (no Rotemberg resource cost was added), that assertion is
an independent check on the whole change. **If E2's assertion fires, something has
leaked into the resource constraint — stop and find it.**

- [ ] **Step 3: Run the experiments' own tests**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest experiments/ code/test_nkpc_blocks.py code/test_eba_calibration.py -v
```

Expected: all pass.

- [ ] **Step 4: Regenerate figures**

Already done by `code/main.py` in Task 14 step 4 (steps 6 and 8 write to
`outputs/`). Confirm the files are newer than the calibration change:

```bash
ls -la outputs/ | head -20
```

- [ ] **Step 5: Commit**

```bash
git add docs/experiments_results.md outputs/ diagnostics/regimes/*.npz docs/PROGRESS.md
git commit -m "regen: rebuild regime cache and regenerate E1-E4 on the sticky model"
```

If the `.npz` cache files are gitignored, drop them from the `git add` — check
`diagnostics/regimes/.gitignore` first.

---

### Task 16: Documentation

The doc hooks require STATE.md, PROGRESS.md and HANDOFF.md on any commit that
stages Python. This task writes them properly rather than the one-liners used
during the phases.

**Files:** Modify `docs/STATE.md`, `docs/PROGRESS.md`, `docs/HANDOFF.md`,
`docs/SPEC.md`, `CLAUDE.md`

- [ ] **Step 1: `docs/STATE.md`**

Add a new section dated 2026-08-05 containing: the new calibration rows (`mu_p`,
`kappa_p`, `mc`, `omega_pi_D`, `i_dep`, tuned `psi_lambda_B`); the three-way impact
table from Task 10 step 3 and Task 13 step 6 (flex / sticky-real-deposits /
sticky-nominal-deposits, for `Y_D[0]`, `C_D[0]`, `I_D[0]`, `n_inter_D[0]`); the
`kappa_p` and `omega_pi_D` sweeps; and the post-change Walras residuals.

State explicitly whether `C_D[0]` changed sign. If it did not, say so plainly and
cross-reference the spec's *"The `C_D[0]` motivation, stated honestly"* section —
Bi-Foerster-Traum get consumption rising on impact too, with a Taylor rule, nominal
debt and a loan-in-advance constraint.

- [ ] **Step 2: `docs/PROGRESS.md`**

Consolidate the per-task one-liners into a single dated changelog entry covering
the refactor, the price NKPCs, the markup rent routing, nominal deposits, the
`psi_lambda_B` re-tune, and the regeneration.

- [ ] **Step 3: `docs/HANDOFF.md`**

Update the incidence paragraph at line 65. It currently says "every quintile's
consumption *rises* on impact ... must be confronted in the draft." Replace with
the post-change finding and note what was tried (sticky prices, then nominal
deposits) and what remains untried (a Sims-Wu loan-in-advance constraint; nominal
sovereign bonds).

- [ ] **Step 4: `docs/SPEC.md`**

Under *Key modelling choices*, add: the price Phillips curve and its
subsidy-neutralised steady state; the union-inflation normalisation as the nominal
anchor and why there is no Taylor rule; nominal deposits against real bonds as a
deliberate asymmetry that must be stated in the paper; and the markup rent's
distribution rule.

- [ ] **Step 5: `CLAUDE.md`**

Update the *Architecture* section to mention `build_block_list()` as the single
model definition. Add `mu_p`, `kappa_p`, `omega_pi_D` and `i_dep` to the key
modelling choices. Update the *Typical iteration* residual list to include the
four new targets. Add `code/test_nkpc_blocks.py` to *Running and testing* as the
fast unit-test entry point:
`/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -v`.

- [ ] **Step 6: Final verification and commit**

```bash
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py experiments/ -v
/opt/anaconda3/envs/ssj/bin/python code/main.py 2>&1 | tail -30
```

Expected: tests pass; pipeline completes with all thresholds held.

```bash
git add docs/ CLAUDE.md
git commit -m "docs: sticky prices and nominal deposits become the baseline

Records the three-way impact comparison (flex / sticky-real / sticky-nominal),
the re-tuned psi_lambda_B, the kappa_p and omega_pi_D sweeps, and states the
C_D[0] outcome plainly."
```

---

## Notes on what is deliberately NOT here

- **Wage rigidity.** Author decision 2026-08-05. `labor_market_{D,F}` is untouched
  and stays in every block list.
- **A Taylor rule.** Author decision 2026-08-05. A policy rule only bites if it
  pins a real rate, which frees both deposit-market conditions and needs either a
  cross-border banking claim (rewriting `external_account_D`) or an ECB reserve
  asset inside the GK incentive constraint.
- **Nominal sovereign bonds.** Deliberate asymmetry; candidate follow-on spec.
- **Rotemberg resource costs.** Quadratic around `pi_ss = 0`, so first-order
  irrelevant under `solve_jacobian`, and including them would inject a nonlinear
  term into `goods_mkt_D` which holds at 1e−14.
- **A Sims-Wu loan-in-advance constraint.** The natural next lever if Tasks 10 and
  13 leave `Y_D[0]` implausibly small; needs its own design pass.
