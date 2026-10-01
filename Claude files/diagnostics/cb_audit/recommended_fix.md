# Central Bank block audit — recommended fixes

**None of these are implemented.** Audit-only; no model source was edited. Findings
and evidence are in `VERDICT.md` and `run_log.md`.

Ordered by what actually needs doing, not by finding number.

---

## R-1 (F-2) — Replace the condition-number pole scan with the exact eigenvalue test, and trim the terminal window before testing

**Two separate defects, one fix each.**

### R-1a. Detect poles exactly, not on a grid

`code/tpi.py:344-353` and `diagnostics/regimes/lottery_math.closed_loop_pole:55`
locate the singularity of `I − γ·A_cb` by scanning 240 grid points for
`cond > 1e4`. The poles of that pencil are known in closed form: for every real
positive eigenvalue `λ` of `A_cb`, `γ = 1/λ` is exactly singular. One `eigvals` call
at `O(T³)` — the same cost as one `cond` evaluation, of which the current code does
240 — replaces the whole scan and cannot step over anything.

```python
# sketch, NOT applied
w = np.linalg.eigvals(A_cb)
poles = sorted(1.0 / x.real for x in w
               if abs(x.imag) < 1e-10 and x.real > 1e-8)
gamma_pole = poles[0] if poles else np.inf
```

This alone would have surfaced the γ = 2.2116 root the current scan misses.

### R-1b. Test the interior operator, not the terminal boundary

`A_cb`'s last few columns are truncation boundary, not model: `||A[:,499]|| = 3.86`
against ~0.0065 for every interior column, and `A[499,499] = +1.080` is the only
positive diagonal entry and the only entry with modulus above 1. Both "poles" —
γ = 2.2116 and γ ≈ 27.4 — come entirely from those columns; dropping five of them
removes every pole below γ = 36 and changes the reported peak spread by nothing at
γ = 2, 5, 10.

So R-1a should be applied to a trimmed operator, e.g. `A_cb[:T-T_trim, :T-T_trim]`
with `T_trim` on the order of 10–20, and the trim must be **stated in the printed
diagnostic** so it can never go silent — the same discipline
`solve_jacobian_padded` already uses for its zero-padded rows.

**Validation before adopting:** re-run the effectiveness curve with the trimmed
detector and confirm `peak_arr` is unchanged to ≤ 0.05 bp at every γ on the existing
grid. The audit already checked this at γ = 2, 5, 10, 20; the fix should check the
whole grid.

**Consequences to follow through, all documentation:**

* the γ = 19.88 cap on the effectiveness curve is imposed for a spurious reason and
  should be removed or re-derived (the curve may legitimately extend further; whether
  it *should* is a linearisation question, not a pole question);
* `CLAUDE.md`'s and `code/tpi.py:332-343`'s claim of "a closed-loop pole at
  gamma ~ 27.3 on the post-GK-refactor calibration" is an artefact reading and must be
  corrected, not softened;
* `docs/STATE.md` and any paper text asserting a stability ceiling on γ derived from
  that pole needs the same treatment.

**Do not** simply tighten the grid. A finer grid still misses poles between points; it
only moves the threshold at which the same failure recurs.

---

## R-2 (F-1) — Separate the realised transfer from the expected loss, in the code's own output

This is a reporting fix, not a model fix. The model is correct; the printed table
invites the wrong reading.

`code/tpi.py:300-323` prints `EL PV`, `prem PV`, `MTM PV` and the per-γ line
"F bears EL PV = … , receives prem PV = …" in one block, with no marker separating
the on-path objects from the off-path ones. `el_pv` is an expectation over a default
event that `writeoff_enabled_D = 0` guarantees never occurs on the traced path; it
appears in no budget constraint. `prem_pv`, `mtm_pv` and `carry_ss_pv` are realised
cash flows that do move `rem_cb_F` and `b_gov_F`.

Recommended:

1. Split the table into an **ON-PATH** group (`purchases_pv`, `prem_pv`,
   `carry_ss_pv`, `mtm_pv` — the flows that actually pass through
   `budget_residual_D_tpi` / `_F_tpi`) and an **OFF-PATH / EXPECTATION** group
   (`el_pv`, and the `loading` ratio built from it), with the gate stated in the
   header: "`writeoff_enabled = 0`: no credit loss is realised on this branch".
2. Add a direct on-path incidence line computed from the conduit itself rather than
   reconstructed — `rem_cb_F` and `rem_cb_D` are already model outputs with
   `cb_buy_D` columns in `G_tpi` (confirmed: `G['rem_cb_F']['cb_buy_D']` exists;
   they have no `shock_def_D` column, since nothing but CB purchases moves them).
   Reporting the discounted sum of `rem_cb_F × size_F × p` next to the P&L table
   makes the two objects visibly different quantities rather than two rows of one
   table.
3. Mirror the split in `experiments/e1_backstop_schedule.py`'s A5-1 reporting and in
   `docs/paper_draft_results.md`'s German-side captions.

**Paper-side rule this implies:** no sentence may net, sum, or trade off the realised
German transfer against the expected loss. A burden-sharing claim must name which of
the two it quantifies. This is the same class of discipline as CLAUDE.md's existing
ban on "x% fundamental / y% non-fundamental".

**Note on scope.** Making the credit loss *actually flow* through the conduit is a
different proposition: it means `writeoff_enabled = 1`, which is the E3
`e3a_realised_writeoff` arm and a change of the paper's S-1 framing, not a fix. Do not
do it as a side effect of R-2.

---

## R-3 (F-4) — Retire or repair the three dangling `diagnostics/` scripts

| file | state | recommended |
|---|---|---|
| `diagnostics/psilam_breakdown_sweep.py` | raises at line 70 | **delete** — it swept `psi_spread`'s linearity in `psi_lambda_B`, and that relationship no longer exists. Its finding (`PSILAM_BREAKDOWN = 15.0`) is already recorded in `diagnostics/regimes/regime_model.py:41-70`. |
| `diagnostics/psilam_moment_sweep.py` | raises at line 59 | **delete** — same reason; it tuned `psi_lambda_B` to the 150bp moment, and CLAUDE.md already records that tuning history as void. |
| `diagnostics/solve_configs.py` | **does not raise**; silently degenerate | **repair or delete.** Lines 167-173 must stop writing `psi_spread_D/F` into `ss0.toplevel`. With `psi_lambda_B = 0` live, the "G0" arm is now identical to its baseline, so the script's central comparison is vacuous and *looks* like it ran. |

Deleting is preferable to repairing all three: they exist to characterise a dial that
is now fixed at 0, and git history keeps them.

**Complementary hardening.** `code/test_nkpc_blocks.py`'s
`test_no_ad_hoc_sovereign_spread_wedge_anywhere` scans `code/*.py` only, which is
structurally why none of these were caught. Extending its glob to `diagnostics/*.py`
and `experiments/*.py` — with the docstring/comment exemption it already implements —
would close the gap. If that is too broad, a narrower rule catching *writes* of banned
names into a steady-state dict would have caught `solve_configs.py:169` specifically.

---

## R-4 (F-5) — Unstale the figure caption

`code/tpi_plots.py:243` hardcodes `[δ_b = 0.10 → insensitive to q_b_D]`. Derive it
from `ss_final['delta_b_D']` / `['delta_b_F']` as the surrounding captions already do
(`code/tpi_plots.py` was reworked for derived captions in commit `231327c`, and this
line was missed). Cosmetic; no computed number depends on it.

---

## R-5 (F-6) — De-duplicate `cb_pnl`

`experiments/e1_backstop_schedule.py:55-99` reimplements `code/tpi.py:272-296`. The
two currently agree, and E1's copy is actually the better one (it uses each leg's own
`delta_b`, with a comment recording that an earlier draft using D's duration on both
legs contaminated `carry_ss_pv`).

`experiments/` deliberately does not import from `code/`, so the fix is not a plain
import. Options, in order of preference:

1. move the function into `experiments/common.py` and have `code/tpi.py` compute its
   P&L by calling into a shared, dependency-free module both can import (a small
   `code/cb_accounting.py` importable from both is the cleanest — it introduces no
   `code/` → `experiments/` dependency in either direction);
2. failing that, add a test asserting the two implementations agree on a fixed
   synthetic IRF, so drift fails loudly.

Leaving it as two copies is the option that produced the `audit_artifacts/` failure.

---

## R-6 (F-7) — Document what the "central bank" is, in `docs/SPEC.md`

No code change. Two modelling choices are currently implicit in `code/tpi.py`'s
comments and should be stated where the paper's modelling choices live:

* **Full per-period pass-through, no CB capital.** The entire net cash flow is
  remitted each period by capital key. There is no retained-earnings buffer, so a loss
  hits the two treasuries in the period it occurs. A real Eurosystem NCB absorbs
  losses against capital and provisions first.
* **No reserve liability, no policy rate.** Purchases are funded by a same-period
  capital call on the treasuries rather than by creating remunerated reserves. This is
  internally consistent — the model has no policy rate, by the design recorded in
  CLAUDE.md's "No policy rate" note — but the seigniorage and
  reserve-remuneration legs of an actual APP/TPI are absent by construction, and a
  reader will expect them.

Both belong next to the existing "Nominal deposits against REAL sovereign bonds"
entry, which is the same kind of declared asymmetry.

---

## Not recommended

* **Adding a CB balance sheet with capital and a capital key** — the brief flags this
  as a possible prerequisite to any creditor-side claim. It is not needed here: the
  capital key and both remittance legs already exist and are correctly signed and
  scaled (F-1). Adding CB *capital* (a retained-earnings stock absorbing losses before
  they reach treasuries) would be a modelling build, not a patch — a new state, a
  loss-allocation rule, and a recapitalisation trigger — and it is not required by
  anything this audit found.
* **Tightening the pole scan's grid** — see R-1a. It reproduces the same failure at a
  smaller scale.
* **Touching the four CB equations** — `domestic_bond_clearing_tpi`,
  `budget_residual_D_tpi`, `budget_residual_F_tpi`, `external_account_D_tpi` are
  correct, mark exclusively at endogenous `q_b`, match `bond_return_D`'s payoff term
  for term, and are exactly SS-neutral. Nothing in this audit calls for a change to
  any of them.
