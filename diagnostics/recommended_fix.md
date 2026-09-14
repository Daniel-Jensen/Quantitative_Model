# Recommended fix (PROPOSED — not implemented)

**Diagnosis:** sovereign default risk transmits only through the friction
parameter `psi_lambda_B`; there is no fundamental expected-loss channel that
survives `psi_lambda_B → 0` (see [`VERDICT.md`](VERDICT.md)). This is a
**modelling decision**, not a mechanical bug — it touches documented choices
(S-1 write-off, F-1 market-value rule), so it needs author sign-off before
implementation. What follows is a concrete, minimal, Walras-safe proposal.

---

## Core idea

Add a **fundamental expected-loss loading** to the bank's sovereign-bond
optimality condition, so that a higher expected default probability raises the
required bond return and therefore **lowers the equilibrium price `q_b` today** —
independent of `psi_lambda_B`. The resulting mark-to-market fall in `q_b`
propagates to bank net worth through the *existing* (endogenous) plumbing
(`k_balance_sheet`, `bank_return`), tightening the IC. Because it changes only a
**required return / equilibrium price** — not any realized cash flow — it books
**no** fiscal relief and no realized bank loss, so it (a) preserves the
risk-premium framing S-1/F-1 require and (b) cannot leak Walras (all agents face
the same `q_b`).

Structurally this is the expected-loss term of the defaultable-bond Euler
equation `q_b = E[Λ·payoff]`; `psi_spread` (the IC-shadow-value premium) then
becomes a pure **amplifier on top of** a fundamental channel — the economically
correct ordering.

## Precise change

Introduce a per-unit expected-loss coefficient `EL_price_D`, `EL_price_F` and add
it, alongside the existing `psi_spread`, to the required spread in the three bond
FOC blocks. It is **zero at the steady state** (it multiplies `def_rate(+1)=0`),
exactly like `psi_spread`.

**1. `code/equations_D.py` — `divert_bond_foc_D` (lines 426–437).** Current:

```python
req_spread = excess_return_bD_D_ss + psi_spread_D * def_rate_D(+1)
```

Proposed (add `EL_price_D` to the signature and the loading):

```python
req_spread = excess_return_bD_D_ss + (EL_price_D + psi_spread_D) * def_rate_D(+1)
```

**2. `code/equations_F.py` — `divert_bond_foc_F` (lines ~394–399).** Symmetric,
with `EL_price_F` on `def_rate_F(+1)`.

**3. `code/equations_global.py` — `divert_portfolio_adj` (lines 79–96).** The
cross-border holdings price the *issuer's* default risk. D holds F-bonds → load
on `def_rate_F(+1)`; F holds D-bonds → load on `def_rate_D(+1)`:

```python
prem_FD = excess_return_F_D_ss + (EL_price_F + psi_spread_D) * def_rate_F(+1)
prem_DF = excess_return_D_F_ss + (EL_price_D + psi_spread_F) * def_rate_D(+1)
```

**4. `code/steady_state.py` — `_apply_ss_anchors` (lines 27–77).** Compute
`EL_price_D/F` as an SS anchor (mirroring how `psi_spread_*` is set), from the
existing loss primitives, and write it into `ss_in.toplevel` + `cal`:

```python
# expected fractional loss per unit default probability, priced by bondholders
# (independent of writeoff_enabled — this is pricing, not fiscal booking)
h_D   = 1.0 - cal['recovery_rate_D']                      # loss-given-default
q_D   = float(ss_in['q_b_D']); db_D = cal['delta_b_D']; z_D = cal['zeta_writeoff_D']
'EL_price_D': h_D * (db_D + z_D * (1.0 - db_D) * q_D) / q_D,   # ~0.10 (z=0) … ~1.0 (z=1)
```

(and symmetric `EL_price_F`). **Magnitude note:** with the current SS
(`q_b ≈ 0.976`, `delta_b = 0.10`, `recovery = 0`) this is `≈ 0.10` if only the
coupon is at risk (`zeta = 0`) and `≈ 1.0` if the continuation value is too
(`zeta = 1`). This is the natural dial for "how much fundamental repricing per
unit default probability."

**5. `code/calibration.py` — (optional).** If you prefer an explicit, directly
calibrated constant over the derived formula, add `EL_price_D/F` to the
calibration dict (e.g. `0.5`) and skip step 4's formula. Either way, add a master
toggle if you want to switch the fundamental channel off for the nested-model
ladder rather than relying on `psi_lambda_B`.

## Simpler alternative (equivalent in spirit)

If you would rather not add a parameter: **decouple the haircut multiplier** by
splitting `writeoff_enabled` into `price_writeoff_enabled` (default **1**, used
by the *pricing* FOC only) and `writeoff_enabled` (default **0**, used by
`budget_residual_*` for fiscal booking). Implement the pricing side via the FOC
required-spread as above — **not** by re-enabling the haircut inside
`bond_return_*` / `rb_actual`, which would make realized bank returns reflect a
loss the government never books and reintroduce the S-2 asymmetry (a Walras
leak). The FOC route keeps the loss in the *price*, where it is symmetric.

## Steady-state side-effects

**None.** Every added term is proportional to `def_rate(+1)`, which is `0` at the
steady state, so the SS is **bit-for-bit unchanged** — no recalibration of
`beta`, `q_b`, portfolio shares, or `Delta` is required. (Same reason
`psi_spread` is SS-neutral.) This is verifiable by asserting the SS solve and all
`ss_checks` in `summary.json` are identical before/after.

## Required re-verification after implementing (author to run)

1. **Symptom resolved:** at `psi_lambda_B = 0`, the 1pp default shock now
   produces a **nonzero, correctly-signed** response — `q_b_D < 0`, `n_inter_D <
   0`, `spread_rb > 0`, `Y_D < 0` on impact (rerun `diagnostics/solve_configs.py`;
   the ψλ=0 column should no longer be all zeros).
2. **Baseline behaviour:** decide whether `EL_price` is *added on top* of
   `psi_lambda_B = 3` (baseline doom loop strengthens — you may want to
   recalibrate `psi_lambda_B` *down*, since it is now a pure amplifier) or whether
   `psi_lambda_B` is reduced so the baseline peak spread is unchanged.
3. **Walras intact:** `code/main.py` — `goods_mkt_F`,
   `ca_res_D` still within thresholds (the fix moves only a common price, so they
   should be).
4. **Fiscal stability:** re-map `phi_lamb` — `q_b` now responds to `def_rate`
   even at low `psi_lambda_B`, which changes the debt→spread feedback gain.
5. **F-1 interaction (critical):** re-run under `mv_rule = 1`. The market-value
   Bohn rule now sees a `q_b` that falls on *every* default shock (not just via
   the friction). Confirm the response is **not** perverse (spread must *widen*,
   `n_inter` and `Y` must *fall*) — i.e. that the fundamental channel dominates
   the mv-rule's counteracting fiscal loosening. This is the one place the fix
   could interact badly with existing choices; test it explicitly.

## Docs to update once accepted

`docs/STATE.md` (retire/relabel the "pure risk-premium loop = psi_lambda_B"
framing; add the EL_price channel), `docs/SPEC.md` (default-risk feedback), and
`CLAUDE.md` (Key modelling choices; S-1/F-1 coupling now has a third lever).
