# VERDICT — why `psi_lambda_B = 0` gives null sovereign-risk results

**The cause is Case 3 (no fundamental channel): in the baseline calibration the
*only* live transmission of sovereign default risk runs through the structural
`psi_lambda_B` terms, so zeroing that one parameter leaves `def_rate` entering
no equation and the shock becomes inert.** The single decisive test — the
response of bank net worth `n_inter_D` to a 1pp default shock at
`psi_lambda_B = 0` — is **exactly `0.0`** (vs **−3.52%** at baseline), and it is
zero *together with* the endogenous bond price `q_b_D` (also exactly `0.0`, vs
−2.60e−2 baseline), while the driver `def_rate_D` is unchanged (`+0.01`, ratio
1.000). Case 1 (stale price) is ruled out; Case 2 (conflation) is a **mislabel**
that nonetheless describes the true effect (see below).

Evidence git commit `f25dcf6`; interpreter `ssj` (Python 3.12.12). Raw numbers in
[`summary.json`](summary.json) and the tables in [`run_log.md`](run_log.md);
figures `01_*`–`03_*`.

---

## The three cases, adjudicated

### Case 3 — No fundamental channel independent of `psi_lambda_B`  ✅ **ROOT CAUSE**

Complete map of every `def_rate` consumer in the dynamic model and its gate
(Step 3 in `run_log.md`). In the baseline calibration (`writeoff_enabled = 0`,
`T1 = 0`) the cash-flow (`bond_return`), fiscal (`budget_residual`) and macropru
(`macro_pru_tax`) channels are already dead, leaving **only** three terms — all
scaled by `psi_lambda_B`:

- `intermediation_IC_D/F`: `Delta_b·_eff = Delta_b· + psi_lambda_B·def_rate(+1)`
- `divert_bond_foc_D/F`: `req_spread = … + psi_spread·def_rate(+1)`, `psi_spread ∝ psi_lambda_B`
- `divert_portfolio_adj`: cross-border `prem = … + psi_spread·def_rate(+1)`

Set `psi_lambda_B = 0` (⇒ `psi_spread = 0`) and `def_rate` enters **no**
structural equation. Empirically, on the 1pp default shock at `psi_lambda_B = 0`
**every** endogenous variable is *identically* `0.0000e+00` (not merely small —
the Jacobian column is structurally zero), while `def_rate_D` still moves by the
full `0.01`:

| var | peak\|·\| baseline | peak\|·\| ψλ=0 | ratio |
|-----|------------------:|--------------:|------:|
| def_rate_D | 1.00e−02 | 1.00e−02 | **1.000** |
| q_b_D | 2.60e−02 | 0.00e+00 | 0.000 |
| spread_rb | 3.31e−03 | 0.00e+00 | 0.000 |
| n_inter_D | 3.52e−02 | 0.00e+00 | 0.000 |
| theta_D | 3.41e−02 | 0.00e+00 | 0.000 |
| Y_D | 2.74e−04 | 0.00e+00 | 0.000 |
| C_D | 1.09e−03 | 0.00e+00 | 0.000 |

The fundamental doom-loop `def_rate → q_b → n_inter → IC(θ) → K/credit → Y`
*is* real at baseline (correct signs: `n_inter_D[0] = −3.52%`, `Y_D[0] < 0`,
`spread` widens) but is powered **entirely** by the single friction parameter
`psi_lambda_B`. There is no expected-loss / risk-premium channel that survives
`psi_lambda_B → 0`.

### Case 1 — Dead mark-to-market / stale price  ❌ **RULED OUT**

Two independent lines of evidence:

1. **Code:** `bond_return_D`, `bank_return_D`, `k_balance_sheet_D` all value the
   bond book at the **endogenous current** `q_b` (and `q_b(-1)`), not a
   fixed/SS price (Step 2 in `run_log.md`). The revaluation plumbing is intact.
2. **Decisive empirics:** on the **TFP control** at `psi_lambda_B = 0`, `q_b_D`
   *does* move (`7.68e−3`) and feeds `n_inter_D` (`0.366`, 89% of baseline) — so
   `q_b` is a live endogenous price wired to net worth. A stale-price bug would
   show `q_b` moving on the *default* shock while `n_inter` stayed flat. We see
   the **opposite**: on the default shock `q_b` and `n_inter` are flat
   **together**, and `q_b` is demonstrably alive elsewhere. The break is
   *upstream* of the price (nothing forces `q_b`), not between price and net
   worth.

### Case 2 — Parameter conflation  ⚠️ **MISLABEL — but the effect is real**

`psi_lambda_B` appears **only in structural equations** (the GK incentive
constraint and the bond FOCs); **no macroprudential policy rule references it**
(Step 1 in `run_log.md`). The model's actual macropru instrument is a *separate*
tax `tau_mp = T0 + T1·def_rate` (`macro_pru_tax_D/F`), which is independently off
(`T1 = 0`). So `psi_lambda_B` is **not** "doing double duty as a policy *and*
structural coefficient" — it is purely structural. But because it is the **sole
structural carrier of sovereign-risk sensitivity**, zeroing it removes
fundamental risk sensitivity — exactly the *symptom* Case 2 predicts. The
premise that `psi_lambda` is "the macroprudential policy-response coefficient"
is incorrect; treat it as the structural risk-premium / divertability
sensitivity.

---

## Fourth framing (the deeper "why"): `writeoff_enabled` conflates pricing with booking

The reason no *expected-loss* channel survives is that the switch which turns off
realized default losses also turns off **priced** default losses. In
`bond_return_*` and `budget_residual_*` the haircut is gated by a single
multiplier `haircut_mult = writeoff_enabled`. Turning it off (deliberate — S-1,
and required by F-1's market-value rule to avoid a perverse response) removes
**both**:

- (i) forward-looking **expected-loss pricing** — a higher default probability
  should lower the bond's expected payoff and hence `q_b` *today*; and
- (ii) ex-post **realized fiscal/bank booking** of a default that actually fires.

Only (ii) needs to be off for the risk-premium framing. (i) is basic defaultable-
bond pricing and should operate regardless. With (i) gated off, the *only*
remaining sensitivity of `q_b` to `def_rate` is the collateral-friction premium
`psi_spread ∝ psi_lambda_B` — hence the total dependence on one parameter.

---

## Sanity checks on the diagnostic itself (Step 5)  ✅

- **Shock reaches the driver:** `def_rate_D` impact `= 0.01` in *both* configs —
  a live shock, not a dead one; the null is a dead *channel*.
- **Model not globally broken:** the TFP control transmits normally at
  `psi_lambda_B = 0` (`Y_D` 98.7% of baseline, `q_b_D`, `n_inter_D`, `w_D` all
  move). The null is specific to the sovereign-risk channel.
- **Steady state clears:** `goods_mkt_D = −3.9e−7`, `deposit_mkt_D = 1e−13`,
  `deposit_mkt_F = −2e−15`, `ca_res_D = 0`, IC/P1 block residuals ≤ 1.8e−15;
  `spread_rb(SS) = 0` exactly (equal duration/price, `q_b_D = q_b_F = 0.9757`).
  We are diagnosing the intended steady state.

---

## Consequence for the research programme

This is a genuine blocker for the nested-model ladder (planned Task 3): stripping
the GK/IC block — where `psi_lambda_B` lives — would necessarily zero **all**
sovereign-risk transmission in the non-GK nested models, which is economically
nonsensical and would sink the "justify the GK block" comparison. A fundamental,
`psi_lambda_B`-independent default channel must exist first. Proposed fix (not
implemented): [`recommended_fix.md`](recommended_fix.md).
