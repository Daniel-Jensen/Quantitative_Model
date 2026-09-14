# Central Bank block audit — VERDICT

**Date:** 2026-08-19 · **Branch:** `gk-structural-foc` · **HEAD:** `ea23e94` **plus
uncommitted working-tree changes** to ten files under `code/` (773 insertions / 239
deletions). The audit is of the working tree.

**Method:** static trace of every CB code path plus a live pipeline solve
(`diagnostics/cb_audit/probe_pipeline.py` — `main.py` stages 1–5 and `run_tpi()`,
figures skipped), a stability/Prony analysis (`probe_stability.py`), and a modal
decomposition of the closed-loop operator. No model source was edited. Full evidence
in `run_log.md`; proposed fixes in `recommended_fix.md`, none implemented.

---

## BASE AUDITED — the well-posed base. No inadmissibility caveat applies.

The brief anticipates that this might sit on the structurally inadmissible EBA
concentration regime (`phi_bD ≈ 2.39` against a well-posedness ceiling of 0.336) and
asks for a prominent flag if so. **It does not.** `EBA_CALIBRATION = True` and
`BANK_SCOPE = "broad"`:

| | D | F | required |
|---|---|---|---|
| `lambda_gk` (GK incentive-constraint multiplier) | **+2.108746** | **+0.641387** | > 0 ✅ |
| `Omega` (banker's marginal value of net worth) | **+10.346180** | **+4.037867** | > 0 ✅ |
| `nu_K` (marginal value of capital per unit net worth) | **+0.103205** | **+0.040278** | > 0 ✅ |
| `phi_own` (own-sovereign book / net worth) | **0.452489** | 0.294078 | < 0.8232 ✅ |
| GK feasibility margin | **+0.296567** | **+0.591774** | > 0 ✅ |

`phi_bD_D = 0.4525`, not 2.39. The 2.39 figure is the retired `BANK_SCOPE="ct1"`
(stress-test-sample) reading that GK-2 replaced on 2026-07-31; at the live broad scope
the own-sovereign concentration sits at 55% of its well-posedness ceiling. All four GK
portfolio FOC legs verify to ≤ 2.1e−13, `assert_gk_well_posed` passes, and the solved
SS reproduces the CLAUDE.md calibration row exactly (`q_b_D = 0.974906`,
`q_b_F = 0.965974`, both yields 80.00 bp, SS spread −0.0000 bp, `EL_load_D = 0.7014`,
peak spread 205.87 bp, `n_inter_D[0] = −11.4073%`, `Y_D[0] = −1.9742%`).

**This audit's results transfer.** There is nothing to revert first.

### One prior correction, because two audit steps are framed on it

The brief describes the refactor as *"`EL_price` is now integrated into the GK
incentive constraint itself"*. That is not what happened. `EL_price_D/F` was
**deleted**, and expected default loss now enters through `rb_exp_D/F` — the
state-contingent expected payoff from `bond_return_D/F` — read by
`intermediation_P1_D/F`, the banker's **Bellman envelope**, not by the IC. The IC does
carry a default channel, but it is a different one: `collateral_quality_D`'s bounded
pledgeability map, whose slope `psi_lambda_B_D/F` is **0** at the preferred baseline,
so that channel is switched off and `Delta_eff ≡ Delta = 0.20` on all four legs. The
brief also anticipates no CB balance sheet and no capital key; both exist
(`kappa_cb_F = 0.929`).

---

## Findings, most severe first

### F-1 — HIGH (interpretation, not a code defect). What Germany actually pays is not the number the P&L table reports it bears.

Answering the brief's decisive question directly: **CB income does NOT remit back to
the Greek treasury, the intervention is NOT self-financing, and a creditor side IS
representable and correctly signed.** The conduit is real:

* **(a) a CB balance sheet holding `cb_buy_D`** — yes. `budget_residual_D_tpi`
  (`code/tpi.py:62-78`) computes `cb_flow_D`, the full net cash flow on the CB's
  D-bond book, in the same coupon/survival form as `bond_return_D`.
* **(b) a capital-key split of CB P&L** — yes. `kappa_cb_F = 0.929` to F,
  `1 − kappa_cb_F = 0.071` to D.
* **(c) remittances to *both* treasuries** — yes. `rem_cb_D` enters `b_gov_res_D`,
  `rem_cb_F` enters `b_gov_res_F`.

Measured aggregate split at impact (F converted at `p` and scaled by
`size_F = 11.6967`): F/D = **13.059** against the capital key's 13.085, the 0.2% gap
being the endogenous terms of trade. Both remittances are **negative on impact** — a
capital call, funded by each treasury through its own fiscal rule at its own sovereign
terms. Downstream, German cumulative taxes rise monotonically with γ
(Σ₁₀₀ `TAX_F`: −1.4470e−03 → −1.2724e−03 from γ=0 to γ=10) and German welfare falls
monotonically (ΔW_F = −0.0195 / −0.0482 / −0.0944 at γ = 2/5/10).

**The finding is what sits on top of that.** `writeoff_enabled_D = 0` (the S-1
framing), so `haircut_mult_D = 0` in `cb_flow_D` and **no default loss ever flows
through the conduit**. On the branch the IRF traces, the German treasury books a pure
gain: at γ=10 it receives `kappa_cb_F × prem_pv = 0.1211%` of quarterly SS `Y_D` plus
its share of `mtm_pv = 0.0869%`, and pays nothing for credit risk. Meanwhile
`run_tpi` prints "F bears EL PV = 0.2513% Y_D" — an **off-path** expectation computed
by hand in `cb_pnl` from `EL_load_D × def_rate_D`, which appears in no budget
constraint anywhere in the model.

Both objects are correct and both are needed. They are not commensurable, and the
paper must never net them, sum them, or present the realised transfer and the expected
loss in one column. Concretely: the loading 0.52 / 0.50 / 0.48 says the CB earns ~50
cents of premium per euro of *expected* loss absorbed — under-compensated in
expectation — while the traced path shows the same position ending in profit. Any
German-tail or burden-sharing sentence has to say which of the two it is quantifying.

### F-2 — MEDIUM. The closed-loop pole guard is unsound, and the pole it documents does not exist.

`code/tpi.py:332-353` locates a "closed-loop pole" by scanning
`np.linspace(0.25, 60.0, 240)` for `cond(I − γ·A_cb) > 1e4`, reports γ = 26.50, and
caps the effectiveness curve at `0.75 × 26.50 = 19.88`. `CLAUDE.md` and `code/tpi.py`
record this as "γ ~ 27.3 on the post-GK-refactor calibration".
`diagnostics/regimes/lottery_math.closed_loop_pole` has the same construction.

Two things are wrong.

**The scan misses a nearer singularity.** `A_cb` has a real eigenvalue +0.452155, so
`I − γ·A_cb` is exactly singular at **γ = 2.2116** — inside the intended range. The
0.25-step grid evaluates 2.00 (cond 6.5e2) and 2.25 (cond 4.5e3), both under
threshold, and steps over it; `sign(det)` flips +1 → −1 across the gap. This is
exactly the failure the guard's own comment warns about.

**But neither pole is economic — both are T=500 terminal-truncation artefacts.** The
resonant eigenvector carries **0.0000 of its mass in the first 100 quarters and 0.9922
in t = 400–499**, peaking at index 499. The terminal columns of `A_cb` are
pathological while every interior column is clean (`||A[:,499]|| = 3.86` vs ~0.0065
for interior columns; `A[499,499] = +1.080`, the only positive diagonal and the only
entry exceeding 1 in magnitude). Dropping the last five rows and columns collapses the
worst conditioning over γ ∈ [0.1, 30] from 5.45e+05 to 8.51e+01, removes every pole
below γ = 36, and **changes the reported peak spread by nothing at γ = 2, 5, 10** and
by 0.03 bp at γ = 20.

**No reported number is wrong** — every TPI statistic is computed on `[:100]`, where
the artefact has no mass. What is wrong is the guard (it would equally miss a genuine
pole falling between grid points), the γ = 19.88 cap (imposed for a spurious reason),
and the documented claim of a pole at γ ≈ 27.3.

### F-3 — MEDIUM (structural interpretation). The bonds the CB buys come off *German* balance sheets, not Greek ones.

Differentiating the clearing identity `b_D_D = b_gov_D − size_F·b_D_F − cb_buy_D`:

| horizon | `d b_D_D / d cb_buy` | `d b_gov_D / d cb_buy` | `size_F · d b_D_F / d cb_buy` |
|---|---|---|---|
| 0 | −0.1268 | +0.0325 | **−0.8408** |
| 1 | −0.0671 | −0.0539 | **−0.9867** |
| 12 | −0.0253 | −0.0202 | **−0.9949** |

84% of a purchase at impact and ~99% from t=1 onward is absorbed by German banks
shedding their Greek book — governed by `psi_bD_F = 0.5`, the cross-border portfolio
adjustment cost in `gk_cross_border_foc` — not by Greek banks. Only 13% at impact and
2–6% thereafter comes off Greek balance sheets.

The spread still compresses, and the reason is worth stating precisely rather than as
"the CB relieves Greek banks":

```
 d b_D_D    /d cb_buy [0,0] = -0.1268     modest quantity relief on D banks
 d theta_D  /d cb_buy [0,0] = -1.1777     required leverage falls, IC slackens
 d n_inter_D/d cb_buy [0,0] = +0.5258     D bank net worth recovers (MTM on q_b_D)
 d K_D      /d cb_buy [0,0] = +0.0170     capital crowded back in
 d q_b_D    /d cb_buy [0,0] = +0.0373     bond price up -> spread compresses
```

The channel is **price support → mark-to-market recapitalisation of Greek banks**,
with quantity relief a minor contributor. This is not a defect; it is a result the
paper should own rather than describe as balance-sheet relief.

### F-4 — LOW. Three `diagnostics/` scripts are wired to the deleted `psi_spread`. The CB block is clean.

Answering Step 3 directly: **zero live `psi_spread` references in the CB block, the
clearing condition, or the residual equations.** The only occurrence in `code/tpi.py`
is a comment at line 313 recording that the old `loading ≈ 1 + psi_spread/EL_price`
closed form no longer exists and is deliberately not replaced. Every `psi_spread`
string in `code/` and `experiments/` is prose. `code/test_nkpc_blocks.py`'s AST
scanner over `code/*.py` passes (22/22 tests).

Outside that scanner's reach, three scripts still execute against the deleted symbol:

* `diagnostics/psilam_breakdown_sweep.py:70` and
  `diagnostics/psilam_moment_sweep.py:59` — bare `float(ss["psi_spread_D"])`; these
  **raise** on the current model. Loud, therefore harmless.
* `diagnostics/solve_configs.py:130/167/169/173` — the read at 130 is inside
  `try/except` and degrades to `MISSING`, but line 169 then **writes**
  `ss0.toplevel["psi_spread_D"] = 0.0`, inventing a symbol no block reads, and line
  173 reads it back successfully. The script does not crash; it silently produces a
  "`psi_lambda_B = 0`" arm identical to its own baseline, since `psi_lambda_B` is
  already 0 at the live calibration. **Silent, therefore worse.**

### F-5 — LOW (cosmetic). Stale hardcoded parameter in a CB figure caption.

`code/tpi_plots.py:243` hardcodes the subtitle `[δ_b = 0.10 → insensitive to q_b_D]`.
Live values are `delta_b_D = 0.0777` and `delta_b_F = 0.0568`. No number is computed
from it.

### F-6 — NOTE. `cb_pnl` is duplicated.

`experiments/e1_backstop_schedule.py:55-99` carries a second copy of `code/tpi.py`'s
`cb_pnl`. It currently agrees (both read `EL_load_D` from the solved model rather than
from an anchor), but duplicated model algebra is the exact drift mechanism that made
the retired `audit_artifacts/` harness test a different model for weeks.

### F-7 — NOTE. What the "central bank" is, stated plainly for the paper.

It is a **capital-key conduit with full per-period pass-through**, not a central bank
with a balance sheet in the institutional sense:

* **no CB capital or retained earnings** — the entire net cash flow is remitted every
  period, so there is no loss-absorption buffer and no delay between a loss and its
  fiscal incidence;
* **no reserve liability and no policy rate** — the asset purchase is funded by a
  same-period capital call on the two treasuries, not by creating remunerated
  reserves. This is consistent with the model having no policy rate at all (the
  union-inflation normalisation is the `phi_pi → ∞` limit), but it means the
  seigniorage/reserve-remuneration leg of a real APP/TPI is absent by construction.

Both are defensible modelling choices and both should be stated as such.

---

## Step-by-step results

| Step | Question | Result |
|---|---|---|
| **0** | Which base? | **Well-posed.** EBA + `BANK_SCOPE="broad"`; `lambda_gk`, `Omega`, `nu_K` all > 0 in both countries; `phi_bD_D = 0.4525` with feasibility margin +0.2966. No caveat. |
| **1** | Where do CB coupons and P&L go? | **Both treasuries, by capital key 0.929/0.071**, measured ratio 13.059 vs 13.085 implied. Not self-financing; creditor side representable and correctly signed (German taxes up, German welfare down, monotone in γ). **But** see F-1: with `writeoff_enabled = 0` no credit loss flows through the conduit, so the realised transfer and the reported expected loss are different objects. |
| **2** | Does the CB mark to the endogenous `q_b`? | **Yes, everywhere.** All four CB equations use endogenous `q_b_D` / `q_b_D(-1)`, and `cb_flow_D` matches `bond_return_D`'s payoff term for term including `zeta_writeoff_D` and `writeoff_enabled_D`. Established by diff as well as inspection: the refactor did not touch these four blocks at all. The SS-price sites in `cb_pnl` are the documented linearisation convention (`cb_buy_ss = 0`, so `q_b_D_ss × cb` is the correct first-order value and `mtm_pv` restores the revaluation separately), not stale prices. **No stale-price reference found.** |
| **3** | Dangling `psi_spread`? | **None in the CB block, clearing, or residual equations.** Three `diagnostics/` scripts are dangling (F-4). CB transmission confirmed to run through `intermediation_IC_D`'s `phi_bD_D` and `k_balance_sheet_D` and through nothing else — no additive spread term anywhere on the path. |
| **4** | Feedback sign and stability | **`A_cb[0,0] = −4.397e−03 < 0`** — purchases compress. Peak spread monotone: 205.87 → 193.25 → 176.75 → 154.36 bp at γ = 0/2/5/10. Prony moduli all < 1 (spread 0.940 → 0.947); `b_gov_D[499]` ≤ 3.7e−04; Walras residuals stay clean and `max|ca_res_D|` actually *falls* with γ. **No breakdown region in [0, 30]** once the terminal-truncation artefact is removed (F-2). |
| **5** | SS neutrality | **Exact.** `cb_buy_ss = 0`; `cb_flow_D`, `rem_cb_D`, `rem_cb_F` all identically 0 at SS; every TPI block output **bit-identical** (difference 0.000e+00) to its non-TPI counterpart — `b_gov_res_D`, `nfa_D`, `ca_res_D`, `b_D_D`, `b_F_F`. Confirmed dynamically by `G_tpi[cb=0]` vs baseline `G`: `max|err| = 0.00e+00`. The refactor did not move the SS through the CB block; prior SS-invariance arguments survive. |

---

## What this audit did not do

* No global/nonlinear check. Every statement here is about the linearised model.
* `psi_lambda_B_D/F = 3.01` (the diagnostic arm) was not audited; at the preferred
  baseline the pledgeability channel is off, and the CB's interaction with it when
  switched on is untested.
* No sensitivity of F-3's 84/99% split to `psi_bD_F`; the number is reported at the
  calibrated `psi_bD_F = 0.5` only.
* The Prony estimator that `docs/STATE.md:2317` refers to is not in the working tree
  (it went with `audit_artifacts/`); a self-contained order-selected reimplementation
  with a passing synthetic self-test is at `diagnostics/cb_audit/prony.py`.
