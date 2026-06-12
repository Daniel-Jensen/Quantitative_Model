# HHBANK v12 — Post-Audit Verification Report

**Date:** 2026-06-11
**Scope:** verify fixes implemented + correct + regression-free; assess calibration. Not a re-audit.
**Codebase verified:** branch `audit` (396cbd9 + 4c810e1), working tree. `main` is intentionally pre-fix.

---

## Finding status table

| Finding | Audit status | Current status | Evidence | Confidence |
|---------|-------------|----------------|----------|------------|
| **W-1** capital timing | Fixed (proven by repair) | **Fixed** | `labor_*` use `K_t`; `capital_adj_*` mpk use `K_t`; `capital_producer_profit_*` add `mpk*(K−K(−1))`. ca_res_D ≤ 5.8e−8, goods_mkt_D ≤ 1e−16 (run4). | High |
| **W-2** F-bank ToT revaluation | Fixed | **Fixed** | `bank_return_F` L265–266: `rb_F_fg=(1+rb_actual_F)*p(-1)/p−1`, `rb_D_fg` likewise; both feed `rn_F`. goods_mkt_F ≤ 8e−10 (run4). | High |
| **W-3** FOC p-conversion | Fixed | **Fixed** | `divert_bond_foc_F` L344: `rb_F_fg_next=(1+rb_actual_F(+1))*p/p(+1)−1`. Cross FOCs in `divert_portfolio_adj` already converted. No other D-return-vs-F-rate mismatches (swept rn/FOC/Bellman). | High |
| **T-2** deposit timing | Fixed (needs phi_lamb retune) | **Fixed** | `deposit_return_*` `Rgross=(1+rdep(-1))*P(-1)/P`; `bank_return_*` rn funding leg uses `rdep(-1)`; P1 Bellman & FOCs use ex-ante `rdep`. Default shock: n_inter_D(0)=−3.5%, Y_D(0)=−2.5e−4 (both **fall** — sign correct). Stationary at phi_lamb=0.15. | High |
| **A-2** m vs Phi/T | Fixed | **Fixed** | `smart_steady_*` `m = n*(1−(1−f)(1+rn))` (no +Phi+T). P2 residual=8.9e−16 at SS (run4). Critical for any chi1≠0 calibration (e.g. bank-cal). | High |
| **TPI-1** CB budget | Fixed | **Fixed** | Notebook cells TPI-1/2: `budget_residual_D_tpi` with `rem_cb_D`; `domestic_bond_clearing_tpi`. ca_res_D γ-invariant ≈5e−8 across γ∈{0,2,5,10} (tpi_log2). | High |
| **C-1** Δ_cross > 1 | Open (author decision) | **Not Fixed** | Notebook cell 12 `_ic_delta` with `ratio_D=ratio_F=2.0` unchanged → Δ_bF_D=Δ_bD_F=1.4545. Still degenerate. | High |
| **S-1** no realized losses | Open (author decision) | **Not Fixed** | Calibration: `writeoff_enabled_D/F=0.0`. recovery/zeta irrelevant. Default still zero cash-flow impact — pure risk-premium loop. | High |

**Verdict:** 6/6 structural/accounting fixes (W-1, W-2, W-3, T-2, A-2, TPI-1) implemented and verified correct. 2 economic-design findings (C-1, S-1) remain open — both explicitly flagged "author decision" in the audit, neither is an accounting bug.

---

## Correctness notes (did fixes introduce inconsistencies?)

- **W-1 rework is self-consistent.** Under `Y=F(K_t)`: banks earn `mpk·K(−1)` via `rk`; capital producer earns `mpk·(K−K(−1))`; sum = `mpk·K = αY` exhausts output. Term vanishes at SS so SS unchanged (verified: same SS values as pre-W-1). The `mpk_D` dependency added to `capital_producer_profit_D/F` signature wires correctly (mpk is an output of `capital_adj_*`).
- **T-2 applied consistently, not over-applied.** SS-only blocks (`smart_steady_*`, `steady_auxilliary_*`) correctly retain plain `rdep` (constant at SS). Dynamic funding legs use `rdep(-1)`; forward Bellman/FOC use ex-ante `rdep` (was `rdep(+1)`). No double-dating.
- **T-2 side-effect (not a regression, a revealed property):** the old state-contingent deposit windfall was an accidental stabilizer. Removing it makes the genuine debt→spread→debt loop locally explosive at the old phi_lamb=0.02. Minimal stable phi_lamb=0.15 at main's amplification (psi_lambda_B=3.0, def_scale=0.25). See bank_cal_review.md — at lower amplification the stable region may include empirically-plausible phi_lamb.
- **No regression in accounting identities:** goods_mkt_D ≤1e−16, deposit_mkt_D/F ≤4e−15, ca_res_D ≤5.8e−8, goods_mkt_F ≤8e−10, global_goods_res ≤2e−9 (run4_kt_log.txt, both TFP and default shocks). **CA = ΔNFA confirmed** (ca_res_D is the CA−ΔNFA residual by construction; at solver tolerance).

---

## Regression search (Phase 5)

| Identity | Pre-fix (main) | Post-fix (audit) | Status |
|----------|---------------|------------------|--------|
| goods_mkt_D (targeted) | ~1e−16 | ~1e−16 | clean |
| goods_mkt_F (untargeted) | **2.0e−2** (TFP) | **8.0e−10** | fixed, no regression |
| ca_res_D = CA−ΔNFA (untargeted) | **1.5e−4** (TFP) | **3.5e−8** | fixed, no regression |
| deposit_mkt_D/F | ~1e−15 | ~1e−15 | clean |
| Jacobian conditioning H_U cond | 8.0e5 | 5.7e7 (phi_lamb=0.02) → finite & stationary at 0.15 | acceptable post-retune |

No new leaks introduced. The condition-number rise is the genuine doom loop becoming visible (previously masked by T-2 windfall), resolved by phi_lamb.

---

## Executive summary

1. **Unresolved audit findings:** C-1 (Δ_cross=1.45>1, degenerate multi-asset IC) and S-1 (writeoff_enabled=0 → default produces no realized losses). Both are economic-design choices the audit flagged for the author, not accounting bugs. Everything accounting/timing-related is fixed.
2. **Correctly implemented:** W-1 (author's K_t variant), W-2, W-3, T-2, A-2, TPI-1 — all verified statically and dynamically.
3. **Is the model internally consistent now?** **Yes** — all market-clearing and stock-flow identities hold to solver tolerance (≤6e−8), CA=ΔNFA holds, doom-loop signs are economically correct. Two *modelling* limitations remain (C-1, S-1) that affect interpretation, not consistency.
4. **Merge bank-cal?** **No — do not merge the branch.** It is 19 commits behind main and predates every structural fix; merging would revert W-1/W-2/W-3/T-2/TPI-1. Instead **port its calibration values** onto `audit` (see bank_cal_review.md), then re-verify stability. bank-cal is materially better on duration (δ_b≈0.037 vs 0.10), GK exit (f=0.03), EBA exposures, and avoids C-1 (hardcoded Δ=0.2/0.4); but its phi_lamb=0.03 was tuned on the buggy model.
5. **Next highest-priority task:** resolve the **phi_lamb tension**. The empirically-defensible Bohn coefficient (0.10–0.15 annual → phi_lamb≈0.03) was tuned for the pre-fix model. The fixed model needs phi_lamb=0.15 (Bohn=0.6, aggressive) *at main's amplification*. Determine whether porting bank-cal's lower amplification (psi_lambda_B=0, def_scale=0.05, recovery=0.40) lets the fixed model run stably at empirically-plausible phi_lamb. (Probe in progress — `bankcal_stability_test.py`; result in bank_cal_review.md §Stability.)
