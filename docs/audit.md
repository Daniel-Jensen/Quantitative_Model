# HHBANK v12 — Forensic Audit Master Log

**Audit start:** 2026-06-10
**Auditor stance:** Hostile referee. Objective is falsification, not validation.
**Sources of truth:** `equations_D.py`, `equations_F.py`, `equations_global.py`, `model_v12.ipynb`. Documentation treated as evidence only; code prevails on conflict.
**Companion files:** `equation_reconstruction.md` (model as implemented), `walras_forensics.md` (leak derivation + fix proof), `calibration_review.md`, `mechanism_validation.md`. Numerical artifacts in `audit_artifacts/` (`run_audit.py`, `fix_test.py`, `tpi_test.py` + JSON/npz results).

---

# Executive Summary

The model solves, converges, and is *internally inconsistent*. Three independent accounting violations and one mechanism inversion were found, each verified numerically; two of the three violations were additionally proven by repair (patched copies drive the residuals to machine tolerance).

**Highest-risk findings, ranked (severity × confidence × impact):**

1. **T-2 (Critical/High): deposits are state-contingent** — the deposit rate is a period-t unknown paid on the t−1 deposit stock (9× quarterly GDP). On a default shock, the equilibrium fall in rdep transfers ~13% of quarterly GDP from households to banks at impact: **bank net worth and GDP *rise* after a sovereign default shock.** The doom loop runs backwards. (This is the same issue behind the reverted commit "Fix perverse GDP increase under default shock".)
2. **W-2 (Critical/High): F-bank's realized return ignores terms-of-trade revaluation** on its D-good-denominated bond book → Walras's law fails; `goods_mkt_F` leaks up to **2% of GDP** on a 1% TFP shock. All cross-country spillover IRFs (C_F, U_F, Y_F, n_inter_F, p) are wrong by 39–124%.
3. **W-1 (Critical/High): capital-timing off-by-one** (`Y_t = F(K_t)` but banks earn `mpk_t` on `K_{t−1}`) → income on newly installed capital accrues to nobody; `goods_mkt_D + ca_res_D = αY·ΔK/K ≠ 0`; **CA = ΔNFA fails at first order** (leak ≡ predicted gap, corr 0.9999998).
4. **TPI-1 (Critical/High): the TPI policy experiment has no central-bank budget constraint** — `cb_buy_D` removes bonds from private hands with nobody paying for them; the closed-loop rule injects unbacked resources exactly during stress, contaminating the headline welfare-gain figures.
5. **S-1 (Major/High): `writeoff_enabled = 0`** — baseline "default" never touches any cash flow. No haircuts, no recovery, no debt relief. The doom loop is a pure risk-premium loop; any text implying realized losses is unsupported.
6. **C-1 (Major/High): back-solved divertabilities Δ_cross = 1.45 > 1** (cross-border bonds "145% divertable"); the multi-asset IC is degenerate at SS (Σ(1−Δᵢ)φᵢ ≡ 0 forced by the single-asset λ_gk formula).
7. **A-2 (Major/latent): m vs Phi/T bookkeeping** — SS and dynamic banker laws disagree whenever adjustment costs or the macropru tax are nonzero (currently both zero).

**What is sound:** household block (EGM/GHH budget exact), sovereign cash-flow accounting (bond return ↔ government budget mutually consistent, no double-counting), trade/CES/NX identities, NFA valuation, deposit and bond market clearing, perpetuity yield formula, inner GK Bellman/IC solve, and the *choice* of which equations to drop from the target system. The Walras failures come from two budget-constraint elements, not from the external block where one would first look.

**Ranked improvement roadmap:**

| # | Action | Effort | Effect | Status (2026-06-11) |
|---|---|---|---|---|
| 1 | W-1 fix (author variant: keep `Y=F(K_t)`, add `mpk·ΔK` to `cap_profit`) | 6 lines | restores CA=ΔNFA | **APPLIED + verified** (ca_res_D ≤ 6e−8; commit 4c810e1) |
| 2 | W-2/W-3 fix: p-conversions in `bank_return_F`, `divert_bond_foc_F` | ~6 lines | restores goods_mkt_F; corrects all spillovers | **APPLIED + verified** (goods_mkt_F ≤ 1e−9) |
| 3 | T-2 fix: re-date deposit rate (`rdep(-1)` on old deposits) | ~6 lines + re-anchor | restores doom-loop sign | **APPLIED + verified**; required phi_lamb 0.02→0.15 (see T-2 follow-on) |
| 4 | TPI: add CB budget (remit purchases/coupons through `budget_residual_D`) | 1 block | makes welfare figures meaningful | **APPLIED + verified** (ca_res γ-invariant at 5e−8); notebook cells TPI-1/2 patched |
| 5 | Re-run **all** IRF/TPI/welfare results after 1–4 | compute only | paper tables/figures | **REQUIRED — notebook must be re-executed**; all stored figures stale. Post-fix TPI: ΔW_D +1.88 / ΔW_F −1.90 at γ=10 (near zero-sum burden shift) |
| 6 | S-1: decide loss regime (writeoff=1, ζ>0, δ_b≈0.04) or re-frame paper as pure risk-premium model | calibration | mechanism honesty | open (author decision) |
| 7 | C-1: multi-asset-consistent λ_gk or drop Δ apparatus | derivation | removes Δ>1 embarrassment | open (author decision) |
| 8 | A-2: align m/Phi/T between SS and dynamics | 3 lines | needed before any chi1/T0/T1 ≠ 0 experiment | **APPLIED** (`+ Phi + T` removed from m in `smart_steady_*`) |
| 9 | X-0/X-1: notebook execution-order hardening, delete dead code/params | hygiene | reproducibility | open (dead blocks still imported by notebook cell 7 — delete jointly) |

Acceptance test for 1–4 (**passed**): max |ca_res_D|, |goods_mkt_F| ≤ 1e−7 on all shocks incl. TPI, n_inter_D and Y_D fall on a default shock, system stationary (b_gov t499 ≈ 2e−5).

---

# Open Hypotheses

| ID | Hypothesis | Status | Notes |
|----|------------|--------|-------|
| H1 | Walras-law leak in external account block | REVISED | External block itself algebraically correct; leak originates in production/capital timing (W-1) and F-bank return units (W-2). ca_res_D is the symptom, not the cause. |
| H2 | D-good/F-good conversion inconsistency (p conversions) | CONFIRMED | W-2 (bank_return_F), W-3 (divert_bond_foc_F). |
| H3 | Timing mismatch in sovereign bond returns | REJECTED | def_rate_t from b_gov_{t−1}; coupons/issuance/pricing internally consistent; perpetuity yield formula correct. New timing problem found in *capital*, not bonds (W-1). |
| H4 | Double-counting of default losses | REJECTED (baseline) / UNTESTED (writeoff=1) | With writeoff_enabled=0 there are no realized losses to double-count (S-1). intermediation_P2_F comment shows writedown terms were already removed once. If writeoff_enabled=1 is ever used, m/Phi/T bookkeeping (A-2) and budget/bank symmetry must be re-audited. |
| H5 | CA ≠ ΔNFA in transition | CONFIRMED (analytical; numeric pending) | goods_mkt_D + ca_res_D = Y − wN − mpk·K(−1) ≠ 0 off SS. See W-1, walras_forensics.md. |
| H6 | Incorrect valuation of foreign bond positions | PARTIAL | NFA valuation itself is market-value consistent; the *return* on F-bank's D-good positions misses ToT revaluation (W-2). |
| H7 | Calibration suppresses doom loop | OPEN | def_scale elasticity, psi_bF_D=0.5 level penalties, writeoff=0 all candidates; quantify in mechanism runs. |
| H8 | Convergence masks structural errors | CONFIRMED | 23 targets imposed ⇒ solver happily produces paths violating the 2 untargeted identities; convergence proves nothing about W-1/W-2. |
| H9 | Deposit/bank funding timing off-by-one | REJECTED | rn uses t−1 shares with t returns; deposit return paid on Dep_{t−1}; consistent. |
| H10 | GHH labour aggregation inconsistency | REJECTED | Uniform N across households; w/P=φN^{1/ν} consistent with firm FOC via 2 unknowns (w,N), 2 targets. |

---

# Findings Log

## W-1 — Capital-timing mismatch destroys Walras's law (CA ≠ ΔNFA at first order)

### Severity
Critical

### Confidence
High (analytical derivation; numerical confirmation in run1_results.json pending)

### Location
`equations_D.py`: `labor_D` (line ~233), `capital_adj_D` (lines ~217–223), `bank_return_D`; F analogues in `equations_F.py`. Symptom appears in `ca_res_D` (`equations_global.py: external_account_D`) and `goods_mkt_F`.

### Description
Production uses **current** capital, `Y_t = Z·K_t^α·N^{1−α}`, `mpk_t = αY_t/K_t`, but banks (capital owners) earn `rk_t = (mpk_t + (1−δ)Q_t)/Q_{t−1} − 1` on capital purchased at **t−1** (`bank_return_D` uses lagged shares). The capital producer receives only `Q·ΔK − I`. Summing every implemented budget constraint in country D yields the exact identity
`goods_mkt_D + ca_res_D = Y_t − w_t N_t − mpk_t K_{t−1} = αY·(K_t − K_{t−1})/K_t`.
The income from the within-period marginal product of newly installed capital accrues to **no agent**. Since `goods_mkt_D` is a solver target (=0), the entire gap is dumped into the *untargeted* `ca_res_D`. CA = ΔNFA fails at first order in any IRF with investment dynamics; `goods_mkt_F` inherits the F-country analogue.

### Evidence
Full derivation in `walras_forensics.md` §2. Standard SSJ examples (KS, two-asset) use `K(−1)` in production precisely to avoid this.

### Proposed Fix
Two equivalent closures:
(a) standard SSJ timing — `K(-1)` in `labor_*` and `capital_adj_*` mpk; or
(b) keep `Y = F(K_t)` and pay the within-period product of new capital to the capital producer: `cap_profit += mpk*(K - K(-1))` (zero at SS).

### Resolution (author decision, 2026-06-11)
Author chose (b): production keeps the current-K_t convention. Implemented on `audit` branch commit 4c810e1: `capital_producer_profit_*` gains `mpk*(K - K(-1))`; `labor_*`/`capital_adj_*` keep `K_t`. Verified equivalent on the Walras criterion: max |ca_res_D| ≤ 5.8e−8, |goods_mkt_F| ≤ 8e−10, stationary, default-shock signs correct (run4_kt_log.txt). Note: under (b) `rk` pays banks the marginal product of *current* K on capital bought at t−1, and the capital producer earns rents on installation — fine for accounting, but state in the paper's model section.

### Status
Confirmed (analytical) — numeric verification pending

---

## W-2 — F-bank realized portfolio return ignores terms-of-trade revaluation

### Severity
Critical (stock-flow consistency of F-bank; second independent Walras leak)

### Confidence
High

### Location
`equations_F.py`: `bank_return_F` (lines ~252–262)

### Description
Both sovereigns' bonds are D-good claims. F-bank accounts are in F-goods. A bond position worth `q(−1)b(−1)/p(−1)` F-goods at t−1 pays `(1+rb_t)q(−1)b(−1)/p_t` F-goods at t, i.e. F-good return `(1+rb_t)·p_{t−1}/p_t − 1`. `bank_return_F` instead applies the raw D-good return `rb_actual` to F-good shares converted at `p(−1)`: the capital gain/loss from Δp on the entire bond book (φ_bF + φ_bD ≈ 0.40 of net worth) vanishes from `n_inter_F` evolution. This is precisely the cross-border valuation channel the model is supposed to study.

### Evidence
Unit table in `equation_reconstruction.md` §0. Contrast with `divert_portfolio_adj`/`domestic_bond_foc_F`, which do convert expected returns with `p/p(+1)`.

### Proposed Fix
In `bank_return_F`:
```python
rb_F_fg = (1 + rb_actual_F) * p(-1) / p - 1
rb_D_fg = (1 + rb_actual_D) * p(-1) / p - 1
rn_F = (kappa_lag_F*(rk_F - rdep_F) + phi_bF_lag_F*(rb_F_fg - rdep_F)
        + phi_bD_lag_F*(rb_D_fg - rdep_F) + rdep_F)
```
No SS change (p constant at SS).

### Status
Confirmed (analytical) — numeric decomposition pending

---

## W-3 — F-bank own-bond FOC drops the p-conversion its predecessor had

### Severity
Major

### Confidence
High (that it is inconsistent); Medium (that it is unintended — but the superseded block `domestic_bond_foc_F` contains the conversion, strongly suggesting it was lost in refactoring)

### Location
`equations_F.py`: `divert_bond_foc_F` (lines ~374–384) vs dead `domestic_bond_foc_F` (lines ~329–337)

### Description
`divert_bond_foc_F` prices q_b_F from `(rb_actual_F(+1) − rdep_F(+1)) = spread`: a D-good return compared to an F-good rate, no `p/p(+1)` conversion. The same bank's FOC for cross-holdings (`divert_portfolio_adj`) and the *replaced* version of this same FOC (`domestic_bond_foc_F`) both convert via `(1+rb)·p/p(+1) − 1`. First-order effect: q_b_F mispriced by the expected ToT depreciation term whenever p moves.

### Proposed Fix
Use `rb_F_fg_next = (1 + rb_actual_F(+1))*p/p(+1) - 1` in `divert_bond_foc_F`, mirroring `domestic_bond_foc_F`.

### Status
Confirmed (inconsistency); Open (intent)

---

## A-1 — SS goods-market residuals only at solver tolerance, not identically zero

### Severity
Minor

### Confidence
High

### Location
SS solve (notebook cells 10–13); diagnostic cell 15

### Description
`goods_mkt_D = goods_mkt_F = −3.9e−7` at the final SS. Walras holds at SS analytically (K constant kills W-1; p constant kills W-2), so the residual is Broyden tolerance from the 3-unknown SS solve, propagated into anchors. Harmless for IRFs (linear deviations), but it biases the *level* of anchored constants (excess_return_ss etc.) at the 1e−7 scale.

### Status
Confirmed

---

## A-2 — Latent bookkeeping inconsistency: m vs Phi/T in SS vs dynamics

### Severity
Major (latent — inert in baseline because Phi = T = 0)

### Confidence
High

### Location
`equations_D.py`: `smart_steady_D` (m_D line ~103), `intermediation_P2_D`, `banker_div_res_D`, `banker_div_D`; F analogues

### Description
SS block: `m = n(1−(1−f)(1+rn)) + Phi + T`. Dynamic net-worth law P2: `n = (1−f)(1+rn)n(−1) + (1−f)cap_profit + m` ⇒ at SS requires `m = n(1−(1−f)(1+rn))` (no Phi+T). With Phi+T ≠ 0: (i) SS is not a rest point of P2 (residual = Phi+T); (ii) `banker_div_res` at SS gives `rn·n − 2(Phi+T)` vs `banker_div`'s `rn·n − (Phi+T)`. Currently chi1=T0=T1=0 so all coincide. Any experiment turning on intermediary adjustment costs or the macroprudential tax (e.g. TPI with T1>0) silently linearizes around an inconsistent point.

### Proposed Fix
Decide who pays Phi and T (continuing bank vs dividends). If dividends (current P2): remove `+ Phi_D + T_D` from `m_D` in `smart_steady_D` and remove the extra `− Phi − T` from one of the two dividend definitions so SS/dynamics agree. Verify P2 residual = 0 at SS with chi1, T0 > 0 as the acceptance test.

### Status
Confirmed (latent)

---

## S-1 — Baseline "default" never causes losses: writeoff_enabled = 0

### Severity
Major (interpretation / paper-claims risk, not a code bug)

### Confidence
High

### Location
Calibration (cell 2): `writeoff_enabled_D/F = 0`, `recovery_rate = 0`, `zeta_writeoff = 0`; `bond_return_*`, `budget_residual_*`

### Description
With the haircut multiplier 0, `rb_actual = [δ_b + (1−δ_b)q_t]/q_{t−1} − 1`: coupons and face value are **never** haircut, regardless of def_rate. Default risk transmits *only* through (i) the FOC spread loading `psi_spread·def_rate(+1)` and (ii) divertability `Δ + ψ_λB·def_rate(+1)` in the IC. The sovereign can default with probability 8%/quarter and bondholders lose nothing ex post; the government's budget also never receives debt relief. Any paper text describing realized defaults, haircuts, or recovery is unsupported by the baseline runs. The "doom loop" is a pure risk-premium loop.

### Status
Confirmed

---

## C-1 — Back-solved divertability Δ_cross ≈ 1.455 > 1; multi-asset IC degenerate at SS

### Severity
Major

### Confidence
High (numerically confirmed: Delta_bD_D = Delta_bF_F = 0.7273, Delta_bF_D = Delta_bD_F = 1.4545)

### Location
Notebook cell 12 (`_ic_delta`), `steady_auxilliary_*` (lambda_gk formula), `intermediation_IC_*`

### Description
`steady_auxilliary` computes λ_gk from the **single-asset** GK condition, which forces franchise value/λ = θ identically. The multi-asset IC then requires Σ(1−Δᵢ)φᵢ = 0, so with Δ_own < 1 the cross-bond Δ must exceed 1 (with ratio 2: Δ_cross = 1.4545): cross-border bonds are "145% divertable" — a bank can abscond with more than the asset's worth. Economically meaningless; also inverts the calibration table's documented intent (comments say Δ_bD = 0.2, Δ_bF = 0.4). The Δ levels are inert at SS (they cancel by construction) and matter only through the ψ_λB perturbation — i.e. the entire Δ apparatus reduces to a re-parameterized spread elasticity.

### Proposed Fix
Either (a) derive λ_gk from the multi-asset IC jointly with chosen Δ ∈ (0,1] (one more SS unknown), or (b) drop the Δ machinery and put the default-sensitivity directly in the FOC spread (which is what effectively remains).

### Status
Confirmed

---

## W-4 — Fix verification: leaks collapse 3–6 orders of magnitude; spillover results first-order contaminated

### Severity
Critical (consequence quantification of W-1/W-2)

### Confidence
High (numerical, `audit_artifacts/fix_test.py`)

### Location
`audit_artifacts/fix_results.json`, summarized in `walras_forensics.md` §7

### Description
With W-1+W-2 patched (in block copies; model files untouched): max |goods_mkt_F| falls 2.0e−2 → 4.2e−8, max |ca_res_D| 1.5e−4 → 4.1e−7 (SS-tolerance floor) on the TFP shock. Pre-fix cross-country IRFs are wrong by 39–124% (C_F, U_F, Y_F, n_inter_F, p); D-side variables move ≤ 7%. The perverse positive Y_D response to a default shock survives the fix (see T-2).

### Status
Confirmed

---

## T-2 — Deposit rate is contemporaneous on predetermined deposits → state-contingent deposits → perverse doom loop

### Severity
Critical (mechanism validity — explains the "perverse GDP increase under default shock" the git history shows was patched and reverted)

### Confidence
High (exact decomposition)

### Location
`deposit_return_D/F` (`Rgross_t = (1+rdep_t)P(−1)/P`), `bank_return_D/F` (`rdep_t` on t−1 balance sheet); `rdep_D/F` are period-t unknowns

### Description
Deposits placed at t−1 pay `rdep_t`, an unknown determined *after* the period-t shock: deposits are implicitly state-contingent. On a 1pp default shock, equilibrium rdep_D(0) = −1.49pp (flight to safety). Paid on the predetermined deposit stock (9× quarterly GDP), this transfers ≈ 0.134 from households to banks at impact — within rounding, the *entire* perverse n_inter_D(0) = +0.139 rise. Decomposition of rn_D(0): bond losses −0.25%, funding-cost windfall (κ leverage on −Δrdep) +5.4%, capital +1.4%. With rdep frozen at SS, rn_D(0) falls from +5.2% to +0.76%. The same mechanism makes n_inter_D(0) = −20.6% on a +1% TFP shock (rdep_D(0) = +9.4pp/q) — implausible volatility.

### Proposed Fix
Re-date: households earn `(1+rdep_D(-1))·P(−1)/P` and banks pay `rdep_D(-1)` on `Dep_{t−1}`; `rdep_t` (chosen at t) then applies to deposits formed at t. One-line change in `deposit_return_*` plus `rdep → rdep(-1)` in `bank_return_*`, and `rdep(+1) → rdep` in the FOC comparisons if the ex-ante rate is meant. This preserves deposit-market clearing (unlike the reverted "fix rdep exogenously" commit, which broke the unknown/target count).

### Follow-on discovery (2026-06-11, after applying the fix)
With T-2 applied, impact responses get the correct sign (Y_D(0) < 0 on a default shock) and Walras still holds at impact, **but the system becomes explosive** (b_gov_D −3.7e4 by t=499; H_U condition number 8e5 → 5.7e7). Interpretation: the state-contingent deposit rate was an *accidental automatic stabilizer* — the bank windfall after adverse shocks recapitalized banks and raised taxable dividends, choking off the debt–spread spiral. Removing it exposes the genuine doom loop, which is locally explosive under `phi_lamb = 0.02`. The notebook's own stability comment (cell 22) prescribes the remedy: raise `phi_lamb`.

**phi_lamb sweep (audit_artifacts/philamb_results.json):** 0.05/0.07/0.10 explosive or near-unit-root; 0.12 borderline (b_gov t499 = 7e−4); **0.15 minimal fully stable** (t499 ≈ 2e−5); impact economics insensitive across 0.12–0.30 (n_inter_D(0) ≈ −1.1%, spread(0) ≈ +30–40bp, Y_D(0) < 0). Applied phi_lamb_D/F = 0.15 in notebook calibration (cell 2, with comment). ⚠ This is a *calibration decision the authors should own*: alternatives are weakening def_scale or psi_lambda_B (weaker doom loop) instead of stronger fiscal feedback. 0.15/quarter is aggressive relative to Bohn-rule estimates — i.e. the post-fix doom loop is strong enough that only aggressive consolidation stabilizes debt. That tension is a *result*, not a nuisance.

**This is a feature, not a bug, of the fix: the doom loop the paper wants now exists, and fiscal policy strength becomes a substantive stabilization question rather than a free parameter.**

### Status
**FIXED in code** (`deposit_return_*`, `bank_return_*`, `intermediation_P1_*`, `divert_bond_foc_*`, `divert_portfolio_adj` re-dated; phi_lamb = 0.15). Verified: leaks ≤ 7e−8, stationary, doom-loop signs correct.

---

## TPI-1 — TPI experiment has no central-bank budget constraint: welfare gains partly bought with unbacked resources

### Severity
Critical (headline policy result)

### Confidence
High (analytical); numerical sizing in progress (tpi_test.py)

### Location
Notebook cells 39–47 (`domestic_bond_clearing_tpi`, `compute_tpi_irfs`, welfare figures)

### Description
`b_D_D = b_gov_D − b_D_F − cb_buy_D` removes bonds from private balance sheets, but no equation records the CB paying `q_D·cb_buy` for them, funding itself, or remitting coupon income. Extending the Walras identity: `goods_mkt_D + ca_res_D = factor_gap + q_D·cb_t − (1+rb_D,t)q_D(−1)cb_{t−1}`. The closed-loop rule cb = γ·spread therefore creates unbacked resource flows during stress, and the discounted-welfare bar chart (TPI-9) is computed on top of them.

### Evidence (tpi_test.py, closed-loop replication of cells 41–43)
| γ | max\|ca_res_D\| | corr(ca_res, predicted CB hole) | ΔW_D vs γ=0 | discounted unbacked flow (same units) |
|---|---|---|---|---|
| 0 | 6.3e−5 | — | — | 0 |
| 2 | 5.0e−3 | −0.9999 | +0.0236 | −0.0092 |
| 5 | 1.26e−2 | −0.99998 | +0.0561 | −0.0222 |
| 10 | **2.61e−2** | −0.99999 | +0.1027 | −0.0420 |

The Walras violation grows **417×** over the no-policy baseline at γ=10 and reaches **2.6% of quarterly GDP per period**; its path matches the analytical CB-hole formula with correlation |1.0000|. The discounted unbacked flow is ~40% of the measured welfare gain (opposite sign). The reported TPI welfare effects (ΔW_D ≈ +0.10, ΔW_F ≈ −0.08 % of SS consumption at γ=10) are therefore not interpretable: the experiment violates the aggregate resource constraint by amounts of the same order as the effects being measured.

### Proposed Fix
Add a CB budget: purchases financed by lump-sum levy on D households (or remittances netted into `budget_residual_D`): subtract `q_b_D·(cb_t − surv(1−δ_b)cb_{t−1}) − δ_b(1−def·h·m)cb_{t−1}` from the government's net issuance need (CB profits remitted), which both closes the accounting and models OMT-style sterilized purchases. Re-run welfare comparison after.

### Status
Confirmed (analytical) / numbers pending

---

## T-MAP — Timing audit map (Phase 4)

Convention: SSJ end-of-period stocks. Verdict per lag/lead operator group:

| Object | Implementation | Verdict |
|---|---|---|
| Bond return | `rb_actual_t = [δ_b(1−def_t·h·m) + (1−δ_b)q_t(1−ζ·def_t·h·m)]/q_{t−1} − 1` | ✓ buy at t−1, default/coupon realized at t |
| Default rate | `def_t = f(b_gov_{t−1}) + shock_t`; priced via `def(+1)` in FOCs/IC | ✓ debt chosen at t determines t+1 risk; consistent everywhere |
| Govt budget | coupon on `b_gov(−1)`, issuance at `q_t·(b_t − surv(1−δ_b)b_{t−1})` | ✓ matches bond cash-flow identity exactly |
| Fiscal rule | `T_ls,t = φ(b_gov_{t−1} − b_ss)` | ✓ predetermined |
| Deposit return | household `Rgross_t = (1+rdep_t)P(−1)/P`; bank cost `rdep_t` on `Dep_{t−1}` | ✓ both sides use same convention. ⚠ note: `rdep_t` (not `rdep_{t−1}`) pays on t−1 deposits ⇒ deposit rate is *not predetermined*; FOCs use `rdep(+1)` consistently, OK |
| Bank portfolio return | t−1 shares (`q(−1)b(−1)/n(−1)`, `θ(−1)`) × t returns | ✓ structure; ✗ F-bank return *units* (W-2) |
| Bellman P1 | `ν_t = SDF_t·Ω_{t+1}·(r_{t+1} − rdep_{t+1})`, `Ω_{t+1} = f + (1−f)λθ_{t+1}` | ✓ GK recursion |
| Net worth P2 | `n_t = (1−f)[(1+rn_t)n_{t−1} + cap_profit_t] + m` | ✓ (m constant: simplification, not error) |
| **Capital/production** | `Y_t = F(K_t, N_t)` but `rk_t` pays `mpk_t` on `K_{t−1}`; `K_t = (1−δ)K_{t−1} + Φ(I_t/K_{t−1})K_{t−1}` | **✗ OFF-BY-ONE (W-1).** Production must use `K(−1)` for income = output |
| Capital producer | `cap_profit_t = Q_t(K_t − (1−δ)K_{t−1}) − I_t` | ✓ given the fix to production timing |
| ToT conversions in FOCs | `divert_portfolio_adj`: `p(+1)/p` and `p/p(+1)` ✓; `divert_bond_foc_F`: **missing** (W-3); `bank_return_F`: **missing `p(−1)/p`** (W-2) | mixed |
| NFA | end-period market value `q_t·b_t`; receipts at t use `q_{t−1}b_{t−1}` | ✓ |

Phase-4 verdict: exactly one off-by-one (capital), two unit errors (F-bank); all sovereign/deposit timing clean.

---

## S-2 — Sovereign-risk accounting under writeoff_enabled = 1 (counterfactual audit)

### Severity
Major (latent — relevant the moment haircuts are switched on)

### Confidence
Medium-High (analytical only)

### Location
`bond_return_*`, `budget_residual_*`, `government_ss_*`, `intermediation_P2_F` comment

### Description
If `writeoff_enabled = 1`: bondholder loses `δ_b·def·h` of coupon and `ζ·def·h` of continuation face; government's budget saves exactly the same coupon and reduces survivors by the same `surv` factor. Bond-return and budget blocks are **mutually consistent** (no double-count: the haircut appears once in each side's cash flow, and the P2_F comment confirms an earlier double-count via separate writedown terms was already removed). Two remaining problems: (i) A-2 (m includes +Phi+T in SS only) is independent; (ii) with ζ = 0 but writeoff = 1, only coupons are haircut while face survives fully — a "default" that costs holders ~δ_b·def·h ≈ 0.1% per 1pp default — likely far weaker than intended; ζ is a hidden leverage on loss-given-default that is currently 0. Recovery_rate = 0 means haircut h = 1 (total loss on the affected coupon flow) — fine.

### Status
Confirmed (consistency); Open (whether ζ=0 is intended)

---

## X-0 — Notebook code audit (Phase 11)

### Severity
Minor–Major (reproducibility)

### Confidence
High

### Findings

1. **Execution-order dependence (Major):** cell 10 contains an *inline copy* of the anchor logic that **omits `psi_spread_D/F`**; these are only injected by `_apply_ss_anchors` (cells 11–13). Running cells 10 → 20 directly (skipping 12–13) fails or silently mismatches; `ss_final` is defined in cell 10 *and* redefined in cell 13 — partial re-runs linearize around stale calibrations. Fix: make cell 10 call `_apply_ss_anchors` (define it before cell 10) and create `ss_final` only once, in cell 13.
2. **Aliasing:** `cali_D = cali_F = ss` (cells 10, 13) — three names for one mutable object. Currently only floats are read; fragile.
3. **Stale diagnostics:** cell 25's annotations ("theta_D RISES = perverse IC") no longer match output (θ falls; n rises instead); cell 22's ρ_b formula hardcodes 0.953/0.95/0.05 unrelated to current calibration; cell-2 comments stale (yield 1.7%, debt 60%, Δ interpretations — see calibration_review.md).
4. **Dead value seeding:** `ss.toplevel['value_D/F']` set with a comment claiming P1 consumes it; no block takes `value_*` as input.
5. **Duplicated logic:** `cap_profit` computed in both `smart_steady_*` and `capital_producer_profit_*`; `banker_div_*` (SS) vs `banker_div_res_*` (dynamic) — consistent today only because Phi=T=0 (A-2).
6. **Repo hygiene:** `model_v11.ipynb` and `OLD models ` (trailing space in dirname) at top level; stray `code/` dir holding only `__pycache__`.
7. `model_notes` contains unrun scratch code with undefined references (`unknowns_tp` loop) — fine as notes, but half of it (ECB backstop, Cochrane perpetuity) describes features *not* in the model; do not cite as documentation.

### Status
Confirmed

---

## X-1 — Dead parameters and dead code

### Severity
Minor

### Confidence
High

### Location
Calibration cell 2; `equations_global.py`; `equations_D/F.py`

### Description
Never referenced by any equation: `lambda_BD_D/F`, `lambda_BF_D/F`, `psi_nu_bD_*`, `psi_nu_bF_*`, `mc_D/F`, `lamb_ss_D/F`; `calibration_hh_D/F` dicts built and never used (incl. hardcoded `beta_D = 0.9920094934` that is *not* the solved beta). Dead blocks: `domestic_bond_foc_D/F`, `portfolio_adj_cost` (imported in cell 7 but not in `ha_full`). Dead branch: `nZ == 19` loading `Discretisation/Outputs/Px_GMAR.txt` (hardcoded path). Household output `uce_*` aggregated but unused.

### Proposed Fix
Delete or mark explicitly. The unused 19-state branch with a file dependency is a reproducibility hazard.

### Status
Confirmed

---
