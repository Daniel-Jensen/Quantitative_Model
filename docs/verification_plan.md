# HHBANK v12 — Post-Audit Verification Plan

**Purpose:** verify the audit's recommended fixes are implemented correctly on disk, find regressions, assess calibration vs the `bank-cal` branch. Not a re-audit.

**State of repo (2026-06-11):** fixes live on branch `audit` (commits 396cbd9 + 4c810e1); `main` deliberately left pre-fix (7e94175). Working tree currently = `audit`. `bank-cal` = separate calibration branch, **13 ahead / 19 behind main**, code under `code/`, predates all structural fixes.

## Outstanding items extracted from audit docs

| ID | Recommended fix (audit) | Acceptance test |
|----|------------------------|-----------------|
| W-1 | Close capital-income gap. Author chose: keep `Y=F(K_t)`, add `mpk·(K−K(−1))` to `cap_profit_*` | max\|ca_res_D\| ≤ 1e−7 all shocks; CA=ΔNFA |
| W-2 | `bank_return_F`: convert bond returns to F-goods via `p(−1)/p` | max\|goods_mkt_F\| ≤ 1e−7 |
| W-3 | `divert_bond_foc_F`: expected return via `p/p(+1)`; sweep for other missing conversions | no contemporaneous D-return vs F-rate mismatches remain |
| T-2 | Deposit rate predetermined: liabilities pay `rdep(−1)`; FOCs use ex-ante `rdep` | n_inter_D & Y_D **fall** on default shock; stationary |
| A-2 | Remove `+Phi+T` from `m` in `smart_steady_*` | P2 residual=0 at SS with chi1/T0/T1≠0 |
| TPI-1 | Add CB budget/remittance to TPI experiment | ca_res_D γ-invariant; welfare interpretable |
| C-1 | Multi-asset-consistent λ_gk OR drop Δ apparatus | Δ_cross ≤ 1 |
| S-1 | Decide loss regime (writeoff=1, ζ>0, δ_b≈0.04) or reframe paper | default produces realized balance-sheet losses |

## Verification method
1. Static: grep each fix on disk (done — all structural fixes present).
2. Dynamic regression: `audit_artifacts/run_audit.py` reads live `equations_*.py` → `run4_kt_log.txt` is the regression test (leaks, signs, stationarity).
3. Calibration: diff `bank-cal` calibration vs current; test whether the *fixed* model is stable at bank-cal-style low-amplification params (`bankcal_stability_test.py`).

Deliverables: `verification_report.md`, `bank_cal_review.md`, executive summary (below in report).
