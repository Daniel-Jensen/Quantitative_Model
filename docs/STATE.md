# Project State

**Branch:** `audit` | **Date:** 2026-06-11 | **Status:** post-forensic-audit baseline

## Current status

Six structural/accounting bugs were found and fixed in the 2026-06-11 forensic audit (W-1, W-2, W-3, T-2, A-2, TPI-1). See `docs/audit.md` for the full ranked finding list and `docs/verification_report.md` for verified fix status. All six are applied on branch `audit`; PR #26 is open for co-author review.

Core equations: `code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py`. Active notebook: `code/model_v12.ipynb`. TPI output figures: `plots/`.

`main` is deliberately left at the pre-fix state to preserve a clean PR diff. Do not use `main` for new work until PR #26 is merged.

## What is complete (post-audit)

- Household deposit choice and GHH preferences for D (Greece) and F (Germany) households.
- Bank steady-state and intermediation blocks: capital, bond returns, fees, GK Bellman (P1) and IC constraint (P3/lambda_gk).
- Production, capital adjustment, and capital producer profit — W-1 fixed: production uses `Y=F(K_t)`, capital producer receives `mpk·(K−K(−1))` so all capital income is allocated.
- Deposit return predetermined correctly: `Rgross = (1+rdep(−1))·P(−1)/P` — T-2 fix. Funding legs in `bank_return_*` and FOCs in `intermediation_P1_*`/`divert_*` use ex-ante deposit rate.
- F-bank bond returns converted to F-goods via `p(−1)/p` in `bank_return_F` — W-2 fix (the dominant leak; was causing ~2% of F GDP goods_mkt_F residual on a 1% TFP shock).
- Cross-border bond FOC in F-bank uses `p/p(+1)` for expected return — W-3 fix (optimality condition; does not affect Walras but required for internally consistent portfolio choice).
- Smart steady-state blocks: `m = n·(1−(1−f)·(1+rn))` without spurious `+Phi+T` — A-2 fix (required for any `chi1≠0` calibration, e.g. bank-cal's chi1=0.5).
- Global goods market, external account, bond clearing, and portfolio adjustment cost blocks.
- Domestic and foreign bond pricing, yields, spreads, and Hatchondo-Martinez geometric-decay perpetuity default mechanics.
- TPI extension (cells TPI-1/TPI-2 in notebook): CB budget closed via `budget_residual_D_tpi` with `rem_cb_D` remittance — TPI-1 fix. Before the fix, unbacked CB flows inflated welfare gains by ~40% at γ=10.
- EBA bilateral sovereign exposures in calibration cell: b_D_D/asset=24.47%, b_F_F/asset=25.79%, b_F_D/asset=0.18%, b_D_F/asset=0.65%.

## Walras accounting (post-fix, verified)

| Residual | 1% TFP-D shock | 1pp default-D shock |
|----------|----------------|---------------------|
| goods_mkt_D (targeted) | ≤1e−16 | ≤1e−16 |
| goods_mkt_F (untargeted) | ≤8e−10 | ≤1e−9 |
| ca_res_D = CA−ΔNFA (untargeted) | ≤5.8e−8 | ≤3.5e−8 |
| deposit_mkt_D/F | ≤4e−15 | ≤4e−15 |

Pre-fix peaks for reference: goods_mkt_F 2.0e−2 (~2% of F GDP); ca_res_D 1.5e−4. All cross-country spillover and welfare results from the pre-fix model are first-order invalid and must be regenerated from `audit` branch.

## IRF summary (post-fix, audit branch, phi_lamb=0.15)

**1pp default shock to D (ρ=0.8):**
- `n_inter_D[0] = −3.5%` (falls), `Y_D[0] = −2.5e−4` (falls) — both signs correct post-T-2-fix; were positive/perverse pre-fix.
- `n_inter_F[0] ≈ −0.33%` — contagion small, sign correct.
- Spread widens on impact; doom loop is live with correct sign.

**TPI (γ=10, post-fix):**
- ΔW_D = +1.88% SS consumption equivalent; ΔW_F = −1.90%. TPI is approximately a zero-sum burden transfer from D to F; spread is not compressed (rises slightly with γ because default is debt-driven). All pre-fix TPI welfare figures in `plots/` are stale until notebook is re-run from `audit` branch.

## Calibration summary (current, audit branch)

| Parameter | Value | Source / note |
|-----------|-------|---------------|
| `phi_lamb_D/F` | 0.15 | Bohn=0.60/yr; min stable at current amplification. Literature: 0.10–0.15/yr (Staehr 2008 EA periphery). Tension: bank-cal's 0.03 was tuned on pre-fix model. Re-map needed (see §Next priorities). |
| `def_scale_D` | 0.25 | Strong amplification. Exceeds GR 2011 crisis peak (0.12–0.23 from spread-debt slope calibration). |
| `delta_b_D/F` | 0.10 | 2.5yr avg maturity. Empirically too short; bank-cal has 0.036/0.038 matching GR/DE 2011 ~7yr/6.5yr. |
| `theta_D/F` | 4.0 | GK leverage; conservative vs 2011 historical 10–25×. |
| `psi_lambda_B_D/F` | 3.0 | State-dependent divertability; primary amplification dial; no direct empirical counterpart. |
| `f_D/F` | 0.12 | Bank exit rate; bank-cal has 0.03 (standard GK range). |
| `Delta_cross` | 1.4545 | Back-solved (`_ic_delta`, ratio=2.0); degenerate >1. See C-1. |
| `recovery_rate_D` | 0.00 | No realized losses; inert while writeoff_enabled=0. |
| `writeoff_enabled_D/F` | 0.0 | Default produces no balance-sheet losses. See S-1. |
| `chi1_D/F` | 0.0 | Intermediation adjustment cost off. A-2 fix makes chi1≠0 safe. |
| `frisch` | 0.5 | Frisch elasticity. |
| nDep_D/F | 500/500 | Household deposit grid points. |
| income rho_z/sigma/nZ | 0.90/0.30/15 | D and F income process (Markov approximation). |

## Open issues

| ID | Description | Status |
|----|-------------|--------|
| C-1 | `Delta_cross=1.45 > 1`: divertable fraction exceeds 1, making the multi-asset IC constraint degenerate. `lambda_gk` absorbs the slack but the theoretical interpretation breaks down. | Author decision. Preferred resolution: hardcode `Delta=0.2/0.4` per bank-cal (avoids back-solve entirely). |
| S-1 | `writeoff_enabled=0`: default shock produces zero realized bank losses. `recovery_rate` and `zeta_writeoff` are set but inert. Model is currently a pure risk-premium loop, not a balance-sheet doom loop. | Author decision. Resolution: set `writeoff_enabled=1` with `recovery=0.40, zeta=1.0` (GR 2012 ~50% haircut → ~0.4–0.5 recovery). |
| X-1 | Dead-code imports in notebook cell 7: blocks no longer in the model remain in the import list. | Minor cleanup; no numerical effect. |

## Next priorities

1. **Port bank-cal calibration values** onto `audit` branch: `delta_b=0.036/0.038`, `f=0.03`, EBA bilateral exposures (verified targets from bank-cal cell `96c6bd50`), hardcode `Delta=0.2/0.4` (resolves C-1), `recovery=0.40`. See `docs/bank_cal_review.md` §Recommendation for the full porting list.
2. **Decide S-1**: set `writeoff_enabled=1` to give default realized losses, or keep pure risk-premium framing and state it explicitly in the paper.
3. **Re-map (phi_lamb, def_scale) stability on the fixed model** with ported duration and amplification. Bank-cal's bifurcation diagram (bifurcation at def_scale≈0.13 at phi_lamb=0.03) is invalid post-T-2-fix — the accidental deposit-windfall stabilizer is gone. Find the lowest empirically-plausible phi_lamb that gives a non-trivial, stable doom loop. This is the gating calibration result.
4. **Re-generate all figures** from `audit` branch. All figures in `plots/` and in the notebook were produced on a pre-fix or mid-fix model state.
