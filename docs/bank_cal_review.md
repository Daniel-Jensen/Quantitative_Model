# HHBANK v12 — bank-cal Branch Review

**Date:** 2026-06-11
**Question:** what is in `bank-cal`, what was lost, and should it be merged after the structural fixes?

## Branch facts
- `bank-cal`: **13 ahead / 19 behind `main`**. Branched before the "Sunday fixes" PR #25 and before the audit. Model code lives under `code/equations_*.py` (it predates the top-level layout). Adds `refs/` (literature PDFs), `requirements.txt`, calibration commits, conference slide figures.
- It therefore **does not contain** any of W-1, W-2, W-3, T-2, A-2, TPI-1. Merging the branch would revert all of them.

## What bank-cal changed (calibration commits)
| Commit | Content |
|--------|---------|
| 61cf2e8 | Empirical def_scale/phi_lamb targets (spread-debt slope ≈ def_scale×430 bps/pp; Bohn mapping) |
| 7529353 | def_scale bifurcation: three regimes; bifurcation at def_scale≈0.305 (phi_lamb=0.40), ≈0.13 (phi_lamb=0.03) **on the old model** |
| ff3afe1 | Bond duration → 2011 maturities: δ_b_D=0.036 (GR 7yr), δ_b_F=0.038 (DE 6.5yr); documents impact magnitude duration-invariant, persistence not |
| 26baed2 | phi_lamb 0.40→0.03 (annual Bohn 1.6→0.12, Staehr 2008 EA periphery 0.10–0.15) |
| 17b5ff3 | SS convergence for long-duration bonds; safe def_scale scenarios |

## Which calibration currently exists where
- **main / audit:** phi_lamb=0.15, def_scale=0.25, δ_b=0.10, f=0.12, psi_lambda_B=3.0, recovery=0, zeta=0, writeoff=0, Δ back-solved (Δ_cross=1.45), chi1=0, frisch=0.5, EBA via cell-12 share targets (b_F_D/asset≈0.15).
- **bank-cal:** phi_lamb=0.03, def_scale=0.05, δ_b=0.036/0.038, f=0.03, psi_lambda_B=0.0, recovery=0.40, zeta=1.0, writeoff=0, **Δ hardcoded 0.2/0.4** (no back-solve → no Δ>1), chi1=0.5, frisch=1.0, EBA bilateral (b_F_D/asset=0.18%, b_D_F/asset=0.65%), income rho_z_F=0.97/sigma=0.8/nZ=7.

## What bank-cal does better (and was lost from main)
1. **Bond duration (δ_b≈0.037 vs 0.10).** Directly addresses the audit's calibration_review concern that δ_b=0.10 (2.5yr) understated the MTM doom-loop channel. bank-cal ties it to GR/DE 2011 maturities. **Port this.**
2. **GK exit f=0.03 vs 0.12.** More standard (GK 0.03–0.10); longer franchise horizon, stronger net-worth persistence. **Port.**
3. **Δ hardcoded 0.2/0.4 — avoids C-1 entirely.** No `_ic_delta` back-solve, so no Δ_cross=1.45. Trade-off: IC may not bind exactly at SS (λ_gk absorbs the slack). Still preferable to a >1 divertability. **Strongly consider porting** as the C-1 resolution.
4. **EBA bilateral exposures** (cross-border 0.18%/0.65% of assets) — empirically grounded, supersedes main's round-number share targets. **Port.**
5. **recovery=0.40, zeta=1.0** set to sensible Greek values — but inert while writeoff_enabled=0 (S-1). Porting these *plus* flipping writeoff_enabled=1 would finally give default realized losses.
6. **Empirical phi_lamb / def_scale mapping** (Bohn coefficient, spread-debt slope) — valuable literature anchoring.

## What bank-cal does NOT fix
- **S-1 unchanged:** writeoff_enabled=0 even in bank-cal. Default still produces zero realized losses despite recovery/zeta being set. Same gap as main.
- All structural accounting/timing bugs (predates them).

## Stability of the fixed model at bank-cal phi_lamb (the crux)

The empirical phi_lamb≈0.03 was tuned **on the buggy model**, where the T-2 deposit windfall accidentally stabilized debt. The fixed model removes that stabilizer.

Probe (`audit_artifacts/bankcal_stability_test.py`, fixed model, psi_lambda_B=0 / def_scale=0.05 / psi_spread=0):
| phi_lamb | stationary? | leaks |
|----------|-------------|-------|
| 0.03 | **yes** | ca_res ≤1.1e−7 |
| 0.05–0.12 | yes | ≤1e−7 |

⚠ **Caveat (important):** zeroing `psi_spread` and `psi_lambda_B` together with writeoff=0 removed *every* default-transmission channel — the default-shock IRFs came back identically 0.0 (default became inert). So this probe proves only that the **fiscal block alone** is stable at phi_lamb=0.03; it does **not** prove the doom loop is stable there, because the probe has no doom loop. The TFP shock still transmits (n_inter_D(0)=+0.35 on 1% Z), confirming the rest of the model is live.

**Honest conclusion:** phi_lamb and amplification strength (def_scale, psi_spread/psi_lambda_B, δ_b) **trade off**. Contrast:
- main amplification (psi_lambda_B=3, def_scale=0.25): needs **phi_lamb=0.15** (Bohn 0.6, too aggressive).
- zero amplification: stable at **phi_lamb=0.03** but default does nothing.
The empirically-defensible regime is in between and **must be re-mapped on the fixed model**: redraw the (phi_lamb, def_scale) bifurcation that bank-cal drew on the old model. bank-cal's "bifurcation at def_scale≈0.13 (phi_lamb=0.03)" is **not valid post-fix** — the T-2 stabilizer is gone, so the stable def_scale ceiling at phi_lamb=0.03 will be lower.

## Recommendation

**Do not merge `bank-cal`.** Instead, on the `audit` branch:
1. Port calibration *values*: δ_b=0.036/0.038, f=0.03, EBA bilateral exposures, Δ hardcoded 0.2/0.4 (resolves C-1), recovery=0.40, frisch (decide 0.5 vs 1.0), income process.
2. Decide S-1: set writeoff_enabled=1 (with recovery=0.40, zeta=1.0) to give default real losses — or keep the pure-risk-premium framing explicitly.
3. **Re-map (phi_lamb, def_scale) stability on the fixed model** with the ported duration/amplification. Find the empirically-plausible corner that is both stationary and has a non-trivial doom loop. This is the gating task.
4. Re-run all IRF/TPI/welfare on the merged calibration.

Porting values (not code) keeps every structural fix and pulls in the empirical work — the best of both branches.

## Phase 4 — Calibration parameter assessment

| Parameter | Main / audit | bank-cal | Literature | Recommendation |
|-----------|-------------|----------|------------|----------------|
| `phi_lamb` (Bohn) | 0.15 (Bohn 0.6) | 0.03 (Bohn 0.12) | 0.025–0.038 (Staehr 2008, EA periphery 0.10–0.15 annual) | **Re-map on fixed model.** bank-cal's value is literature-correct but pre-fix; 0.15 is the only stationary value at main's strong amplification. Target: lowest amplification giving a non-trivial doom loop, then minimal stationary phi_lamb (likely 0.05–0.10). |
| `def_scale` | 0.25 | 0.05 | 0.007–0.047 (normal); 0.12–0.23 (GR 2011 peak) | **0.05–0.10.** main's 0.25 exceeds even crisis-peak; bank-cal's 0.05 ≈ moderate. Re-check bifurcation post-fix (old ceiling 0.13 at phi_lamb=0.03 no longer valid). |
| `delta_b` (duration) | 0.10 (2.5yr) | 0.036/0.038 (7/6.5yr) | GR 2011 ~7yr, DE ~6.5yr | **Port bank-cal.** Realistic duration; lengthens MTM doom-loop persistence (audit flagged 0.10 as too short). |
| `f` (bank exit) | 0.12 | 0.03 | GK 0.03–0.10 | **Port 0.03.** Standard GK; main's 0.12 shortens franchise horizon. |
| `psi_lambda_B` (state-dep divertability) | 3.0 | 0.0 | no empirical counterpart | **Reduce toward bank-cal.** This is a free amplification dial; 3.0 drives most of main's spread response and forces phi_lamb=0.15. Set low, source amplification from def_scale + duration instead. |
| `Delta_cross` | 1.4545 (back-solved) | 0.4 (hardcoded) | must be ≤1 (divertable fraction) | **Port bank-cal's hardcoded 0.2/0.4.** Resolves C-1. Accept small SS IC slack (λ_gk absorbs). |
| `recovery_rate` | 0.00 | 0.40 | Greece 2012 ~0.5 haircut → recovery ~0.4–0.5 | **Port 0.40** — but inert unless writeoff_enabled=1. |
| `zeta_writeoff` | 0.0 | 1.0 | — | Port 1.0 (full write-off) with writeoff_enabled=1. |
| `writeoff_enabled` | 0.0 | 0.0 | — | **Set 1.0** to resolve S-1 (give default realized losses), or keep 0 and reframe paper as pure risk-premium. |
| `theta` (leverage) | 4 | 4 | GK 4–6; 2011 historical 10–25 | Keep 4 as conservative baseline; note historical is higher. |
| `frisch` | 0.5 | 1.0 | 0.5–1.0 | Either defensible; bank-cal=1.0 is common in NK. Author choice. |
| `chi1` (interm. adj cost) | 0.0 | 0.5 | — | If porting chi1=0.5, **A-2 fix is required** (present on audit branch); would have been silently inconsistent on bank-cal's own code. |
| income `rho_z_F/sigma/nZ` | 0.90/0.3/15 | 0.97/0.8/7 | — | bank-cal's higher persistence/dispersion → fatter wealth distribution; fewer grid points (7) is coarse. Keep main's nZ=15 with bank-cal's rho/sigma if porting. |

**Cross-finding:** bank-cal's chi1=0.5 (adjustment costs on) would have triggered the audit's A-2 latent bug on bank-cal's own (pre-fix) code — the SS `m` would not be a rest point of P2. The audit branch's A-2 fix makes porting chi1=0.5 safe. This is an example of a structural fix being a prerequisite for a bank-cal calibration choice.
