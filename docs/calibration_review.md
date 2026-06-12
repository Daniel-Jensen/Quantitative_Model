# HHBANK v12 — Calibration Review

**Status:** complete (focused per prioritisation: 5% budget). Values below are the *final* post-pipeline values (after notebook cells 12–13 overwrite the cell-2 entries), from the audited re-run (`audit_artifacts/run1_results.json`).

## Headline inconsistencies between stated calibration and effective calibration

| Parameter | Stated (cell 2) | Effective (final SS) | Issue |
|---|---|---|---|
| `lambda_gk` | 0.2 | **0.3459** | Overwritten by `steady_auxilliary` single-asset formula. Comment value never used in dynamics. |
| `Delta_bD_D`/`Delta_bF_F` | 0.2 ("preferred collateral") | **0.7273** | Back-solved in cell 12. |
| `Delta_bF_D`/`Delta_bD_F` | 0.4 | **1.4545** | **> 1: cross-border bonds "145% divertable" — economically meaningless.** Forced by single-asset λ_gk (value/λ ≡ θ ⇒ Σ(1−Δᵢ)φᵢ = 0). See audit.md C-1. |
| `delta` (depreciation) | 0.025 | 0.02241 | Recalibrated to hit rk = 1%/q. Fine, but cell-2 comment stale. |
| Debt/GDP | "≈60% of annual GDP" | **30.7%** (b_gov = 1.230 face, q = 0.9757, Y_ann = 4) | Cell-12 portfolio-share targeting silently halves the debt stock. |
| Bond yield | "β_inter → 1.7%/yr" | **rb ≈ 0.249%/q ≈ 1.0%/yr** | Comment stale. Spread over deposits (rdep = 0) = 1.0%/yr. |
| `beta_D` (hh) | guess 0.985 / dead `calibration_hh` value 0.9920 | **0.99946** | Near-unity; forced by DEP/Y = 9 (225% of annual GDP in deposits) at rdep = 0. |

## Parameter assessment table

| Param | Value | Interpretation | Literature range | Assessment |
|---|---|---|---|---|
| eis 0.5, frisch 0.5 | — | standard | 0.5–1 / 0.5–1 | OK |
| α 0.35, K/4Y = 2.7 | — | standard | ✓ | OK |
| δ_b 0.10 | avg debt maturity 1/δ_b = 10q = 2.5y | EA sovereigns ~6–7y (δ_b ≈ 0.04) | **Short.** Halves the duration channel: capital loss on q_b per unit spread is ~duration-proportional, so doom-loop amplification via bond MTM losses is mechanically understated by factor ~2–3. Works *against* the paper's mechanism (conservative), but inconsistent with "long-term bonds" framing. |
| f 0.12 | banker exit 12%/q | GK: 1−σ ≈ 0.03–0.10/q | Slightly high; shortens franchise horizon, weakens net-worth channel persistence. |
| θ 4 | bank leverage | GK 4–6, EA banks ~10+ on sovereign books | Defensible for "broad" intermediaries. |
| φ_bD+φ_bF = 0.40 | bonds = 40% of bank NW = 10% of bank assets | EA periphery banks: sovereign ≈ 5–10% of assets | OK in assets terms. |
| `def_scale` 0.25, curv 0.5, offset 0.05 | d(def)/d(debt ratio) = 0.25·0.5·(1.28)^−0.5 ≈ **0.11/q** | — | A 10% rise in debt (Δratio 0.123) ⇒ +1.36pp *quarterly* default prob ⇒ +0.778·1.36 ≈ +1.06pp/q ≈ +4.2pp annual spread. **Strong**, not suppressive. But see S-1: zero cash-flow consequence. |
| `psi_spread` (derived) 0.778 | spread per unit E[def] | risk-neutral pricing with full loss ⇒ ≈ 1·haircut | <1 with haircut=1: mildly *under*-prices risk. Internally derived from λψ_λB/βΩ — at least consistent. |
| `psi_lambda_B` 3.0 | Δ_eff sensitivity | no direct empirical counterpart | Free amplification dial; the *only* live role of the entire Δ apparatus. |
| `psi_bF_D` = `psi_bD_F` = 0.5 | level penalty on cross-border face positions | no counterpart; chosen by half-life grid loop (model_notes) | **Suppressive by design**: pins b_F_D, b_D_F near SS ⇒ portfolio-rebalancing channel is heavily damped (IRF: b_F_D moves ~0.1% of its SS on a 1pp default shock). The "contagion via portfolio reallocation" mechanism is therefore quantitatively shut by calibration. |
| `phi_lamb` 0.02 | Bohn coefficient | Bohn estimates ~0.02–0.1 | OK; debt half-life ~35q. |
| `writeoff_enabled` 0 | no default losses ever | — | **The single most consequential calibration switch** (S-1): kills bondholder losses, recovery, and debt relief. Doom loop = pure risk-premium loop. |
| ω 0.85, ε_trade 1.5 | home bias, Armington | 0.7–0.9 / 1–2 | OK |
| G/Y 8.2%, TAX/Y 8.5% | residual fiscal scale | data ~20–45% | Small government: understates the fiscal-capacity stakes of default risk. |
| rdep = 0 (SS) | deposit rate fixed at SS | — | Hardwires a 1%/yr SS spread; spread is a *choice*, not an outcome of the friction (excess-return anchors absorb whatever SS implies). |

## Mechanisms suppressed or hardwired by calibration (H7 verdict)

1. **Suppressed:** portfolio reallocation (ψ = 0.5 level penalties chosen for IRF half-life cosmetics, not data).
2. **Suppressed:** MTM bond-loss amplification (δ_b = 0.1 short duration; and writeoff = 0 removes face losses entirely).
3. **Hardwired:** SS spread (β_inter chosen + anchors), so IC friction does not *determine* the spread level, only its dynamics.
4. **Not suppressed:** endogenous default sensitivity (def_scale strong) and divertability spread loading (ψ_λB = 3).

H7 status: **PARTIALLY CONFIRMED** — see audit.md.
