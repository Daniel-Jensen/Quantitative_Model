# Project State

## Current status

- The core two-country HH-bank model is implemented in `model_v12.ipynb`.
- The model is structured around domestic and foreign blocks in `equations_D.py`, `equations_F.py`, and `equations_global.py`.
- `model_v12` adds New Keynesian labour-side features and a portfolio adjustment cost on cross-border bond holdings.
- The code can solve for a steady state and is designed to generate impulse response functions.

## What is complete

- Household deposit choice and GHH preferences for D and F households.
- Bank steady-state and intermediation blocks, including capital, bond returns, and fees.
- Global trade and external account clearing.
- Domestic and foreign bond pricing, yields, and spreads.
- Portfolio adjustment cost blocks anchoring cross-border bond positions.
- **Corrected EBA bilateral sovereign exposure calibration** (bond-to-assets, system aggregates):
  - b_D_D/asset_D = 24.47%, b_F_F/asset_F = 25.79%, b_F_D/asset_D = 0.18%, b_D_F/asset_F = 0.65%
- **Portfolio frictions reduced** `psi_bF_D = psi_bD_F = 0.5` (from 1.5) to avoid over-penalising small cross-border positions.
- **Steady state validated** (2026-06-08): all IC, P1, and labour-market residuals < 1e-15; goods-market residuals < 1e-7.
- **IRF validated** (1pp default shock to D, ρ=0.8):
  - `n_inter_D[0] = −5.23%`, `n_inter_F[0] = −0.33%` — domestic net worth falls ~16× more than foreign ✓
  - `q_b_D[0] = −70.7 bps`, `q_b_F[0] = −5.7 bps` — contagion to F is small ✓
  - `spread_rb[0] = +167 bps` — spread widens sharply on impact ✓
  - `goods_mkt_D[0] ≈ 0` (targeted); `goods_mkt_F[0] = 7.2e-4` — Walras residual on untargeted F market, decays quickly, within linearisation tolerance ✓
  - `b_F_D ↑`, `b_D_F ↑` on impact: flight-to-quality (D banks) and opportunity carry (F banks) — economically correct, small absolute magnitudes ✓

## Calibration summary (SS)

| Parameter | Value | Note |
|-----------|-------|------|
| `n_inter_D = n_inter_F` | 3.0 | net worth |
| `theta_D = theta_F` | 4.0 | leverage (conservative vs 2011 historical 10–25×) |
| `asset_D = asset_F` | 12.0 | total bank assets |
| B_supply_D / Y_D (annualised) | 75.4% face, 46.9% market | vs GR 2011 ~150% historical |
| B_supply_F / Y_F (annualised) | 77.9% face, 47.3% market | vs DE 2011 ~80% ✓ |

## Known limitations

- `theta = 4` is 2.5–6× more conservative than 2011 historical leverage (10–25×). Baseline "stress" scenario; raise to 10–15 to match history if needed.
- D annualised face debt/GDP = 75% vs GR historical 150%: model output normalised to 1 so debt in face-value units is under-scaled. Model-design artifact.
- `goods_mkt_F` not in `targets_tp` (clears by Walras). Linearisation residual ~7e-4 at impact is acceptable but worth monitoring if shocks are large.
- ECB/ESM backstop and macroprudential extensions are conceptual only.

## def_scale calibration: empirical targets (2026-06-08)

### Analytical slopes

Two structural slopes tie `def_scale` and `phi_lamb` to the literature:

**1. Spread-debt slope** = d(ex-ante spread)/d(annual B/Y):
- Endogenous channel: `def_scale × f'(0) × 4 × d(spread)/d(def_rate)`
- `f'(0) = curv × offset^(curv−1) = 0.5 × 0.05^(−0.5) = 2.236`
- `d(spread)/d(def_rate) ≈ delta_b × (1−recovery) / q_b = 0.05×0.60/0.623 ≈ 0.048 per pp`
- **Net: ≈ def_scale × 430 bps per pp of annual B/Y**

| Calibration scenario | Lit target | Implied def_scale |
|---|---|---|
| Pre-crisis periphery | 3–5 bps/pp (Laubach 2009) | 0.007–0.012 |
| Crisis regime | 10–20 bps/pp (Haugh et al 2009) | 0.023–0.047 |
| GR 2011 peak | 50–100 bps/pp | 0.12–0.23 |

**2. Bohn coefficient** = d(PB/Y)/d(annual B/Y):
- `phi_lamb` acts on quarterly b_gov: `d(PB/Y)/d(ann B/Y) = phi_lamb × 4 = 0.4 × 4 = 1.6`
- Literature: 0.03–0.07 (Bohn 1998; Checherita-Westphal 2012), 0.10–0.15 (EA periphery, Staehr 2008)
- **Current phi_lamb = 0.4 implies Bohn = 1.6 — 10–50× above literature range**
- phi_lamb was calibrated for numerical stability, not fiscal reaction function
- To match EA periphery (0.10–0.15): phi_lamb_target = 0.10 / 4 = 0.025 to 0.15/4 = 0.038
- **Implication**: over-strong Bohn rule suppresses the fiscal-financial feedback; for papers that study debt-default spirals, phi_lamb should be reduced

### Paper scenarios (recommended)

| Scenario | def_scale | phi_lamb | Spread-debt | Bohn coef |
|---|---|---|---|---|
| Baseline (exogenous) | 0.0 | 0.40 | 0 endogenous | 1.6 (over-stabilised) |
| Moderate endogenous | 0.05 | 0.03 | ~22 bps/pp | 0.12 |
| Strong endogenous | 0.10 | 0.03 | ~43 bps/pp | 0.12 |
| Crisis | 0.20 | 0.03 | ~86 bps/pp | 0.12 |

Note: phi_lamb = 0.03 gives annual Bohn = 0.12 — top of EA periphery range. Reduce further to 0.015 for core-country baseline.

⚠️ Keep def_scale ≤ 0.20 — bifurcation at def_scale ≈ 0.305 (see below).

## def_scale probe: endogenous default risk (2026-06-08)

`def_scale_D` controls how strongly debt-to-output feeds back into default risk:
`def_rate_D = shock_def_D + def_scale_D * ((debt_gap + offset)^curv - offset^curv)`

At SS `debt_gap = 0` so SS is identical for all `def_scale` values. Only the Jacobian changes.

**Three regimes** (1pp shock_def_D, ρ=0.8, T=100):

| Regime | def_scale | spread[0] bps | n_inter_D[0] | lamb_D peak | b_gov peak | Notes |
|--------|-----------|--------------|--------------|-------------|------------|-------|
| Baseline | 0.0 | +167 | −5.2% | −0.13 pp | +0.10% | exogenous only |
| Stable amplification | 0.1 | +204 | −6.6% | −0.16 pp | +0.13% | +22% more spread |
| | 0.2 | +288 | −9.5% | −0.24 pp | +0.20% | +72% more spread |
| Near-critical | 0.3 | +2189 | −77% | −2.27 pp | +1.89% | explosive; austerity trap |
| **Bifurcation** | **≈0.305** | — | — | — | — | **Jacobian eigenvalue = 1; model indeterminate** |
| Sign-flip zone | 0.31–0.35 | sign reversed | sign reversed | — | — | linearisation unreliable |
| Re-stabilised | 0.5–2.0 | +102–162 | −2.6 to −4.8% | −0.05 to −0.09 pp | weak | Bohn rule dominates |

**Mechanism**: higher `def_scale` → default risk rises with debt → spreads blow up → bank net worth collapses → output falls → debt/GDP rises → more default risk (self-reinforcing). Bohn rule (φ_lamb = 0.4) stabilises via fiscal tightening but is overwhelmed near the critical point.

**Recommended calibration range**: `def_scale ∈ [0.0, 0.20]` — monotonic, well-behaved amplification. `def_scale = 0.1` gives +22% spread amplification and +26% bank loss amplification vs pure exogenous. `def_scale > 0.25` approaches instability and linearization breaks down.

## Next focus

- Decide whether to raise `theta` to match 2011 leverage.
- Consider adding `goods_mkt_F` as an explicit check (not target) in the IRF diagnostic cell.
- Calibrate `def_scale` for the paper: 0.0 (baseline), 0.1 (moderate endogenous), 0.2 (strong endogenous) as three scenarios.
- Policy extensions: ECB backstop, macroprudential bond tax.
