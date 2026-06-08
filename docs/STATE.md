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

## Next focus

- Decide whether to raise `theta` to match 2011 leverage.
- Consider adding `goods_mkt_F` as an explicit check (not target) in the IRF diagnostic cell.
- Policy extensions: ECB backstop, macroprudential bond tax.
