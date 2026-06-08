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
- **Corrected EBA bilateral sovereign exposure calibration**: Greek banks (D) and German banks (F) are anchored to observed sovereign asset shares (bond-to-assets, not bond-to-equity; system aggregates, not single banks).
  - All four shares now match EBA exactly: b_D_D = 0.2447, b_F_F = 0.2579, b_F_D = 0.0018, b_D_F = 0.0065

## Known limitations

- Bond supplies dropped 47-63% after correcting for unit mismatch; debt-to-output ratio changes significantly.
- Bank leverage `theta_D = 4` may need adjustment now that sovereign-to-assets shares are realistic.
- Solver calibration is fragile for portfolio adjustment cost parameters `psi_bF_D` and `psi_bD_F`, though corrected calibration shows improved stability.
- External habit formation is noted in the research log but not yet integrated.
- ECB/ESM backstop concepts are present in notes only and not incorporated into the active model.
- Macroprudential and long-term bond extensions are conceptual, not fully implemented.

## Next focus

- Verify that corrected bond supplies (B_supply_D ≈ 3.0, B_supply_F ≈ 3.1) reflect realistic sovereign debt-to-output ratios.
- Check whether bank leverage `theta_D = 4` is still appropriate given the lower sovereign-to-assets shares.
- Re-calibrate portfolio adjustment costs `psi_bF_D` and `psi_bD_F` for the smaller cross-border positions.
- Verify impulse responses for sovereign spread and cross-border portfolio shocks with the corrected calibration.
- Document the exact meaning of key calibration scalars and policy parameters.
- Consider whether the corrected calibration stabilizes the solver fragility noted in earlier iterations.
