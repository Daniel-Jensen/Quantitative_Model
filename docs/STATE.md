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

## Known limitations

- Solver calibration is fragile for portfolio adjustment cost parameters `psi_bF_D` and `psi_bD_F`.
- External habit formation is noted in the research log but not yet integrated.
- ECB/ESM backstop concepts are present in notes only and not incorporated into the active model.
- Macroprudential and long-term bond extensions are conceptual, not fully implemented.

## Next focus

- Stabilize the steady-state solve for `model_v12`.
- Verify impulse responses for sovereign spread and cross-border portfolio shocks.
- Document the exact meaning of key calibration scalars and policy parameters.
