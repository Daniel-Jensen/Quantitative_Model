# psi_lambda_B breakdown sweep — broad-scope EBA calibration

Generated 2026-07-31 16:49:53. Steady state solved once (psi_lambda_B enters only via the psi_spread anchor, which is exactly linear in it); Jacobian re-solved per point with BOTH dials moved together.

`n_inter_D_ss = 2.1379`, `lambda_gk_D = 2.2129`, `Omega_D = 10.8515`.

| psi_lambda_B | psi_spread_D | peak spread (bp ann) | b_gov_D[T-1] | n_inter_D[0] (%) | Y_D[0] (%) | verdict |
|---:|---:|---:|---:|---:|---:|:--|
| 8.5 | 1.738 | 150.3 | +1.39e-05 | -3.380 | -0.0149 | ok |
| 10.0 | 2.044 | 172.1 | +4.48e-06 | -3.803 | -0.0164 | ok |
| 12.0 | 2.453 | 196.3 | +1.41e-06 | -4.195 | -0.0170 | ok |
| 14.0 | 2.862 | 223.8 | +1.21e-06 | -4.823 | -0.0193 | ok |
| 15.0 | 3.067 | 225.8 | +2.34e-06 | -4.633 | -0.0171 | ok |
| 16.0 | 3.271 | 228.1 | +2.65e-06 | -4.499 | -0.0153 | ok |
| 17.0 | 3.475 | 229.2 | +2.47e-06 | -4.336 | -0.0135 | ok |
| 18.0 | 3.680 | 234.6 | +1.72e-06 | -4.329 | -0.0128 | ok |
| 20.0 | 4.089 | 273.6 | -9.60e-08 | -5.288 | -0.0180 | ok |
| 22.0 | 4.498 | 343.1 | -1.00e-06 | -7.221 | -0.0298 | ok |
| 25.0 | 5.111 | 625.3 | -2.32e-06 | -15.467 | -0.0826 | ok |
| 26.0 | 5.315 | 1034.5 | -4.08e-06 | -27.716 | -0.1627 | **BREAKDOWN** |
| 27.0 | 5.520 | 8903.8 | -3.85e-05 | -264.922 | -1.7213 | **BREAKDOWN** |
| 28.0 | 5.724 | 41.0 | +4.31e-06 | +29.866 | +0.2164 | **BREAKDOWN** |
| 29.0 | 5.929 | 22.4 | +1.87e-06 | +12.940 | +0.1054 | **BREAKDOWN** |
| 30.0 | 6.133 | 16.6 | +1.12e-06 | +7.702 | +0.0713 | **BREAKDOWN** |

- peak spread monotone increasing up to psi_lambda_B = **27.0**
- first psi_lambda_B at which n_inter_D[0] REVERSES (shrinks while the spread still rises): **15.0**
- first breakdown row (A7 / stationarity / sign): **26.0**

Breakdown criteria: peak spread > 1000bp (A7), |b_gov_D[T-1]| > 1e-2, or a sign flip in n_inter_D[0] / Y_D[0] on a default shock.

`regime_model.PSILAM_BREAKDOWN` is set from the FIRST pathology (the n_inter_D[0] reversal), not from the pole, so the guard has real margin rather than sitting on the edge of the singularity.
