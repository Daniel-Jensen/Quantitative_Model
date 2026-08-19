# Step 4 — CB feedback sign and closed-loop stability

A_cb = d(spread_rb)/d(cb_buy_D), T = 500

## 4a. Sign of the feedback

| lag h | A_cb[h,0]  (impulse at t=0) | A_cb[h,h] (diagonal) |
|---|---|---|
| 0 | -4.397083e-03 | -4.397083e-03 |
| 1 | -1.533000e-03 | -4.169005e-03 |
| 2 | -1.518925e-03 | -3.944466e-03 |
| 3 | -1.437022e-03 | -3.724130e-03 |
| 4 | -1.325046e-03 | -3.519005e-03 |
| 5 | -1.202372e-03 | -3.336355e-03 |
| 6 | -1.079358e-03 | -3.178659e-03 |
| 7 | -9.616160e-04 | -3.045240e-03 |

- `A_cb[0,0] = -4.397083e-03`  -> impact effect of a unit purchase on the impact spread is **NEGATIVE (compresses)**
- column sum `sum_h A_cb[h,0] = -1.597028e-02`  (cumulative spread response to a one-period purchase at t=0)
- fraction of the t=0 column that is negative: 16.4%
- most positive entry in the t=0 column: +1.849e-05 at lag 111
- spectral radius of A_cb = 0.452155; max real eigenvalue = +0.452155; min real eigenvalue = -0.023262
- 1 / max real eigenvalue = +2.2116 (the gamma at which I - gamma*A_cb becomes singular, if the max real eigenvalue is the binding one)

## 4b. Closed-loop pole and margin

- condition-number scan: {'gamma_pole': 26.5, 'gamma_safe_max': 19.875}
- intended gamma range: [0, 10] (code/tpi.py gamma_values = [0, 2, 5, 10])
- margin: pole / gamma_max = 2.65x; the 0.75-safety cap sits at gamma = 19.88

## 4c. Per-gamma stability

| gamma | peak spread (bp ann) | b_gov_D[499] | n_inter_D[0] | Y_D[0] | |lam| spread | |lam| b_gov_D | |lam| n_inter_D |
|---|---|---|---|---|---|---|---|
| 0 | 205.87 | +3.024e-05 | -2.4388e-01 | -1.9742e-02 | 0.939968 | 0.828311 | 0.934720 |
| 2 | 193.25 | +3.694e-04 | -2.1289e-01 | -1.5702e-02 | 0.941372 | 0.841631 | 0.937532 |
| 5 | 176.75 | -5.615e-05 | -1.7260e-01 | -1.0428e-02 | 0.943455 | 0.844431 | 0.941139 |
| 10 | 154.36 | -1.506e-04 | -1.1837e-01 | -3.2741e-03 | 0.946993 | 0.923171 | 0.946025 |

All moduli must be < 1 for a stationary closed loop.

## 4d. Breakdown scan over gamma (open-loop grid, up to the pole)

| gamma | peak spread (bp ann) | compression vs g=0 | cond(I - g A_cb) | |lam| spread |
|---|---|---|---|---|
| 0 | 205.87 | +0.0% | 1.000e+00 | 0.939968 |
| 1 | 199.38 | +3.2% | 3.004e+01 | 0.940676 |
| 2 | 193.25 | +6.1% | 6.461e+02 | 0.941372 |
| 5 | 176.75 | +14.1% | 3.090e+02 | 0.943455 |
| 10 | 154.36 | +25.0% | 4.803e+02 | 0.946993 |
| 15 | 136.73 | +33.6% | 7.833e+02 | 0.951232 |
| 19.88 | 122.90 | +40.3% | 1.406e+03 | 0.958447 |
| 22 | 117.75 | +42.8% | 2.007e+03 | 0.978901 |
| 25 | 111.41 | +45.9% | 4.656e+03 | 0.997033 |
| 26.5 | 109.50 | +46.8% | 1.274e+04 | nan |

## 4e. On-path fiscal incidence of the conduit

conduit series present in the dump: ['G_cb_flow_D__cb_buy_D', 'G_rem_cb_D__cb_buy_D', 'G_rem_cb_F__cb_buy_D', 'irf_g0_cb_flow_D', 'irf_g0_rem_cb_D', 'irf_g0_rem_cb_F', 'irf_g10_cb_flow_D', 'irf_g10_rem_cb_D', 'irf_g10_rem_cb_F', 'irf_g2_cb_flow_D', 'irf_g2_rem_cb_D', 'irf_g2_rem_cb_F', 'irf_g5_cb_flow_D', 'irf_g5_rem_cb_D', 'irf_g5_rem_cb_F']
G columns that do not exist (pure-CB objects have no shock_def_D loading): ['rem_cb_D__shock_def_D', 'rem_cb_F__shock_def_D', 'cb_flow_D__shock_def_D']
- g=2  rem_cb_D  : t0 = -6.688102e-04, t1 = +3.165945e-05, min = -6.688102e-04, max = +3.636726e-05, undiscounted sum(0:100) = +7.724208e-06
- g=2  rem_cb_F  : t0 = -7.466957e-04, t1 = +3.534631e-05, min = -7.466957e-04, max = +4.060237e-05, undiscounted sum(0:100) = +8.623722e-06
- g=2  cb_flow_D : t0 = -9.419861e-03, t1 = +4.459077e-04, min = -9.419861e-03, max = +5.122150e-04, undiscounted sum(0:100) = +1.087917e-04
- g=5  rem_cb_D  : t0 = -1.529293e-03, t1 = +6.941594e-05, min = -1.529293e-03, max = +8.148337e-05, undiscounted sum(0:100) = +9.616795e-06
- g=5  rem_cb_F  : t0 = -1.707385e-03, t1 = +7.749969e-05, min = -1.707385e-03, max = +9.097242e-05, undiscounted sum(0:100) = +1.073671e-05
- g=5  cb_flow_D : t0 = -2.153934e-02, t1 = +9.776893e-04, min = -2.153934e-02, max = +1.147653e-03, undiscounted sum(0:100) = +1.354478e-04
- g=10 rem_cb_D  : t0 = -2.671150e-03, t1 = +1.131731e-04, min = -2.671150e-03, max = +1.378937e-04, undiscounted sum(0:100) = -1.127434e-05
- g=10 rem_cb_F  : t0 = -2.982216e-03, t1 = +1.263526e-04, min = -2.982216e-03, max = +1.539519e-04, undiscounted sum(0:100) = -1.258728e-05
- g=10 cb_flow_D : t0 = -3.762183e-02, t1 = +1.593988e-03, min = -3.762183e-02, max = +1.942164e-03, undiscounted sum(0:100) = -1.587936e-04
- g=2  TAX_D     : t0 = -1.057603e-02, peak|.| = +1.057603e-02
- g=2  TAX_F     : t0 = +4.620592e-04, peak|.| = +4.620592e-04
- g=5  TAX_D     : t0 = -8.136094e-03, peak|.| = +8.136094e-03
- g=5  TAX_F     : t0 = +2.489603e-04, peak|.| = +3.360112e-04
- g=10 TAX_D     : t0 = -4.840994e-03, peak|.| = +4.840994e-03
- g=10 TAX_F     : t0 = -3.985060e-05, peak|.| = +4.535629e-04

## 4f. Walras residuals along the closed loop

| gamma | max|ca_res_D| | max|goods_mkt_D| | max|goods_mkt_F| |
|---|---|---|---|
| 0 | 7.21e-08 | 5.86e-17 | 1.89e-10 |
| 2 | 5.55e-08 | 5.47e-17 | 1.96e-10 |
| 5 | 3.38e-08 | 5.59e-17 | 2.04e-10 |
| 10 | 1.47e-08 | 5.84e-17 | 2.15e-10 |
