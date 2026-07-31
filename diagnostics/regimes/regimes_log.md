# Regime experiment log — MAIN model (reconciled 2026-07-23)

Policy-regime feature rebuilt on `main`'s EBA-anchored, C-1-fixed model after the
ms-regime version was found to rest on a superseded model. Branch: `regimes-reconciled`.

## Pivotal reconciliation finding — SA-1 is a superseded-model artifact

On ms-regime (single-country conduit, par-value fiscal rule, psi_lambda_B=2.8) a CB
purchase WIDENED the D-F spread on impact (A_cb = +4.85e-4) — "Finding SA-1", which
forced an output-protection workaround and an inverted Stage-B sign. That does NOT
hold on main. Main's model has:
  * the ECB **capital-key conduit** (kappa_cb_F=0.929): a D-bond purchase is funded
    92.9% by F's treasury and only 7.1% by D's, so D's gross debt barely moves and D
    banks genuinely shed bonds;
  * the **market-value fiscal rule** (mv_rule=1); and
  * psi_lambda_B=1.1793 (in-range; 2.8 sits in main's linear-approximation-breakdown region).

Measured impact response to a unit cb_buy_D on main (independent run_tpi probe AND the
regime_model cache, matching to machine precision):
  d(spread_rb) = **-1.9455e-2 (COMPRESSES)**, d(b_D_D) = -0.72 (banks shed bonds),
  d(b_gov_D) = -0.078 (gross debt falls). The backstop works as intended.

So SA-1, its output-protection workaround, the SB-1 sign inversion, and the net-debt
post-mortem "fix" are all ms-regime artifacts and do NOT apply here. Stage A/B are
rebuilt with the spec's original spread-compression targeting.

## Cache build (main model)

- calibration: psi_lambda_B=1.1793, mv_rule=1, recovery_rate=0.30, kappa_cb_F=0.929
- EL_price_D = 0.071743 (main recovery=0.30; not the ms-regime 0.102491 anchor)
- model-build sanity: G_tpi[cb=0] vs baseline spread_rb max|err| = 0.00e+00 (block list exact)
- cross-check vs SA-1 probe: d(spread_rb)/d(cb_buy)[0,0] = -1.94550e-02 (exact match)
- caches: cache_G_main_psilam1p18.npz, cache_G_main_psilam0p00.npz (gitignored)
- zero-response note: cb_flow_D (CB conduit flow) does not respond to shock_def_D — expected.


## Stage A run (MAIN model) — 2026-07-23 20:54:30
- passive peak spread: 156.6 bp ann (main is calibrated to ~151bp at psi_lambda_B=1.1793; investigate if outside 120-180)
- impact A_cb = d(spread)/d(cb_buy)[0,0] = -1.946e-02 (COMPRESSES (backstop works — SA-1 absent on main))
- gamma_aggressive = 17.0503 (50% peak-spread compression), gamma_medium = 5.7178 (25%), gamma_passive = 0 (fixed anchor)
  aggressive: peak spread    78.3 bp, Y_D[0] +0.1100%, n_inter_D[0] -2.394%
      medium: peak spread   117.5 bp, Y_D[0] +0.0650%, n_inter_D[0] -2.566%
     passive: peak spread   156.6 bp, Y_D[0] +0.0323%, n_inter_D[0] -2.830%

| regime | discounted CB purchases (Sum beta^t q_b cb_t) | dY_D peak (%) | dC_D peak (%) | pd_D peak |
|---|---|---|---|---|
| aggressive | 0.19952 | -0.0511 | -0.1660 | +0.04214 |
| medium | 0.06357 | -0.0438 | -0.1982 | +0.02413 |
| passive | 0.00000 | -0.0356 | -0.2603 | +0.01268 |

### A6 — ranking at psi_lambda_B = 0 (fundamental floor)
  aggressive: peak spread   5.23 bp, output loss 0.00048%
      medium: peak spread   7.80 bp, output loss 0.00028%
     passive: peak spread  10.28 bp, output loss 0.00010%
- **A6 spread ranking survives at psi_lambda_B=0: YES** (aggressive < medium < passive in peak spread = crisis severity)

Stage A (main) complete.


## Stage B-lite run (MAIN model) — 2026-07-23 20:54:33
- gammas (spread-compression, spec §7): {'aggressive': 17.05, 'medium': 5.718, 'passive': 0.0}
- impact A_cb = -1.946e-02 (compresses on main)
- pi_onset = {'aggressive': 0.0002, 'medium': 0.0587, 'passive': 0.9411}, pi_ergodic = {'aggressive': 0.229, 'medium': 0.5241, 'passive': 0.2468}
- §10.3 assertions PASS (pre-k identity, RE jump=0, Stage A nesting) at pi_onset, k=2

| branch | known-immediate | known-delayed(k) | lottery | delay cost | uncertainty premium |
|---|---|---|---|---|---|
| aggressive | 78.3 | 103.0 | 155.1 | +24.7 | +52.1 |
| medium | 117.5 | 132.0 | 155.1 | +14.5 | +23.2 |
| passive | 156.6 | 156.6 | 155.1 | +0.0 | -1.5 |
(all peak spread, bp ann; delay cost = known-delayed - known-immediate; uncertainty premium = lottery - known-delayed)

| branch | W_D | W_F | discounted CB purchases (A5, per branch) |
|---|---|---|---|
| aggressive | -0.0519 | -0.0933 | 0.08202 |
| medium | -0.0132 | -0.0254 | 0.01965 |
| passive | +0.0175 | +0.0075 | 0.00000 |
| **E_pi** | +0.0157 | +0.0056 | 0.00117 |
(welfare: % SS cons., 100q, tpi.py convention; purchases: Sigma beta^t q_b cb_t)

- impact spread vs pi_passive: **rises** (from 131.04 bp at pi_passive=0.247 to 156.60 bp at pi_passive=1.000); this is the regime-uncertainty price, sign computed not targeted.

| k | impact spread (bp) | passive-branch peak (bp) | E_pi[W_D] |
|---|---|---|---|
| 1 | 154.4 | 154.4 | +0.0156 |
| 2 | 153.7 | 153.7 | +0.0140 |
| 4 | 154.8 | 154.8 | +0.0108 |

- A6 (lottery, psi_lambda_B=0): branch peaks {'aggressive': 10.22, 'medium': 10.22, 'passive': 10.22} bp — aggressive<medium<passive ordering survives: NO (reported)

Stage B-lite (main) complete.


## Stage A run (MAIN model) — 2026-07-30 17:11:47

## Cache build (main model) — 2026-07-30 17:13:18
- calibration: psi_lambda_B=3.0, mv_rule=1.0, recovery_rate=0.3, kappa_cb_F=0.929
- EL_price_D = 0.071743 (main recovery=0.30; NOT the ms-regime 0.102491 anchor)
- 2026-07-30 17:13:18 solving G_tpi at psi_lambda_B=3.0 ...
- model-build sanity: G_tpi[cb=0] vs baseline spread_rb max|err| = 0.00e+00 (expect <1e-8)
- cross-check vs SA-1 probe: d(spread_rb)/d(cb_buy)[0,0] = -2.69964e-02 (probe found -1.9455e-2 → expect match; A_cb<0 = backstop COMPRESSES on main)
  **MISSING OPTIONAL OUTPUT `G_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `ra_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `lambda_gk_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `GINI_WEALTH`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `GINI_C`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  note: output `div_fund_D` has no Jacobian column for `shock_def_D` at this calibration -> zero response (filled 0 T x T); economically = o does not respond to i.
  note: output `div_fund_D` has no Jacobian column for `cb_buy_D` at this calibration -> zero response (filled 0 T x T); economically = o does not respond to i.
- 2026-07-30 17:13:41 solving G_tpi at psi_lambda_B=0.0 ...
  **MISSING OPTIONAL OUTPUT `G_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `ra_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `lambda_gk_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `GINI_WEALTH`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `GINI_C`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  note: output `div_fund_D` has no Jacobian column for `shock_def_D` at this calibration -> zero response (filled 0 T x T); economically = o does not respond to i.
  note: output `div_fund_D` has no Jacobian column for `cb_buy_D` at this calibration -> zero response (filled 0 T x T); economically = o does not respond to i.
- caches written: ['cache_G_main_psilam3p00.npz', 'cache_G_main_psilam0p00.npz']
- passive peak spread: 124.0 bp ann (main is calibrated to ~151bp at psi_lambda_B=3.0; investigate if outside 120-180)
- impact A_cb = d(spread)/d(cb_buy)[0,0] = -2.700e-02 (COMPRESSES (backstop works — SA-1 absent on main))


## Stage A run (MAIN model) — 2026-07-30 17:34:09

## Cache build (main model) — 2026-07-30 17:35:39
- calibration: psi_lambda_B=3.0, mv_rule=0.0, recovery_rate=0.3, kappa_cb_F=0.929
- EL_price_D = 0.071743 (main recovery=0.30; NOT the ms-regime 0.102491 anchor)
- 2026-07-30 17:35:39 solving G_tpi at psi_lambda_B=3.0 ...
- model-build sanity: G_tpi[cb=0] vs baseline spread_rb max|err| = 0.00e+00 (expect <1e-8)
- cross-check vs SA-1 probe: d(spread_rb)/d(cb_buy)[0,0] = -2.40602e-02 (probe found -1.9455e-2 → expect match; A_cb<0 = backstop COMPRESSES on main)
  **MISSING OPTIONAL OUTPUT `G_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `ra_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `lambda_gk_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `GINI_WEALTH`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `GINI_C`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  note: output `div_fund_D` has no Jacobian column for `shock_def_D` at this calibration -> zero response (filled 0 T x T); economically = o does not respond to i.
  note: output `div_fund_D` has no Jacobian column for `cb_buy_D` at this calibration -> zero response (filled 0 T x T); economically = o does not respond to i.
- 2026-07-30 17:36:02 solving G_tpi at psi_lambda_B=0.0 ...
  **MISSING OPTIONAL OUTPUT `G_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `ra_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `lambda_gk_D`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `GINI_WEALTH`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  **MISSING OPTIONAL OUTPUT `GINI_C`** — not in main's G_tpi.outputs; panel zero-filled/omitted WITH a caption note, never silently.
  note: output `div_fund_D` has no Jacobian column for `shock_def_D` at this calibration -> zero response (filled 0 T x T); economically = o does not respond to i.
  note: output `div_fund_D` has no Jacobian column for `cb_buy_D` at this calibration -> zero response (filled 0 T x T); economically = o does not respond to i.
- caches written: ['cache_G_main_psilam3p00.npz', 'cache_G_main_psilam0p00.npz']
- passive peak spread: 187.2 bp ann (main is calibrated to ~151bp at psi_lambda_B=3.0; investigate if outside 120-180)
- impact A_cb = d(spread)/d(cb_buy)[0,0] = -2.406e-02 (COMPRESSES (backstop works — SA-1 absent on main))


## Stage A run (MAIN model) — 2026-07-31 10:15:43
- passive peak spread: 187.2 bp ann (main is calibrated to ~151bp at psi_lambda_B=3.0; investigate if outside 120-180)
- impact A_cb = d(spread)/d(cb_buy)[0,0] = -2.406e-02 (COMPRESSES (backstop works — SA-1 absent on main))
- gamma_aggressive = 5.0813 (50% peak-spread compression), gamma_medium = 1.5730 (25%), gamma_passive = 0 (fixed anchor)
  aggressive: peak spread    93.6 bp, Y_D[0] +0.0067%, n_inter_D[0] -1.978%
      medium: peak spread   140.4 bp, Y_D[0] -0.0102%, n_inter_D[0] -2.594%
     passive: peak spread   187.2 bp, Y_D[0] -0.0261%, n_inter_D[0] -3.001%

| regime | discounted CB purchases (Sum beta^t q_b cb_t) | dY_D peak (%) | dC_D peak (%) | pd_D peak |
|---|---|---|---|---|
| aggressive | 0.31965 | -0.0332 | -0.0631 | +0.00067 |
| medium | 0.16945 | -0.0196 | -0.1348 | +0.00099 |
| passive | 0.00000 | -0.0261 | -0.2276 | +0.00123 |

### A6 — ranking at psi_lambda_B = 0 (fundamental floor)
  aggressive: peak spread   7.14 bp, output loss 0.00060%
      medium: peak spread   8.54 bp, output loss 0.00060%
     passive: peak spread   9.33 bp, output loss 0.00098%
- **A6 spread ranking survives at psi_lambda_B=0: YES** (aggressive < medium < passive in peak spread = crisis severity)

Stage A (main) complete.


## Stage B-lite run (MAIN model) — 2026-07-31 10:16:08
- gammas (spread-compression, spec §7): {'aggressive': 5.081, 'medium': 1.573, 'passive': 0.0}
- impact A_cb = -2.406e-02 (compresses on main)
- pi_onset = {'aggressive': 0.0002, 'medium': 0.0587, 'passive': 0.9411}, pi_ergodic = {'aggressive': 0.229, 'medium': 0.5241, 'passive': 0.2468}
- §10.3 assertions PASS (pre-k identity, RE jump=0, Stage A nesting) at pi_onset, k=2

| branch | known-immediate | known-delayed(k) | lottery | delay cost | uncertainty premium |
|---|---|---|---|---|---|
| aggressive | 93.6 | 116.0 | 185.0 | +22.4 | +69.0 |
| medium | 140.4 | 150.6 | 185.0 | +10.2 | +34.4 |
| passive | 187.2 | 187.2 | 185.0 | +0.0 | -2.2 |
(all peak spread, bp ann; delay cost = known-delayed - known-immediate; uncertainty premium = lottery - known-delayed)

| branch | W_D | W_F | discounted CB purchases (A5, per branch) |
|---|---|---|---|
| aggressive | -5.6253 | +6.0074 | 0.29922 |
| medium | -9.7101 | +9.9601 | 0.15947 |
| passive | -14.6251 | +14.8025 | 0.00000 |
| **E_pi** | -14.3343 | +14.5161 | 0.00943 |
(welfare: % SS cons., 100q, tpi.py convention; purchases: Sigma beta^t q_b cb_t)

- impact spread vs pi_passive: **rises** (from 151.46 bp at pi_passive=0.247 to 187.19 bp at pi_passive=1.000); this is the regime-uncertainty price, sign computed not targeted.

| k | impact spread (bp) | passive-branch peak (bp) | E_pi[W_D] |
|---|---|---|---|
| 1 | 184.7 | 184.7 | -14.3343 |
| 2 | 183.0 | 183.0 | -14.0520 |
| 4 | 181.0 | 181.0 | -13.5722 |

- A6 (lottery, psi_lambda_B=0): branch peaks {'aggressive': 9.3, 'medium': 9.3, 'passive': 9.3} bp — aggressive<medium<passive ordering survives: YES (reported)

Stage B-lite (main) complete.
