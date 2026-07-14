# Diagnostic v2 — substitution vs deleveraging (macro-pru-fix)

**Objective (shifted from v1).** The `psi_lambda_B` sweep already shows peak ΔY_D < 0
at every dial value, so the perverse "GDP rises in a crisis" symptom is resolved on
this branch. This diagnostic instead **confirms deleveraging dominates and
characterises whether substitution is still present but dominated**, at the
data-disciplined dial `psi_lambda_B = 2.8` (targets ~150bp 2010 GR spread). The
decisive classifier is the **sign of the capital LEVEL `K_D`**, not the portfolio
share.

**Branch/commit:** `macro-pru-fix` `c6f5707`. `EL_price` fix present
(`req_spread = excess_return + (EL_price_D + psi_spread_D)·def_rate(+1)`).
Calibration: `writeoff_enabled=0`, `mv_rule=0`, `recovery=0`, `zeta=0`
⇒ `EL_price_D = δ_b/q_b ≈ 0.1025`. Diagnose only; no `code/` edits.

**Configs solved (both at `psi_lambda_B = 2.8`, `psi_spread = 0.777747·2.8/3`):**
- `ELon`  — `EL_price` anchored (0.102491): the fixed model.
- `ELoff` — `EL_price = 0`: removes the fundamental expected-loss/net-worth channel
  the fix added; retains the pre-fix agency (`psi_spread`) + IC-collateral channels.
  `K_D(ELon) − K_D(ELoff)` = the deleveraging pull the fix supplies; `K_D(ELoff)` is
  the (impure) substitution-leaning residual.

---

## Step 0 — solve — 2026-07-14 09:19:06
- 2026-07-14 09:25:34 anchors: EL_price_D(anchor)=0.102491; psi_spread_D(base@3)=0.777747; psi_spread@2.8=0.725897
- SS levels: n_inter=3.0000 K_D=10.8000 b_D_D=0.7687 q_b_D=0.9757 theta=4.000 Y_D=1.0000 I_D=0.2420 C_D=0.6763
- SS market clearing: goods_mkt_D=-3.90e-07  goods_mkt_F=-3.90e-07  ca_res_D=0.00e+00  deposit_mkt_D=1.05e-13

### Step 0/1 checks
- 2026-07-14 09:25:34 solving [ELon] psi_lambda_B=2.8, psi_spread=0.725897, EL_price=0.102491
    [ELon] shock reaches q_b: spread_rb impact=3.6805e-03 peak(bp,ann)=147.2; q_b_D impact=-2.8788e-02
    [ELon] Y_D impact=-2.8267e-04 (-0.0283%SS)  peakabs=-3.2282e-04
    [ELon] K_D impact=-1.1450e-03 (-0.0106%SS)  peakabs=-7.0943e-03 (-0.0657%SS)
    [ELon] n_inter_D impact=-3.7815e-02 (-1.2605%SS); b_D_D impact=5.7579e-03 (+0.7491%SS); I_D impact=-1.1450e-03
- 2026-07-14 09:26:21 solving [ELoff] psi_lambda_B=2.8, psi_spread=0.725897, EL_price=0.000000
    [ELoff] shock reaches q_b: spread_rb impact=2.9863e-03 peak(bp,ann)=119.5; q_b_D impact=-2.3690e-02
    [ELoff] Y_D impact=-2.2501e-04 (-0.0225%SS)  peakabs=-2.3154e-04
    [ELoff] K_D impact=-9.7332e-04 (-0.0090%SS)  peakabs=-5.0935e-03 (-0.0472%SS)
    [ELoff] n_inter_D impact=-3.1873e-02 (-1.0624%SS); b_D_D impact=4.7636e-03 (+0.6197%SS); I_D impact=-9.7332e-04

- 2026-07-14 09:27:08 SOLVE v2 COMPLETE. saved irfs_2p8_ELon.npz, irfs_2p8_ELoff.npz, ss_values.json

## Step 1 — aggregate output sign at psi_lambda_B=2.8 (ELon)
- Y_D impact = -2.8267e-04 (-0.0283%SS); extremum(100q) = -3.2282e-04 (-0.0323%SS) at t=26
- Sweep implied ΔY at 2.8 ≈ interpolate(2.6:-0.0271%, 3.0:-0.0389%) ≈ -0.033%SS. Consistent: YES

## Step 2 — bank balance sheet (LEVELS), classify by sign of capital LEVEL K_D
- K_D (capital LEVEL): impact -0.0106%SS, extremum -0.0657%SS at t=36
- b_D_D (Greek sov QUANTITY): impact +0.7491%SS, extremum +1.3519%SS at t=4
- n_inter_D (net worth): impact -1.2605%SS, extremum -1.2605%SS at t=0
- total assets (θ·n): impact -0.3443%SS, min -0.3443%SS
- bond MV share of book: impact Δ -0.1161pp, extremum -0.1161pp (>0 = tilt TOWARD bonds)
- REGIME: DELEVERAGING-DOMINANT: capital LEVEL falls (K_D ≤ 0 throughout)

## Step 3a — ΔY_D decomposition (goods-market identity)
- contribution I_D       : impact -1.1450e-03  extremum -1.1450e-03 at t=0
- contribution P_CES·C_D : impact -3.2524e-04  extremum -8.4937e-04 at t=4
- contribution NX_D      : impact +1.1876e-03  extremum +1.1876e-03 at t=0
- contribution G_D       : impact +0.0000e+00  extremum +0.0000e+00 at t=0
- identity residual max|·| = 5.35e-17 (should be ~0)

## Step 3b — capital: substitution push vs deleveraging pull (EL_price on/off)
- net K_D (ELon)                    : impact -0.0106%SS  extremum -0.0657%SS at t=36
- substitution-leaning (ELoff)      : impact -0.0090%SS  extremum -0.0472%SS at t=32
- deleveraging pull (ELon-ELoff)    : impact -0.0016%SS  extremum -0.0196%SS at t=49

## Step 4 — net-worth-to-spread pass-through at 2.8 (validation moment)
- peak spread = 147.2bp; peak Δn_inter = -1.261%SS; pass-through = -0.856%/100bp  (sweep flag ≈ -0.85)

## Step 6 — empirical-prior consistency
- capital/credit: contracting; Greek-sov quantity b_D_D: mixed/up; bond MV share tilt: toward capital
- Model omits renationalisation motive (always-binding IC). Check it does NOT produce the doubly-counterfactual bonds↓/capital↑ WITH credit expanding.
- doubly-counterfactual (bonds↓ & capital↑ & credit↑)? NO

- analysis complete; figures v2_01/02/03 written.

---

## Conclusions

- **Y_D falls** at psi_lambda_B=2.8 (−0.028% impact / −0.032% trough) — matches sweep. ✓
- **Deleveraging-dominant by the decisive classifier:** capital LEVEL K_D < 0 throughout
  (trough −0.066%SS), net worth −1.26%, total assets −0.34%. Capital never rises. ✓
- **Substitution present but dominated, and it points toward SOVEREIGNS not capital:**
  Greek-sov *quantity* b_D_D ↑ (+0.75→+1.35%) while capital falls; but bond *market-value*
  share ↓ (q_b −2.88% crash dominates), so it is deleveraging in MV with a quantity tilt.
- **ΔY is investment/deleveraging, NOT terms-of-trade:** I_D −1.145e-3 and P·C −0.33e-3
  drive it; NX_D is +1.19e-3 (cushions). The small aggregate masks a larger credit
  contraction offset by rising net exports.
- **EL_price deepens deleveraging** (capital −0.047% off → −0.066% on); substitution toward
  capital never emerges even with EL off. EL_price is additive, not double-counted.
- **Pass-through −0.856%/100bp** at 2.8 (validation moment, matches sweep ≈−0.85). ✓
- **Step 5 not triggered** (capital does not rise); no recommended_fix.md.

## Deliverables (diagnostics/substitution_v2/)
`VERDICT.md`, `run_log.md`, `ss_values.json`, `irfs_2p8_ELon.npz`, `irfs_2p8_ELoff.npz`,
`solve_v2.py`, `analyze_v2.py`, figures `v2_01_Y_decomposition.png`,
`v2_02_balance_sheet.png`, `v2_03_EL_on_off_capital.png`, `env.txt`. No code/ modified.
