# Diagnostic Run Log — `psi_lambda = 0` null-result investigation

**Purpose.** Determine why a sovereign default shock produces null IRFs when
`psi_lambda_B = 0`, and classify the cause as Case 1 (dead mark-to-market /
stale price), Case 2 (parameter conflation), Case 3 (no fundamental channel), or
a documented alternative. **Diagnose only — no fixes applied to `code/`.**

Model git commit: `f25dcf6` (branch `market-value-debt`). Full env in
[`env.txt`](env.txt). Interpreter `/opt/anaconda3/envs/ssj/bin/python`
(Python 3.12.12, numpy 1.26.4, scipy 1.17.1, sequence_jacobian import OK).

---

## Step 1 — Static classification of `psi_lambda` (targets Case 2)

Every occurrence of `psi_lambda_B` / its derived anchor `psi_spread` in `code/`
(`grep -rn "psi_lambda_B\|psi_spread"`), with role classification:

| File:line | Context | Role |
|-----------|---------|------|
| `calibration.py:38` | `psi_lambda_B_D/F = 3.0` | **parameter definition** |
| `steady_state.py:37-40` | `psi_spread_{D,F} = lambda_gk·psi_lambda_B/(beta_inter·Omega)` | **derived anchor** (converts `psi_lambda_B` → a required-spread loading; STRUCTURAL) |
| `equations_D.py:282-283` | `intermediation_IC_D`: `Delta_bD_eff = Delta_bD_D + psi_lambda_B_D·def_rate_D(+1)` | **STRUCTURAL** — bank IC / divertability |
| `equations_D.py:432` | `divert_bond_foc_D`: `req_spread = excess_return_bD_D_ss + psi_spread_D·def_rate_D(+1)` | **STRUCTURAL** — bank bond valuation FOC |
| `equations_F.py:247-248` | `intermediation_IC_F` (symmetric) | **STRUCTURAL** — bank IC |
| `equations_F.py:398` | `divert_bond_foc_F` (symmetric) | **STRUCTURAL** — bank bond FOC |
| `equations_global.py:86,93` | `divert_portfolio_adj`: `prem_FD/prem_DF = excess_return_ss + psi_spread·def_rate(+1)` | **STRUCTURAL** — cross-border bond FOC |

**Finding (Step 1).** `psi_lambda_B` occurs **only inside structural equations**
(the GK incentive constraint and the sovereign-bond optimality conditions).
There is **no macroprudential policy rule that references `psi_lambda_B`.** The
model's actual macroprudential instrument is a *separate* bond tax
`tau_mp = T0 + T1·def_rate` (`macro_pru_tax_D/F`), which uses `T0/T1`, not
`psi_lambda_B`, and is independently off (`T1 = 0`).

⇒ The premise "`psi_lambda` is the macroprudential policy-response coefficient"
is a **mislabel**. `psi_lambda_B` is the *structural* state-dependent
divertability / risk-premium-sensitivity parameter. Because it is the **sole**
carrier of sovereign-risk sensitivity in the baseline calibration (see Step 3),
zeroing it removes fundamental risk sensitivity — the *effect* Case 2 predicts,
but via a structural, not a policy, coefficient. Detailed reconciliation in
`VERDICT.md`.

---

## Step 3 (static) — Does the spread / default rate enter any structural equation? (targets Case 3)

Every consumer of `def_rate_{D,F}` in the **dynamic** model (`ha_full`, built in
`full_model.py`) — from `grep -rn "def_rate" code/*.py` — with the switch that
gates it and whether it survives the baseline calibration
(`writeoff_enabled = 0`, `T1 = 0`) and `psi_lambda_B = 0`:

| Channel | Block | Gate | Baseline (wo=0,T1=0) | + `psi_lambda_B=0` |
|---------|-------|------|:---:|:---:|
| Bond cash-flow (coupon haircut) | `bond_return_D/F` | `writeoff_enabled` | **DEAD** | dead |
| Gov budget realized loss | `budget_residual_D/F` | `writeoff_enabled` | **DEAD** | dead |
| Macropru tax (the real policy) | `macro_pru_tax_D/F` | `T1` | **DEAD** | dead |
| IC collateral divertability | `intermediation_IC_D/F` | `psi_lambda_B` | live | **DEAD** |
| Bond FOC required spread | `divert_bond_foc_D/F` | `psi_spread ∝ psi_lambda_B` | live | **DEAD** |
| Cross-border FOC premium | `divert_portfolio_adj` | `psi_spread ∝ psi_lambda_B` | live | **DEAD** |

(Steady-state-only consumers `bond_price_ss_*`, `government_ss_*`, `smart_steady_*`
are irrelevant to IRFs and see `def_rate = 0` at SS.)

**Finding (Step 3, static).** In the baseline calibration the *only* live
default-risk channels are the three `psi_lambda_B`-scaled structural terms.
Setting `psi_lambda_B = 0` zeroes all three simultaneously, after which
`def_rate` (and hence `shock_def_D`) enters **no** structural equation ⇒ the
sovereign default shock is inert. This is the Case 3 signature, established
statically; the empirical solve (Step 0/2/5 below) confirms it.

## Step 2 (static, Case 1 pre-check) — Is the mark-to-market plumbing endogenous or stale?

Inspection of the net-worth / bond-return blocks:

- `bond_return_D` (`equations_D.py:211-218`): `rb_actual_D = (coupon + (1-δ_b)·q_b_D)/q_b_D(-1) - 1` — uses the **endogenous current** price `q_b_D` and its lag. Not stale.
- `bank_return_D` (`equations_D.py:293-303`): `rn_D` built from `phi_*_lag = q_b_*(-1)·b(-1)/n(-1)` and `rb_actual_*` — **endogenous** revaluation.
- `k_balance_sheet_D` (`equations_D.py:320-322`): `K_res = Q·K + q_b_D·b_D_D + q_b_F·b_F_D - θ·n` — net worth valued at **current endogenous** `q_b`.

⇒ The revaluation channel is wired with **endogenous** prices, so a move in
`q_b` *would* move `n_inter`. Case 1 (stale/SS price feeding net worth) is
**not** supported by code inspection. The decisive empirical discriminator is
whether `q_b` itself moves at `psi_lambda_B = 0` (Case 3: it will not; Case 1: it
would). Recorded by the solve below.

---

---

## Solve run — 2026-07-13 22:43:52
- `2026-07-13 22:43:52` Importing pipeline modules from code/ (unmodified).
- `2026-07-13 22:45:32` Step 1: get_calibration()
- `2026-07-13 22:45:32` psi_lambda_B_D (baseline) = 3.0, writeoff_enabled_D = 0.0, T1_D = 0.0, mv_rule_D = 0.0
- `2026-07-13 22:45:32` Step 2: solve_steady_state()
- `2026-07-13 22:46:40` Step 3: calibrate_ic_delta()
- `2026-07-13 22:46:40` Step 4: calibrate_depreciation()  (final SS re-solve)

### SS sanity checks (market clearing + key SS values)
    goods_mkt_D      = -3.897868541408167e-07
    goods_mkt_F      = -3.8978685477919495e-07
    ca_res_D         = 0.0
    deposit_mkt_D    = 1.0480505352461478e-13
    deposit_mkt_F    = -1.7763568394002505e-15
    psi_lambda_B_D   = 3.0
    psi_lambda_B_F   = 3.0
    psi_spread_D     = 0.777747203781219
    psi_spread_F     = 0.777747203781219
    q_b_D            = 0.9756984802086146
    q_b_F            = 0.9756984802086146
    rb_D             = 0.0024906792707301896
    rb_F             = 0.0024906792707301896
    n_inter_D        = 3.0
    theta_D          = 4.0
    beta_D           = 0.9994621763351645
    beta_F           = 0.9994621763351645
    p                = 0.9999999999999996
    spread_rb (SS)   = 0.0   (expect ~0)

- `2026-07-13 22:47:49` Building baseline Jacobian G (psi_lambda_B = 3.0) via build_and_solve()...
- `2026-07-13 22:48:40` Baseline G computed. T=500. shock peaks: def=0.01, TFP=0.01

- `2026-07-13 22:48:40` Building psi_lambda_B = 0 Jacobian G0 (zero psi_lambda_B_* and psi_spread_*)...
- `2026-07-13 22:48:40` Confirming ss0 overrides: psi_lambda_B_D=0.0, psi_lambda_B_F=0.0, psi_spread_D=0.0, psi_spread_F=0.0
- `2026-07-13 22:49:38` G0 computed.
- `2026-07-13 22:49:39` Saved irfs_baseline.npz (42 series), irfs_psilam0.npz (42 series).
- `2026-07-13 22:49:39` Wrote summary.json

### Default shock (1pp): baseline vs psi_lambda_B=0

| var | impact (base) | peak|·| (base) | impact (ψλ=0) | peak|·| (ψλ=0) | peak ratio ψλ0/base |
|-----|--------------:|---------------:|--------------:|---------------:|--------------------:|
| def_rate_D | 1.0000e-02 | 1.0000e-02 | 1.0000e-02 | 1.0000e-02 | 1.000e+00 |
| q_b_D | -2.5980e-02 | 2.5980e-02 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| q_b_F | 5.5558e-03 | 5.5558e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| rb_actual_D | -2.3965e-02 | 2.3965e-02 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| rb_actual_F | 5.1248e-03 | 5.1248e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| spread_rb | 3.3127e-03 | 3.3127e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| n_inter_D | -3.5198e-02 | 3.5198e-02 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| n_inter_F | 5.2398e-03 | 7.5425e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| theta_D | 3.4081e-02 | 3.4081e-02 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| theta_F | -5.2042e-03 | 8.1125e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| K_D | -1.0864e-03 | 6.0273e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| I_D | -1.0864e-03 | 1.0864e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| Y_D | -2.5400e-04 | 2.7369e-04 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| Y_F | 2.4404e-04 | 3.5603e-04 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| C_D | -7.2250e-04 | 1.0852e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| C_F | 9.5321e-04 | 9.5321e-04 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| w_D | 5.3693e-05 | 1.0307e-04 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| b_gov_D | 4.7122e-03 | 9.5926e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| TAX_D | -1.4663e-03 | 1.7578e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| rdep_D | -1.2756e-03 | 1.2756e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |
| div_D | -4.7997e-03 | 4.7997e-03 | 0.0000e+00 | 0.0000e+00 | 0.000e+00 |

### TFP shock (1%, CONTROL): baseline vs psi_lambda_B=0

| var | impact (base) | peak|·| (base) | impact (ψλ=0) | peak|·| (ψλ=0) | peak ratio ψλ0/base |
|-----|--------------:|---------------:|--------------:|---------------:|--------------------:|
| def_rate_D | 0.0000e+00 | 6.1742e-03 | 0.0000e+00 | 4.2016e-03 | 6.805e-01 |
| q_b_D | 3.9591e-02 | 4.4693e-02 | 7.6810e-03 | 9.3316e-03 | 2.088e-01 |
| q_b_F | 3.3835e-04 | 6.0795e-03 | 4.4677e-03 | 5.5280e-03 | 9.093e-01 |
| rb_actual_D | 3.6519e-02 | 3.6519e-02 | 7.0851e-03 | 7.0851e-03 | 1.940e-01 |
| rb_actual_F | 3.1210e-04 | 1.2998e-03 | 4.1211e-03 | 4.1211e-03 | 3.171e+00 |
| spread_rb | -4.1232e-03 | 4.4972e-03 | -3.3753e-04 | 4.0629e-04 | 9.034e-02 |
| n_inter_D | 4.1124e-01 | 4.1124e-01 | 3.6622e-01 | 3.6622e-01 | 8.905e-01 |
| n_inter_F | -1.2872e-02 | 1.2872e-02 | -3.6783e-03 | 3.6783e-03 | 2.858e-01 |
| theta_D | -3.9704e-01 | 3.9704e-01 | -3.5328e-01 | 3.5328e-01 | 8.898e-01 |
| theta_F | 1.2660e-02 | 1.2660e-02 | 3.5172e-03 | 4.2485e-03 | 3.356e-01 |
| K_D | 1.9236e-02 | 7.0512e-02 | 1.7862e-02 | 6.4213e-02 | 9.107e-01 |
| I_D | 1.9236e-02 | 1.9236e-02 | 1.7862e-02 | 1.7862e-02 | 9.286e-01 |
| Y_D | 3.0152e-02 | 3.0152e-02 | 2.9755e-02 | 2.9755e-02 | 9.868e-01 |
| Y_F | -5.3908e-05 | 8.4791e-04 | 3.4282e-04 | 3.6013e-04 | 4.247e-01 |
| C_D | 1.2932e-02 | 1.2932e-02 | 1.1298e-02 | 1.1298e-02 | 8.736e-01 |
| C_F | -6.7412e-04 | 1.3751e-03 | 9.6461e-04 | 1.0424e-03 | 7.581e-01 |
| w_D | 1.3069e-02 | 1.3069e-02 | 1.3163e-02 | 1.3163e-02 | 1.007e+00 |
| b_gov_D | -2.5629e-02 | 5.5880e-02 | -1.9722e-02 | 3.8027e-02 | 6.805e-01 |
| TAX_D | 2.0135e-02 | 2.0135e-02 | 1.8193e-02 | 1.8193e-02 | 9.035e-01 |
| rdep_D | 1.8351e-03 | 1.8351e-03 | 6.6663e-04 | 1.8015e-03 | 9.817e-01 |
| div_D | 5.6079e-02 | 5.6079e-02 | 4.9940e-02 | 4.9940e-02 | 8.905e-01 |

- `2026-07-13 22:49:39` SOLVE RUN COMPLETE.

---

## Conclusions (Steps 0/2/4/5)

- **Step 0 — reproduction:** CONFIRMED. Default shock at `psi_lambda_B=0` →
  every endogenous variable identically `0`; baseline is nonzero and
  correctly-signed (`n_inter_D[0]=−3.52%`, `Y_D[0]<0`, spread widens).
- **Step 2 — decisive `n_t` test:** `n_inter_D` response at `psi_lambda_B=0` is
  **exactly 0** (baseline −3.52%), *together with* `q_b_D` (exactly 0). The bond
  price is a live endogenous variable (moves 7.68e−3 on the TFP control) → the
  break is *upstream* of the price (no forcing), i.e. **Case 3, not Case 1**.
- **Step 4 — disambiguation:** the revaluation term connecting `def_rate`→`q_b`
  →net worth is present and endogenous, but receives no `def_rate` forcing when
  `psi_lambda_B=0`; there is **no** fundamental-loss term feeding it. Case 3.
- **Step 5 — sanity:** shock reaches the driver (`def_rate_D`=0.01 both configs);
  TFP control transmits normally at `psi_lambda_B=0` (`Y_D` 98.7% of baseline);
  SS market clearing holds (`goods_mkt_D`=−3.9e−7, `ca_res_D`=0, deposits ~1e−13);
  `spread_rb(SS)=0`.

**Verdict:** Case 3 (no fundamental channel independent of `psi_lambda_B`), with
Case 2 reframed as a mislabel (`psi_lambda_B` is structural, not policy) and
Case 1 ruled out. Deeper cause: `writeoff_enabled` conflates expected-loss
*pricing* with realized *booking*. Full write-up: [`VERDICT.md`](VERDICT.md).

## Deliverables

| File | Contents |
|------|----------|
| [`VERDICT.md`](VERDICT.md) | Headline verdict + case-by-case adjudication with evidence |
| [`recommended_fix.md`](recommended_fix.md) | Proposed (unimplemented) fix: fundamental expected-loss loading in the bond FOC; exact equations/lines, SS side-effects, re-verification checklist |
| [`run_log.md`](run_log.md) | This timestamped log (static + empirical) |
| [`summary.json`](summary.json) | Machine-readable SS checks + impact/peak stats, both configs, both shocks |
| [`irfs_baseline.npz`](irfs_baseline.npz), [`irfs_psilam0.npz`](irfs_psilam0.npz) | Full IRF series (42 each) for re-analysis without re-solve |
| `01_default_transmission_chain.png` | Chain grid, default shock, baseline vs ψλ=0 |
| `02_decisive_qb_networth.png` | Decisive `q_b`/`n_inter` panel (Case 3 vs Case 1) |
| `03_tfp_control.png` | TFP control — model transmits normally at ψλ=0 |
| [`solve_configs.py`](solve_configs.py), [`make_figures.py`](make_figures.py) | Diagnostic scripts (import unmodified `code/`; write only to `diagnostics/`) |
| [`env.txt`](env.txt) | Environment + git commit `f25dcf6` |

**No files under `code/`, `routines/`, or the notebook were modified.**

---

