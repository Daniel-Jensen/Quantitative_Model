# Central Bank block audit — run log

Audit-only. No model source was edited. All artefacts under `diagnostics/cb_audit/`.

Environment: `/opt/anaconda3/envs/ssj/bin/python`.
Repo state at audit start: branch `gk-structural-foc`, HEAD `ea23e94`
("GK refactor stage 1: bounded, exported pledgeability map"), with **uncommitted
working-tree changes** across `code/calibration.py`, `code/depreciation_calibration.py`,
`code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py`,
`code/full_model.py`, `code/steady_state.py`, `code/test_nkpc_blocks.py`,
`code/tpi.py`, `code/tpi_plots.py` (773 insertions / 239 deletions). The audit is of
the **working tree**, not of `ea23e94`.

---

## 2026-08-19 08:12 — Probe A: brief correction to the audit's stated premise

The task describes the refactor as *"`EL_price` is now integrated into the GK incentive
constraint itself"* and as possibly leaving *"no CB balance sheet and no capital key"*.
Both are out of date; recorded here because two of the audit steps are framed around
them and the framing has to be re-pointed before the findings mean anything.

1. **`EL_price_D/F` (the anchored expected-loss loading) is DELETED, not moved into the
   IC.** Expected default loss now enters through `rb_exp_D` / `rb_exp_F` — the
   expected payoff `bond_return_D/F` computes over the default distribution — which is
   read by `intermediation_P1_D/F`, the banker's **Bellman envelope**, not by
   `intermediation_IC_D/F`. Source: `code/equations_D.py:281` (`bond_return_D`
   docstring: "the model's SINGLE SOURCE OF TRUTH for default loss") and
   `code/equations_D.py:576` (`intermediation_P1_D`: "EXPECTED DEFAULT LOSS ENTERS
   HERE, AND ONLY HERE").
2. **The IC does carry a separate, differently-motivated default channel** —
   `Delta_bD_eff_D` / `Delta_bF_eff_D` from the new `collateral_quality_D`
   (`code/equations_D.py:468`), a bounded pledgeability map with local slope
   `psi_lambda_B_D`. That is a *collateral haircut* channel, not `EL_price`. At the
   preferred baseline `psi_lambda_B_D/F = 0` (`code/calibration.py:288-289`) so this
   channel is **switched off** and `Delta_eff == Delta` identically.
3. **There is a CB balance sheet and there is a capital key.** `code/tpi.py:62-97`
   implements an explicit ECB conduit: `cb_flow_D` (net cash flow on the CB's D-bond
   book) split `kappa_cb_F` to the F treasury and `1 - kappa_cb_F` to the D treasury,
   with `kappa_cb_F = 0.929` (`code/calibration.py:438`).

Steps 1–3 are therefore run against what the code actually does, and the divergence
from the brief is itself logged as a finding in `VERDICT.md`.

---

## 2026-08-19 08:14 — Probe B: `psi_spread` reference sweep (Step 3)

```
grep -rn "psi_spread" --include="*.py" .
```

Classified by whether the hit is **live code** or **prose**:

| File | Lines | Kind | Live? |
|---|---|---|---|
| `code/equations_D.py` | 299, 703, 706 | docstring (records the deleted block's algebra) | no |
| `code/equations_F.py` | 511 | docstring | no |
| `code/equations_global.py` | 119, 127 | docstring | no |
| `code/steady_state.py` | 202, 310 | comment | no |
| `code/calibration.py` | 212–271, 398–402 | comment (tuning history, marked void) | no |
| `code/tpi.py` | 313 | comment (records the removed closed form) | no |
| `code/test_nkpc_blocks.py` | 249, 261, 270, 292 | the AST ban-list itself | n/a |
| `experiments/e3_writeoff_s1.py` | 36, 273 | docstring / printed prose | no |
| `experiments/paper_outputs.py` | 703, 705 | docstring | no |
| **`diagnostics/solve_configs.py`** | **130, 167, 169, 173** | **executable** | **YES** |
| **`diagnostics/psilam_breakdown_sweep.py`** | **15, 66–88, 123–149** | **executable** | **YES** |
| **`diagnostics/psilam_moment_sweep.py`** | **7, 59, 65, 73–86** | **executable** | **YES** |

**Result for the audit's actual Step-3 target — the CB block, the clearing condition
and the residual equations (`code/tpi.py`, `equations_global.domestic_bond_clearing`,
`budget_residual_D/F`): ZERO live references.** The only `psi_spread` string in
`code/tpi.py` is a comment at line 313 explaining that the old
`loading ≈ 1 + psi_spread/EL_price` closed form no longer exists and is not being
replaced by another closed form.

The three `diagnostics/` scripts are dangling but are outside the CB block; they are
recorded in `VERDICT.md` under a separate heading. Their failure modes differ:

* `psilam_breakdown_sweep.py:70` and `psilam_moment_sweep.py:59` do a bare
  `float(ss["psi_spread_D"])` — these will **raise** on the current model.
* `solve_configs.py:130` is inside `try/except`, so it degrades to `MISSING (...)`;
  but line 169 then **writes** `ss0.toplevel["psi_spread_D"] = 0.0`, inventing a
  symbol no block reads, and line 173 reads it back successfully. That script does
  not crash — it silently produces a `psi_lambda_B = 0` arm that is now identical to
  its own baseline, because `psi_lambda_B` is already `0` at the live calibration.
  Silent, not loud: worse.

Corroborating automated check — the AST scanner that enforces the deletion:

```
/opt/anaconda3/envs/ssj/bin/python -m pytest code/test_nkpc_blocks.py -q
22 passed in 1.16s
```

`test_no_ad_hoc_sovereign_spread_wedge_anywhere` (`code/test_nkpc_blocks.py:261`)
scans `code/*.py` for `psi_spread`, `EL_price`, `divert_bond_foc`,
`divert_portfolio_adj`, `excess_return_` and passes. Note its scope is `code/` only —
it is structurally incapable of catching the three `diagnostics/` hits above.

Companion sweep for the other deleted names (`excess_return`, `bond_price_ss`,
`divert_bond_foc`, `divert_portfolio_adj`, `portfolio_adj_cost`): all remaining hits
are docstrings. `tau_mp_D/F` still exists as live code (`equations_D.py:630`,
`equations_F.py:441`) but is now a **lump-sum macro-prudential tax**, not a term in
the bond FOC, and is calibrated to zero (`T0_* = T1_* = 0.0`,
`code/calibration.py:446-447`).

---

## 2026-08-19 08:16 — Probe C: what the GK refactor actually changed in the CB layer

`git diff HEAD -- code/tpi.py` (48 lines) shows the refactor touched **only**:

* three import renames (`divert_bond_foc_D/F` → `gk_bond_foc_D/F`,
  `divert_portfolio_adj` → `gk_cross_border_foc`);
* `cb_pnl`'s expected-loss coefficient, `EL_price_D` → `EL_load_D`;
* the removal of the `loading ≈ 1 + psi_spread/EL_price` printed claim;
* a new closed-loop-pole guard capping the effectiveness grid at `0.75 × pole`.

**The four CB equations themselves — `domestic_bond_clearing_tpi`,
`budget_residual_D_tpi`, `budget_residual_F_tpi`, `external_account_D_tpi` — are
byte-identical to their pre-refactor versions.** They were never written against the
old FOC's price decomposition (they were written against the *payoff*, in
coupon/survival form), which is why the refactor did not break them. This is the
answer to Step 2's central worry, established by diff rather than by inspection.

---

## 2026-08-19 08:18 — Probe D: line-by-line price-marking trace (Step 2)

Every place a CB quantity meets a price, and which price it meets:

| Site | Expression | Price used | Verdict |
|---|---|---|---|
| `tpi.py:50` `domestic_bond_clearing_tpi` | `b_D_D = b_gov_D - size_F*b_D_F - cb_buy_D` | none (quantities) | n/a |
| `tpi.py:73` `cb_flow_D` coupon leg | `delta_b_D*(1 - def_rate_D*h*mult)*cb_buy_D(-1)` | none | matches `bond_return_D` |
| `tpi.py:74` `cb_flow_D` continuation leg | `q_b_D * surv_cont_D * (1-delta_b_D) * cb_buy_D(-1)` | **endogenous `q_b_D`** | OK |
| `tpi.py:75` `cb_flow_D` purchase leg | `- q_b_D * cb_buy_D` | **endogenous `q_b_D`** | OK |
| `tpi.py:108` `external_account_D_tpi` payments | `(1+rb_actual_D)*q_b_D(-1)*(size_F*b_D_F(-1) + kappa_cb_F*cb_buy_D(-1))` | **endogenous `q_b_D(-1)`** | OK |
| `tpi.py:110` `nfa_D` | `- q_b_D*(size_F*b_D_F + kappa_cb_F*cb_buy_D)` | **endogenous `q_b_D`** | OK |
| `tpi.py:290-295` `cb_pnl` | `q_b_D_ss * cb`, plus `mtm_pv` in `dq` | SS price + separate MTM | see below |
| `tpi.py:366` cost curve | `(cb_buy_D[:100] * q_b_D_ss).sum()` | SS price | see below |

Algebraic cross-check of the coupon/survival form against `bond_return_D`
(`code/equations_D.py:325-331`):

```
bond_return_D : payoff = delta_b*(1 - def*h*mult) + (1-delta_b)*q_b*(1 - zeta*def*h*mult)
cb_flow_D     : inflow = delta_b*(1 - def*h*mult)*cb(-1)
                       + q_b*(1 - zeta*def*h*mult)*(1-delta_b)*cb(-1)
```

Identical term for term, including the `zeta_writeoff_D` continuation write-down and
the `writeoff_enabled_D` realisation gate. **No stale price, no hardcoded price
decomposition, no surviving reference to the old FOC anywhere in the CB equations.**

The two SS-price sites are **not** stale-price bugs — they are the deliberate
linearisation convention. `cb_buy_ss = 0`, so `cb` is a pure deviation; `q_b_D_ss * cb`
is the correct first-order market value and `mtm_pv` restores the
(second-order) `dq × cb` revaluation as a separate line, exactly as
`tpi.py:252-255` documents for `el_pv` and `prem_pv`. Flagged in `VERDICT.md` as a
*documented convention*, not a defect.

One genuinely stale string was found in the CB plotting layer:
`code/tpi_plots.py:243` hardcodes the axis subtitle
`[δ_b = 0.10 → insensitive to q_b_D]`. Live values are `delta_b_D = 0.0777` and
`delta_b_F = 0.0568`. Cosmetic; it is a figure caption, no number is computed from it.

---

## 2026-08-19 08:20 — Probe E: live pipeline solve

`diagnostics/cb_audit/probe_pipeline.py` — runs `main.py` stages 1–5 plus `run_tpi()`,
skipping figure generation, and dumps Steps 0/4/5 quantities to
`probe_pipeline.json` / `.npz`. Started 08:20:41.

`EBA_CALIBRATION=True`, `BANK_SCOPE="broad"` read at 08:20:41.

## 2026-08-19 08:19 — Probe F: Prony estimator

`docs/STATE.md:2317` and `docs/PROGRESS.md:1499` refer to a "validated order-selected
Prony" eigenvalue estimator. **No such file exists in the working tree** (`grep -rn
-i prony` matches only those two prose lines) — it went with the retired
`audit_artifacts/` harness. Reimplemented self-contained as
`diagnostics/cb_audit/prony.py`, with order selection (the property
`docs/STATE.md:2270` records as necessary: fixed-order fits overfit and manufacture
spurious near-unit-circle roots) and a synthetic self-test:

```
  single real 0.90         -> |lam|=0.900000000 (truth 0.90, err 1.11e-16, order 1, R2 1.000000000)
  two real 0.95/0.60       -> |lam|=0.950000000 (truth 0.95, err 1.22e-15, order 2, R2 1.000000000)
  damped osc r=0.97        -> |lam|=0.970000000 (truth 0.97, err 5.55e-16, order 2, R2 1.000000000)
  selftest PASS
```


---

## 2026-08-19 08:27 — Probe E results: Step 0 (which base) and Step 5 (SS neutrality)

`probe_pipeline.py` completed 08:26:54 (first attempt 08:23:27 died in the dump
stage on a `KeyError: 'shock_def_D'` — `rem_cb_D`, `rem_cb_F` and `cb_flow_D` are
outputs with a `cb_buy_D` column and no `shock_def_D` column, because nothing but the
CB's own purchases moves them. The guard was widened and the probe re-run; run 1's
stdout is kept as `probe_pipeline_stdout_run1.txt`. The two runs agree on every
number.)

**Step 0 — the base is ADMISSIBLE.** `EBA_CALIBRATION = True`, `BANK_SCOPE = "broad"`.

| object | D | F |
|---|---|---|
| `lambda_gk` (IC multiplier) | **+2.108746** | **+0.641387** |
| `Omega` (banker marginal value of net worth) | **+10.346180** | **+4.037867** |
| `nu_K` (marginal value of capital per unit net worth) | **+0.103205** | **+0.040278** |
| `phi_own` = `q_b·b_own/n_inter` | **0.452489** | 0.294078 |
| `phi_cross` | 0.003406 | 0.010710 |
| GK feasibility margin `f·theta − [(1−Δ_own)φ_own + (1−Δ_cross)φ_cross]` | **+0.296567** | **+0.591774** |

All strictly positive; `assert_gk_well_posed` passes. `theta_D = 5.5107`, `f_D = 0.12`,
so at `Delta = 0.20` the well-posedness ceiling on `phi_own_D` is
`(f·theta − 0.8·phi_cross)/0.8 = 0.8232` and the live value 0.4525 sits at 55% of it.

The GK portfolio FOC table on the final SS:

```
  bank  leg       nu_i/nu_K   Delta_eff     residual   status
  D     own        0.200000    0.200000   -1.624e-14       OK
  D     cross      0.200000    0.200000    2.058e-13       OK
  F     own        0.200000    0.200000   -4.746e-15       OK
  F     cross      0.200000    0.200000   -2.268e-13       OK
  cross-border SS wedges:  F-in-D = +8.2e-11 bp/yr   D-in-F = -9.1e-11 bp/yr
  psi_lambda_B = 0.0/0.0   zeta_writeoff = 1.0/1.0   writeoff_enabled = 0.0/0.0
```

`q_b_D = 0.974906`, `q_b_F = 0.965974`, both yields 80.00 bp annualised, SS spread
−0.0000 bp. `EL_load_D = 0.7014`. All reproduce the CLAUDE.md calibration row exactly.
Baseline default-shock impact also reproduces exactly: peak spread 205.87 bp,
`n_inter_D[0] = −11.4073%`, `Y_D[0] = −1.9742%`.

**The audit brief's Step-0 premise does not describe this tree.** It anticipates
`phi_bD ≈ 2.39` against a ceiling of 0.336. `2.39` is the `BANK_SCOPE="ct1"`
(stress-test-sample) number that GK-2 retired on 2026-07-31; the live broad-scope
value is 0.4525 and it is inside the ceiling. There is no inadmissible-base caveat to
attach to this audit.

**Step 5 — SS neutrality is EXACT, not approximate.** Each TPI block was evaluated at
the solved steady state with `cb_buy_D = 0`, alongside its non-TPI counterpart:

| block output | TPI version | non-TPI version | difference |
|---|---|---|---|
| `b_gov_res_D` | −4.163336e−17 | −4.163336e−17 | **0.000e+00** |
| `nfa_D` | −1.338351e−01 | −1.338351e−01 | **0.000e+00** |
| `ca_res_D` | +1.665335e−16 | +1.665335e−16 | **0.000e+00** |
| `b_D_D` | +9.922837e−01 | +9.922837e−01 | **0.000e+00** |
| `b_F_F` | +4.953501e−01 | +4.953501e−01 | **0.000e+00** |

and the conduit's own objects are identically zero:
`cb_flow_D = 0.000e+00`, `rem_cb_D = 0.000e+00`, `rem_cb_F = 0.000e+00`,
`b_gov_res_F = −2.776e−17`. Bit-identical, so the refactor cannot have moved the
steady state through the CB block and every prior SS-invariance argument survives.

Corroborated dynamically by `run_tpi`'s own gate:
`Sanity check G_tpi[cb=0] vs baseline G: max |err| = 0.00e+00`.

---

## 2026-08-19 08:29 — Probe G: Step 4, feedback sign

`probe_stability.py` (output `probe_stability.md`).

```
A_cb[0,0] = -4.397083e-03      NEGATIVE: a unit purchase compresses the impact spread
sum_h A_cb[h,0] = -1.597028e-02   cumulative response to a one-period purchase
```

Closed-loop peak spread, monotone and declining:

| γ | peak spread (bp ann) | compression | `b_gov_D[499]` | `n_inter_D[0]` | `Y_D[0]` |
|---|---|---|---|---|---|
| 0 | 205.87 | — | +3.02e−05 | −2.4388e−01 | −1.9742e−02 |
| 2 | 193.25 | 6.1% | +3.69e−04 | −2.1289e−01 | −1.5702e−02 |
| 5 | 176.75 | 14.1% | −5.62e−05 | −1.7260e−01 | −1.0428e−02 |
| 10 | 154.36 | 25.0% | −1.51e−04 | −1.1837e−01 | −3.2741e−03 |

Prony dominant moduli (estimator validated above), all < 1 and barely moving with γ:
spread 0.9400 → 0.9470, `b_gov_D` 0.8283 → 0.9232, `n_inter_D` 0.9347 → 0.9460.
Walras residuals stay clean along the whole loop: `max|ca_res_D|` *falls* 7.2e−08 →
1.5e−08 as γ rises, `max|goods_mkt_D|` ~6e−17, `max|goods_mkt_F|` ~2e−10.

The mechanism, read straight off the Jacobian columns at impact:

```
 d b_D_D   /d cb_buy [0,0] = -0.126779     bonds off D banks' books
 d theta_D /d cb_buy [0,0] = -1.177721     required leverage falls (IC slackens)
 d n_inter_D/d cb_buy[0,0] = +0.525814     bank net worth recovers
 d K_D     /d cb_buy [0,0] = +0.017015     capital crowded back in
 d q_b_D   /d cb_buy [0,0] = +0.037316     bond price up -> spread compresses
```

**Step 3's second half confirmed structurally, not just by grep.** `cb_buy_D` reaches
the price through exactly two doors — `intermediation_IC_D`'s `phi_bD_D` and
`k_balance_sheet_D` — and through nothing else. There is no term of the form
`spread += parameter × def_rate` anywhere in the path.

Note further that with `psi_lambda_B = 0` the own-leg FOC composed with
`intermediation_P1_D` reduces *exactly* to

```
rb_exp_D(+1) - rdep_D = 0.20 * (rk_D(+1) - rdep_D)
```

— `SDF_banker` and `Omega_p1` cancel in the ratio `nu_bD/nu_K`. Verified at the SS:
`rb_exp_D = 0.002000`, `rk_D = 0.010000`, `rdep_D = 0`, and `0.20 × 0.01 = 0.002`
exactly. So the CB compresses the sovereign spread **only** by lowering `rk_D(+1)`,
i.e. by crowding capital back in. TPI's spread effect and its investment effect are
the same effect in this model; they cannot be reported as separate channels.

---

## 2026-08-19 08:33 — Probe H: the "closed-loop pole" is a truncation artefact

`run_tpi` prints `closed-loop pole at gamma = 26.50`. Cross-checking that against the
spectrum of `A_cb` did not agree, so the discrepancy was chased down.

Eigenvalues of `A_cb` (T = 500):

```
   +0.452155        |lam|=0.452155   pole gamma = 2.2116
   +0.036524        |lam|=0.036524   pole gamma = 27.3793
   -0.023262±0.013576j                 (complex, no pole)
```

Direct scan of `I - gamma*A_cb`:

```
gamma   smin           cond        sign(det)   peak bp
2.1000  6.1119e-03  1.349e+03     +1          192.65
2.2000  6.0812e-04  1.419e+04     +1          192.06
2.2116  1.6529e-06  5.248e+06     +1          189.76
2.2500  1.9616e-03  4.498e+03     -1          191.77
```

**There is a genuine singularity of the T=500 system at γ = 2.2116, inside the
intended γ range, and the guard in `code/tpi.py:332-353` misses it.** The guard scans
`np.linspace(0.25, 60.0, 240)` — step 0.25 — so it evaluates 2.00 (cond 6.5e2) and
2.25 (cond 4.5e3), both under the 1e4 threshold, and steps straight over. What it
then reports at 26.50 is the *second* eigenvalue's shoulder. This is precisely the
failure its own comment warns about ("a coarse scan steps over a pole this narrow and
sees nothing").

**But the pole is not economic.** Modal decomposition of the resonant eigenvector `v`:

```
|v| peaks at index 499 (the LAST period of the truncation), max|v| = 0.2779
mass in the first 100 quarters = 0.0000
mass in t = 400..499        = 0.9922
projection of the forcing on the resonant left eigenvector = 6.0e-04 (relative)
```

and the terminal columns of `A_cb` are pathological while every interior column is
not:

```
||A[:,  0]|| = 5.88e-03      A[  0,  0] = -4.397e-03
||A[:,200]|| = 6.51e-03      A[200,200] = -2.176e-03
||A[:,400]|| = 6.43e-03      A[400,400] = -2.184e-03
||A[:,495]|| = 2.14e-01      A[495,495] = -5.204e-02
||A[:,499]|| = 3.86e+00      A[499,499] = +1.080e+00   <- only positive diagonal, only |.|>1
```

The dominant eigenvalue only exists once the last few columns are included
(leading-block spectra: T≤450 → |λ|max ≈ 0.027 complex, no positive real root at all;
T=495 → −0.416; T=500 → +0.452).

Decisive test — drop the last k rows and columns and re-solve the closed loop:

| k dropped | smallest positive-real pole γ | max cond over γ∈[0.1,30] | peak bp @ γ=2 / 5 / 10 / 20 |
|---|---|---|---|
| 0 | **2.2116** | 5.45e+05 | 193.25 / 176.75 / 154.36 / 122.59 |
| 5 | 36.81 | **8.51e+01** | 193.25 / 176.75 / 154.36 / 122.56 |
| 20 | 38.18 | 1.87e+01 | 193.25 / 176.75 / 154.36 / 122.56 |
| 100 | 40.96 | 7.76e+00 | 193.25 / 176.75 / 154.36 / 122.55 |

Dropping five columns collapses the worst conditioning by four orders of magnitude,
removes every pole below γ = 36, and **changes the reported peak spread by nothing at
γ = 2, 5, 10 and by 0.03 bp at γ = 20**. The residual "pole" that remains keeps
migrating with k, so it is boundary-driven too.

Conclusions: (a) no reported TPI number is affected — every statistic is computed on
`[:100]`, where the artefact has zero mass; (b) the γ = 19.88 cap the effectiveness
curve now carries is imposed for a spurious reason; (c) the claim in `CLAUDE.md`
and `code/tpi.py` that there is a closed-loop pole at γ ≈ 27.3 on this calibration is
wrong — that is the T=500 terminal condition, not the model.

---

## 2026-08-19 08:36 — Probe I: Step 1, where the money actually goes

Impact remittances (D-goods; F converted at `p` and per F capita in the block, scaled
back to an aggregate here with `size_F = 11.6967`):

| γ | `cb_flow_D[0]` | `rem_cb_D[0]` (D aggregate) | `rem_cb_F[0]` × `size_F` (F aggregate) | ratio F/D |
|---|---|---|---|---|
| 2 | −9.4199e−03 | −6.6881e−04 | −8.7338e−03 | 13.059 |
| 5 | −2.1539e−02 | −1.5293e−03 | −1.9971e−02 | 13.059 |
| 10 | −3.7622e−02 | −2.6711e−03 | −3.4882e−02 | 13.059 |

Capital key implies `0.929/0.071 = 13.085`; the 0.2% gap is the endogenous
terms-of-trade conversion (`p = 1.00197`). The split is correct and it is not a
self-financing loop.

Both remittances are **negative at impact** — a capital call. Each treasury funds its
share through its own fiscal rule at its own sovereign terms.

Downstream fiscal incidence:

| γ | `TAX_D[0]` | `TAX_F[0]` | Σ₁₀₀ `TAX_D` | Σ₁₀₀ `TAX_F` | `C_D[0]` | `C_F[0]` | ΔW_D | ΔW_F |
|---|---|---|---|---|---|---|---|---|
| 0 | −1.2449e−02 | +6.2537e−04 | +9.4021e−02 | −1.4470e−03 | −1.6959e−02 | +1.5843e−03 | — | — |
| 2 | −1.0576e−02 | +4.6206e−04 | +9.0618e−02 | −1.4255e−03 | −1.2681e−02 | +1.2794e−03 | +0.1956 | −0.0195 |
| 5 | −8.1361e−03 | +2.4896e−04 | +8.6147e−02 | −1.3804e−03 | −7.0811e−03 | +8.8080e−04 | +0.5155 | −0.0482 |
| 10 | −4.8410e−03 | −3.9851e−05 | +8.0123e−02 | −1.2724e−03 | +5.4953e−04 | +3.3869e−04 | +1.1290 | −0.0944 |

German cumulative taxes rise with γ (Σ₁₀₀ `TAX_F` moves −1.447e−03 → −1.272e−03),
German consumption gain shrinks, German welfare falls monotonically. **The creditor
side is real, signed correctly, and quantitatively representable.**

Who actually gives up the bonds — differentiating the clearing identity
`b_D_D = b_gov_D − size_F·b_D_F − cb_buy_D` with respect to `cb_buy_D`:

| horizon | `d b_D_D` | `d b_gov_D` | `size_F · d b_D_F` |
|---|---|---|---|
| 0 | −0.1268 | +0.0325 | **−0.8408** |
| 1 | −0.0671 | −0.0539 | **−0.9867** |
| 4 | −0.0589 | −0.0488 | **−0.9899** |
| 12 | −0.0253 | −0.0202 | **−0.9949** |

84% of a CB purchase at impact, and ~99% from t=1 on, is absorbed by **German** banks
shedding their Greek book, not by Greek banks. Only 13% at impact and 2–6% thereafter
comes off Greek balance sheets. Recorded as F-3 in `VERDICT.md`; `psi_bD_F = 0.5` is
the elasticity that governs it.

---

## 2026-08-19 10:09 — Probe J: the 2×2 sovereign-holdings matrix (decisive diagnostic)

`probe_portfolio.py`, full output `portfolio_matrix.md` / `.json` / `.npz`. Required
its own solve: `probe_pipeline.py` dumped `b_D_D` and `b_D_F` only, and the four legs
carry different per-capita conventions, so partial data cannot be rescaled after the
fact.

**Units, which are the point of the diagnostic.** `b_D_D` and `b_F_D` are per D capita
(D size = 1); `b_D_F`, `b_F_F` and `b_gov_F` are PER F CAPITA and take `size_F =
11.696651`; `cb_buy_D` is already a D aggregate and takes no weight. `q_b_D` and
`q_b_F` are both D-good prices, so market values need no `p` conversion. Market-value
deviations are `q·db + b·dq`, dropping the second-order cross term, matching
`cb_pnl`'s convention.

**Both clearing identities hold to ≤ 1.2e−15** at the SS and at t = 0, 4, 20 for every
γ ∈ {0, 2, 5, 10} — 26 checks, all machine-zero. The matrix is internally consistent.

Steady state (aggregate market value, D goods):

```
holder      D paper    F paper      total    share of D issue
D banks    0.967383   0.007282   0.974666      87.27%
F banks    0.141118   5.596796   5.737913      12.73%
CB         0.000000   0.000000   0.000000       0.00%
issued     1.108501   5.604078   6.712579
```

γ = 0 (crisis, no TPI) → γ = 10, at impact:

```
             D paper              share of D issue
           g=0       g=10        g=0      g=10
D banks  0.930292  0.926341     87.76%   86.86%
F banks  0.129756  0.102520     12.24%    9.61%
CB       0.000000  0.037622      0.00%    3.53%
issued   1.060048  1.066483
```

Sourcing of the CB book, Δ(γ=10 − γ=0):

```
F banks / D paper   -0.027236   72.4% of the book
D govt new issue    +0.006435   17.1%
D banks / D paper   -0.003951   10.5%
CB                  +0.037622  100.0%
identity: -0.003951 - 0.027236 + 0.037622 = +0.006435 = D govt issue change
```

Concentration `phi = q_b·b/n_inter`, what `intermediation_IC_D/F` reads:

```
                          SS       g=0       g=2       g=5      g=10
phi_bD_D (D bank, own)  0.452489  0.486757  0.479736  0.470614  0.458345
phi_bF_D (D bank, F)    0.003406  0.003936  0.004167  0.004463  0.004853
phi_bD_F (F bank, D)    0.007415  0.006774  0.006428  0.005983  0.005397
phi_bF_F (F bank, own)  0.294078  0.293864  0.294447  0.295200  0.296203
```

Three results:

1. The shock raises `phi_bD_D` by +0.034268; γ = 10 removes 0.028412 of that —
   **82.9% of the crisis concentration spike undone**.
2. It is a **denominator** effect. Numerator `q_b_D·b_D_D` moves −0.4% between γ = 0
   and γ = 10 (0.930292 → 0.926341); `n_inter_D` recovers +6.6%. The CB barely takes
   bonds off Greek banks, and does not need to.
3. `phi_bD_F` falls monotonically in γ, to 27% below its SS level at γ = 10, while
   `phi_bF_D` rises 42% above SS. **TPI accelerates cross-border retrenchment in both
   directions**; it does not reverse it.

Cross-check against Probe I: in market value F banks supply 72.4% of the book, in pure
quantity (differentiating the clearing identity, price effect stripped out) 84% at
impact and ~99% from t = 1. The gap is the price support inflating every holder's
mark. Both figures are correct; the quantity figure is the one to quote for "who sold".

---

## 2026-08-19 10:15 — Canonicalisation

Audit promoted to `docs/cb_mechanism.md` (**canonical**), with an index row added to
CLAUDE.md's docs table. No model source modified at any point in this audit; the
`code/` entries in `git status` are the pre-existing uncommitted GK refactor.
