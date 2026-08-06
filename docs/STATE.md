# Project State

**Branch:** `eba-recalibration` | **Date:** 2026-07-31 | **Status:** **EBA calibration rebuilt, identified, and LIVE** (`EBA_CALIBRATION=True`, `BANK_SCOPE="broad"`). Y-1 and RK-1 resolved; spread on target at 150.4bp; TPI loading declining.

## Nominal rigidities (`add-nkpc`): Task 11 — nominal deposit contracts split into `i_dep`/`rdep`/`rdep_expost`

`rdep_D/F` is split into three names to avoid a wide, risky rename. `deposit_rates_D/F`
(new `@simple` block) takes the nominal unknown `i_dep_D/F` and `pi_D/F` and returns two
outputs: `rdep_D/F` (unchanged name, unchanged meaning — the ex-ante real rate for
t → t+1, locked at t) and `rdep_expost_D/F` (new — the realised real rate at t on
deposits placed at t−1, `= (1+i_dep(-1))/(1+pi) - 1`, carrying the inflation surprise).
`deposit_return_D/F` now takes `i_dep`, `P_CES`, `pi` instead of `rdep`, keeping the T-2
timing (rate locked at t−1) and adding the period-t CPI deflator via `pi`.

Because `rdep_D/F` keeps its name and its ex-ante meaning, `intermediation_P1_D/F`,
`divert_bond_foc_D/F` and `divert_portfolio_adj` (global) are **untouched** — verified by
introspecting `.inputs` after the change (`rdep_D` present, `rdep_expost_D` absent).

At SS `pi = 0`, so `rdep_D = rdep_expost_D = i_dep_D` and `Rgross_D = 1 + i_dep_D`
exactly — this is what keeps the steady state bit-identical once `i_dep` is wired in as
the solver unknown (Task 13). `deposit_rates_D/F` and the new `deposit_return_D/F` are
not yet wired into `full_model.py`'s block list or the calibration/steady-state — that is
Task 13.

`code/test_nkpc_blocks.py` gained 4 tests (11 → 15 passed): SS collapse, deflation raises
the realised real rate (sign check for the Fisher channel), `Rgross` unchanged at
`pi=0`, and D/F symmetry.

## Nominal rigidities (`add-nkpc`): Task 12 — `bank_return_D/F`/`capital_fund_D/F` now pay the ex-post real rate

The Fisher channel: `bank_return_D/F` and `capital_fund_D/F` pay for deposits taken on at
t−1, so they must use the *realised* real rate on that funding, not the ex-ante one.
Signature changed `rdep_D/F` → `rdep_expost_D/F`; all four (three, for `capital_fund`)
occurrences of `rdep_D(-1)`/`rdep_F(-1)` in the bodies became bare `rdep_expost_D/F`
(`rdep_expost` already carries its own `(-1)` internally — writing `rdep_expost_D(-1)`
would double-lag it). `intermediation_P1_D/F` and `divert_bond_foc_D/F` are untouched —
they are ex-ante (t → t+1) and still correctly read `rdep_D/F`. F's W-2 `p(-1)/p`
terms-of-trade conversion in `bank_return_F` is unrelated and was left alone.

Sign logic: D banks are net nominal debtors on deposits. When D deflates,
`rdep_expost_D` rises (the real value of what the bank owes goes up), so `bank_return_D`
falls further — this is what Task 13's gate checks (`n_inter_D[0]` more negative, not
less). Getting the substitution backwards would silently invert the whole channel.

`code/test_nkpc_blocks.py` gained 2 tests (15 → 17 passed): a signature check that
`bank_return_D`/`capital_fund_D` take `rdep_expost_D` and not `rdep_D`, and a check that
`intermediation_P1_D`/`divert_bond_foc_D` still take `rdep_D` and not `rdep_expost_D`.

## Nominal rigidities (`add-nkpc`): Task 13 — nominal deposits WIRED IN; Fisher gate PASSED

The model solves end-to-end again. Three edits, no equation changed:

- `code/calibration.py` — `'rdep_D'/'rdep_F': 0.000` renamed to `'i_dep_D'/'i_dep_F': 0.000`.
  `rdep_D/F` is no longer a calibration parameter; it is an output of `deposit_rates_D/F`.
  Every `ss['rdep_D']` read elsewhere (`steady_state.py` excess-return anchors and
  `Rgross` seeding, `depreciation_calibration.py` printout, `irf_plots.py`,
  `diagnostics/solve_configs.py`) still works unchanged — the object exists, it is just
  solved rather than supplied.
- `code/steady_state.py` and `code/full_model.py` — `deposit_rates_D/F` imported and
  placed immediately before `deposit_return_D/F` in both block lists.
- `code/full_model.py` — `unknowns_tp`: `rdep_D → i_dep_D`, `rdep_F → i_dep_F`.
  `targets_tp` **unchanged**; the system is still 27×27.

**SS gate: PASSED, bit-identical.** Every numeric steady-state line diffs clean against
the pre-change baseline: `goods_mkt_D = -4.2493506589857954e-07`,
`goods_mkt_F = -4.1914559989475464e-07`, `ca_res_D = 6.852157730108388e-17`,
`IC_D: θ − θ_tgt = 1.776357e-15`, all residuals `< 1e-8 ✓`, `ρ_b = 0.8451`,
`b_gov_D[499] = 1.6e-05`. And the three rates collapse exactly at `pi = 0`:
`i_dep_D = rdep_D = rdep_expost_D = 0.000000000000`.

**Fisher gate: PASSED.** Impact response to the D default shock, sticky prices in both
columns; the only difference is real vs nominal deposit contracts.

| | sticky, REAL deposits (Task 10) | + NOMINAL deposits (Task 13) | change |
|---|---|---|---|
| `Y_D[0]`       | −0.4923% | **−0.5449%** | −0.0525 |
| `C_D[0]`       | −0.4904% | **−0.5499%** | −0.0595 |
| `I_D[0]`       | −0.9907% | **−1.0849%** | −0.0942 |
| `n_inter_D[0]` | −4.0140% | **−4.6155%** | −0.6016 |

The sign is right and the ordering is right: the amplification lands hardest on bank net
worth (−0.60pp, ~15% deeper), and only reaches output through the intermediary. That is
the Fisher channel doing exactly what it should — D banks are net nominal debtors on
deposits, the crisis deflates D, `rdep_expost_D` rises, and the realised real funding
cost eats further into `n_inter_D`. Had the Task 12 ex-post/ex-ante substitution been
backwards, `n_inter_D[0]` would have gone *less* negative; it does not.

First 8 quarters, % of SS:

```
Y_D real dep   : -0.4923 +0.0083 +0.1006 +0.0881 +0.0580 +0.0329 +0.0161 +0.0063
Y_D nominal dep: -0.5449 -0.0011 +0.1021 +0.0912 +0.0604 +0.0342 +0.0165 +0.0061
C_D real dep   : -0.4904 +0.1141 +0.1553 +0.0790 +0.0046 -0.0440 -0.0685 -0.0761
C_D nominal dep: -0.5499 +0.1144 +0.1676 +0.0903 +0.0124 -0.0395 -0.0665 -0.0758
```

The **one-quarter-spike caveat from Task 10 survives unchanged** and must not be dropped
from the write-up: nominal deposits deepen the impact quarter but do not fix the sign
reversal from quarter 2 on — `Y_D` is already back to ≈0 at q1 and positive by q2 in both
columns. Nominal contracts amplify the impact; they do not lengthen the recession.
`C_D[1]` is essentially unmoved (+0.1141 → +0.1144), so the whole Fisher effect is an
impact-quarter effect.

`code/test_nkpc_blocks.py`: **17 passed, 0 failed** (unchanged — Task 13 added no
equations, so no new tests).

Caveat carried forward: the spread moment has drifted with the nominal side and has not
been re-tuned. That is Task 14.

**Next: Task 14** — re-tune `psi_lambda_B` to the 150bp spread moment on the full
sticky-price / nominal-deposit model.

## Nominal rigidities (`add-nkpc`): Task 10 — price-stickiness-only result at the calibrated slope

Deposits are still **real** here; this isolates what price stickiness alone does.
`kappa_p = 0.0871` (Calvo theta=0.75 at beta=0.985). Steady state byte-identical
to pre-change (`goods_mkt_D = -4.2493506589857954e-07`, `IC_D = 1.776357e-15`,
`rho_b = 0.8451` all unchanged); both impact signs negative; all four TPI gammas
inside 1e-7.

### Impact responses to the 1pp default shock (% of own SS level)

| | flex (pre-change) | sticky `kappa_p=0.0871` | change |
|---|---|---|---|
| `Y_D[0]` | −0.0149 | **−0.4923** | 33x larger |
| `C_D[0]` | **+0.2164** | **−0.4904** | **SIGN FLIP** |
| `I_D[0]` | −0.7718 | −0.9907 | 1.28x |
| `n_inter_D[0]` | −3.3804 | −4.0140 | 1.19x |

The flex column reproduces this file's own E1 passive row to every printed digit,
which is the check that the SS divisors are right.

**`C_D[0]` DOES change sign.** The spec warned it might not, and that Bi-Foerster-Traum
get consumption rising on impact even with a Taylor rule and nominal debt. Here
price stickiness alone flips it.

**But read the path before claiming the consumption-bust problem is solved.**

```
Y_D sticky:  -0.4923  +0.0083  +0.1006  +0.0881  +0.0580  +0.0329  +0.0161  +0.0063
C_D sticky:  -0.4904  +0.1141  +0.1553  +0.0790  +0.0046  -0.0440  -0.0685  -0.0761
C_D flex:    +0.2164  +0.0702  -0.0161  -0.0624  -0.0829  -0.0873  -0.0822  -0.0721
```

This is a **one-quarter spike, not a downturn**: output and consumption are both
positive from quarter 1. Bi-Foerster-Traum's output falls 0.6% and stays negative
for ~20 quarters. Flex consumption is in fact *more* persistently negative from
quarter 2 on. So the honest claim is: sticky prices fix the **impact quarter**
and arguably worsen the medium-run path. Do not write "resolves the investment-bust
counterfactual" without this qualification.

**The extra output decline is not an investment story.** `I_D` moved 1.28x while
`Y_D` moved 33x — it is the markup wedge shifting labour demand directly, which is
the mechanism the design intended.

### `kappa_p` sweep — monotone and stable throughout

| `kappa_p` | `Y_D[0]` | `C_D[0]` | `I_D[0]` | `n_inter_D[0]` | `b_gov_D[499]` |
|---|---|---|---|---|---|
| 0.03 (stickier) | −0.7280 | −0.8401 | −1.0987 | −4.3244 | 1.50e−05 |
| **0.0871** (calibrated) | **−0.4923** | **−0.4904** | **−0.9907** | **−4.0140** | 1.4e−05 |
| 0.2 (more flexible) | −0.3105 | −0.2211 | −0.9077 | −3.7736 | 1.42e−05 |
| flex limit | −0.0149 | +0.2164 | −0.7718 | −3.3804 | 1.4e−05 |

Every point stable (`b_gov_D[499]` ~1.5e−05) and monotone in `kappa_p`. **The
consumption sign flip is not knife-edge** — `C_D[0]` is negative across the whole
sticky range and only positive in the flexible limit.

### Benchmark context

`Y_D[0] = -0.4923%` against Bi-Foerster-Traum (FRBSF WP 2025-10) at −0.6%. The
pre-change model was two orders of magnitude below that and the spec flagged the
gap as the main risk to the exercise. Price stickiness closes it **without** their
Sims-Wu loan-in-advance constraint.

## Nominal rigidities (`add-nkpc`, in progress): Task 9 — 27×27 sticky-price system LIVE, flex-price equivalence gate PASSED

The dynamic model is now sticky-price. `full_model.build_block_list()` carries
`firm_profit_D/F`, `price_nkpc_D/F`, `terms_of_trade`, `union_inflation`, and
`build_and_solve` solves **27×27** (was 23×23): `mc_D, pi_D, mc_F, pi_F` added to
`unknowns_tp`; `nkpc_p_res_D/F, tot_res, union_pi_res` added to `targets_tp`.
Nothing existing was renamed or dropped; `labor_mkt_res_D/F` is unchanged because
wages stay flexible.

**Equivalence gate result (a real number, quote it).** As `kappa_p -> inf` the
NKPC forces `mu_p*mc -> 1` and labour demand collapses to competitive, so the
27×27 system must reproduce the pre-change 23×23 IRFs. It does, with **textbook
O(1/kappa_p) convergence**:

| `kappa_p` | worst relative IRF deviation vs pre-change baseline |
|---|---|
| 1e4 | `2.925e-03` |
| **1e5** | **`2.925e-04`** ← gate threshold 1e-3, **PASSED** |
| 1e6 | `2.925e-05` |

Every one of the 30 IRF series shrinks by a ratio of **exactly 10.00 per decade
of `kappa_p`**, and the 4 SS levels are bit-identical (`0.000e+00`) at every
`kappa_p`. This is a *limit*, not an identity — a finite `kappa_p` cannot match
to machine precision, and the clean 1/kappa_p rate (rather than a stuck floor) is
the evidence that the blocks are wired correctly. The binding series are
`irfs_def_D__w_D` (2.925e-04 at 1e5) and `irfs_def_D__N_D` (1.018e-04) — the two
objects the markup wedge acts on directly. Comparison harness:
`code/dump_irfs.py`.

**Open item carried into Task 10 — `full_model.solve_jacobian_padded()`.** SSJ
1.0.0's `Block.solve_jacobian` cannot solve this system as-is.
`CombinedBlock._jacobian` ends with
`total_Js[original_outputs & total_Js.outputs, :]` and only visits a block whose
inputs intersect the shock list, so a target that is a pure function of the
solver's *own unknowns* never enters H_Z. All four new targets are exactly that
(no `Z_*`/`shock_def_*` symbol appears anywhere in them), so SSJ returned a
23-row H_Z against a 27×27 H_U and numpy raised
`size 11500 is different from 13500`. The new helper restores those rows as
zeros — **exact, not an approximation**, since dH/dZ at fixed unknowns is
identically zero when the shock does not appear in the equation — and otherwise
reproduces `Block.solve_jacobian` line-for-line. `code/tpi.py` routes through it
too. It prints the padded row names on every solve so the padding can never go
silent.

**Task 9b (2026-08-06) — every Jacobian call site now routes through it.**
The seven remaining raw `Block.solve_jacobian` calls
(`diagnostics/regimes/regime_model.py:177`, `experiments/e4_distribution.py:255`,
`diagnostics/solve_configs.py:176`, `diagnostics/psilam_moment_sweep.py:76`,
`diagnostics/psilam_breakdown_sweep.py:83`,
`diagnostics/substitution_v2/solve_v2.py:106`,
`diagnostics/substitution_v2/exp_psilam0.py:64`) were converted to
`full_model.solve_jacobian_padded`, imported locally per-function (matching
each file's existing local-import style). `grep -rn "\.solve_jacobian("
--include="*.py" code experiments diagnostics | grep -v solve_jacobian_padded`
now returns nothing — no call site in the repo can hit the 23-vs-27 H_Z/H_U
shape mismatch (`size 11500 is different from 13500`) any more.
`code/test_nkpc_blocks.py` unaffected, still 11/11. The two call sites on the
critical path for `experiments/run_all.py` — `regime_model.py` (the E1-E4
Jacobian cache) and `e4_distribution.py` (quintile incidence) — are now safe to
re-run against the sticky-price system; the regime cache itself has not yet
been rebuilt (`diagnostics/regimes/regime_model.py --force`), that's for
whoever next needs fresh E1-E4 numbers.

**Full-pipeline regression at the flex limit** (`code/main.py`, `kappa_p=1e5`)
reproduces the pre-change baseline on every monitored number:
`goods_mkt_D=-4.2493506589857954e-07`, `goods_mkt_F=-4.1914559989475464e-07`,
`ca_res_D=6.852157730108388e-17`, `All residuals < 1e-8 ✓`,
`b_gov_D[499]=0.000014`, `ρ_b=0.8451`, `n_inter_D[0]=-3.3804%` of SS,
`Y_D[0]=-0.0149%` of SS (both negative), TPI `max|ca_res_D| ≤ 7.55e-08`,
`max|goods_mkt_F| ≤ 2.44e-09`, both TPI sanity checks exactly `0.00e+00`.
`code/test_nkpc_blocks.py`: 11 passed, 0 failed.

**Next: Task 10** — dial `kappa_p` to the calibrated `0.0871` and record the
price-stickiness result. That is the first run whose IRFs are *supposed* to
differ from the baseline.

## Nominal rigidities (`add-nkpc`): Task 8 — blocks seeded into the steady state

`code/steady_state.py`'s `create_model([...])` list (inside `solve_steady_state`)
gains `firm_profit_D/F`, `price_nkpc_D/F`, `terms_of_trade`, `union_inflation` —
the same blocks Tasks 2, 3 and 5 added to the dynamic model, now also carried
through the SS solve (they were previously absent from `steady_state.py`'s
import lists and block list, so `mc`, `pi`, `profit`, `tot_res`, `union_pi_res`
were computed dynamically but not at SS). `labor_demand_D/F` deliberately
stays out of the SS block list (SS still uses `labor_ss_D/F`), so Task 4's
markup wedge in labour demand cannot touch the SS solve. Verified: `profit_D`,
`profit_F`, `nkpc_p_res_D`, `nkpc_p_res_F`, `tot_res`, `union_pi_res` are all
exactly `0.000000e+00` at the solved SS, and the full initial-SS→IC-delta→
depreciation-calibration chain reproduces the pre-change baseline byte-for-byte
(portfolio shares, betas, p, IC residuals, delta, rk, final betas, rb/rdep/q_b,
goods_mkt_D/F, ca_res_D all identical). `code/test_nkpc_blocks.py`: still 11
passed, 0 failed. Next: Task 9, wire the 27×27 dynamic system and run the
`kappa_p = 1e4` (near-rigid-price) equivalence gate against the flex-price
Jacobian.

## Model plumbing (2026-08-05): single block-list definition

The model's block list now has **one definition**, `full_model.build_block_list()`.
`code/tpi.py` and `diagnostics/regimes/regime_model.py` previously each hardcoded
their own copy of the `sj.create_model([...])` list; both now call the shared
factory, with the TPI layer supplying its four `_tpi` swaps
(`budget_residual_D/F`, `external_account_D`, `domestic_bond_clearing`) via
`tpi.tpi_overrides()`. Verified bit-identical `main.py` output before/after
(no-op). Done ahead of Task 1 of `docs/superpowers/plans/2026-08-05-nominal-rigidities.md`,
which adds six new blocks and would otherwise need editing in three places.

## Nominal rigidities (`add-nkpc`, in progress): Task 7 — calibration parameters

`code/calibration.py` gains the parameters Tasks 2–6 already reference:
`mu_p_D/F = 1.20` (gross price markup, `epsilon_p = 6`; free to first order
under subsidy neutralisation, so the level needs no defending yet),
`kappa_p_D/F = 0.0871` (Rotemberg slope from Calvo `theta_p = 0.75` at
`beta = 0.985`, i.e. `(1-theta)(1-beta*theta)/theta`; matches Bi-Foerster-Traum's
implied 0.0846 to within 3%, consistent with the ~4-quarter euro-area IPN
median price duration), `pi_D/F = 0.0` (SS PPI inflation, exact), and
`omega_pi_D = 0.071` (union-inflation weight on D = 1 − the renormalised
two-country capital key, BuBa 26.1 / BoG 2.0 — deliberately *not* model GDP
weights, which would split the terms-of-trade move ~50/50 and erase the
93/7 internal-devaluation pattern Task 5 derived). **`mc_D`/`mc_F` retargeted**
from a dead placeholder `1.0` (nothing read them before Task 2) to
`1.0/mu_p = 0.8333`, the subsidy-neutralised SS real marginal cost that makes
labour demand collapse to the competitive `w=(1-alpha)Y/N` and keeps the SS
bit-identical to flex-price. Verified: `mu_p_D*mc_D == 1` and `mu_p_F*mc_F
== 1` exactly; no duplicate dict keys (checked by AST walk, given `mc_D/F`
already existed); `kappa_p` arithmetic reproduces 0.0871 to 5e-5.
`code/test_nkpc_blocks.py`: still 11 passed, 0 failed (calibration values
only, no new tests). Next: Task 8, seed `mc`, `pi`, `profit` into the steady
state.

## Nominal rigidities (`add-nkpc`, in progress): Task 6 — markup rent reaches households

`income_D`/`income_F` now add `profit_D*e_grid_D`/`profit_F*e_grid_F` to the
household's labour income (alongside `div`/`div_fund`), distributed in
proportion to productivity `e` rather than lump-sum. Because `w*N*e +
profit*e = mu_p*mc*(1-alpha)*Y*e + (1-mu_p*mc)*(1-alpha)*Y*e =
(1-alpha)*Y*e` exactly, household labour-plus-profit income is identical to
the flex-price model — the markup wedge (Task 4) bites only on the firm's
hiring decision, never on household income — so `labor_market_D/F` (labour
supply) stays valid and untouched. Zero at SS (`mu_p*mc=1`), so the steady
state is still bit-identical. `income_D/F` are hetinputs, not `@simple`
blocks, so adding `profit_D/F` to the function signature was sufficient to
wire them into `hh_extended_D/F.inputs`. `code/test_nkpc_blocks.py`: 11
passed, 0 failed (unchanged — no new tests needed for this task).

## Nominal rigidities (`add-nkpc`, in progress): Task 5 — global closure, no policy rate

`firm_profit_D`/`firm_profit_F` (appended after `labor_demand_D/F` in
`equations_D.py`/`equations_F.py`) now exist: `profit = (1 - mu_p*mc)*(1-alpha)*Y`,
the residual left over once labour demand pays `mu_p*mc*(1-alpha)*Y` (Task 4)
while capital keeps earning `alpha*Y` unchanged. Distributed to
households in proportion to productivity `e`, not lump-sum — this is what lets
`labor_market_D/F` stay untouched. Zero at SS (`mu_p*mc=1`), so current SS is
bit-identical.

`price_nkpc_D`/`price_nkpc_F` (appended after `firm_profit_D/F`) add the
Rotemberg curves themselves: `pi = beta*pi(+1) + kappa_p*(mu_p*mc - 1)`. The
gap is the ratio `mu_p*mc - 1`, not a level difference — unit-free, and
linearises to exactly `mc_hat` for any `mu_p`, so published Calvo slopes map
straight onto `kappa_p` with no SS rescaling. Subsidy-neutralised
(`mc_ss = 1/mu_p`) so `pi_ss = 0` and the steady state is still bit-identical
to the flex-price model. Discounted at constant `beta_D/F`, not the SDF —
immaterial to first order since `pi_ss = 0`. `kappa_p -> inf` recovers
flexible prices (`mu_p*mc = 1`), the equivalence gate for later tasks.

`labor_demand_D`/`labor_demand_F` now solve `w = mu_p*mc*(1-alpha)*Y/N` in
place of the competitive `w = (1-alpha)*Y/N`. Employment is no longer purely
supply-determined: previously, combining flex labour supply
(`w/P_CES = vphi*N^(1/frisch)`) with the old competitive demand condition
eliminated `Y` entirely, leaving `N` pinned by `Z`, `K`, `P_CES` alone — nothing
about aggregate demand entered. With the wedge, a change in `mc` shifts labour
demand, so `N` and `Y` move together. SS is unchanged bit-for-bit (`mu_p*mc=1`
at `mc_ss=1/mu_p` collapses the condition identically to the old one).
`labor_market_D/F` (labour supply) deliberately untouched, per the routing
argument above. `code/test_nkpc_blocks.py`: 8/8 tests green (2 `firm_profit`
+ 4 `price_nkpc` + `test_labor_demand_collapses_to_competitive_at_ss_markup`,
which turns the previously-red factor-exhaustion test green).

**Task 5 (`terms_of_trade` + `union_inflation`, appended to
`equations_global.py` after `bond_yield`) closes the nominal side with no
policy rate anywhere in the model.** `terms_of_trade` turns the identity
`p/p(-1) = (1+pi_F)/(1+pi_D)` (fixed nominal exchange rate in a monetary union,
so terms-of-trade movement *is* the inflation differential) into a residual on
`p`, an unknown that already exists (`goods_mkt_D`'s target). `union_inflation`
supplies the missing level normalisation — `omega_pi_D*pi_D +
(1-omega_pi_D)*pi_F = 0`, the `phi_pi -> inf` limit of an ECB Taylor rule on
union-wide PPI inflation, stated as an abstraction and not a modelled rule; no
Fisher relation is needed since no financial contract in the model carries a
policy rate. Solved together: `pi_D = -(1-omega_pi_D)*dlog p`, `pi_F =
omega_pi_D*dlog p`. At the capital-key `omega_pi_D = 0.071`, **93% of any
terms-of-trade move is Greek PPI deflation, 7% German inflation** — the
2010-12 internal-devaluation pattern; GDP weights would give ~0.5/0.5
(counterfactual) since the model normalises `Y_D_ss ~ Y_F_ss ~ 1`. Both
residuals zero at SS (`p/p(-1)=1`, `pi_D=pi_F=0`). `code/test_nkpc_blocks.py`:
11/11 green. One test-authoring bug found and fixed in review: the added
`test_closure_puts_93pct_of_tot_move_into_D_deflation` originally asserted
`log((1+pi_F)/(1+pi_D)) == dlog_p`, which only holds to first order — the
closed-form pi's are a linear approximation to the *exact* nonlinear
`tot_res`, and the O(dlog_p^2) truncation gap is `0.429*dlog_p` in relative
terms (derived: the gap is `-dlog_p*(pi_F+pi_D)/2`, and `pi_F+pi_D =
(2*omega-1)*dlog_p = -0.858*dlog_p` at `omega=0.071`), which swamped the
test's own `rel=1e-6` at `dlog_p=1e-4`. Fixed by asserting on the net rates
instead — `pi_F - pi_D == dlog_p` and `share_D = |pi_D|/(|pi_D|+|pi_F|) ==
1-omega` — both exact arithmetic with no linearisation, so the test is
robust to any `omega_pi_D` (a parameter swept later in this plan) rather than
tuned to dodge a truncation error. `terms_of_trade`/`union_inflation`
themselves were correct throughout; only the test needed fixing. Next: Task 6,
route the markup rent into `income_D/F`.

## Policy experiments (`experiments/`) — **COMPLETE 2026-08-03**

E1, E2 and E3 all land, with the orchestrator writing
`docs/experiments_results.md`. **Two results change how the paper must be
written** — E2's ΔY decomposition (report the decomposition, never the headline)
and E3's writeoff finding (full writeoff inverts Live Claim 1). Two items need an
author decision before anything is quoted: **S-1** (now quantified) and **A5-1's
misnamed third object**.

**Production regression re-run after all of it, `code/main.py`, exit 0, 18m46s —
bit-identical to the pre-work baseline:** `goods_mkt_D = -4.2493506589857954e-07`,
`ca_res_D = 6.852157730108388e-17`, `K_D = 10.800` / `K_F = 10.832`,
`n_inter_D[0] = -3.3804%`, `Y_D[0] = -0.0149%`, `b_gov_D[499] = 1.4e-05`, peak
spread +0.376 pp, loading 4.35/4.01/3.44. `git diff main -- code/` is **empty** —
the package has no side effects on production, by construction.

Full test suite: **31 passed** (`experiments/` 13, `diagnostics/regimes/` 8,
`code/test_eba_calibration.py` 10).

### S-1 RESOLVED (2026-08-04): `writeoff_enabled=0` — pure risk-premium framing

Author decision, taken with E3's numbers in hand. The paper commits to the
risk-premium framing: `def_rate` is a genuine probability, agents price the
expected loss, the IRF traces the **no-default branch**. Standard
risk-premium-shock device — *not* "default is impossible", and the paper must say
so explicitly.

E3 stays as an **appendix robustness result**, and its content is a caveat the
paper owes the reader: the over-compensation claim is **conditional on no realised
principal writedown**. With `zeta_writeoff=1` the loading falls to 0.37/0.28,
below 1. State the conditionality; do not bury it.

### E4 — distributional incidence (2026-08-05). DIST-1 addressed; Ginis dropped

`experiments/e4_distribution.py` adds per-bin consumption hetoutputs to the D
household block and re-solves the Jacobian (~5 min). **No Gini is computed** —
per DIST-1 it is the wrong statistic for the Greek crisis. Author decision
2026-08-05: report quantile responses instead.

**A methodological trap, found and avoided — read before adding any other
distributional cut.** Binning on *wealth* (steady-state deposits) with fixed
boundaries produces a per-capita consumption response that is **overwhelmingly
composition, not behaviour**. Measured, bottom decile, PV over 40q: the
consumption term is **−41.6** and the mass term **−44.4**, netting **+2.8**. The
whole deposit distribution shifts across fixed thresholds, so bin membership
churns and the net is a small residue of two large nearly-cancelling terms. The
arithmetic is exact and the object is well defined, but it must never be described
as "how poor households responded".

**The fix is to bin on the exogenous income state.** The marginal distribution
over productivity is the stationary distribution of the exogenous Markov chain, so
each bin's mass is **invariant** to the shock — households move between states but
the share in each state does not. The per-capita response is then purely
behavioural. **Income quintiles are the cut to report**; wealth deciles are
retained in the tables with the caveat attached and are deliberately not plotted.

### Paper first-draft outputs (2026-08-04)

`experiments/paper_outputs.py` → five captioned figures in `experiments/paper/`
and three tables in `docs/paper_draft_results.md`. Every number derived live; only
source citations are literal text.

**Corrected here and in CLAUDE.md: the default-loading split is 3.1% fundamental /
96.9% collateral friction** (`EL_price_D=0.056134`, `psi_spread_D=1.737724`), not
the 10.9%/89% recorded at the pre-EBA calibration. This is Live Claim 3's
quantitative core and 96.9% is a materially stronger version of it.

**New finding from the transmission figure — the backstop does not shift the whole
spread path down.** Its cushioning is concentrated at impact; by roughly quarter
four the net-worth and investment paths across regimes have converged, and the
spread ordering **reverses** (t=8: passive +0.5bp, medium +12.6, aggressive
+15.5; t=12: passive −10.1bp, aggressive +5.9). Passive overshoots downward later
while intervention decays monotonically. The right description is that the
backstop damps the *oscillation* — the initial spike and the later undershoot —
not that it uniformly lowers the spread. Any claim about intervening-regime paths
beyond the first few quarters must reckon with this.

Figure palette re-validated with the dataviz six-checks validator: the project's
previous red `#8C1515` **failed** the lightness band and the navy `#002147` failed
both the band and the chroma floor (it reads gray). The paper set uses
`#1B6CA8 / #A62B22 / #c87941 / #1a6e3a`, all passing in light mode.

### Details, newest first

New `experiments/` package on branch `experiments`, building the paper's standard
results set on top of `diagnostics/regimes/regime_model.py`'s solve/cache layer.
`code/main.py` is untouched and remains the regression path. Design spec:
`docs/superpowers/specs/2026-08-01-policy-experiments-design.md`.

**Cache schema v2 (this commit).** `regime_model.cache_path` computed the
calibration fingerprint at **import** time. A later experiment (E3) runs the model
under a calibration override applied at *run* time, so its cache would have been
written to the baseline filename and silently clobbered it — the same
stale-cache-passing-every-check failure the fingerprint was introduced to prevent,
reintroduced one level up. The fingerprint is now computed at call time, stamped
into the `.npz`, and asserted on load; `CACHE_SCHEMA` is in the filename so an
output-list change can never be masked by an unchanged calibration. New cached
outputs: `Phi_D` and `def_rate_D` (`T_D` is OPTIONAL — identically zero at
`T0=T1=0`, so zero-fill is the correct value, not a hole).

Schema-1 `.npz` files are now unused rather than overwritten; the cache must be
rebuilt before any experiment runs.

**Hardened after code review (2026-08-03).** `common.calibration_override` now
raises `KeyError` on an unrecognised override key instead of silently adding it as
junk while the intended (typo'd) parameter stays at its default — the same
wrong-but-plausible-number failure mode this package exists to prevent, one layer
up. `common.write_results` now refuses to serialise `NaN`/`Infinity`
(`allow_nan=False`) rather than writing a token most JSON parsers reject.

**Cache rebuilt under schema 2 (2026-08-03), 7m27s.** `G_tpi[cb=0]` vs the
baseline Jacobian `max|err| = 0.00e+00`; `A_cb[0,0] = -1.88891e-02` (backstop
compresses); `n_inter_D[0] = -3.3804% of SS`, `Y_D[0] = -0.0149% of SS`,
`b_gov_D[499] = 1.4e-05`, `ρ_b = 0.8451` — all identical to the committed
`code/main.py` baseline. SS block residuals all `< 1e-8`; SS goods residuals
`goods_mkt_D = -4.25e-07`, `goods_mkt_F = -4.19e-07`, `ca_res_D = 6.85e-17`,
matching the baseline run exactly.

> Two outputs are **zero-filled, and verified so rather than assumed**: `Phi_D`
> has no Jacobian column (the portfolio adjustment cost is quadratic about its
> anchor, so its level deviation is second-order and vanishes to first order), and
> `G_D` is absent from `G_tpi.outputs` entirely (government spending is constant).
> Both are genuinely zero, which is why E2's identity still closes.

> **Stale number, noted not fixed.** `EL_price_D` at the live calibration is
> **0.056134**, not the `0.0717` quoted in this file's older sections and in
> CLAUDE.md. 0.0717 was computed at the pre-EBA `delta_b=0.10` and `q_b=0.83`;
> the live EBA values are `delta_b_D=0.0777006` and `q_b_D=0.968941`, giving
> `(1-0.30)*0.0777006/0.968941 = 0.056134`. Anywhere `EL_price` is quoted as a
> denominator (the TPI loading), re-derive it rather than reusing 0.0717.

### E2 — ΔY decomposition (2026-08-03) — **DONE**

`experiments/e2_dy_decomposition.py`. Decomposes the output response against the
`market_clearing_D` identity `Y_D = P_CES_D·C_D + I_D + G_D + Φ_D + T_D + NX_D`.
Self-verifying, since `goods_mkt_D` is a targeted residual: **closure achieved at
5.8e−17 / 1.1e−16 / 1.5e−16** across the three regimes, against a 1e−7 assertion.

Impact (t=0) level deviations, by regime:

| component | passive | medium | aggressive |
|---|---|---|---|
| consumption (quantity) | +1.462e−03 | +2.055e−03 | +2.606e−03 |
| consumption (price) | +1.763e−04 | −3.416e−04 | −7.955e−04 |
| **investment** | **−1.868e−03** | **−7.025e−04** | **+2.945e−04** |
| **net exports** | **+7.955e−05** | **−9.005e−04** | **−1.767e−03** |
| government / portfolio cost / macro-pru tax | 0 | 0 | 0 |
| **dY[0] total** | **−1.494e−04** | **+1.107e−04** | **+3.377e−04** |

γ = 0 / 5.0798 / 12.7260. `dY[0]` as % of SS: −0.0149 / +0.0111 / +0.0338;
trough −0.0149 / +0.0065 / +0.0043.

**The watch item is answered, and the answer is the SPEC's warning made
quantitative.** Going passive → aggressive, the headline `dY[0]` moves by
**+4.87e−04**, but that is the residue of an investment swing of **+2.16e−03** and
a net-export swing of **−1.85e−03** — each roughly **4× the headline**, pointing
opposite ways. Consumption behaves the same way (quantity +1.14e−03 against price
−0.97e−03). `docs/SPEC.md`'s "a small headline output number can be *only* small
because two large channels are netting out — and they land on different
households" is confirmed in-model. **Do not report ΔY as a headline; report this
decomposition.**

On whether the positive intervening-regime output is linear-rule overshoot: the
proximate driver is **consumption quantity, not investment**. At `medium`,
investment is still negative on impact (−7.025e−04) and `dY[0]` is positive only
because consumption (+2.055e−03) outweighs it and NX. Investment does not turn
positive on impact until `aggressive`. Consumption is *already* positive on impact
at `passive` (+1.462e−03) with no backstop at all, so the consumption response is
not created by the policy rule — which argues against pure overshoot. But the
magnitudes are 0.01–0.03% of SS, i.e. the small difference of much larger
offsetting terms, so the headline is not robust and should not be leaned on.

### E3 — S-1 writeoff (2026-08-03) — **DONE. This one changes a paper claim.**

`experiments/e3_writeoff_s1.py`. Two nested variants, because the two switches do
different things. All rows evaluated at the **baseline's** γ (0 / 5.0798 /
12.7260) held fixed, so the model changes without the policy also changing.

| setting | `writeoff_enabled` | `zeta_writeoff` | `EL_price_D` | peak spread @γ=0 | loading (medium) | loading (aggressive) |
|---|---|---|---|---|---|---|
| baseline | 0 | 0.0 | 0.056134 | 150.3 bp | 4.00 | 3.17 |
| **E3a** coupon-only | 1 | 0.0 | 0.056134 | 149.1 bp | 3.93 | 3.13 |
| **E3b** full | 1 | 1.0 | **0.701743** | **168.9 bp** | **0.37** | **0.28** |

**E3a is economically negligible.** Haircutting the *coupon* alone moves the peak
spread by −1.2 bp and the loading by −0.04. The coupon is only `delta_b ≈ 7.8%` of
the bond, so writing it down barely registers. **Its SS drift is exactly
0.000e+00**, confirming `writeoff_enabled` is strictly steady-state-neutral as
predicted.

**E3b inverts SPEC Live Claim 1.** The loading falls from 4.00/3.17 to
**0.37/0.28 — below 1**. Under full writeoff the CB is *under*-compensated,
receiving roughly 30% of the actuarially fair expected loss, where the paper's
central claim is over-compensation. **"The monetary-financing objection fails on
the model's own terms" does not survive `zeta_writeoff = 1`.**

**The mechanism is clean and attributable to the denominator alone.** At `medium`,
premium income barely moves (`prem_PV` 0.00909 → 0.00988, +9%) while the priced
expected loss explodes (`EL_PV` 0.00228 → 0.02696, **×11.8**). The inversion is a
*repricing of the expected loss*, not a change in what the CB earns.

**`psi_lambda_B = 8.5` no longer hits its target.** Peak spread goes 150.3 →
168.9 bp (+12.4% over the 150 bp anchor) under E3b; E3a stays on target at
149.1 bp. Reported, **not** re-tuned away — whether to re-tune is a separate
author decision this result informs.

> **Correction to the design spec's prediction.** The spec said E3b "moves the
> steady state via `EL_price`". More precisely: `EL_price` (a derived anchor)
> changes value 12.5×, but **no steady-state allocation moves** — measured SS
> drift is exactly 0.000e+00 across `q_b`, `Y`, `C`, `I`, `NX`, `n_inter`, `K`,
> `TAX`, `P_CES`, `b_gov`, `b_D_D`. `EL_price` multiplies `def_rate`, which is 0
> at SS, so it is SS-neutral in *allocation* while still changing the linearised
> bond FOC and hence every dynamic result.

**New finding, not anticipated: under E3b the named-regime construction itself
breaks.** Peak spread stops being monotone in γ, so `gamma_for_compression`'s
bisection is invalid and "25% / 50% compression" no longer identifies a unique γ.
Two violations over a 40-point grid on γ ∈ [0,15]: a negligible one at γ≈0.385
(+0.7 bp) and a **large spike at γ≈3.46 (144.4 → 166.6 bp, +22 bp)**, after which
the curve resumes falling monotonically to 82.7 bp at γ=15. The isolated spike
sits where the closed-loop matrix `I − γ·A_cb` is likely near-singular, so treat
it as a **linear-algebra pathology rather than economics** until confirmed — but
it does mean compression-targeted regimes cannot be defined under full writeoff.

`A_cb` impact stays negative throughout (−1.889e−02 / −1.894e−02 / −2.158e−02),
so the backstop still compresses on impact in every setting.

### Orchestrator (2026-08-03)

`experiments/run_all.py` → `docs/experiments_results.md` (generated; do not
hand-edit). `--skip-e3` avoids E3's two re-solves; `--render-only` rebuilds the
document from results already on disk. Every table carries a live provenance
stamp including a **working-tree-dirty flag**, so a document produced with
uncommitted edits cannot be mistaken for a clean run at that SHA.

### E1 — backstop schedule (2026-08-03) — **DONE**

`experiments/e1_backstop_schedule.py`. Named regimes canonical, γ **solved** for
0/25/50% peak-spread compression.

| | passive | medium | aggressive |
|---|---|---|---|
| γ | 0 | 5.0798 | 12.7260 |
| peak spread (bp ann) | **150.3** | 112.7 | 75.2 |
| `Y_D[0]` (% SS) | −0.0149 | +0.0111 | +0.0338 |
| `C_D[0]` (% SS) | +0.2164 | +0.3040 | +0.3855 |
| `I_D[0]` (% SS) | −0.7718 | −0.2903 | +0.1217 |
| `n_inter_D[0]` (% SS) | −3.380 | −2.167 | −1.099 |
| **loading** | n/a | **4.00** | **3.17** |

**Cross-checks against the independent production pipeline all pass.** `code/main.py`
gives loading 4.01 at γ=5 (E1's `medium` is γ=5.08 → 4.00) and 3.44 at γ=10 (E1's
`aggressive` is γ=12.73 → 3.17, correctly below). Peak spreads hit the compression
targets exactly (150.3 → 112.7 = 25%, → 75.2 = 50%), and `n_inter_D[0]` /
`Y_D[0]` reproduce this file's own regime table to every printed digit.

**Live Claim 5 (self-extinguishing premium) confirmed on a fine grid.** The
loading schedule is **monotone decreasing at all 59 finite grid points**, 4.51 at
γ=0.51 down to 2.07 at γ=30. This is stronger than the three-point evidence
previously on record.

**A5-1, three separate objects (never summed):**

| regime | exposure PV (% Y) | expected loss PV (% Y) | `pd_D` differential PV |
|---|---|---|---|
| medium | 0.4265 | 0.00228 | −0.001541 |
| aggressive | 1.1062 | 0.00391 | −0.004651 |

> **The third object's sign needs an author decision, and "Greek fiscal saving" is
> the wrong name for what is computed.** The code reports
> `Σ β^t (pd_passive − pd_intervention)`, which is **negative** because the
> backstop lets Greece run a *larger* primary deficit — i.e. it relaxes the
> required austerity. So a negative number here means Greece is better off, which
> is the opposite of what "saving" implies. Either flip the sign convention or
> rename it (e.g. "austerity relief, PV"). **Do not report this number under its
> current label.** The magnitudes (0.0015 / 0.0047 PV) are unaffected.

`carry_ss_pv` is 1.2e−16 / 3.2e−16 — numerically zero, confirming SS yields are
equalised and that the `delta_b_F` fix below is correct.

**Welfare (SECONDARY — SPEC says do not lead with it):** `W_D` +0.0399 / +0.0646 /
+0.1056, `W_F` +0.0403 / +0.0189 / −0.0310. Near-exactly zero-sum, as before.

> **Cache schema 3 (2026-08-03).** `delta_b_F_ss` added. E1's `cb_pnl` computes
> each carry leg's SS yield as `delta_b·(1/q_b_ss − 1)`, and the first draft used
> `delta_b_D` on **both** legs — but `delta_b_F = 0.056779 ≠ delta_b_D = 0.077701`
> (the two countries' bank books have different measured maturity ladders). That
> put the SS spread at −9.2e−04 instead of its true ~1e−17 and would have
> contaminated `carry_ss_pv`. Caught by an assertion written into `cb_pnl` before
> the code was first run. Cache rebuilt; E2 re-ran identically.

## EBA REBUILD (2026-07-31) — read this first

The EBA 2011 moment set was rebuilt from scratch to be identified rather than
back-solved (`code/eba_calibration.py`, `data/eba_moments.json`,
`docs/eba_calibration.md`). **The rebuild succeeded; feeding the result into the
model did not, and that is the finding.**

**The EBA calibration is now LIVE** — see *RESOLVED* below for the final state
and verification table. The two sections immediately following record the journey
(what was measured, and the two structural blockers found and fixed on the way);
they are kept because the failure modes are instructive, but the *RESOLVED*
section is what describes the committed model. Setting `EBA_CALIBRATION=False`
still reproduces the pre-EBA placeholder calibration bit-exactly for regression
comparison.

### What is now measured (was assumed or absent)

| Parameter | Was | Measured | Source |
|---|---|---|---|
| `delta_b_D/F` | 0.10 (no EBA counterpart) | **0.0777 / 0.0568** | maturity ladder repriced at the end-2010 market yield |
| `theta_D/F` | 4.0 assumed | **5.51 / 6.94** | (corp ex-CRE + CRE + sovereign) EAD / CT1 |
| `omega_K_D/F` | back-solved plug | **0.117 / 0.067** | corp+CRE EAD ÷ K, `K/Y_ann=2.7` |
| `n_inter`, `phi_*`, `B_supply` | placeholders | 0.4075/0.1748, 2.390/2.758/0.069/0.018, 1.116/0.483 | CT1 and sovereign books |

`delta_b` measurement **retires the F-1 duration blocker**: the standing
"port 0.036/0.038 (7y)" target measured the *sovereign's whole outstanding
stock*, not the *bank-held book at the yields banks faced*. The bank book's
modified duration is 3.12y (GGB) / 4.22y (Bund) — near the old 0.10, so
`mv_rule=1` + `phi_lamb=0.60` are not needed.

`omega_K` is **kept, not dropped**: banks fund ~12%/7% of the capital stock, so
`omega_K=1` is counterfactual. It is now measured, and `K` becomes an output —
the balance sheet delivers `K_D = 10.800` against the `K/Y=2.7` target of 10.8.

### The blocking finding: GK feasibility

The GK block is well-posed only if

```
f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross
```

At the measured moments this is **violated by −1.26 (D) / −1.42 (F)**, giving
`lambda_gk_D = -0.087`, `Omega_D = -0.301` — negative IC multiplier and negative
banker franchise value. **The solver still converged, every Walras residual was
machine-zero, the IC-δ check passed exactly, stability passed, and the TPI
loading still declined.** All of it meaningless. This is the C-1 failure mode
reborn: silent degeneracy passing every check the pipeline ran.

The constraint bounds `Delta_own > ~0.73`. What had pinned it at 0.2 was an
undocumented convention inside `ic_delta_calibration` —
`ratio = Delta_cross/Delta_own = 2.0` — which with `Delta_cross <= 1` capped
`Delta_own <= 0.5`. That convention is a back-solve closure masquerading as a
consistency check, and it is exactly why `0.2/0.4` "passed".

**RESOLVED — collateral mapping fixed.** The ratio convention is gone; `Delta` is
now free and `ic_delta_calibration` checks the IC *residual* directly (tol 1e-8,
plus a positive-divertable-leverage check). Sweeping:

| `Delta_own`/`Delta_cross` | `lambda_gk_D` | `lambda_gk_F` | `Omega_D` | `Omega_F` |
|---|---|---|---|---|
| **0.85 / 0.90** | **+0.927** | **+0.960** | +4.62 | +5.98 |
| 0.90 / 0.95 | +0.488 | +0.456 | +2.49 | +2.91 |

`lambda_gk_D = +0.927` is essentially the pre-EBA `+0.923`, so the amplification
block keeps its strength while concentration becomes measured. `K_D=10.80`,
`K_F=10.65` vs the 10.8 target. This is a correction, not a fudge: leverage of
5.5× on a 43%-sovereign book is inconsistent with sovereigns being good
collateral — if they were, the bank would lever further.

**Guarded.** `steady_state.assert_gk_well_posed` runs inside `_apply_ss_anchors`,
on every solved SS in both `steady_state.py` and `depreciation_calibration.py`.

### RESOLVED — `omega_K` and the `n_inter` scope; **the EBA calibration is LIVE**

Two fixes, in order:

1. **`omega_K` as a fixed share** made the passive fund hold `(1-omega_K)*K`
   always, so it mirrored bank deleveraging: `dK/dN = theta/omega_K` (47.1).
   New **`fund_rule_D/F = 1`** — fund holds a constant `K_fund`, bank is the
   marginal holder, `dK/dN = theta` (5.5). Steady state **identical** under both
   rules; purely dynamic. Helped 17× but did not stabilise alone.
2. **`n_inter` scope was the real blocker.** EBA CT1 is the capital of the
   *stress-test sample*; the model's `n_inter` is the net worth of the agent
   intermediating the *whole* capital stock. New **`BANK_SCOPE = "broad"`**:
   `n_inter = (Q*K + sovereign)/theta`, so **`omega_K = 1` by construction** and
   the fund device disappears entirely.

| | CT1 scope | **broad (live)** | pre-EBA placeholder |
|---|---|---|---|
| `n_inter_D/F` | 0.408 / 0.175 | **2.138 / 1.627** | 3.0 / 3.0 |
| `phi_own_D/F` | 2.390 / 2.758 | **0.456 / 0.296** | 0.25 / 0.25 |
| `omega_K_D/F` | 0.117 / 0.067 | **1.0 / 1.0** | 1.0 |
| `theta_D/F` | 5.51 / 6.94 (measured) | **5.51 / 6.94** | 4.0 |

Kept measured: `theta`, the sovereign book, `delta_b` (ladder), `K/Y`. Given up:
`n_inter` as observed CT1, and `phi_own = 2.39` as a *model* parameter (2.39 is
concentration *within the stress-tested slice*, not within the whole
capital-funding sector — only the latter is what the model's `phi_own` means).
**Load-bearing assumption:** applying the EBA sample's `theta` to the whole
sector. If non-stress-tested capital funding levers less, `N_broad` is larger and
`phi_own` smaller. This is now the *only* such assumption left in the bank block.

`Delta` returns to **0.2/0.4** (the 0.85 bound was a CT1 artifact; at
`phi_own=0.456`, `f*theta = 0.661 > 0.367`), the fiscal rule to `phi_lamb=0.15`,
`mv_rule=0`, and `psi_lambda_B` retunes to **8.5** for the 150bp target.

**Verified end-to-end** (`code/main.py`, exit 0):

| Check | Result |
|---|---|
| over-identifying `K` | `K_D = 10.800`, `K_F = 10.832` (target 10.8) |
| IC residual | −8.9e−16 (D) / 0.0 (F) at `Delta = 0.2/0.4` |
| Walras | `ca_res_D = 6.9e−17`; all block residuals < 1e−8 |
| stability | `b_gov_D[499] = +1.4e−05`, `ρ_b = 0.845` |
| `n_inter_D[0]` | **−3.380% of SS** ✓ (level dev −7.227 — see the units fix below) |
| `Y_D[0]` | **−0.0149%** ✓ — **Y-1 RESOLVED** |
| `rk_D`, `rk_F` | both exactly 0.010000 — **RK-1 RESOLVED** |
| peak spread (γ=0) | 0.376pp = **150.4bp**, on target |
| TPI loading | **4.35 / 4.01 / 3.44** at γ=2/5/10, declining |

The paper's self-extinguishing-premium claim survives on a properly identified
calibration.

> **Sweep-method caveat.** `psi_spread_D` is derived from `psi_lambda_B` inside
> `_apply_ss_anchors`, so a `psi_lambda_B` sweep **must re-solve the SS per
> point**; patching the flag on a solved SS leaves `psi_spread` stale and
> inverts the apparent sign of the spread response.

### Secondary results

- **Mechanical MTM channel measured** (Acharya–Steffen construction; the 2011
  adverse-scenario CT1 depletion is *deliberately rejected* — it excluded
  banking-book sovereign default): **−5.73%/100bp** of CT1 for GR banks,
  −39.8% of CT1 realised over 2011 (GGB 10y 12.01%→21.14%). The pre-EBA
  calibration generates only **−0.61%/100bp**, ~10× too weak — essentially all
  of the gap is `phi_own` (0.25 vs 2.39), not duration. This is what
  `psi_lambda_B=3.0` had been standing in for, and it reframes S-1's "89%
  collateral friction" split as a calibration artifact.
- **`rk_F` (RK-1) improves** to 0.009896 vs target 0.0100 under the measured
  moments (was 0.0133 in the 2026-07-22 build).
- **Regimes staleness fixed**: cache filenames now carry a calibration
  fingerprint (they keyed only on `psi_lambda_B`, so the `psilam=0` cache would
  have been silently reused across different models); the hardcoded
  `"market-value rule"` figure suptitle now reads the live `mv_rule`.
- **Units bug FIXED (2026-07-31)**: SSJ IRFs are *level* deviations, so `×100` is
  a percentage only where the SS level is ≈1. `Y_D_ss≈1` passes;
  `n_inter_D_ss=2.138` and `K_D_ss=10.8` do not. `main.py` printed
  `n_inter_D[0]×100` as `%` (a 2.1× overstatement — the widely-quoted "−7.227%"
  is the level deviation, the true impact is **−3.380% of SS**), and the Stage A
  figure titled two panels `(%)` on the same basis (2.1× on net worth, **10×** on
  capital). Both now divide by the SS level and print/plot `% of SS`, with the
  level deviation retained alongside for continuity with the historical logs.

### `psi_lambda_B` breakdown ceiling — re-derived 2026-07-31

`diagnostics/psilam_breakdown_sweep.py`, 16-point grid. The SS is solved once
(`psi_spread` is exactly linear in `psi_lambda_B`) and the Jacobian re-solved per
point with **both** dials moved together.

| `psi_lambda_B` | 8.5 | 14 | 15 | 18 | 20 | 25 | 26 | 27 | 28 |
|---|---|---|---|---|---|---|---|---|---|
| peak spread (bp) | 150.3 | 223.8 | 225.8 | 234.6 | 273.6 | 625.3 | 1034.5 | 8903.8 | 41.0 |
| `n_inter_D[0]` (% SS) | −3.38 | −4.82 | −4.63 | −4.33 | −5.29 | −15.47 | −27.72 | −264.9 | **+29.87** |

There is a **pole between 27 and 28** — the amplification denominator crosses
zero, the response runs up superlinearly and returns sign-flipped. The A7
>1000bp flag first fires at **26**. The *first* pathology is earlier and milder:
over `[14,18]` the net-worth response **shrinks while the peak spread keeps
rising**, an internally inconsistent doom loop.

`PSILAM_BREAKDOWN` is therefore set to **15.0** — from that first pathology, not
from the pole, so the guard has real margin (55% of the pole, 1.76× the live
8.5). It was 2.5, derived for CT1-scope net worth, and would have blocked the
live calibration outright. `run_regimes.py`'s A7 flag remains the empirical
backstop.

### Not done

- The `omega_K` cross-check against ECB BSI bank-credit statistics is unpulled.
- `beliefs.json` still dates from 2026-07-23 (estimated MS chain on the FRED
  peripheral–Bund composite; calibration-independent, so not stale in principle).
- `PT-1` (bank-net-worth-to-spread pass-through vs Acharya–Drechsler–Schnabl)
  was computed from the mis-scaled net-worth number and should be recomputed:
  at the live calibration it is **−2.25%/100bp** (−3.380% at 150.3bp), still
  inside the literature-implied −1.8% to −8.6%/100bp range but now at the low
  end rather than mid-range.

## Historical (2026-07-22 onward) — superseded by the section above

## CURRENT CALIBRATION (2026-07-30) — supersedes every table below

**The calibration was reverted to its pre-EBA values.** Everything below this section
describes the EBA-anchored calibration and is now **historical**. Structural fixes
(C-1, W-1/W-2/W-3, T-2, A-2, TPI-1, `omega_K`, capital-key conduit) are all retained —
only parameter values moved.

| Parameter | Value | Note |
|---|---|---|
| `psi_lambda_B_D/F` | 3.0 | pre-EBA default; gives 187.2bp ann vs 150bp target (~25% over) |
| `n_inter_D/F` | 3.0 (=0.75×4) | bank net worth, 0.75 of annual GDP |
| `omega_K_D/F` | 1.0 | capital fund empty (`div_fund=0`) — pre-EBA balance sheet exactly |
| `phi_lamb_D/F` | 0.15 | minimum stabilising value under the par rule |
| `mv_rule_D/F` | 0 (par) | see Finding F-1 below — **1 is unusable at `phi_lamb=0.15`** |
| `recovery_rate_D/F` | 0.30 | EL-1 **retained**; live only via `EL_price` while `writeoff_enabled=0` |
| `delta_b_D/F` | 0.10 | 2.5yr; empirical 7yr/6.5yr still not ported (needs `mv_rule=1`) |
| portfolio targets | 0.25/0.15/0.15/0.25 | `steady_state.py`; EBA's 2.39/0.018/0.069/2.76 retired |

**Verified end-to-end** (`code/main.py`, exit 0): `n_inter_D[0]=-3.0009%`,
`Y_D[0]=-0.0261%`; `b_gov_D[499]=-1.3e-5`; `rho_b=0.8451`; IC-δ exact;
`max|ca_res_D|=6.3e-8`, `max|goods_mkt_F|=1.1e-9`; TPI monotone in γ.

**Finding F-1 sharpened.** The `phi_lamb∈[0.15,0.18]` zone identified under `mv_rule=1`
is a **hard break**, not the "narrow, mild zone" described below. Measured directly:
`mv_rule=1` + `phi_lamb=0.15` → `n_inter_D[0]=-1554%`, `Y_D[0]=+0.170%` (perverse sign),
`b_gov_D[499]=1.6e-2`. `mv_rule=1` needs `phi_lamb=0.60` (→ `-5.89%` / `-0.024%` / `0.0`).

**Default-loading split.** `EL_price_D=0.0717` vs `psi_spread_D=0.8385` → fundamental
expected loss is **10.9%** of the default loading, GK collateral friction **89%**.

> **SUPERSEDED 2026-08-04.** Those are pre-EBA numbers. At the live calibration
> `EL_price_D=0.056134` and `psi_spread_D=1.737724`, so the split is **3.1%
> fundamental / 96.9% collateral friction** — a materially stronger version of the
> constrained-seller claim. Do not quote 10.9%/89%.

**Units.** `spread_rb` is a *quarterly* rate deviation; annualise ×4×1e4 for comparison
with the 150bp target.

**Policy regimes — RE-RUN 2026-07-31 at the broad-scope EBA calibration.** Caches
rebuilt from scratch (the calibration fingerprint forced it), Stage A, Stage B-lite,
the certainty-equivalence decomposition and 18/18 unit tests all pass (exit 0). All
three figures regenerated.

| | aggressive | medium | passive |
|---|---|---|---|
| `gamma` | 12.7260 | 5.0798 | 0 |
| peak spread (bp ann) | 75.2 | 112.7 | **150.3** |
| `Y_D[0]` (% SS) | +0.0338 | +0.0111 | −0.0149 |
| `n_inter_D[0]` (% SS) | −1.099 | −2.167 | −3.380 |
| A6 peak at `psi_lambda_B=0` | 3.78 | 4.67 | 5.44 |

`A_cb = -1.889e-2` — the backstop still **compresses**, so the ms-regime SA-1
pathology remains absent. The passive peak is 150.3bp, so `run_regimes.py`'s
120–180bp sanity band now **passes** (it did not at the pre-EBA 187.2bp). The
`gamma` values roughly doubled (12.73/5.08 vs 5.08/1.57) because `A_def` scaled with
`psi_lambda_B` while `A_cb` barely moved — more gain is needed per unit of
compression.

> **Watch item.** `Y_D[0]` is **positive** under both intervening regimes, and the
> A5 table's `dY_D` trough is `+0.0000` (aggressive) / `+0.0024` (medium) — output
> never falls at all. A backstop that fully offsets a sovereign shock on impact is
> a strong claim; at `gamma_aggressive=12.7` this is more likely linear-rule
> overshoot than economics. Do not report the intervening-regime output paths
> without checking this.

**A6 at the lottery stage — fixed 2026-07-31, and the invariance is real.** The old
check ranked the *full-sample* peak, which is the common pre-`k` spread (no branch has
acted yet), so it was identical across branches by construction and the strict
inequality compared floating-point noise: it read YES at `psi_lambda_B=0` on a
`+2e-14`bp gap and NO at the calibrated `psi_lambda_B` on a `-5e-13`bp gap — same
expression, opposite verdicts. Now ranked on the **post-revelation window `t>=k`** and
checked at **both** amplifier settings, with a `1e-3`bp separation margin so noise can
never produce a verdict:

| `psi_lambda_B` | aggressive | medium | passive | ordered |
|---|---|---|---|---|
| as calibrated (8.5) | 51.12 | 72.77 | 91.59 | YES |
| 0 (fundamental floor) | 2.51 | 3.09 | 3.59 | YES |

Separation at `psi_lambda_B=0` is 0.508bp, ~500x the margin. **A6 amplifier invariance
holds in the lottery as well as in deterministic Stage A** — the previous "false pass"
verdict was a measurement-window error, not a failure of the economics.

### What regime uncertainty actually is under first-order certainty equivalence

`diagnostics/regimes/certainty_equivalence.py` (new, 2026-07-31). The model prices
the first moment exactly and the second moment not at all, so it is fair to ask
whether the Stage B lottery is degenerate — i.e. observationally identical to one
*known* CB with the belief-weighted coefficient `gamma_bar`. It is not, but the
non-degeneracy is smaller than the Stage B figure suggests, and the two must not be
conflated.

Three constructions, all with the same silence-until-`k=2` timing so delay is not
doing the work: **CE** = one known CB at `gamma_bar`; **MIX** = belief-weighted
mixture of known-type economies; **LOT** = the actual lottery.

| beliefs | `gamma_bar` | CE | MIX | LOT | LOT − CE |
|---|---|---|---|---|---|
| onset (crisis-conditional) | 0.582 | 148.562 | 148.410 | 148.526 | **−0.036** |
| ergodic (unconditional) | 5.577 | 132.089 | 133.132 | 133.567 | **+1.477** |

(impact spread, bp ann.) The wedge is *exactly* decomposable — pre-revelation both
economies see the same `A_def@eps`, so
`LOT_0 − CE_0 = A_cb[0,:] @ (cb^e − cb_CE)`, verified to 6e−15 bp. At ergodic
beliefs the horizon split is `[2,5): +1.828`, `[5,12): +0.228`, `[12,40): −0.538`,
`[40,∞): −0.041`. **Mechanism: regime uncertainty back-loads the expected backstop**
(the lottery expects less buying in the first three quarters, more later), and the
date-0 spread weights near-term purchases most heavily, so the impact spread rises.

The decisive non-degeneracy check is that the impact spread is **non-affine in
beliefs**: 0.289bp of curvature over the medium→passive path. An affine profile
would make the lottery observationally equivalent to a known-`gamma` CB.

**Reporting discipline.** The 16.6bp move in `fig_stageB_premium_vs_pi.png`
(133.7→150.3 as beliefs go ergodic→pure-passive) is ~91% a shift in the
**conditional mean** of the backstop path — which linearisation prices correctly —
and only ~9% (1.48bp) the genuine uncertainty-vs-equivalent-certainty wedge.
Neither component is a risk premium. Call it the *expected-backstop* or
*anticipation* channel, never a "regime-uncertainty risk premium".

> A first version of this script asserted a negative `Cov_pi(gamma_s, spread^s)`
> and read the sign off the purchase *totals*. Both are wrong: the covariance is
> positive here, and the totals point the opposite way to the wedge because
> `A_cb[0,:]` decays in the horizon, so the two economies differ in **timing**
> while nearly agreeing in total. The exact identity above replaces the argument;
> the covariance is retained as a reported diagnostic only.

> **Note:** `audit_artifacts/` was removed 2026-07-30 — the harness carried its own
> hardcoded copy of the calibration instead of importing `get_calibration()`, so it
> silently tested a different model than `code/main.py`. Every `audit_artifacts/*`
> path cited below is **historical provenance** (what was run at the time), not a
> live command. Scripts recoverable from git history at `0c99013`.

## EBA calibration + C-1 fix (2026-07-22) -- current state supersedes the sections below

Everything from here to "What is complete (post-audit)" describes the state as of
the 2026-06-22 forensic audit. Since then (2026-07-13 to 2026-07-22), the
calibration was re-anchored to EBA 2011 stress-test data (see `docs/eba_calibration.md`
for the full derivation) and a structural bug (C-1) was found and fixed. Summary:

- **Calibration is now EBA-anchored, not the placeholder values below.**
  `n_inter_D/F`, `phi_bD_D_ss` etc., bank-held debt stocks, and bank concentration
  ratios all come from the EBA 2011 disclosure (`data/DATA_DISCLOSURE.CSV`), not
  the round-number placeholders in the table further down. See `docs/eba_calibration.md`
  for the parameter -> moment map.
- **A passive capital-intermediation fund (`omega_K`) was added**: banks hold only
  `omega_K` of the physical capital stock (`omega_K_D=0.0601`, `omega_K_F=0.0190`);
  the rest sits in a deposit-funded passive fund that rebates its spread to
  households (`capital_fund_D/F`, `div_fund_D/F`). Needed because EBA's thin bank
  net worth (`n_inter_D=0.408`) can't plausibly intermediate the entire capital
  stock. Verified Walras-neutral and a no-op at `omega_K=1`.
- **C-1 (`Delta_cross > 1`, degenerate IC) is fixed at its root, not calibrated
  around.** The single-asset `lambda_gk` formula in `steady_auxilliary_D/F`
  implicitly assumed full (`Delta=1`) bond divertability; substituted into the
  back-solve, this forced `Delta_cross > 1` whenever `Delta_own > 0.5` -- i.e.
  essentially always at realistic (EBA) concentration. `steady_auxilliary_D/F`
  now solve `lambda_gk` from the multi-asset IC directly (closed form, since
  `value` is linear in `lambda_gk`), taking `Delta_bD_D=0.2`/`Delta_bF_D=0.4` (D)
  and `Delta_bF_F=0.2`/`Delta_bD_F=0.4` (F) as genuine hardcoded calibration
  inputs instead of back-solved outputs. Verified to recover those exact values
  (to floating-point) via `code/ic_delta_calibration.py`'s consistency check and
  independently via `audit_artifacts/run_audit.py`. See `docs/eba_calibration.md`
  "Why C-1 forces Delta->1" for the full derivation.
- **The "explosive EBA doom loop" finding in `docs/eba_calibration.md` was
  superseded by the C-1 fix, not just improved.** With C-1 fixed and *no other
  calibration change* (`psi_lambda_B=0.31`, `mv_rule=1`, `phi_lamb=0.60` all
  unchanged from the pre-fix explosive runs), `b_gov_D[499]` is `~1e-5` on both
  shocks (was 3.4-79.5 across the variants previously tested) -- near-stationary,
  not merely damped. The explosive root was the C-1 degeneracy amplifying itself
  through the leverage/IC loop, not an intrinsic feature of the 2.39x EBA
  concentration ratio as previously concluded.
- **A second, independent bug was found and fixed alongside C-1**: `code/tpi.py`
  builds its own copy of the dynamic model and had not been updated for the
  `omega_K` capital-fund commits -- it was missing `capital_fund_D/F`, so
  `div_fund_D/F` was silently frozen at its SS level instead of responding to
  shocks in the TPI experiment. This produced a real (not machine-precision)
  accounting leak on every `gamma` (`max|ca_res_D|` up to 6.2e-2, comparable in
  size to the reported welfare gains). Fixed by adding the missing blocks; the
  `G_tpi[cb=0]` vs. baseline-Jacobian sanity check now matches to `0.00e+00`
  (was `2.54e-3`). **Any TPI welfare number from before 2026-07-22 is superseded.**
- **`audit_artifacts/run_audit.py` had drifted out of sync** with the shared
  equation files across three separate additions (`omega_K`, `mv_rule`,
  `EL_price` from the macro-pru-fix) and could not run at all without patching.
  Fixed: its calibration fixture now has the missing (neutral-default) params,
  its duplicated `_apply_ss_anchors` now imports the one true definition from
  `steady_state.py`, and it gained real regression assertions (Walras thresholds
  from this file's own acceptance table, the default-shock sign check from
  "Typical iteration" step 4, and a C-1 consistency check) that fail loudly with
  a collected summary instead of silently passing or crashing on unrelated
  `KeyError`s.
- **Two smaller items are flagged, not yet resolved** (see `docs/eba_calibration.md`
  "Remaining open items"): a `rk_F` depreciation-calibration miss (0.0133 vs.
  target 0.01, likely a pre-existing one-shot-calibration imprecision) and a small
  (two orders of magnitude below bank net worth) positive `Y_D[0]` on the default
  shock, plausibly a portfolio-substitution effect rather than a bug.
- **`psi_lambda_B` recalibrated 0.31 → 1.1284 (2026-07-22, later same day).**
  0.31 was never a calibration target -- it was chosen only to dodge the (now
  fixed) C-1 degeneracy and undershot the paper's external anchor (2010 GR-DE
  spread ~150bp on a 1pp default shock) by more than 3x (44bp). Re-ran the
  moment-matching exercise on today's model (`audit_artifacts/psilam_moment_sweep_postC1.py`)
  and found the spread-vs-`psi_lambda_B` response is smooth/monotonic only up
  to about 1.5-2.0, then turns wildly non-monotonic (219bp at 2.0, 97bp at 2.6,
  853bp at 2.8, 353bp at 3.0, 11478bp at 5.0) -- a linear-approximation
  breakdown, not real moments. **Both the pre-fix literature value (2.8, per
  `docs/SPEC.md`'s theoretical framing) and the original round-number default
  (3.0) now sit inside that broken region on this model and must not be
  restored as-is.** `psi_lambda_B=1.1284` verified directly to hit 151.3bp
  (smooth neighbourhood: 147bp at 1.10, 154bp at 1.15) and is now the
  committed value. Side effect: this also resolves a "loading <1" concern
  raised in a same-day paper-direction hostile review (the model's TPI
  compensation appeared to have inverted below the actuarially-fair level at
  the pre-recalibration `psi_lambda_B=0.31`) -- loading is back above 1
  (2.54/2.14/1.74 at γ=2/5/10) and still declining, restoring the paper's
  central over-compensation claim, at roughly a third the magnitude of the
  stale "~7x" figure in `docs/SPEC.md`'s theoretical framing, which has been
  updated accordingly.

## Policy regimes (exogenous backstop aggressiveness) — added 2026-07-23

Three exogenous ECB backstop regimes over the TPI feedback coefficient γ on the
endogenous D–F yield spread (`TPI_t = γ·(spread_t − spread_ss)`, dormant at SS),
plus a Stage-B lottery in which banks price an ex-ante-unknown CB type revealed at
date k. All post-Jacobian numpy on cached `G_tpi` response matrices — production
`main.py` untouched. Code/outputs in `diagnostics/regimes/` (`run_regimes.py` =
Stage A, `uncertain_regime.py` = Stage B, `lottery_math.py`+tests = the math core,
`regime_model.py` = the cache layer; see `regimes_log.md`). Beliefs disciplined by
the 3-state Markov-switching intervention-regime estimation in `Empirics/` (via
`beliefs_from_empirics.py`).

**The backstop compresses spreads on this model, as intended.** A unit CB purchase
gives `d(spread_rb)/d(cb_buy) = −1.95e-2` (compression): the ECB **capital-key
conduit** (`kappa_cb_F=0.929`) funds a D-bond purchase 92.9% through F's treasury,
so D banks shed bonds (`d(b_D_D)=−0.72`) and the periphery–core spread narrows.
Stage A (spec §7 spread-compression targeting): aggressive/medium γ = 17.05/5.72
hit 78/117 bp peak spread vs 156 bp passive, cushioning bank net worth; the A6
amplifier-invariance ranking **survives** at `psi_lambda_B=0`. Stage B: the impact
spread rises with the belief weight on a passive CB (fear of no backstop is priced
before any policy acts), with a positive regime-uncertainty premium (+52/+23 bp on
the aggressive/medium branches); the §10.3 structural identities hold exactly.

**Provenance caveat (important).** An earlier build of this feature on the retired
`ms-regime` branch reported the OPPOSITE — that CB purchases *widen* the impact
spread ("Finding SA-1"), forcing an output-protection workaround and an inverted
Stage-B sign. That result was an artifact of `ms-regime`'s **superseded** model
(single-country conduit, par-value fiscal rule, `psi_lambda_B=2.8` — now known to
sit in this model's linear-approximation-breakdown region). It does **not** hold on
main; the capital-key conduit + market-value rule resolve it. Any SA-1/SB-1
"spread-widening" or "output-protection" statement from `ms-regime` is void here.

## Historical state (as of 2026-06-22 forensic audit) -- see above for what changed since

## Current status

Six structural/accounting bugs were found and fixed in the 2026-06-11 forensic audit (W-1, W-2, W-3, T-2, A-2, TPI-1). See `docs/audit.md` for the full ranked finding list and `docs/verification_report.md` for verified fix status. All six are **merged into `main`** via PR #27 (`AB-audit`, merged 2026-06-11). PR #26 (`audit` branch) was closed as superseded. `main` was subsequently reorganised into modular Python files via PR #28 (merged 2026-06-22).

Core equations: `code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py`. `code/main.py` is the production entry point (added by PR #28; orchestrates `code/full_model.py`, `code/steady_state.py`, `code/tpi.py`, …). The legacy notebook `code/model_v12.ipynb` has been removed. TPI/IRF figures are produced by `code/tpi_plots.py` and `code/irf_plots.py` (no committed `plots/` dir).

`main` now contains all six fixes (PR #27) plus the modular-file reorganisation (PR #28). Use `main` — or a feature branch off it — for all new work. The `audit` / `AB-audit` branches are historical and should not be reused.

## What is complete (post-audit)

- Household deposit choice and GHH preferences for D (Greece) and F (Germany) households.
- Bank steady-state and intermediation blocks: capital, bond returns, fees, GK Bellman (P1) and IC constraint (P3/lambda_gk).
- Production, capital adjustment, and capital producer profit — W-1 fixed: production uses `Y=F(K_t)`, capital producer receives `mpk·(K−K(−1))` so all capital income is allocated.
- Deposit return predetermined correctly: `Rgross = (1+rdep(−1))·P(−1)/P` — T-2 fix. Funding legs in `bank_return_*` and FOCs in `intermediation_P1_*`/`divert_*` use ex-ante deposit rate.
- F-bank bond returns converted to F-goods via `p(−1)/p` in `bank_return_F` — W-2 fix (the dominant leak; was causing ~2% of F GDP goods_mkt_F residual on a 1% TFP shock).
- Cross-border bond FOC in F-bank uses `p/p(+1)` for expected return — W-3 fix (optimality condition; does not affect Walras but required for internally consistent portfolio choice).
- Smart steady-state blocks: `m = n·(1−(1−f)·(1+rn))` without spurious `+Phi+T` — A-2 fix (required for any `chi1≠0` calibration, e.g. bank-cal's chi1=0.5).
- Global goods market, external account, bond clearing, and portfolio adjustment cost blocks.
- Domestic and foreign bond pricing, yields, spreads, and Hatchondo-Martinez geometric-decay perpetuity default mechanics.
- TPI extension: CB budget closed via `budget_residual_D_tpi` remittance — TPI-1 fix. Before the fix, unbacked CB flows inflated welfare gains by ~40% at γ=10.
- **ECB balance sheet as a capital-key conduit** (ecb-balance-sheet branch): the TPI-1 wiring remitted the *entire* CB cash flow to the Greek treasury (self-financing; no German side). Now `cb_flow_D = (1+rb_actual_D)·q_b(-1)·cb_buy(-1) − q_b·cb_buy` splits `kappa_cb_F=0.929` to the F treasury (`budget_residual_F_tpi`, converted `/p`) and the rest to D; the F share of the CB book enters the external account like `b_D_F` (`external_account_D_tpi`), so the CB hole cancels identically and `ca_res_D`/`goods_mkt_F` stay at baseline leak levels at every γ. Carry/credit legs and off-path expected loss are post-processed in `cb_pnl` (`code/tpi.py`), discounted at `beta_F`; `audit_artifacts/tpi_test.py` imports the production blocks and asserts `CB_CONDUIT_TEST: PASS` at ≤1e−7.
- EBA bilateral sovereign exposures in calibration cell: b_D_D/asset=24.47%, b_F_F/asset=25.79%, b_F_D/asset=0.18%, b_D_F/asset=0.65%.

## Walras accounting (post-fix, verified)

| Residual | 1% TFP-D shock | 1pp default-D shock |
|----------|----------------|---------------------|
| goods_mkt_D (targeted) | ≤1e−16 | ≤1e−16 |
| goods_mkt_F (untargeted) | ≤8e−10 | ≤1e−9 |
| ca_res_D = CA−ΔNFA (untargeted) | ≤5.8e−8 | ≤3.5e−8 |
| deposit_mkt_D/F | ≤4e−15 | ≤4e−15 |

Pre-fix peaks for reference: goods_mkt_F 2.0e−2 (~2% of F GDP); ca_res_D 1.5e−4. All cross-country spillover and welfare results from the pre-fix model are first-order invalid and must be regenerated from `main`.

## IRF summary (historical, pre-EBA — phi_lamb=0.15 era; SUPERSEDED)

> **Superseded (2026-07-24).** The committed fiscal-feedback coefficient is
> `phi_lamb=0.60` (see the C-1-fix section at line ~37 and `code/calibration.py`),
> **not** the `0.15` this section's title assumed — that was a stale table value,
> now corrected. The impact numbers below also predate the current EBA /
> `psi_lambda_B` / `recovery_rate` calibration; the 2026-07-24 production run gives
> `n_inter_D[0]=−2.83%`, `Y_D[0]=+0.032%`, peak D–F spread `+0.392pp`,
> `b_gov_D[499]≈2e−6`.

**1pp default shock to D (ρ=0.8):**
- `n_inter_D[0] = −3.5%` (falls), `Y_D[0] = −2.5e−4` (falls) — both signs correct post-T-2-fix; were positive/perverse pre-fix.
- `n_inter_F[0] ≈ −0.33%` — contagion small, sign correct.
- Spread widens on impact; doom loop is live with correct sign.

**TPI (ECB capital-key conduit, ecb-balance-sheet branch, 2026-07-16):**
- ΔW_D = +1.38/+2.75/+4.13, ΔW_F = −1.36/−2.70/−4.03 (% quarterly SS consumption, 100q) at γ=2/5/10. Near-transfer, slightly positive-sum (+0.10 at γ=10). Spread now compresses with γ: peak 0.409pp (γ=0) → 0.161pp (γ=10).
- ECB P&L (PV at β_F over 100q, % of quarterly SS Y_D, γ=10): peak exposure 1.57%, purchases PV 4.20%, expected-loss leg 0.0070%, default-premium leg 0.0227%, MTM leg +0.0253%; SS-carry 0 (SS yields equalised).
- **Loading (premium PV / EL PV) declines monotonically in γ: 4.86 → 4.06 → 3.22** — the self-extinguishing premium (`docs/SPEC.md` "Theoretical framing," Live Claim 5) confirmed in-model at this (now superseded) calibration; `loading_arr` over `gammas_fine` in `tpi_results` is the key-figure schedule. **Superseded 2026-07-22**: current committed calibration gives 2.54/2.14/1.74 at γ=2/5/10 (see the EBA/C-1-fix section at the top of this file).
- Germany (κ=0.929): bears EL PV 0.0065% Y_D, receives premium PV 0.0210% at γ=10 (memo at full EA key 26.1%: EL 0.0018%).
- Conduit residuals at every γ: `ca_res_D` ≤ 7.4e−8, `goods_mkt_F` ≤ 6.2e−10. The pre-conduit "zero-sum transfer, no spread compression" numbers (ΔW_D +1.88/ΔW_F −1.90) were from the notebook-era model and are superseded.

## Calibration summary (current, main)

| Parameter | Value | Source / note |
|-----------|-------|---------------|
| `phi_lamb_D/F` | 0.60 | ~Bohn (1998) fiscal-feedback magnitude. Hand-set for stationarity under EBA doom-loop amplification (not moment-matched); enters only `tax_rule_D/F`, and under `mv_rule=1` reacts to the market-value debt gap (so it couples to `q_b`). Well above the F-1 near-unit-root zone `[0.15,0.18]`; stability confirmed (2026-07-24 run: `b_gov_D[499]≈2e−6`, `ρ_b=0.373`). Governs the debt/fiscal mode only — **not** the ~25q financial-accelerator ring (that mode is set by the PAC `psi_bF_D/psi_bD_F`). **Corrected 2026-07-24 from a stale `0.15` entry.** |
| `def_scale_D` | 0.25 | Strong amplification. Exceeds GR 2011 crisis peak (0.12–0.23 from spread-debt slope calibration). |
| `delta_b_D/F` | 0.10 | 2.5yr avg maturity. Empirically too short; bank-cal has 0.036/0.038 matching GR/DE 2011 ~7yr/6.5yr. |
| `theta_D/F` | 4.0 | GK leverage; conservative vs 2011 historical 10–25×. |
| `psi_lambda_B_D/F` | 3.0 | State-dependent divertability; primary amplification dial; no direct empirical counterpart. |
| `f_D/F` | 0.12 | Bank exit rate; bank-cal has 0.03 (standard GK range). |
| `kappa_cb_F` | 0.929 | F share of CB conduit cash flows: two-country renormalised euro-area capital key (BuBa 26.1 / BoG 2.0). TPI layer only; SS-neutral. Memo reporting also at full EA key 26.1%. |
| `Delta_cross` | 1.4545 | Back-solved (`_ic_delta`, ratio=2.0); degenerate >1. See C-1. |
| `recovery_rate_D` | 0.00 | No realized losses; inert while writeoff_enabled=0. |
| `writeoff_enabled_D/F` | 0.0 | Default produces no balance-sheet losses. See S-1. |
| `chi1_D/F` | 0.0 | Intermediation adjustment cost off. A-2 fix makes chi1≠0 safe. |
| `frisch` | 0.5 | Frisch elasticity. |
| nDep_D/F | 500/500 | Household deposit grid points. |
| income rho_z/sigma/nZ | 0.90/0.30/15 | D and F income process (Markov approximation). |

## Open issues

| ID | Description | Status |
|----|-------------|--------|
| C-1 | ~~`Delta_cross=1.45 > 1`: divertable fraction exceeds 1, making the multi-asset IC constraint degenerate.~~ | **RESOLVED (2026-07-22).** `lambda_gk` now solved from the multi-asset IC directly; `Delta_bD_D/F=0.2/0.4` are genuine hardcoded inputs, verified to bind exactly. See "EBA calibration + C-1 fix" above and `docs/eba_calibration.md`. |
| S-1 | `writeoff_enabled=0`: default shock produces zero realized bank losses. `recovery_rate` and `zeta_writeoff` are set but inert. Model is currently a pure risk-premium loop, not a balance-sheet doom loop. | Author decision. Resolution: set `writeoff_enabled=1` with `recovery=0.40, zeta=1.0` (GR 2012 ~50% haircut → ~0.4–0.5 recovery). |
| RK-1 | `rk_F` depreciation calibration misses its target (0.0133 vs. 0.0100) post-EBA; the one-shot `delta_F` calibration in `depreciation_calibration.py` isn't iterated to convergence against the re-solved SS. `rk_D` hits its target exactly; only F misses. | Open (2026-07-22). Not yet investigated against a pre-EBA baseline. |
| Y-1 | Small positive `Y_D[0]` on the default shock (`+0.01%` to `+0.03%` across calibrations) even though `n_inter_D[0]` is correctly negative. Plausibly a portfolio-substitution effect (bonds are now good collateral at `Delta=0.4` and become relatively worse on impact) partially offsetting the direct hit — a tension already flagged in an earlier "substitution vs deleveraging" diagnostic. | Open (2026-07-22). Small magnitude; worth checking before reporting `Y_D` impact signs. |
| EL-1 | ~~`recovery_rate_D/F=0.00` is a placeholder~~ **RESOLVED (2026-07-22).** Set to `0.30` (NPV-recovery framing for the actual March 2012 Greek PSI): Zettelmeyer, Trebesch & Gulati ("The Greek Debt Restructuring: An Autopsy", PIIE WP13-8) find actual investor NPV losses of 59-65%, considerably below the ~75% commonly quoted; contemporary bank estimates (Credit Suisse, Morgan Stanley) put NPV haircuts at 73-78%. 0.30 recovery (70% haircut) sits centrally in this range. NPV framing chosen over the face-value framing (53.5% face cut -> ~46.5% recovery) because `EL_price` multiplies a payoff already priced at market terms (`q_b`), matching the NPV concept. `EL_price` fell from 0.1025 to 0.0717 as a result; `psi_lambda_B` re-tuned from 1.1284 to 1.1793 to restore the 150bp target exactly (verified: 150.02bp). Full pipeline reverified clean at the new calibration. | **Resolved.** Side effect: TPI loading rose to 3.59/3.03/2.47 at γ=2/5/10 (was 2.54/2.14/1.74) — a smaller `EL_price` denominator makes the same premium income look like a larger multiple of fair compensation. Update `docs/SPEC.md`'s theoretical framing numbers accordingly. |
| PT-1 | ~~Bank-net-worth-to-spread pass-through unvalidated~~ **CHECKED (2026-07-22).** Model gives ≈−4.5%/100bp at the calibrated point. Acharya-Drechsler-Schnabl (2014 JF, Table 6 col. 9, post-bailout, full controls) report bank equity return regressed on Δlog(sovereign CDS): coefficient −0.096 (s.e. 0.017, 1% sig.). Converting their log-CDS elasticity to a per-100bp-of-level moment requires a baseline CDS level (log vs. level don't translate 1:1): using their post-bailout sample mean (108.5bp) gives ≈−6.3%/100bp; using their reported typical weekly move (11.3% log-change at that mean, linearly rescaled to 100bp) gives ≈−8.6%/100bp; using a higher, more Greece-specific crisis baseline (~500bp) gives ≈−1.8%/100bp. **The model's −4.5%/100bp sits within this literature-implied range (−1.8% to −8.6%/100bp) under every reasonable baseline choice** — same order of magnitude, not a stray number. Altavilla-Pagano-Simonelli (2017, Rev. Finance) corroborates the *mechanism* (bank sovereign exposure amplifies transmission of sovereign stress to bank outcomes) but its headline results are on lending, not equity/net-worth, so it wasn't used for a second point estimate. | **Resolved as "within plausible range," not exact-match** — appropriate given the model has no free parameter tuned to hit this moment; it emerges from the calibrated collateral/leverage block. Re-check if `psi_lambda_B` or the bank-capital calibration changes again. |
| DIST-1 | Distributional resolution is Gini-only this iteration; per-decile deferred. `GINI_C` is specifically the wrong statistic for the Greek crisis — measured inequality barely moved (already highly unequal pre-crisis; worst-affected households dropped out of surveys). | Open. Any incidence claim in the paper must state its distributional resolution explicitly. |
| A5-1 | The TPI's German fiscal cost/benefit should be reported as three distinct objects — **exposure** (discounted purchases, what capital-key sharing acts on), **priced expected loss** (hand-computed off-path per `docs/SPEC.md`'s "structural constraints"), and **Greek fiscal saving** (`pd_D` differential vs. a passive counterfactual, the cleanest headline, no pricing assumption needed) — not conflated into a single "implicit transfer." Discount rate should be German risk-free for a cross-border object (currently `beta_F`, the creditor side — check this is deliberate and justified in the paper, not just a default). | Open; check `code/tpi.py`'s `cb_pnl` reports all three cleanly before writing this section of the paper. |
| X-1 | Dead-code imports in notebook cell 7: blocks no longer in the model remain in the import list. | Minor cleanup; no numerical effect. |

## Next priorities

**Superseded by the 2026-07-22 EBA calibration + C-1 fix (see top of file):**
EBA bilateral exposures are ported and live (item 1's exposures/Delta=0.2/0.4
portion is done); `delta_b`/`f` are **not yet** ported to the bank-cal empirical
values (still `0.10`/`0.12`) — that part of item 1 remains open. Figures were
regenerated 2026-07-22 under the fixed, EBA-anchored calibration (item 4 done for
now, but re-do again after any further calibration change).

1. **Port remaining bank-cal calibration values**: `delta_b=0.036/0.038` (empirical
   7yr/6.5yr duration), `f=0.03` (bank exit rate). Given C-1 is now fixed and the
   baseline doom loop is stable without it, re-test Finding F-1's explosive-root
   concern (below) under the *fixed* model before assuming it still applies —
   it was diagnosed on the pre-C-1-fix, non-EBA calibration.
2. **Decide S-1**: set `writeoff_enabled=1` to give default realized losses, or
   keep pure risk-premium framing and state it explicitly in the paper.
3. ~~Re-explore `psi_lambda_B`~~ **DONE (2026-07-22).** Recalibrated to 1.1284
   (150bp spread target), replacing the never-validated 0.31. See the entry
   above this table for the full rationale and the discovered non-monotonic
   breakdown region above `psi_lambda_B~2`.
4. **Investigate RK-1 and Y-1** (open issues table above) before reporting `rk_F`
   or `Y_D` impact-sign results in the paper.
5. **Re-generate all figures** from `main` after any of the above — the current
   `outputs/` figures (2026-07-22) reflect the EBA calibration + C-1 fix + TPI
   conduit-leak fix, but not any further calibration change.

## Finding F-1: market-value fiscal rule (duration ↔ fiscal-stability tension)

**Historical diagnosis below predates the C-1 fix** (kept verbatim; the doom loop's
amplification was resting on a degenerate collateral value at the time). **Re-tested
2026-07-22** on the now-fixed model, using today's actual committed calibration
(EBA-anchored, `psi_lambda_B=0.31`, `def_scale=0.25`, `writeoff=0`) with only
`delta_b` (duration), `phi_lamb`, and `mv_rule` varied. This went through three
rounds before landing on a trustworthy answer -- worth keeping all three, since the
methodological correction is as important as the result:

**Round 1 (coarse grid, energy-ratio proxy, `audit_artifacts/philamb_sweep_postC1.py`).**
`mv_rule=0` (par value): explosive at every `phi_lamb` in `{0.05,...,0.60}` tested,
but mildly (modulus 1.002-1.005, vs. 1.005-1.015 pre-fix; `bgov[499]` ~1e-6, vs.
9.5-79.5 pre-fix). `mv_rule=1` (market value): STABLE at 0.05, EXPLOSIVE at 0.10 and
0.15, STABLE again at 0.25-0.60 -- a non-monotonic "explosive island" sitting exactly
where the pre-fix finding's stable plateau (`phi_lamb` in `[0.07,0.12]`) was reported.

**Round 2 (finer grid, naive Prony/eigenvalue estimator, `philamb_sweep_postC1_fine.py`).**
The two-window energy-ratio proxy used in round 1 is known to be foolable by
oscillatory near-unit-circle modes, so a proper linear-prediction (Prony) estimate
was built: fit an order-`p` AR recursion to the IRF tail (t=150-500), take the
eigenvalues of the resulting companion matrix. Validated first against synthetic
signals (recovered known poles to ~1e-15) and an adversarial two-mode case (correctly
resolved a true dominant pole an energy-ratio proxy would have underestimated by 7x)
-- then applied at `p=8` (with `p=4` as a cross-check) to the real IRFs. **This gave
an implausible modulus of ~3.2** at every grid point, flatly contradicting the
observed `bgov[499]` levels (~1e-6) a true modulus of 3.2 could never produce over
500 periods.

**Root cause (`audit_artifacts/philamb_order_selection_check.py`).** R^2 of the AR
fit already saturates to machine precision (`1-1e-9`) at order `p=2-3`; every order
beyond that fits pure numerical noise into spurious high-modulus "ghost" poles that
contribute negligibly to the actual signal but have huge companion-matrix
eigenvalues -- textbook AR/Prony overfitting once model order exceeds what the data
supports. `p=4` (round 2's "conservative" cross-check) was *already* past this onset
for at least one grid point -- it wasn't a safe fallback either.

**Round 3 (corrected: automatic order selection, `philamb_sweep_postC1_fine_v2.py`).**
Refit every cell at `p=1..4`, taking the modulus at the *smallest* order whose R^2
clears `0.99999` (falling back to `p=3`, flagged, if none does -- never triggered).
Result, `mv_rule=1`, fine grid over `phi_lamb` in `{0.05, 0.075, ..., 0.25}`:

```
phi_lamb:  0.05    0.075   0.10    0.125   0.15    0.175   0.20    0.225   0.25
trusted:   0.997   0.982   0.998   0.996   0.999   1.003   0.990   0.990   0.990
verdict:   STABLE  STABLE  STABLE  STABLE  ~1.000  mild-X  STABLE  STABLE  STABLE
```

**Verdict: the "explosive island" was mostly an artifact of the energy-ratio
proxy.** The properly order-selected estimate is comfortably stable (modulus
0.98-0.997) across nearly the *entire* `phi_lamb ∈ [0.05, 0.25]` range under
`mv_rule=1` -- including 0.10 and 0.125, which round 1's proxy called explosive.
There is a real (not noise-floor -- this system's Jacobian conditioning implies a
~1e-10 floor, far below the ~0.003 excess seen here), but narrow and mild, near-unit-
root zone right around `phi_lamb ≈ 0.15-0.175` (modulus 0.999-1.003). That is much
smaller and narrower than either the round-1 or round-2 estimates suggested, and it
sits at a different location than the original pre-C-1-fix "stable plateau"
(`[0.07,0.12]`) claimed, but "avoid `phi_lamb` in roughly `[0.15,0.18]`, most other
values in `[0.05,0.25]` are fine" is a defensible, much more precise statement than
either "there's a stable plateau at 0.10" (pre-fix finding) or "there's an explosive
island at 0.10-0.15" (round-1 proxy). `mv_rule=0` (par value) was not re-swept at
finer resolution or order-selected precision -- round 1's mild-explosive-everywhere
reading for it should be treated with the same "proxy could be off" caveat pending
that redo. See `audit_artifacts/philamb_sweep_postC1_fine_v2_results.json` for the
full per-cell R^2/order/modulus breakdown, and `philamb_order_selection_check.py` for
the overfitting diagnostic.

*Historical diagnosis (pre-C-1-fix, kept for record):*

> The Bohn rule responds to the **par/face-value** lagged debt gap (`tax_rule_*`). With empirical long-duration bonds (`delta_b=0.036/0.038`, 7yr/6.5yr) this is **explosive at every `phi_lamb ∈ [0.02, 0.50]`** under both the balance-sheet (write-off ON) and risk-premium channels — the debt dynamics collapse to a near-unit-root ~250-quarter cycle (dominant modulus ≈ 1.005–1.015) that fiscal feedback cannot damp. The original model is stable at `phi_lamb=0.15` only *because* its duration is short (`delta_b=0.10`). So empirical duration, a literature `phi_lamb`, and a stable live doom loop are jointly infeasible with the par-value rule. (Sweeps: `audit_artifacts/philamb_sweep*.py`.)
>
> **Resolution — market-value rule.** Reacting to the mark-to-market debt gap `q_b·b_gov(-1) − q_b_ss·b_gov_ss` (it sees the current spread) opens a stable plateau at `phi_lamb ∈ [0.07, 0.12]` (modulus down to 0.983 at 0.12; `phi_lamb≈0.10` robustly interior) with empirical duration and a live, correctly-signed doom loop — **but only in the risk-premium framing** (`psi_lambda_B=1.0`, `def_scale=0.10`, write-off OFF). With write-off ON it is a *false victory*: `|λ|<1` but the default shock is perverse (spread narrows, bank net worth and output rise). So adopting the market-value rule **forces the risk-premium framing — it couples to the S-1 decision** (keep write-off OFF).

**Status.** Implemented as a switchable option: `mv_rule_D/F` in calibration (`0`=par, default, behaviour unchanged; `1`=market value). `mv_gov_ss_D/F` is set from the solved SS in `build_and_solve`. Adopting `delta_b` (empirical duration) as the baseline calibration remains an author decision. Per the round-3 (order-selected) re-test above, `phi_lamb` in roughly `[0.05, 0.125]` or `[0.20, 0.25]` reads as genuinely stable under `mv_rule=1`; avoid the narrow `[0.15, 0.18]` zone. Do not reuse the pre-fix `[0.07, 0.12]` "stable plateau" language as if it were re-confirmed — it wasn't specifically re-tested at the old paper's exact framing (risk-premium, `psi_lambda_B=1.0`, `def_scale=0.10`), only at today's committed EBA amplification (`psi_lambda_B=0.31`, `def_scale=0.25`), which happens to also land in a stable region there but for different reasons. **Committed value (2026-07-24): `phi_lamb=0.60`** — above the entire fine-swept range `[0.05,0.25]`. Round-1 (coarse) read `mv_rule=1` as stable at `0.25–0.60`, and the 2026-07-24 production run confirms it (`b_gov_D[499]≈2e−6`; debt mode well-damped, `ρ_b=0.373`). The `[0.15,0.18]` caveat therefore does **not** bind the committed calibration — it matters only if `phi_lamb` is ever lowered toward the literature range. Note this is the *debt/fiscal* dominant mode; the ~25q financial-accelerator ring in the asset-price IRFs is a separate, faster complex pair governed by the PAC (`psi_bF_D/psi_bD_F`), not `phi_lamb` — see `audit_artifacts/pac_sweep.py`.

## Finding F-2: financial-accelerator ring (2026-07-24)

The baseline IRFs carry a damped oscillatory "ring", concentrated in the
asset-price / bank-net-worth variables (`n_inter_D`, `q_b_D`, `q_b_F`, `C_D`) and
**absent** from the real block (`Y_D`, `w_D` decay monotonically). Quantified at
the committed calibration with a Prony / companion-eigenvalue extractor
(`audit_artifacts/pac_sweep.py`):

- **Complex pair: |λ| = 0.954, period ≈ 25q (~6.25 yr), half-life ≈ 14.6q
  (~3.6 yr).** R²=1.00000; consistent across `n_inter_D`/`q_b_F`/`q_b_D`, Prony
  orders 6 and 8, and the empirical peak-spacing (~24q). Comfortably damped —
  **not** a stability concern (each swing falls to ~30% of the last; ~99% gone by
  q100).
- **It is the GK IC/leverage financial accelerator, not a portfolio-friction
  artifact.** A PAC sweep (`psi_bF_D=psi_bD_F` over 0.05→5.0, a 100× range) moves
  |λ| by <0.006 (0.954→0.960) and leaves the period fixed at ~25q; the dependence
  is mildly U-shaped with the committed `PAC=0.5` sitting at the modulus
  *minimum*. So the ring is intrinsic to the amplification block
  (`psi_lambda_B`/`def_scale`/`theta`/`f`), which is pinned to the 150bp spread
  target — **there is no free cosmetic fix.** Three modes, three owners:
  debt/fiscal → `phi_lamb` (fast-damped, `ρ_b=0.373`); cross-border position →
  PAC; the ~6yr financial cycle → the accelerator. Neither the fiscal nor the
  portfolio dial touches the third.
- **Takeaway:** treat the ring as a structurally-identified ~6yr financial cycle
  to *describe* in the paper's propagation discussion, not a bug to patch. Full
  per-cell R²/order/modulus table: `audit_artifacts/pac_sweep_results.json`.
