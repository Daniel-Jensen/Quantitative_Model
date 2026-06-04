# Brief — Auclert intermediary cap-adj + macroprudential bond tax refactor

## Repo
- `/Users/Adam/Documents/Uni/PhD/Research/QUANTITATIVE_MODEL`
- Stack: 2-country monetary-union DSGE on `sequence_jacobian` 1.0.0. Gertler-Karadi banks with three ν's (`ν_K`, `ν_bD`, `ν_bF`) per country. Studying sovereign default crises.
- Notebook: `model_v12.ipynb` (single source of truth — calibration, SS solve, off-SS imports, `ha_full` assembly, IRFs, diagnostics).
- Equation modules: `equations_D.py`, `equations_F.py` (mirror), `equations_global.py`.
- Current branch: `main` at `08e1019` (clean). The refactor below is in progress, partly committed, partly uncommitted — read the files; do not assume the latest commit reflects the new state.

## What this refactor changes (manual edits already partly in tree)
The user is replacing the old `risk_weight × mp_wedge` intermediary friction with two cleaner pieces and tightening the SS anchors. Code-review and continue debugging until SS solves and impulse responses are sensible.

1. **Auclert (2019) intermediary capital adjustment cost `Phi`.**
   - `cap_adj_cost_inter_D` ([equations_D.py:272-276](equations_D.py#L272-L276)) and its `_F` twin:
     `arg = (K − (1+r^k)·K(-1)) / (K(-1) + chi0)`; `Phi = (chi1/chi2)·(arg²)^(chi2/2)·(K(-1) + chi0)`.
   - Resource cost on the bank's portfolio of physical capital. Enters net worth law of motion directly in `intermediation_P2_D` (not via `rn`).
   - `chi0` is a regularisation constant (avoids div-zero at SS where `arg=0`); `chi1` scales, `chi2=2` is quadratic.
   - SS spec inside `smart_steady_D` ([equations_D.py:80-81](equations_D.py#L80-L81)) currently uses `arg = -rk·K/(K + chi0)`. **Check this matches the off-SS form at `K_t = K_{t-1} = K_ss` exactly** — it should evaluate to `arg = -rk_ss·K_ss/(K_ss + chi0)` from `(K − (1+rk)·K)/(K+chi0)`. Looks consistent but verify with a print at SS.

2. **Macroprudential bond tax `tau_mp = T0 + T1·def_rate` replaces `risk_weight × mp_wedge`.**
   - New block `macro_pru_tax_D` ([equations_D.py:279-285](equations_D.py#L279-L285)): `T = (T0 + T1·def_rate)·(b_D_D + b_F_D)`. Both bond holdings taxed at the same marginal rate.
   - Marginal rate `tau_mp_D` is fed into the bond Euler residuals (`domestic_bond_foc_D`, `portfolio_adj_cost`) as the FOC wedge. Total `T_D` is fed into `intermediation_P2_D` as the realised net-worth drain.
   - Default-sensitivity (`T1·def_rate`) is the macroprudential channel — when the gov is more likely to default, the regulator charges banks more to hold sovereign debt.

3. **`q_b_D` / `q_b_F` are the primary model variables; `rb` derived as interpretation-only.**
   - `bond_yield` ([equations_global.py:33-38](equations_global.py#L33-L38)) defines `rb = 1/q_b − 1`. Comment explicitly says "rb is derived here and NOT used elsewhere."
   - **VERIFY**: grep that `rb_D` / `rb_F` (the derived yields, distinct from `rb_actual_D/F`) are not consumed by any residual block. `rb_actual` is the realised return after default haircut (still used).

4. **Bank releveraging friction `rho_theta` removed.**
   - Old partial-adj `theta_D = ρ·theta_D(-1) + (1-ρ)·theta_tgt_D` is gone. `intermediation_IC_D` ([equations_D.py:216-225](equations_D.py#L216-L225)) now hardwires `ic_res = theta − theta_tgt` (full re-leveraging).
   - **VERIFY**: no `rho_theta` references remain in `equations_*.py` (grep confirmed clean). Notebook may still have stale `rho_theta` in calibration dict — strip if present.

5. **Market-value balance sheets throughout.**
   - `intermediation_IC_D`: `kappa = Q·K/n`, `phi_bD = q_b·b/n`, `phi_bF = q_b·b/n` ([equations_D.py:220-223](equations_D.py#L220-L223)).
   - `bank_return_D`: portfolio shares at t-1 use `q_b(-1)·b(-1)/n(-1)` ([equations_D.py:232-234](equations_D.py#L232-L234)).
   - `portfolio_adj_cost`: `phi_bF_D = q_b_F·b_F_D/n_inter_D` and twin ([equations_global.py:54,59](equations_global.py#L54-L59)).
   - `k_balance_sheet_D`: `Q·K + q_b·b_D + q_b·b_F = θ·n` ([equations_D.py:259-262](equations_D.py#L259-L262)).
   - **Audit**: any leftover face-value `b/n` share that should be `q_b·b/n`? Walk every block that touches portfolio shares.

6. **`excess_return_*_ss` SS anchors subtract `T0`.**
   - Notebook off-SS imports cell — see [model_v12.ipynb:281-284](model_v12.ipynb#L281-L284):
     ```
     excess_return_bD_D_ss = rb_actual_D − rdep_D − T0_D
     excess_return_bF_F_ss = rb_actual_F − rdep_F − T0_F
     excess_return_F_D_ss  = rb_actual_F − rdep_D − T0_D
     excess_return_D_F_ss  = rb_actual_D − rdep_F − T0_F
     ```
   - Reason: at SS, `def_rate = def_rate_ss`, so `tau_mp_ss = T0 + T1·def_rate_ss`. But the bond FOC `(rb − rdep) − excess_ss − ψ·(phi − phi_ss) − tau_mp = 0` at SS reduces to `excess_actual − excess_ss = tau_mp_ss`. **The current code subtracts only `T0`, not `T0 + T1·def_rate_ss`. Check whether `def_rate_ss = 0` (then the two are equal) or non-zero (then the anchor still has a residual `T1·def_rate_ss` mismatch).** Look at default calibration and `government_default_D` at SS.

7. **`phi_bF_D_ss` / `phi_bD_F_ss` use solved SS bond prices, not initial guesses.**
   - [model_v12.ipynb:277-280](model_v12.ipynb#L277-L280):
     ```
     phi_bF_D_ss = ss['q_b_F'] · calibration_start['b_F_D'] / calibration_start['n_inter_D']
     phi_bD_F_ss = ss['q_b_D'] · calibration_start['b_D_F'] / calibration_start['n_inter_F']
     phi_bD_D_ss = ss['q_b_D'] · ss['b_D_D'] / ss['n_inter_D']
     phi_bF_F_ss = ss['q_b_F'] · ss['b_F_F'] / ss['n_inter_F']
     ```
   - **Inconsistency to check**: cross-border shares mix `ss['q_b']` (solved) with `calibration_start['b']` and `calibration_start['n_inter']` (guesses). Domestic shares use `ss[...]` throughout. Is `b_F_D` / `b_D_F` solved by `smart_steady` or held at the guess? If solved, switch to `ss[...]` for consistency.

8. **Diagnostics cell uses market-value `phi` shares and subtracts `tau_mp`.**
   - [model_v12.ipynb:430-440](model_v12.ipynb#L430-L440) computes cross-border and domestic excess-return diagnostics. Verify the formula matches the residual: `(rb − rdep) − excess_ss − ψ·(phi_now − phi_ss) − tau_mp` should print ~0 at SS.

## What to do
1. **Code-review pass.** Walk `equations_D.py`, `equations_F.py`, `equations_global.py`, and every notebook cell that touches calibration/SS/`ha_full`. Confirm `_F` mirrors `_D` exactly. Check:
   - `Phi`/`T` enter `intermediation_P2` not `bank_return`/`rn` (so `Omega/ν/η` SS in `smart_steady` and `intermediation_P1` are consistent).
   - No stale `rho_theta`, `risk_weight`, `mp_wedge` references anywhere.
   - Every share that should be market-value is market-value.
   - The SS-only `bond_price_ss_D` ([equations_D.py:317-323](equations_D.py#L317-L323)) is dropped from `ha_full` (it's a SS-only no-arb anchor; the transition pin is `domestic_bond_foc_D`).
2. **Verify SS solves and residuals print zero.** Run notebook cells through the SS solve + anchors block. Any `*_res` > 1e-8 is a bug.
3. **Verify the `excess_return_*_ss − T0` anchor.** Print `tau_mp_D_ss = T0_D + T1_D·def_rate_D_ss` and compare to `T0_D`. If `def_rate_ss ≠ 0`, the SS anchors need `T0 + T1·def_rate_ss`, not just `T0`.
4. **Run IRFs** on the default shock (`shock_def_D`) and report on/off-SS behaviour. Prior baseline (pre-refactor) had a Y_D overshoot driven by terms-of-trade swing; that issue is *not* this refactor's target but watch for regressions.

## Critical context — do NOT regress these earlier fixes
- `firm_profit_D = (1 − mc_D)·Y_D` (no Rotemberg pricing on main).
- `income_D`: `y_pre = (w·N·e + div)/P_CES`. HH income channel from default flows only via labour + dividend.
- `budget_residual_D` ([equations_D.py:355-362](equations_D.py#L355-L362)): `zeta_writeoff=1` makes default a transfer from bank to gov; setting `zeta=0` makes the gov pay full face value.
- `terms_of_trade` ([equations_global.py:67-69](equations_global.py#L67-L69)) is in the model but **NOT** used as a target on this branch; `p` is pinned by `goods_mkt_D`.
- `add-nkwpc` branch (open as PR #16, awaiting Daniel-Jensen) added a wage NKPC — **do not merge or pull from it** into this refactor unless explicitly asked.

## Files to read first
1. `equations_D.py` ([smart_steady_D L71-101](equations_D.py#L71-L101), [intermediation_IC_D L216-225](equations_D.py#L216-L225), [bank_return_D L228-239](equations_D.py#L228-L239), [intermediation_P1_D L242-255](equations_D.py#L242-L255), [k_balance_sheet_D L258-262](equations_D.py#L258-L262), [cap_adj_cost_inter_D L271-276](equations_D.py#L271-L276), [macro_pru_tax_D L279-285](equations_D.py#L279-L285), [intermediation_P2_D L288-300](equations_D.py#L288-L300), [domestic_bond_foc_D L326-334](equations_D.py#L326-L334)).
2. `equations_F.py` — sanity-check the mirror.
3. `equations_global.py` (loaded above — short file).
4. `model_v12.ipynb`: calibration cell (search `chi0_D`), off-SS imports cell (search `phi_bF_D_ss`), `ha_full` block list, diagnostics cell.

## Style / guardrails
- This is a personal research repo, single user (Adam, PhD). Terse responses preferred. No trailing summaries.
- Don't add error handling, fallbacks, or speculative abstractions. Trust the framework.
- Default to no comments. Existing equation-block comments document non-obvious modelling choices — preserve them.
- Make `_F` edits a strict mirror of `_D` unless the user signals otherwise.
- Do NOT commit unless asked.
