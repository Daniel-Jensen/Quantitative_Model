# Nominal rigidities: price and wage Phillips curves

**Date:** 2026-08-05
**Branch:** `add-nkpc`
**Status:** design approved, not yet implemented

## Goal

Add genuine nominal price and wage rigidities to the two-country monetary-union
HANK model, and make the sticky model the paper's new baseline. The current model
is fully real — `sj.create_model` is named *"Flex Price & Wage, No CB"*
(`code/full_model.py:91`) and there is no inflation variable, nominal rate, or
policy rule anywhere in `code/`.

The rigidities must deliver three things:

1. **Demand-determined output** — a sovereign-risk shock contracts activity
   through demand, not only through the bank/supply side.
2. **Damped labour-side response** — flexible wages let `N` track the real wage
   one-for-one off the labour-supply curve. Sticky wages push adjustment off the
   wage and onto hours. **Note this is *not* a wealth-effect story:** under GHH
   there is no wealth effect on labour supply at all (marginal utility cancels —
   `code/equations_D.py:236`), so the `add-nkwpc` branch's stated motive does not
   apply to this model. See *Relationship to the `add-nkwpc` branch*.
3. **Distributional incidence** — sticky wages shift adjustment from wages onto
   hours, changing who bears the shock across the E4 income quintiles.

## Explicitly out of scope

**The union monetary channel.** No ECB policy rate, no Taylor rule, no Fisher
equation, no change to the deposit-market architecture. **All financial contracts
stay real** — deposits, sovereign bonds, coupons, and the T-2 predetermined
deposit rate are untouched.

This is a deliberate boundary, and it is the seam a referee will probe, so the
paper must state it. The consequence is recorded under *Limitations* below.

Also out of scope: live steady-state markups and the profit-income channel (see
*Steady-state neutrality*). That is a candidate follow-on spec.

## Why the closure works without a policy rate

In a monetary union the nominal exchange rate is fixed at 1, so the terms of
trade **is** the accumulated inflation differential:

```
p = P_F / P_D        =>        p_t / p_{t-1} = (1 + pi_F) / (1 + pi_D)
```

`p` is already an unknown in the model (target `goods_mkt_D`). That identity pins
the inflation *differential* off an object that already exists. One further
normalisation pins the *level*:

```
omega_pi_D * pi_D + (1 - omega_pi_D) * pi_F = 0
```

i.e. the ECB stabilises union-wide producer-price inflation — the `phi_pi -> inf`
limit of a Taylor rule. This is stated as an abstraction, not a modelled rule.
Because financial contracts remain real, no policy rate or Fisher relation is
needed anywhere.

Given the `p` path, these two equations determine `pi_D` and `pi_F`; the two price
Phillips curves then determine `mc_D` and `mc_F`; labour demand determines `w`;
the wage Phillips curves determine `N`; and `goods_mkt_D` determines `p`. No
redundancy.

## Equations

For each country `i` in `{D, F}`.

### Price Phillips curve (new, `equations_{D,F}.py`)

```
nkpc_p_res_i = pi_i - beta_i * pi_i(+1) - kappa_p_i * (mu_p_i * mc_i - 1)
```

Rotemberg form in producer-price inflation. The gap is written as a **ratio**, so
it is unit-free and linearises to exactly `mc_hat`; published Calvo slopes are
therefore directly usable for `kappa_p` with no steady-state rescaling.

### Labour demand (changed, `labor_demand_{D,F}`)

```
w_res_i = w_i - mu_p_i * mc_i * (1 - alpha_i) * Y_i / N_i
```

The `mu_p_i` factor **is** the production subsidy `tau_s = 1 - 1/mu_p`. At
`mc_i = 1/mu_p_i` this collapses to today's competitive condition
`w = (1-alpha) Y / N` identically.

### Wage Phillips curve (new; replaces `labor_market_{D,F}` as a target)

```
pi_w_i     = (w_i / w_i(-1))   * (1 + pi_i)    - 1
pi_w_i_p1  = (w_i(+1) / w_i)   * (1 + pi_i(+1)) - 1
gap_w_i    = vphi_i * N_i ** (1 / frisch_i) / (w_i / P_CES_i) - 1
nkpc_w_res_i = pi_w_i - beta_i * pi_w_i_p1 - kappa_w_i * gap_w_i
```

Nominal wage `W = w * P_i`, hence `1 + pi_w = (w/w(-1)) * (1 + pi)`.

Under GHH the MRS is `vphi * N^(1/frisch)` with marginal utility cancelling —
`labor_ss_D` states this explicitly (`code/equations_D.py:236`). So `gap_w_i` is
the existing `labor_mkt_res_i` in ratio form, and heterogeneity requires no
marginal-utility-weighted aggregation.

**Implementation constraint:** SSJ `@simple` blocks apply lead/lag operators only
to *block inputs*, never to locally computed variables. `pi_w_i_p1` must therefore
be built explicitly from `w_i(+1)` and `pi_i(+1)` as shown — `pi_w_i(+1)` will not
work.

**Disposition of `labor_market_{D,F}`:** the function stays in
`equations_{D,F}.py` (it documents the flexible-wage condition, of which
`gap_w_i` is the ratio form) but is **removed from all three dynamic
`create_model` lists**, replaced by `wage_nkpc_{D,F}`. It is not in the
steady-state block list (`code/steady_state.py:151`), so the SS solve is
unaffected.

### Global (new, `equations_global.py`)

```
tot_res      = p / p(-1) - (1 + pi_F) / (1 + pi_D)
union_pi_res = omega_pi_D * pi_D + (1 - omega_pi_D) * pi_F
```

### Discounting

The curves discount at the constant `beta_i` rather than at `SDF_i`. Because
`pi_ss = 0`, the SDF deviation multiplies a zero and the two are **identical to
first order** — and the model is solved by linearised `solve_jacobian`. The choice
is immaterial; `beta_i` is the textbook form.

## System: 23x23 -> 27x27

New unknowns: `mc_D, pi_D, mc_F, pi_F`.
New targets: `nkpc_p_res_D, nkpc_p_res_F, tot_res, union_pi_res`.
Renamed target: `labor_mkt_res_{D,F}` -> `nkpc_w_res_{D,F}`.

```python
unknowns_tp = [
    'K_D','n_inter_D','div_D','I_D','Q_D','b_gov_D','N_D','b_F_D','w_D','rdep_D','mc_D','pi_D',
    'K_F','n_inter_F','div_F','I_F','Q_F','b_gov_F','N_F','b_D_F','w_F','rdep_F','mc_F','pi_F',
    'p','q_b_D','q_b_F',
]
targets_tp = [
    'deposit_mkt_D','K_res_D','n_inter_val_D','div_res_D','capital_res_D','q_res_D',
    'b_gov_res_D','b_F_D_res','nkpc_w_res_D','w_res_D','nkpc_p_res_D',
    'deposit_mkt_F','K_res_F','n_inter_val_F','div_res_F','capital_res_F','q_res_F',
    'b_gov_res_F','b_D_F_res','nkpc_w_res_F','w_res_F','nkpc_p_res_F',
    'goods_mkt_D','rb_D_res','rb_F_res','tot_res','union_pi_res',
]
```

## Steady-state neutrality

Markups are neutralised by subsidy, so **the steady state is bit-identical to
today's**. With `pi_i_ss = 0` and `mc_i_ss = 1/mu_p_i`, all five new residuals are
*exactly* zero at the existing steady state:

| Residual | At SS |
|---|---|
| `nkpc_p_res_i` | `0 - 0 - kappa_p*(mu_p*(1/mu_p) - 1) = 0` |
| `w_res_i` | `w - 1*(1-alpha)Y/N = 0` (today's condition) |
| `pi_w_i` | `(1)(1) - 1 = 0` |
| `nkpc_w_res_i` | `0 - 0 - kappa_w*0 = 0` (gap zero by `labor_ss_i`) |
| `tot_res` | `1 - 1 = 0` |
| `union_pi_res` | `0` |

`K`, `rk`, `w`, `N`, the spread, `EL_price_D`, the IC-delta consistency check,
`assert_gk_well_posed`, and every Walras residual are therefore unchanged.
`mu_w` stays at **1.0** — it is the steady-state neutralising device, while a
finite `epsilon_w = 11` governs the wage-curve slope. `labor_ss_{D,F}`'s `vphi`
calibration is untouched.

`steady_state.py` needs only `mc_{D,F} = 1/mu_p_{D,F}` and `pi_{D,F} = 0.0` seeded
into `calibration_start` so `ss_final` carries them into the dynamic solve.

## The flex model is the exact `kappa -> inf` limit

Dividing each Phillips-curve residual by `kappa` and letting `kappa -> inf` gives
back `gap_w_i = 0` (today's `labor_mkt_res_i`) and `mu_p*mc = 1` (today's
`w_res_i`) — identically, not approximately.

Two consequences:

- **No `STICKY` switch is needed.** The flex robustness run is a calibration
  override `kappa_p = kappa_w = 1e4` through
  `experiments/common.calibration_override`. No branching inside the equations
  and no second code path to drift — the failure mode CLAUDE.md records for the
  retired `audit_artifacts/` harness.
- **It supplies the regression gate** (rollout step 1). Very large `kappa` may be
  ill-conditioned; if `1e4` fails to converge, step down and record the largest
  value that does.

## The het block is unchanged

Hours become demand-determined and are rationed in proportion to productivity
`e` — already exactly how `income_D` distributes labour income
(`code/equations_D.py:66`) — and the GHH disutility already takes aggregate `N_D`.
That is the standard union-wage-setting assumption. Nothing in the EGM changes.

## Rotemberg resource costs are omitted

`(phi/2) * pi^2 * Y` is quadratic around `pi_ss = 0`, so its derivative there is
zero and `solve_jacobian` linearises — it **cannot** move any IRF. Including it
would inject a nonlinear term into `goods_mkt_D`, which currently holds at 1e-14.
It would matter only for a second-order welfare exercise; E1's welfare is
first-order off `U_D = X_D / C_D_ss`.

A useful side effect: the `market_clearing_D` identity is untouched, so E2's
self-verifying dY decomposition must still close at 1e-7 — an independent check
on the whole change.

## Calibration

New parameters in `code/calibration.py`. The slope formulae are evaluated at
`beta = 0.985`; `beta_D` and `beta_F` are separately solved SS unknowns near that
value, and the slopes are fixed calibration constants, not functions of the
solved betas.

| Parameter | Value | Basis |
|---|---|---|
| `mu_p_D/F` | 1.20 | `epsilon_p = 6`, standard. **Free to first order** — see below |
| `kappa_p_D/F` | 0.0871 | Calvo `theta_p = 0.75`, `(1-theta)(1-beta*theta)/theta`. Euro-area IPN median price duration ~4 quarters (Alvarez et al. 2006; Dhyne et al. 2006) |
| `kappa_w_D/F` | 0.0038 | `theta_w = 0.75`, `epsilon_w = 11`, `phi = 1/frisch = 2`: `(1-theta)(1-beta*theta)/(theta*(1+epsilon_w*phi))` |
| `omega_pi_D` | 0.071 | `1 - kappa_cb_F`, the already-documented renormalised capital key (BuBa 26.1 / BoG 2.0 of the euro-area key) |

### `mu_p` does not matter to first order

Under subsidy-neutralisation the gap `mu_p*mc - 1` linearises to exactly `mc_hat`
irrespective of `mu_p` (since `mu_p * mc_ss = 1`), and labour demand's steady
state is `w = (1-alpha)Y/N` either way. So `mu_p` is a free normalisation in this
build and needs no defending. It becomes a genuine calibration choice only if the
live-markup / profit-income follow-on is taken up.

### Do not use model GDP weights for `omega_pi_D`

The model normalises `Y_D_ss ~ Y_F_ss ~ 1`, so GDP weights would give
`omega_pi_D ~ 0.5`, implying a Greek deflation forces German inflation up nearly
one-for-one — the opposite of 2010-12. The capital key (~0.071, and ~0.076 on
actual 2011 GR/DE nominal GDP) is the correct weight and is already in the
codebase with a citation.

Working the closure through makes the stake concrete. Combining the two global
equations gives `pi_D = -(1 - omega_pi_D) * dlog p` and
`pi_F = omega_pi_D * dlog p`. At `omega_pi_D = 0.071`, **93% of any
terms-of-trade adjustment appears as D (Greek) producer-price deflation and 7% as
F (German) inflation** — the 2010-12 internal-devaluation pattern. At 0.5 it
splits evenly, which is counterfactual for GR/DE.

Bi, Foerster and Traum (2026) use 0.5/0.5, but for Italy/Germany as two
comparably-sized blocs, and inside a *Taylor rule* where the weight sets a policy
response rather than allocating a given differential. Report the choice and note
the 0.5 alternative; do not adopt it.

## Rollout

**Step 0 — refactor first, as a pure no-op.** `full_model.py:69`, `tpi.py:145`
and `diagnostics/regimes/regime_model.py:160` each hardcode the `create_model`
block list. Extract one `build_block_list()` in `full_model.py` and point the
other two at it. **Verify `code/main.py` output is bit-identical before adding
anything.** Doing this after would mean validating two changes at once, and
CLAUDE.md records that a drifting duplicate model is what invalidated
`audit_artifacts/`.

**Step 1 — structural-equivalence gate.** Add the six blocks, wire the 27x27
system, run at `kappa = 1e4`. Must reproduce the current baseline IRFs to solver
tolerance and hold every threshold: `goods_mkt_D <= 1e-14`, `goods_mkt_F` and
`ca_res_D <= 1e-7`, `deposit_mkt_D/F <= 1e-13`. If this fails the wiring is wrong
and nothing downstream is worth debugging.

**Step 2 — dial `kappa` to calibrated values.** Check residuals, doom-loop signs
(`n_inter_D[0] < 0`, `Y_D[0] < 0`), stability (`b_gov_D[499] ~ 0`), the IC-delta
consistency check, and `assert_gk_well_posed`.

**Step 3 — re-tune `psi_lambda_B`** to the 150bp-per-1pp-default-shock target. It
will move, because spread transmission now runs through a sticky terms of trade.

**Step 4 — regenerate.** `diagnostics/regimes/regime_model.py --force`, then
`experiments/run_all.py`, then figures. E1, E2, E3, E4 and the declining-loading
key figure are all rebuilt on the sticky model; flex becomes an appendix
comparison.

**Step 5 — docs.** STATE.md, PROGRESS.md, HANDOFF.md (hook-enforced), plus
SPEC.md and CLAUDE.md.

## Verification

`code/main.py` remains the structural regression test. On top of it:

- **Equivalence:** at `kappa_p = kappa_w = 1e4`, IRFs match the pre-change
  baseline to solver tolerance.
- **Steady state:** every solved SS object is bit-identical to `main`'s.
- **Residual thresholds:** unchanged, as listed in step 1.
- **Sign checks:** `n_inter_D[0]` and `Y_D[0]` both negative on the default shock.
- **E2 closure:** dY decomposition still asserts at 1e-7.
- **Slope sweep:** `kappa_w` in {0.002, 0.0038, 0.01, 0.03} and `kappa_p` in
  {0.03, 0.087, 0.2}; report the stable region. This is a required robustness
  table regardless.

## Risks

**`kappa_w = 0.0038` is flat — wages are quite rigid.** Against a model with
documented near-unit-root and explosive regions (F-1, GK-2), that may
destabilise. Mitigated by the slope sweep above.

**The `psi_lambda_B` re-tune may land in a breakdown region.** CLAUDE.md puts the
documented breakdown around 4-5 at `n_inter = 3.0`. Step 3 must re-verify
stability at whatever value it lands on, not merely hit the moment.

**Solve time grows.** 23 -> 27 unknowns at T=500; expect ~3 min to become 4-5.
E3 does two re-solves, so `experiments/run_all.py --skip-e3` matters more during
iteration.

## Limitations to state in the paper

**The demand channel is relative-price, not intertemporal.** Excluding the union
channel means there is no real-rate route from stickiness to spending. Output is
demand-determined here through the **sluggish terms of trade** — a Greek
contraction cannot depreciate `p` quickly, so net exports do not cushion it — plus
the markup wedge in labour demand. That is a defensible mechanism and *is* the
internal-devaluation story for 2010-12, but it is not the textbook monetary
channel, and the paper must say which one it claims.

**Steady-state markups are subsidised away**, so the profit-income channel is
absent. Live markups (`mu_p ~ 1.2`) would put pure profits at ~17% of output as a
new, highly concentrated household income stream — which would change E4 incidence
on its own. Deferred to a follow-on spec.

## Benchmark: Bi, Foerster and Traum (2026)

*"Asset Purchases in a Monetary Union With Default and Liquidity Risks", FRBSF
Working Paper 2025-10, https://doi.org/10.24148/wp2025-10.* The closest published
analogue to this model: two-country monetary union, Gertler-Karadi
intermediaries, endogenous sovereign default, cross-border sovereign holdings,
targeted ECB asset purchases — calibrated to Italy/Germany 2012 rather than
Greece/Germany 2010-12.

| | Bi-Foerster-Traum | This spec |
|---|---|---|
| Price rigidity | Rotemberg, exact nonlinear (their eq. 2.14) | Rotemberg, linear-equivalent |
| Wage rigidity | **none** — flexible, `chi*L^sigma_l = U_c*w` (A.9) | wage NKPC, both countries |
| SS markup | **live**, `theta^c = 11` -> `mc_ss = 10/11` | subsidy-neutralised, SS bit-identical |
| Firm profits | routed to representative household | absent (neutralised) |
| Nominal anchor | Taylor rule, `phi_pi=1.6, phi_y=0.07, phi_r=0.85` | union-inflation normalisation, no rate |
| Union weights | 0.5 / 0.5 | `omega_pi_D = 0.071` |
| ToT identity | `rer_t/rer_{t-1} = pi*_t/pi_t` (A.80), CPI form | same relation, PPI form |
| Financial contracts | nominal (deposits, debt, net worth all deflated) | real |
| Solution | 2nd-order perturbation (endogenous regime switching) | 1st-order SSJ |

**The price slope agrees.** Their `psi` maps to a Calvo-equivalent slope
`(1-xi)(1-beta*xi)/xi = 0.0846` at `xi_p = 0.75, beta = 0.995`. This spec's
`kappa_p = 0.0871` at `beta = 0.985` is the same number to within 3%; the gap is
entirely the discount factor. No change required.

**The wage slope has no benchmark in this model class.** They have no wage curve.
Note also that `kappa_w`'s flatness is driven substantially by the *preference*
parameter, not by any wage-rigidity estimate: the `(1 + epsilon_w/frisch)`
denominator is 23 at `frisch = 0.5` versus 12 at their `sigma_l = 1`. This is why
the slope sweep is load-bearing rather than cosmetic.

**They never report a flexible-price counterfactual.** Their Table 1 decomposes
over the liquidity-risk channel, the fiscal-limit shift, and the debt change —
never over price stickiness. The `kappa -> inf` gate in this spec produces exactly
that counterfactual as a by-product, so it is reportable output, not merely a
regression test.

**Their nominal side does little propagation work, and that is a live risk here.**
Inflation moves +/-0.1% while investment moves 9% and output 0.6% (their Figure
3); Tables 2 and 3 show inflation entries of 0.00-0.02 against investment at 0.53.
Their Section 4.1 mechanism narrative treats inflation as an *outcome* of the
relative-price move, never as a channel. What drives their output contraction is
the **loan-in-advance constraint** (`eta^I = 0.65/0.75` of investment must be
debt-financed, their eq. 2.10) — a real financial friction this model does not
have. Against their `-0.6%` output impact, this model's `Y_D[0] = -0.0149%` is two
orders of magnitude smaller, and a Phillips curve is unlikely to close that gap.
A Sims-Wu working-capital constraint is the natural candidate if the sticky build
leaves `Y_D[0]` implausibly small; it is **out of scope here** and would need its
own design pass.

## Relationship to the `add-nkwpc` branch

`add-nkwpc` (commit `2377f79`, off `08e1010`, pre-reorganisation) is a single
26-line commit adding `wage_setting_{D,F}` only. It was never wired into the model
list and has no `kappa_w` calibration. It is a **template, not code to port**:

- It is a **real-wage** Rotemberg curve (`pi_w = w/w(-1) - 1` on the *real* wage,
  explicitly "no CB needed") — a real adjustment friction, not a nominal rigidity.
  There is no price-side counterpart at all.
- It divides the MRS by `UCE_D`, which is correct under separable preferences but
  **wrong under the GHH preferences this model uses** — cf. `labor_ss_D` and
  `labor_market_D`, neither of which contains `UCE`.
- It uses `w_D` where the existing labour condition uses `w_D / P_CES_D`, dropping
  the CES bundle deflator.
- **Its stated motivation is wrong for this model.** The comment says sticky wages
  stop "the household wealth effect" translating into `N`. GHH preferences have no
  wealth effect on labour supply by construction — that is the point of GHH, and
  both `labor_market_D` and `labor_ss_D` show marginal utility cancelling out of
  the intratemporal condition. What sticky wages actually buy here is that
  adjustment moves off the real wage and onto hours, which is what the
  distributional goal needs (hours are rationed proportional to `e`), but the
  channel must be described correctly in the paper.
