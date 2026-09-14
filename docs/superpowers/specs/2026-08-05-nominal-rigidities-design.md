# Nominal rigidities: price Phillips curves and nominal deposits

**Date:** 2026-08-05
**Branch:** `add-nkpc`
**Status:** design approved, not yet implemented

## Goal

Add nominal price rigidity and nominal deposit contracts to the two-country
monetary-union HANK model, and make the sticky model the paper's new baseline.
The current model is fully real — `sj.create_model` is named *"Flex Price & Wage,
No CB"* (`code/full_model.py:91`) and there is no inflation variable, nominal
rate, or policy rule anywhere in `code/`.

1. **Demand-determined output.** Today `N` is purely supply-determined: combining
   the two flexible labour conditions gives `(1-alpha)Y/(N*P_CES) =
   vphi*N^(1/frisch)`, a function of `Z`, `K` and `P_CES` alone. A markup wedge in
   labour demand breaks that.
2. **Fisher debt deflation on bank balance sheets.** Nominal deposits against real
   assets make banks net nominal debtors, so a deflation at impact raises the real
   value of their liabilities and deepens the net-worth loss. This targets the
   fact that `Y_D[0] = -0.0149%` is two orders of magnitude below the benchmark's
   `-0.6%`.
3. **Distributional incidence** — the changed paths of output, wages, the terms of
   trade and the real deposit return alter how the shock lands across the E4
   income quintiles.

### The `C_D[0]` motivation, stated honestly

At the live calibration `C_D[0] = +0.2164%` against `Y_D[0] = -0.0149%` and
`I_D[0] = -0.7718%` (`docs/STATE.md`, E1 passive). `docs/HANDOFF.md` already flags
this: the model's crisis is an investment bust, not a consumption bust, which is
counterfactual for Greece 2010-13.

**This spec is not guaranteed to flip that sign, and must not be sold as if it
were.** Two forces pull opposite ways once deposits are nominal:

- *Intertemporal substitution* — expected deflation raises the ex-ante real
  deposit rate; at `eis = 0.5` households substitute out of current consumption.
  Pushes `C_D[0]` **down**.
- *Impact revaluation* — households are net nominal creditors, so the same
  deflation is a windfall on the existing deposit stock. Pushes `C_D[0]` **up**.

Which dominates is quantitative. What is certain is that the channel becomes
connected at all: the model currently generates `pi_D = -0.93 * dlog p` and lets
it affect nothing. For context, Bi, Foerster and Traum (2026) get consumption
rising on impact too — with a Taylor rule, nominal debt *and* a loan-in-advance
constraint (their §4.1: "Consumption rises briefly on impact, but quickly turns
negative"). The impact sign is a property of this model class.

## Explicitly out of scope

**Wage rigidity (author decision, 2026-08-05).** Price stickiness only; the labour
market keeps today's flexible GHH condition `labor_market_{D,F}` unchanged. This
matches Bi-Foerster-Traum. Consequence: no hours-rationing channel, so
distributional incidence works only through the output / wage / terms-of-trade /
deposit-return path.

**A Taylor rule (author decision, 2026-08-05).** The blocker is structural, not
effort: a policy rule only has traction if it pins a real rate, and pinning
`rdep_D`/`rdep_F` frees *both* deposit-market conditions, which then need an
absorber — either a zero-net-supply cross-border banking claim (rewriting
`external_account_D`, currently at 1e-7) or an ECB reserve asset inside the GK
incentive constraint. A deposit spread over the policy rate does not work, since
a free spread absorbs the rate one-for-one and the rule does nothing.

**Note this is *not* what makes deposits nominal expensive.** The absorber
problem belongs to the policy rule alone. With no rule, the deposit rate stays a
free unknown and simply becomes nominal — see below.

**Nominal sovereign bonds (author decision, 2026-08-05).** Bonds, coupons and the
external account stay real. This is a **deliberate asymmetry** that maximises
banks' Fisher exposure, and the paper must label it as such: euro-area sovereign
debt is nominal, Bi-Foerster-Traum deflate debt and net worth alike, and a referee
will ask. Making bonds nominal would touch `bond_return`, `bond_price_ss`,
`budget_residual`, `government_ss`, both bond FOCs, and `external_account_D`
(1e-7, and it carries the W-2 `p`-conversion). Candidate follow-on spec.

Also out of scope: live steady-state markups, and a Sims-Wu loan-in-advance
constraint (see *Benchmark*).

## Why the closure works without a policy rate

In a monetary union the nominal exchange rate is fixed at 1, so the terms of trade
**is** the accumulated inflation differential:

```
p = P_F / P_D        =>        p_t / p_{t-1} = (1 + pi_F) / (1 + pi_D)
```

`p` is already an unknown (target `goods_mkt_D`), so the identity pins the
inflation *differential* off an existing object. One normalisation pins the level:

```
omega_pi_D * pi_D + (1 - omega_pi_D) * pi_F = 0
```

i.e. the ECB stabilises union-wide producer-price inflation — the `phi_pi -> inf`
limit of a Taylor rule, stated as an abstraction rather than a modelled rule.
Given the `p` path these two determine `pi_D` and `pi_F`; the Phillips curves
determine `mc_D` and `mc_F`; labour demand determines `w`; the unchanged labour-
supply condition determines `N`; the deposit markets determine `i_dep_{D,F}`; and
`goods_mkt_D` determines `p`. No redundancy.

## Equations

For each country `i` in `{D, F}`.

### Price Phillips curve (new, `equations_{D,F}.py`)

```
nkpc_p_res_i = pi_i - beta_i * pi_i(+1) - kappa_p_i * (mu_p_i * mc_i - 1)
```

Rotemberg form in producer-price inflation. The gap is a **ratio**, so it is
unit-free and linearises to exactly `mc_hat`; published Calvo slopes are directly
usable for `kappa_p` with no steady-state rescaling.

### Labour demand (changed, `labor_demand_{D,F}`)

```
w_res_i = w_i - mu_p_i * mc_i * (1 - alpha_i) * Y_i / N_i
```

The `mu_p_i` factor **is** the production subsidy `tau_s = 1 - 1/mu_p`. At
`mc_i = 1/mu_p_i` this collapses to today's competitive condition identically.

### Markup rent (new, `firm_profit_{D,F}`) — required, not optional

With the markup in labour demand only, factor payments stop exhausting output:
`w*N = mu_p*mc*(1-alpha)Y` while the capital return is untouched
(`capital_adj_D` keeps `mpk = alpha*Z*K^(alpha-1)*N^(1-alpha)`). Off steady state
`mu_p*mc != 1` leaves an unrouted rent, which is a Walras leak of the W-1 / W-2
class. It **must** be routed.

```
profit_i = (1 - mu_p_i * mc_i) * (1 - alpha_i) * Y_i
```

Distributed to households **in proportion to productivity `e`** (Auclert-Rognlie-
Straub), so `income_i`'s numerator gains `profit_i * e_grid_i`. Household labour
plus profit income is then `w*N*e + profit*e = (1-alpha)*Y*e` — *identical to the
flexible model* — and factor payments exhaust output exactly.

Three properties make this the right rule rather than the textbook lump-sum one:

- Markups are **countercyclical** (`mc` falls in a downturn), so a lump-sum profit
  rebate would hand households rising income exactly when output falls, pushing
  `C_D[0]` further up — the opposite of the motivation.
- The wedge then affects the **firm's hiring decision only**, which is the channel
  goal 1 wants, and leaves household income unchanged.
- The distribution is proportional to the household's *type* `e`, not to hours, so
  the marginal wage is still `w_i` and **`labor_market_{D,F}` is unchanged**.
  `labor_ss_{D,F}`'s `vphi` calibration is untouched (`profit = 0` at SS).

### Nominal deposits (changed, `deposit_return_{D,F}` + substitutions)

`rdep_{D,F}` is reinterpreted as a **nominal** rate `i_dep_{D,F}`. It remains the
free unknown clearing `deposit_mkt_{D,F}` — no absorber, no cross-border claim,
`external_account_D` untouched, count unchanged.

```
Rgross_i = (1 + i_dep_i(-1)) * P_CES_i(-1) / P_CES_i / (1 + pi_i)
```

`deposit_return_i` already has exactly this shape — a predetermined rate times a
period-t deflator — so this is a one-line change. **T-2 is not reopened:** the
rate stays locked at t-1; only the deflator becomes period-t, which the block
already does. T-2 was about using a period-t *unknown rate*.

Two derived real rates (`@simple` outputs, **not** unknowns):

```
rdep_expost_i = (1 + i_dep_i(-1)) / (1 + pi_i)     - 1    # realised at t on t-1 deposits
rdep_exante_i = (1 + i_dep_i)     / (1 + pi_i(+1)) - 1    # locked at t, for t -> t+1
```

Substitutions, following the existing timing convention exactly:

| Block | Today | Becomes |
|---|---|---|
| `bank_return_i` | `rdep_i(-1)` | `rdep_expost_i` |
| `capital_fund_i` | `rdep_i(-1)` | `rdep_expost_i` |
| `intermediation_P1_i` | `rdep_i` | `rdep_exante_i` |
| `divert_bond_foc_i` | `rdep_i` | `rdep_exante_i` |
| `divert_portfolio_adj` | `rdep_D`, `rdep_F` | `rdep_exante_{D,F}` |
| `steady_auxilliary_i`, `smart_steady_i` | `rdep_i` | `i_dep_i` (SS only, `pi = 0`) |

At the steady state `pi = 0` and all three rates collapse to today's `rdep`.

Because banks hold **real** assets against **nominal** liabilities, they are net
nominal debtors: deflation raises the real value of deposits and deepens the
net-worth loss. That is the Fisher-Bernanke channel goal 2 is buying.

### Labour supply — unchanged

`labor_market_{D,F}` (`code/equations_D.py:287`) stays exactly as it is, remains
in all three dynamic block lists, and keeps `labor_mkt_res_{D,F}` as the target
for `N_{D,F}`.

### Global (new, `equations_global.py`)

```
tot_res      = p / p(-1) - (1 + pi_F) / (1 + pi_D)
union_pi_res = omega_pi_D * pi_D + (1 - omega_pi_D) * pi_F
```

### Discounting

The Phillips curve discounts at constant `beta_i` rather than `SDF_i`. Because
`pi_ss = 0` the SDF deviation multiplies a zero, so the two are **identical to
first order** and the model is solved by linearised `solve_jacobian`. Immaterial;
`beta_i` is the textbook form.

## System: 23x23 -> 27x27

New unknowns: `mc_D, pi_D, mc_F, pi_F`. New targets: `nkpc_p_res_D,
nkpc_p_res_F, tot_res, union_pi_res`. `rdep_{D,F}` is renamed `i_dep_{D,F}`; no
targets are renamed or removed.

```python
unknowns_tp = [
    'K_D','n_inter_D','div_D','I_D','Q_D','b_gov_D','N_D','b_F_D','w_D','i_dep_D','mc_D','pi_D',
    'K_F','n_inter_F','div_F','I_F','Q_F','b_gov_F','N_F','b_D_F','w_F','i_dep_F','mc_F','pi_F',
    'p','q_b_D','q_b_F',
]
targets_tp = [
    'deposit_mkt_D','K_res_D','n_inter_val_D','div_res_D','capital_res_D','q_res_D',
    'b_gov_res_D','b_F_D_res','labor_mkt_res_D','w_res_D','nkpc_p_res_D',
    'deposit_mkt_F','K_res_F','n_inter_val_F','div_res_F','capital_res_F','q_res_F',
    'b_gov_res_F','b_D_F_res','labor_mkt_res_F','w_res_F','nkpc_p_res_F',
    'goods_mkt_D','rb_D_res','rb_F_res','tot_res','union_pi_res',
]
```

## Steady-state neutrality

Markups are subsidy-neutralised and `pi_ss = 0`, so **the steady state is
bit-identical to today's**. All new residuals are *exactly* zero at it:

| Residual | At SS |
|---|---|
| `nkpc_p_res_i` | `0 - 0 - kappa_p*(mu_p*(1/mu_p) - 1) = 0` |
| `w_res_i` | `w - 1*(1-alpha)Y/N = 0` (today's condition) |
| `profit_i` | `(1 - 1)*(1-alpha)Y = 0` |
| `tot_res` | `1 - 1 = 0` |
| `union_pi_res` | `0` |
| `rdep_expost_i`, `rdep_exante_i` | both `= i_dep_i = ` today's `rdep_i` |

`labor_mkt_res_{D,F}` is unchanged and already zero. `K`, `rk`, `w`, `N`, the
spread, `EL_price_D`, the IC-delta check, `assert_gk_well_posed` and every Walras
residual are unchanged. `steady_state.py` needs only `mc_{D,F} = 1/mu_p_{D,F}` and
`pi_{D,F} = 0.0` seeded into `calibration_start`.

## The flex model is the exact `kappa_p -> inf` limit

Dividing the Phillips-curve residual by `kappa_p` and letting `kappa_p -> inf`
gives `mu_p*mc = 1`, hence today's `w_res_i` and `profit_i = 0` — identically, not
approximately. With `profit = 0`, `income_i` reverts exactly, and with `pi` driven
only by the (then flexible) terms of trade the nominal deposit terms still bite,
so the limit recovers **flex prices with nominal deposits**, not the current
`main`. See the rollout for how the gate is staged around this.

No `STICKY` switch is needed: the flexible-price run is a calibration override
`kappa_p = 1e4` through `experiments/common.calibration_override`. No branching
inside the equations and no second code path to drift — the failure mode CLAUDE.md
records for the retired `audit_artifacts/` harness.

## Rotemberg resource costs are omitted

`(phi/2)*pi^2*Y` is quadratic around `pi_ss = 0`, so its derivative there is zero
and `solve_jacobian` linearises — it **cannot** move any IRF. Including it would
inject a nonlinear term into `goods_mkt_D`, which holds at 1e-14. It would matter
only for a second-order welfare exercise; E1's welfare is first-order off
`U_D = X_D / C_D_ss`. Side effect: `market_clearing_D` is untouched, so E2's
self-verifying dY decomposition must still close at 1e-7 — an independent check on
the whole change.

## Calibration

New parameters in `code/calibration.py`. The slope is evaluated at `beta = 0.985`;
`beta_D` and `beta_F` are separately solved SS unknowns near that value, and the
slope is a fixed constant, not a function of the solved betas.

| Parameter | Value | Basis |
|---|---|---|
| `mu_p_D/F` | 1.20 | `epsilon_p = 6`, standard. **Free to first order** — see below |
| `kappa_p_D/F` | 0.0871 | Calvo `theta_p = 0.75`, `(1-theta)(1-beta*theta)/theta`. Euro-area IPN median price duration ~4 quarters (Alvarez et al. 2006; Dhyne et al. 2006) |
| `omega_pi_D` | 0.071 | `1 - kappa_cb_F`, the documented renormalised capital key (BuBa 26.1 / BoG 2.0) |

### `mu_p` does not matter to first order

Under subsidy-neutralisation the gap `mu_p*mc - 1` linearises to exactly `mc_hat`
irrespective of `mu_p` (since `mu_p*mc_ss = 1`), and labour demand's steady state
is `w = (1-alpha)Y/N` either way. `mu_p` is a free normalisation here and needs no
defending. It becomes a genuine choice only if the live-markup follow-on is taken.

### Do not use model GDP weights for `omega_pi_D`

The model normalises `Y_D_ss ~ Y_F_ss ~ 1`, so GDP weights would give
`omega_pi_D ~ 0.5`, implying Greek deflation forces German inflation up nearly
one-for-one — the opposite of 2010-12. Combining the two global equations gives
`pi_D = -(1 - omega_pi_D)*dlog p` and `pi_F = omega_pi_D*dlog p`, so at
`omega_pi_D = 0.071` **93% of any terms-of-trade adjustment appears as D
producer-price deflation and 7% as F inflation** — the internal-devaluation
pattern. At 0.5 it splits evenly, which is counterfactual for GR/DE.

Bi-Foerster-Traum use 0.5/0.5, but for Italy/Germany as comparably-sized blocs and
inside a *Taylor rule*, where the weight sets a policy response rather than
allocating a given differential. Report the choice; do not adopt it.

**This parameter is now load-bearing twice over.** It scales `pi_D`, which drives
the Fisher revaluation on bank balance sheets. Include it in the sweep.

## Rollout

**Step 0 — refactor first, as a pure no-op.** `full_model.py:69`, `tpi.py:145` and
`diagnostics/regimes/regime_model.py:160` each hardcode the `create_model` block
list. Extract one `build_block_list()` in `full_model.py` and point the other two
at it. **Verify `code/main.py` output is bit-identical before adding anything.**
CLAUDE.md records that a drifting duplicate model is what invalidated
`audit_artifacts/`.

**Step 1 — sticky prices only, real deposits, equivalence gate.** Add the price
NKPCs, `firm_profit`, the two global blocks; wire 27x27; run at `kappa_p = 1e4`
with the deposit blocks still real. Must reproduce the current baseline IRFs to
solver tolerance and hold every threshold: `goods_mkt_D <= 1e-14`, `goods_mkt_F`
and `ca_res_D <= 1e-7`, `deposit_mkt_D/F <= 1e-13`. If this fails the wiring is
wrong and nothing downstream is worth debugging. Very large `kappa_p` may be
ill-conditioned; if `1e4` fails to converge, step down and record the largest
value that does.

**Step 2 — dial `kappa_p` to 0.0871**, deposits still real. Check residuals,
doom-loop signs (`n_inter_D[0] < 0`, `Y_D[0] < 0`), stability
(`b_gov_D[499] ~ 0`), the IC-delta check, `assert_gk_well_posed`. **Record
`Y_D[0]`, `C_D[0]`, `I_D[0]`, `n_inter_D[0]` here** — this is the clean measure of
what price stickiness alone does, and it is a reportable result.

**Step 3 — switch deposits to nominal.** One change at a time: this is where the
Fisher channel arrives, and it must be attributable. Re-run the same checks.
`n_inter_D[0]` must fall by *more* than at step 2; if it does not, the sign is
wrong somewhere in the ex-post/ex-ante substitution table.

**Step 4 — re-tune `psi_lambda_B`** to the 150bp-per-1pp-default-shock target.
Expect a large move — see *Risks*.

**Step 5 — regenerate.** `diagnostics/regimes/regime_model.py --force`, then
`experiments/run_all.py`, then figures. E1-E4 and the declining-loading key figure
all rebuilt; flex-price becomes an appendix comparison.

**Step 6 — docs.** STATE.md, PROGRESS.md, HANDOFF.md (hook-enforced), plus
SPEC.md and CLAUDE.md.

## Verification

`code/main.py` remains the structural regression test. On top of it:

- **Equivalence (step 1):** at `kappa_p = 1e4` with real deposits, IRFs match the
  pre-change baseline to solver tolerance.
- **Steady state:** every solved SS object bit-identical to `main`'s, at every step.
- **Residual thresholds:** unchanged throughout, as listed in step 1.
- **Sign checks:** `n_inter_D[0]` and `Y_D[0]` both negative on the default shock.
- **Fisher sign (step 3):** `n_inter_D[0]` strictly more negative than at step 2.
- **E2 closure:** dY decomposition still asserts at 1e-7.
- **Sweep:** `kappa_p` in {0.03, 0.087, 0.2} and `omega_pi_D` in {0.071, 0.2, 0.5};
  report the stable region. Required robustness table regardless.
- **Report `C_D[0]` explicitly** at steps 2 and 3, whatever its sign.

## Risks

**The Fisher channel may dominate rather than supplement the doom loop.**
`D_supply ~ (theta-1)*n_inter ~ 9` against `n_inter ~ 3`, so a price-level surprise
is levered roughly 3x onto net worth. If `pi_D[0]` lands near -0.1%, that is on the
order of a 9% net-worth hit against today's -3.38%. Given F-1 and GK-2, treat
instability as a live possibility and be ready to report the `omega_pi_D` sweep as
the containing parameter.

**`psi_lambda_B` will move a long way**, since spread transmission now runs through
both a sticky terms of trade and a Fisher revaluation. CLAUDE.md puts the
documented breakdown around 4-5 at `n_inter = 3.0`; step 4 must re-verify stability
at whatever value it lands on, not merely hit the moment.

**`C_D[0]` may not change sign, or may rise further.** See *The `C_D[0]`
motivation*. Step 2 and step 3 both report it; do not commit to a claim about it in
the draft before those numbers exist.

**Solve time grows.** 23 -> 27 unknowns at T=500; expect ~3 min to become 4-5. E3
does two re-solves, so `--skip-e3` matters more during iteration.

## Limitations to state in the paper

**No policy rate.** There is no Taylor rule, so no conventional monetary
transmission and no "the ECB did not respond to Greece" experiment. The nominal
anchor is an assumed union-inflation stabilisation, not a modelled rule.

**Sovereign bonds stay real while deposits are nominal.** A deliberate asymmetry
that maximises banks' Fisher exposure. Euro-area sovereign debt is nominal; say so.

**Wages are flexible**, so adjustment is not shifted from wages onto hours and the
model is silent on that component of distributional incidence.

**Steady-state markups are subsidised away**, so there is no profit-income level
effect — the markup rent is distributed proportional to `e` and nets out of
household income by construction.

## Benchmark: Bi, Foerster and Traum (2026)

*"Asset Purchases in a Monetary Union With Default and Liquidity Risks", FRBSF
Working Paper 2025-10, https://doi.org/10.24148/wp2025-10.* Closest published
analogue: two-country monetary union, Gertler-Karadi intermediaries, endogenous
sovereign default, cross-border sovereign holdings, targeted ECB asset purchases —
calibrated Italy/Germany 2012.

| | Bi-Foerster-Traum | This spec |
|---|---|---|
| Price rigidity | Rotemberg, exact nonlinear (their 2.14) | Rotemberg, linear-equivalent |
| Wage rigidity | none — flexible, `chi*L^sigma_l = U_c*w` (A.9) | **none** (same) |
| SS markup | live, `theta^c = 11` -> `mc_ss = 10/11` | subsidy-neutralised, SS bit-identical |
| Markup rent | lump-sum `Pi^f` to representative household | proportional to `e`, nets out of income |
| Nominal anchor | Taylor rule, `phi_pi=1.6, phi_y=0.07, phi_r=0.85` | union-inflation normalisation, no rate |
| Union weights | 0.5 / 0.5 | `omega_pi_D = 0.071` |
| ToT identity | `rer_t/rer_{t-1} = pi*_t/pi_t` (A.80), CPI form | same relation, PPI form |
| Deposits | nominal | **nominal** (same) |
| Sovereign bonds | nominal | real (deliberate asymmetry) |
| Solution | 2nd-order perturbation (endogenous regime switching) | 1st-order SSJ |

**The price slope agrees.** Their `psi` maps to a Calvo-equivalent slope
`(1-xi)(1-beta*xi)/xi = 0.0846` at `xi_p = 0.75, beta = 0.995`. This spec's
`kappa_p = 0.0871` at `beta = 0.985` is the same number to within 3%; the gap is
entirely the discount factor.

**They never report a flexible-price counterfactual.** Their Table 1 decomposes
over the liquidity-risk channel, the fiscal-limit shift and the debt change — never
over price stickiness. The staged rollout here produces both that counterfactual
and a clean split between the price-stickiness and Fisher contributions, so those
are reportable output rather than merely regression tests.

**Their nominal side does little propagation work.** Inflation moves +/-0.1% while
investment moves 9% and output 0.6% (their Figure 3); Tables 2 and 3 show inflation
at 0.00-0.02 against investment at 0.53. Their §4.1 treats inflation as an *outcome*
of the relative-price move, never as a channel. What drives their output contraction
is the **loan-in-advance constraint** (`eta^I = 0.65/0.75` of investment must be
debt-financed, their 2.10) — a real financial friction this model does not have.
Against their `-0.6%` output impact, this model's `Y_D[0] = -0.0149%` is two orders
of magnitude smaller. If steps 2 and 3 leave `Y_D[0]` implausibly small, a Sims-Wu
working-capital constraint is the natural next lever; it is **out of scope here**
and would need its own design pass.

## Relationship to the `add-nkwpc` branch

`add-nkwpc` (commit `2377f79`, off `08e1010`, pre-reorganisation) is a single
26-line commit adding `wage_setting_{D,F}` only. It was never wired into the model
list and has no `kappa_w` calibration. **Nothing from it is used** — this spec is
price-side only and the branch has no price-side content. Recorded so it is not
revisited without knowing what is in it:

- It is a **real-wage** Rotemberg curve (`pi_w = w/w(-1) - 1` on the *real* wage,
  explicitly "no CB needed") — a real adjustment friction, not a nominal rigidity.
- It divides the MRS by `UCE_D`, correct under separable preferences but **wrong
  under the GHH preferences this model uses** — cf. `labor_ss_D` and
  `labor_market_D`, neither of which contains `UCE`.
- It uses `w_D` where the existing labour condition uses `w_D / P_CES_D`, dropping
  the CES bundle deflator.
- **Its stated motivation is wrong for this model.** The comment says sticky wages
  stop "the household wealth effect" translating into `N`. GHH preferences have no
  wealth effect on labour supply by construction.
