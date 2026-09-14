# Runtime Function Reference — `code/global/`

> **STALE — documents machinery that no longer exists.** This reference was written for
> the perfect-foresight pipeline (`transition.py`, `risk_branch.py`, `solvers.py`), which
> was **deleted in commit `b1f0b81`**, and for the mechanical CB price floor
> (`_cb_price_floor`, `psi_cb_D`, `cb_buy_D`, `Q_floor_D`, `Q_bD_free`, `recap_D`,
> `recap_path`), removed on 2026-08-30 with the TPI rework. Anything below describing
> those objects is history, not current behaviour.
>
> For what actually runs, read `CLAUDE.md` and the module headers under
> `code/global/solver_recursive/`. In particular the TPI backstop is now a **one-sided
> yield peg with real purchases** (`phi_tpi`, `Q_peg_D`, the `x_cb` unknown, the `b_cb`
> state and the `rem_cb_D` remittance in `point_map.py`), not a portfolio-balance price
> floor. The sections below on `Q_floor_D` / `cb_buy_D` / `psi_cb_D` describe the
> superseded design.


**Scope.** This document is a practical reference for the principal functions on the
main execution path of the two-country HANK–GK monetary-union model
(`python3 code/global/main.py`). It covers the `transition`, `risk_branch`, `bank`,
`household`, `firms`, `government`, `capital`, and `distribution` modules and, for
each core function, records its interface, internal logic, economic purpose, and
computational properties. The trade block (`trade.py`), steady-state driver
(`steady_state.py`), Newton solver (`solvers.py`), numba kernels (`fast_kernels.py`)
and plotting (`plots.py`) are referenced where they interact with these modules but
are not documented in detail here.

**Status.** Reflects the code as of 2026-07-21 (branch `bocola-rewrite`), including
the TPI (Transmission Protection Instrument) extension and the 2026-07-21 clarity
cleanup. Parameter values cited are the current entries in `calibration.py` and may
drift; treat parameter *names* as authoritative and re-check values there.

**Superseded content.** Everything in this document that referenced Cole-Kehoe
self-fulfilling crisis zones, the `sunspot_*` shock, `solve_transition_ck`/
`solve_transition_ck_risk`, `ck_default_prob`, `chi_tilt`, or a selectable
`cal["sdf_mode"]` (`"income"`/`"empirical"`/`"model"`) described the **pre-rewrite**
architecture (branch `global`) and has been removed. The 2026-07-21 cleanup
additionally deleted the branch scarring flags (`def_output_cost_D`,
`def_output_rho_D`, `def_capital_quality_D` and its `quality0` plumbing,
`recap_share_D`), the `pin_rdep` deposit-rate diagnostic, the `hybr_factor`
argument, `capital_branch_summary`, and `hm_bond_price_ss`/`hm_bond_return_ss` —
all of which were switched off or unused at the baseline. The 2026-07-16 "Bocola-faithful
rewrite" replaced the CK crisis-zone wrapper with an exogenous priced-default-
probability path (`def_price_D`, never a function of debt) and the single always-
recomputed two-branch (now three-branch, with TPI) kernel in `bank.py`. See
`CLAUDE.md` and git history on `bocola-rewrite` for the full rationale.

---

## Table of contents

1. [Execution overview](#1-execution-overview)
2. [Notation and conventions](#2-notation-and-conventions)
3. [Module `transition`](#3-module-transition)
   - [`_inner_economy`](#_inner_economy)
   - [`make_residual`](#make_residual)
   - [`solve_transition`](#solve_transition)
   - [`market_residuals`](#market_residuals)
4. [Module `bank`](#4-module-bank)
   - [`steady_state_bank`](#steady_state_bank)
   - [`bank_backward`](#bank_backward)
   - [`bank_forward`](#bank_forward)
5. [Module `risk_branch`](#5-module-risk_branch)
   - [`extract_init_state`](#extract_init_state)
   - [`solve_default_branch`](#solve_default_branch)
   - [`make_risk_inputs`](#make_risk_inputs)
   - [`solve_tpi_branch`](#solve_tpi_branch)
   - [`make_tpi_inputs`](#make_tpi_inputs)
   - [`solve_transition_risk`](#solve_transition_risk)
   - [`bond_decomposition`](#bond_decomposition)
6. [Module `government`](#6-module-government)
   - [`govt_steady_state`](#govt_steady_state)
   - [`govt_transition`](#govt_transition)
7. [Module `capital`](#7-module-capital)
   - [`gamma_params`](#gamma_params)
   - [`capital_demand`](#capital_demand)
   - [`solve_capital_path`](#solve_capital_path)
8. [Module `firms`](#8-module-firms)
   - [`markup_ss`](#markup_ss)
   - [`steady_state_firm`](#steady_state_firm)
   - [`solve_firm_path`](#solve_firm_path)
9. [Module `household`](#9-module-household)
   - [`make_asset_grid`](#make_asset_grid)
   - [`egm_step`](#egm_step)
   - [`solve_steady_state_household`](#solve_steady_state_household)
   - [`solve_backward_transition`](#solve_backward_transition)
10. [Module `distribution`](#10-module-distribution)
    - [`get_lottery_weights`](#get_lottery_weights)
    - [`forward_iterate`](#forward_iterate)
    - [`forward_paths`](#forward_paths)
    - [`stationary_distribution`](#stationary_distribution)
    - [`aggregate_assets` / `aggregate_consumption`](#aggregate_assets--aggregate_consumption)
11. [Call-graph summary](#11-call-graph-summary)

---

## 1. Execution overview

`main.py` runs four sections in order:

1. **Steady state (always runs).** `steady_state.solve_steady_state(cal)` solves the
   symmetric two-country steady state (two-stage: capital markets + current
   account, then deposit markets). Every transition experiment starts from and
   terminates at this steady state.
2. **TFP shock (`RUN_TFP`).** A 1% AR(1) TFP shock in country D
   (`Z_D[t] = Z_ss · exp(0.01 · 0.9^t)`) is fed to [`solve_transition`](#solve_transition)
   with no default risk — the baseline real-shock IRF used to validate the
   perfect-foresight machinery.
3. **Bocola sovereign-risk pass-through (`RUN_RISK`, centerpiece).** An exogenous
   *priced* default-probability path `π_t = 0.01·0.95^t` (Bocola's s-shock analog,
   eqs. 11–12 — an input path, **never** a function of debt) is fed to
   [`solve_transition_risk`](#solve_transition_risk), which solves the fixed point
   between the no-default base path and a representative post-default branch;
   [`bond_decomposition`](#bond_decomposition) then splits the sovereign spread into
   default compensation, risk premium, and liquidity premium (plotted as a
   standalone figure). Default is priced but never realized (`def_real ≡ 0`):
   pure pass-through.
4. **TPI backstop (`RUN_TPI`).** The *same* sovereign-risk shock as experiment 3,
   plus a Markov-switching central-bank backstop on D-sovereign bonds: a priced
   probability path that the backstop holds (`pi_tpi_D_path`) weights a THIRD
   representative branch (the backstop reneging) inside
   [`solve_transition_risk`](#solve_transition_risk), and a realized activation
   path (`s_tpi_D_path`) drives a mechanical price-floor override inside
   [`bank_backward`](#bank_backward). `prints.print_tpi_table` reports the
   intervention size and the compression relative to the no-TPI run at the same
   shock.

All console output — the steady-state table, the market-clearing residual checks,
and both experiment diagnostic tables — is formatted in `prints.py`; `main.py` only
orchestrates the solves and the figures.

The solver hierarchy, from outermost to innermost:

```
solve_transition_risk         (risk_branch)  base ↔ default-branch ↔ TPI-branch
 │                                            fixed point (whichever are "live")
 ├─ solve_default_branch      (risk_branch)  ONE representative default event
 ├─ solve_tpi_branch          (risk_branch)  ONE representative "backstop reneged" event
 └─ solve_transition          (transition)   7T-unknown Newton (damped, jac_cache-
     └─ residual (nested in                  reused) → hybr fallback → Newton polish
        make_residual)        (transition)   one full economy per evaluation
         └─ _inner_economy    (transition)   firms → capital → bank_backward → govt →
                                              [bond clearing + CB remittance] →
                                              bank_forward → households → distribution
                                              → trade
```

Every Newton residual evaluation solves a *complete* general-equilibrium economy
given the 7T guessed paths — there are no inner fixed points besides the household
backward/forward passes, which are direct (non-iterative) given prices. The only
*outer* fixed points left in the model are the branch ↔ base-path loops in
`solve_transition_risk` (default risk, and now TPI); there is no crisis-zone
indicator loop (removed with Cole-Kehoe).

## 2. Notation and conventions

| Symbol / suffix | Meaning |
|---|---|
| `D`, `F` | Country suffixes: D = domestic/periphery (Greece), F = foreign/core (Germany). |
| `T` | Transition horizon in quarters (`cal["T"]`, currently 200). |
| `p` | Relative price of the F good in D goods (terms of trade / real exchange rate); `p_ss = 1` in the symmetric steady state. |
| `P_CES` | CES consumption-basket price index in units of the home good. |
| `N`, `Kap` | Aggregate employment and end-of-period capital stock. |
| `rdep` | Deposit rate **set at t, paid at t+1** (predetermined: the rate received at t was locked at t−1). |
| `Q`, `rk` | Price of capital (Jermann) and realized return on capital claims. |
| `Q_bD`, `Q_bF` | Hatchondo-Martinez perpetuity prices for D- and F-government bonds. `Q_bD_free` is the pre-TPI-floor price; `Q_bD` (used everywhere downstream) is the post-floor, actually-traded price. |
| `Q_floor_D` | The TPI "fundamental-only" price floor: prices default compensation alone (zeroing the risk premium and the liquidity/IC spread) at the steady-state IC spread. |
| `def_price`, `def_real` | *Priced* default probability (enters `Q` and expected-return FOCs) vs *realized* haircut indicator (enters realized returns and government flows). Only D is default-risky. |
| `pi_tpi_D`, `s_tpi_D` | *Priced* probability the CB backstop remains active (enters the three-branch kernel; adverse "reneged" weight is `1−pi_tpi_D`) vs the *realized* mechanical-activation indicator (drives the price-floor override). Mirrors the `def_price`/`def_real` split exactly. |
| `cb_buy_D` | Realized CB purchase quantity of D-bonds (closed-form, not a Newton unknown). `rem_cb_D` is the CB's own net cash flow each period (coupon + continuation value of last period's holding, minus this period's purchase cost), rebated lump-sum to households by SS-GDP share. |
| `psi_cb_D` | Portfolio-balance elasticity translating a `Q_floor − Q_bD_free` price gap into `cb_buy_D`. |
| `b_gov`, `b_gov_eop` | Government bond stock at beginning / end of period. |
| `n`, `alpha`, `mu` | Bank net worth, franchise value per unit of net worth (V/n), and incentive-constraint multiplier. |
| `vN` | GHH labour disutility `χ·N^(1+1/frisch)/(1+1/frisch)`; the GHH composite is `x = c − vN`. |
| `ss`, `cal` | Steady-state dict from `solve_steady_state`; calibration dict from `get_calibration`. |

Bond denomination convention (from `calibration.py`): D-bonds are D-good claims
priced off `rdep_D`; F-bonds are F-good claims priced off `rdep_F`. Cross-border
positions convert at `p` (e.g. the D-bank's F-bond leg in D-goods is `p·Q_bF·b_F_D`).
The CB rebate's F-share is a genuine cross-border transfer and is converted the same
way (`rebate_F = share_F·rem_cb_D / p`).

---

## 3. Module `transition`

**File:** `code/global/transition.py`

The nonlinear perfect-foresight (MIT-shock) transition solver. Stacks 7 unknown
paths of length T into a single vector and solves 7T market-clearing residuals
with a damped Newton method (`solvers.newton_solve`, explicit finite-difference
Jacobian, `hybr` fallback). Three design decisions define the module:

- **Endogenous debt inside every residual call.** Bond prices come from bank
  marginal conditions alone (`bank_backward`); the debt stock is then
  forward-integrated from the government budget identity
  ([`govt_transition`](#govt_transition)), and banks clear the bond market against the
  *true* end-of-period stock, net of any CB purchase. This keeps Walras exact when
  beliefs (or the CB) move the debt stock held by banks — clearing against a fixed
  `B_gov_ss` re-opens a leak of ~0.5% of GDP per 5% debt deviation.
- **Predetermined capital (Bocola eq. 6).** The stock producing at t is the one
  carried INTO t (`Kap_prod[t] = Kap[t-1]`), so impact output moves through hours
  alone — see [`solve_capital_path`](#solve_capital_path) / [`solve_firm_path`](#solve_firm_path).
- **Walras redundancy.** The F goods market and the current account are *dropped*
  from the residual system and only monitored as diagnostics (thresholds:
  goods_D ≤ 5e−9 imposed, goods_F ≤ 2e−6 diagnostic — both including a CB-rebate
  cross-border transfer term when TPI is active, see [`_inner_economy`](#_inner_economy)).

### `_inner_economy`

```python
_inner_economy(N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p_path,
               Z_D_path, Z_F_path, ss, cal,
               def_price_D=None, def_real_D=None,
               init=None, risk_D=None,
               tpi_D=None, s_tpi_D=None) -> dict
```

| Input | Type | Role |
|---|---|---|
| `N_*, Kap_*, rdep_*, p_path` | `(T,)` arrays | The 7 guessed unknown paths. |
| `Z_*_path` | `(T,)` arrays | Exogenous TFP paths. |
| `def_price_D`, `def_real_D` | `(T,)` or `None` | Priced default probability / realized haircut paths. **D only** — F never defaults. |
| `init` | dict or `None` | Mid-crisis initial conditions (lagged states) for default/TPI branches and policy runs; `None` ⇒ start from steady state. Gains `cb_buy_D_lag0` (the CB's carried-over D-bond holding entering the launch date) when TPI is in play — see [`extract_init_state`](#extract_init_state). |
| `risk_D` | dict or `None` | Bocola two-branch risk inputs for `bank_backward` (see [`make_risk_inputs`](#make_risk_inputs)); `None` ⇒ risk-neutral default pricing. |
| `tpi_D` | dict or `None` | Priced TPI-reneging inputs for `bank_backward` (see [`make_tpi_inputs`](#make_tpi_inputs)); adds a third branch to the kernel. Requires `risk_D` also be set. |
| `s_tpi_D` | `(T,)` or `None` | REALIZED mechanical-purchase activation path — independent of `tpi_D`; drives the price-floor override in `bank_backward` on its own. |

**Returns** a dict of all endogenous block outputs: firm/capital/bank/government
sub-dicts, dividends, CES price indices, household policies and income paths,
aggregate `C`/`A` paths, trade flows, `rem_cb_D` (CB net cash flow), and the
sequence of start-of-period cross-sectional distributions `D_start_*` (shape
`(T+1, n_a, n_e)`).

**Logic** — one full economy per call, evaluated block by block in dependency order:

1. **Firms.** [`solve_firm_path`](#solve_firm_path) maps `(N, Kap_prod, Z)` into
   `Y`, `w` (frictionless), `mpk` — all on the **predetermined** capital vintage
   `Kap_prod[t] = Kap_lag` at t=0, `Kap[t-1]` thereafter.
2. **Capital.** [`solve_capital_path`](#solve_capital_path) inverts the Jermann
   accumulation technology on the guessed `Kap` path (bought at t, producing at
   t+1) to obtain investment `I`, the capital price `Q`, the realized return `rk`,
   and capital-producer profit.
3. **CES price indices.** `trade.ces_price(p)` per period, per country.
4. **Bank backward pass.** [`bank_backward`](#bank_backward) computes, from
   expected-return FOCs under the priced default probability (and, if `risk_D`/
   `tpi_D` are set, the two/three-branch expectations, plus the TPI price floor if
   `s_tpi_D` is set): bond prices `Q_bD` (post-floor), `Q_bF`, `Q_bD_free`,
   `Q_floor_D`, `cb_buy_D`, franchise values `alpha`, IC multipliers `mu`, discount
   factors `Omega`, and the *cross-border* bond holdings `b_D_F`, `b_F_D` from
   portfolio FOCs.
5. **Working-capital wedge (Neumeyer-Perri).** With the IC multiplier now known,
   form `r_wc_t = rdep_{t−1} + λμ_t/Ω̃_t` and divide the firm wage by
   `1 + ζ_wc·r_wc_t` (`ζ_wc = zeta_wc_*`; `ζ_wc = 0` is a no-op). The lowered
   wage feeds the labour-market residual and household income — the
   spread→output transmission channel, and (with predetermined capital) the ONLY
   channel from spreads into impact output. The financing income is accumulated
   for step 8 (routed to households as dividends, not onto the bank balance sheet).
6. **Government.** [`govt_transition`](#govt_transition) forward-integrates the
   debt stock under the Bohn tax rule at the just-computed (post-floor) bond
   prices, applying any *realized* haircuts and default-branch `recap_D_path`
   outlays. Optional `init` keys `b_gov0_*`, `b_anchor_*`, `recap_D_path` support
   mid-path and post-default starts.
7. **CB remittance.** `rem_cb_D = delta_b_D·surv_cb_D·cb_buy_D_lag +
   Q_bD·surv_cb_D·(1−delta_b_D)·cb_buy_D_lag − Q_bD·cb_buy_D` — the CB's own net
   cash flow (coupon plus the surviving continuation value of last period's
   holding, minus this period's purchase cost), using the same
   `def_real_D`/`recovery_rate_D` survival convention as bank-held bonds. Zero
   identically when `cb_buy_D ≡ 0` (TPI off). Ports the historical TPI-1 audit
   fix (`docs/audit.md`): the pre-rewrite prototype omitted exactly this term and
   leaked ~2.6% of quarterly GDP per period.
8. **Bond-market clearing against the true stock.** Domestic banks are the
   residual holders, net of both the foreign cross-border leg and any CB purchase:

   ```
   b_D_D[t] = b_gov_D_eop[t] − b_D_F[t] − cb_buy_D[t]
   b_F_F[t] = b_gov_F_eop[t] − b_F_D[t]
   ```

   A `RuntimeError` is raised if either residual holding turns non-positive —
   caught by the outer solver and converted into a penalty.
9. **Bank forward pass.** [`bank_forward`](#bank_forward) rolls net worth forward
   from *realized* returns (marked-to-market bond and capital revaluations,
   realized haircuts via `def_real_D`, plus any `recap_D` equity injection),
   producing `n`, `n_IC`, dividends `div`, and deposit supply `Dep_supply`.
10. **Dividends to households.** `Div = (1 − mc)·Y + cap_profit + div_bank +
    wc_income` (flexible-price markup is constant).
11. **CB rebate to households, split by SS-GDP share.** `rem_cb_D` (D-goods
    denominated) is split `share_D = Y_ss_D/(Y_ss_D+Y_ss_F)` to D and `share_F` to
    F. The D-share stays a purely domestic financial flow (no separate goods-
    market term needed, exactly like `Tax`/`coupon`/`net_issuance`). The F-share
    is a genuine cross-border real transfer: it is converted to F-goods via `p`
    before entering F household income, AND appears as an explicit compensating
    term in both `goods_D_resid` (see [`make_residual`](#make_residual)) and the
    `goods_F` diagnostic — omitting either reopens the TPI-1 class of leak, this
    time via a terms-of-trade conversion bug.
12. **Household income.** In composite-good units, per idiosyncratic state `e`:

    ```
    y_D_t(e) = (w_D_t/P_CES_D_t)·N_D_t·e + (Div_D_t − Tax_D_t + rebate_D_t)/P_CES_D_t
    ```

    (symmetric for F), with the GHH disutility `vN_t` passed separately.
13. **Real deposit returns (Fisher equation, predetermined rate).**
    `r_real[t] = (1+rdep[t−1])·P_CES[t−1]/P_CES[t] − 1`, built as a length-`T+1`
    array; period −1 anchors (`rdep_prev`, `P_lag`) come from `init` when the path
    starts mid-crisis.
14. **Household EGM backward.** [`solve_backward_transition`](#solve_backward_transition)
    for each country: policies `c[t]`, `a'[t]` by backward induction from the
    steady-state terminal condition (numba kernel or numpy fallback).
15. **Distribution forward.** [`forward_paths`](#forward_paths) from `init["D_*"]`
    or the stationary distribution; all T+1 start-of-period distributions are
    stored (`D_start_*`) so a default or TPI branch can be launched from any base
    date.
16. **Trade.** `trade.import_demand` per period and `trade.trade_balance` give
    `IM` and `NX` for both countries.

**Economic purpose.** This is the model's general-equilibrium map: given prices
and quantities on the 7 guessed paths, it produces every other endogenous object
consistently with agent optimization (banks, households, firms), government
policy, and (when active) the central bank's backstop. The block ordering embodies
the model's causal structure under perfect foresight: prices from marginal
conditions (backward passes) precede flows and stocks (forward passes).

**Computational aspects.**
- No internal iteration: every block is a direct computation given its inputs, so
  cost is linear in T. The dominant costs are the household EGM and the
  distribution forward pass (both numba-JITed when available, `cal["use_numba"]`).
- Errors from infeasible guesses (negative Jermann bracket, non-positive bond
  holdings, `mu ≥ 1`, NaN powers) are raised as exceptions, not returned as NaN, so
  the outer solver can penalize immediately.
- The distinction between the bank *backward* pass (expected returns, FOCs, priced
  probabilities) and *forward* pass (realized flows, net worth, realized haircuts)
  is what implements the PRICED vs REALIZED default (and TPI) split.

### `make_residual`

```python
make_residual(spec, verbose=False) -> residual  # residual: (7T,) -> (7T,)
```

Builds the stacked-residual closure from a **picklable** `spec` dict (`ss, cal,
Z_D_path, Z_F_path, def_price_D, def_real_D, init, risk_D, tpi_D, s_tpi_D`) — the
same dict multiprocessing Jacobian workers unpickle to rebuild an identical
residual for parallel finite-difference columns (`solvers.fd_jacobian`).

Maps the stacked unknown vector `y` (ordering `[N_D | N_F | Kap_D | Kap_F | rdep_D
| rdep_F | p]`, each block length T) to 7T residuals:

| # | Residual (per period, normalized) | Pins |
|---|---|---|
| 1 | Capital-IC complementarity D: Fischer-Burmeister `φ(μ_D/μ_ss, slack_D/n_ss)` | `Kap_D` |
| 2 | Capital-IC complementarity F: same, country F | `Kap_F` |
| 3 | Labour market D: `(χ_D·N_D^(1/frisch) − w_D/P_CES_D) / (w_D/P_CES_D)` | `N_D` |
| 4 | Labour market F: same, country F | `N_F` |
| 5 | Union deposit clearing (D-good units, both countries' imbalances netted) | `rdep_D` |
| 6 | Deposit-UIP real-rate parity `(1+rdep_D) = (1+rdep_F)·p'/p` | `rdep_F` |
| 7 | Goods market D (incl. the CB-rebate cross-border transfer term when TPI is active) | `p` |

Notes on the residuals:

- The IC is **occasionally binding** (Bocola): μ (from the capital FOC, asset-
  agnostic under the single-λ assumption) and slack (`α·(n−n_IC)`) satisfy
  `0 ≤ μ ⊥ slack ≥ 0`, imposed via the smooth Fischer-Burmeister function
  `φ(a,b)=a+b−√(a²+b²)` rather than an always-binding equality — this is the
  general pattern any new kinked/occasionally-binding mechanism (e.g. the TPI
  price floor) should follow to keep the FD Jacobian well-behaved; see
  [`bank_backward`](#bank_backward)'s `_CB_SMOOTH_EPS`.
- Deposits are a **union-wide** market (deposit-UIP integration, replacing two
  national clearings): own-good claims at national rates, cleared once in D-good
  units, with the cross-border deposit position (`nfa_dep_D`) as the new
  absorption margin, plus real-rate parity as the second condition.
- Goods market F and the current account are *not* imposed (Walras).

Robustness devices inside `residual`:

- **Domain guard:** `p ≤ 0.05`, `N ≤ 0.01`, `Kap ≤ 0.1` return a flat penalty
  vector `np.full(7T, 10.0)` before any computation.
- **Uniform wall height:** every failure path (guard, exception in
  `_inner_economy`, non-finite residuals) returns the *same* penalty height 10.0,
  so `hybr`'s finite-difference gradient isn't biased toward one wall.

### `solve_transition`

```python
solve_transition(ss, cal, Z_D_path, Z_F_path,
                 def_price_D=None, def_real_D=None,
                 verbose=True, maxiter=300, y0=None,
                 init=None, risk_D=None, jac_cache=None, accept_tol=None,
                 tpi_D=None, s_tpi_D=None) -> dict
```

| Parameter | Default | Role |
|---|---|---|
| `y0` | flat SS paths | Initial guess for the stacked 7T unknown vector (warm start). |
| `init` | `None` | Mid-crisis initial state (passed through to `_inner_economy`). |
| `risk_D` | `None` | Bocola two-branch risk inputs (risk-neutral if `None`). |
| `tpi_D` | `None` | Priced TPI-reneging inputs (no third branch if `None`; requires `risk_D`). |
| `s_tpi_D` | `None` | Realized mechanical-purchase path (floor inactive if `None`). |
| `jac_cache` | `None` | Caller-owned dict carrying the Jacobian **across** solves — the source of the branch/risk warm-resolve speedup. `None` ⇒ a fresh local cache. |
| `accept_tol` | `None` | Max-abs residual accepted; `None` ⇒ `cal["tol_transition"]` (1e−10). Branch probes pass `1e-9`. |
| `maxiter` | `300` | Scales `maxfev` for the `hybr` fallback. |

**Returns** a flat dict of all solved paths: the 7 unknowns, all firm/capital/bank/
government/household/trade outputs (suffixed `_D`/`_F`), `rem_cb_D`, `s_tpi_D`,
`nfa_dep_D`, the default paths actually used, `mu_min_D`/`mu_min_F`/`slack_min_D`/
`slack_min_F` (the complementarity monitor), all bank-block outputs including
`cb_buy_D`/`Q_bD_free`/`Q_floor_D` (spread from `**out["bk"]`), and `y_vec` — the
solved unknown vector, used as a warm start by every outer loop.

**Logic** (solver in `solvers.py`).

1. Build the default initial guess: all seven paths flat at their steady-state
   values.
2. **Damped Newton** (`newton_solve`) on an explicit finite-difference Jacobian
   (`fd_jacobian`, built in parallel via multiprocessing), with Broyden updates
   and stall-triggered rebuilds. Reused across solves via `jac_cache`.
3. **Fallback.** If Newton stalls above `accept_tol`, fall back to
   `scipy.optimize.root(method="hybr")`, keep it only if it improves, then run a
   final Newton **polish** with a fresh Jacobian (`hybr` alone plateaus near
   `max|resid| ≈ 5e-11` on xtol, not the residual).
4. Raise `RuntimeError` if the final residual exceeds `accept_tol`; otherwise
   re-evaluate `_inner_economy` at the solution and assemble the output dict.
5. **Complementarity monitor.** After a successful solve, check
   `min(mu_D), min(mu_F), min(slack_D), min(slack_F)`; print a warning if any is
   negative beyond tolerance (a spurious corner the Fischer-Burmeister residual
   shouldn't have accepted). Checked here, *not* inside the residual.

**Economic purpose.** This is the model's equilibrium concept: a perfect-foresight
path on which the occasionally-binding bank incentive constraint, the GHH labour
FOC, union deposit-market clearing with real-rate parity, and the D goods market
all hold every period, with government debt (and, when TPI is active, the CB's own
bond position) endogenous throughout.

**Computational aspects.**
- Problem size 7T (1,400 unknowns at T=200). Each Newton iteration costs `O(7T)`
  inner-economy evaluations built in parallel; the TFP experiment runs in ~10s,
  the risk-channel experiment (with its outer branch loop) in ~40s, the TPI
  experiment (base + full three-branch fixed point) in ~100–160s.
- `y_vec` in the output enables warm starting: all outer loops
  (`solve_transition_risk`, homotopies) restart the Newton from the previous
  solution, which usually converges in a handful of iterations.
- Zero-shock regression: with flat `Z` paths the solver must stay at the steady
  state to ≤ 1e−5 (`tests/test_transition_walras.py`).

### `market_residuals`

```python
market_residuals(out, cal, ss=None) -> dict
# goods_D, goods_F, dep_union, uip, cap_D, cap_F,
# slack_min_D/F, mu_min_D/F, nfa_dep_D
```

The single definition of the market-clearing diagnostics for a solved path, shared
by `prints.print_transition_residuals` and every regression test — there is no
second copy to drift. `cap_*` is the complementarity product `μ·slack` (zero at an
exact solution), not the old always-binding gap. Pass `ss` whenever `out` carries a
nonzero `rem_cb_D` (TPI active) so the CB rebate's cross-border transfer enters
both goods residuals exactly as it does in the imposed `goods_D_resid`; omitting it
reports a Walras leak that is not in the solve.

---

## 4. Module `bank`

**File:** `code/global/bank.py`

The Gertler-Karadi/Bocola financial-intermediary block: each bank holds capital
plus domestic and foreign bonds, subject to a single-λ (asset-agnostic)
divertability constraint (Bocola 2016 eq. 3). This is the module carrying the
model's entire risk and policy pass-through logic — the two/three-branch pricing
kernel, the occasionally-binding IC, and (new) the TPI price floor all live here.

### `steady_state_bank`

```python
steady_state_bank(cal, rk_ss, Kap_ss, Q_bD_ss, Q_bF_ss,
                  b_dom_ss, b_for_ss, p_ss, country="D") -> dict
```

Solves the steady-state bank block given prices: the franchise-value fixed point
`α = Ω(1+rdep)/(1−μ)` (via `_alpha_ss_fixed_point`), the IC-implied net worth
`n_ss_IC`, the accumulation-implied net worth `n_ss_ACCUM` (must equal `n_ss_IC` at
the SS — a solved identity, not imposed), leverage `theta_ss`, and deposit supply.
`calibrate_bank_targets` (used by `steady_state.py`) inverts this fixed point to
solve for the single λ and entrant transfer `ω_ent` that hit target leverage and
credit-spread moments — note the **fold** in this fixed point at high leverage/low
spread combinations (see the calibration.py header comment on `leverage_target`).

### `bank_backward`

```python
bank_backward(rk_D, rk_F, rdep_D, rdep_F, p_path,
              cal, ss_bk_D, ss_bk_F,
              def_price_D=None, risk_D=None,
              tpi_D=None, s_tpi_D=None) -> dict
```

**Returns** `(T,)` arrays: `alpha_D/F`, `mu_D/F`, `Omega_D/F`, `Q_bD`, `Q_bF`,
`b_F_D`, `b_D_F` (cross-border FOC holdings), `ic_spread_bD_D`, `ic_spread_bF_F`,
`Q_bD_ss_val`, `Q_bF_ss_val`, and (new) `cb_buy_D`, `Q_bD_free` (pre-floor price),
`Q_floor_D` (the fundamental-only price). `cb_buy_D`/the `Q_bD`–`Q_bD_free` gap are
always present, zero/no-op when `s_tpi_D` is off, so callers never need to branch
on key existence.

**Logic.** A backward pass, `t = T−1 … 0`, carrying `alpha_next`/`Q_b*_next` as
continuation values from the SS terminal condition. At each t, prices come from
marginal conditions only (no stocks) — what lets debt be forward-integrated
afterwards. Three nested modes, selected by which of `risk_D`/`tpi_D` are `None`:

- **Risk-neutral** (`risk_D is None`): `def_price_D` enters linearly,
  `surv_D_price = 1 − defp_D_next·(1−recovery_rate_D)`,
  `Q_bD_free = surv_D_price·payoff_D_nd / (1 + rdep_D[t] + ic_spread_bD_D)`.
- **Two-branch (Bocola)** (`risk_D` set, `tpi_D is None`): with probability
  `pi_def1 = def_price_D[t+1]` the D-default event hits at t+1 and pricing/
  returns jump to the branch values (`Om_d_D[t]`, `rk_d_D`, `Q_bD_d`), discounted
  at the branch-specific kernel weight `Ω^d` instead of `Ω^nd`. `pi_def1 ≡ 0`
  collapses this bit-for-bit to the risk-neutral formula.
- **Three-branch (Bocola + TPI)** (`risk_D` and `tpi_D` both set): a THIRD branch
  — "the backstop reneged despite being priced active" — is nested inside the
  no-default continuation. Weights at each t (`pi_tpi1 = tpi_D["pi"][t+1]`):

  ```
  w_nd  = (1 − pi_def1)·pi_tpi1          # no default, backstop holds
  w_def = pi_def1                        # existing default branch, unchanged
  w_tpi = (1 − pi_def1)·(1 − pi_tpi1)    # backstop reneged
  ```

  (`w_nd+w_def+w_tpi ≡ 1` for any `pi_def1, pi_tpi1 ∈ [0,1]`.) Every base object
  (`Ω̃`, `μ`, `α`, `Q_b`, the cross-border FOC returns) becomes a three-term
  weighted mixture. `pi_tpi1 ≡ 1` (backstop never doubted) forces `w_tpi ≡ 0`
  **exactly** in floating point (`1.0−1.0=0.0`, `0.0·x=0.0`), collapsing the
  three-branch formulas to the two-branch ones bit-for-bit — this is how
  `tpi_D=None` nests exactly without a separate code path (the caller-facing
  `tpi_D is None` case builds inert placeholders — `pi_tpi_path=ones(T)`,
  `Om_tpi_D=zeros(T)` — internally and reuses the same formulas). `tpi_D`
  requires `risk_D` also be set (raises `ValueError` otherwise) — TPI pricing
  rides on the two-branch kernel, it isn't defined standalone. F-bonds are safe
  from the D-default haircut in every mode, but their price still jumps in each
  branch (safe-haven repricing); the D-default branch applies
  `surv_d = recovery_rate_D` (a genuine haircut), the TPI-reneged branch applies
  `surv_tpi_d = 1.0` (a repricing, never a haircut — reneging doesn't destroy
  the claim).

**TPI price-floor override** (independent of which pricing mode above is active —
even the risk-neutral mode supports it): if `s_tpi_D is not None and
s_tpi_D[t] > 0`, compute the "fundamental-only" price

```
Epay_def_only = (1−pi_def1)·payoff_D_nd + pi_def1·payoff_D_d      # plain probability
                                                                    # weight, no Ω
                                                                    # covariance premium
Q_floor       = Epay_def_only / (1 + rdep_D[t] + ic_spread_bD_ss)  # SS liquidity
                                                                    # spread, not the
                                                                    # actual (stressed) one
gap           = smooth_max(0, Q_floor − Q_bD_free)                 # see below
cb_buy_D_t    = gap / psi_cb_D
Q_bD          = Q_bD_free + psi_cb_D·cb_buy_D_t   # ≈ max(Q_bD_free, Q_floor)
```

i.e. the floor zeros out BOTH the risk premium (Ω-covariance effect) and the
liquidity/IC-spread wedge, keeping only default compensation — exactly the split
`bond_decomposition` already computes post-hoc. `cb_buy_D` is a closed-form
function of already-available objects, not a new Newton unknown (the stacked
system stays at 7T). The override is computed *before* the cross-border FOC (D-in-
F leg), so the F-bank sees the same post-floor price it would actually transact
at. `smooth_max(0,x) = 0.5·(x+√(x²+ε²))`, `ε = _CB_SMOOTH_EPS = 1e-5` (module
constant): a raw `max(0,·)` is non-differentiable exactly where the floor turns
on/off, which stalls the FD-Jacobian Newton solver — the same reason the IC
complementarity in `transition.py` uses Fischer-Burmeister rather than a hard
`max`. This smoothing has a systematic bias of `ε/2` at the boundary (not exactly
0), which **compounds** through the backward recursion (`Q_bD_next` carries each
period's bias forward) to an economically negligible but non-infinitesimal ~1e-4
in price level over a 200-period horizon — size any new test tolerances against
this compounded magnitude, not the raw `ε`.

`psi_cb_D` (calibration.py) is a portfolio-balance elasticity translating a
*price-level* gap into a purchase quantity — **not** the same magnitude as the
cross-border elasticities `psi_bF_D`/`psi_bD_F`, which govern much smaller
*return-differential* gaps; naively reusing `0.05` implies purchases of 25–50% of
the entire bond stock and destabilizes the Newton solve (verified). `0.5` converges
cleanly with `cb_buy_D` peaking at a plausible sub-1% of `B_gov_D_ss`.

**Economic purpose.** The pricing engine for every asset the banks hold: expected-
return FOCs under priced default risk (and, when active, priced TPI-reneging
risk), the occasionally-binding IC multiplier, and — new — a mechanical central-
bank price support that only ever raises `Q_bD`, never lowers it, and is a pure
no-op away from the floor.

**Computational aspects.** `O(T)` backward recursion, vectorized nowhere (each
period's `alpha`/`Q_b` depends on the next period's), so this is a Python loop —
still cheap relative to the household EGM. The `w_tpi≡0`/`Om_tpi_D=0` inertness
trick (rather than branching the whole formula block on `tpi_mode`) is what keeps
the existing two-branch and risk-neutral code paths byte-for-byte unchanged when
TPI is off (regression-tested in `tests/test_tpi.py::test_tpi_off_nests_baseline`
and `test_pi_tpi_one_nests_two_branch`).

### `bank_forward`

```python
bank_forward(Kap_D, Kap_F, Q_D, Q_F, rk_D, rk_F, rdep_D, rdep_F, p_path,
            b_D_D_path, b_F_F_path, bwd, cal, ss_bk_D, ss_bk_F,
            def_real_D=None, init_D=None, init_F=None,
            Q_bD_lag0=None, Q_bF_lag0=None, p_lag0=None,
            recap_D=None) -> dict
```

**Returns** `(T,)` arrays: `n_IC_D/F` (IC-implied net worth), `n_D/F` (accumulated
net worth), `rn_D/F` (realized portfolio return), `div_D/F`, `theta_D/F`
(leverage), `Dep_supply_D/F`, `rb_D/F` (realized bond returns), `b_D_D`, `b_F_F`
(echoed holdings).

**Logic.** A forward pass, `t = 0 … T−1`, rolling net worth from *realized*
returns on positions bought at t−1: bond returns `rb_D_path` use REALIZED survival
(`def_real_D`, not `def_price_D`) and the actual (post-floor) `Q_bD` path with its
own one-period lag; capital returns use `rk`. Portfolio shares (`kappa`,
`phi_bdom`, `phi_bfor`) are recomputed each period on ACTUAL net worth, not the
IC-implied level — the two coincide only when the IC binds exactly. `recap_D`
(default-branch government equity injection) adds directly to retained net worth,
not gross income. `b_D_D_path`/`b_F_F_path` are the bond quantities the caller
computed in `_inner_economy`'s clearing step (already net of any CB purchase for
D) — `bank_forward` is agnostic to *why* the bank holds less than the government
issued, it just uses the actual holding.

**Economic purpose.** Where the two-branch/three-branch *pricing* kernel's
consequences actually hit bank balance sheets: mark-to-market gains/losses on the
realized (not expected) price path determine net worth, dividends, and hence the
whole downstream real economy.

**Computational aspects.** `O(T)` forward recursion, one Python loop per country;
cheap relative to the household EGM.

---

## 5. Module `risk_branch`

**File:** `code/global/risk_branch.py`

Implements the **Bocola (2016) risk channel**, and (new) the TPI policy-risk
channel, on top of the risk-neutral base solver. Bankers discount with the
household SDF and hold state-contingent continuation values: news that a feared
event (default, or a reneged CB backstop) is more likely pairs marginal valuations
in that state with the state's asset payoffs — a covariance premium on bonds and
capital, i.e. precautionary deleveraging even when funding is cheap.

The perfect-foresight implementation ("R2-lite") averages branches at each base
date t. **One representative branch per risk source** (default at τ* = 1; TPI
reneging at τ* = 1) is solved and its period-0 objects are reused at every base
date — deliberately, to keep the model tractable (a full state-dependent branch
tree, or a 2×2 default×TPI cross, is explicitly out of scope; see CLAUDE.md).
Documented approximations: `Λ^nd ≡ beta_inter` on the base path, branch state-
dependence across base dates treated as second order, and household π-blindness
(the deposit Euler never weights any branch — risk pricing lives entirely in the
bank block, faithful to Bocola where household deposits are riskless too). Setting
`π_def ≡ 0` and `π_tpi ≡ 1`/`s_tpi ≡ 0` nests the risk-neutral model exactly
(regression-tested in `tests/test_risk_channel.py` and `tests/test_tpi.py`).

### `extract_init_state`

```python
extract_init_state(out, ss, cal, tau) -> dict
```

**Logic.** Builds the `init` dict needed to launch a transition at base period
`tau`, using only period `tau − 1` objects of the solved base path `out`:

- household cross-sectional distributions `D_D`, `D_F` (= `D_start[tau]`);
- bank states per country: lagged net worth `n_prev`, portfolio shares
  `kappa_prev = Q·Kap/n`, `phi_bdom_prev`, `phi_bfor_prev` (cross-border leg,
  converted at `p`), and the predetermined deposit rate `rdep_prev`;
- initial government stocks `b_gov0_*` (end-of-period at `tau − 1`);
- all price/stock lags the inner economy needs: `Kap_lag`, `Q_lag`, `Q_bD_lag`,
  `Q_bF_lag`, `p_lag`, `P_lag_*`;
- `cb_buy_D_lag0 = out["cb_buy_D"][tau-1]` — the CB's carried-over D-bond holding
  entering `tau`. Without this, any branch launched from a state where the base
  path already had an active CB position would silently drop the CB's legacy
  coupon/continuation income at the branch's own h=0 (its `rem_cb_D[0]` computed
  as if the CB started from nothing) — a real, if typically small, understatement
  of branch resources; found in review and regression-tested
  (`tests/test_tpi.py::test_extract_init_state_carries_cb_buy_lag`).

**Economic purpose.** The state vector of the economy at a point mid-path — the
mechanism by which default branches, TPI-reneging branches, and (in future work)
other policy interventions can start from a crisis state rather than the steady
state.

**Computational aspects.** Pure indexing; portfolio shares are normalized by
ACTUAL net worth (matching `bank_forward`'s convention, not the IC-implied level).
Two small helpers support branch launches: `_shift_path(x, tau, T)` (shift a path
forward by `tau`, padding with the final value) and `_shifted_y0(out, tau, T,
)` (the shifted base solution, used as a Newton warm start).

### `solve_default_branch`

```python
solve_default_branch(out, ss, cal, tau=1, verbose=False, y0=None,
                     jac_cache=None) -> dict
```

**Returns** a full `solve_transition` output dict for the post-default economy,
plus `recap_D_path` (always zeros — see step 4).

**One deterministic solve of a single fixed event**: a full write-down to
`recovery_rate_D` (Greek PSI, 0.45) on the WHOLE claim, realized at branch period 0
(`def_real_D[0] = 1.0`) — not a scale search over haircut sizes. The default-state
recession must arise endogenously via bank balance sheets: the model carries no
output-cost, capital-quality or recap scarring flags (they were all 0 at the
Bocola-pure baseline and were deleted in the 2026-07-21 cleanup).

**Logic.**

1. Extract the base-path state entering period `tau` via
   [`extract_init_state`](#extract_init_state).
2. **Re-anchor the Bohn rule** to the post-haircut stock:
   `b_anchor_D = b_gov0_D · recovery_rate_D`. Keeping the steady-state anchor
   would turn the haircut into a large tax-cut windfall, making default
   expansionary and flipping the risk premium's sign.
3. Solve the branch transition once with `solve_transition` (`accept_tol=1e-9`),
   warm-started from the previous round's branch solution if available, else
   `_shifted_y0`.
4. If that direct solve stalls, a **recap-share continuation ladder**
   (`_RECAP_LADDER`, 0.5 → 0.05) chains warm starts through intermediate (larger,
   easier) recap shares and then re-solves at zero recap. The intermediate solves
   are purely numerical scaffolding and are discarded; the returned branch always
   has `recap_D_path = 0`. If even that is infeasible, raise `RuntimeError`.

**Economic purpose.** The representative "feared" state: a pure-haircut sovereign
default. Its period-0 asset returns and valuations are the *default-branch
payoffs* that the base-path bankers price under `risk_D`; the wedge between them
and base-path payoffs is the entire Bocola risk channel.

**Computational aspects.** Each branch solve is a full 7T Newton problem — a
~54% net-worth wipeout is far outside the shifted-base Newton basin, hence the
continuation ladder. Round-1 cold start is the dominant cost of a risk round;
warm-starting from the previous round's branch solution plus `jac_cache` reuse
makes later rounds cost seconds.

### `make_risk_inputs`

```python
make_risk_inputs(branch, base, ss, cal) -> dict
```

**Returns** the `risk_D` dict consumed by `bank_backward` (without the
probability `pi`, which the outer loop attaches): default-branch discount loadings
`Omega_d_D`, `Omega_d_F` (shape `(T,)`), branch period-0 prices/returns `rk_d_D`,
`rk_d_F`, `Q_bD_d`, `Q_bF_d`, `p_d` (scalars, "one representative branch reused at
every date"), and `surv_d = recovery_rate_D`.

**Logic.** The banker's default-state discount factor is

```
Ω^d_X[t] = Λ^d_X[t] · [f_X + (1−f_X)·α^d_X(0)],     X ∈ {D, F}
```

with `α^d(0)` the branch's period-0 franchise value and `f` the banker exit share.
The default-state SDF `Λ^d` is an income-based proxy:
`Λ^d = beta_inter·(Y^d(0)/Y^nd_{t+1})^(−σ)` — Euler-consistent loading on
aggregate output, a robust rep-agent stand-in for the HA economy. **Sign-gated**:
if the branch isn't a recession anywhere (`Y^d(0) ≥ Y^nd_{t+1}` for some t), the
proxy falls back to `Λ^d = β` (no SDF-side premium) with a printed warning — default
must never be priced as a good state.

**Economic purpose.** Translates the solved default branch into the two numbers
per asset the bank pricing FOCs need: how much the banker *values* payoffs in the
default state (`Ω^d`) and what those payoffs *are* (`rk^d`, `Q^d`, `surv_d`).
Because `Ω^d > Ω^nd` while default-state payoffs are low, expected-discounted
returns fall and bond/capital prices are depressed beyond actuarially fair default
compensation — the risk premium.

**Computational aspects.** Cheap (array arithmetic).

### `solve_tpi_branch`

```python
solve_tpi_branch(out, ss, cal, tau=1, verbose=False, y0=None,
                 jac_cache=None) -> dict
```

**Returns** a full `solve_transition` output dict for the "backstop reneged"
economy.

**Logic.** Mirrors [`solve_default_branch`](#solve_default_branch) structurally
(same `extract_init_state`/`_shift_path`/`_shifted_y0` helpers), but is **not** a
haircut: `def_real_D` stays zero throughout. The branch is literally a plain
`solve_transition` call (`tpi_D=None, s_tpi_D=None` — no further priced TPI or
mechanical purchases *inside* the branch, exactly as the default branch carries no
further priced default risk) launched from the base path's state at `tau` —
representing "the market believed the backstop was active, but from `tau` onward
it wasn't." No recap-style continuation ladder — reneging removes a price support
rather than destroying net worth via a realized haircut, a milder GE perturbation
than the default branch's ~54% wipeout.

**Economic purpose.** The representative feared TPI event: what the economy looks
like if the CB backstop, having been priced as likely to hold, fails to
materialize. Its period-0 objects are the branch payoffs base-path bankers price
under `tpi_D`.

**Computational aspects.** A full 7T Newton problem per solve, same order of cost
as the default branch, but typically converges more easily (milder perturbation)
so rarely needs the default branch's continuation ladder.

### `make_tpi_inputs`

```python
make_tpi_inputs(branch, base, ss, cal) -> dict
```

**Returns** the `tpi_D` dict consumed by `bank_backward` (without `pi`, attached
by the outer loop): `Omega_tpi_D`, `Omega_tpi_F` (shape `(T,)`), branch period-0
`rk_tpi_d_D`, `rk_tpi_d_F`, `Q_bD_tpi_d`, `Q_bF_tpi_d`, `p_tpi_d` (scalars), and
`surv_tpi_d = 1.0` (reneging repriced the claim; it never haircuts it, unlike the
default branch's `surv_d = recovery_rate_D`).

**Logic.** Identical construction to [`make_risk_inputs`](#make_risk_inputs),
including the same income-SDF sign gate (falls back to `Λ^tpi = β` with a warning
if the reneging branch isn't a recession relative to the base continuation — in
practice this gate trips often for TPI, since reneging removes a price support
rather than destroying resources, so the branch is frequently *not* clearly
recessionary).

**Economic purpose.** Translates the solved TPI-reneged branch into the two
numbers per asset the three-branch kernel needs, the same role
`make_risk_inputs` plays for the default branch.

**Computational aspects.** Cheap (array arithmetic).

### `solve_transition_risk`

```python
solve_transition_risk(ss, cal, Z_D_path, Z_F_path, pi_D_path=None,
                      pi_tpi_D_path=None, s_tpi_D_path=None,
                      verbose=True, max_rounds=12, damp=0.5, tol=1e-3,
                      y0=None) -> dict
```

**Returns** the base-path `solve_transition` dict plus `branch` (default-branch
dict or `None`), `risk_D_inputs`, `pi_D`, `risk_converged`, and (new)
`tpi_branch`, `tpi_D_inputs`, `tpi_converged`.

**Logic.** The outer fixed point among the base path and up to two representative
branches, with three **independent, decoupled** off-states, each nesting exactly:

- `pi_D_path ≡ 0` (or `None`) — no priced default risk;
- `pi_tpi_D_path ≡ 1` (or `None`) — backstop never doubted, no TPI-reneged
  branch;
- `s_tpi_D_path ≡ 0` (or `None`) — no mechanical CB purchases.

Two cheap short-circuits before any branch machinery: if all three are off,
returns the plain risk-neutral solve; if only `s_tpi_D_path` is live (mechanical
purchases with a certain, undoubted backstop), returns a single `solve_transition`
call with `s_tpi_D` set — **no branch solve is needed at all**, since the price
floor in `bank_backward` reads only `s_tpi_D`, independent of the priced `tpi_D`
channel. This decoupling means "CB always buys, market never doubts it" costs the
same as the plain risk-neutral path.

Otherwise, one `jac_cache` per LIVE system kind (`jc_base` always; `jc_branch` if
default risk is priced; `jc_tpi` if TPI reneging is priced), and each round:

```
Round 0:  risk-neutral (or mechanical-TPI-only) base path
Round k = 1 … max_rounds:
  1. If def_live:  solve_default_branch → make_risk_inputs → damp risk_in
                    (damp = 0.5); conv_def from (Q_bD_d, rk_d_D, Omega_d_D)
  2. If tpi_priced_live:  solve_tpi_branch → make_tpi_inputs → damp tpi_in;
                    conv_tpi analogously
  3. conv = max(conv_def, conv_tpi)   (0 for any channel not live)
  4. Build risk_D: the damped risk_in dict if def_live, else — IF
     tpi_priced_live but NOT def_live — an INERT placeholder dict (pi≡0, so the
     default-branch term contributes exactly zero weight regardless of its
     values) purely to satisfy bank_backward's "tpi_D requires risk_D" precondition,
     else None.
  5. Build tpi_D: the damped tpi_in dict if tpi_priced_live, else None.
  6. Re-solve the base path with risk_D, tpi_D, and s_tpi_D_path (unconditional,
     REALIZED — never damped/iterated, passed through every round unchanged).
  7. Stop when conv < tol AND k ≥ 2.
```

Failure handling: if a branch re-solve fails after round 1, the previous round's
inputs for that branch are kept and `fixed_point_ok` is cleared (rather than
aborting — the OTHER branch, if live, keeps iterating).

**Economic purpose.** The centerpiece experiment solver. The fixed point is
economically necessary for each live branch: it is launched from the *base*
impact state (which depends on pricing), while base pricing depends on branch
objects. At convergence, bankers' expectations (over default, and — new — over
the backstop's persistence) are consistent with the states they fear.

**Computational aspects.** Cost per round ≈ one default-branch solve (if live) +
one TPI-branch solve (if live) + one warm-started base solve. Damping at 0.5
stabilizes each branch↔base map independently; `rd ≥ 2` prevents spurious one-
round "convergence" before damping has acted. The full three-way fixed point
(both branches live) typically needs more rounds than either alone (~10–12 vs
~6–9) since the two branches converge jointly, not independently.

### `bond_decomposition`

```python
bond_decomposition(out, ss, cal) -> dict
```

**Returns** `(T,)` arrays in annualized basis points (deviations from steady
state where applicable): `total_yield`, `defcomp`, `risk`, `liquidity`,
`promised_excess`.

**Logic.** An *exact* per-period identity splitting the D-bond promised excess
return:

```
payoff^nd/Q − 1 − rdep  =  (payoff^nd − E[payoff]) / Q          [default compensation]
                         + (E[payoff]/Q − 1 − rdep − λμ/Ω̃)      [risk premium]
                         + λμ/Ω̃                                 [liquidity premium]
```

where `Q = out["Q_bD"]` is the **actual, post-TPI-floor** price. In risk-neutral
mode the middle term is 0 *by the pricing equation*; with the Bocola channel it is
positive because `Ω^d > Ω^nd` depresses `Q`. **Known limitation:** when TPI is
mechanically active, the `risk` term as computed here implicitly also absorbs the
mechanical price-support wedge (`Q_bD` vs `Q_bD_free`) — the identity is still
exact, but `risk` no longer cleanly isolates only the Ω-covariance premium in that
case. A dedicated `tpi_support` fourth component (computable from the now-
available `out["Q_bD_free"]`/`out["cb_buy_D"]`) was scoped but not implemented;
flagged as a follow-up, not a correctness bug (`print_tpi_table`'s "Sov spread
peak" figure should be read as the *net* spread, support included).

**Economic purpose.** The paper's key diagnostic: how much of the sovereign
spread is actuarial default compensation vs a true risk premium vs bank
balance-sheet tightness — Bocola's own decomposition (risk channel share of the
lending-spread response, his benchmark "up to 45%").

**Computational aspects.** Pure post-processing; the identity is regression-tested
(`tests/test_risk_channel.py`) so any change to bank pricing that breaks the
decomposition fails loudly. Annualization: quarterly rate × 4e4 → bps/yr.

---

## 6. Module `government`

**File:** `code/global/government.py`

Two ingredients: **Hatchondo-Martinez (2009)** geometric-decay perpetuity bonds
and the **Bohn (1998)** tax rule. Default risk is EXOGENOUS (Bocola 2016 eqs.
11–12): the priced default probability `π_t` is an input path to the transition
solver, never a function of the debt stock — there is no Cole-Kehoe crisis-zone
machinery in this module (removed in the 2026-07-16 rewrite).

### `govt_steady_state`

```python
govt_steady_state(cal, rdep_ss, country) -> dict
# Q_B_ss, Tax_ss, b_gov_ss
```

The HM perpetuity pays coupon `δ_b` per unit of face value and the stock decays
at rate `1 − δ_b`; Macaulay duration ≈ `1/δ_b` quarters (`δ_b=0.036` ⇒ ~7y, the
Greek anchor — long duration is what converts a *priced* default probability into
large mark-to-market losses on bank balance sheets, and what makes the TPI price
floor's job non-trivial). At the steady state `Q_B_ss = δ_b/(rdep_ss + δ_b)`
discounts the coupon stream at `rdep_ss`, and the realized return equals
`rdep_ss` by no-arbitrage (checked in `tests/test_ss_identities.py`).

With a constant stock `B_gov_ss` and no default risk, maturing principal is
rolled over: issuance of `δ_b·B_gov` new bonds at price `Q_B_ss` each period,
giving the balanced-budget tax `Tax_ss = G + δ_b·B_gov·(1 − Q_B_ss)`. Used by
`steady_state.py` and as the anchor for the transition Bohn rule. (The former
standalone `hm_bond_price_ss`/`hm_bond_return_ss` helpers were inlined and
removed in the 2026-07-21 cleanup.)

### `govt_transition`

```python
govt_transition(cal, gs, Q_B_path, def_real_path, country,
                b_gov0=None, b_anchor=None, recap_path=None) -> dict
# Tax, coupon, net_issuance, b_gov (bop), b_gov_eop   — all (T,)
```

| Parameter | Default | Role |
|---|---|---|
| `gs` | — | Country's `govt_steady_state` dict (anchors). |
| `Q_B_path` | — | Bond price path from `bank_backward` — the ACTUAL (post-TPI-floor, if active) price; the government pays/issues at this price uniformly regardless of who holds the bonds (bank or CB). |
| `def_real_path` | zeros | Realized haircut indicator per period. `country="F"` callers pass `None` (F never defaults; `recovery_rate` defaults to 1.0). |
| `b_gov0` | `b_gov_ss` | Initial stock (mid-path starts). |
| `b_anchor` | `b_gov_ss` | Bohn-rule anchor; **must** be re-set to the post-haircut stock on default branches. |
| `recap_path` | zeros | Bank-recapitalization outlays, extra government spending financed by issuance. Used ONLY by `solve_default_branch`'s numerical continuation ladder; the returned branch always carries zeros. |

**Logic.** One forward pass over the budget identity; for each t:

```
surv_t     = 1 − def_real_t·(1 − recovery_rate)
Tax_t      = Tax_base + φ·(b_t·surv_t − b_anchor)             (Bohn rule on the
                                                              SURVIVING stock)
coupon_t   = δ_b · b_t · surv_t
new_bonds  = (G + recap_t + coupon_t − Tax_t) / Q_t
b_{t+1}    = (1 − δ_b)·b_t·surv_t + new_bonds
```

The Bohn rule responds to the **post-haircut** stock (a no-op when `def_real=0`);
responding to the pre-haircut stock produced a ~31%-of-GDP one-quarter tax spike
that alone made full-event default branches infeasible.

**Economic purpose.** Fiscal policy and debt dynamics. Two properties are
load-bearing: **endogenous issuance at market prices** is the amplification leg
(when beliefs — or a TPI floor — move `Q`, financing a given deficit requires
correspondingly less/more face value); **anchor re-basing on default branches**
prevents the haircut from becoming a tax-cut windfall.

**No changes for TPI.** `govt_transition` is agnostic to who holds the bonds — it
does not need (and was not given) a `cb_buy_D` argument; the CB's own separate
cash-flow ledger is computed in `transition._inner_economy` using the same
`Q_B_path`/`def_real_path`/`recovery_rate` this function uses internally, so both
sides of the CB's position use one shared, consistent survival convention.

**Computational aspects.** A single `O(T)` scalar recursion — deliberately with
*no feedback from taxes to prices inside the block*; all price feedback runs
through the outer Newton.

---

## 7. Module `capital`

**File:** `code/global/capital.py`

Capital producers with the **Jermann (1998)** concave accumulation technology.

### `gamma_params`

```python
gamma_params(cal, country="D") -> (gamma0, gamma1)
```

Pins the adjustment-cost function so that at the steady state (`ι = δ`):
`Γ(δ) = δ` (no adjustment cost) and `Γ'(δ) = 1` (hence `Q_ss = 1`). The curvature
`ξ = cal["ksi_*"]` controls the elasticity of `Q` to investment.

### `capital_demand`

```python
capital_demand(rk_ss, mc_ss, cal, country="D") -> K_ss
```

Inverts the steady-state capital FOC `mpk_ss = rk_ss + δ` to get the capital
stock consistent with a guessed required return — used by the steady-state
solver's capital-market stage.

### `solve_capital_path`

```python
solve_capital_path(Kap_path, Kap_lag0, Q_lag0, mpk_path, cal, country="D",
                   Kap_lag_path=None) -> dict
# iota, Q, rk, I, cap_profit   — all (T,)
```

**Predetermined-capital timing (Bocola eq. 6).** `Kap_path[t]` is the stock
**bought/priced at t** and producing at t+1; `Kap_lag_path[t]` is the stock
**carried INTO t** — it is both the Jermann rebuilding base AND the production
stock (`mpk_path` must be computed on it by the caller — see
[`solve_firm_path`](#solve_firm_path)). If `Kap_lag_path` is `None` it defaults to
`[Kap_lag0, Kap_path[:-1]]` (a one-period-lagged version of the guessed path
itself); `_inner_economy` always passes it explicitly.

**Logic.**

1. Invert the accumulation technology on the guessed capital path:
   `bracket_t = (Kap_t/Kap_lag_t − (1−δ) − γ1)/γ0`, `ι_t = bracket_t^(1/(1−ξ))`. A
   negative bracket raises `ValueError` immediately (fractional powers of
   negatives return NaN silently otherwise).
2. Price of capital from the capital producer's FOC: `Q_t = 1/Γ'(ι_t)`.
3. Realized return on capital claims held from t−1 to t:
   `rk_t = (mpk_t + (1−δ)·Q_t)/Q_{t−1} − 1`, with `Q_{−1} = Q_lag0`.
4. Investment `I_t = ι_t·Kap_lag_t` and capital-producer profit
   `cap_profit_t = Q_t·(Kap_t − (1−δ)·Kap_lag_t) − I_t` — **no** additional mpk-
   reconciliation term (unlike the pre-rewrite contemporaneous-capital timing):
   firms rent exactly the bank-held vintage, so there is nothing left to
   reconcile.

**Economic purpose.** `rk` is the return on the claims banks hold against firms —
the lending-spread object. With capital predetermined, mark-to-market
revaluations of capital claims still hit bank net worth, but IMPACT OUTPUT can
only move through hours (`N`), not through a contemporaneous investment boom —
this is what restores the correct-signed comovement under a sovereign-risk shock
(reverting to contemporaneous timing re-opens the comovement problem; see
CLAUDE.md "Known limitations").

**Computational aspects.** Fully vectorized `O(T)`; the only failure mode is the
negative-bracket domain error. Note the *inversion* structure: the solver guesses
`Kap` and this block backs out `I` and `Q`, rather than integrating `Kap` forward
from `I`.

---

## 8. Module `firms`

**File:** `code/global/firms.py`

Cobb-Douglas production with monopolistic competition and *fully flexible*
prices — a deliberate benchmark choice: the markup is constant, so the block is
purely static and contemporaneous *in the state it's evaluated at* (which is the
predetermined capital vintage, not the period's own guessed `Kap`).

### `markup_ss`

```python
markup_ss(cal, country="D") -> float
```

Real marginal cost under flexible prices: `mc = (ε − 1)/ε` with
`ε = cal["epsilon_*"]` the demand elasticity (ε = 6 ⇒ mc = 5/6, a 20% markup).

### `steady_state_firm`

```python
steady_state_firm(cal, Kap_ss, country="D") -> dict
```

With `N_ss = 1` normalized: `Y_ss = Z_ss·K_ss^α`,
`w_ss = mc·(1−α)·Y_ss / (1+ζ_wc·r_wc_ss)` (working-capital wedge, a constant at
SS), `mpk_ss = mc·α·Y_ss/K_ss`, `I_ss = δ·K_ss`, `C_ss = Y_ss − I_ss − G`, and
`chi = w_ss/N_ss^(1/frisch)` (GHH static FOC, pins `N_ss=1`). The returned `chi`
overwrites the calibration warm start; `Z_ss` is separately rescaled upstream to
normalize `Y_ss = 1`.

### `solve_firm_path`

```python
solve_firm_path(N_path, Kap_prod_path, Z_path, cal, country="D") -> dict
# returns Y, w, mpk (each (T,)) and mc (scalar)
```

**Logic.** Vectorized contemporaneous evaluation ON THE PREDETERMINED VINTAGE:

```
Y_t   = Z_t · Kap_prod_t^α · N_t^(1−α)
w_t   = mc·(1−α)·Y_t / N_t
mpk_t = mc·α·Y_t / Kap_prod_t
```

`Kap_prod_path[t]` is the stock PRODUCING at t — under the predetermined-capital
timing (Bocola eq. 6), the caller (`_inner_economy`) passes the LAGGED guessed
stock (`Kap[t-1]`, quality-scaled at t=0), so `mpk_t` is the marginal product of
the vintage banks bought at t−1, not of the period's own `Kap_t` guess. The wage
`w_t` returned here is the frictionless FOC wage; the **working-capital wedge**
is applied downstream in `_inner_economy` (divides by `1+ζ_wc·r_wc_t`, needs the
IC multiplier from the bank backward pass, so can't be applied here).

**Economic purpose.** Supplies output (goods-market resource constraint), the
wage (labour-market residual and household income), and the marginal product of
capital (input to the Jermann return calculation in `capital`). This is the
channel by which predetermined capital mutes the impact-output response to any
shock that would otherwise want to move `Kap_t` contemporaneously.

**Computational aspects.** Pure NumPy arithmetic, no state.

---

## 9. Module `household`

**File:** `code/global/household.py`

One-asset incomplete-markets consumption-savings block with GHH preferences,
solved by the endogenous grid method (EGM, Carroll 2006). Utility is
`u(x) = x^(1−σ)/(1−σ)` over the composite `x = c − v(N)`,
`v(N) = χ·N^(1+1/frisch)/(1+1/frisch)`. GHH kills the wealth effect on labour
supply: the labour FOC is static and handled in the transition solver, so the
household block only chooses consumption/savings. Bocola's own calibration is
log utility (`σ_D=σ_F=1.0`) — NOT Epstein-Zin, per his explicit statement.

### `make_asset_grid`

```python
make_asset_grid(cal, country="D") -> (n_a,) array
```

Power-spaced grid `a_min + (a_max − a_min)·linspace(0,1,n_a)^curve` — with
`curve > 1`, points concentrate near the borrowing constraint. Current
calibration: 250 points on `[0, 87.2]`, curvature 2.

### `egm_step`

```python
egm_step(c_next, a_grid, Pi, r_today, r_next, y_e, beta, sigma, a_min,
         vN_today=0.0, vN_next=0.0) -> (c_today, a_pol_today)
```

One backward EGM step: invert the Euler equation on the *savings* grid
(`x_endo = (β(1+r')·E[u'(x')])^(−1/σ)`), recover consumption and the endogenous
asset grid from the budget constraint, then linearly interpolate `(a_endo →
c_endo)` back onto the fixed grid per productivity state; grid points below
`a_endo[0,e]` are borrowing-constrained (`a'=a_min`, consumption absorbs the
rest). The `1e−11` floor on the GHH composite `x'` guards the fractional power
near the constraint.

**Economic purpose.** The household Euler equation under incomplete markets:
precautionary savings against idiosyncratic income risk (via the expectation over
`Π`) generates the wealth distribution and the aggregate deposit supply banks
intermediate.

**Computational aspects.** No root-finding: the Euler equation is inverted
analytically; the only numerical operation is a 1-D interpolation per income
state. Cost `O(n_a · n_e)` per step.

### `solve_steady_state_household`

```python
solve_steady_state_household(a_grid, Pi, r_ss, y_e, beta, sigma, a_min, tol,
                             maxiter=10_000, vN_ss=0.0) -> (c, a_pol)
```

Time-iterates [`egm_step`](#egm_step) with constant `(r, y, vN)` until
`max|c_new − c| < tol` (`cal["tol_hh"] = 1e−12`); raises `RuntimeError` on
non-convergence. Used inside the steady-state deposit-market clearing stage, and
as the *terminal condition* of every transition backward pass.

### `solve_backward_transition`

```python
solve_backward_transition(a_grid, Pi, r_path, y_path, c_ss, beta, sigma, a_min,
                          vN_path=None, use_fast=True) -> (c_path, a_pol_path)
# each (T, n_a, n_e)
```

| Input | Shape | Note |
|---|---|---|
| `r_path` | `(T+1,)` | Real returns; entry `t+1` is relevant for the period-t Euler equation. |
| `y_path` | `(T, n_e)` | Income by period and productivity state. |
| `c_ss` | `(n_a, n_e)` | Terminal condition: steady-state consumption policy. |
| `vN_path` | `(T,)` | GHH disutility path (zeros if `None`). |
| `use_fast` | `bool` | Dispatch to the numba kernel (`fast_kernels.hh_backward`) when available; else a pure-numpy loop. Equivalence between the two is regression-tested (`tests/test_fast_kernels.py`). |

**Logic.** Backward induction `t = T−1, …, 0`, each step one
[`egm_step`](#egm_step) with `(r_today, r_next) = (r_path[t], r_path[t+1])`. The
terminal GHH disutility uses `vN_path[-1]` — period T is permanently at steady
state.

**Economic purpose.** Household expectations under perfect foresight: policies at
every t are consistent with the entire future path of returns and incomes. Note
the documented **π-blindness**: the deposit Euler never weights any branch (default
or TPI), so there is no precautionary-savings response to state-contingent income
risk from either channel — risk pricing lives entirely in the bank block.

**Computational aspects.** `O(T · n_a · n_e)`; the numba path (`cal["use_numba"]`,
default `True`) is materially faster and is what makes a full risk/TPI outer
fixed point (many inner-economy evaluations per branch solve) tractable.

---

## 10. Module `distribution`

**File:** `code/global/distribution.py`

Non-stochastic simulation of the household cross-sectional distribution over
`(a, e)` using the **Young (2010) lottery method**: off-grid savings choices are
split across the two neighbouring grid points with weights that preserve the
mean exactly.

### `get_lottery_weights`

```python
get_lottery_weights(a_pol, a_grid) -> (idx_lo, idx_hi, w_lo, w_hi)
```

Clip the policy to the grid range; find bracketing indices with `searchsorted`
(clipped to `[1, n_a−1]`); set `w_hi = (a'−a_grid[lo])/(a_grid[hi]−a_grid[lo])`,
`w_lo = 1−w_hi` (guarding zero-width brackets). By construction
`w_lo·a_grid[lo] + w_hi·a_grid[hi] = a'`, so aggregate assets are preserved
without histogram approximation error.

### `forward_iterate`

```python
forward_iterate(D, a_pol, a_grid, Pi) -> D_next    # (n_a, n_e)
```

One period of the distributional law of motion: a **single `np.bincount`** over
flattened `(a, e)` indices scatters mass onto the lottery-weighted neighbouring
grid points for all income states at once (~4x faster than a per-income-state
`np.add.at` loop), then post-multiplies by the Markov matrix (`pre @ Pi`).

**Economic purpose.** The Kolmogorov-forward step: given today's distribution and
policies, tomorrow's distribution.

**Computational aspects.** `O(n_a · n_e)` per period plus an `(n_a×n_e)·(n_e×n_e)`
matmul. Mass is conserved exactly.

### `forward_paths`

```python
forward_paths(D0, a_pol_path, c_path, a_grid, Pi, use_fast=True) -> (A_path, C_path, D_start)
# A_path, C_path: (T,);  D_start: (T+1, n_a, n_e)
```

**Logic.** The function `_inner_economy` actually calls for the full transition
(not a bare loop of [`forward_iterate`](#forward_iterate)): dispatches to the
numba kernel `fast_kernels.dist_forward` when available, else a numpy loop.
Timing convention: `C_t` is aggregated on the **start-of-period** distribution
(the population that consumes at t), `A_t` on the **end-of-period** one (deposits
carried into t+1, matching the bank's funding leg). `D_start[t]` (the distribution
entering period t) is stored for every t, so a default or TPI branch can be
launched from any base date via `extract_init_state`.

**Computational aspects.** `O(T · n_a · n_e)`; numba-JITed alongside the
household backward pass for the same tractability reason.

### `stationary_distribution`

```python
stationary_distribution(a_pol, a_grid, Pi, pi_e_stationary, tol,
                        maxiter=100_000) -> D    # (n_a, n_e)
```

Fixed-point iteration of [`forward_iterate`](#forward_iterate) under the
*steady-state* policy, initialized as a point mass at `a_min` distributed by the
ergodic income distribution `pi_e_stationary` (Rouwenhorst). Converges when
`max|D_new − D| < tol` (`cal["tol_dist"] = 1e−12`).

**Economic purpose.** The invariant wealth distribution — the initial condition
`D_start[0]` of every transition.

### `aggregate_assets` / `aggregate_consumption`

```python
aggregate_assets(D, a_grid)     = Σ_{i,e} D[i,e] · a_grid[i]
aggregate_consumption(D, c_pol) = Σ_{i,e} D[i,e] · c_pol[i,e]
```

Distribution-weighted sums.

---

## 11. Call-graph summary

```
main.py
├── get_calibration()                          calibration
├── solve_steady_state(cal)                    steady_state
│     [uses: firms.steady_state_firm, capital.capital_demand,
│      government.govt_steady_state, household.make_asset_grid /
│      solve_steady_state_household,
│      distribution.stationary_distribution, bank.steady_state_bank /
│      calibrate_bank_targets]
├── run_tfp:   solve_transition(...)           transition
├── run_risk:  solve_transition_risk(..., pi_D_path)              risk_branch
│     ├── solve_default_branch(...)            risk_branch [iterated: rounds]
│     │     ├── extract_init_state(...)
│     │     └── solve_transition(..., init, def_real_D)
│     ├── make_risk_inputs(...)                risk_branch
│     ├── solve_transition(..., risk_D)        transition  (two-branch re-solve)
│     └── bond_decomposition(...)              risk_branch (spread-decomposition
│                                              figure: plots.plot_bond_decomposition)
└── run_tpi:   solve_transition_risk(..., pi_D_path, pi_tpi_D_path,
               s_tpi_D_path)                    risk_branch  (base ↔ default ↔
        │                                                     TPI-reneged fixed point)
        ├── solve_default_branch(...)          risk_branch  [iterated, if def_live]
        ├── solve_tpi_branch(...)              risk_branch  [iterated, if tpi_priced_live]
        │     └── extract_init_state(...)      (now also threads cb_buy_D_lag0)
        ├── make_risk_inputs(...) / make_tpi_inputs(...)         risk_branch
        ├── solve_transition(..., risk_D, tpi_D, s_tpi_D)        transition
        └── bond_decomposition(...) / prints.print_tpi_table(...)  risk_branch / prints

prints.py (all console output; called only from main.py)
├── banner, print_ss_table, print_transition_residuals
├── print_risk_table, print_tpi_table
└── lending_spread_bps                        (also reused by plots.py)

solve_transition (every residual evaluation):
  _inner_economy
  ├── firms.solve_firm_path            (×2 countries; predetermined Kap_prod)
  ├── capital.solve_capital_path       (×2; predetermined timing, Kap_lag_path)
  ├── trade.ces_price
  ├── bank.bank_backward               (prices, FOC holdings; def_price, risk_D,
  │                                     tpi_D, s_tpi_D — TPI price floor here)
  ├── government.govt_transition       (×2; debt forward-integrated, Bohn tax;
  │                                     ACTUAL post-floor Q_B_path)
  ├── [CB remittance: rem_cb_D]
  ├── [bond clearing: b_dom = b_gov_eop − b_foreign (− cb_buy_D for D)]
  ├── bank.bank_forward                (net worth, dividends, deposits; def_real)
  ├── [CB rebate split by SS-GDP share, F-share p-converted]
  ├── household.solve_backward_transition  (×2; EGM, numba or numpy)
  ├── distribution.forward_paths           (×2, numba or numpy)
  └── trade.import_demand / trade_balance
```

**Regression anchors** (run before and after touching any of these functions):
`tests/test_ss_identities.py`, `tests/test_bank_block.py`,
`tests/test_fast_kernels.py`, `tests/test_transition_walras.py`,
`tests/test_signs_bocola.py`, `tests/test_risk_channel.py`, `tests/test_tpi.py`.
Acceptance thresholds are listed in `CLAUDE.md`.
