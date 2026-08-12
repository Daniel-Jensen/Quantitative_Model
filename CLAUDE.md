# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Two-country heterogeneous-agent model of a monetary union with Gertler-Karadi
financial intermediaries and sovereign default risk, calibrated to the
2010–2012 Greek sovereign debt crisis. The default mechanism follows
**Bocola (2016, JPE) "The Pass-Through of Sovereign Risk"**: an EXOGENOUS
rise in the *priced* probability of default π_t (his s-shock, eqs. 11–12 —
an input path, never a function of debt) makes bond prices fall, banks take
mark-to-market losses, the single-λ occasionally-binding incentive
constraint tightens, lending spreads rise and output falls — with no default
ever realized. Only D is default-risky; F bonds are safe. Application: ECB
asset purchases (TPI). Primary output is a research paper (Overleaf:
https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c).
(The 2026-07-16 Bocola-faithful rewrite replaced the earlier Cole-Kehoe
crisis-zone wrapper, the always-binding IC, and the patched default branch;
see git history on branch `bocola-rewrite`.)

## Environment

Plain `python3` (numpy/scipy/matplotlib). **Do not use the old
`/opt/anaconda3/envs/ssj` environment or the `sequence_jacobian` library** —
that was the previous implementation (see "History" below); the path no
longer exists.

## Model code (`code/global/`)

**Package layout.** The modules are grouped into subpackages; imports are
absolute from the `code/global/` root (`from blocks.bank import …`). `main.py`
sits at the root (run `python3 main.py`) and is the CHEBYSHEV-SMOLYAK PROJECTION
driver — there is NO perfect-foresight / representative-branch machinery
(`solver_pf/` was deleted 2026-08-11; git history preserves it).
- `main.py` — projection driver: SS → TFP → risk pass-through → OMT/TPI (all recursive)
- `config/` — `calibration.py`, `steady_state.py`
- `blocks/` — economic blocks (solver-agnostic): `bank.py`, `government.py`,
  `household.py`, `distribution.py`, `rouwenhorst.py`, `fast_kernels.py`,
  `firms.py`, `capital.py`, `trade.py`
- `solver_recursive/` — the ONLY solver: recursive global solution (Smolyak
  time iteration over a 7-state grid): `state_grid.py`, `decision_rules.py`,
  `point_map.py`, `expectations.py`, `recursive_main.py`, `recursive_residual.py`,
  `recursive_experiment.py` (risk + TFP), `tpi_recursive_experiment.py` (OMT/TPI)
- `reporting/` — `prints.py` (SS table), `plots.py` (activation-IRF figure)
- `tests/` — regression suite

The model is solved GLOBALLY as recursive decision rules on a Smolyak sparse
grid (Chebyshev interpolation), over the 7-state vector
`[K_D, K_F, P_D, P_F, B_D, s, Z_D]` — two capital stocks, two banks' gross
deposit obligations, the risky D-debt stock, the sovereign-risk factor s, and
the TFP state Z_D (deterministic AR(1); the TFP experiment reads the IRF along a
Z-decay path). At each grid point the SEVEN market-clearing unknowns
`[N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p]` are solved (the per-period image of
the old stacked system) with Bocola's closed-form occasionally-binding μ; α/Q_b
and the household aggregates are read off the recursions. Expectations are
genuine multi-branch Gauss-Hermite quadrature over the s-innovation × the default
fork (no-default / default-honoured / default-reneged). Driver: time iteration
(`recursive_main.time_iteration`), converging at μ=1 (15 grid points; μ=2 does
NOT converge — methods note §8). Hot kernels (household EGM backward, distribution
forward) are numba-JITed with an exact pure-numpy fallback (`cal["use_numba"]`).

| File | Contents |
|------|----------|
| `calibration.py` | All parameters. Single λ per bank (Bocola IC); Bocola/Greece anchors documented inline. f = exit/payout share; Ω = β·[f + (1−f)α′] (Bocola's ψ = 1−f survival weight on the franchise value). |
| `steady_state.py` | Two-stage SS solve: {rk_D, rk_F, p} on capital markets + current account, then {β_D, β_F} on deposit markets. Symmetric SS required (see docstring). |
| `bank.py` | GK/Bocola bank block. `bank_backward` (α, μ, bond prices, cross-border FOC holdings), `bank_forward` (net worth, dividends, deposit supply; portfolio shares on ACTUAL net worth). PRICED (`def_price_D`) vs REALIZED (`def_real_D`) default split; only D is risky, F bonds are safe. |
| `government.py` | HM perpetuity bonds, Bohn rule. `govt_transition` forward-integrates the debt stock in one pass. Default risk is exogenous (no crisis zones). |
| `solver_recursive/point_map.py` | The per-point period map (image of the old stacked system): 7 market-clearing residuals at one grid point given the frozen continuation rules. Bocola closed-form μ; THREE-branch default quadrature (no-default / honoured / reneged); Z_D read from the state; TPI via `tpi_activation`/`recovery_tpi_D`/`tpi_real_shield`. |
| `solver_recursive/state_grid.py` | Smolyak sparse grid + Chebyshev basis; `build_state_box` (7-state box incl. Z_D), `default_prob`, `s_process_params` (s-process + deterministic Z-process). |
| `solver_recursive/recursive_main.py`, `recursive_experiment.py`, `tpi_recursive_experiment.py` | Time-iteration driver; the risk + TFP experiments; the OMT/TPI activation comparison. |
| `fast_kernels.py` | numba kernels for EGM backward + distribution forward; exact numpy fallback when numba is absent (`cal["use_numba"]`). |
| `household.py`, `distribution.py` | EGM with GHH utility; stationary distribution and forward iteration. |
| `firms.py`, `capital.py`, `trade.py` | Flexible-price production with the Neumeyer-Perri working-capital wedge (w ÷ (1+ζ·r_wc), the spread→output channel; ζ=0 nests exactly — Bocola §V.C's own open-economy fix), Jermann adjustment costs, CES/Armington trade. |
| `prints.py` | Console reporting: `banner`, `print_ss_table` (the recursive experiments print their own IRF tables). |
| `plots.py` | `plot_activation_irf` (the OMT/TPI activation overlay), written to `output/`. |
| `main.py` | Projection driver: SS → TFP → risk pass-through → OMT/TPI, each a full time-iteration solve. Heavy by design (~20–30 min). |
| `tests/` | Regression suite (see below). |

## Running and testing

```bash
cd code/global
python3 main.py                                       # full projection pipeline (SS+TFP+risk+TPI), ~20-30 min
python3 -m solver_recursive.recursive_experiment      # risk pass-through only
python3 -m solver_recursive.tpi_recursive_experiment  # OMT/TPI activation only (0/50/100%)
python3 tests/test_ss_identities.py          # SS theory identities (fast)
python3 tests/test_bank_block.py             # bank FOC/no-arbitrage identities (fast)
python3 tests/test_fast_kernels.py           # numba/numpy kernel equivalence (fast)
python3 tests/test_state_grid.py             # Smolyak grid exactness (fast)
python3 tests/test_recursive_nesting.py      # recursive SS rest point + pi=0 nesting
```

**Comment convention** (enforced across `code/global/`): every module and every
function carries exactly ONE leading ALL-CAPS comment saying what it is; any
further explanation is lowercase `#` comments attached to the specific hard
line. No docstrings, no bold markers, no prose blocks inside function bodies.
Console output lives in `prints.py`, never inside the model blocks.

**Acceptance thresholds** (all enforced in tests):
- goods_D (imposed) ≤ 1e−9; goods_F (Walras-redundant diagnostic) ≤ 2e−6 —
  including when the debt stock moves.  (The Newton solver typically lands
  goods_D near 1e−13; acceptance is `tol_transition` = 1e−10 normalized.)
- Zero-shock transition stays at SS to ≤ 1e−5.
- Risk-only shock (exogenous π): Q_bD↓, n_D↓, n_F↓, Y_D[0]↓, C_D[0]↓,
  lending spread↑, b_gov↑, Tax↑ (a positive Y or n response to sovereign
  risk = bug).
- Complementarity on every solved path: μ ≥ 0, slack = αn − λ·assets ≥ 0,
  μ·slack ≈ 0 (`out["mu_D/F"]`, `out["slack_D/F"]` from point_map.py).
  Known open item: risk-on n_D[0] can sit above
  risk-off (M1 deposit-rate channel; test warning, not assert) and
  post-impact Y_D runs mildly positive — both die with the union deposit
  market (docs/sunspot_transition_study.md §8).

## Key modelling choices — do not "fix" without checking docs/SPEC.md

- **Single λ (Bocola 2016 eq. 3):** all three asset classes carry the same
  divertability. Diverging them re-opens the portfolio-substitution margin
  that made sovereign risk *expansionary* pre-rework.
- **Priced vs realized default:** `def_price` enters bond pricing and
  expected-return FOCs; `def_real` enters realized returns and government
  flows. The baseline experiment prices risk but never realizes it
  (Bocola's pass-through design); a realized-default variant just passes
  `def_real ≠ 0`.
- **Endogenous debt in clearing:** the government's end-of-period stock is
  forward-integrated inside every residual evaluation and absorbed by banks.
  Clearing against a fixed `B_gov_ss` instead re-opens a Walras leak of
  ~0.5% of GDP per 5% debt deviation.
- **Symmetric steady state:** country asymmetries enter through shocks only.
  An asymmetric SS (e.g. δ_b_D ≠ δ_b_F) shifts p_ss off 1 and opens an
  O(1e−4) SS goods-market wedge (p is weakly identified by external balance
  at trade elasticity 0.5; see steady_state.py docstring).
- **Occasionally-binding IC (Bocola):** the leverage constraint enters the
  stacked system as the Fischer-Burmeister complementarity between μ (from
  the capital FOC, valid in both regimes) and slack = αn − λ·assets, scaled
  by μ_ss and n_ss (FB's zero set is scaling-invariant). At the SS the
  constraint binds (μ_ss ≈ 0.02, slack = 0), where FB is smooth. Portfolio
  shares and branch initial conditions divide by ACTUAL net worth, not n_IC.
- **Ω-kernel weights (Bocola Prop. 1):** f = exit/payout share, so
  Ω = β·[f + (1−f)·α′] — weight 1−f ≈ 0.95 on the franchise value α′
  (Bocola's survival ψ). beta_inter ≈ β_hh ≈ 0.99 proxies the household SDF;
  values ≪ 1/(1+rdep) drive α_ss below 1 and mute the franchise channel
  (the pre-rewrite code had the weights swapped AND beta_inter = 0.96).
- **Risk channel = genuine multi-branch quadrature (solver_recursive/), NOT a
  representative branch.** The default fork enters `point_map.py`'s banker FOCs
  as a real probability-weighted integral: Gauss-Hermite over the s-innovation ×
  the default realization d′∈{0,1} weighted by π_t (EXOGENOUS input path), where
  the default state is the SAME fitted decision rules evaluated at a reachable
  next-period point — never a frozen stand-in economy. The premium is endogenous
  (Ω^d > Ω^nd on the low default payoffs). `pi ≡ 0` nests the risk-neutral model
  exactly (test_recursive_nesting). The earlier perfect-foresight
  representative-branch pricing got the sign wrong (expansionary); the entire PF
  stack (`solver_pf/`: transition + solvers + risk_branch) was deleted 2026-08-11
  and the Chebyshev-Smolyak projection solver is now the ONLY machinery (TFP is a
  deterministic 7th state Z_D).
- **OMT/TPI = THIRD quadrature branch (solver_recursive/):** the default fork
  splits into backstop-HONOURED and RENEGED, weighted by the priced activation
  probability `cal["tpi_activation"]` (per-experiment scalar, not a state).
  Honoured redeems the D-bond at `recovery_tpi_D` and averts a fraction
  `tpi_real_shield` of the recession (continuation blended toward d′=0); at both
  = 1 the effective default prob is π·(1−a) (OMT removes the premium).
  `tpi_activation = 0` (or `tpi_real_shield = 0`) nests the two-branch solve.
  Monotonic on the financial channels; `tpi_recursive_experiment.py` compares
  0/50/100%. Full activation (a=1, shield=1) trips the μ=1 slack-slip → shield
  defaults to 0.5.
- **Predetermined deposit rate:** the rate paid at t was locked at t−1
  throughout (bank funding legs, household EGM returns, μ timing).
- **Predetermined capital (Bocola eq. 6):** the stock producing at t was
  bought at t−1 (`Kap_prod[t] = Kap[t−1]`); mpk is the marginal product of
  the bank-held vintage, so impact output moves through hours alone. The
  old contemporaneous timing let the sovereign-risk investment boom raise
  Y_0 directly — reverting it re-opens the comovement problem.
- **Union deposit market (deposit-UIP):** deposits are own-good claims at
  national rates; a frictionless union interbank replaces the two national
  clearings with ONE union-wide clearing (D-good units) plus real-rate
  parity (1+rdep_D) = (1+rdep_F)·p′/p — the flexible-price image of one
  nominal union rate + national inflation differentials (BKK/Baxter-Crucini
  single-traded-bond margin). UIP makes the interbank pass-through
  zero-profit → no Walras leak; the cross-border deposit position
  (`out["nfa_dep_D"]`) is the absorption margin that broke the national
  S=I trap (the M1 comovement mechanism). A literal rdep_D=rdep_F with
  own-good legs is WRONG (unassigned RER valuation profit → Walras leak).
  Stage-2 SS imposes β_F = β_D (symmetric-SS doctrine).
- **Hatchondo-Martinez perpetuity:** stock decays at rate 1−δ_b; duration
  ≈ 1/δ_b quarters (0.036 ⇒ ~7y). Long duration is what makes priced risk
  generate large MTM losses — an interlude with δ_b=0.25/recovery=0.80 cut
  the repricing ~6x and made the risk channel expansionary (study §8).
- **Default branch = ONE pure-haircut feared event (Bocola):** the branch
  solves a single deterministic event — a full write-down to recovery
  `recovery_rate_D` = 0.45 (Greek PSI; Bocola's D = 0.55) on the whole
  claim, with the default-state recession arising endogenously through bank
  balance sheets. There are no scarring add-ons: the old Arellano output
  cost, GK ξ_K capital-quality loss and HFSF recap FLAGS were all 0 at the
  Bocola-pure baseline and were deleted in the 2026-07-21 cleanup (see git
  history if a variant needs them back). The recap machinery survives only
  as `_RECAP_LADDER`, a warm-start continuation used when the direct branch
  solve stalls. If the event is infeasible after that, the branch RAISES.
  Bohn taxes respond to the SURVIVING stock (taxing the pre-haircut stock at
  t=0 was a ~31%-of-GDP artifact).
- **Working capital (Neumeyer-Perri):** ζ_wc=1 × wage bill pre-financed at
  r_wc = rdep(−1) + λμ/Ω̃; the wedge is the only channel from spreads into
  impact output (without it Y_D moved −0.2% even at Q_bD −30%, n_D −20%).
  Financing income passes through to household dividends (intra-period,
  never on the bank balance sheet); routing it through bank equity would
  break the closed-form leverage/spread calibration.
- **Walras redundancy:** goods_F and the current account are *dropped* from
  the residual system and monitored as diagnostics.
- **Policy rules present:** the Bohn tax, and the TPI backstop (a Markov-
  switching CB price floor on the D-sovereign, `psi_cb_D`, with an explicit
  CB-budget rebate `rem_cb_D`). No macroprudential policy, by design.

## Known limitations (documented, next thesis phases)

- Comovement problem RESOLVED (2026-07-18): predetermined capital + the
  union deposit market restored the impact contraction at the headline
  shock (Y_D[0] and I_D[0] both negative at π = 1%·0.95^t). Remaining
  next-phase dial: NK/union nominal block (needed for the TPI
  application); real interest parity currently plays the role of the
  single policy rate.
- Risk channel (recursive) approximations: Λ^nd ≡ beta_inter, rep-agent
  income-SDF proxy for Λ^d, household-side π-blindness (the deposit Euler never
  weights the default branch — faithful to Bocola, where household deposits are
  riskless too). Magnitudes are indicative at μ=1 (13-pt grid); μ=2 does not
  converge (the near-unit-root d′=1 corners — see the methods note §8). The
  weak output pass-through at the standard calibration is why the OMT/TPI
  activation comparison cushions the financial channels cleanly but barely
  moves output.

## Branch convention

- `bocola-rewrite` — current working branch (Bocola-faithful trim, 2026-07-16).
- `global` — pre-rewrite snapshot (CK zones, always-binding IC).
- `main` — merge target.
- `audit`, `bank-cal` — historical SSJ-era branches; do not use for new work.

## History

The previous implementation used the `sequence_jacobian` (SSJ) library
(`code/model_v12.ipynb`, `equations_*.py`, `audit_artifacts/`) — superseded
by the standalone `code/global/` model in July 2026. The SSJ-era audit trail
(six structural fixes W-1…TPI-1, Walras forensics) lives in `docs/audit.md`,
`docs/walras_forensics.md`, `docs/verification_report.md` and git history.
`docs/STATE.md` records the current model state and calibration.
