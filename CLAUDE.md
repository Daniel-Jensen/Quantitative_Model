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

The model is solved with global nonlinear methods over 7T stacked unknowns
`[N_D, N_F, Kap_D, Kap_F, rdep_D, rdep_F, p]` under perfect foresight
(MIT shocks), T=200.  Solver (`solvers.py`, 2026-07): damped Newton on an
explicit finite-difference Jacobian — built in parallel (multiprocessing)
and REUSED across warm re-solves via caller-owned `jac_cache` dicts, which
is what makes the CK/risk fixed points fast — with scipy hybr as the
cold-start/fallback path and a Newton polish to `tol_transition` (1e-10).
Do not tighten `tol_transition` toward 1e-12 casually: hybr alone plateaus
near 5e-11 (it stops on xtol, not the residual), and 1e-12 previously made
every solve grind through futile restarts and then raise.  Hot kernels
(household EGM backward, distribution forward) are numba-JITed with an
exact pure-numpy fallback (`cal["use_numba"]`; equivalence is tested).

| File | Contents |
|------|----------|
| `calibration.py` | All parameters. Single λ per bank (Bocola IC); Bocola/Greece anchors documented inline. f = exit/payout share; Ω = β·[f + (1−f)α′] (Bocola's ψ = 1−f survival weight on the franchise value). |
| `steady_state.py` | Two-stage SS solve: {rk_D, rk_F, p} on capital markets + current account, then {β_D, β_F} on deposit markets. Symmetric SS required (see docstring). |
| `bank.py` | GK/Bocola bank block. `bank_backward` (α, μ, bond prices, cross-border FOC holdings), `bank_forward` (net worth, dividends, deposit supply; portfolio shares on ACTUAL net worth). PRICED (`def_price_D`) vs REALIZED (`def_real_D`) default split; only D is risky, F bonds are safe. |
| `government.py` | HM perpetuity bonds, Bohn rule. `govt_transition` forward-integrates the debt stock in one pass. Default risk is exogenous (no crisis zones). |
| `transition.py` | 7T stacked system: 2 FB complementarities, 2 labour, UNION deposit clearing + deposit-UIP, goods-D. IC is OCCASIONALLY BINDING: the capital-market block imposes the Fischer-Burmeister complementarity 0 ≤ μ ⊥ slack ≥ 0 (slack = αn − λ·assets); μ from the capital FOC is valid in both regimes. Capital is PREDETERMINED (Bocola eq. 6): the production stock at t is Kap[t−1] (impact output moves through hours alone). Debt is endogenous inside every residual call; banks clear bonds against the true end-of-period stock (`b_D_D = b_gov_eop − b_D_F`). `make_residual` builds the residual from a picklable spec (shared with Jacobian workers). Supports mid-crisis initial conditions (`init=`) for default branches and policy runs. |
| `solvers.py` | Parallel FD Jacobian + damped Newton with Broyden updates, stall-triggered rebuilds, and cross-solve Jacobian reuse (`jac_cache`). |
| `fast_kernels.py` | numba kernels for EGM backward + distribution forward; exact numpy fallback when numba is absent (`cal["use_numba"]`). |
| `risk_branch.py` | **Bocola risk channel**: representative post-default branch — ONE deterministic solve of the PURE-HAIRCUT feared event (def_real_D[0]=1, recovery 0.45 = Greek PSI). Diagnostic flags off at baseline: `def_output_cost_D`, `def_capital_quality_D`, `recap_share_D`. Two-branch risk inputs for `bank_backward`, `solve_transition_risk` outer loop (base ↔ branch fixed point at exogenous π), and `bond_decomposition` (default comp. + risk premium + liquidity premium, exact identity). |
| `household.py`, `distribution.py` | EGM with GHH utility; stationary distribution and forward iteration. |
| `firms.py`, `capital.py`, `trade.py` | Flexible-price production with the Neumeyer-Perri working-capital wedge (w ÷ (1+ζ·r_wc), the spread→output channel; ζ=0 nests exactly — Bocola §V.C's own open-economy fix), Jermann adjustment costs (+ branch capital-quality hook `quality0`), CES/Armington trade. |
| `main.py` | End-to-end run: SS → TFP IRF → Bocola pass-through experiment (exogenous π path 1%·0.95^t). Full pipeline ≈ 1–2 min (branch cold-start ≈ 30s incl. the recap-share warm-start continuation; later rounds ~1s via `jac_cache`). |
| `tests/` | Regression suite (see below). |

## Running and testing

```bash
cd code/global
python3 main.py                              # full pipeline + figures (~1.5 min)
python3 tests/test_ss_identities.py          # SS theory identities (fast)
python3 tests/test_bank_block.py             # bank FOC/no-arbitrage identities (fast)
python3 tests/test_fast_kernels.py           # numba/numpy kernel equivalence (fast)
python3 tests/test_transition_walras.py      # fixed point + Walras with moving debt (~30s)
python3 tests/test_signs_bocola.py           # sign acceptance criteria (~20s)
python3 tests/test_risk_channel.py           # risk-channel nesting/identity/signs (~2 min)
```

**Acceptance thresholds** (all enforced in tests):
- goods_D (imposed) ≤ 1e−9; goods_F (Walras-redundant diagnostic) ≤ 2e−6 —
  including when the debt stock moves.  (The Newton solver typically lands
  goods_D near 1e−13; acceptance is `tol_transition` = 1e−10 normalized.)
- Zero-shock transition stays at SS to ≤ 1e−5.
- Risk-only shock (exogenous π): Q_bD↓, n_D↓, n_F↓, Y_D[0]↓, C_D[0]↓,
  lending spread↑, b_gov↑, Tax↑ (a positive Y or n response to sovereign
  risk = bug).
- Complementarity on every solved path: μ ≥ 0, slack = αn − λ·assets ≥ 0,
  μ·slack ≈ 0 (`out["mu_min_D/F"]`, `out["slack_min_D/F"]`; transition.py
  warns on violations).  Known open item: risk-on n_D[0] can sit above
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
- **Risk channel (Bocola) = two-branch expectations, not a wedge:** bankers
  discount with the household SDF (Λ = β·u_c′/u_c — Bocola uses log utility,
  NOT Epstein-Zin) and weight a post-default branch by the priced default
  probability π_t (EXOGENOUS input path). The premium is endogenous:
  Ω^d > Ω^nd multiplies the low default-branch payoffs. Approximations
  (documented in risk_branch.py): Λ^nd ≡ beta_inter on the base path, ONE
  representative branch reused across dates, aggregate-income SDF as the HA
  rep-agent proxy. `pi ≡ 0` nests the risk-neutral model exactly —
  regression-tested.
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
  balance sheets. Diagnostic flags, all 0 at baseline: `def_output_cost_D`
  (Arellano), `def_capital_quality_D` (GK ξ_K), `recap_share_D` (HFSF-style
  equity injection financed by branch issuance). If the event is infeasible
  the branch RAISES. Bohn taxes respond to the SURVIVING stock (taxing the
  pre-haircut stock at t=0 was a ~31%-of-GDP artifact).
- **Working capital (Neumeyer-Perri):** ζ_wc=1 × wage bill pre-financed at
  r_wc = rdep(−1) + λμ/Ω̃; the wedge is the only channel from spreads into
  impact output (without it Y_D moved −0.2% even at Q_bD −30%, n_D −20%).
  Financing income passes through to household dividends (intra-period,
  never on the bank balance sheet); routing it through bank equity would
  break the closed-form leverage/spread calibration.
- **Walras redundancy:** goods_F and the current account are *dropped* from
  the residual system and monitored as diagnostics.
- **No macroprudential policy** (by design, current phase). The only policy
  rule is the Bohn tax (the default-branch recap flag is a diagnostic, off
  at baseline).

## Known limitations (documented, next thesis phases)

- Comovement problem RESOLVED (2026-07-18): predetermined capital + the
  union deposit market restored the impact contraction at the headline
  shock (Y_D[0] and I_D[0] both negative at π = 1%·0.95^t). Remaining
  next-phase dial: NK/union nominal block (needed for the TPI
  application); real interest parity currently plays the role of the
  single policy rate.
- Risk channel approximations: single representative default branch,
  Λ^nd ≡ beta_inter, rep-agent income-SDF proxy for Λ^d, household-side
  π-blindness (the deposit Euler never weights the default branch —
  faithful to Bocola, where household deposits are riskless too; see
  risk_branch.py docstring). Validation moment: risk-channel share of the
  lending-spread response vs Bocola's "up to 45%".
- On a deterministic path the occasionally-binding IC regime is a
  perfect-foresight switch, not a distribution over binding states —
  the precautionary motive still comes only through the two-branch
  expectations, not through path uncertainty.

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
