# PROGRESS — Changelog

Reverse-chronological development log for the two-country MU-HANK sovereign-risk
model (Greece 2010–12 / ECB TPI). **Newest first.** Reverse-engineered from the
git history (135 commits, first commit 2026-05-14), `docs/STATE.md`,
`docs/audit.md`, `docs/eba_calibration.md`, and `docs/HANDOFF.md` — the detailed
derivations live in those docs; this file is the timeline.

Convention: `[hash]` is a commit on `main`; dates are commit dates. Doc-sync
policy: every commit touching model/code (`code/**`, `*.py`) gets an entry here.
Enforced by `.claude/hooks/require-docs-before-commit.sh` (Claude Code commits)
and `.githooks/pre-commit` (terminal commits; enable with
`git config core.hooksPath .githooks`).

---

## 2026-08-24 — Paper introduction rewritten; motivation figures added (`gk-structural-foc`)

**No model source changed.** `code/`, the calibration and every solved object are
untouched; this entry covers `Empirics/` and the Overleaf project.

- **`Empirics/motivation_figures.py`** (new). Emits two paper figures with captions
  baked into the PNGs, per the repo's figure convention. Outputs to the gitignored
  `Empirics/outputs/`; the tracked copies live in the Overleaf project at
  `VIVA/figures/`.
  - `fig_greece_motivation` — 2×2: Greek general government debt (€bn), debt/GDP,
    the ten-year yield against the Bund with the spread shaded, and real GDP against
    real investment indexed to 2007Q4.
  - `fig_euro_yield_decoupling` — ITA/ESP/PRT/GRC ten-year yields against the Bund,
    1995–2026, with ECB/OMT/TPI markers.
  - Sources: Eurostat `gov_10q_ggdebt` (S13, `na_item=GD`) and `namq_10_gdp` (B1GQ,
    P51G, CLV10_MEUR, SCA); the FRED yield panel already cached by
    `Empirics/graph_spreads.py`.
  - Palette is the existing Okabe–Ito order, validated for colour-vision separation
    (worst adjacent pair ΔE 11.0 deutan, 15.6 normal). Germany is drawn in ink and
    dashed rather than given a categorical hue, since it is the benchmark.

- **Two facts the figures establish, both checked against the series.** Real
  investment fell 69.7% from 2007Q4 to the 2015Q3 trough against a 27.4% fall in
  GDP — the asymmetry that motivates an intermediary-constraint channel rather than
  a demand or labour-wedge channel. And the March 2012 PSI cut €75bn and 33.7 points
  of debt ratio, with the ratio back at its 2011Q4 level by 2013Q2 and the stock not
  until 2021Q2.

- **Introduction rewritten** (Overleaf `866a537`, `a0aeb1c`). Reframed from a
  counterfactual-history question to three mechanism questions; TPI named as the
  modelled instrument via the three design features the model uses; Bi–Foerster–Traum
  separation expanded to three substantive points (CES portfolio aggregator vs.
  equilibrium retrenchment; their reduced-form financing friction vs. no additive
  component here; their joint fundamental/non-fundamental scope vs. our deliberate
  narrowing). Distributional incidence promoted to the second reported finding.

- **Two errors in the previous draft, corrected.** The "87% of the outstanding
  stock" motivating statistic was a model share of model-issued paper and wrong as an
  empirical claim by roughly five times; replaced with the EBA figures (€54bn against
  €23bn of Core Tier 1, 2.4× equity). And `Delta` was glossed as *pledgeability*,
  which is backwards relative to the proportionality identity — from
  `intermediation_IC_D` it is the **divertable** share, so `Delta = 0.20` makes
  sovereign paper better collateral than capital. **`CLAUDE.md` still carries the old
  wording** ("makes Greek paper worse collateral than capital") and contradicts both
  the code comment and the IC algebra; flagged, not yet changed.

- **`docs/referee_report_2.md`** (new, `8fc0786`) — hostile second-referee report on
  the paper's motivation, ten findings. The introduction has been rewritten against
  it; Sections 2–4, the abstract and the appendices are unaddressed.

---

## 2026-08-19 — Central-bank block audit; `docs/cb_mechanism.md` promoted to canonical (`gk-structural-foc`)

Diagnose-and-report audit of the CB block against the refactored pricing. **No model source
changed.** Evidence in `diagnostics/cb_audit/` (`run_log.md`, `VERDICT.md`,
`recommended_fix.md`, four probes, a self-tested Prony estimator).

- **The CB block survived the refactor untouched and is correct.** The four TPI blocks
  (`domestic_bond_clearing_tpi`, `budget_residual_D/F_tpi`, `external_account_D_tpi`) are
  byte-identical to their pre-refactor versions — they were written against the *payoff* in
  coupon/survival form, never against the old FOC's price decomposition. All mark exclusively
  at endogenous `q_b_D`; `cb_flow_D` matches `bond_return_D` term for term including
  `zeta_writeoff_D` and `writeoff_enabled_D`. Zero live `psi_spread` references in the CB
  block, clearing condition or residual equations.
- **SS neutrality is exact, not approximate.** Every TPI block output is bit-identical
  (difference `0.000e+00`) to its non-TPI counterpart at the steady state, and
  `cb_flow_D = rem_cb_D = rem_cb_F = 0`.
- **The mechanism, stated canonically.** With `psi_lambda_B = 0`, `SDF_banker` and `Omega_p1`
  cancel in `nu_bD_D/nu_K_D`, leaving `rb_exp_D(+1) - rdep_D = 0.20*(rk_D(+1) - rdep_D)`
  exactly. The CB has **no direct lever on the spread**; it compresses only by lowering
  `rk_D`. TPI's spread effect and its investment effect are therefore ONE effect and must not
  be reported as two.
- **The 2x2 sovereign-holdings matrix (new decisive diagnostic, `probe_portfolio.py`).** Both
  clearing identities close to <=1.2e-15 at 26 checkpoints. ~72% of the CB book (84% in pure
  quantity, ~99% from t=1) is bought from **German** banks, 17% is new Greek issuance, only
  10% comes off Greek banks. TPI at gamma=10 undoes 82.9% of the crisis rise in `phi_bD_D`,
  and does it through the DENOMINATOR — numerator moves -0.4%, `n_inter_D` recovers +6.6%.
  German banks' Greek exposure falls to 27% below SS. The ECB buys German banks out; it does
  not share the exposure.
- **The "closed-loop pole at gamma ~ 27.3" is a T=500 terminal-truncation artefact.** The
  resonant eigenvector carries 0.0000 of its mass in the first 100 quarters and 0.9922 in
  t=400-499; `||A_cb[:,499]|| = 3.86` against ~0.0065 for every interior column and
  `A_cb[499,499] = +1.080` is the only positive diagonal in the matrix. Dropping five columns
  removes every pole below gamma=36 and changes the reported peak spread by **nothing** at
  gamma=2/5/10. The condition-number scan in `code/tpi.py` and
  `lottery_math.closed_loop_pole` also steps clean over a nearer apparent singularity at
  gamma=2.2116. No reported number is affected; the guard and the documented claim are wrong.
  Fix proposed (R-1), **not implemented**.
- **Reporting hazard recorded (F-1).** `writeoff_enabled = 0` means no credit loss ever flows
  through the conduit, so the realised German transfer and the printed "F bears EL PV" are
  different objects living in different places. They must never be netted.
- Three `diagnostics/` scripts still execute against the deleted `psi_spread`;
  `solve_configs.py` is the dangerous one — it does not crash, it silently produces a
  `psi_lambda_B = 0` arm identical to its own baseline.

---

## 2026-08-18 — GK structural refactor stages 2–5: no sovereign spread wedge (`gk-structural-foc`)

The sovereign spread is now generated by the bond's state-contingent payoff inside the genuine
GK portfolio FOC. Full detail and numbers in `docs/STATE.md` -> *GK structural refactor*.

- **Root cause was the PAYOFF, not the pricing block.** `zeta_writeoff_D = 0` wrote down only
  the current coupon on default and left the perpetuity's continuation value whole,
  understating the loss on a 12.9-quarter claim by `[delta_b + (1-delta_b)q_b]/delta_b = 12.6x`
  (`EL 0.0561` against the contract's `0.7014`). `psi_spread_D = 0.615` was standing in for
  almost exactly that gap. The stage-1 conclusion that the GK mechanism was too weak to
  generate the spread was wrong — it was being fed the wrong payoff.
- **Deleted, not recalibrated.** `psi_spread_D/F`; `EL_price_D/F` as a pricing wedge;
  `divert_bond_foc_D/F`; `divert_portfolio_adj`; `bond_price_ss_D/F`; `domestic_bond_foc_D/F`;
  `portfolio_adj_cost`; `excess_return_bD_D_ss`, `excess_return_bF_F_ss`,
  `excess_return_F_D_ss`, `excess_return_D_F_ss`.
- **`bond_return_D/F` is the single source of truth** for the payoff and emits three things:
  `rb_exp` (expected — the only return the pricing equations read), `rb_actual` (realised
  branch, still gated by `writeoff_enabled = 0`, so S-1's pure risk-premium framing stands),
  and `EL_load` (diagnostic; read by `code/tpi.py`'s CB P&L and nothing else).
  `zeta_writeoff_D/F = 1`.
- **`gk_bond_foc_D/F`** impose `nu_own = Delta_own_eff * nu_K`, which with
  `intermediation_P1_D/F` is `rb_exp(+1) - rdep = Delta_eff*(rk(+1) - rdep)`. `q_b_D`/`q_b_F`
  became SS UNKNOWNS with `rb_D_res`/`rb_F_res` as their targets.
- **`gk_cross_border_foc`** (in `equations_global.py`) states the same FOC on the two
  cross-border legs plus the `psi_bF_D`/`psi_bD_F` stock cost, divided through by
  `SDF_banker*Omega_p1` so those keep their calibrated units. This also removed a genuine
  DOUBLE COUNT: from the 2026-08-17 draft the expected loss sat in both `intermediation_P1`
  and `divert_portfolio_adj`, so the cross-border legs netted it twice while the own legs
  netted it once. `Omega_p1_D/F` are now exported from `P1` rather than duplicated.
- **`psi_lambda_B_D/F = 0` is the preferred baseline.** No independent Greek observable
  identifies a sovereign-specific haircut *elasticity*. `3.01` is retained as a diagnostic arm
  only.
- **`Delta_bF_D`, `Delta_bD_F`: 0.40 -> 0.20, forced not fitted.** With the own legs pinning
  both bond prices, and `rk_D = rk_F`, `rdep_D = rdep_F = 0` at a riskless SS, the cross-border
  ratio `nu_cross/nu_K` is no longer free. Holding 0.40 leaves a constant 80bp/yr cross-border
  wedge. Measured wedges at the live calibration: `8e-11` / `-9e-11` bp/yr.
- **Pre-existing units bug fixed in `intermediation_P1_F`.** `q_b_D` and `q_b_F` are both
  D-good prices, but `P1_F` compared an unconverted `rb_actual_F(+1)` with the F-good `rdep_F`
  while the old cross-border block applied `p/p(+1)` to the same return. The conversion now
  lives once, in `P1_F`. SS-neutral (`p` constant at SS); first-order relevant off it.
- **Depreciation calibration is now ITERATED** to a fixed point (5 passes, `5e-14`). One pass
  left `rk_D = 0.009981` once `q_b` joined the SS unknowns and the `delta -> K -> rk -> q_b -> K`
  loop closed. `rk_D = rk_F = 0.010000` exactly, as RK-1 requires.
- **New guard `steady_state.report_gk_steady_state`** prints the §13 diagnostic table and
  RAISES if any of the four portfolio FOCs is violated or any `Delta_*_eff` leaves `[0,1]`.
  Runs on every solved SS. Two new tests in `code/test_nkpc_blocks.py`: an AST scan that fails
  if any deleted name reappears in live `code/*.py`, and a structural double-counting check
  that no portfolio condition takes both a `nu` and a `def_rate`/`EL` object. 38 fast tests pass.
- **Results, 1pp shock (nothing tuned to any moment).** Peak spread **205.9 bp**; German yield
  **-16.2 bp** on impact (endogenous flight to quality, no F-side wedge); `b_DD` **+2.12%**
  while `K_D` **-0.067%** (the intended balance-sheet crowding-out); `b_DF` **-2.10%** — German
  banks now DO retrench, reversing the open problem in HANDOFF.md; `n_inter_D` **-11.41%**,
  `Y_D` **-1.974%**, `C_D` **-2.511%**. `b_gov_D[499] = 3.0e-05`.
- **TPI loading is 0.520 / 0.504 / 0.482 at gamma = 2/5/10 — BELOW 1.** The old 3.82/2.90
  over-compensation headline is reversed, entirely through the denominator: the CB earns the
  same premium but absorbs a 12.6x larger expected loss. The decline in gamma (the
  self-extinguishing premium) survives. E3's 2026-08-06 `zeta_writeoff = 1` finding is no
  longer a robustness variant; it is the baseline.
- **Documentation fix, not a regression:** `goods_mkt_D` is `-4.23e-07`, and a clean worktree
  at `91ac778` prints `-4.2493e-07` on the same pipeline. CLAUDE.md's `<= 1e-14` threshold for
  that residual was never met on this calibration and has been corrected to `1e-6`.
  `ca_res_D` (`1.7e-16`) is the residual that genuinely reaches machine zero.
- **`experiments/` rebased.** `fig04_spread_decomposition` no longer draws an
  `EL_price`/`psi_spread` share bar — that framing is forbidden in a linearised model and the
  split was mostly a calibration artefact. It now plots the Greek yield under direct
  expected-loss pricing (required return frozen at SS, bond price solved forward on the
  model's own recursion) against the equilibrium Greek yield and the GR-DE spread.
  **Finding: the intermediary channel is a QUANTITY amplifier, not a price amplifier** --
  direct pricing alone gives 213.8bp against an equilibrium 189.7bp Greek yield, and the
  spread exceeds the Greek yield only because the German leg falls 16.2bp. A first draft
  plotted `EL_load_D*def_rate_D` (a one-period capital-loss rate) against `spread_rb` (a
  coupon-equivalent yield); not commensurate, and it produced a spurious 0.62
  "amplification factor". E3's variants
  are rebased on `zeta = 1` (`e3a_realised_writeoff`, `e3b_coupon_only_pricing`). Cache key
  `EL_price_D` -> `EL_load_D` in `diagnostics/regimes/regime_model.py` and
  `experiments/e1_backstop_schedule.py`.
- **The "aggressive" named regime is no longer 50% compression — it is 40.3%.** The closed
  loop has a POLE at `gamma ~ 27.3`; max compression below it is 46.6% and the 50% target is
  met only on the far branch. New `lottery_math.closed_loop_pole` locates it by CONDITION
  NUMBER (a 61-point monotonicity scan of [0,40] steps straight over a pole this narrow and
  reports a spurious non-monotonicity — exactly how the first regeneration attempt failed) and
  `CompressionInfeasible` is raised. **The fallback is 0.75 x pole, not 0.98**: measured, the
  loading schedule is monotone in gamma only up to ~0.85 x pole, and at 0.98 x pole the
  discounted consumption gains hit +11..+12.4% of SS consumption and Greek output goes +1.15%
  on impact — the singularity, not the policy. New shared constant
  `lottery_math.POLE_SAFETY_FRACTION = 0.75`; `common.named_regime_gammas` falls back to
  `gamma = 19.875` (40.3%), `e1.loading_schedule` and `code/tpi.py`'s effectiveness curve cap
  their grids the same way. `medium` is unaffected at `gamma = 9.989` (25.0%).
  `experiments/paper_outputs.py`'s second copy of the gamma solve now routes through
  `common.named_regime_gammas` — one definition, not two. **Paper prose calling the aggressive
  regime "50% compression" must be corrected to 40.3%.**
- **`code/tpi_plots.py` figure-8 text was stale in four places and is now derived.** It
  hard-coded "rho=0.8" (rho_def has been 0.9408 since 2026-08-06), a panel title
  "Premium Peaks (gamma~26)" that no longer describes the shape, annotations written for
  a schedule that STARTS above ell=1 ("timid intervention -> high loading (SMP-type)"),
  and a caption asserting the loading falls "toward the fair-insurance limit ell=1" while
  printing "from 0.5x to 0.5x". All four now read off the data; the caption states plainly
  that the schedule sits BELOW 1 throughout. The 'actuarially fair' label also moved off
  the subplot title it was overprinting.
- **New impact-sign table on BOTH shocks in `build_and_solve`.** Added because
  `fig_irf_overview_macro.png` looks like Y_D collapses under TFP; it does not — that is
  the default shock's line. Measured: +1% TFP gives Y_D **-0.073%** (flat) with N_D -3.71%,
  w_D -7.26%, I_D +5.00% — the standard sticky-price contractionary-technology result
  under a phi_pi -> infinity normalisation and GHH labour supply. Not a defect; the guard
  fires only below -0.5%.
- **Two figure defects the pole caused, both fixed.** `fig02` plotted the loading spiking to
  1.17 and collapsing to 0.38 across two grid points, and its caption's own two-branch test —
  written for the old world where the loading STARTED above 1 — read that artefact as
  "crossing below the actuarially fair benchmark of 1". `fig05` showed German exposure
  plunging to -50% of `Y_D`. The caption now has a third branch and states plainly that the
  loading stays BELOW 1 throughout, so over-compensation must not be asserted.
- **STALE:** every E1–E4 artefact and `experiments/paper/fig0*.png`. The SS moved
  (`q_b_D 0.968941 -> 0.974906`); rebuild the regime cache before the experiments, not after.
- Old `diagnostics/psilam_*` and `diagnostics/substitution_v2/` scripts still reference
  `psi_spread`/`EL_price` and will fail. They are one-off historical probes of the deleted
  specification; superseded, not ported.

## 2026-08-17 — GK structural refactor stage 1: bounded pledgeability (`gk-structural-foc`)

Audit of the sovereign-risk-to-bank-financing block, then the first of five stages. Full
detail in `docs/STATE.md` -> *GK structural refactor*.

- **Audit finding.** The chain `p_def -> Delta_bD_eff -> IC -> lambda_gk/Omega -> P1 -> q_b_D`
  is broken at the third arrow. `Delta_bD_eff` moves only `theta_D`; the Greek spread comes
  entirely from `divert_bond_foc_D`, which touches no endogenous GK object and carries the
  frozen `psi_spread_D`. `bond_price_ss_D`, `steady_auxilliary_D`, `smart_steady_D` are all
  SS-only and absent from `build_block_list()`.
- **GK portfolio optimality is violated at the SS.** `nu_bD_D/nu_K_D = 0.2491` against
  `Delta_bD_D = 0.20`, and `nu_bD_D == nu_bF_D` bit-identically (0.02696043) while
  `Delta_bD_D = 0.20` vs `Delta_bF_D = 0.40`. `steady_auxilliary_D` defines the marginal
  values from returns and never restricts them; the portfolio FOCs are imposed nowhere. This
  is *why* the wedges exist.
- **Stage 1 (this commit).** New `collateral_quality_D/F` export the four `Delta_*_eff_*`
  under a bounded map `Delta + (1-Delta)*z/(1+z)`, `z = psi_lambda_B*def_rate(+1)/(1-Delta)`.
  Local slope is `psi_lambda_B` exactly (SSJ Jacobian: 3.0100000000), so IRFs are unchanged;
  range `[Delta,1)` closes the domain hole at `def_rate(+1) > 0.266` where the old linear form
  drove `1-Delta_eff` negative.
- **New SSJ gotcha recorded.** `np.exp` in a `@simple` block raises
  `TypeError: ... AccumulatedDerivative`. Simple blocks differentiate through a dual-number
  type supporting arithmetic operators only. Hence the rational rather than exponential
  saturation.
- **Doc drift corrected.** `docs/eba_calibration.md` ledger said `Delta_own` committed at 0.80
  and cross at 0.90; CLAUDE.md's GK-1 row said `Delta=0.85/0.90 -> lambda_gk_D=+0.927`. Live
  values are **0.20 / 0.40** with `lambda_gk_D = 2.2129`. Both docs described the CT1-scope
  world; GK-2's broad scope cut `phi_own` 2.39 -> 0.456, which satisfies feasibility at the
  inherited `Delta`, so the raise was never adopted. ~4x error for anyone computing the
  collateral channel from the old numbers.
- **Verification.** `code/main.py` exit 0, bit-identical: `n_inter_D[0] = -6.7366%`,
  `Y_D[0] = -0.8521%`, peak spread +0.375 pp, `goods_mkt_D = -4.2493163257550925e-07`,
  `max abs(goods_mkt_F)` 2.06e-10..2.12e-10 across the gamma grid. 35 fast tests pass.
- **Rejected en route.** A prior `writeoff-test` branch flipped `zeta_writeoff`/
  `writeoff_enabled` to 1 to test whether realising the default loss cures the `n_inter_D`
  overshoot (+3.14% at t=8). It does not: everything scales 3-4x (peak spread 538.5 bp,
  `Y_D[0]` -3.68%) with timing untouched — `n_inter_D` still turns positive at t=4 and peaks
  *higher*, +5.53%. Relative overshoot halves (peak/|trough| 0.465 -> 0.234), so it bites on
  the right margin but nowhere near enough. Branch deleted, S-1 stands.

---

## 2026-08-07 — Fiscal rule and fiscal limit audited; `Empirics/fiscal_limit.py` added

No model changes. Two existing mechanisms audited to see whether the fiscal block could
generate foreign retrenchment without new wedges or shocks. It cannot, and the audit
changed what is claimable about both parameters. Full detail in `docs/STATE.md` ->
*Fiscal rule and fiscal limit: what is identified*.

- **`phi_lamb_D` sweep {0.15, 0.10, 0.07, 0.05}.** Stability boundary is between **0.10
  and 0.07**, not the 0.05 the pipeline's printed `rho_b` gate predicts — that gate is
  partial-equilibrium and omits `def_scale_D`, so it is optimistic enough to land a user
  on a divergent calibration. Retrenchment (`b_D_F` < 0) appears ONLY in the divergent
  region, so it is an artefact: **no stationary calibration of the Bohn rule produces
  retrenchment.** `Y_D[0]` is insensitive across the stable range (-0.852% to -0.871%).
  `phi_lamb_D` = 0.15 keeps ~1.5x margin over the true floor and stays.
- **`mv_rule_D` = 0 justified.** Par and market-value debt gaps move in opposite
  directions in a crisis (par positive 39/40 quarters, market-value negative 40/40, the
  latter 2.88x larger); the market-value rule would CUT taxes 0.75% of quarterly GDP at
  impact, reading a wider spread as a windfall. Maastricht debt is nominal face value, so
  the par rule is institutionally correct. The "market-value rule REQUIRED" comment was
  stale from the CT1 scope (`phi_bD_D` = 2.39 vs 0.456 under the broad scope) and is retired.
- **Fiscal limit estimated** (`Empirics/fiscal_limit.py`, new): BFT's logistic on Eurostat
  `gov_10q_ggdebt` plus the repo's Greek-Bund spreads. Preferred pre-OMT sample
  `eta0 = -14.80 (0.55)`, `eta_s = 7.67 (0.47)`, R2 = 0.849, n = 50 — same family as BFT's
  Italian `-10.70 / 5.25`. Post-2012 data must be excluded because OMT severed the
  debt-spread link (debt 152% -> 181%, spread 13.4pp -> 6.3pp -> 1.4pp), which would build
  the studied policy into the parameter.
- **`def_scale_D` stays at 0.25**, now with provenance: it sits inside the estimated range
  (0.04 full / 0.19 crisis / 0.63 pre-OMT). The best-fitting 0.633 is unusable — at that
  value `psi_lambda_B` is not continuously calibratable (divergences at 2.00 and 2.34
  bracketing a marginal island at 2.10-2.20; the 150bp target sits within 0.09 of a
  blow-up). Curvature is qualitatively wrong (0.5 concave vs 3.9-10.9 convex estimated)
  but is not identified at first order.
- **Benchmark check:** BFT set `phi_T` = 3 *"to ensure stability of the debt path"* — they
  do not calibrate their fiscal rule either. What they estimate is the fiscal limit. Their
  rule can be weaker because default is realised (`Delta_t = delta_b` writes debt down) and
  their default probability is logistic/bounded; with `writeoff_enabled_D` = 0 the tax rule
  is this model's only stabilising device.

---

## 2026-08-07 — Country-size asymmetry: F is 11.7x D (`fix-cross-border-units`)

**The defect.** The model normalised `Y_D_ss = Y_F_ss = 1` — Greece and Germany the same
size — while every EBA moment is a ratio to its **own** country's net worth. Cross-border
stocks built as `phi * n_holder / q` therefore landed in the **holder's** units. The model
could match the portfolio-composition moment (`phi_bD_F = 0.0075`, DE banks' Greek book /
DE bank net worth) **or** the market-structure moment (foreigners hold 12.72% of the
bank-held Greek stock), never both: joint consistency needs `n_F/n_D = 8.85` against the
model's 0.761, and the gap is exactly the Germany/Greece GDP ratio. Matching composition,
as the model did, put the foreign share at **1.25%** against 12.72% in the data, and
symmetrically overstated Greek banks' share of the bank-held Bund stock at 1.50% vs 0.13%.

**The fix.** `size_F = 11.697` (Eurostat 2010 annual GDP, `data/eba_moments.json`
`raw_EURm`, exposed by `calibration.load_eba_size_ratio`). Convention: **every F variable
is per F capita and O(1); every D variable is a D aggregate** (`size_D == 1`). The weight
appears in exactly the blocks where the two countries meet — `trade_balance`,
`external_account_D`, `global_goods_mkt`, `domestic_bond_clearing`, plus the three `_tpi`
overrides. Nothing inside the F household, bank or production blocks changes, and no grid
is rescaled.

**Home bias split.** A single shared `omega` is inconsistent with size asymmetry: at
`omega_F = omega_D` the larger country's imports from the smaller come out `size_F` times
too large. Symmetric bilateral trade intensity pins the pair,
`size_F*(1-omega_F) = (1-omega_D)`, giving `omega_D = 0.85` (unchanged) and
`omega_F = 0.98717`.

**Bug this exposed — every TPI result before today carried it.** `budget_residual_F_tpi`
paid `rem_cb_F = kappa_cb_F * cb_flow_D / p`, a **D-aggregate** ECB cash flow, into a
**per-F-capita** budget. At equal country size the missing weight was exactly 1.0, so it
was invisible: `goods_mkt_F` sat at 2e-10 for the whole history of the block. Under
`size_F` it leaked up to **1.98e-2 of F GDP at gamma=10** while gamma=0 stayed clean —
the signature of a conduit-only units error. Fixed by `/ size_F`; `max|goods_mkt_F|` is
now 2.06e-10..2.12e-10 across the whole gamma grid.

**Recalibration.** `size_F` makes a Greek shock a much smaller shock to F, damping
cross-border amplification: peak spread fell to 145.20 bp at the incumbent
`psi_lambda_B = 2.92`. Re-bisected on the same 150.14 bp moment (SS re-solved per point):
`2.92 -> 145.20`, **`3.01 -> 149.93` adopted**, local slope 52.6 bp/unit.

**Verification** (`code/main.py`, exit 0): foreign shares **0.1274** (EBA 0.1272) and
**0.001298** (EBA 0.001301) — both moments now hold jointly; `phi_bD_F` exact;
`K_D = 10.800`, `K_F = 10.824` against the 10.8 over-identifying check;
`goods_mkt_D/F ~ 4.2e-07`; `ca_res_D = -2.8e-17`; IC residuals machine-zero; GK
well-posed; `n_inter_D[0] = -6.7366%`, `Y_D[0] = -0.8521%` (both correct sign);
peak spread 150.0 bp. 48/48 fast tests pass, including a new
`code/test_cross_border_units.py` that locks the per-capita/aggregate convention and
asserts the direction of the weight so it cannot be silently inverted.

**Reported outputs that moved.** TPI loading schedule 5.60 / 5.43 / 5.18 at
gamma = 2/5/10. Spread compression at gamma=10 is 22.4%. German-side responses are now
an order of magnitude smaller, which is the point: `rdep_F` falls 2.5 bp on the default
shock where it fell 17.1 bp before, and `n_inter_F` rises 0.14% where it rose 1.20%.
Cross-border absorption correspondingly matters far more for Greece — freezing `b_D_F`
now costs `Y_D` -0.83% -> -1.32% and `C_D` -0.77% -> -1.58%.

**Superseded within the same branch.** An interim patch first scaled the two cross-border
stocks directly into issuer units. That matched market structure but broke composition
(`phi_bD_F` 0.0075 -> 0.086) and pushed `K_F` to 10.672 — it traded one moment for the
other rather than satisfying both. Replaced by `size_F`.

---

## 2026-08-06 — Regeneration on the MS-disciplined shock

- Rebuilt regime cache, E1-E4 and all paper figures at `rho_def=0.9408`, `psi_lambda_B=2.92`. E2 closes at 1.1e-16.
- **Live Claim 5 weakened:** loading schedule 4.43->1.49 becomes 5.65->4.59. Still monotone, but the premium no longer approaches extinction. Live Claim 1 correspondingly stronger (floor 1.49 -> 4.59).

## 2026-08-06 — `rho_def` promoted to the calibration and disciplined at 0.9408; `psi_lambda_B` re-tuned 7.85 → 2.92 (`add-nkpc`)

**Problem.** The sovereign-risk shock's persistence was **hardcoded at `rho_def = 0.80` in
`code/full_model.py:220`**, next to `rho_Z_D = 0.8`. It was therefore neither stated nor
defended as a calibration choice, and it implied a **14-month** crisis
(`0.80^(1/3) = 0.9283` monthly → 13.95 months). The repo's own estimation contradicts that.

**The estimate.** `Empirics/outputs/ms_regime_GRC.npz` fits three Markov-switching states to
monthly Greek–Bund spreads (348 obs, 1997-06 to 2026-06). The crisis state has mean spread
**9.63pp** and monthly persistence **0.9798499**, an expected duration of **49.6 months**; the
realised episode ran 2010-04 to 2017-12, **92 months**. Quarterly equivalent
`0.9798499^3 = 0.94076` → **`rho_def = 0.9408`** (16.9 quarters).

**Changes.**
- `code/full_model.py` — `rho_Z_D` and `rho_def_D` now read from `calibration_start` with the
  old literals as a fallback, and the resolved values are printed.
- `code/calibration.py` — new *Shock processes* block: `rho_def_D/F = 0.9408`,
  `rho_Z_D/F = 0.80`. **`rho_Z` deliberately unchanged** — the MS estimate is about sovereign
  spreads, not TFP.
- `code/calibration.py` — `psi_lambda_B_D/F` **7.85 → 2.92** (`EBA_CALIBRATION` branch only;
  the `else 3.0` branch is untouched).

**Re-tune.** Peak spread is monotone increasing in `psi_lambda_B`; at the new persistence the
old 7.85 gave **470.62 bp** against the 150bp GR–DE moment. Full pipeline re-solve per point:
7.85 → 470.62, 2.73 → 139.60, 2.8909 → 148.50, 2.9181 → 149.99, **2.92 → 150.09 (adopted)**.
Harness sanity anchor at `7.85 / rho=0.80` reproduced the recorded baseline bit-for-bit
(150.14bp, `Y_D[0] = −0.5064`, `C_D[0] = −0.5103`).

**Method finding worth keeping.** Sweeping `psi_lambda_B` by patching it and `psi_spread` onto
an already-solved SS and re-solving only the Jacobian is **wrong**, even though the SS really
is `psi_lambda_B`-neutral (bit-identical `goods_mkt_D`/`K_D`/`beta_D` at every value). It
predicted 150.33bp at `psi = 2.73` where the pipeline gives 139.60 — the
`intermediation_IC_D` `Delta_bD_eff` collateral channel does not pick the patch up, only
`divert_bond_foc_D`'s `psi_spread` does. The first bisection was discarded and redone.

**Results.** Like-for-like at 150bp: `Y_D[0]` −0.5064 → **−0.7502**, `C_D[0]` −0.5103 →
**−0.7014**, `I_D[0]` −1.0114 → **−1.7107**, `n_inter_D[0]` −4.2962 → **−6.2710**. Cumulative
40-quarter `Y` −0.0492 → **−2.5420** (51.7×); negative-`Y` quarters in the first 40: 5 → **37**;
spread above half-peak q3 → **q11**. All four impact signs stay negative, so the `add-nkpc`
consumption sign flip survives.

**Paper-level consequence.** `psi_spread_D` is linear in the dial, so it falls 1.604839 →
**0.596959** against an unchanged `EL_price_D = 0.056134`. The default-loading split moves
**3.38% / 96.62% → 8.60% / 91.40%** fundamental / collateral friction — a friction:fundamental
ratio of **10.63:1**, down from 28.59:1. The constrained-seller claim survives in direction but
"essentially all of it" must become "roughly nine tenths of it".
`experiments/paper_outputs.py`'s `fig04_spread_decomposition` prose needs re-deriving again.

**Issue I-1 substantially resolved** — and by the shock process, not by a capital friction, so
the earlier conclusion that `chi1`/`omega_I` were the wrong hypotheses holds. Residual defect:
`Y_D` still blips marginally positive at q2–q4 (all under +0.03% of SS) before going negative
from q5 and staying there.

**Stability — passes, and moves away from the risk.** SS bit-identical
(`goods_mkt_D = -4.2493506589857954e-07`, `IC_D = 1.776357e-15`, `All residuals < 1e-8 ✓`);
`b_gov_D[499]` on the default shock **fell** 4.63e−05 → 2.04e−05; `ρ_b = 0.8451 < 0.95`; no
`assert_gk_well_posed` failure; all four TPI gammas converge (`max|ca_res_D| ≤ 6.39e−08`).
The re-tune lowers `psi_lambda_B`, i.e. away from the high-`psi_lambda_B` breakdown region.
TPI loading 5.55 / 5.37 / 5.13 at γ = 2/5/10 — monotone decreasing and above 1, so Live
Claims 1 and 5 both survive.

**Tests:** 42 passed (`code/test_nkpc_blocks.py code/test_eba_calibration.py experiments/`),
unchanged from HEAD.

**STALE:** E1–E4, `docs/experiments_results.md`, `docs/paper_draft_results.md` and the eight
tracked `experiments/paper/fig0*.png` all predate this calibration. Regenerate in order:
`diagnostics/regimes/regime_model.py --force` → `experiments/run_all.py` →
`experiments/e4_distribution.py` → `experiments/paper_outputs.py`.

---

## 2026-08-06 — Investment-flow adjustment cost `S(I/I(-1))`, added inactive at `omega_I = 0` (`add-nkpc`)

**Problem.** On the default shock output falls for exactly **one quarter** and then turns
positive: `Y_D` = −0.5064, −0.0026, **+0.0929**, +0.0829, +0.0548, … The comparable published
model (Bi–Foerster–Traum) keeps output negative for ~20 quarters. Investment is the driver —
`I_D` = −1.0114, −0.2671, then a sustained **boom** peaking +0.3324 at q5 that drags `Y` up
with it. The diagnosis: the model had **no adjustment cost on the flow of investment**.
Nothing penalised `I/I(-1)`; the only capital friction was the Q-based cost on the `I/K`
ratio, which lets investment jump down and snap straight back.

**What was added.** `capital_adj_D/F` now carry `S(x) = (omega_I/2)(x-1)^2` with `x = I/I(-1)`.
Effective investment is `I_eff = (1-S)*I` and the installation technology runs on `I_eff`, so
the investment FOC becomes

```
1 = Q*mpi*[(1-S) - S'*(I/I(-1))] + beta*Q(+1)*mpi(+1)*S'(+1)*(I(+1)/I)^2
```

with `mpi = gamma0*(1-ksi)*iota^(-ksi)`. New calibration entry `omega_I_D/F`, **committed at
0.0**, so the committed model is provably unchanged.

**Why it is exactly SS-neutral.** `S(1) = S'(1) = 0`. At the steady state the FOC collapses to
`Q*mpi = 1`, which is the old `q_res = Q - 1/mpi` — the *same root*. The two forms differ by
the factor `mpi`, and since the old residual is zero at SS, `dr = mpi_ss * dr̃` exactly: a
constant row scaling of the target system, which `-H_U^{-1} H_Z` is invariant to. So the
linearised solution is invariant too, not merely the steady state.

**Discounted at `beta`, not the SDF — and this is not an approximation.** `S'(1) = 0` means
the intertemporal term multiplies a factor that is *zero at SS*, so linearising it uses only
`SDF_ss = beta`; the SDF's own deviation contributes nothing to first order. This is the
identical argument `price_nkpc_D/F` already uses (`pi_ss = 0` there). It is also **required**:
taking `SDF_D` makes SSJ's topological sort fail outright with
`hh_D -> capital_fund_D -> capital_adj_D -> sdf_D -> ghh_composite_D -> hh_D`. Locked by an
assertion in `test_flow_adjustment_cost_vanishes_at_steady_state` that `SDF_*` is **not** an
input and `beta_*` is.

**Equivalence gate.** `omega_I = 0` reproduces the pre-change model to **1.08e-13** worst
relative deviation across all 45 dumped arrays (`dump_irfs.py` run at `231327c` immediately
before the edit). Note the older `/tmp/nkpc_irfs_nominal.npz` is **stale** — it predates the
`psi_lambda_B` 8.5 → 7.85 re-tune and differs by 1.56; do not use it as a reference.

**Sweep — the hypothesis is NOT supported.** SS bit-identical at every value
(`K_D = 10.8000000000`, `beta_D = 0.999534992056`), confirming SS-neutrality.

| `omega_I` | `Y_D[0]` % | `Y_D` trough % | contiguous neg. quarters | cum. `Y_D` (40q) | `I_D[0]` % | `I_D` peak boom % | peak spread |
|---|---|---|---|---|---|---|---|
| **0** | −0.5064 | −0.5064 | 2 | −0.0492 | −1.0114 | +0.3324 | 150.1 bp |
| 2 | −0.0287 | −0.0483 | 3 | **+0.2097** | −0.3786 | +0.2573 | 163.7 bp |
| 5 | **+0.0486** | −0.0078 | **0** | +0.3001 | −0.2201 | +0.2082 | 167.2 bp |
| 10 | **+0.0854** | −0.0015 | **0** | +0.3697 | −0.1290 | +0.1748 | 168.4 bp |

The cost does smooth investment — the impact drop shrinks monotonically from −1.01 to −0.13
and the q5 boom from +0.33 to +0.17 — but it **does not convert the V into a sustained U**.
It shrinks the whole contraction toward zero. `omega_I = 2` buys one extra negative quarter
(3 vs 2) at the price of an impact trough **18× shallower** and a cumulative 40-quarter `Y`
response that flips **positive**. At `omega_I >= 5` `Y_D[0]` itself goes positive, tripping
the sign check in CLAUDE.md's *Typical iteration* step 4.

**Mechanism, and why it echoes the `chi1` result.** Making investment sluggish frees the
household budget rather than the economy's resources: `C_D[0]` moves from −0.5103 (at 0) to
+0.1092 (at 2) to +0.2276 (at 10). This is the *same failure mode* as the earlier rejected
`chi1` diagnostic — penalising a capital/investment margin just shifts the burden between `I`
and `C` instead of deepening the aggregate contraction. The persistence problem is therefore
**not** a missing investment friction, and the next hypothesis should look elsewhere (the
`n_inter` rebound at +3.6% by q5 is the more likely engine).

**Left at `omega_I = 0`** pending an author decision. Peak spread drifting 150.1 → 163–168 bp
off the 150 bp target is a second reason not to adopt a positive value without re-tuning
`psi_lambda_B`.

Tests: `code/test_nkpc_blocks.py` **19 passed**; full suite **50 passed**.

---

## 2026-08-06 — Paper figure captions derive from results (`add-nkpc`, Task 17)

**Why.** `experiments/paper_outputs.py` carried a module-level `CAPTIONS` dict of literal
prose written against the flexible-price model. The sticky-price conversion and the
`psi_lambda_B` 8.5 → 7.85 re-tune (Tasks 1–16) left every caption stale and three of them
*inverted*, and because captions are baked into the PNGs, the repo was shipping eight tracked
figures and a generated `docs/paper_draft_results.md` whose prose contradicted its own tables.
The clearest case: `fig08_deciles` claimed the lowest quintile "gains 0.95%" and the highest
loses 0.59%, against a Table 4 *in the same file* reading **+0.4250** and **−0.9073**. This is
the identical hazard Task 15 fixed inside `run_all.py`.

**Fix (structural, not a substitution).** `CAPTIONS` is now empty at import and filled at run
time. `save(fig, name, caption)` takes the caption as a required argument and registers it;
each figure builds it from the arrays it just plotted, via a new `_caption_figNN` helper.
Directional claims are selected from the data by `_monotone`, `_first_quarter` and sign
tests, so a flip rewrites the sentence rather than lying in it — e.g. `fig02` will print "does
NOT fall — the self-extinguishing-premium claim fails at this calibration" if the loading
schedule ever stops declining, and `fig05` will refuse the German-ledger reading if exposure
and loading stop moving in opposite directions. `main()` gained a prose-vs-table assertion:
fig01's caption and Table 3 must agree on impact bank net worth via their two independent
routes (cache vs `e1.run()` payload).

**What the captions were wrong about.**

| figure | was | now (derived) |
|---|---|---|
| `fig01` | "net worth 3.4%", "investment −0.77%" | −4.3% / −1.0%; adds the reversal quarters (net worth q5, spread q8) |
| `fig02` | "4.5× … 2.1×" | 4.43× at γ=0.51 → 1.49× at γ=30, monotone, above 1 throughout |
| `fig03` | "each roughly four times the headline" | **inverted**: consumption carries 0.99× the headline, investment +0.25×, NX −0.21× |
| `fig04` | "3% / 97%" | **3.4% / 96.6%**, re-derived at `psi_lambda_B = 7.85` |
| `fig05` | qualitative only | endpoints: exposure 0 → 0.92% of quarterly `Y_D`, loading 4.43× → 1.49× |
| `fig06` | net path smaller "at every horizon" | **false in the impact quarter**; true in 14 of the first 16 |
| `fig07` | 23/52/25 hardcoded | read from the npz (22.9/52.4/24.7), hawk span 2010–2014 derived; genuinely model-independent |
| `fig08` | "lowest gains 0.95%, highest 0.59%"; "consumption rises on impact" | Q1 +0.4250 / Q5 −0.9073, backstop gain +2.01 / +1.34; consumption **falls** ~0.51% in every quintile on impact |

**`fig04` derivation.** Loading per unit of default probability = `EL_price_D + psi_spread_D`
(bond-pricing FOC, `equations_D.py:566`). `EL_price_D = (1−0.30)·0.0777006/0.968941 =
0.056134`, invariant to `psi_lambda_B`. `psi_spread_D = lambda_gk_D·psi_lambda_B_D /
(beta_inter_D·Omega_D)` (`steady_state.py:104`) is linear in `psi_lambda_B`, so 8.5 → 7.85
took it 1.737724 → 1.604839. Split 0.056134/1.660973 = **3.4% fundamental / 96.6% friction**.

**Verification.** 8 figures + `docs/paper_draft_results.md` regenerated (~95s, no Jacobian
re-solve — runs off the existing regime cache). Every caption cross-checked against the
corresponding table: no contradictions. `pytest code/test_nkpc_blocks.py
code/test_eba_calibration.py experiments/` → **40 passed**. No model, calibration or equation
change; results tables are numerically identical to the pre-Task-17 run.

---

## 2026-08-05/06 — Nominal rigidities: sticky prices + nominal deposit contracts become the baseline (`add-nkpc`, Tasks 1–16)

*One entry for the whole workstream (16 commits, `2015edd`…`120dcf6` plus this doc pass).
Per-task detail is in the commits; the consolidated state is `docs/STATE.md`.*

**Why.** The flexible-price model's response to a 1pp default shock was
`Y_D[0] = −0.0149%` and `C_D[0] = +0.2164%` — two orders of magnitude below
Bi-Foerster-Traum's −0.6% and with consumption *rising* in a crisis. With flexible labour
supply and competitive labour demand, `Y` drops out of the labour block entirely and `N` is
pinned by `Z`, `K`, `P_CES` alone: there was nothing for aggregate demand to act on.

**Structural changes** (`code/equations_D.py`, `equations_F.py`, `equations_global.py`):

- **Task 1** — extracted `full_model.build_block_list()`, now the single model definition,
  shared by `full_model.py`, `tpi.py` (via a new `tpi_overrides()` for its four `_tpi` swaps)
  and `diagnostics/regimes/regime_model.py`. Pure no-op, verified `main.py` output
  byte-identical (8172 bytes, empty diff).
- **Task 2** — `firm_profit_D/F`: `profit = (1 − mu_p*mc)*(1−alpha)*Y`, the markup rent left
  once labour is paid `mu_p*mc*(1−alpha)*Y` and capital keeps `alpha*Y`. Unrouted this is a
  Walras leak of the W-1/W-2 class.
- **Task 3** — `price_nkpc_D/F`: Rotemberg curves `pi = beta*pi(+1) + kappa_p*(mu_p*mc − 1)`.
  The gap is the *ratio*, so it is unit-free, linearises to exactly `mc_hat` for any `mu_p`,
  and published Calvo slopes map onto `kappa_p` with no SS rescaling.
- **Task 4** — markup wedge in `labor_demand_D/F`: `w = mu_p*mc*(1−alpha)*Y/N`. Employment is
  no longer purely supply-determined. `labor_market_D/F` (labour supply) deliberately
  untouched — wages stay flexible.
- **Task 5** — `terms_of_trade` + `union_inflation` close the nominal side with **no policy
  rate**. `p/p(-1) = (1+pi_F)/(1+pi_D)` pins the inflation differential off the existing
  unknown `p`; `omega_pi_D*pi_D + (1−omega_pi_D)*pi_F = 0` pins the level (the `phi_pi → ∞`
  limit of an ECB rule on union PPI, stated as an abstraction). At `omega_pi_D = 0.071`,
  93% of any terms-of-trade move is Greek deflation, 7% German inflation.
- **Task 6** — markup rent reaches households through `income_D/F`, in proportion to
  productivity `e` rather than lump-sum. `w*N*e + profit*e = (1−alpha)*Y*e` exactly (max abs
  diff 2.2e-16), so the wedge bites only on the firm's hiring decision and household income
  is untouched. `income_D/F` are hetinputs, so a signature change was sufficient to wire
  `profit_D/F` into `hh_extended_D/F.inputs`.
- **Tasks 11–12** — nominal deposit contracts. New `deposit_rates_D/F(i_dep, pi)` emit
  `rdep` (**unchanged name, unchanged ex-ante meaning**) and `rdep_expost` (realised real
  rate, carrying the inflation surprise); `deposit_return_D/F` takes `i_dep`, `P_CES`, `pi`.
  `bank_return_D/F` and `capital_fund_D/F` switch to `rdep_expost_D/F` — the Fisher channel.
  Keeping the `rdep` name meant `intermediation_P1_D/F`, `divert_bond_foc_D/F` and
  `divert_portfolio_adj` needed **zero changes** and remain correctly ex-ante (verified by
  `.inputs` introspection in both directions). T-2 not reopened: the rate is still locked at
  `i_dep(-1)`; only the deflator is period-t. Sovereign bonds stay **real** — a deliberate
  asymmetry that maximises banks' Fisher exposure.

**Calibration** (`code/calibration.py`, Tasks 7 and 14): `mu_p_D/F = 1.20`,
`mc_D/F = 1/1.20` (retargeted from a dead placeholder `1.0`, the subsidy neutralisation that
keeps the SS bit-identical), `kappa_p_D/F = 0.0871` (Calvo θ=0.75 at β=0.985; Bi-Foerster-
Traum's implied 0.0846 to within 3%), `pi_D/F = 0.0`, `omega_pi_D = 0.071` (renormalised
two-country capital key — deliberately *not* GDP weights, which would erase the 93/7 split),
`rdep_D/F → i_dep_D/F`. And **`psi_lambda_B_D/F` 8.5 → 7.85**: stickiness plus Fisher had
pushed peak spread to 162.0bp, an 8% overshoot of the paper's 150bp moment. Re-bisected
(8.5 → 162.14bp, 7.0 → 136.21bp, 7.8 → 149.16bp, **7.85 → 150.14bp adopted**), `b_gov_D[499]`
in the ~1e−5..1e−4 band throughout. `EBA_CALIBRATION` branch only; the pre-EBA `else 3.0`
branch untouched. Bisection table recorded in a comment at the parameter.

**Solver system 23×23 → 27×27** (Tasks 8, 9, 13). `+mc_D, pi_D, mc_F, pi_F` to `unknowns_tp`
(and `rdep_D/F → i_dep_D/F`); `+nkpc_p_res_D/F, tot_res, union_pi_res` to `targets_tp`.
`steady_state.py` carries the six new blocks too (`labor_demand_D/F` deliberately excluded —
SS still uses `labor_ss_D/F`).

**The steady state is bit-identical throughout.** Subsidy-neutralised markups
(`mu_p*mc = 1`, `profit_ss = 0`) and `pi_ss = 0` make every new SS residual exactly
`0.000000e+00`. `goods_mkt_D = -4.2493506589857954e-07`,
`goods_mkt_F = -4.1914559989475464e-07`, `ca_res_D = 6.852157730108388e-17`,
`IC_D: θ − θ_tgt = 1.776357e-15`, `ρ_b = 0.8451` — unchanged at every task, including after
the `psi_lambda_B` re-tune (the dial only touches dynamics). At `pi = 0`,
`i_dep_D = rdep_D = rdep_expost_D = 0.0` exactly.

**Gates passed.**

- *Flex-price equivalence* (Task 9): as `kappa_p → ∞` the 27×27 system reproduces the
  pre-change 23×23 IRFs with textbook O(1/`kappa_p`) convergence — worst relative deviation
  2.925e−03 / **2.925e−04** / 2.925e−05 at `kappa_p` = 1e4/1e5/1e6, gate threshold 1e−3.
  Every one of the 30 IRF series shrinks by exactly 10.00× per decade and the SS levels are
  bit-identical. Binding series `w_D`, `N_D` — the two objects the wedge acts on. Harness:
  `code/dump_irfs.py`.
- *Fisher sign* (Task 13): real → nominal deposits deepens `n_inter_D[0]` −4.0140% →
  −4.6155% (~15% deeper, ~11× the effect on output), reaching output only through the
  intermediary. Had the Task 12 ex-post/ex-ante substitution been backwards, net worth would
  have gone *less* negative.

**SSJ 1.0.0 defect found and worked around** (Tasks 9 and 9b). `CombinedBlock._jacobian`
seeds from the shock list and returns `total_Js[original_outputs & total_Js.outputs, :]`, so
a target reachable from no shock is silently dropped from H_Z; `Block.solve_jacobian` then
hands mismatched shapes to `np.linalg.solve` (`size 11500 is different from 13500`). All four
new targets are pure functions of the solver's own unknowns, so H_Z came back with 23 rows
against a 27×27 H_U. New **`full_model.solve_jacobian_padded()`** restores the rows as zeros
— **exact, not an approximation**, since `dH/dZ` at fixed unknowns is identically zero when
the shock never appears in the equation — and otherwise mirrors `Block.solve_jacobian`
line-for-line, printing the padded row names on every solve. All nine call sites across
`code/`, `experiments/` and `diagnostics/` were converted;
`grep -rn "\.solve_jacobian(" --include="*.py" code experiments diagnostics | grep -v
solve_jacobian_padded` must stay empty. A 25×25 rewrite was considered and rejected — it
would hit the identical defect with smaller numbers.

**Results.** Impact on the 1pp default shock, % of own SS level, both columns on the same
150bp moment:

| | flex, real deposits (8.5) | sticky + nominal (7.85) |
|---|---|---|
| peak spread | 150.4 bp | 150.0 bp |
| `Y_D[0]` | −0.0149 | **−0.5064** |
| `C_D[0]` | **+0.2164** | **−0.5103** |
| `I_D[0]` | −0.7718 | −1.0114 |
| `n_inter_D[0]` | −3.3804 | −4.2962 |

Price stickiness alone (Task 10, deposits still real) does most of it: `Y_D[0]` −0.0149 →
−0.4923 (33×), `C_D[0]` +0.2164 → −0.4904 (**sign flip**), `I_D[0]` only 1.28× — so the extra
output decline is the markup wedge shifting labour demand, not an investment story. The
`kappa_p` sweep {0.03, 0.0871, 0.2} is monotone and stable (`b_gov_D[499]` ~1.5e−05
throughout) and `C_D[0]` is negative across the whole sticky range, so the sign flip is not
knife-edge.

**Caveat, recorded and not to be dropped: this is a one-quarter spike, not a downturn.**
Output and consumption are both positive from quarter 1, and flexible-price consumption is
*more* persistently negative from quarter 2 on. `C_D[1]` is essentially unmoved by nominal
deposits (+0.1141 → +0.1144) — the entire Fisher effect is an impact-quarter effect.
Bi-Foerster-Traum's output stays negative ~20 quarters. The honest claim is that the model
fixes the **impact quarter**, not that it resolves the investment-bust counterfactual.

**E1–E4 regenerated** (Task 15). Regime cache rebuilt with
`diagnostics/regimes/regime_model.py --force` **first**, `experiments/run_all.py` second —
the ordering is load-bearing, because `experiments/` never re-solves the model and would
otherwise have silently re-reported flex-price numbers. New caches tagged
`psilam7p85_cal685f7838`; confirmed consumed via every provenance stamp in
`docs/experiments_results.md`.

- **E1**: γ for the same 0/25/50% compression falls 5.0798 → 3.2515 (medium) and 12.7260 →
  9.0163 (aggressive) — the backstop is more powerful per unit under sticky prices. Loading
  4.00/3.17 → **3.82/2.90**, still monotone decreasing, 4.43 → 1.49 over 59 grid points and
  above 1 throughout. **Live Claims 1 and 5 both survive.**
- **E2**: `market_clearing_D` closes at 3.5e−17 / 1.1e−16 / 2.2e−16 against the 1e−07
  assertion — no Rotemberg resource cost leaked into the resource constraint. But the
  headline-vs-channels finding **reverses**: the largest single channel is now 0.25× the
  headline, where under flex prices it was ~4×. `docs/SPEC.md`'s ΔY caution was restated
  accordingly — it now rests on the channels *cancelling*, not on the headline being the
  smaller object.
- **E3**: `writeoff_enabled=1` alone still negligible and SS-neutral (drift 0.000e+00). Full
  writeoff now inverts Live Claim 1 **only at aggressive** (loading 0.26) — medium holds at
  2.46, where the flex model had both below 1 (0.37/0.28). The appendix robustness claim was
  narrowed accordingly in `CLAUDE.md` and `docs/STATE.md`.
- **E4 + paper artefacts**: `experiments/cache_e4_deciles.npz` was stale and E4 is **not**
  wired into `run_all.py` — `e4_distribution.py` is a separate entry point feeding
  `paper_outputs.py`. Both rebuilt, re-emitting all 8 tracked `experiments/paper/fig0*.png`
  and `docs/paper_draft_results.md`.

**NEW WATCH ITEM.** `n_inter_D[0]` is now **positive (+0.924)** under the aggressive
backstop, where it was −1.099. With `Y_D[0] = +0.8721` and `C_D[0] = +1.5143`, the aggressive
backstop produces an impact *boom* in the crisis country rather than merely cushioning the
bust. The old watch item (`Y_D[0]` positive under intervening regimes) survives and is an
order of magnitude larger.

**Two generated-document hazards found in Task 15**: `run_all.py` carried two prose captions
with flex-price numbers hardcoded as string literals. E3's `psi_lambda_B = 8.5` was merely
stale; E2's "each roughly 4× the headline" was **asserting the opposite of the table printed
immediately above it**. Both now compute from provenance/results. No test covers agreement
between rendered prose and rendered tables — still open.

**Also fixed along the way.** `code/dump_irfs.py` now stores SS levels for every dumped
series with an assertion — `I_D` had no `ss__` entry, so a consumer fell back to a divisor of
1.0 and reported a level deviation as a percentage (the mislabelling class `CLAUDE.md`
already records for `n_inter` and `K`). And one test-authoring bug in Task 5:
`test_closure_puts_93pct_of_tot_move_into_D_deflation` originally asserted a first-order log
identity against the exact nonlinear `tot_res`, whose O(dlog_p²) truncation is `0.429*dlog_p`
in relative terms and swamped its own `rel=1e-6`; replaced with an exact net-rate-split
assertion (`pi_F − pi_D == dlog_p`, `share_D == 1−omega`) that is robust to any
`omega_pi_D`. The blocks were correct throughout; only the test needed fixing.

**Tests:** `code/test_nkpc_blocks.py` 17 passed (~1s); full suite
`code/test_nkpc_blocks.py code/test_eba_calibration.py experiments/` → **40 passed**.

**Task 16 (this commit)** — documentation pass: `docs/STATE.md`, `docs/PROGRESS.md`,
`docs/HANDOFF.md` consolidated from sixteen per-task appendices into one section each;
`docs/SPEC.md` gains the four new modelling choices and its restated ΔY caution;
`CLAUDE.md` gains `build_block_list()`, `solve_jacobian_padded()` and its grep invariant,
the new calibration, the four new residuals in the iteration checklist, and the
`test_nkpc_blocks.py` entry point, plus the corrected S-1 and `experiments/` descriptions.

---

## 2026-08-05 — E4 distributional incidence; net-effects and MS-regime figures  [this commit]

Three additions to the first-draft set, at the author's request.

**`fig06_net_effects`** — contributions to ΔY quarter by quarter with the net path
overlaid. Makes visually what E2 made numerically: the net path is roughly an
order of magnitude smaller than the components generating it, at every horizon.

**`fig07_ms_regimes`** — the empirical scenario chart. A three-state
Markov-switching model on peripheral–Bund spreads dates the intervention stance:
dove through the 1998–2008 convergence era, **hawk across 2010–14**, base
thereafter; ergodic shares 23/52/25%, durations 105/73/49 months. This is the
discipline behind the Stage-B beliefs. The pre-1999 stretch is marked on the
figure as predating the ECB — those regimes are EMU convergence, not a policy
stance, and the shading would otherwise imply one that could not have existed.

**`fig08_deciles` + `experiments/e4_distribution.py`** — DIST-1 addressed. Adds
per-bin consumption hetoutputs to the D household block and re-solves the Jacobian.
**No Gini computed** (author decision; DIST-1 says it is the wrong statistic here).

**The methodological finding, which changed what gets reported.** Binning on
*wealth* with fixed boundaries yields a per-capita consumption response that is
overwhelmingly composition rather than behaviour: bottom decile, PV over 40q, the
consumption term is −41.6 against a mass term of −44.4, netting +2.8. The deposit
distribution shifts across fixed thresholds and membership churns. Binning on the
exogenous **income** state instead makes each bin's mass invariant (it is the
stationary distribution of the exogenous Markov chain), so the measure is purely
behavioural. Income quintiles are now the reported cut; wealth deciles are kept in
the tables with the caveat and are deliberately not plotted.

Three SSJ hetoutput constraints were hit and are documented in the module, since
each fails deep inside the het-block Jacobian with an unhelpful error: the
function's source must be introspectable (no `exec`, no closure), the return must
be a bare single line (a parenthesised tuple yields the name `(cdec01_D`), and it
must carry no trailing comment (which corrupts the last name). `regime_model.
build_tpi_model_main` gained optional `hh_D`/`hh_F` overrides so the augmented
block substitutes into the one canonical block list rather than a second copy.

## 2026-08-04 — S-1 resolved (`writeoff_enabled=0`); first-draft figures and tables

**S-1 resolved by author decision: the paper keeps `writeoff_enabled=0`**, the
pure risk-premium framing. E3 is retained as an appendix robustness result, and
its content becomes a stated caveat rather than an open question — the
over-compensation claim is *conditional on no realised principal writedown*.

`experiments/paper_outputs.py` emits five captioned figures (`experiments/paper/`)
and three tables (`docs/paper_draft_results.md`): a calibration/identification
ledger separating measured from targeted from free parameters, a moment-match
table, and the main results table. Every number is derived live from the solved SS
or the cached response matrices; only source citations are literal text. Captions
are baked into each PNG rather than living in the LaTeX, so they survive reuse.

**Corrected: the default-loading split is 3.1% fundamental expected loss / 96.9%
collateral friction** (`EL_price_D=0.056134`, `psi_spread_D=1.737724`), not the
10.9%/89% recorded in CLAUDE.md and STATE.md at the pre-EBA calibration. This is
the quantitative core of the constrained-seller argument and 96.9% is a materially
stronger version of it.

**New finding — the backstop damps the oscillation, it does not shift the spread
path down.** Building the transmission figure surfaced that regime differences are
concentrated at impact: by quarter four the net-worth and investment paths have
converged, and the spread ordering *reverses* (t=8: passive +0.5bp vs aggressive
+15.5bp; t=12: passive −10.1bp vs aggressive +5.9bp). Passive overshoots downward
later; intervention decays monotonically. An earlier caption claimed the backstop
"cushions each link roughly proportionally", which is true only on impact — it was
corrected rather than shipped.

Figure palette re-validated with the dataviz validator: the project's `#8C1515`
failed the lightness band and `#002147` failed the chroma floor (reads gray). The
paper set uses `#1B6CA8 / #A62B22 / #c87941 / #1a6e3a`, all passing.

## 2026-08-03 — `experiments/` package complete; production regression bit-identical

Consolidation of the E1/E2/E3 work. `CLAUDE.md` gains an `experiments/`
architecture section, three doc-reference rows, and a rewritten **S-1** row
carrying E3's numbers — S-1 is no longer an open flag but a quantified decision
that determines whether the paper's central claim holds.

**`CLAUDE.md`'s `EL_price_D = 0.0717` was stale and is corrected to 0.056134.**
That figure predates the EBA `delta_b = 0.0777` / `q_b = 0.969`. It is the TPI
loading's denominator, so the error propagated into every loading figure quoted
from memory rather than re-derived. E3 made it load-bearing.

**Production regression re-run after all experiment work** (`code/main.py`,
exit 0, 18m46s): **bit-identical to the pre-work baseline** —
`goods_mkt_D = -4.2493506589857954e-07`, `ca_res_D = 6.852157730108388e-17`,
`K_D = 10.800`/`K_F = 10.832`, `n_inter_D[0] = -3.3804%`, `Y_D[0] = -0.0149%`,
`b_gov_D[499] = 1.4e-05`, loading 4.35/4.01/3.44. `git diff main -- code/` is
empty. The package reads the production equation files but writes nothing into
them, so `code/main.py` remains the regression path. Full suite: 31 passed.

**Independent cross-validation of E1.** Two separately written code paths agree:
`code/main.py` gives loading 4.01 at γ=5, E1 gives 4.00 at γ=5.0798; production
3.44 at γ=10, E1 3.17 at γ=12.726 (correctly lower). Peak spreads hit the
compression targets exactly.

**Two items now blocking paper text, both author decisions:**
1. **S-1** — `zeta_writeoff=1` inverts Live Claim 1 (loading 0.37/0.28, below 1).
2. **A5-1's third object is misnamed** — the reported
   `Σ β^t (pd_passive − pd_intervention)` is negative *because* the backstop
   relaxes austerity, so negative means Greece is better off, the opposite of what
   "fiscal saving" implies. Flip the sign or rename it before it is quoted.

## 2026-08-03 — E3: S-1 writeoff. Full writeoff INVERTS Live Claim 1

`experiments/e3_writeoff_s1.py`. S-1 resolved into two nested variants, because
`writeoff_enabled` and `zeta_writeoff` do different things and only the first is
steady-state-neutral.

**E3a (coupon-only, `zeta=0`) is negligible**: peak spread 150.3 → 149.1 bp,
loading 4.00 → 3.93. Measured SS drift exactly **0.000e+00**, confirming
`writeoff_enabled` is strictly SS-neutral. The coupon is only ~7.8% of the bond,
so haircutting it alone barely registers.

**E3b (full, `zeta=1`) inverts SPEC Live Claim 1.** `EL_price_D` goes 0.056134 →
0.701743 (**12.5×**, matching the closed form to 1e−12) and the loading collapses
from 4.00/3.17 to **0.37/0.28 — below 1**. The CB becomes *under*-compensated,
receiving ~30% of the actuarially fair expected loss, where the paper's central
claim is over-compensation. **The "monetary-financing objection fails on the
model's own terms" argument does not survive `zeta_writeoff = 1`.**

The mechanism is attributable to the denominator alone: premium income barely
moves (`prem_PV` +9%) while priced expected loss goes ×11.8. It is a repricing of
the expected loss, not a change in what the CB earns.

`psi_lambda_B = 8.5` also stops hitting its 150 bp anchor under E3b (168.9 bp,
+12.4%). Reported, not re-tuned away.

**Correction to the design spec.** It predicted E3b "moves the steady state via
`EL_price`". More precisely: `EL_price` changes value 12.5× but **no SS allocation
moves** (drift 0.000e+00 across eleven quantities) — it multiplies `def_rate`,
which is 0 at SS, so it is allocation-neutral while still changing the linearised
bond FOC and hence every dynamic result.

**Unanticipated: under E3b the named-regime construction breaks.** Peak spread
stops being monotone in γ, so compression targeting has no unique solution. Two
violations on a 40-point grid over γ∈[0,15]: a trivial one at γ≈0.385 and a large
spike at γ≈3.46 (144.4 → 166.6 bp), after which the curve resumes falling to
82.7 bp at γ=15. The isolated spike sits where `I − γ·A_cb` is plausibly
near-singular, so it reads as a linear-algebra pathology rather than economics —
but it means compression-targeted regimes are undefined under full writeoff. E3
therefore evaluates every variant at the **baseline's** γ held fixed, so the model
changes without the policy also changing.

Two bugs were caught by assertions written before the code was first run. The
`gamma_for_compression` monotonicity guard surfaced the finding above rather than
silently bisecting to a meaningless γ. And `expected_EL_price` was initially handed
a stale calibration because `run()` did `from calibration import get_calibration`
at the top, binding the original function *before* the override context opened —
the exact footgun documented in `calibration_override`'s docstring after the
earlier code review. The closed-form check caught it; the fix is to import the
module and resolve at use time.

## 2026-08-03 — orchestrator: `run_all.py` → `docs/experiments_results.md`

`experiments/run_all.py` runs every experiment and renders the generated results
document. E2 runs first because it is self-verifying, so it validates the cache
before anything else reports numbers off it. `--skip-e3` re-renders without paying
for E3's two model re-solves; `--render-only` rebuilds the document from results
already on disk.

Every table carries a provenance stamp read **live** from the calibration —
including a **working-tree-dirty flag**, so a document generated with uncommitted
edits can no longer be mistaken for a clean run at that SHA. The stamp exists
because `run_regimes.py` once shipped a hardcoded "market-value rule" caption while
actually running the par rule.

The document is generated and carries a do-not-hand-edit warning; it also states
plainly when an experiment's results are missing rather than rendering a partial
table silently.

## 2026-08-03 — E1: backstop schedule; cache schema 3 (`delta_b_F`)

`experiments/e1_backstop_schedule.py`. Named regimes canonical, γ **solved** for
0/25/50% peak-spread compression (0 / 5.0798 / 12.7260). Reports the regime table,
A5-1's three German objects separately, the loading schedule, and welfare labelled
secondary. Figure: `experiments/figures/fig_e1_loading_schedule.png`.

**Every cross-check against the independent `code/main.py` pipeline passes.**
Loading 4.00 at `medium` (γ=5.08) vs production's 4.01 at γ=5; 3.17 at
`aggressive` (γ=12.73), correctly below production's 3.44 at γ=10. Peak spreads
hit the compression targets exactly (150.3 → 112.7 → 75.2 bp). `n_inter_D[0]`
(−3.380 / −2.167 / −1.099 % SS) and `Y_D[0]` (−0.0149 / +0.0111 / +0.0338)
reproduce `docs/STATE.md`'s regime table to every printed digit.

**Live Claim 5 confirmed on a fine grid.** The loading schedule is monotone
decreasing at **all 59 finite grid points**, 4.51 (γ=0.51) → 2.07 (γ=30) —
stronger evidence for the self-extinguishing premium than the three points
previously on record.

**Cache schema 3: `delta_b_F_ss` added, and the reason is a bug this caught.**
`cb_pnl` computes each carry leg's SS yield as `delta_b·(1/q_b_ss − 1)`. The first
draft used `delta_b_D` on *both* legs, but `delta_b_F = 0.056779 ≠ delta_b_D =
0.077701` — the two countries' bank books have different measured maturity
ladders. That put the SS spread at −9.2e−04 rather than its true ~1e−17 and would
have silently contaminated `carry_ss_pv`. An assertion written into `cb_pnl`
before the code was ever run caught it. `carry_ss_pv` now comes out 1.2e−16 /
3.2e−16, i.e. numerically zero as it should be. Cache rebuilt (schema 3); E2
re-ran with identical γ and `dY` values.

**Open, needs an author decision: A5-1's third object is misnamed.** The code
reports `Σ β^t (pd_passive − pd_intervention)`, which comes out **negative**
(−0.0015 medium, −0.0047 aggressive) because the backstop lets Greece run a
*larger* primary deficit — it relaxes required austerity. A negative number
therefore means Greece is better off, the opposite of what "Greek fiscal saving"
implies. Either flip the sign or rename it ("austerity relief, PV"). Magnitudes
are unaffected; the label must not ship as-is.

## 2026-08-03 — E2: ΔY decomposition against the `market_clearing_D` identity

`experiments/e2_dy_decomposition.py` + `experiments/test_e2_identity.py`, and five
cache/IRF helpers appended to `experiments/common.py` (`load_cache`,
`cache_outputs`, `irf_from_cache`, `named_regime_gammas`, `regime_irfs`). The
named regimes are canonical: γ is **solved** for 0/25/50% peak-spread compression
(0 / 5.0798 / 12.7260), not chosen as round numbers, so the regimes keep their
meaning across recalibrations.

**Self-verifying by construction.** The decomposition is the linearised
`market_clearing_D` identity `dY = P_ss·dC + C_ss·dP + dI + dG + dΦ + dT + dNX`,
and `goods_mkt_D` is a *targeted* solver residual, so the components must sum to
`dY` to solver tolerance. Achieved **5.8e−17 / 1.1e−16 / 1.5e−16** against a 1e−7
assertion that halts rather than warns.

**Finding — the headline output number is the residue of two much larger
offsetting channels.** Passive → aggressive, `dY[0]` moves +4.87e−04 while
investment moves +2.16e−03 and net exports −1.85e−03, each ~4× the headline and
opposite in sign. This confirms `docs/SPEC.md`'s standing caution ("a small
headline output number can be *only* small because two large channels are netting
out — and they land on different households") **as a measured property of this
calibration**, and settles the gate it placed on the trade-channel claim.

**The `Y_D[0] > 0` watch item is answered.** `dY[0]` = −0.0149 / +0.0111 /
+0.0338 % of SS. The proximate driver of the sign flip is **consumption quantity,
not investment**: at `medium`, investment is still negative on impact and output
is positive only because consumption outweighs it. Consumption is already positive
on impact at `passive` with no backstop, so it is not manufactured by the policy
rule — which argues against pure linear-rule overshoot. But the magnitudes
(0.01–0.03% of SS) are small differences of much larger terms and should not be
leaned on. Full table in `docs/STATE.md`.

`Phi_D` and `G_D` are **verified** zero rather than merely uncached: `Phi_D` has
no Jacobian column (the portfolio adjustment cost is quadratic about its anchor,
so its level deviation is second-order), and `G_D` is absent from `G_tpi.outputs`
because government spending is constant. The identity closes *because* both are
genuinely zero. Cache rebuilt under schema 2 (7m27s); `G_tpi[cb=0]` vs baseline
`max|err| = 0.00e+00`.

Also noted: `EL_price_D` is **0.056134** at the live calibration, not the `0.0717`
still quoted in older doc sections and CLAUDE.md — that figure predates the EBA
`delta_b`/`q_b`. It is the TPI loading's denominator, so it must be re-derived
wherever quoted, not copied. Not fixed in this commit.

## 2026-08-03 — `experiments/` package: cache schema v2 (call-time fingerprint)

First commit of the new `experiments/` package (branch `experiments`), which will
produce the paper's standard policy results on top of the regimes cache layer.
`code/main.py` untouched.

**The fix that had to come first.** `regime_model.cache_path` built its filename
from `CAL_FINGERPRINT`, computed at module **import**. Experiment E3 needs to solve
the model under a calibration override applied at run time; with an import-time
fingerprint its cache would have been written to the baseline filename and silently
overwritten it. Now computed at call time, stamped into the `.npz` as
`cal_fingerprint`, and asserted against the live calibration on load, with a clear
`FileNotFoundError` naming the rebuild command when no cache matches.

`CACHE_SCHEMA = 2` now appears in the filename: the calibration fingerprint alone
cannot detect a change to the cached *output list*, so without it an old cache would
reload under an unchanged name missing the new keys — invisibly, since `irf_all`
discovers outputs by scanning cache keys. Added `Phi_D` and `def_rate_D` to
`REQUIRED` (`Phi_D` closes the `market_clearing_D` identity for E2; `def_rate_D` is
the off-path expected-loss leg for E1). `T_D` went to `OPTIONAL` deliberately —
`T0=T1=0` makes it identically zero, so zero-filling is correct rather than a silent
hole, and E2's closure assertion catches it if that ever changes.

Also `build_caches` now reads `psi_lambda_B` live rather than from the import-time
`PSILAM_MAIN` constant, so an override wins there too.

New: `experiments/common.py` (calibration-override context manager, unit helpers,
provenance stamp, results writer) and `experiments/test_common.py` (regression
guards including one that the override changes the cache filename and does not
leak out of the context manager).

**Hardened after code review (same day).** `calibration_override` now raises
`KeyError` on an override key not present in the calibration dict — previously a
typo (e.g. `psi_lambda_b_D` for `psi_lambda_B_D`) would silently add a junk key
while leaving the real parameter at its default, producing a wrong-but-plausible
number without any error: exactly the failure mode this whole package exists to
close off, one level up from the cache-fingerprint fix above. `write_results` now
serialises numpy arrays/scalars properly (`json.dump(default=float)` raised on any
multi-element array — the normal shape of an IRF payload) and passes
`allow_nan=False`, so a `NaN` in a result is a loud `ValueError` at write time
rather than a token that travels silently into a table and that strict JSON
parsers reject anyway. Also fixed: the exception-restore test asserted on a
name bound at import (`from calibration import get_calibration`), which cannot
observe the module attribute the context manager patches and so passed
regardless of whether the restore worked — rewritten to assert on
`calibration.get_calibration` directly, and verified to fail when the `finally:`
restore is removed. `cache_path`'s `fingerprint` parameter is now actually used by
`load_cache` instead of being computed and discarded twice. 17/17 tests pass
(`experiments/` + `diagnostics/regimes/test_lottery_math.py`).

## 2026-07-31 — regimes re-run at the broad scope; `PSILAM_BREAKDOWN` re-derived; units bug fixed

Follow-up to the broad-scope commit below: the policy-regime diagnostics had never
been run at the live calibration, and their hard guard blocked it.

**`PSILAM_BREAKDOWN` re-derived: 2.5 → 15.0.** New
`diagnostics/psilam_breakdown_sweep.py` sweeps 16 points above the live 8.5,
solving the SS once (`psi_spread` is exactly linear in `psi_lambda_B`) and
re-solving the Jacobian per point with both dials moved together. Peak spread
150.3 → 223.8 → 273.6 → 625.3 → 1034.5 → 8903.8 bp at 8.5 / 14 / 20 / 25 / 26 / 27,
then **sign-flips at 28** — a pole between 27 and 28. The A7 >1000bp flag first
fires at 26; the *first* pathology is earlier, `n_inter_D[0]` shrinking over
`[14,18]` while the spread still rises. The guard is set from that, not from the
pole, so it keeps real margin. The old 2.5 was CT1-scope and would have blocked
the live calibration outright.

**Regimes re-run, all three figures regenerated.** Caches rebuilt (the calibration
fingerprint forced it). `A_cb=-1.889e-2` — backstop still compresses, SA-1 absent.
`gamma_aggressive=12.726` / `gamma_medium=5.080`; peak spread 75.2 / 112.7 /
**150.3** bp, so `run_regimes.py`'s 120–180bp band now passes. A6 amplifier
invariance holds in both Stage A (3.78/4.67/5.44 at `psi_lambda_B=0`) and the
lottery (2.51/3.09/3.59). 18/18 tests pass.

**New: `diagnostics/regimes/certainty_equivalence.py`** — answers whether the
regime lottery is degenerate under first-order certainty equivalence. It is not,
but the non-degeneracy is ~9% of what the Stage B figure shows: comparing the
lottery to one *known* CB at `gamma_bar` (same silence-until-`k` timing) gives
`LOT − CE = +1.477bp` at ergodic beliefs against a 16.6bp total belief-shift
effect. The wedge is exactly `A_cb[0,:] @ (cb^e − cb_CE)` (verified to 6e−15 bp);
mechanism is that uncertainty **back-loads** the expected backstop while the
date-0 spread weights near-term purchases most. See `docs/STATE.md`.

**Units bug fixed (pre-existing, affected published numbers).** SSJ IRFs are
*level* deviations, so `×100` is a percent only where the SS level is ≈1.
`Y_D_ss≈1` passes; `n_inter_D_ss=2.138` and `K_D_ss=10.8` do not. `main.py`
printed `n_inter_D[0]×100` as `%` — the widely-quoted **−7.227%** is the level
deviation and the true impact is **−3.380% of SS** — and the Stage A figure
titled two panels `(%)` on the same basis (2.1× on net worth, **10×** on
capital). Both fixed; level deviations retained alongside for continuity.
Consequence: `PT-1`'s pass-through moment is **−2.25%/100bp**, not ≈−4.5%; still
inside the Acharya–Drechsler–Schnabl-implied range but at its low end.

**Watch item recorded:** `Y_D[0]` is positive under both intervening regimes and
the A5 `dY_D` trough never goes negative — output never falls at all under the
backstop. At `gamma_aggressive=12.7` that is more likely linear-rule overshoot
than economics; flagged in STATE.md, not reported as a result.

Full `code/main.py` re-verified (exit 0): residuals, `b_gov_D[499]=+1.4e−05`,
`ρ_b=0.845`, peak spread 0.376pp, TPI loading 4.35/4.01/3.44 all unchanged.

---

## 2026-07-31 — broad capital-funding-sector scope: the EBA calibration goes LIVE

The last blocker was the **scope of `n_inter`**, not any parameter. EBA CT1 is the
capital of the *sovereign-exposed stress-test sample*; the model's `n_inter` is
the net worth of the agent intermediating the **whole capital stock**. Using CT1
pins `omega_K` tiny and makes the accelerator gain ~`1/n_inter`.

**New `BANK_SCOPE` in `code/calibration.py`** (`"broad"` live, `"ct1"` kept for
comparison). Under `"broad"` the intermediary is the entire capital-funding
sector and its net worth follows from the measured leverage and the balance
sheet, `N = (Q*K + sovereign)/theta`, so **`omega_K = 1` by construction** and the
passive-fund device disappears entirely.

| | CT1 scope | **broad (live)** | pre-EBA placeholder |
|---|---|---|---|
| `n_inter_D/F` | 0.408 / 0.175 | **2.138 / 1.627** | 3.0 / 3.0 |
| `phi_own_D/F` | 2.390 / 2.758 | **0.456 / 0.296** | 0.25 / 0.25 |
| `omega_K_D/F` | 0.117 / 0.067 | **1.0 / 1.0** | 1.0 |
| `theta_D/F` | 5.51 / 6.94 (measured) | **5.51 / 6.94** | 4.0 |

Kept measured: `theta`, the sovereign book, `delta_b` (ladder), `K/Y`. Given up:
`n_inter` as observed CT1, and `phi_own = 2.39` as a *model* parameter — 2.39 is
concentration *within the stress-tested slice*, not within the whole
capital-funding sector, and only the latter is what the model's `phi_own` means.
**Load-bearing assumption:** applying the EBA sample's `theta` to the whole
sector. This is now the only such assumption left in the bank block.

Consequences: `Delta` returns to **0.2/0.4** (the 0.85 bound was a CT1 artifact;
at `phi_own=0.456`, `f*theta = 0.661 > 0.367`), the fiscal rule to
`phi_lamb=0.15` / `mv_rule=0`, and `psi_lambda_B` retunes to **8.5**.

**Verified end-to-end** (`code/main.py`, exit 0):

| Check | Result |
|---|---|
| over-identifying `K` | `K_D = 10.800`, `K_F = 10.832` (target 10.8) |
| IC residual | −8.9e−16 (D) / 0.0 (F) |
| Walras | `ca_res_D = 6.9e−17`; all block residuals < 1e−8 |
| stability | `b_gov_D[499] = +1.4e−05`, `ρ_b = 0.845` |
| `n_inter_D[0]` | **−7.227%** ✓ |
| `Y_D[0]` | **−0.0149%** ✓ — **Y-1 RESOLVED** |
| `rk_D`, `rk_F` | both exactly 0.010000 — **RK-1 RESOLVED** |
| peak spread (γ=0) | 0.376pp = **150.4bp**, on target |
| TPI loading | **4.35 / 4.01 / 3.44** at γ=2/5/10, declining |

The paper's self-extinguishing-premium claim survives on a properly identified
calibration.

`psi_lambda_B = 8.5` is higher than the historical 0.31 / 1.18 / 3.0 because the
broad scope's `phi_own = 0.456` is well below CT1's 2.39, so more of the default
loading comes from the friction. The mapping is smooth and monotone with **no
breakdown region** — 5.4 / 23.1 / 58.6 / 111.0 / 142.5 / 157.8 bp at
0 / 1 / 3 / 6 / 8 / 9 — so the old "breakdown above ~1.5-2.0" warning (specific to
CT1-thin net worth) does not apply.

**Method bug found and fixed.** `psi_spread_D` is derived from `psi_lambda_B`
inside `_apply_ss_anchors`, so a sweep **must re-solve the SS per point**.
Patching the flag onto an already-solved SS leaves `psi_spread` stale and
*inverts the apparent sign* of the spread response — an earlier sweep did exactly
that and reported the spread falling in `psi_lambda_B`. Those numbers are void;
all figures above come from full per-point re-solves.

---

## 2026-07-31 — fix `omega_K`: capital fund holds a fixed QUANTITY, not a fixed share  [301ffd2]

`omega_K` as a **fixed share** was the defect. The passive fund held
`(1-omega_K)*K` at all times, so it mechanically *mirrored* bank deleveraging —
a 1% fall in the bank's book dragged the other ~88% of the capital stock with
it. That is the `1/omega_K` lever, and it is an assumption nobody would defend
stated plainly: non-bank capital holders do not shrink in lockstep with bank
equity.

**Fix: `fund_rule_D/F`** (`code/calibration.py`, committed at **1**).

| `fund_rule` | fund holds | bank holds | `dK/dN` |
|---|---|---|---|
| 0 (legacy) | `(1-omega_K)*K` | `omega_K*K` | `theta/omega_K` |
| **1** | **constant `K_fund`** | `K - K_fund` | **`theta`** |

With `K_fund = (1-omega_K)*K_ss` the rules are **identical in steady state** —
the change is purely dynamic. Verified: `lambda_gk_D=+0.9271`, `Omega_D=+4.62`,
`K_D=10.800` under both. Gain 47.1 → 5.5. At `omega_K=1` (pre-EBA), `K_fund=0`
and the rules coincide exactly, so the placeholder calibration is untouched by
construction — confirmed by a full `main.py` run reproducing
`n_inter_D[0]=-3.0009%`, `Y_D[0]=-0.0261%`, `ρ_b=0.8451`, `b_gov_D[499]=-1.3e-5`,
IC residual `+8.9e-16`.

Touched four sites per country (`smart_steady`, `intermediation_IC` kappa,
`capital_fund`, `k_balance_sheet`) plus the `steady_state.py` over-identifying
print, which is now exact for both countries (`K_D = K_F = 10.800`).

**Helps substantially, does not stabilise alone:** `b_gov[499]` +3.98e2 → +2.25e2
at `psi_lambda_B=0`, and −2.08e4 → −1.23e3 at `psi_lambda_B=1` (17×).

**Three compounding amplifiers** (at `fund_rule=1`, `psi_lambda_B=0`):

| bank block | `phi_own` | `b_gov[499]` |
|---|---|---|
| pre-EBA (`n=3.0, theta=4.0, omega_K=1`) | 0.25 | **−4.4e−08 stable** |
| EBA (`n=0.41, theta=5.51`) | 0.25 | +1.85 |
| EBA | 2.39 (measured) | +5.32e+02 |
| pre-EBA | 2.39 | solver failed to converge |

1. the fixed-share fund (`1/omega_K` ≈ 8.5×) — **now fixed**;
2. thin measured net worth (`n_inter` 3.0 → 0.41): ~1e−8 → 1.85;
3. measured concentration (`phi_own` 0.25 → 2.39): 1.85 → 532.

(2) and (3) are both *measured*, so neither can be tuned away. `EBA_CALIBRATION`
stays **False**.

**Also ruled out this session** (all with evidence in `docs/eba_calibration.md`):
the sovereign-risk schedule — `def_scale` 0.25 → **0.00** makes it *worse*
(−2.1e4 → −1.6e5), so flattening it (e.g. a bounded `tanh`) cannot help; note
also that at first order only the local slope `a*b` matters, so a bounded form is
observationally equivalent in the linearised IRFs. Separately, measuring
`B_supply` (2.4 → 1.116) steepened that schedule's local slope ~45% (0.0799 →
0.1158) as a side effect, since lower debt sits further up the concave curve —
real, but not the driver.

---

## 2026-07-31 — fix the collateral mapping; EBA steady state now correct, dynamics still explosive  [c6230a2]

Follow-up to the rebuild below. The previous entry concluded the feasible `Delta`
set was empty. **That was wrong, and the reason is instructive:** the emptiness
came entirely from an undocumented convention, not from the economics.

`ic_delta_calibration._ic_delta` back-solved `Delta` from one equation in two
unknowns and closed the system with a hardcoded
`ratio = Delta_cross/Delta_own = 2.0`. With `Delta_cross <= 1` that caps
`Delta_own <= 0.5`, against the `> ~0.73` GK well-posedness requires. It was a
back-solve closure masquerading as a consistency check — which is exactly why the
inherited `0.2/0.4` pair "passed" it.

**Fixed.** The convention is gone. `Delta_own`/`Delta_cross` are free structural
parameters and the module now checks the IC **residual** directly
(`ic_residual`, tol 1e-8, plus a positive-divertable-leverage check). Verified on
the pre-EBA calibration: residual `+8.9e-16` (D) / `+1.8e-15` (F).

**`Delta = 0.85/0.90` resolves the steady state:**

| `Delta_own`/`Delta_cross` | `lambda_gk_D` | `lambda_gk_F` | `Omega_D` | `Omega_F` | `K_D` | `K_F` |
|---|---|---|---|---|---|---|
| **0.85 / 0.90** | **+0.927** | **+0.960** | +4.62 | +5.98 | 10.80 | 10.65 |
| 0.90 / 0.95 | +0.488 | +0.456 | +2.49 | +2.91 | 10.80 | 10.65 |

`lambda_gk_D = +0.927` is essentially the pre-EBA `+0.923`, so the amplification
block keeps its strength with measured concentration. Not a fudge: leverage of
5.5× on a 43%-sovereign book is inconsistent with sovereigns being good
collateral — if they were, the bank would lever further and `theta=5.5` would not
bind. Also the right story for 2010–12 Greece (collapsing GGB eligibility, rising
ECB haircuts). Note `Delta=0.80/0.90` lands just past the `lambda_gk` **pole**
(`lambda_gk_F=-12.45`): the closed-form bound ignores endogenous `rn`, so the real
frontier is a fixed point — sweep, don't solve.

**STILL BLOCKING — dynamic instability (GK-2).** With the SS correct,
`b_gov_D[499] ~ 1e2–1e3` vs a ~1e-5 target. Measured amplification is
`theta*phi_own = 13.17` against `4.0*0.25 = 1.0` for the placeholder. Diagnosed:

| Test | Result |
|---|---|
| `psi_lambda_B = 0` | still explosive, `b_gov[499] = -2038` |
| `phi_lamb` 0.6 → 25 | **flat** (peak ~1.1e7bp at 0.6, 1.5 *and* 25) — not the fiscal mode |
| `mv_rule = 1` | does not fix it |
| `chi1` 0 → 0.5 | peak spread **1.1e7bp → 6.0bp**, `b_gov[499]` −2038 → +70 |
| `chi1` ∈ [0.2, 5.0] | `b_gov[499]` stays 70–560 — amplitude damped, root not removed |

`chi1` (Auclert intermediary adjustment cost, currently 0) is the strongest lever
and restores a sane spread response, but no tested value removes the root. Three
open routes in `docs/eba_calibration.md` → *Dynamic instability*.

**Also fixed:** `steady_state.py` did `from calibration import EBA_CALIBRATION`,
binding the flag at import time so a sweep flipping it would silently keep the old
portfolio targets — the same stale-binding trap as the regimes cache key and
`PSILAM_MAIN`. Now imports the module and resolves at call time. The GK guard also
reports both countries' `lambda_gk`/`Omega` on failure, not just the failing one.

`EBA_CALIBRATION` stays **False**. Full `main.py` re-verified: `n_inter_D[0]=-3.0009%`,
`Y_D[0]=-0.0261%`, `ρ_b=0.8451`, `b_gov_D[499]=-1.3e-5`, all residuals < 1e-8.

---

## 2026-07-31 — EBA calibration rebuilt and identified; measured moments found infeasible in the GK block  [050e54a]

Rebuilt the EBA 2011 moment set to be **identified rather than back-solved**, then
found that the measured moments are structurally incompatible with the model's
collateral constraint. The rebuild is sound; the incompatibility is the finding.

**Live calibration unchanged.** `EBA_CALIBRATION = False` (`code/calibration.py`)
selects the pre-EBA values bit-exactly — verified parameter-by-parameter and by a
full `main.py` run reproducing `n_inter_D[0]=-3.0009%`, `Y_D[0]=-0.0261%`,
`ρ_b=0.8451`, `b_gov_D[499]=-1.3e-5`, TPI peaks 0.468/0.330/0.236/0.163 pp.

**Newly measured** (`code/eba_calibration.py` → `data/eba_moments.json`):

| Parameter | Was | Now | From |
|---|---|---|---|
| `delta_b_D/F` | 0.10, no EBA counterpart | 0.0777 / 0.0568 | maturity ladder (`MATURITY_CODE` 125–155, previously unread) repriced at the end-2010 market yield |
| `theta_D/F` | 4.0 assumed | 5.51 / 6.94 | (corp ex-CRE + CRE + sovereign) EAD / CT1 |
| `omega_K_D/F` | back-solved plug | 0.117 / 0.067 | corp+CRE EAD ÷ K at `K/Y_ann=2.7` |

- **`delta_b` retires the F-1 duration blocker.** The old "port 0.036/0.038 (7y)"
  target measured the sovereign's *whole outstanding stock*; `delta_b` governs the
  *bank-held book at the yields banks faced*, whose modified duration is 3.12y
  (GGB) / 4.22y (Bund). Close to the old 0.10 — `mv_rule=1` + `phi_lamb=0.60` not
  needed.
- **`omega_K` kept, not dropped.** Banks fund ~12%/7% of the capital stock, so
  `omega_K=1` is counterfactual. Now measured, with `K` an output: the balance
  sheet delivers `K_D = 10.800` vs the 10.8 target.
- **Amplification moment is Acharya–Steffen, not CT1 depletion.** Ladder ×
  observed 2011 yield moves (GGB 12.01%→21.14%) gives mechanical
  **−5.73%/100bp** of CT1, −39.8% realised over 2011. The pre-EBA calibration
  generates −0.61%/100bp, **~10× too weak** — almost entirely `phi_own`
  (0.25 vs 2.39), not duration. This is what `psi_lambda_B=3.0` was standing in
  for, and it reframes S-1's "89% collateral friction" as a calibration artifact.
  The 2011 adverse scenario is *deliberately rejected* (it excluded banking-book
  sovereign default) and guarded by a test.

**BLOCKING FINDING — GK feasibility.** The block is well-posed only if
`f*theta > (1-Delta_own)*phi_own + (1-Delta_cross)*phi_cross`. At the measured
moments this is violated by −1.26 (D) / −1.42 (F):

```
lambda_gk_D = -0.0869   Omega_D = -0.3013
lambda_gk_F = -0.0723   Omega_F = -0.3217
```

Negative IC multiplier and negative banker franchise value — yet the solver
converged, `ca_res_D = -7.7e-17`, the IC-δ check passed exactly, stability passed,
and the TPI loading still declined. Every IRF meaningless. C-1's failure mode
reborn.

The bound is `Delta_own > ~0.73`, but `_ic_delta` hardcodes
`ratio = Delta_cross/Delta_own = 2.0`, which with `Delta_cross <= 1` caps
`Delta_own <= 0.5`. **Feasible set empty.** `f` would need > 0.349 (literature
0.03–0.12), `theta` > 16.03 (measured 5.51). Tested: `Delta = 0.80/0.90` →
`lambda_gk_F = -12.45`, `Omega_F = -75.91`, just past the pole. Escaping is a
modelling decision; clearing the pole needs `Delta_own ≈ 0.85–0.95`, i.e.
sovereign bonds nearly worthless as collateral — removing the channel the doom
loop runs on.

**New guard: `steady_state.assert_gk_well_posed`**, called from
`_apply_ss_anchors`, so it fires on every solved SS in both `steady_state.py` and
`depreciation_calibration.py`. Exports `gk_feasibility_margin` / `min_Delta_own`.
The most valuable artifact here — it makes this failure class impossible to commit
silently.

**Also fixed** (all pre-existing, surfaced while working):
- `diagnostics/regimes/` cache filenames keyed only on `psi_lambda_B`, so the
  `psilam=0` cache would be **silently reused across different models**. Now
  carries a SHA of the whole calibration. Stale pre-EBA caches deleted.
- `PSILAM_BREAKDOWN` 4.0 → 2.5, with its net-worth dependence documented.
- `run_regimes.py`'s figure suptitle hardcoded `"market-value rule"` while the
  committed calibration is the par rule; now reads live `mv_rule`. (The earlier
  fix covered the JSON provenance string only.)
- `data/README.md`: full worksheet-4/5 code map and the maturity-bucket decode.

**Tests:** `code/test_eba_calibration.py` 10/10 pass, including ladder-exhausts-total,
duration-below-maturity, `delta_b` round-trip, the balance-sheet identity, the
mechanical-MTM magnitude, and a guard that the adverse scenario is never used.

**NOT done:** the `psi_lambda_B` retune to 150bp. A sweep was run but on the
degenerate model, so it is **void and discarded**; redo it after the GK
feasibility question is settled.

---

## 2026-07-31 - fix A6 lottery invariance test (measurement-window error)

The A6 amplifier-invariance check in `uncertain_regime.py` ranked the **full-sample**
peak spread across the three lottery branches. That peak is the common pre-`k` spread —
no branch has acted before revelation (the same pre-`k` identity asserted in sec 10.3) —
so it is identical across branches *by construction*, and `pk0[0] < pk0[1] < pk0[2]` was
comparing floating-point noise. Measured both ways:

| `psi_lambda_B` | full-sample peaks | strict `<` | gap |
|---|---|---|---|
| 3.0 (calibrated) | all `184.992662` bp | **False** | `-5.1e-13` bp |
| 0 | all `9.302930` bp | **True** | `+2.3e-14` bp |

Same expression, opposite verdicts, decided by the last ULP — the reported "YES" was
luck, not a result.

**Fix:** rank on the **post-revelation window `t>=k`**, where branches have actually
diverged; check **both** amplifier settings (that is what makes it an invariance test
rather than a single reading); and require separation above a `1e-3` bp margin so
numerical noise can never produce a verdict.

| `psi_lambda_B` | aggressive | medium | passive | ordered |
|---|---|---|---|---|
| as calibrated | 76.98 | 117.87 | 160.19 | YES |
| 0 | 4.91 | 5.89 | 6.43 | YES |

**A6 invariance genuinely holds in the lottery** (separation 0.548bp at `psi_lambda_B=0`,
~550x the margin) — the earlier "false pass" was a measurement-window error, not a
failure of the economics. The deterministic Stage A A6 (9.33/8.54/7.14bp) is unaffected
and independently confirms it.

No other Stage B output changed: `E_pi[W_D]=-14.3343`, impact spread 151.46 -> 187.19bp,
the `k` sweep and the welfare table are all identical. Full suite re-verified, exit 0:
`test_lottery_math.py`, Stage A, Stage B-lite.

Not touched: the Stage B table's `lottery` column is 185.0bp for all three branches, and
that is *correct and deliberate* — before revelation nobody can tell the branches apart,
which is exactly what the uncertainty premium (lottery minus known-delayed) prices.

---

## 2026-07-31 - doc-sync enforced by a git hook; regimes feature unblocked and verified  [this commit]

**Doc-sync is now actually enforced.** Correction to yesterday's entry: the policy *was*
already implemented as a Claude Code PreToolUse hook
(`.claude/hooks/require-docs-before-commit.sh`, tracked); the "not installed" claim came
from checking `.git/hooks/`. Two real gaps found:

- A PreToolUse hook cannot see terminal commits at all. Added **`.githooks/pre-commit`**,
  the git-native twin and now the primary enforcement: enable once per clone with
  `git config core.hooksPath .githooks`. It inspects the git *index*, so unlike the
  command-string-matching PreToolUse hook it cannot false-positive. Verified: exits 1
  with the block message when code is staged without the docs.
- The PreToolUse hook's `if` filter is prefix-matched, so it never fired on compound
  commands (`cd X; git commit ...`). Removing the filter was tried and **reverted** - it
  makes the script scan every Bash command string and it then denied an unrelated
  documentation edit. The filter stays; the git hook is the real gate, and the
  PreToolUse hook is early feedback for the simple-command case.
- Stale `PROCESS.md` reference in the script header corrected to `PROGRESS.md`.

**`run_regimes.py` provenance string no longer hardcoded.** It read `"mv_rule=1"`
regardless of calibration, silently mislabelling every number in
`regimes_calibration.json`. Now reads `psi_lambda_B`, `mv_rule`, `recovery_rate` live.

**Policy-regime feature now runs end-to-end at the pre-EBA calibration** - it did not
before. `gamma_for_compression`'s scan range narrowed 60 -> 25: bisection's validity
condition is monotonicity on the *bracketing interval*, and `peak(gamma)` falls
monotonically 187.2 -> 34.2bp on [0,25] then ticks up 1.1bp at gamma=30 (saturation at
81% compression, not economics), which aborted the whole run.

Verified, all exit 0: `test_lottery_math.py`, Stage A `run_regimes.py`, Stage B-lite
`uncertain_regime.py`.

- Stage A: `A_cb=-2.406e-2` (compresses; SA-1 absent), passive peak 187.2bp,
  `gamma_aggressive=5.0813` (50% compression), `gamma_medium=1.5730` (25%). Peak spread
  187.2 / 140.4 / 93.6 bp passive / medium / aggressive.
- Stage A A6 (deterministic, `psi_lambda_B=0`): 9.33 / 8.54 / 7.14 bp - amplifier
  invariance genuinely holds, clean separation.
- Stage B-lite: sec 10.3 assertions PASS; uncertainty premium +69.0 / +34.4 / -2.2 bp
  (aggressive / medium / passive); impact spread rises with the passive belief weight
  (151.5 -> 187.2bp), sign computed not targeted.

**Open - A6 invariance at the LOTTERY stage is a false pass.** All three branch peaks are
numerically identical (`9.302980` bp, gap `0.000000`), yet
`pk0[0] < pk0[1] < pk0[2]` returns True on last-ULP floating-point noise. Structural, not
a coding slip: with `k=2` the peak spread falls in the common pre-`k` window where no
branch has acted, so the branches coincide there by construction. Stage A's deterministic
A6 is the informative test; the lottery A6 line should be reworded or dropped rather than
reported as a result.

**Also open, none EBA-related:** `psi_lambda_B=3.0` gives 187.2bp, outside
`run_regimes.py`'s own 120-180bp sanity band (it logs "investigate"); `delta_b=0.10`
still needs `mv_rule=1` AND `phi_lamb=0.60` jointly; `beliefs.json` dates from
2026-07-23, before the calibration revert.

---

## 2026-07-30 — revert calibration to pre-EBA; drop audit_artifacts; F-1 hard-break measured  [this commit]

**Calibration reverted to the pre-EBA values** in force at `abcbb6e` (the last commit
before the EBA work began at `eade414`). Verified parameter-by-parameter against that
commit: the only additions are `omega_K_D/F=1.0`, which is the structural no-op
(capital fund empty, `div_fund=0`, so the pre-EBA balance sheet is recovered exactly).

| Parameter | EBA | now |
|---|---|---|
| `psi_lambda_B_D/F` | 1.1793 | **3.0** |
| `n_inter_D/F` | 0.408 / 0.175 | **3.0 / 3.0** |
| `omega_K_D/F` | 0.0602 / 0.0190 | **1.0 / 1.0** |
| `B_supply`/`b_gov`/`b_gov_ss` | 1.19 / 0.591 | **2.4** |
| `phi_lamb_D/F` | 0.60 | **0.15** |
| `mv_rule_D/F` | 1.0 | **0.0** (par) |
| `phi_bF_D_ss` / `phi_bD_F_ss` | 0.018 / 0.069 | **0.25 / 0.25** |
| portfolio targets (`steady_state.py`) | 2.39/0.018/0.069/2.76 | **0.25/0.15/0.15/0.25** |

**All structural fixes retained** — C-1 (multi-asset `lambda_gk`), W-1/W-2/W-3, T-2, A-2,
TPI-1, the `omega_K` generalisation and the capital-key conduit are untouched. Only
parameter values moved. `recovery_rate_D/F=0.30` (EL-1) was **kept**, not reverted: 0.00
asserts a counterfactual 100% loss-given-default.

**Verified end-to-end** (`code/main.py`, exit 0): `n_inter_D[0]=-3.0009%`,
`Y_D[0]=-0.0261%` (both negative ✓); `b_gov_D[499]=-1.3e-5` (default) / `+7.8e-5` (TFP);
`rho_b=0.8451`; IC-δ exact at 0.2000/0.4000; `max|ca_res_D|=6.3e-8`,
`max|goods_mkt_F|=1.1e-9`. TPI peak spread monotone in γ: 0.468→0.330→0.236→0.163 pp.

**Finding F-1 sharpened — the `[0.15,0.18]` zone is a hard break, not a mild one.**
Measured directly by switching `mv_rule` on the otherwise-pre-EBA calibration:

| `phi_lamb` | `mv_rule` | `n_inter_D[0]` | `Y_D[0]` | `b_gov[499]` | spread |
|---|---|---|---|---|---|
| 0.15 | 1 | **-1554.0%** | **+0.170%** ✗ | 1.6e-2 | 124.0bp |
| 0.60 | 1 | -5.89% | -0.024% ✓ | 0.0 | 219.8bp |
| 0.15 | 0 | -3.00% | -0.026% ✓ | -1.3e-5 | 187.2bp |

`mv_rule=1` and `phi_lamb=0.15` are **not a usable pair**. Porting empirical duration
(`delta_b=0.036/0.038`) is therefore a two-parameter move, not one.

**Default-loading decomposition measured.** `EL_price_D=0.0717` vs `psi_spread_D=0.8385`
(= `lambda_gk*psi_lambda_B/(beta_inter*Omega)` = 1.8031*3.0/(0.9975*6.4670)): the
fundamental expected loss is only **10.9%** of the total default loading; the GK
collateral friction is **89%**. Quantifies S-1 — the model is an amplification story,
not a credit-loss story. Consequence: `recovery_rate` 0.00→0.30 moves the loading by
3.3% partial / 6.0% in GE (peak spread 199.2→187.2bp ann).

**Units correction:** `spread_rb` is a *quarterly* rate deviation; the 150bp target is
annualised. Annualise with ×4×1e4 (`BP_ANN` in `run_regimes.py`). At `psi_lambda_B=3.0`
the model gives **187.2bp ann**, i.e. ~25% *over* target — not under.

**`audit_artifacts/` deleted (30 files).** `run_audit.py` carried its own hardcoded copy
of the calibration instead of importing `get_calibration()`, so it silently tested a
different model than `code/main.py` — its `ca_res_D` "regression failure" (1.479e-7 vs
1e-7) is reproduced exactly by the committed baseline artifact and predates any change
here. The `ms-regime` branch had already deleted it for the same reason (`2fa1b55`).
`code/main.py` is now the only regression path. Doc references updated; historical
provenance citations in finding write-ups left verbatim under a deprecation note.

**`diagnostics/regimes/` — two fixes and one open blocker.** `PSILAM_MAIN` was hardcoded
to `1.1793` while serving as **both** the guard and the **cache filename key**, so it
pointed at a stale cache built under a different model; it now reads the live
calibration. The `psi_lambda_B<1.5` guard was EBA-specific (thin net worth) and is now
`PSILAM_BREAKDOWN=4.0` for pre-EBA net worth. Open: `gamma_for_compression` scans
`linspace(0,60,61)` and requires *global* monotonicity; at the current calibration
`peak(γ)` is monotone on `[0,25]` and ticks up 1.1bp at γ=30 (81% compression,
saturation), so the run still aborts. Targets sit at γ≈1.6 (25%) and γ≈5.1 (50%) — far
inside the monotone region; narrowing `hi` to 25 would fix it (bisection only needs
monotonicity on its bracketing interval). **`run_regimes.py:75` hardcodes `"mv_rule=1"`
in its JSON provenance string — now false.** Note the regimes feature has *no*
structural `mv_rule` dependence; `A_cb=-2.406e-2` (compression) holds at `mv_rule=0`, so
SA-1 does not recur.

**`docs/PROCESS.md` deleted** — superseded by PROGRESS.md at `0c99013`; dangling
references in CLAUDE.md/HANDOFF.md/PROGRESS.md fixed. The doc-sync pre-commit hook is
documented but **not installed** in `.git/hooks`.

---

## 2026-07-24 — phi_lamb doc correction, Finding F-2 (ring), doc-sync hook  [803ecd2, +this commit]
- **Corrected a stale `phi_lamb`**: the committed value is **0.60** (~Bohn), not the
  `0.15` STATE.md's calibration table/IRF header claimed nor the "0.30" in
  `calibration.py`'s comment. STATE.md line 37 and the code itself already had 0.60
  (the run's `ρ_b=0.373` inverts to 0.60). `phi_lamb` governs the debt/fiscal mode
  only — well above the F-1 near-unit-root zone.
- **Finding F-2 (financial-accelerator ring).** The damped ~25q oscillation in the
  asset-price / bank-net-worth IRFs (`n_inter_D`, `q_b_D`, `q_b_F`, `C_D`; absent
  from `Y_D`/`w_D`) is the GK IC/leverage accelerator: **|λ|=0.954, period ≈25q
  (~6.25 yr), half-life ≈14.6q (~3.6 yr), R²=1.0.** A 100× PAC sweep
  (`psi_bF_D=psi_bD_F`, 0.05→5.0) moves |λ| by <0.006 and leaves the period fixed;
  committed `PAC=0.5` sits at the modulus *minimum*. The ring is intrinsic to the
  amplification block (pinned to the 150bp spread target) → **no free cosmetic
  fix**; a structural financial cycle to describe, not patch. New:
  `audit_artifacts/pac_sweep.py`, `pac_sweep_results.json`.
- Re-ran `code/main.py` end-to-end (clean: SS `goods_mkt_D≈-4.8e-7`,
  `ca_res_D≈2.3e-16`; `b_gov_D[499]≈2e-6`). Confirmed default-shock
  `n_inter_D[0]=-2.83%`, peak spread `+0.392pp`, TPI loading 3.59/3.03/2.47 and
  spread compression 0.392→0.244pp (−38%) over γ=0→10.
- **Added the pre-commit doc-sync hook** (`.claude/settings.json`,
  `.claude/hooks/require-docs-before-commit.sh`) and **this changelog**
  (`PROGRESS.md`); the hook's required-doc set is `STATE.md / PROGRESS.md /
  HANDOFF.md`.

## 2026-07-23 — Policy regimes (exogenous backstop aggressiveness)  [55f2031, 9e123a8]
- Rebuilt the policy-regime feature on `main`: three exogenous ECB backstop regimes
  over the TPI feedback coefficient γ (`TPI_t = γ·(spread−spread_ss)`), plus a
  Stage-B belief lottery over an ex-ante-unknown CB type. All post-Jacobian numpy on
  cached `G_tpi`; production `main.py` untouched (code in `diagnostics/regimes/`).
- **The capital-key backstop compresses spreads** (`d(spread_rb)/d(cb_buy)=−1.95e-2`);
  aggressive/medium γ hit 78/117 bp peak vs 156 bp passive.
- Provenance: the opposite "SA-1 spread-*widening*" result from the retired
  `ms-regime` branch was an artifact of a superseded model (single-country conduit,
  par rule, `psi_lambda_B=2.8` in the breakdown region) — **void on `main`.**

## 2026-07-22 — C-1 fixed at root; psi_lambda_B → 1.1793; EL-1/PT-1 resolved  [53e1783, e73659d, af02685, 0c70882]
- **C-1 fixed at its root** [53e1783]: `steady_auxilliary_D/F` now solve `lambda_gk`
  from the multi-asset IC directly; `Delta_bD_D/F=0.2/0.4` are genuine hardcoded
  inputs (was a degenerate back-solve, `Δ_cross=1.45>1`). The "explosive EBA doom
  loop" was superseded, not merely improved — `b_gov_D[499]`≈1e-5 near-stationary.
  Fixed 3 downstream drift bugs found while verifying (stale audit harness, a
  diagnostic sign/scale error, a real TPI conduit accounting leak). Re-tested
  Finding F-1 with a validated **order-selected Prony** eigenvalue estimator (the
  earlier energy-ratio proxy was itself an overfitting artifact).
- **psi_lambda_B recalibrated** [e73659d]: 0.31 → **1.1284** to hit the 2010 GR-DE
  ~150 bp spread target (0.31 was only a bug-dodge, undershooting ~3×). Found a
  non-monotonic linear-approximation breakdown region above ψ≈2 — the old literature
  values 2.8/3.0 now sit inside it and must not be restored.
- **EL-1 + PT-1 resolved** [0c70882]: `recovery_rate` 0.00 → **0.30** (NPV Greek-PSI
  framing; Zettelmeyer–Trebesch–Gulati). `psi_lambda_B` re-tuned 1.1284 → **1.1793**
  to restore the 150 bp target after `EL_price` shrank. Pass-through (PT-1) validated
  at ≈−4.5%/100 bp, inside the Acharya–Drechsler–Schnabl range. **TPI loading rose to
  3.59/3.03/2.47** at γ=2/5/10 (over-compensated, declining — self-extinguishing
  premium holds).
- Retired `docs/FRAMING_HANDOFF.md` into SPEC/STATE/HANDOFF/CLAUDE [af02685].

## 2026-07-21 — EBA-2011 calibration goes live  [eade414, e003284, f43f3a7, dd704eb, 1ce5312]
- Ported EBA 2011 stress-test bank-sovereign concentration (`phi_bD_D=2.39`,
  `phi_bF_F=2.76`; thin bank net worth `n_inter_D=0.408`) and added a passive
  capital-intermediation fund (`omega_K_D=0.0601`, `omega_K_F=0.0190`) to reconcile
  thin EBA net worth with the aggregate capital stock (Walras-neutral; no-op at
  `omega_K=1`).
- Correct-signed doom loop achieved at `psi_lambda_B=0.31`, `mv_rule=1`,
  `phi_lamb=0.6` [e003284]. Documented the endogenous-`omega_K` negative result and
  the C-1 root cause analytically (ruling out EBA leverage as `theta`).

## 2026-07-13–17 — Fundamental EL channel, ECB conduit, EBA spec  [c6f5707, 995a957, 1e82e22, 738311d, abcbb6e, b25cd1f, c981f19, a6961ab, 7d560ed, a129238]
- Added the fundamental expected-loss channel (`EL_price`) to the sovereign bond FOC
  (macro-pru-fix) [c6f5707], independent of `psi_lambda_B`; disciplined `psi_lambda_B`
  by a moment sweep and a `psi_lambda_B=0` output experiment [995a957, 1e82e22,
  738311d].
- **ECB balance sheet as a capital-key conduit** [abcbb6e]: CB net cash flow split
  `kappa_cb_F=0.929` to the F treasury / rest to D; two-leg (carry + credit) P&L;
  off-path expected loss hand-computed (never read off the linear path). Figure 8
  insurance-loading schedule + captions baked into all TPI figures [b25cd1f]
  (loading then 4.86→4.06→3.22 at the ψ=0.31 calibration).
- EBA-calibration design spec + a moments module decoding the 2011 stress-test CSV
  into calibration targets [a6961ab, c981f19]; 2010 base-year decision and the
  `omega_K` fork resolution [7d560ed].
- **Retired `model_v12.ipynb`; designated `code/main.py` as the production
  pipeline** [a129238].

## 2026-06-23–29 — Market-value fiscal rule (Finding F-1)  [d4fa259, f25dcf6]
- Added the switchable **market-value fiscal rule** (`mv_rule`): the Bohn rule can
  react to the mark-to-market debt gap `q_b·b_gov(-1) − mv_gov_ss` instead of par —
  it "sees" the current spread, restoring stationarity with empirical long-duration
  bonds under the risk-premium framing. Documented as Finding F-1 across nav docs.

## 2026-06-22 — Modular reorganization (PR #28)  [7037d72, ab52722, e6affd0]
- Reorganized the monolithic `model_v12.ipynb` into modular Python files
  (`code/equations_{D,F,global}.py`, `main.py`, `steady_state.py`,
  `ic_delta_calibration.py`, `depreciation_calibration.py`, `full_model.py`,
  `tpi.py`, plot modules); removed stale notebooks/folders. `main.py` becomes the
  single source of truth. Fixed stale branch/PR references after the audit merge.

## 2026-06-10–12 — Forensic audit + six structural fixes (PR #27)  [396cbd9, 4c810e1, 52f17d5, a0ddc18, 1e68440]
- Hostile-referee forensic audit: the model converged but was **internally
  inconsistent**. Six fixes, each verified numerically (`docs/audit.md`):
  - **T-2** (critical) — deposit rate re-dated to `Rgross=(1+rdep(-1))·P(-1)/P` (was a
    period-t unknown on the t−1 deposit stock, ~9× quarterly GDP). The doom loop had
    run *backwards*: bank net worth and GDP *rose* after a sovereign-default shock.
    Required `phi_lamb` 0.02→0.15.
  - **W-2/W-3** — p-conversion of F-bank's D-good-denominated bond book
    (`bank_return_F`, `divert_bond_foc_F`); `goods_mkt_F` had leaked up to **2% of
    GDP**, corrupting all cross-country spillovers by 39–124%.
  - **W-1** — capital-timing: keep `Y=F(K_t)`, pay `mpk·ΔK` to the capital producer so
    **CA=ΔNFA** holds at first order [4c810e1].
  - **TPI-1** — added the CB budget constraint (remittance); unbacked closed-loop flows
    had inflated welfare gains ~40% at γ=10.
  - **A-2** — aligned `m` vs Φ/T bookkeeping between SS and dynamics (needed before any
    `chi1≠0` experiment).
- Acceptance (**passed**): max|ca_res_D|,|goods_mkt_F| ≤1e-7 on all shocks incl. TPI;
  `n_inter_D` and `Y_D` fall on a default shock; system stationary. Post-fix
  notebook-era TPI: ΔW_D +1.88 / ΔW_F −1.90 at γ=10 (later superseded). Left open at
  audit: **S-1** (writeoff regime), **C-1** (`Δ_cross>1`).

## 2026-06-03–10 — TPI experiment, Walras-leak fixes, bond-pricing inversion  [f810776, 6aaa9ec, 25caad4, 754ea86, ecbae65, 09cf30e]
- Added the **TPI closed-loop spread rule** with welfare figures (PR #25) [754ea86];
  endogenous default rule + lump-sum tax fiscal rule [ecbae65].
- Fixed a Walras leak by adding `cap_profit` to the resource constraint [f810776];
  **inverted bond pricing** so `q_b` is the forward-looking PV and `rb` the implied
  YTM [6aaa9ec]. Precomputed the Jacobian `G` and got IRFs via matrix multiply
  [25caad4]. Added `spread_rb`/welfare aggregates; "sovereign risk without actual
  write-off" via `writeoff_enabled` [09cf30e].
- Discovered and reverted a "perverse GDP rise under default" fix [985b0b6→670c0ce] —
  the precursor symptom later diagnosed as T-2.

## 2026-05-21–28 — Multi-asset GK IC, endogenous default, macroprudential charge  [353c053, 5d7227a, a0c0828, 19df483, 26a223d, afd3bf8, b15811d, 6214e84]
- Built the **multi-asset Gertler-Karadi incentive constraint** (three-ν,
  asset-specific diversion rates) and a canonical-GK steady state (Ω, `f` convention,
  `theta` retune, `lambda_gk` calibration).
- Added **endogenous default** via a smooth-power debt-gap response (`def_scale`), an
  endogenous sovereign risk weight (lagged default + rating shock), and a Basel-style
  `mp_wedge` sovereign capital charge. Merged several topic branches (PRs #5–14).
- Forked `model_v12` as the no-central-bank baseline.

## 2026-05-14–20 — Foundations  [4aa2155, PR #1–3, 2327f4c, 19df483, 9385c05]
- Initial repo, README, and two-country model assembled from merged predecessor
  notebooks (v10/v11).
- Government default (PR #2); bond pricing; cross-border excess-return anchors; a
  spread-miscalibration fix [9385c05] (had targeted `rdep_F−rdep_D` not `rb−rdep`,
  linearizing around a non-stationary point → ~60/12 bp errors).
- nbstripout/nbdime notebook hygiene (PR #3); `zeta_writeoff` dial; endogenous-default
  channel via `def_scale`.

---

## Predecessor models (pre-git / `OLD models/`)
- `model_v11.ipynb` — predecessor with free bond trade between intermediaries.
- `model_v12` — added NK labour, portfolio adjustment costs, and the TPI extension;
  reorganized into `code/main.py` (PR #28) and then removed.

*(This section folded in from the former `docs/PROCESS.md` "Version history".)*

---

**Sources:** `git log`, `docs/STATE.md`, `docs/audit.md`, `docs/eba_calibration.md`,
`docs/HANDOFF.md`.
