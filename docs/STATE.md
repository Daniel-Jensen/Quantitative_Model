# Project State

**Branch:** `file-reorganisation` | **Date:** 2026-07-07 | **Status:** Bocola (2016) / Cole-Kehoe sovereign-risk mechanism implemented and verified; **risk channel added** via two-branch default-branch pricing (standalone `code/global/` model)

## Risk channel (added 2026-07-07, `risk_branch.py`)

Bocola's second transmission channel — precautionary deleveraging — is now
implemented. Bankers discount with the household SDF Λ = β·u_c′/u_c (Bocola
uses log utility, **not** Epstein-Zin; verified from the paper) and weight a
**representative post-default branch** by the priced default probability in
every backward-pass expectation (Ω̃, μ, α, bond prices, cross-border FOCs):

    Ω̃_{t+1} = (1−π)Ω^nd + π·Ω^d,   Ω^d = Λ^d·[(1−f)+f·α^d(0)]

The branch is a full PF transition with the haircut realized at its period 0,
launched from the base-path impact state (new `init=` support in
solve_transition — also the scaffolding for future mid-crisis TPI runs) and
absorbing (post-default debt exits the CK crisis zone). Outer fixed point
base ↔ branch, damped, warm-started. `pi ≡ 0` nests the risk-neutral model
exactly (regression-tested). Diagnostics: `bond_decomposition` splits the
sovereign spread into default compensation + **risk premium** + liquidity
premium via an exact per-period identity. Approximations (documented):
single representative branch, Λ^nd ≡ β_inter, aggregate-composite SDF as the
rep-agent proxy for the HA household, and household-side π-blindness — the
deposit Euler never weights the default branch by π (deposits are safe in
rate terms even in the branch, so what is omitted is only the
precautionary-savings response to default-state income risk; risk pricing
lives entirely in the bank block). Validation moment: risk-channel share
of the lending-spread response vs Bocola's "up to 45%".

Motivation (analysis of the flexible-rate caveat): without risk pricing,
risk-neutral banks re-lever into capital when the deposit rate collapses in
the crisis — the model's investment boomed. This is Bocola's own
"comovement problem" (his §VI); per author decision the comovement fixes
(union deposit market, working capital) are deferred, while the risk channel
disciplines banks' expansion through the covariance premium on capital.

**Default-state specification** (required for a correctly-signed premium —
discovered numerically): with the plain Bohn rule, a default's 55% haircut
became φ·(b−b_ss) ≈ −1.05 of *tax cuts* per quarter — default was
expansionary for households (branch Y(0) = +1.1%) and the risk premium came
out ≈ 0. Two canonical ingredients fix the default state:
1. **Fiscal re-anchoring** (`b_anchor` in govt_transition): post-default the
   Bohn rule anchors to the post-haircut stock — debt relief is not handed
   to households as windfall transfers (Greek post-PSI reality).
2. **Output cost of default** (Arellano 2008 tradition): branch TFP =
   Z·(1 − 0.05·0.9^h) — conservative next to the Greek 2012 collapse
   (`def_output_cost_D`, `def_output_rho_D` in calibration.py).

## Validity of the default-pricing mechanism in the two-country HANK setting

Reviewed 2026-07-08 (logic review vs Bocola 2016; FOC algebra, priced/realized
split, HM haircut conventions, timing and p-conversions all verified correct).
Assessment of why the Bocola pricing mechanism remains valid when embedded in
a two-country heterogeneous-agent economy:

1. **Pricing is bank-side by construction — no HA aggregation problem
   contaminates Q.** Households hold only safe one-asset deposits
   (household.py: non-contingent rate, locked one period ahead); they never
   hold sovereign bonds. Every pricing equation lives in the bank block,
   where the pricer is a representative banker with a well-defined
   objective. Heterogeneity reaches Q only through equilibrium objects, all
   inside the residual system: rdep (deposit-market clearing against the HA
   wealth distribution → the discount in Q), rk (MPC-weighted goods demand →
   MPK → μ → the liquidity spread), and the debt path (Bohn taxes →
   consumption → output → the crisis-zone indicator). This also matches the
   bank-centric holding structure of euro-area periphery debt.

2. **Banker SDF is an assumption, not an approximation.** GK/Bocola's
   "bankers inside a representative family" has no HANK analogue — there is
   no single household SDF to inherit (the constrained household's and the
   wealthy saver's differ enormously, and any aggregate-composite proxy is
   an arbitrary aggregation rule, additionally wrong-signed here via the
   comovement problem). Λ^d = β_inter·κ_d is therefore a banker-specific
   discount with an empirically disciplined default-state loading — a
   stated assumption (see calibration.py), disciplined by the ≈45%
   risk-channel share target.

3. **Household π-blindness is bounded.** Deposits stay risk-free in rate
   terms even in the default branch (the haircut feasibility ladder exists
   precisely to rule out equity wipeout / deposit impairment), so what the
   household Euler misses by not weighting the branch is only default-state
   *income* risk (wages, dividends, Bohn taxes). Second-order for pricing
   (π enters Q only via rdep); potentially first-order for distributional /
   welfare statements — paper caveat, not a pricing defect. A consistent fix
   (HA problem under two-branch expectations) is a research extension.

4. **What HANK genuinely adds is correctly wired and priced.** The fiscal
   amplification loop is MPC-weighted: depressed Q → rollover at low prices
   → debt ↑ → lump-sum Bohn tax ↑ → constrained households cut C
   one-for-one → Y ↓ → feeds rk, μ and the zone indicator back into Q,
   entirely inside the fixed point. Dividend-cut incidence is uniform
   per-capita (stated assumption; mildly amplifying through constrained
   households). The default branch launches from the actual HA distribution
   snapshot (`extract_init_state` passes `D_start`), so α^d(0) is
   distribution-consistent; representative-branch reuse ignores only
   distribution drift across pricing dates (second-order, documented).

5. **Residual design notes.** CK zone thresholds condition on b/Y_ss, not
   current Y (avoids another fixed-point layer; understates zone-deepening
   in deep recessions). Cross-border risk sharing flows only through banks
   (households cannot hold foreign assets) — the right restriction for a
   bank-centric monetary-union crisis model. In risk mode, F's own default
   risk is priced risk-neutrally and independently of the D-event (survival
   factor on both F-bond branch payoffs; no F default branch — see
   bank_backward docstring).

---

## Current model (`code/global/`) — Bocola–Cole-Kehoe rework (2026-07-07)

The sovereign-default mechanism was rebuilt to follow Bocola (2016, JPE)
"The Pass-Through of Sovereign Risk" embedded in Cole-Kehoe crisis zones,
after verification showed the previous ψ_bd reduced form produced an
*expansionary* default-risk shock in general equilibrium.

### Structural changes

| Change | Where | Rationale |
|--------|-------|-----------|
| Single λ per bank (λ_K = λ_bD = λ_bF = 0.22) | `bank.py`, `calibration.py` | Bocola eq. (3): banker diverts a fraction of TOTAL assets. Asset-specific λ let banks substitute into capital when bond IC tightened → wrong GE sign. |
| ψ_bd·ξ IC-tightening removed | `bank.py` | Replaced by expected-haircut pricing: sunspot = priced default probability. |
| PRICED vs REALIZED default split (`def_price` / `def_real`) | `bank.py`, `government.py`, `transition.py` | Bocola experiment: news of default is priced (MTM losses) but default never happens. Realized-default variant = pass `def_real ≠ 0`. |
| Endogenous debt in bond-market clearing | `transition.py`, `government.py` | Debt is forward-integrated (Bohn tax inside the recursion) within every residual evaluation; banks hold the true end-of-period stock. Closes the Walras leak (pre-fix: 0.47% of Y_F per 5% debt deviation) and restores the issuance-absorption amplification. Removes the old CK/BD outer debt loops. |
| Split `solve_bank_paths` → `bank_backward` + `bank_forward` | `bank.py` | Prices come from marginal conditions only, so debt can be integrated between the passes. |
| SS external balance = current account (not NX=0) | `steady_state.py` | With 20% cross-border bond books, net foreign income ≠ 0 matters. |
| SS household income `/P_CES` conversion fix | `steady_state.py` | Was missing (invisible while p_ss=1 exactly); with any asymmetry it caused a 1.3e-4 goods wedge. |
| Symmetric SS enforced (δ_b_F = δ_b_D) | `calibration.py` | p is weakly identified by external balance (η=0.5 ⇒ NX ∝ p^½·(P_F^½C_F − P_D^½C_D)); asymmetric SS opens an O(1e-4) wedge. Asymmetries enter through shocks. |
| BD solver (`solve_transition_bd`), `verify_mechanism.py` deleted | `transition.py` | Superseded by the CK–Bocola design and `tests/`. Git history preserves them. |

### Calibration (current)

| Parameter | Value | Target / source |
|-----------|-------|-----------------|
| λ (single, both banks) | 0.22 | Leverage θ_ss = 4.45 (GK11/Bocola range 4–6) |
| B_gov_ss (both) | 12.80 | D-bank sovereign exposure Q·b/n = 0.89 (Bocola: GIPS domestic sov holdings ≈ 93% of bank equity, 2009); face debt = 93% of annual GDP |
| b_D_F_ss = b_F_D_ss | 2.56 | Foreign bank holds 20% of each bond supply (union contagion leg) |
| δ_b (both) | 0.036 | ~7y duration (GR/DE pre-crisis average maturity) |
| recovery_rate (both) | 0.45 | Haircut 0.55, Greek PSI 2012 (Zettelmeyer-Trebesch-Gulati; used by Bocola) |
| b_ck_low_D / b_ck_high_D | 3.0 / 6.0 | SS b/Y_ss = 3.72 sits inside the crisis zone; fundamental default unreachable in the risk-only experiment |
| φ_lamb (Bohn) | 0.15 | Sweep over {0.02, 0.05, 0.15} changes C/I shapes little (C crash is driven by the deposit-rate collapse, not taxes) |
| f, ω_ent, β_inter | 0.028, 0.002, 0.96 | Unchanged; jointly give rk_ss − rdep = 1.77% ann (= SS sovereign spread under single λ) |
| ψ_bF_D = ψ_bD_F | 0.01 | Cross-border portfolio adjustment cost |
| a_max, n_a | 300, 250 | Deposit demand ≈ 35.3 per country (A = Dep_supply) |

### Residuals (verified 2026-07-07, T=100)

| Check | Value |
|-------|-------|
| SS: IC resid, Bellman, bond-FOC identities | ≤ 1e-12 (machine) |
| SS: goods market (Y − C − I) | 8.8e-07 (household-grid floor) |
| Zero-shock transition: max deviation from SS | ≤ 1.4e-06; goods_F 6.1e-07 |
| TFP shock (1%, ρ=0.8): goods_D / goods_F | 2.0e-11 / 6.0e-07 |
| CK risk-only shock (ξ₀=7%, ρ=0.95): goods_D / goods_F | 4.8e-10 / 5.9e-07 |
| Bond FOC E[rb]−rdep = λμ/Ω along shocked paths | ≤ 1e-12 per period |

The old "W-G1 structural limitation" (goods_F ≈ 1.7e-2 on TFP shocks) is
**gone** — it was accounting error, not structure; the diagnostic now sits at
the household-grid floor (~6e-7) on all shocks including moving debt.

### Centerpiece experiment (main.py): CK sunspot, risk only

> **2026-07-15 update.** The centerpiece is now ξ₀ = 1%, ρ = 0.95 (small
> persistent shock; the numbers just below are the older ξ₀ = 7% run, kept
> for reference).  Since then the default branch gained a GK
> capital-quality loss (`def_capital_quality_D = 0.05`) and a contingent
> government recap, a Neumeyer-Perri **working-capital** wedge
> (`zeta_wc = 1`) was added as the spread→output channel, and the δ_b/rec
> interlude values (0.25/0.80) were reverted to the documented 0.036/0.45.
> At the 1% sunspot: Q_bD[0] −5.7%, n_D[0] −3.4%, n_F[0] −2.4%, lending
> spread +344bp ann, risk premium +74bp; 12/13 sign criteria pass (risk-on
> n_D[0] above risk-off and a mild post-impact Y boom survive via the M1
> deposit-rate channel — killed by the deferred union deposit market).
> A μ-monitor now warns if the always-binding IC is violated on a solved
> path.  Full detail: `docs/sunspot_transition_study.md` §8.

ξ₀ = 7% quarterly default probability, ρ = 0.95, priced in the crisis zone,
never realized. Solved with a 3-step homotopy (~30s). Results:

- Q_bD −32% on impact (MTM repricing of 7y bonds)
- n_D −23% (no default!), n_F −5.6% (contagion via 20% cross-holding)
- Sovereign yield spread +770bps ann at peak; lending spread +1931bps at impact; pass-through 0.33 at peak
- Y_D −0.33% trough; C_D −3.9% trough
- b_gov +5% (rollover at depressed prices), Bohn taxes +0.097 peak
- Banks recapitalize in ~8-12 quarters via ex-post excess bond returns
  (bought at 60% of par, repaid in full — the fiscal cost of the belief shock)

### Known limitations (deliberate, next phases)

1. **Flexible prices / no union nominal rate**: rdep collapses in the crisis
   (−350bps ann on impact), cushioning banks and pushing the contraction into
   consumption while investment rises (real-model crowding-in). The
   monetary-union nominal block is the next major layer and is required for
   the TPI application.
2. ~~**No risk channel**~~ — added 2026-07-07 (`risk_branch.py`); see the
   section at the top of this file.
3. **IC always binding** (Bocola's binds occasionally, μ≈0 in calm times).
   As of 2026-07-15 a monitor warns when μ goes negative on a solved path;
   the 1% centerpiece keeps μ_D > 0, but stress shocks (10%) still trip it
   and the branch-side μ_F can be negative — occasionally-binding IC is the
   structural fix.
4. CK zones use Y_ss (no output feedback into the crisis zone).
5. Debt/annual GDP = 93% (face) is below Greek crisis peaks; the binding
   anchor is bank exposure/net worth ≈ 0.9 (banks hold 80% of supply).

### Next priorities

1. Nominal rigidities + single union policy rate (kills the rdep escape
   valve; expected to flip investment response and deepen the recession).
2. TPI/asset-purchase experiment: central bank buys D-bonds in the crisis
   zone (Bocola's LTRO experiment is the template).
3. Realized-default comparison run (`def_real` = PSI event at a chosen date).
4. Optional: occasionally-binding IC; risk-channel proxy via exogenous SDF wedge.

---

## Historical: SSJ model (`code/model_v12.ipynb`, branch `audit`) — superseded

The SSJ-era state (six structural fixes W-1, W-2, W-3, T-2, A-2, TPI-1;
Walras forensics; TPI welfare results; C-1/S-1 open issues) is preserved in
`docs/audit.md`, `docs/verification_report.md`, `docs/walras_forensics.md`
and in this file's git history (pre-2026-07-07 versions). The
`/opt/anaconda3/envs/ssj` environment no longer exists; do not use the
`audit`/`bank-cal` branches for new work.
