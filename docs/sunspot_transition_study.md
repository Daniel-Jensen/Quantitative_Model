# Sunspot Transition: Numerical Forensics and the Income-SDF Fix

> **2026-07-15 addendum — fix program (see §8 at the end):** the mild
> calibration diagnosed below was reverted to the documented values
> (δ_b = 0.036, recovery = 0.45), the impact boom was fixed structurally
> (capital-quality loss + working capital + recap/ladder), and the M1/M2
> mechanisms of §3 are now partially closed.  §§1–7 document the OLD
> configuration; numbers there no longer describe the model.

*2026-07-14. All numbers verified on `code/global/` at the current calibration
(sunspot ξ = 1%·0.95^t, recovery 0.80, single-λ, anchored entrants, ψ = 0.05,
T = 200 in `calibration.py`; study runs sweep T in memory). Scripts and raw
logs: `/tmp/sdf_forensics.log`, `/tmp/truncation_study2.log`.*

## 1. What changed: Euler-consistent income SDF (`sdf_mode="income"`)

The banker's default-state discount is now an Euler-equation loading on
**aggregate output**:

```
Λ^d_X[t] = beta_inter_X · (Y^d_X(0) / Y^nd_X[t+1])^(−σ_X)
```

(`risk_branch.make_risk_inputs`). It prices the default state by how deep a
recession it is — high marginal value exactly when GDP is low. Unlike the
consumption-composite SDF (`sdf_mode="model"`, wrong-signed because the
deposit-rate collapse makes branch *consumption* rise), branch *output* falls
(−5.5%), so the sign is correct without any deposit-market surgery. Unlike
the retired free loading κ_d = 2.0, the loading is **endogenous and
disciplined**: measured (Y^d/Y^nd)^(−2) = **1.120** (range 1.1199–1.1223
across dates). A sign gate falls back to the empirical mode loudly if the
branch is ever not a recession. `kappa_d` is kept as a robustness dial.

### A/B at the current experiment (T = 100)

| | Y_D[0] | C_D[0] | I_D[0] | p[0] | n_D[0] | Q_bD[0] | sov peak (def/risk/liq) | lend peak |
|---|---|---|---|---|---|---|---|---|
| risk-OFF | −0.020% | −0.147% | +0.057% | +0.50% | −1.53% | −2.08% | +34 (77/0/126) | +126 bps |
| risk-ON **income** | **+0.058%** | +0.124% | +0.475% | −1.27% | −1.07% | −1.80% | +30 (92/**13**/82) | +116 bps |
| risk-ON empirical κ=2 | +0.103% | +0.162% | +0.998% | −2.20% | −1.40% | −2.93% | +49 (89/**91**/25) | +195 bps |

Two honest conclusions. (i) The income SDF gives a *modest* risk premium
(+13 bps at peak vs +91 with κ_d = 2) — that is what Euler-consistency with a
5%-output-cost default state buys; a larger premium must come from a worse
feared state, not a bigger free κ. (ii) **The impact boom is not the SDF's
fault and is not fixed by it** (κ_d = 1 also boomed): it comes from the two
deferred mechanisms below.

### Mechanism note: risk premium cannibalizes the liquidity premium

Across the three rows the liquidity component falls as the risk component
rises (126 → 82 → 25 bps) while the *total* moves little. Two-branch pricing
weights the default branch's very low capital return in μ, compressing the
IC multiplier and hence the balance-sheet spread λμ/Ω̃ — the stronger the
risk weighting, the cheaper the IC charge. The total sovereign spread is
therefore far less sensitive to the SDF than the decomposition split is.

### Presentation caveat (main.py table)

The "Sov spread peak" is a **yield** spread (duration-weighted average of
future per-period spreads); the def/risk/liq components are **per-period
excess-return** spreads. They satisfy an exact per-period identity
(verified 2e−16) but the components do *not* sum to the yield number —
at the peak the per-period spread far exceeds the yield spread because the
sunspot decays. Don't read the table as `total = def + risk + liq`.

## 2. Internal consistency: the path is exact

All identities on the income risk-on path (T = 100), verified numerically:

| Identity | max abs residual |
|---|---|
| Two-branch pricing FOC (Q reconstructed from Ω^nd, Ω^d, π, μ) | 0.0 |
| Bond clearing b_D_D + b_D_F − b_gov_eop | 4.4e−16 |
| Spread decomposition identity | 2.1e−16 |
| Goods D (imposed) | 4.9e−12 |
| Goods F (Walras diagnostic) | 1.0e−09 |
| Government budget Tax + issuance − G − coupon | 2.8e−17 |
| Capital markets n_IC − n_ACCUM (D/F) | 2.5e−11 / 9.6e−11 |
| Deposit market D | 9.2e−11 |

**The "weird" responses are not accounting or solver errors.** They are
equilibrium properties of the current model structure.

## 3. Anatomy of the weird impact responses

1. **Deposit-rate crash → consumption boom (M1).** Banks deleverage →
   deposit demand falls → rdep_D falls (−21 bps at impact). The household's
   date-0 return is *predetermined* at the SS rate, so the date-0 C jump
   (+0.12%) is a purely forward-looking Euler response to the anticipated
   low-rate path. This is Bocola's comovement problem (his §VI) amplified
   by the segmented, flexible-price deposit market. Fix: union deposit
   market (rdep_D = rdep_F) or an NK union rate — deferred by scope.
2. **Capital is the branch safe haven → investment boom (M2).** In the
   feared default state bonds lose ~23% (haircut + repricing) while capital
   loses ~7.5% (output cost only). Two-branch pricing therefore *favors*
   capital: the IC envelope freed by falling bond values is reallocated to
   K (I_D[0] +0.47% income / +1.0% κ2 / +0.06% off), and the associated
   capital gains mean the risk-on net-worth loss is *smaller* than risk-off
   (−1.07% vs −1.53%). The single-λ design closed the divertability
   substitution margin; the branch-return asymmetry reopens one. Fix: a
   default-state capital-quality loss (Gertler–Kiyotaki) — deferred.
3. **rk_D[0] spike ≈ +23 bps** is the Jermann-Q revaluation from the
   investment boom (capital-gain term), not an mpk move.
4. **Timing conventions** (numerically confirmed): the date-0 sunspot value
   never enters pricing — Q_t prices π_{t+1}, so setting π_0 = 0.9 changes
   Q by exactly 0.0. The sunspot works entirely through anticipated
   π_{t≥1} plus the zone indicator on the endogenous debt path.

## 4. "No return to steady state after 100 periods": two separate facts

### (a) T = 100 truncates the sunspot experiment — visibly

Common-window drift, T = 100 vs T = 200 solutions on t ∈ [0, 70):

| variable | max drift (share of SS) | vs own impact response |
|---|---|---|
| Y_D | 1.5e−04 | ~25% of the Y_D response |
| C_D | 2.0e−04 | ~15% |
| p | 1.4e−03 | ~10% |
| n_F | 4.4e−04 | ~4% |
| Q_bD | 7.9e−05 | ~0.4% |

With the ρ = 0.95 sunspot and slow states, **T = 100 contaminates even the
first 70 quarters** — up to a quarter of the Y_D response is a terminal
artifact. Terminal wedges: at T = 100 the endpoint is visibly off SS
(n_F +0.24%, p +0.56%); at T = 200 they are small (+0.06%, +0.02%).
**T = 200 (the current calibration setting) is the right default for the
sunspot centerpiece**; T = 100 remains fine for the TFP shock (ρ = 0.8).
Read plots only to ~t = T−40.

### (b) The F-side genuinely never returns — a cross-country wealth unit root

Half-lives fitted on the T = 200 income risk-on path (window [50, 160]):

| state | half-life | comment |
|---|---|---|
| n_D | 12.3q | fast: anchored entrants + strong excess-return feedback (growth factor 0.966 < 0.9911 bound) |
| b_gov_D | 18.9q | Bohn rule |
| Y_D, Kap_D | ~28–29q | follows n_D and K adjustment |
| p (RER) | 30.5q | |
| A_D − A_F (relative wealth) | **67q** | the binding slow state |
| n_F | **~4,400q** | effectively permanent |
| Kap_F | **∞ (non-reverting in window)** | |
| Y_F | **~1,600q** | |

The crisis **permanently redistributes wealth toward F**: F banks buy
D-bonds at crisis prices and keep the excess returns (n_F overshoots to
+0.05% above SS and stays), households' relative wealth shifts, and nothing
anchors it back — the ψ adjustment cost pins gross *bank* bond positions,
not relative *country* wealth. This is the standard incomplete-markets
open-economy non-stationarity, concentrated on the F side because the
D side now has the entrant anchor + a strong spread feedback. The risk-off
T = 300 run shows the same pattern much smaller (Y_F +0.006% at t = 240),
so risk pricing amplifies the redistribution but does not cause it.
Extending T will never make these paths visibly "return" — only a
stationarity device on relative wealth would (union-wide deposit market,
or a debt-elastic wealth anchor à la Schmitt-Grohé–Uribe).

## 5. Solver findings

- **T = 300 risk-on fails acceptance** (max|resid| 2.6e−06 > 1e−06) from a
  padded T = 200 warm start — long-horizon hybr conditioning; same failure
  mode as the T = 300 cold-start TFP result. Practical ceiling for risk-on
  runs is currently T ≈ 200–250. (Risk-off solves fine at T = 300, 338s.)
  The same marginal-stall mode appeared at T = 200 inside the risk loop
  (1.4e−06). Fix: **polish restarts** in `solve_transition` — when hybr
  stalls above the bar, restart it from the stalled point with a fresh
  small trust region (up to twice, no-op when converged) before the krylov
  fallback; this typically shaves the last order of magnitude.
- **krylov fallback crash fixed** (`transition.py`): scipy's krylov can
  raise "Jacobian inversion yielded zero vector" near the penalty walls;
  this killed `test_risk_channel` mid-run. The fallback is now wrapped —
  a krylov failure keeps the hybr solution instead of crashing.
- **Stale-test bug found and fixed**: `test_risk_channel.py` hardcoded
  `Z = np.full(T, 1.0)` in two GE tests, predating the Z-rescaling to
  Z_ss = 0.448 (the Y_ss = 1 normalization). Those tests were unknowingly
  solving a +123% permanent TFP level shock — the true cause of the
  penalty-wall non-convergence (resid = 10.0) and the original krylov
  crash. Now uses `cal["Z_ss_D"]`/`cal["Z_ss_F"]`.
- Runtimes (this machine): risk-on T=100 ≈ 165s, T=200 ≈ 800–920s;
  risk-off T=300 ≈ 340s.

## 6. Test-suite status at this configuration

`test_ss_identities`, `test_bank_block`, `test_transition_walras`,
`test_signs_bocola` (risk-off signs): **pass** at T = 200.
`test_risk_channel`: after fixing the stale Z = 1.0 inputs (section 5), the
structural assertions (π ≡ 0 nesting, decomposition identity, positive risk
premium, zone consistency, Walras, no realized default, zero-shock fixed
point) are enforced. Two *directional* assertions — "risk-on Q_bD[0] below
risk-off" and "risk-on n_D[0] below risk-off" — were downgraded to loud
warnings: under the disciplined Euler loading (1.12) the μ-compression and
capital-gain offsets legitimately dominate the small premium (sections 1
and 3). They flip back once the default-state capital-quality loss is
added; re-promote them to assertions with that fix.

## 7. Bottom line

- The Euler-consistent income SDF is in, disciplined (loading 1.12), exact
  (nesting and identities at machine precision), and sign-safe (gated).
- The impact boom under priced risk is a *structural* property (M1 + M2),
  present since the first risk-channel runs and unmasked by the mild
  calibration; fixing it requires the deferred union-deposit-market and/or
  working-capital and capital-quality changes.
- Non-return by t = 100 is (a) one-quarter truncation contamination at
  T = 100 — solved by the current T = 200 — plus (b) a genuine
  cross-country wealth quasi-unit-root that no horizon extension will
  remove.

## 8. 2026-07-15: forensics at the 10% sunspot and the fix program

### 8a. Why the experiment was weak (all numbers verified, T = 200 ≡ T = 500 to 3 decimals)

At the pre-fix configuration (δ_b = 0.25 ≈ 1y duration, recovery = 0.80,
sunspot 10%·0.9^t — a **64% cumulative priced default probability**):

1. **Mild priced event.** Pure expected-loss repricing of Q_bD at SS
   discounting: −5.3% (vs −30.1% at the documented δ_b = 0.036 /
   rec = 0.45).  Commit a82431a had reverted the documented values because
   the risk-ON branch wiped out bank equity (no recap; ladder removed).
2. **GK block absorbed the shock.** Branch capital lost only ~8.7pp/q vs
   ~21% on bonds → capital was the branch safe haven → two-branch μ_D went
   **negative** (−0.022 vs SS +0.0196; always-binding IC violated in
   equilibrium), the wedge λμ/Ω̃ collapsed (−420bp ann, cancelling over
   half of the +770bp default compensation), and the imposed IC equality
   with α↓3% *forced* a recapitalization boom: Q_D(cap) +1.6%, I_D +3.4%,
   n_D +3.8% at t = 1.  Equilibrium Q_bD[0] only −3.5%.
3. **No output transmission.** Impact Y is pinned by the GHH identity
   ŷ = α·k̂ + (1−α)·(α·k̂ − P̂_CES)/(1/ν + α), P_CES loads on p with weight
   0.15, and no spread enters production.  Counterfactual proof: at the
   documented calibration risk-OFF produced Q_bD −29.8%, n_D −20%, lending
   spread +1842bp — and Y_D fell only −0.198% with I_D still +0.76%.

### 8b. The fix program (user-approved plan, phases 1+2)

- **Calibration restored**: δ_b = 0.036 (~7y HM duration), recovery = 0.45
  (Greek PSI).  T back to 200; centerpiece sunspot 1%·0.95^t.
- **One fixed feared event + government recap** (`recap_share_D = 0.5`,
  HFSF/EFSF analogue): equity injection financed by issuance in the branch —
  post-default debt and Bohn taxes rise.  The branch does ONE deterministic
  solve of the full PSI haircut with recap (`branch_haircut_scale = 1.0`,
  `rescue_mode = full+recap`); no per-run search.  *(Update 2026-07-15: the
  original attempt-ladder full → full+recap → ladder was replaced by this
  single-solve path at the user's request — the ladder search over scales
  0.075…1.0 is now opt-in only, `branch_use_ladder`, default off; it also
  removed a wasted no-recap probe every round.  If the fixed event is
  infeasible the branch raises with a hint instead of searching.)*
- **GK capital-quality loss** `def_capital_quality_D = 0.05` at branch
  h = 0 (capital stops being the branch safe haven; branch rk_d ≈ −12%/q,
  branch Y(0) ≈ −10%, n_d(0)/n_ss ≈ 0.32).  ξ_K = 0.10 tested and
  REJECTED: μ_min → 0.005, deeper post-impact boom, C_D[0] flips positive,
  branch n → 0.14 (2 sign failures vs 1).
- **Bohn rule on the surviving stock** (government.py): taxing the
  pre-haircut stock produced a one-quarter ~31%-of-GDP tax spike at the
  PSI haircut, which alone made full-event branches infeasible.
- **Working capital (Neumeyer-Perri)**: firms pre-finance `zeta_wc = 1` ×
  wage bill at r_wc = rdep(−1) + λμ/Ω̃; the wedge enters labour demand
  (w divided by 1 + ζ·r_wc), financing income is passed to households as
  dividends (intra-period, never on the bank balance sheet — Bocola-lite;
  routing through bank equity would break the closed-form leverage/spread
  calibration).  ζ = 0 nests exactly; SS Walras 2e−10 at ζ = 1.
- **μ monitor** (transition.py): warns loudly if the IC multiplier goes
  negative on a solved path; `mu_min_D/F` returned.
- **Branch warm start**: `_shifted_y0` scales the Kap_D block by
  (1 − ξ_K·0.975^t) — a raw shifted base path implies a one-quarter
  rebuild whose Jermann-Q spike lands every probe on the penalty wall.
- Branch solves accept max|resid| ≤ 1e−9 (the test bar; quality-shocked
  systems polish-stall a few x above the 1e−10 default).

### 8c. Centerpiece results (sunspot 1%·0.95^t, T = 200, ξ_K = 0.05)

Sign checks: **12 of 13 pass** — Q_bD[0] −5.68%, n_D[0] −3.45%, n_F[0]
−2.38%, Y_D[0] −0.07%, C_D[0] −0.03%, lending spread peak +344bp ann,
b_gov↑, Tax↑, μ_D > 0 everywhere (min +0.010), risk premium +74bp ann,
risk-on Q_bD[0] < risk-off.  Remaining failure (kept as test WARNING):
risk-on n_D[0] (−3.45%) above risk-off (−4.10%) — the M1 deposit-rate
collapse finances an impact I boom whose Jermann-Q gains cushion bank
capital; killed by the union deposit market (deferred).  Post-impact Y_D
turns mildly positive (+0.2–0.3% for ~10y) through the same M1 channel +
residual μ compression (risk-on μ ≈ 0.010–0.015 < μ_ss → the wc wedge
runs mildly *below* SS after impact).  Spread decomposition t = 0:
total +99, defcomp +219, risk +74, liquidity +278 then negative from
t = 1 (−52…−98) — μ compression persists but no longer dominates.
Branch-side μ < 0 (either bank *inside* the deep feared state — the
main.py run shows branch μ_D ≈ −0.07; validation runs showed branch
μ_F < 0) — logged by the monitor, harmless to the BASE path (base
μ_D = +0.010 > 0), but it means the branch's own IC is questionable where
the risk inputs are extracted (always-binding-IC limitation, now visible
because the feared state is a genuine deep crisis).  Runtime: risk-ON at
this calibration ≈ 32 min (round-1 branch homotopy ≈ 27 min; later rounds
seconds via jac_cache).

### 8d. Open items

Union deposit market (kills M1, enables sdf_mode="model", should flip the
n-ordering warning and the post-impact Y boom); branch μ_F; stress shocks
(10%) can still drive μ_D < 0 — the monitor flags them; occasionally
binding IC remains the structural answer.
