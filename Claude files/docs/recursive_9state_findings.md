# Recursive solver — 9-state rework: findings

**Date:** 2026-08-26
**Scope:** diagnosis and repair of the sovereign-risk IRF in `code/global/`, plus the
structural rework requested on 2026-08-25 (endogenous D-sovereign market, union deposit
market, calibration on the 8 bp spread). Measurement record, not a re-audit.
**Branch:** `bocola-rewrite`, working tree.

> **SUPERSEDED IN PART, 2026-08-28.** The Bocola-replication audit and the collocation
> rework that followed overturn three claims below. (i) **The benchmark in §4 is wrong**:
> Table 5's −1.05/−1.44/−1.53 is a cumulated quarterly *growth* gap ×400 over an
> 8-quarter estimated shock sequence, whose output *level* equivalent is
> −0.26/−0.36/−0.38%; the like-for-like single-shock targets are −0.157% (his §V.C open
> economy) and −0.222% (his closed benchmark). (ii) **O-3's "the response is inside the
> solution's own error" does not survive**: re-solving instead of reading the μ=1 fit
> moves output by 8%, the closed-form labour FOC reproduces the solved hours to four
> decimals, and the degree-2 basis *over*-states the bond channel by 23% on a controlled
> test. The response was small for an economic reason: `r_wc = rdep + λμ/E[Ω]`, and in
> the symmetric union the deposit rate fell 44.8 bp/yr against a 57.6 bp/yr credit-spread
> rise, cancelling 78% of the wedge. (iii) **O-1's 8 bp discussion stands, but the
> conclusion that the constraint is "barely binding" reads differently against Bocola**:
> computed exactly from his solved coefficients, μ = 0 along his entire benchmark IRF
> except one quarter, and binds on 1.2% of his ergodic set — this model's constraint is
> *more* active than his, not less.
> The solver described here (damped time iteration) is no longer the solver; see
> `solver_recursive/collocation.py` and CLAUDE.md.

**Method:** every claim below is a number produced by a run of the shipped code. Where a
mechanism was hypothesised and then measured, both the hypothesis and its outcome are
recorded, including where the hypothesis was wrong.

---

## 1. Finding status table

| # | Finding | Status | Evidence | Confidence |
|---|---------|--------|----------|-----------|
| **B-1** | Risk IRF's "trough at q7 then flat recovery" was the collocation **box wall**, not a mechanism | Fixed | `B_D` pinned at exactly +8.00% (`b_band`) from q7 to q24. Unclipped, `Y_D` runs to −2.26% at q24 and is still falling. `dynamic_irf` now reports escapes instead of clipping silently. | High |
| **B-2** | Bohn rule did not stabilise debt: root 0.9929 (half-life 97 q) | Fixed | Implemented as unit-**elasticity** `Tax_ss·(B/B_ss)^1` giving `dTax/dB = 0.00303`, because `G_D = 0` leaves `Tax_ss = 0.00297` against `B_ss = 0.98`. Replaced by a **linear level** rule with `gamma_tau` solved from a target debt root (`debt_root = 0.93`). Raising the elasticity instead was tested and rejected: φ=15 swings taxes 0.29×–3.17× across the B band and to 6e−6 at the default node, and the aliasing moved the SS spread 129→250 bp. | High |
| **B-3** | Time iteration was reporting an **unconverged** rest point as the answer | Fixed | Binding mode is the franchise-value recursion, slope `β(1−f)R/(1−μ) = 0.961`, i.e. 0.990 per sweep at `damp=0.25` — **235 sweeps per decade**, not the 0.75 the code comment assumed. Joint stage read 205 bp at 100 sweeps, 331 at 200, settling at 288 only past ~490. Budget raised to 400/300/800; every stage now exits on the rule-change test. | High |
| **B-4** | Default branch was a **fiscal windfall**, making the feared event expansionary | Fixed | `point_map` used a fixed tax anchor in both regimes while `government.py` already re-anchored its branches ("else the haircut becomes a tax-cut windfall → default expansionary"). Measured: tax at the default node −3.44% of Y against +0.30% at d=0, a transfer worth 12.6× the SS tax level on impact. After re-anchoring, d=1 vs d=0: **Y −2.60%, C −1.54%, N −3.70%, n_D −32.11%**. | High |
| **B-5** | `euler_F` error 22× the D side, from deriving `W_F` off the union identity | Fixed | `corr(|euler_F|, |W_F gap|) = +0.973`. A 0.06% error in `W_F` became a 1.1% error in `euler_F`, because `C_F = W_F/P_CES + inc − A_F` is a 0.79 difference of ~8-sized terms. Replaced by carrying `V_dep = W_D − P_D`, so `W_D = P_D + V`, `W_F = P_F − V/p` and the identity holds by construction. `euler_F` mean −2.03 → **−3.86** (68×); `dep_clear` −3.32 → −4.90; `ALL` −3.31 → −3.71. | High |
| **B-6** | `sigma_s` used Bocola's **reported** parameter against his **effective** process | Fixed | `gh_nodes` uses numpy `hermegauss`, the probabilists' rule (nodes already in sd units); Bocola's `GaussHermite.m` is the physicists' rule under the same `σ·node` map, so his solved model behaves as if `σ = 0.63/√2 = 0.4455`. `calibration.py`'s own comment stated this on 2026-08-15; the value was never changed. Measured cost: unconditional sd of s 2.02 vs 1.43, ergodic E[p^d] 0.66% vs 0.27%, D-bond 0.9057 vs 0.9178, SS spread 32.8 vs 25.5 bp. | High |
| **B-7** | s-box coverage was hard-wired in **absolute** units, not sd | Fixed | `s_halfwidth = 4.35` is Bocola's ±2.16 unconditional sd only at `σ = 0.63`; at the corrected σ the same literal is ±3.05 sd. Now computed from the process. | High |
| **B-8** | Several private copies of the state convention | Fixed | `recursive_main._sweep` held a duplicated 7-entry `x_ss` literal; `eds.py` defined its own `IS, IZ = 5, 6`; `accuracy.py` had `default_prob(states[:, 5])`, which after the state change reported `logistic(b_DF) = 54.9%` as the ergodic default probability. All single-sourced from `state_grid`. | High |
| **B-9** | TFP grid had 1/15 points frozen at `|F| = 1.7e−01` | Fixed | `solve_tfp` called `build_state_box` with bare defaults (`p_band = 0.25`); the period map has no solution at +25% wealth. Both experiments now share one `BOX_KW`. Pre-existing — present in every run before this work. | High |
| **O-1** | **8 bp is not attainable at the stochastic rest point** | Open (documented) | See §3. | High |
| **O-2** | Bond FOCs are the accuracy floor | Open | `bondFOC_D` −2.78, `bondFOC_F` −2.89 against −3.7…−5.2 elsewhere. `corr(|bondFOC_D|, |s−s*|) = +0.885`; fitted `Q_bD` overshoots to 0.981–0.986 and turns non-monotone where `p^d ≈ 0` and the true price is the risk-free ~0.946. | High |
| **O-3** | Risk-IRF output response is inside the solution's own error | Open | See §4. | Medium-High |

---

## 2. Structural rework

States 7 → 9: `[K_D, K_F, P_D, P_F, b_DD, b_DF, V_dep, s, Z_D]`.
Unknowns 7 → 11: added `Q_bD`, `b_DF`, `A_D`, `A_F`.
Residuals added: both banks' D-bond FOCs, `euler_F`, union deposit clearing.

- **D-sovereign market.** Previously `b_D_D = (1−shareF)·B'` at a fixed SS share, with the
  price read off the D bank's Euler at that imposed quantity — no demand schedule, no
  market clearing. Now both banks' FOCs are residuals, `b_DD = B' − b_DF` clears, and
  `Q_bD` is the price that does it. `psi_bD_F` (previously dead code) became load-bearing:
  at 0.05 a 0.5% price wedge supported a 54% position swing and the split was numerically
  indeterminate; set to 2.0, the split holds to 4 decimals.
- **Union deposit market.** `euler_F` was computed and discarded, and `A_D = dep_D/P_CES_D`
  force-fed each household its own bank's funding need — national clearing, not the union
  market the model is documented to run. Consequence: `C_D` was the bookkeeping residual of
  the bank balance sheet, with income contributing ~2% of its movement against ~42% from
  each gross leg. Both savings are now solved with clearing explicit.
- **Verification.** SS rest point exact across all 11 residuals (`≤1e−10`, `mu_D = 0.001001`,
  `C_D` and `A_D` at their SS values); `test_recursive_nesting` N2 probe 9.04 → 9.7e−15.

Independent confirmation that the union market was the right change: the **TFP** experiment's
consumption response, which had the wrong sign, corrected.

| | national clearing | union market |
|---|---|---|
| `Y_D` q0 | +1.147% | +1.494% |
| `C_D` q0 | **−0.084%** | **+0.608%** |
| `I_D` q0 | +1.80% | +2.87% |

---

## 3. The 8 bp target (O-1)

`calibrate_bank_targets` solves λ and ω_ent analytically from the leverage and spread
targets at the **deterministic** SS. The solved model does not rest there: ergodic
E[p^d] is 0.27% against 0.10% at `s*`, so the bank permanently holds the sovereign at
~0.92 rather than 0.946 and its divertable base is smaller at any λ.

Each instrument was measured with a full solve:

| instrument | measured | interpretation |
|---|---|---|
| `credit_spread_target` 8.04 → 0.04 bp | spread 32.8 → 24.6 bp, **λ = 0.199980 throughout** | disconnected: `λ = α/θ` and `α → Ω(1+rdep)` as `μ → 0`, so λ is pinned by *leverage* |
| `leverage_target` 5.0 → 5.5 | spread 32.8 → **39.8** bp (μ 0.0042 → 0.0056) | **wrong sign**: more assets per unit net worth binds harder, dominating the λ effect |
| `f` 0.02 / 0.04 / 0.08 | 17.6 / 25.5 / 62.5 bp; leverage 5.087 / 5.089 / 5.114 | clean and orthogonal to leverage, but convex (525 bp/unit near 0.02, 743 near 0.08), so `f → 0` floors at **12–14 bp** |

Reaching 8 bp requires `f < 0`. **Conclusion:** the deterministic 8 bp and the stochastic
8 bp are different objects; the gap is the priced risk premium. Leverage 5.0 and f = 0.04
retained; the model rests at ~25 bp.

---

## 4. Output response (O-3)

Impact `Y_D` = −0.028%, trough −0.028%, against Bocola Table 5's −1.05 / −1.44 / −1.53%.

The model's own labour FOC implies hours −0.145% on impact; solved hours move +0.003%.
The 0.148 pp gap corresponds to a `lab_D` residual of ~7e−4, which is that equation's
own 90th-percentile error.

Two hypotheses were tested and one survived:

- **"The feared event is expansionary"** — confirmed (B-4), fixed. But after the fix the
  IRF impact moved only +0.002% → −0.029%, i.e. flipping the feared event from
  expansionary to a −2.60% recession barely moved the IRF.
- **"`P_CES` falls enough to offset the working-capital wedge, and the household smooths
  through the cross-border position"** — **rejected**. `P_CES` offsets 0.0185 of the 0.0900
  wedge (≈20%, not the near-cancellation required), and `nfa_dep` is −0.0010 at impact.

Accuracy improved substantially across this work (`ALL` −2.96 → −3.83) while the impact
response did not trend with it (−0.029%, −0.041%, −0.028%). Candidate binding constraints
on magnitude, in order of suspicion and **not yet tested**: (i) the IC goes slack at q4,
after which there is no pass-through at any grid resolution; (ii) exposure 7.6% × leverage
5 caps the mark-to-market loss; (iii) the Neumeyer-Perri wedge is the only spread→output
channel. Grid refinement addresses none of these directly.

---

## 5. Calibration changes made

| parameter | from | to | basis |
|---|---|---|---|
| `debt_root_D/F` | — (`phi_lamb = 1.0`) | 0.95 → 0.93 | debt stationarity; root is now the calibrated object |
| `sigma_s` | 0.63 | 0.4455 | quadrature convention (B-6), not a recalibration |
| `psi_bD_F` | 0.05 | 2.0 | pins the bond split, previously indeterminate |
| `kappa_nfa` | 0.0 | 0.01 | Bocola's own foreign friction; keyed to the real position, not a `P_D−P_F` proxy |
| `S_COVER_SD` | (4.35 absolute) | 2.75 sd | see below |
| P bands | 0.12 / 0.20 | 0.04 / 0.04 | old width assumed `P_X` was both obligation and claim |
| `leverage_target`, `f`, `credit_spread_target` | 5.0, 0.04, 8 bp | **unchanged** | see §3 |

`S_COVER_SD = 2.75`, not Bocola's 2.16: the headline shock (p^d 0.10% → 1.98%) is itself
2.11 sd, so at 2.16 coverage it sits at 97% of the box half-width. The box cannot be
shifted instead, because `s*` must remain the box centre or it stops being a collocation
node and the exact SS rest point is lost. This is a fidelity-versus-usability trade-off
resolved in favour of containing the experiment.

---

## 6. Reproduction

```bash
cd code/global
python3 main.py                                    # RUN_TPI = False in this record
python3 tests/test_recursive_nesting.py            # SS rest point, all 11 residuals
python3 -m solver_recursive.calibrate_stochastic   # §3, the instrument sweep
```

`test_ss_identities` fails on `sov exposure/net worth = 0.36` against a 0.7–1.1 assert.
This is **pre-existing** (it fails identically at HEAD) and is the open "exposure test
threshold" item, not a consequence of this work.
