# A stochastic LTRO backstop: implementation plan

**Date:** 2026-08-31
**Status:** IMPLEMENTED 2026-08-31 (steps 1-3 of S11); calibration and experiments in
progress. Gates passing: N1 (SS rest point, 13 residuals ~1e-10, facility not drawn),
N3 (`phi_ltro = 0` nests the no-backstop model **bit-for-bit**; `ltro = 0` to 1 ULP --
see below), N4 (the facility touches ONLY the constraint), N2 (grid-wide pi = 0 solve,
1.3e-10), plus the five fast suites. Delivered dimensions match the plan exactly:
**10 states, 13 unknowns, 19 stored rules, 4 regimes**.

*One refinement the tests forced.* The two off switches hold to DIFFERENT tolerances and
the difference is not sloppiness. At `phi = 0` the facility regimes carry identically zero
weight, `_regime_weights` aliases them to regime 0, `np.dot` contributes an exact `0.0`
and the sum is bit-identical to the two-regime model. At `phi > 0` with `ltro = 0` the
economy is the same but the arithmetic is not: the same expectation is accumulated as
`(1-phi)a + phi*a` instead of `a`, which differs by one ULP. Demanding bit-identity there
would be demanding that floating-point addition be associative. `test_recursive_nesting`
now encodes both tolerances explicitly.
**Question.** With per-period probability φ the central bank stands ready to conduct a
Bocola-style LTRO — collateralised lending that relaxes banks' incentive constraint.
Agents internalise this. Because they fear the crisis states less, the economy is
stabilised **even along the realised path on which the facility never fires**. This is
the OMT fact — announced September 2012, never used, ~200 bp of compression — expressed
through the instrument the ECB actually deployed.

---

## 1. Why this instrument and not bond purchases

The completed yield-peg experiment established the constraint that motivates this plan.
The D bank's own first-order condition is

```
E[Ω·payD]  =  Q_bD · ( E[Ω]·R + λ_bD·μ )
```

A **bond purchase** raises `Q_bD` only by shrinking the divertable base and pushing `μ`
down, and `μ` is floored at zero. Holding the continuation fixed, no quantity can lift
the price above `E[Ω·payD]/(E[Ω]·R)` — the same claim with the **liquidity premium
removed and nothing else**. Measured (`peg_feasibility_report`): that premium is
**0.2–0.7%** of the price across the grid and **0.63%** at the crisis corner, against a
22.2% gap to a peg set at the φ=0 rest-point price. The converged peg walk reached only
Q\* ≈ 0.731 (824 bp/yr) against a 277 bp target, at 27–33% of the outstanding stock, and
credibility barely moved it (0.7311 at φ=1 vs 0.7318 at φ=0.5).

Bocola ships the alternative in his own replication package
(`Model/Solution Files/residual_model_ltro_firstperiod.m`, `ltro_policies.m`, Figure 8).
It is a different margin entirely:

```
baseline    μ_ratio =  N'      / ( λ · ( Q·K' + q·B'     ) )
LTRO        μ_ratio = (N' + m) / ( λ · ( Q·K' + q·B' − m ) )
```

Central-bank credit of size `m` does **two** things where a purchase does one: the assets
it funds leave the divertable base **and** the funding counts as equity in the constraint.
Measured on our calibration (n = 1.9338, leverage 5.0, λ = 0.2363, assets 9.6692):

| operation | n/(λA) | relief |
|---|---|---|
| none | 0.84625 | — |
| CB **buys** bonds worth m = 0.4 | 0.88276 | +4.3% |
| CB **lends** m = 0.4 against collateral | 1.06536 | **+25.9%** |

The 6.0× is not a coincidence. To first order

```
Δ(LTRO) / Δ(purchase)  =  1 + A/n  =  1 + leverage  =  6
```

at Bocola's leverage of 5. **The `+m` in the numerator is a margin no amount of
bond-buying can reach**, and it is unavailable to any policy that operates through the
bond market alone.

---

## 2. The mechanism is a ONE-EQUATION change

This is the plan's central claim and it deserves the argument in full.

Let the facility be drawn at the deposit rate, `r_ltro = rdep`. The bank's balance sheet
becomes `assets = deposits + m + n` instead of `assets = deposits + n`, so

```
deposits  =  assets − n − m
P'        =  (1 + rdep)·deposits + (1 + r_ltro)·m − (1 + r_wc)·L_wc
          =  (1 + rdep)·(assets − n)             − (1 + r_wc)·L_wc      [r_ltro = rdep]
```

which is **exactly** the existing `Pp_D`. On the household side the claim is unchanged in
total: a euro of bank deposit is replaced by a euro of central-bank claim at the same
rate, so `save_union`, `dep_union`, `nfa_dep_D` and `Vp_dep` are all untouched. The CB
lends at the rate it pays and bears no credit risk (net worth is floored), so its carry
is identically zero and **no remittance identity is needed** — the Walras leak that
dogged the purchase design cannot arise here.

> **The LTRO changes the COMPOSITION of the bank's funding — divertable deposits for
> non-divertable central-bank credit — at an unchanged rate. Every budget identity in the
> model is therefore unchanged, and the entire effect passes through the diversion
> constraint.**

Concretely, in `point_map.py` only these lines move:

```python
m_D      = cal["ltro_D"] if m_reg else 0.0     # facility drawn this period
n_IC_D   = n_D + m_D                            # counts as equity in the constraint
lev_IC_D = max(lev_D - lbDD * m_D, 1e-6)        # and leaves the divertable base
mu_D     = _smin(max(1.0 - E_Om_D*(1.0+rdep_D)*n_IC_D/lev_IC_D, 0.0), _MU_CAP, _GUARD_EPS)
slack_D  = alpha_D_cur * n_IC_D - lev_IC_D
```

`dep_D = assets_D − n_D` keeps **actual** net worth, per the standing convention
(CLAUDE.md: *"Portfolio shares and branch initial conditions divide by ACTUAL net worth,
not n_IC"*). The distinction already exists in the codebase (`bank.py`'s `n_ss_IC`,
`n_IC_D`), so this extends an established object rather than inventing one. `alpha`,
`r_wc = rdep + λ_K·μ/E[Ω]` and every downstream residual pick up the relieved `μ`
automatically.

Under the single-λ doctrine all three λ are equal (0.2363), so `lev_D = λ·assets_D` and
`lev_D − λ·m` is Bocola's `λ(A − m)` exactly. Writing it as `lbDD·m` keeps it correct if
the λ ever diverge, with the collateral read as D-sovereigns.

---

## 3. The economics to be measured — three channels, one of which runs backwards

The plan's value is that the sign is **not obvious ex ante**. Three channels operate.

**(a) Liquidity premium, direct and stabilising.** In the relieved regime `μ` falls, so
`λ_bD·μ` falls out of the denominator and `Q_bD` rises. In the same regime
`r_wc = rdep + λ_K·μ/E[Ω]` falls, so — through the Neumeyer-Perri wedge, the model's only
channel from the financial block into output — hours and output rise directly.

**(b) Risk premium, indirect and stabilising — THIS IS THE USER'S MECHANISM.** The
premium is the covariance leg `E[Ω·payD] / (E[Ω]·E[payD])`. `Ω` is high exactly in the
states where `payD` is low. The facility lowers `μ'`, hence `α'`, hence `Ω'`, **most in
the states where the constraint is tightest** — which are the same states where `payD` is
lowest. The covariance shrinks, the risk premium falls, and today's price rises **in
every state, including those where the CB is absent**. Nothing is imposed; it falls out
of the quadrature. This is why the never-fired path is stabilised, and it is the reason
the default regime must carry the facility (§4).

**(c) Franchise value, indirect and DESTABILISING.** This is the channel that could
reverse the result and it must be confronted rather than hoped away.
`α = E[Ω]R/(1−μ)` and `Ω = β·[f + (1−f)α']`. Lowering future `μ'` lowers `α'`, lowers
`Ω'`, and therefore lowers `E[Ω]` — which **raises** today's

```
μ = max( 1 − E[Ω]·R·n/lev , 0 )
```

The bank's charter value *is* its collateral in a Gertler-Karadi economy: make the future
safer and the bank has less to lose, so the constraint binds harder today. The channel is
not second-order here — `μ` is a small difference of numbers near one, so a 1% fall in
`E[Ω]` moves `μ` by roughly `0.99 × 0.01 = 0.0099`, which against `μ_ss = 0.01231` is a
**~80% increase**. At the algebraic limit (`βR = 1`, `μ → 0`) `α` collapses from 1.1817
to exactly 1.

Whether (a)+(b) beat (c) is a general-equilibrium question that only the solve answers.
**That is the research content of the experiment, and the model is already instrumented
to decompose it**: `output_decomposition.decompose_bond_price` splits `log Q_bD`
additively into continuation/duration, expected loss, risk premium and liquidity premium.
Running it across φ measures which channel does the work — and if (c) dominates, that is a
publishable negative result about standing liquidity backstops, not a failed experiment.

---

## 4. Regimes: the facility must be available in the default state

The compound regime index `(d′, m′)` and its table already exist
(`decision_rules.regime_table`). The peg used three sets because a central bank does not
peg a defaulted bond. **The LTRO is different: it is liquidity support to BANKS, and the
default state is precisely when banks need it.**

More importantly, channel (b) *requires* it. The default branch is where `payD` is lowest
(a 55% haircut) and `Ω` highest; it carries little probability mass but a very large
payoff deviation, so it is disproportionately important for the covariance. Withdrawing
the facility there removes the largest single term in the risk-premium channel.

**Baseline: four regimes**, `(d,m) ∈ {(0,0),(0,1),(1,0),(1,1)}`, m orthogonal to d.

`_regime_weights` needs generalising so the m-probabilities are conditional on d — one
function, correct for every table:

```python
def _regime_weights(wq, pd, phi, reg):
    out = []
    for d_n, m_n in reg:
        p_d = pd if d_n else (1.0 - pd)
        rows = [m for dd, m in reg if dd == d_n]
        p_m  = 1.0 if len(rows) == 1 else (phi if m_n else 1.0 - phi)
        out.append(wq * p_d * p_m)
    return out
```

This reproduces the 2-regime and 3-regime tables exactly (each `d` with a single row takes
all of that `d`'s mass), so the existing nesting gates are unaffected.

*Documented variant, 25% cheaper:* three regimes with no facility in default. This is
**historically accurate** — the ECB suspended the collateral eligibility of Greek
government bonds in February 2012 during the PSI and again in February 2015 — and it
produces the collateral cliff, an amplification mechanism worth reporting in its own
right. It is a different experiment, not a cheaper version of this one.

---

## 5. Calibrating the facility — and the trap in Bocola's own size

Size matters more than it looks, for a reason that connects to the 2026-08-29 kink
finding. Backing `E[Ω]R = 1.16714` out of `μ_ss = 0.01231` and solving for the facility
that drives `μ` to zero:

| facility m | % of quarterly GDP | μ at the SS | μ in the crisis state |
|---|---|---|---|
| 0.000 | 0 | 0.01231 | 0.02339 |
| 0.005 | 0.5 | 0.00924 | 0.02032 |
| 0.010 | 1.0 | 0.00617 | 0.01725 |
| **0.020** | **2.0** | **0.00003** | 0.01110 |
| **0.034** | **3.4** | 0 | **0.00248** |
| 0.100 | 10 | 0 | 0 |
| **0.400** | **40 — Bocola's own** | **0** | **0** |

A facility of **2.0% of quarterly GDP unbinds the constraint at the steady state**, and
**3.4% unbinds it in the crisis state**. Bocola's 40% is roughly **twelve times** the size
that fully neutralises the crisis.

**This is a trap, not a detail.** At his size `μ = 0` with enormous margin in every
relieved state, so the entire `m=1` coefficient set sits **on the KKT kink** — exactly the
region where `μ = max(·,0)` is C0, a Chebyshev interpolant returns `μ > 0` where the truth
is 0, and the fitted-versus-exact read disagrees by more than the response being measured.
That is the pathology the 100 bp spread recalibration was adopted to escape.

**Calibration rule: size `m` so the facility RELIEVES the constraint without unbinding
it.** A target such as *"halve the crisis-state `μ`"* gives `m ≈ 0.010–0.015`, i.e.
**1.0–1.5% of quarterly GDP**, which keeps `μ > 0` in both regimes and keeps the solution
off the kink. Ship that as the baseline; run Bocola's 40% as a documented upper-bound
variant with the identification caveat attached, not as the headline.

Two further sizing notes. The facility is **per country** (`ltro_D`, `ltro_F`), with
`ltro_F = 0` for the targeted experiment and `ltro_F = ltro_D` for the union-wide one —
the actual LTROs were euro-area-wide. And a natural extension makes the envelope
**collateral-linked**, `m = ltro_frac · Q_bD · B_D'`, which is smooth (no `min` needed if
the envelope always binds), procyclical, and reproduces the real doom loop: the facility
shrinks exactly as the collateral it is secured against loses value.

---

## 6. Accounting and cost

The LTRO needs **no new state and no new unknown**. `b_cb` and `x_cb`, added for the peg,
are not required and should be retired for this experiment — which takes the model back to
the pre-TPI dimensions and makes this design *cheaper per regime* than the one it replaces.

| | pre-TPI | yield peg (built) | **LTRO backstop** |
|---|---|---|---|
| states | 10 | 11 (+`b_cb`) | **10** |
| unknowns / point / regime | 13 | 14 (+`x_cb`) | **13** |
| stored rules | 19 | 20 | **19** |
| regimes | 2 | 3 | **4** |
| coarse μ=1 points | 21 | 23 | **21** |
| s-refined (m=5) points | 95 | 105 | **95** |
| total collocation unknowns | 3,610 | 6,300 | **7,220** |
| coarse Jacobian | ~2 min | ~5.6 min | **~5.6 min** |
| s-refined Jacobian | 29 min | ~83 min | **~114 min** |

The refined solve is ~7–8 h for four Newton steps, 10–14 h with the ladder. The coarse
grid is 5.6 min a Jacobian, so **every sign test, nesting gate and channel decomposition
below runs on the coarse grid** and only the final numbers need the refined one. The
three-regime variant is 5,415 unknowns and ~64 min a Jacobian if the budget binds.

`s_refine = 9` is 12,996 unknowns at ~370 min a Jacobian — record it as out of reach and
do not plan around it.

---

## 7. Changes by file

**`solver_recursive/decision_rules.py`** — add `4: ((0,0),(0,1),(1,0),(1,1))` to
`_REG_TABLE`. Nothing else; `RuleSet` is already regime-count-generic.

**`solver_recursive/point_map.py`** — the substantive work, and it is small.
1. `_regime_weights` → the conditional form in §4.
2. The four lines of §2 in the multiplier block (`m_D`, `n_IC_D`, `lev_IC_D`, `mu_D`), and
   the same for F. `slack_D/F` follow `n_IC`/`lev_IC`.
3. `out` gains `m_ltro_D/F`, `n_IC_D/F` and `lev_IC_D/F` for the diagnostics.
4. `no_cb` (already plumbed) forces `phi = 0`, so it doubles as the LTRO's off switch.
5. Retire `x_cb` from `SOLVE`, `b_cb` from the state, and the CB purchase/remittance block
   — **or** keep them behind `phi_tpi`/`Q_peg_D` if both policies are to be compared. Note
   that carrying the peg's state and unknown while running the LTRO costs ~20% for nothing.

**`solver_recursive/state_grid.py`** — revert `STATE_NAMES` to 10 if `b_cb` is retired.

**`config/calibration.py`** — `phi_ltro` (default 0.0, the nesting value), `ltro_D`,
`ltro_F`, and `ltro_frac` for the collateral-linked variant. Document the §5 sizing table
at the parameter.

**`solver_recursive/recursive_experiment.py`** — the solve ladder needs **no peg walk and
no adaptive homotopy**: the LTRO regime is not a complementarity, has no free quantity and
no fixed point to bootstrap, so it is seeded from `d=0` and solved directly. Expect the
ladder to be *simpler* than the current one: coarse d=0 → haircut homotopy → joint. If a
homotopy is needed at all it is on `m` from 0, which is trivially slack at 0.

**`solver_recursive/tpi_recursive_experiment.py`** — the driver already loops over an
activation scalar and reads IRFs off `dynamic_irf`'s rest-point baseline; it needs the
parameter renamed and the CB-footprint panels changed from purchases to facility draw.

---

## 8. Experiments

**E1 — the never-fired path (the headline).** Solve at φ ∈ {0, 0.25, 0.5, 0.75, 1}. Read
the IRF along the regime-(0,0) path throughout: the realisation on which the facility is
announced and **never drawn**. Report `Y_D`, `C_D`, `I_D`, the sovereign spread and
`Q_bD` against φ = 0. The whole difference is the announcement effect. `read_at` and
`dynamic_irf` already read regime 0 by default, so this needs no new machinery.

**E2 — channel decomposition (the diagnostic that decides §3).** Run
`decompose_bond_price` at each φ and report the four legs. The prediction is that the
**risk-premium leg** carries the compression on the never-fired path while the
**liquidity leg** is confined to the drawn regime. If instead the liquidity leg dominates,
the mechanism is not the one claimed and the result is about (a), not (b).

**E3 — the franchise-value counter-test.** Track `E[Ω]`, `α_D` and `μ_D` at the
stochastic rest point as functions of φ. If `μ_D` **rises** with φ on the never-fired
path, channel (c) is offsetting and the size of the offset is itself the finding. Report
it whichever way it goes.

**E4 — realisation band.** Monte Carlo over `m`-draws with **common random numbers**
across the shocked and unshocked paths, plus the two deterministic bounds (never drawn,
always drawn). The gap between the mean and the never-drawn path is the balance-sheet
channel; the gap between never-drawn and φ = 0 is the pure announcement.

**E5 — take-up and cost.** Expected drawdown as % of GDP, against Bocola's 40% and the
actual 3-year LTROs (~EUR 1tn, ~10% of euro-area GDP). Since the carry is zero by
construction, the fiscal cost is zero and the policy's cost is entirely the moral-hazard
margin in E6.

**E6 — the endogenous cost.** The model delivers this for free and it should be reported:
banks that expect relief **lever up ex ante**, because the constraint is looser in
expectation. Measure leverage and `slack` at the rest point against φ. This is the
charter-value/risk-taking cost of a standing backstop, and it is a genuine welfare offset
rather than an artefact.

---

## 9. Verification

Cheap gates first; all on the coarse grid unless noted.

1. **φ = 0 nests exactly.** The existing N3 gate, extended to four regimes: the regime-0
   residuals at φ = 0 must equal the two-regime model's to **0.0**, not to tolerance.
2. **`ltro_D = 0` nests exactly**, independently of φ. Two separate off switches, both
   exact, is what makes any measured effect attributable.
3. **Budget closure — the strong test.** §2 claims no flow changes. Therefore `goods_F`,
   the Walras-redundant diagnostic, must be **unchanged to machine precision** between
   `ltro = 0` and `ltro > 0` at the same policies. If any budget identity moved by
   accident, this catches it. Nothing else in the suite would.
4. **Constraint algebra.** Point-wise on the grid, in the drawn regime: `μ` must fall,
   `slack` must rise, `r_wc` must fall, monotonically in `m`.
5. **Collocation floor.** `sum(F²) ≤ m·(1e-9)²` at every ladder stage over 19 equations ×
   points × 4 regimes.
6. **The kink bracket (§5).** Report `read_at` (fitted) against `read_exact` (period map
   cleared exactly) at every reported point, **for each φ**. The bracket is expected to
   *widen* with φ as more points approach `μ = 0`. If it widens past the response, reduce
   `m` per §5 and say so. This is the accuracy limit of the whole exercise and it must be
   reported next to the headline number, exactly as `impact_table` already does.
7. **Signs.** Higher φ ⇒ `Q_bD` ↑, sovereign spread ↓, credit spread ↓, `Y_D[0]` less
   negative, monotone in φ. Non-monotonicity is a bug **unless** E3 shows channel (c)
   dominating, in which case it is the result — the two must be distinguished before
   either is reported.
8. Existing suite unchanged.

---

## 10. Risks, and what would falsify the claim

| risk | how it shows up | response |
|---|---|---|
| **Franchise value dominates (§3c)** | `μ` rises with φ on the never-fired path (E3); spreads widen | Report it. This is a real result about standing backstops, not a failure |
| **Facility oversized → μ = 0 everywhere (§5)** | fitted/exact bracket blows up at high φ | Resize `m` to halve rather than eliminate crisis `μ` |
| **Four regimes unaffordable** | refined Jacobian ~114 min | Fall back to three regimes, and report that the default-state leg of channel (b) is then missing |
| **Take-up implausible** | E5 far from ~10% of GDP | Recalibrate; the §5 table shows the constraint binds at ~1–3% of quarterly GDP, so realistic sizes are comfortably above what is needed |
| **Time consistency** | not modelled — the CB commits | State the limitation. A CB that reneges is a separate regime and a separate paper |

**The claim is falsified if** E1 shows no material stabilisation on the never-fired path
at φ > 0, **or** if E2 attributes the compression to the liquidity leg rather than the
risk-premium leg. Either outcome is reportable and neither is a solver failure.

---

## 11. Sequencing

1. `_REG_TABLE[4]` + the conditional `_regime_weights`. **No economics.** Gate: the
   2- and 3-regime tables still reproduce their residuals bit-for-bit.
2. The four lines in the multiplier block, behind `ltro_D = 0`. Gate: tests 1–3.
3. Retire `b_cb`/`x_cb` (or gate them behind the peg parameters). Gate: full fast suite.
4. Size `m` per §5 on the coarse grid — solve at three sizes, read crisis `μ`, pick the
   one that halves it.
5. Coarse-grid E1–E3 across φ. **This is where the paper's answer appears**; if channel
   (c) dominates, stop and report rather than proceeding to the refined solve.
6. E4–E6, coarse.
7. Refined ladder at `s_refine = 5`. Final numbers, with the test-6 bracket beside each.

---

## 12. First results, and two corrections the run forced (2026-08-31)

**The mechanism works, and E3 is the evidence.** At an envelope of 2.0% of quarterly GDP,
read at the model's own stochastic rest point with the facility **never drawn**:

| φ | μ at the rest point | credit spread | E[Ω] | α | drawn |
|---|---|---|---|---|---|
| 0% | 0.00983 | 80.0 bp/yr | 1.161834 | 1.17693 | — |
| 50% | 0.00627 | 51.1 bp/yr | 1.160439 | 1.17106 | **0.00%** |
| 100% | 0.00000 | 0.0 bp/yr | 1.160399 | 1.16330 | **0.00%** |

`d(μ)/d(φ) = −0.0098`: **the relief channels dominate the charter-value channel.** The
latter is real and measurable — E[Ω] and α both fall with φ, exactly as §3(c) predicted —
but it does not reverse the sign. Agents price a facility that is never used and the
constraint is looser for it, which is the claim the experiment was built to test.

The four-regime collocation gate shows the two channels separately: α falls **1.18%** in
the no-default facility regime and **2.45%** in the default one. The *level* fall is the
charter-value channel; the *asymmetry* — Ω compressed most where `payD` is worst — is the
risk-premium channel of §3(b).

### Correction 1: the envelope was sized against the wrong multiplier

§5 sized `m` against μ at the **grid centre** (0.01979). The object that must stay off the
KKT kink is the **stochastic rest point**, where μ is 0.00983 — half. At 2.0% the facility
therefore put the ergodic point *on* μ = 0 at full credibility, and the identification
went with it: the fitted-versus-exact output bracket at φ = 1 was
**[−0.1524%, +0.0148%]** — wider than the response and straddling zero — and that solve
stopped at max|F| = 1.5e-04 rather than converging.

Sizing against the rest point instead, with the measured slope
`d(μ_rest)/d(φ·m) = −0.4915`:

```
mu_rest(phi = 1)  =  0.00983 - 0.4915 * m      ->   m = 1.0%  gives  0.0049
```

**1.0% of quarterly GDP ships.** φ and m enter only as a product, so full credibility is
the worst case and sizing there covers every φ. A larger facility is not more policy; it
is less identification.

### Correction 2: the IRF differences away the thing being measured

`dynamic_irf` differences each φ's shocked path against **its own** unshocked path. That
is right for a single experiment and wrong for a comparison across policy regimes: the
backstop's main effect is to **move the ergodic point**, and differencing against that
moved point removes it. E1 as originally specified therefore reports the shock response
*conditional on the regime*, not the stabilisation.

Added **E3b**: every rest point read against the φ = 0 rest point — output, consumption,
investment, the bond price and both spreads, with nothing ever drawn. That is the level
shift, and it is what "the economy is stabilised even when the facility never fires"
actually means. E1 is kept and relabelled as the conditional response.

---

## 13. The result (coarse grid, envelope 1.0% of quarterly GDP)

### E3b — where the economy rests, with the facility NEVER DRAWN

| φ | Y_D | C_D | I_D | Q_bD | sovereign spread | credit spread | drawn |
|---|---|---|---|---|---|---|---|
| 0% | — | — | — | — | — | — | — |
| 50% | **+0.173%** | +0.024% | **+0.653%** | +0.216% | −3.3 bp | **−29.5 bp/yr** | **0.00%** |
| 100% | +0.315% | +0.084% | +1.131% | +0.357% | −3.2 bp | −80.0 bp/yr | 0.00% |

**The claim the experiment was built to test is confirmed.** Agents price a facility that
never fires, and the economy rests at higher output and investment with a materially lower
lending spread. `d(μ)/d(φ) = −0.0098`: the relief channels dominate the charter-value
channel of §3(c), which is present and measured (α falls 1.18% in the no-default facility
regime, 2.45% in the default one) but does not reverse the sign.

### The nuance that decides how this is written up

**The backstop compresses the CREDIT spread by 80 bp and the SOVEREIGN spread by 3.**
That is the mechanism being internally consistent, not a defect. The facility acts on the
*bank's* constraint; the sovereign bond price can only move through the liquidity premium,
which `liquidity_ceiling_report` independently bounds at 0.2–0.7% of the price — and
Q_bD duly moves +0.36%. So the finding is:

> A Bocola-style facility stabilises the real economy **through bank balance sheets**,
> not by making sovereign debt safer.

That is arguably the better reading of what the 3-year LTROs did, but it is a *different*
claim from "it stabilises spreads" and the paper must say which. It also explains why the
yield-peg design failed: it was trying to move the object this instrument cannot move.

### E1 — the conditional shock response

| φ | Y_D impact | C_D | I_D |
|---|---|---|---|
| 0% | −0.1021% | −0.0814% | −0.629% |
| 50% | −0.1098% | −0.0697% | −0.749% |
| 100% | −0.1135% | −0.0504% | −0.875% |

Conditional on the regime, the same shock does *slightly more* damage under the backstop —
monotone and small. Coherent rather than anomalous: from a less-constrained rest point the
constraint has more room to tighten. The total effect is E3b + E1, strongly positive.

### What is and is not identified

- **φ = 0 and φ = 0.5: converged and usable.** Fitted-vs-exact brackets 0.021 and 0.026 pp
  against responses of ~0.11 pp — 20–24%, the model's documented identification state.
- **φ = 1: NOT converged** (stopped at max|F| = 4.0e-05), and its rest point has μ = 0 on
  *both* the fitted and the exact read, so output there is unidentified (bracket
  [−0.1135%, +0.0185%], straddling zero). Report it as indicative or not at all.

The resize from 2.0% to 1.0% did not rescue φ = 1, because **the relief is strongly convex
in φ**: at 1.0% the reduction in μ_rest at φ = 1 is more than 2.7× the reduction at
φ = 0.5, not 2×. Full credibility genuinely unbinds the constraint at the ergodic point —
economic content, not an artefact — which puts the solver on the KKT kink. The honest
presentation is the identified range φ ∈ [0, 0.5] plus the statement that full credibility
unbinds the constraint, rather than a converged-looking number that is not one.

---

## 14. The activation curve (coarse grid, envelope 1.0% of quarterly GDP)

Assembled from two runs at identical calibration; φ = 1 appears in both and agrees
(μ_rest = 0 in each), which cross-checks them.

| φ | μ at the rest point | credit spread | Y_D vs φ=0 | I_D vs φ=0 | converged | IRF bracket |
|---|---|---|---|---|---|---|
| 0 | 0.00983 | 80.0 bp/yr | — | — | yes (1.6e-10) | 0.021 pp (20%) |
| 0.50 | 0.00619 | 50.5 bp/yr | **+0.173%** | **+0.653%** | yes (1.4e-09) | 0.026 pp (24%) |
| 0.75 | 0.00288 | 23.4 bp/yr | ≈ +0.250% | ≈ +0.864% | yes | 0.067 pp (61%) |
| 1.00 | 0.00000 | 0.0 bp/yr | +0.315% | +1.131% | **NO** (7.5e-05) | straddles zero |

The credit spread falls almost linearly in φ and reaches zero at full credibility: a
facility of 1% of quarterly GDP, believed with certainty, **unbinds the intermediary's
constraint at the ergodic point** — with nothing ever drawn.

**Two things this table settles.**

*The three-rung facility ladder was worth adding.* φ = 0.75 stopped at 4.0e-05 with a
single rung and converges with three. φ = 1 still does not (7.5e-05), and that is now
attributable to the economics rather than the basin: at φ = 1 the rest point has μ = 0 on
BOTH the fitted and the exact read, so the KKT max{·,0} is active at the ergodic point and
the Newton is resolving a kink, not a badly-seeded smooth problem.

*The LEVEL results are better identified than the SHOCK results, and by construction.*
E3/E3b read the rest point, where μ > 0 for every φ ≤ 0.75. E1 reads the shocked state,
where the fitted interpolant under-reads μ at every φ (0.00844 against an exact 0.01962 at
φ = 0; 0 against 0.00492 at φ = 0.75) — the documented Gibbs behaviour of a C0 multiplier,
worsening as the facility pushes μ toward the kink. Report E3b as the headline and E1 with
its bracket attached, never the reverse.

**Identified range: φ ∈ [0, 0.75].** φ = 1 is reportable as the limiting statement "full
credibility unbinds the constraint", not as a number.
