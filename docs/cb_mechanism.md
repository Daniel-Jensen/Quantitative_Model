# The central-bank mechanism — canonical statement

**Status: canonical.** This is the reference description of what the TPI central bank
is, how it transmits, and what may and may not be claimed about it. Established by the
CB-block audit of 2026-08-19 (`diagnostics/cb_audit/`) against the working tree of
`gk-structural-foc` after the GK structural refactor. Every number below is measured,
not asserted; the probes that produced them are named at each point.

Supersedes ad-hoc descriptions of the CB in `docs/STATE.md` and in figure captions.
When this document and a caption disagree, this document is right and the caption is
stale.

---

## 1. What the object is — and what it is not

The model's "ECB" is a **capital-key conduit with full per-period pass-through**. It
is not a central bank in the institutional sense, and the difference matters for what
the paper can claim.

**It has:**

* a bond book, `cb_buy_D` (CB holdings of D sovereign paper), a D aggregate;
* a net cash flow on that book, `cb_flow_D` (CB net cash flow, D goods), computed in
  the same coupon/survival form as the private payoff;
* a capital key, `kappa_cb_F` (F share of CB profit and loss) `= 0.929`, splitting
  that flow between the two treasuries;
* remittances to **both** treasuries — `rem_cb_D` and `rem_cb_F` (D and F treasury
  receipts from the CB) — that enter their respective budget constraints.

**It does not have:**

* **capital or retained earnings.** The entire net cash flow is remitted every period.
  There is no loss-absorption buffer, so a loss reaches the two treasuries in the
  period it occurs. A real Eurosystem national central bank absorbs losses against
  capital and provisions first.
* **a reserve liability, or a policy rate.** Purchases are funded by a same-period
  capital call on the two treasuries, not by creating remunerated reserves. This is
  internally consistent — the model has no policy rate at all, the union-inflation
  normalisation being the `phi_pi` (Taylor coefficient on union inflation) → ∞ limit —
  but the seigniorage and reserve-remuneration legs of an actual APP or TPI are absent
  by construction.

Both absences are deliberate modelling choices and both must be stated as such in the
paper. They are not defects and they are not to be "fixed" without a design decision.

## 2. The four equations

`code/tpi.py` supplies four blocks through `tpi_overrides()`, swapped into the shared
`full_model.build_block_list()`. They are the whole central bank; there is no fifth
place where `cb_buy_D` appears.

| block | replaces | what it does |
|---|---|---|
| `domestic_bond_clearing_tpi` | `domestic_bond_clearing` | `b_D_D = b_gov_D − size_F·b_D_F − cb_buy_D` |
| `budget_residual_D_tpi` | `budget_residual_D` | computes `cb_flow_D`; remits `(1−kappa_cb_F)` of it to D |
| `budget_residual_F_tpi` | `budget_residual_F` | remits `kappa_cb_F` of it to F, converted `/p/size_F` |
| `external_account_D_tpi` | `external_account_D` | books the F-owned share `kappa_cb_F·cb_buy_D` as an external claim on D |

**Units.** `cb_buy_D` is already a D aggregate — the ECB buys a quantity of D debt —
so unlike `b_D_F` (F-bank holdings of D paper, per F capita) it takes no `size_F`
(F size relative to D) weight. The F remittance takes **two** conversions, `/p` for
the good and `/size_F` for the population; omitting the second leaked up to 2e−02 of F
GDP through `goods_mkt_F` (F goods-market residual) at γ = 10 while γ = 0 stayed clean
— the signature of a conduit-only units error.

**Pricing.** All four mark exclusively at the endogenous `q_b_D` (D sovereign bond
price). `cb_flow_D` matches `bond_return_D`'s state-contingent payoff term for term,
including the `zeta_writeoff_D` (continuation-value write-down switch) continuation
leg and the `writeoff_enabled_D` (realisation gate) multiplier. Audited line by line
(`run_log.md` Probe D); no stale price, no hardcoded price decomposition, no surviving
reference to the deleted pre-refactor FOC.

The GK structural refactor did **not** touch these four blocks. They were written
against the *payoff*, in coupon/survival form, never against the old FOC's price
decomposition, which is why they survived unchanged (`run_log.md` Probe C).

## 3. Steady-state neutrality is exact

`cb_buy_ss = 0` — TPI is dormant at rest. Evaluated at the solved steady state, every
TPI block output is **bit-identical** to its non-TPI counterpart:

| output | TPI | non-TPI | difference |
|---|---|---|---|
| `b_gov_res_D` (D budget residual) | −4.163336e−17 | −4.163336e−17 | **0.000e+00** |
| `nfa_D` (D net foreign assets) | −1.338351e−01 | −1.338351e−01 | **0.000e+00** |
| `ca_res_D` (D current-account residual) | +1.665335e−16 | +1.665335e−16 | **0.000e+00** |
| `b_D_D` (D-bank holdings of own paper) | +9.922837e−01 | +9.922837e−01 | **0.000e+00** |
| `b_F_F` (F-bank holdings of own paper) | +4.953501e−01 | +4.953501e−01 | **0.000e+00** |

and `cb_flow_D`, `rem_cb_D`, `rem_cb_F` are identically zero. Confirmed dynamically by
`run_tpi`'s own gate: `G_tpi[cb=0]` versus the baseline Jacobian, `max|err| = 0.00e+00`.

**Consequence:** no steady-state result and no SS-invariance argument anywhere in the
project can be affected by the CB block. The refactor did not move the SS through it.

## 4. How the spread is generated, and therefore what the CB can do to it

At the preferred baseline `psi_lambda_B_D/F` (collateral-friction amplification dial)
`= 0`, so `Delta_bD_eff_D` (effective pledgeability of D paper to D banks) `≡ 0.20`.
Compose `gk_bond_foc_D` (`nu_bD_D/nu_K_D = Delta_bD_eff_D`) with `intermediation_P1_D`
and the intermediary stochastic discount factor `SDF_banker_D` and franchise value
`Omega_p1_D` **cancel in the ratio**, leaving

```
rb_exp_D(+1) − rdep_D = 0.20 · ( rk_D(+1) − rdep_D )
```

exactly. Verified at the steady state: `rb_exp_D` (expected D bond return)
`= 0.002000`, `rk_D` (D return on capital) `= 0.010000`, `rdep_D` (D ex-ante real
deposit rate) `= 0`, and `0.20 × 0.01 = 0.002` to machine precision; the four-leg FOC
table verifies to ≤ 2.1e−13 on every solved SS.

**This is the single most important structural fact about the mechanism.** The
sovereign risk premium is pinned as a fixed fraction of the capital premium. The
central bank has no direct lever on the spread — there is no term of the form
`spread += parameter × def_rate` anywhere on its path, and never was one to revive.
**The only way the CB can compress the spread is by lowering `rk_D(+1)`, that is, by
crowding capital back in.**

It follows that in this model **TPI's spread effect and its investment effect are the
same effect**. They cannot be reported as two channels, decomposed against each other,
or traded off. Any sentence of the form "TPI compresses the spread, and separately
supports investment" is wrong.

## 5. The transmission chain

`cb_buy_D` reaches the price through exactly two doors — `intermediation_IC_D`'s
`phi_bD_D` (D-bank sovereign concentration) and `k_balance_sheet_D` — and through
nothing else. Measured off the TPI Jacobian at impact (`probe_stability.py`):

```
 d b_D_D    / d cb_buy [0,0] = -0.126779     modest quantity relief on D banks
 d theta_D  / d cb_buy [0,0] = -1.177721     required leverage falls, IC slackens
 d n_inter_D/ d cb_buy [0,0] = +0.525814     D-bank net worth recovers (MTM on q_b_D)
 d K_D      / d cb_buy [0,0] = +0.017015     capital crowded back in
 d q_b_D    / d cb_buy [0,0] = +0.037316     bond price up -> spread compresses
```

Read in the order the mechanism actually runs: the purchase supports `q_b_D`, the
price support recapitalises D banks through mark-to-market on the book they still
hold, the recapitalisation slackens the incentive constraint, capital is crowded back
in, `rk_D` falls, and §4's identity then requires `rb_exp_D` to fall — which *is* the
spread compression. The quantity relief is a minor contributor; see §6.

Aggregate effect, closed loop: `A_cb[0,0] = d(spread_rb)/d(cb_buy_D)[0,0] =
−4.397083e−03`, negative as required, and peak spread falls monotonically
205.87 → 193.25 → 176.75 → 154.36 bp annualised at γ = 0/2/5/10.

## 6. The decisive diagnostic — the 2×2 sovereign-holdings matrix

`diagnostics/cb_audit/probe_portfolio.py`, full output in `portfolio_matrix.md`.
Aggregate market value in D goods, `q_b × quantity`, with the per-F-capita legs scaled
by `size_F = 11.696651`. Both clearing identities close to ≤ 1.2e−15 at every point
reported, at the steady state and at t = 0, 4, 20 for every γ.

### 6.1 Steady state — a home-biased small system beside a large indifferent one

| holder | D paper | F paper | total | share of D issue |
|---|---|---|---|---|
| **D banks** | 0.967383 | 0.007282 | 0.974666 | **87.27%** |
| **F banks** | 0.141118 | 5.596796 | 5.737913 | **12.73%** |
| **CB** | 0.000000 | 0.000000 | 0.000000 | 0.00% |
| issued | 1.108501 | 5.604078 | 6.712579 | |

The D banking system holds 87% of its own sovereign and essentially nothing else
(F paper is 0.75% of its book). The F system is 5.9× larger and holds D paper worth
2.5% of its own book. **The exposure is radically asymmetric, and that asymmetry is
what the intervention acts on.**

### 6.2 The crisis, before TPI (γ = 0, 1pp default shock, impact)

| holder | D paper | F paper | share of D issue |
|---|---|---|---|
| **D banks** | 0.930292 | 0.007584 | **87.76%** |
| **F banks** | 0.129756 | 5.625571 | **12.24%** |
| issued | 1.060048 | 5.633155 | |

The D-paper stock loses 4.4% of market value (1.1085 → 1.0600), almost entirely price.
F banks shed D paper, D banks buy Bunds — and **the D banks' *share* of their own
sovereign rises, 87.27% → 87.76%.** Retrenchment concentrates the risk on the balance
sheet least able to carry it. This is the doom loop stated as a portfolio fact.

### 6.3 After TPI (γ = 10, same shock, impact)

| holder | D paper | F paper | share of D issue |
|---|---|---|---|
| **D banks** | 0.926341 | 0.009973 | **86.86%** |
| **F banks** | 0.102520 | 5.629347 | **9.61%** |
| **CB** | 0.037622 | 0.000000 | **3.53%** |
| issued | 1.066483 | 5.639320 | |

### 6.4 Who actually sells to the central bank

Holdings at t = 0 relative to the γ = 0 counterfactual, same shock:

| leg | Δ(γ=10 − γ=0) | % of SS D issue | share of the CB book |
|---|---|---|---|
| **F banks / D paper** | **−0.027236** | −2.457% | **72.4%** |
| D banks / D paper | −0.003951 | −0.356% | 10.5% |
| D govt issue (new supply) | +0.006435 | +0.581% | 17.1% |
| **CB / D paper** | **+0.037622** | +3.394% | 100% |
| D banks / F paper | +0.002389 | +0.216% | — |
| F banks / F paper | +0.003777 | +0.341% | — |

The clearing identity closes exactly:
`Δ(D banks) + Δ(F banks) + Δ(CB) = −0.003951 − 0.027236 + 0.037622 = +0.006435 =
Δ(D govt issue)`, to machine precision, and the three sourcing shares sum to 100.0%.

**Roughly three-quarters of the central bank's book is bought from German banks, one
sixth is new Greek issuance, and only one tenth comes off Greek banks.** In pure
quantity terms — differentiating the clearing identity, which strips the price effect
out — the German share is higher still: 84% at impact and ~99% from t = 1 onward. The
elasticity governing it is `psi_bD_F` (F-bank cross-border portfolio adjustment cost)
`= 0.5` in `gk_cross_border_foc`.

**The ECB is not sharing Greek exposure with German banks. It is buying them out of
it, at a price its own purchases are raising.**

### 6.5 What the incentive constraint sees

`phi = q_b·b / n_inter` is the object `intermediation_IC_D/F` reads, so this is where
the portfolio shift becomes a constraint effect:

| ratio | SS | γ=0 | γ=2 | γ=5 | γ=10 |
|---|---|---|---|---|---|
| `phi_bD_D` (D bank, own paper) | 0.452489 | **0.486757** | 0.479736 | 0.470614 | **0.458345** |
| `phi_bF_D` (D bank, F paper) | 0.003406 | 0.003936 | 0.004167 | 0.004463 | 0.004853 |
| `phi_bD_F` (F bank, D paper) | 0.007415 | 0.006774 | 0.006428 | 0.005983 | **0.005397** |
| `phi_bF_F` (F bank, own paper) | 0.294078 | 0.293864 | 0.294447 | 0.295200 | 0.296203 |

Three readings, and they are the core result of this document:

1. **TPI at γ = 10 undoes 82.9% of the crisis-induced rise in Greek banks' sovereign
   concentration** — the shock raises `phi_bD_D` by 0.034268, and TPI removes 0.028412
   of that.
2. **It does so almost entirely through the denominator.** Between γ = 0 and γ = 10
   the numerator `q_b_D·b_D_D` moves −0.4% (0.930292 → 0.926341) while net worth
   `n_inter_D` (D-bank net worth) recovers +6.6%. The concentration relief is a
   net-worth effect, not a quantity effect. §6.4 and §5 say the same thing from two
   directions: the CB barely takes bonds off Greek banks, and it does not need to.
3. **German banks' Greek exposure falls further the more the CB does.** `phi_bD_F`
   goes 0.007415 → 0.006774 under the shock alone and on to 0.005397 at γ = 10, i.e.
   27% below its steady-state level. Meanwhile Greek banks buy *more* Bunds
   (`phi_bF_D` up 42% from SS at γ = 10). TPI accelerates cross-border retrenchment in
   both directions rather than reversing it.

## 7. Fiscal incidence and the capital key

The conduit is real and correctly scaled. Impact remittances, F converted at `p`
(terms of trade) and scaled back to an aggregate by `size_F`:

| γ | `cb_flow_D[0]` | `rem_cb_D[0]` (D agg.) | `rem_cb_F[0]` × `size_F` | ratio F/D |
|---|---|---|---|---|
| 2 | −9.4199e−03 | −6.6881e−04 | −8.7338e−03 | 13.059 |
| 5 | −2.1539e−02 | −1.5293e−03 | −1.9971e−02 | 13.059 |
| 10 | −3.7622e−02 | −2.6711e−03 | −3.4882e−02 | 13.059 |

against `0.929/0.071 = 13.085` implied by the key; the 0.2% gap is the endogenous
terms of trade (`p = 1.00197`). Both remittances are **negative at impact** — a
capital call, which each treasury funds through its own fiscal rule at its own
sovereign terms.

Downstream, German cumulative taxes rise monotonically with γ (Σ₁₀₀ `TAX_F`, F tax
revenue: −1.4470e−03 → −1.2724e−03 from γ = 0 to γ = 10), German consumption gain
shrinks, and German welfare falls monotonically (ΔW_F = −0.0195 / −0.0482 / −0.0944 at
γ = 2/5/10) against Greek gains of +0.1956 / +0.5155 / +1.1290.

**The intervention is not self-financing, the creditor side exists, and it is signed
and scaled correctly.** Statements about German exposure are computable here.

## 8. The two profit-and-loss objects, and the rule against netting them

`writeoff_enabled_D = 0` — the S-1 framing — so `haircut_mult_D = 0` in `cb_flow_D`
and **no credit loss ever flows through the conduit**. On the branch the IRF traces,
the German treasury books a pure gain. Two different objects therefore appear in
`run_tpi`'s output and they are not commensurable:

| object | what it is | where it lives |
|---|---|---|
| `prem_pv`, `carry_ss_pv`, `mtm_pv`, `purchases_pv` | **on-path**, realised cash flows | actually pass through `budget_residual_D/F_tpi` |
| `el_pv` and the `loading` ratio built from it | **off-path**, an expectation over a default that never occurs on this branch | appear in **no** budget constraint |

At γ = 10 the printed line reads "F bears EL PV = 0.2513% `Y_D`, receives prem PV =
0.1211%". The first number is an expectation; the second is a transfer. The loading
0.52 / 0.50 / 0.48 says the CB earns roughly fifty cents of premium per euro of
*expected* loss absorbed — under-compensated in expectation — while the traced path
shows the same position ending in profit. **Both are true and they are about different
things.**

**Rule.** No sentence may net, sum, or trade off the realised German transfer against
the expected loss. Any burden-sharing claim must name which of the two it quantifies.
Same class of discipline as the standing ban on "x% fundamental / y% non-fundamental".

## 9. Stability — and a claim not to re-derive

`A_cb` = `d(spread_rb)/d(cb_buy_D)`. Prony dominant moduli are all below 1 and barely
move with γ (spread 0.9400 → 0.9470, `b_gov_D` 0.8283 → 0.9232, `n_inter_D` 0.9347 →
0.9460); `b_gov_D[499]` stays ≤ 3.7e−04; Walras residuals stay clean and
`max|ca_res_D|` actually *falls* with γ, 7.2e−08 → 1.5e−08. The closed loop is well
behaved across the whole intended range.

**There is no closed-loop pole.** `code/tpi.py` prints one at γ = 26.50 and
`CLAUDE.md` records "γ ~ 27.3"; both are wrong, and there is a nearer apparent
singularity at γ = 2.2116 that the guard's 0.25-step condition-number scan steps over
entirely. All of it is a T = 500 terminal-truncation artefact: the resonant eigenvector
carries 0.0000 of its mass in the first 100 quarters and 0.9922 in t = 400–499,
`||A_cb[:,499]|| = 3.86` against ~0.0065 for every interior column, and
`A_cb[499,499] = +1.080` is the only positive diagonal entry in the matrix. Dropping
five columns removes every pole below γ = 36 and changes the reported peak spread by
**nothing** at γ = 2, 5, 10.

No reported number is affected — every statistic is computed on `[:100]`, where the
artefact has no mass. But the γ = 19.88 cap on the effectiveness curve is imposed for a
spurious reason, and the pole should not be cited as a stability ceiling on γ. Fix
proposed at `diagnostics/cb_audit/recommended_fix.md` R-1; not implemented.

## 10. Reporting rules

Consolidated, for the paper and for anything generated into `docs/`:

1. **Never** describe TPI's spread effect and its investment effect as separate
   channels. §4 — they are one effect.
2. **Never** net the realised German transfer against the expected loss. §8.
3. Describe the relief channel as **price support and mark-to-market recapitalisation**,
   not as balance-sheet or quantity relief. §5, §6.5.
4. State that roughly three-quarters of the CB book is bought from **German** banks.
   §6.4. Do not write that the ECB absorbs bonds from Greek banks.
5. State the two absences — no CB capital, no reserve liability or policy rate — as
   modelling choices. §1.
6. Do not cite a stability ceiling on γ derived from the "closed-loop pole". §9.
7. The standing bans carry over: no "x% fundamental / y% non-fundamental", and the
   loading schedule is a measured outcome, never a target.

## 11. Provenance

| artefact | what it establishes |
|---|---|
| `diagnostics/cb_audit/run_log.md` | timestamped probe-by-probe evidence |
| `diagnostics/cb_audit/VERDICT.md` | audit findings, most severe first |
| `diagnostics/cb_audit/recommended_fix.md` | proposed fixes; **none implemented** |
| `diagnostics/cb_audit/probe_pipeline.py` | live solve: Steps 0/4/5 quantities |
| `diagnostics/cb_audit/probe_stability.py` | feedback sign, Prony moduli, fiscal incidence |
| `diagnostics/cb_audit/probe_portfolio.py` | §6, the 2×2 matrix |
| `diagnostics/cb_audit/prony.py` | order-selected Prony estimator, self-test passing |

Base audited: `EBA_CALIBRATION = True`, `BANK_SCOPE = "broad"`, all four `Delta = 0.20`,
`psi_lambda_B = 0`, `zeta_writeoff = 1`, `writeoff_enabled = 0`, `kappa_cb_F = 0.929`,
`recovery_rate = 0.30`. `lambda_gk_D/F` (GK incentive-constraint multiplier)
`= +2.1087 / +0.6414`, `Omega_D/F` (banker marginal value of net worth)
`= +10.3462 / +4.0379`, `phi_bD_D = 0.4525` against a well-posedness ceiling of 0.8232.
The GK block is well posed; this is not the inadmissible CT1-scope regime, in which
`phi_bD_D` was 2.39.
