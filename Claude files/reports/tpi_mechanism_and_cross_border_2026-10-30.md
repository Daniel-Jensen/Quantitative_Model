# Why the TPI cannot hold its cap, what purchases do, and the cross-border test (2026-10-30)

Global solver only (`global_projection/`), branch `calibration`, the final appendix-B solution: fine grid, s box to
12% p^d, b_band 0.25, k_band 0.04, f_D 0.08 (the rules `python3 run.py global` writes to
`results/GLOBAL/data/rules_risk.pkl`).

**Method.**
- Every number below is an **exact read**: the period map is cleared at the state, with the other states at the
  stochastic rest point.
- **Purchase tests** are one-off purchases with the decision rules held fixed. A credible, ongoing programme would also
  ease expected future constraints, so these numbers are lower bounds on the full TPI equilibrium. That equilibrium is
  not solved (`RUN_TPI = False`).
- **Notation:**
  - p^d: the priced default probability per quarter.
  - `mu_D`: the multiplier on the D bank's incentive constraint (IC). It is 0 when the IC is slack.
  - `E[Ω]`: the D bank's expected value of a marginal euro of equity.
  - `lambda_K`: the divertable share of loans.
  - Q_bD: the D bond price.
  - B: the D debt stock.

## 1. Why the D bank's constraint is slack in normal times

- **Calibration.** Appendix B targets a 5 bp lending spread at rest (pre-crisis Euribor–Eurepo). That puts `mu_D` at
  about 0.001 in the deterministic steady state, so the bank is barely constrained there.
- **Anticipation.** At the stochastic rest point `E[Ω]` is 1.117, about 12% above its no-risk value:
  - Banks expect future crisis states, in which equity is scarce and valuable.
  - A more valuable bank has more to lose by diverting, so it can borrow more.
  - As a result, the IC is slack by roughly 10% at rest.
- **Erosion.** The IC binds only once default risk has eroded equity enough:

| p^d | rest (0.1%) | 1% | 2% | 3% | 4% | 5% | 8% | 12% |
|---|---|---|---|---|---|---|---|---|
| D bank equity / SS | 1.01 | 0.95 | 0.90 | 0.87 | 0.84 | 0.80 | 0.71 | 0.63 |
| `mu_D` | 0 | 0 | 0 | 0 | 0.003 | 0.041 | 0.114 | 0.173 |
| `E[Ω]` | 1.117 | 1.122 | 1.139 | 1.155 | 1.170 | 1.174 | 1.200 | 1.249 |

The threshold is p^d 3.93%, and it is box-dependent (1.3–3.9% across the boxes tried). It is not grid-converged with
3 nodes per bank stock.

## 2. What the D bond's yield premium is made of

The D bond's promised yield over the deposit rate, in bp/yr (`tpi_mechanism_check.py`):
- **Expected default loss:** p^d × the 55% haircut.
- **Risk and liquidity premia:** Bocola's (2016) split of the expected excess return E[R_b] − R:
  - risk premium = −Cov(Ω, R_b)/E[Ω];
  - liquidity premium = λμ/E[Ω].

| p^d | total | expected default loss | risk premium | liquidity premium | liquidity share |
|---|---|---|---|---|---|
| rest | 35 | 22 | 13 | 0 | 0% |
| 1% | 349 | 220 | 128 | 0 | 0% |
| 2% | 687 | 439 | 248 | 0 | 0% |
| 4% | 1357 | 881 | 462 | 14 | 1% |
| 5% | 1855 | 1108 | 574 | 173 | 9% |
| 8% | 3154 | 1809 | 872 | 473 | 15% |
| 12% | 4663 | 2787 | 1183 | 693 | 15% |

The risk premium at p^d 2% (248 bp) is the +235 bp rise the pipeline's bond decomposition reports, plus the 13 bp at rest.

## 3. Why the TPI cannot hold a 200 bp cap

- **What purchases can remove.** In this model a purchase can remove only the liquidity premium: it eases the IC.
- **What they cannot touch.** Expected default loss and the risk premium are fixed by the exogenous p^d, and they are
  85–100% of the premium.
- **Price effect.** A 5% purchase raises the price by:

| p^d | rest | 2% | 5% | 8% | 12% |
|---|---|---|---|---|---|
| price rise | +0.07% | +0.41% | +0.66% | +0.97% | +1.32% |
| fall in `mu_D` | 0 | 0 | −0.009 | −0.029 | −0.047 |
| fall in the liquidity premium | 0 | 0 | −39 bp | −125 bp | −199 bp |

- **Binding vs slack.** Purchases work about three times better when the IC binds. When it is slack, they still raise
  the price a little, because they move default risk off the banks (a smaller risk premium).
- **The gap at 12%.** The 200 bp cap needs Q_bD 0.875; at 12% the price is 0.670, about 23% short. Removing the
  liquidity premium entirely would raise the price at most about 1.7%. The rest of the gap could only be closed by the
  Eurosystem absorbing the default risk, which in the model means buying most of the debt.
- **Solver evidence:**
  - On the coarse grid the cap walk solves 1500 → 1150 bp and stalls at 1100–1125 bp (`tpi_check.py`).
  - On the fine grid it barely moves at 700 bp.
- **The wall itself is not diagnosed.** `tpi_wall.py` reads the 1100 bp system at the 1150 bp solution: max|F| is 0.12,
  and the worst residuals (euler_D) sit at nodes where nothing is bought and the price is far above the floor. The
  earlier guess (the TPI book runs past its box) is therefore unconfirmed.

## 4. The lending side: Bocola's capital premia

E[r^k′] − r = λ_K·μ/E[Ω] (liquidity) − Cov(Ω, r^k′)/E[Ω] (risk), in bp/yr (`tpi_purchase_tests.py`):

| p^d | total | risk premium | liquidity premium | Y_D vs rest |
|---|---|---|---|---|
| rest | 0.1 | 0.1 | 0 | 0 |
| 2% | 1.5 | 1.5 | 0 | −0.031% |
| 5% | 173.9 | 0.7 | 173.2 | −0.135% |
| 12% | 695.6 | 2.6 | 693.0 | −0.409% |

- **No risk channel on lending.** Capital returns barely co-move with the default state, so the crisis lending spread is
  all liquidity premium. Bocola's risk channel on lending (his Table 4 object) is essentially absent here.
- **Real effects.** A purchase reaches output only through `mu_D`.

## 5. Purchases: price and real effects, and the Gertler–Karadi variant

`tpi_purchase_tests.py` compares two designs:
- **Current:** the TPI claim swaps for the bond and counts against the IC at the bond's λ (`point_map.py`, `lev_D`
  includes `lbDD·Z_cb`).
- **GK variant:** central-bank reserves are not divertable, as in Gertler–Karadi. It is applied as a diagnostic patch of
  `_incentive_constraint`; no model edit.

A one-off purchase of 5% or 15% of B, with the rules held fixed:

| p^d | current, 5% | current, 15% | GK, 5% |
|---|---|---|---|
| 2% | price +0.41%, spread 0, Y −0.002%, I +0.19% | price +1.86%, spread 0, Y −0.004%, I +0.80% | as current |
| 5% | price +0.66%, spread −38 bp, Y +0.014% | price +3.03%, spread −173 bp (all), Y +0.070% | price +0.86%, spread −173 bp, Y +0.085% |
| 12% | price +1.32%, spread −198 bp, Y +0.093% | price +5.46%, spread −693 bp (all), Y +0.321%, I +1.72% | price +1.44%, spread −272 bp, Y +0.132% |

**Real effects exist only where `mu_D` > 0.**
- **At p^d 2%** the output loss (−0.03%) does not come from the constraint, so purchases cannot repair it.
- **At 12%** the risk costs 0.41% of D output. A 15% purchase drives `mu_D` to 0 and wins back 0.32%, about 78%.
- **The ceiling.** The most purchases can ever do is remove the constraint; past that, more buying adds nothing real.

**The GK variant:**
- **Per euro:** it makes purchases work up to about 6× better on output (at 5%: +0.085% against +0.014% for a 5%
  purchase).
- **Price:** it barely moves the price (+0.2 pp).
- **Ceiling:** unchanged; at 15% both designs coincide, because `mu_D` hits 0 in both.

## 6. Options (author decisions; nothing implemented)

**To give purchases real effects**, ranked by usefulness:
1. **Crisis-sized shocks.** Evaluate the TPI at crisis p^d (Greek p^d averaged about 9%/qtr in 2010–11, inside the
   binding region). No model change. The 2% headline shock sits in the slack region.
2. **GK reserves not divertable.** One term in the IC, about 6× more output per euro; a modelling choice.
3. **A binding calibration.** A higher pre-crisis lending spread (the old 100 bp) would put 2% in the binding region.
   This conflicts with appendix B's 5 bp target.
4. **Endogenous / panic default risk.** Split p^d into a fundamental and a panic part, and let a credible backstop
   remove the panic part. This has the largest potential and is the ECB's own rationale; SPEC.md describes it for the
   SSJ model only. It is a sizeable model change.
5. **A stronger credit channel.** Today only the wage bill is bank-financed, so even a 693 bp lending spread costs about
   0.3–0.4% of output.

**To make the TPI experiment well posed:**
- **Cap the liquidity premium** instead of the whole spread, or set the cap above the actuarial spread. This matches
  the ECB's "unwarranted spreads" framing.
- **Defend the cap only up to some p^d.** This matches the ECB's eligibility conditions.
- **Use a higher cap** (400–500 bp).

**My recommendation:** options 1 + 2. They are small, consistent with the data and the literature, and give purchases a
clear real effect exactly when banks are constrained.

## 7. The cross-border bond adjustment: theory against data

The full write-up is in `Claude files/diagnostics/cross_border_test/README.md`. The model side is in `model_side.py`;
the data side (BIS CBS, ECB BSI, FRED) is in `compare.py`.

- **The term.** `psi_bD_F` (2) prices F banks' deviations from their reference D-bond book:
  (b_DF − b_ref)/B = (F banks' excess return on D bonds)/psi.
- **Prediction.** D banks set the price and value D bonds less in a D crisis (risk and liquidity premia), so F banks
  see a positive excess return and **buy**. Their book is +2.3% at p^d 2% and +18% at 12% (excess return +219 / +1812
  bp/yr).
- **Data, 2010Q1 → 2010Q4** (spread 314 → 910 bp; model p^d 2.6% → 10.9%):

| | model | data |
|---|---|---|
| German holdings, quantity | +13.5% | about −21% |
| German holdings, value | −7% | −36% |
| Greek banks' holdings, value | −18% | +14% |
| German share of German + Greek holdings | 11.6% → 13.0% | 29.6% → 19.3% |

- **Verdict.** The sign is wrong: foreign banks sold and Greek banks absorbed. A finite psi only scales the
  wrong-signed inflow.
- **What could match the data** (all modelling changes):
  - gambling for resurrection by D banks;
  - moral suasion;
  - ECB collateral value of domestic bonds;
  - a risk-dependent cross-border cost (redenomination fear).
- **Caveats:**
  - The SMP bought from foreign banks from May 2010.
  - The BIS official sector is broader than bonds.
  - The 2011 spreads lie beyond the box.
- **Not verified.** A lower psi should let F banks absorb more and soften the price fall. That needs a re-solve.

## 8. Reproduce

All scripts are read-only on the solved rules. They need `results/GLOBAL/data/rules_risk.pkl` from `python3 run.py global`.

| script | log | output |
|---|---|---|
| `Claude files/diagnostics/tex_calibration/tpi_mechanism_check.py` | `runs/tpi_mechanism_check.log` | §1–3: `mu_D`, `E[Ω]`, the yield-premium split, the 5% purchase |
| `Claude files/diagnostics/tex_calibration/tpi_purchase_tests.py` | `runs/tpi_purchase_tests.log` | §4–5: both designs, capital premia, output (~12 s) |
| `Claude files/diagnostics/tex_calibration/tpi_check.py`, `tpi_wall.py` | `runs/tpi_check.log`, `runs/tpi_wall.log` | §3: the cap walk and the wall |
| `Claude files/diagnostics/cross_border_test/model_side.py` | `model_side.log` | §7, model |
| `Claude files/diagnostics/cross_border_test/compare.py` | `compare.log` | §7, data (from `raw/`) |
