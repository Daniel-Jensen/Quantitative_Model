# The Bi-Foerster-Traum purchase rule as an announced commitment (2026-11-12)

Global solver, PR #38 calibration and solution (`results/GLOBAL_2026-11-12_pr38` in the main repo), branch
`new-rules` (renamed from `bft-rule` 2026-11-13) off 0c80116. Output: `results/GLOBAL_2026-11-12_pr38_bft` (main repo). Not quotable in the paper.

## The rule, checked against BFT (FRBSF WP 2025-10, eq. 2.29, A.73-A.75, Table A.1)

- BFT: `Q_t b_t^cb = phi_cb (ER_t - ER_ss)`, a STOCK at market value; `ER_t = E_t R^b_{t+1} - R^d_t` is the
  quarterly expected excess return net of the haircut; `phi_cb = 16 Q_ss b_ss`. Holdings in % of the SS market
  value of debt = 16 x the quarterly excess-return gap (100 bp/yr = 4% of the debt).
- Here: `Q_bD M' = k Q_bD_ss B_ss (ER - ER_ref)`, `ER = E[pay_D']/Q_bD - 1 - rdep_D` (physical measure, haircut
  in default), `ER_ref` = the no-policy rest point's 15.38 bp/yr. Purchases `m = max(target - runoff, 0)`: the
  Eurosystem buys up to the target and never sells (BFT's linear rule sells when ER falls); `m` is capped at the
  privately held stock carried into the quarter, `(1-delta)(b_DD + b_DF)`.
- Financing differs from BFT by instruction: BFT rebates `R^b Q_{t-1} b_{t-1} - Q_t b_t` to households; here the
  existing treatment is kept (a non-divertible claim on the Eurosystem held by D banks at the D deposit rate, P&L
  remitted by the capital key, 7.1% to D).
- Clearing `B' = b_DD + b_DF + M'` holds to 1e-17, the claim equals the book's value exactly, and the F goods
  market (Walras) closes to 3e-10 on the crisis paths.

## Solves (pre-specified k = 2, 4, 8, 16; continuation in k from the no-policy rules)

- k = 2: converged with all complementarities hard, max|F| 1.2e-8 (44 min).
- k = 4: three attempts (from 2, 3 and 3.5, which all solved) stall at max|F| 3.2e-3 in the D franchise-value
  identity. Run stopped after the third failure; 8 and 16 not attempted. (The halving rule as coded restarts its
  failure count after each successful intermediate step; the stop at three failures was imposed during the run.)
- Cause (numerical domain): at the riskiest node (p^d 12%/qtr) the rule asks for 14.5% of the debt at k = 2,
  25.4% at 3.5, 29.0% at 4, 58% at 8 and 116% at 16, against a box band for the Eurosystem book of 27.9%. The
  failure starts where the book first leaves the box; beyond it the continuation reads the excess as bank-held
  bonds, so no consistent fixed point exists. Reaching k >= 4 needs a wider book band (a box setting, with the
  no-policy baseline re-solved on the wider grid); k = 16 also needs the available-stock bound to bind.

## k = 2: results (exact one-quarter clearing; differences are policy minus no policy)

Rest point (announcement): sovereign spread 61.5 -> 59.3 bp; D bond price +0.09%; D output +0.04%, investment
+0.10%, D consumption +0.03%, F consumption +0.001%; D net worth -0.13%; D constraint slack in both; the book holds
0.02% of the debt (rest-point ER 15.9 bp, 0.5 bp above the trigger).

1.5% shock, impact: purchase 1.06% of the debt (none later; the book runs off); spread 222 -> 214 bp; bond price
-6.2% -> -6.0%; D net worth -10.0% -> -8.0%; mu_D 0.015 -> 0; lending spread 70 -> 0 bp; investment -1.49% ->
-0.87%; output -0.076% -> -0.023%; D deposit rate 89 -> 114 bp. Over 25 quarters: summed output deviation -0.94 ->
-0.54, D consumption -0.36 -> -0.61, F consumption +0.08 -> +0.07.
CORRECTED 2026-11-12 (forensic section below): the earlier "not identified" verdict compared the policy's FITTED read
with the no-policy EXACT read. Like for like the lending-spread cut is robust (fitted 146 -> 66 bp, exact 70 -> 0,
consistent 49 -> 0, exact continuation 98 -> 25); only whether mu reaches exactly zero is not identified.

6% shock, impact: purchase 5.0%; spread 573 -> 532 bp; mu_D 0.099 -> 0.039 (fitted 0.052); lending spread 426 ->
169 bp; investment -3.35% -> -2.14%; output -0.29% -> -0.14% (fitted -0.23%); D net worth -25.8% -> -21.2%. D
consumption summed over 25 quarters -1.45 -> -3.76 and F +0.20 -> +0.33, but the no-policy path hits the box wall
in quarters 1-4, so only the impact comparison is clean.

Accuracy (ergodic set, 1,200 periods): Euler errors match the no-policy rules equation by equation (mean
log10 -3.1 to -3.9) except the fitted purchase rule (mean 0.16%, max 2.4% of the debt) and the D franchise value
(max 2.6% against 1.0%). The book never leaves its band at k = 2; other states leave the box next quarter at 11
nodes (25 without the policy); simulated box escapes 0.7% of periods (1.5% without).

## Forensic decomposition of k = 2 (2026-11-12; no solve, one-quarter clearings on the saved rules)

Scripts in `forensic/`, output in `results/GLOBAL_2026-11-12_pr38_bft/forensic/` (main repo). Common state = the
no-policy rest point with s shocked; A no policy, B the same purchase as a surprise (no-policy expectations), C the
announced rule with no purchase this quarter, D the rule.

- Lending spread (bp), common state: 1.5%: A 70, B 72, C 0, D 0. 6%: A 426, B 400, C 169, D 121; own rest points
  D 169 (the policy rest point adds +48: a deposit bill 0.9% of net worth larger). The credit effect is the
  announcement, not the purchase.
- Shapley split of the exact identity mu = max(1 - E[Om](1+r)n/lev, 0), common state, 1.5% / 6% (bp): capital-price
  revaluation of the legacy capital stock -68 / -218; sovereign capital gain -7 / -58; swap into the non-divertible
  claim -3 / -24; E[Om] (SDF + franchise) +3 / -26; deposit rate -1.5 / -11.5; larger divertible base +8 / +36.
- Bond price (%), common state, 1.5% / 6%: +0.35 / +1.50 = continuation price +0.23 / +0.89, GK premium +0.17 / +0.76,
  covariance +0.01 / +0.13, deposit rate -0.07 / -0.29, default-warranted value 0 / 0 (pi exogenous).
- Net worth (%), common state: +2.6 / +7.2 = capital price +2.24 / +5.64 (Q_K +0.35% / +0.73%), sovereign gain +0.32
  / +1.50, swap 0 (exact).
- One-quarter equilibrium unique (capital Euler monotone in K'); mu falls ~0.12 per 0.1% of K'.
- Divertible claim (lambda_cb = lambda_bD) at D's allocation: 1.5% still slack (psi -0.26% -> -0.11%); 6% mu 0.039
  -> 0.045, spread 169 -> 194 bp (210 re-cleared).
- Accuracy: next-quarter franchise value alpha' fitted with 2-5% error near the kink; replacing the continuation by
  exactly cleared values moves the levels (1.5%: 98 -> 25 bp; 6%: 524 -> 294) but not the size of the cut. Fitted
  today's lending rate is off 77-110 bp at the 1.5% node. No current or next state off the box, except the no-policy
  6% next quarter (deposit bill, 1.8% of the width). Swap, claim, clearing and balance-sheet identities exact.
- Fixed-investment counterfactual inconclusive: holding K' forces a 130-170 bp fall in the D deposit rate.

## Why the price of capital moves (2026-11-12, part 3; one- and joint two-quarter clearings, no solve)

- Exact pricing: Q_K = E[mpk' + (1-delta)Q_K'] / (1+r) x cov x E[Om](1+r)/(E[Om](1+r) + lambda_K mu); identity holds
  to 1e-15. Common state, 1.5% / 6%: Q_K +0.35% / +0.73% = resale value +0.24 / +0.25, constraint premium +0.17 /
  +0.76, deposit rate -0.07 / -0.29, dividends 0 / 0, covariance 0 / +0.01.
- Source (continuation values switched to the policy rules one group at a time, today re-cleared, Shapley), Q_K
  points 1.5% / 6%: next-quarter capital accumulation +0.14 / +0.25, D terms of trade (via the deposit rate) +0.12 /
  +0.14, D consumption (SDF) +0.05 / +0.14, franchise value 0 / +0.07, D bond price +0.01 / +0.05, today's purchase
  +0.01 / +0.07. Impulse at the no-policy allocation +0.27% / +0.53%; equilibrium amplification 1.28 / 1.38.
- Revaluation: 0.349% x (1-f)(1-delta)Q_K K/n = 6.41 -> 2.234% of n; mpk +0.010, entry +0.025 separate.
- dmu/dlnQ_K = -5.6 / -6.3 (accounting, binding side); dlnQ_K/dlnK' = xi Q_K K'/I = 20. The "-0.12 per 0.1% of K'"
  is a secant over [-0.1%, 0] at the no-policy 1.5% state; 0.1% of K' is 4% of investment and 7.5x the distance
  to the kink (0.013%). Loop gain of the capital demand price in itself 0.65.
- Joint two-quarter clearing: Q_K effect 1.5% +0.26 / +0.30% (own / common), 6% +0.57 / +0.68%; spread cut 1.5% 70.5
  -> 9.5 bp, 6% 516 -> 290 bp; next-quarter franchise value 1-3% lower than fitted.

## Viva draft and the two open calibration issues (2026-11-13)

The k = 2 capital-price audit is treated as provisionally complete. A paste-ready draft subsection
for 05-results, presenting the announced proportional rule as an extension and robustness check on
the 200 bp cap, is in `05-results_proportional_rule_DRAFT.tex` (not pasted; review round open).
One correction to the part-3 report: the bond-market contribution to the 1.5% capital-price rise is
0.024 points (bond continuation 0.009, franchise 0.002, impact purchase 0.013), reported as 0.03.

**Bohn 0.15: minimum work.** The no-policy rules at Bohn 0.15 already exist, solved on the same
economy (12% box, recapitalisation floor 0.5) as a rung of the 0.06 -> 0.25 continuation of
2026-11-06: `Claude files/diagnostics/tpi_root_cause/runs/rules_12pct_bohn0.150.pkl`, max|F|
1.1e-13, in git. Debt root at 0.15 (gamma = 0.0375): 0.944 + (0.056 - 0.0375)/0.947 = 0.9635,
half-life about 19 quarters (0.25: 0.937, 11 q; 0.06: 0.987, about 60 q); the rule is explosive
only below q = 1 - gamma/delta_b = 0.33. Minimum test, no code change (calibration passed as
`dict(cal, bohn_annual=0.15)` in a diagnostic script):
1. on the saved rules (minutes, no global solve): rest point, the default probability at which the
   D constraint starts to bind, impact at 1.5% and 6%, exact one- and two-quarter clearing;
2. surprise 200 bp cap on the 0.15 no-policy rules (minutes);
3. announced k = 2 rule at 0.15: one policy solve, continued 0.25 -> 0.20 -> 0.15 from the k = 2
   rules (k = 2 took 44 min cold; budget 1-2 h);
4. rest-point targets will move (already missed at 0.25: leverage 7.96 vs 8, D exposure 13.0% vs
   14.4%); minimum is to report them; re-targeting is the multi-solve rest-point loop.

**Deposit clearing: minimum work.** The code is internally consistent: with the debt-elastic
premium the zero cross-border deposit position is the only rest point, so national clearing in the
steady state and union clearing with deposit parity in the dynamics are one model. The mismatch is
with the paper, which clears deposits nationally every period. Options:
- keep union clearing and rewrite the paper's deposit market (no computation);
- approximate national clearing without a code change by walking kappa_nfa from 0.01 to a
  prohibitive value (0.1, 1, 10), each rung a warm no-policy solve, then re-solve the policies.
  Numerical risk: the position V_dep collapses towards zero and sits in the aligned box block with
  K_D, P_D and b_DD;
- exact national clearing: a one-equation change (deposit parity replaced by nfa_dep_D = 0, behind
  a switch like `union_nominal_rate`), needs approval; keep V_dep in the state so warm starts survive.

**Which is the greater risk to the main mechanism: the Bohn coefficient.** Holding everything
else fixed, moving it from 0.06 to 0.25 moved the default probability at which the D constraint
binds from about 5% to 1.2-1.5% a quarter (PROGRESS 2026-11-07). The 1.5% headline shock engages
the constraint, and every headline number in 05-results and in the draft depends on that, only
under the faster rule, which is five to eight times the euro-area estimates (0.03-0.05). The
channel through which the fiscal rule moves the threshold has not been decomposed. Deposit
clearing changes the strength of transmission (the deposit-rate crowding that limits the capital
accelerator, and the terms-of-trade term worth 0.12 of the 0.35-point capital-price rise at 1.5%)
and the cross-country incidence, but nothing yet suggests it decides whether the constraint binds.
The two interact: at Bohn 0.06 the slowest modes of the exact law of motion load on the
cross-border deposit position (persistence_check.log; relative units overweight a position near
zero at rest, so this is suggestive only).

## Bohn 0.15 against 0.25 on the saved no-policy rules (2026-11-13)

Scripts `bohn015/bohn_compare.py`, `bohn_why.py`, `bohn_curve.py`; output
`results/GLOBAL_2026-11-12_pr38_bft/bohn015/` (main repo). No global solve; rules of the 2026-11-06 continuation
(0.15: max|F| 1.1e-13). The 0.25 rung reproduces the PR #38 solution (rest n 1.4814, mu 0; 1.5%: 69.7 bp; 6%: 425.6 bp).

- Binding threshold (exact one-quarter / two-quarter / fitted), quarterly p^d: 0.25 1.11 / 1.12 / 0.70%;
  0.15 1.31 / 1.30 / 0.78%. Along the continuation (exact): 0.06 4.43%, 0.10 1.88%, 0.15 1.31%, 0.20 1.16%, 0.25 1.11%.
- Rest point (exact), 0.25 vs 0.15: debt 1.897 vs 1.928 (47.4% vs 48.2% of annual GDP), D-bank bonds 1.661 vs 1.692,
  exposure 13.0% vs 13.2%, n 1.481 vs 1.485, leverage 7.96 vs 7.95, mu 0 both, psi -0.038 vs -0.044, E[Omega] 1.030 vs 1.034.
- Impact, exact (two-quarter): 1.5%: lending spread 69.7 (70.5) vs 34.9 (37.2) bp; output -0.076 (-0.071) vs -0.055
  (-0.052)%; investment -1.49 vs -1.36%; bond price -6.22 vs -6.23%. 6%: 426 (516) vs 438 (510) bp; output -0.290 vs
  -0.292%. Fitted at 1.5%: 146 vs 133 bp (factor 2 and 4 off).
- Surprise 200 bp cap at 1.5%: 16.7% vs 17.1% of the debt; spread to 20 vs 0 bp; output -0.076 -> -0.045 vs -0.055 ->
  -0.032%. At 6% the cap is infeasible in both (expected-loss spread above 200 bp).
- Off-grid accuracy (1200 simulated quarters): mean log10 error -4.37 vs -4.34, max -1.98 vs -1.87; the simulated set
  never binds (p^d p95 0.44%), so crisis accuracy rests on the one- vs two-quarter gap. Franchise-value continuation
  error 2.0-2.4% at 1.5%, 3.3-3.5% at 6%. The 6% exact IRF sits on the box edge for 8 (0.25) and 6 (0.15) quarters.
- Why: at a common state (1.5%), 85% of the fall in psi comes from the rules (expectations, tax response), 15% from the
  rest-point balance sheets. The rules effect is all in E[Omega] (+0.42%), all from next quarter's franchise value
  E[alpha'] (+0.50%); the household discount factor moves -0.02%. Slower rule: taxes 0.5% lower and issuance 1% higher
  in the shock quarter, debt more persistent (root 0.9635 vs 0.937), banks expect to carry a larger book of bonds
  earning a 2.6-2.9%/yr excess return; their continuation value rises and the constraint binds later. Rest-point
  E[Omega] rises convexly as the root approaches one (1.030, 1.031, 1.034, 1.046, 1.111 from 0.25 to 0.06).

## The announced k = 2 rule at Bohn 0.15, headline shock 2% (2026-11-13)

Proposed baseline (not yet in the calibration file): Bohn 0.15, headline p^d 2%/qtr, 1.5% kept as the lower-stress
case. Scripts in `bohn015/k2/`; output `results/GLOBAL_2026-11-12_pr38_bft/bohn015/k2/` (main repo).

- Shock size (the code's own logistic and AR(1), rho_s 0.9, sigma_s 0.41, rest p^d 0.1%): 2% is 3.2 unconditional sd
  above rest (7.4 innovation sd in one quarter), stationary tail 0.067% (once in ~370 years); 1.5% is 2.9 sd (6.6
  innovation sd), tail 0.19% (~130 years).
- Solve: from the saved Bohn 0.15 no-policy rules, k 0 -> 1 -> 2 plus hard polish, 39 min, max|F| 1.3e-8 at the
  nodes; book peaks at 26.4% of the debt against the 27.9% band (as at 0.25). Rule reference ER 15.45 bp (15.38 at 0.25).
- ACCOUNTING/NUMERICAL TRAP: the fitted rest-point iteration took the Eurosystem book to -3.4% of the debt, where the
  rule switches itself off (tpi requires M_lag >= 0), and the exact polish kept that point. The first comparison
  (`compare_k2_VOID_negative_book_rest.*`) is void for every own-rest-point number. The true rest point is the exact
  fixed point from a zero book (book 0.022% of the debt); the exact no-shock path converges to it.
- Announcement at rest (exact = 2q): sovereign spread -1.7 bp, output +0.04%, investment +0.10%; mu 0 either way.
- Impact, exact (two-quarter), Bohn 0.15, from each economy's own rest point:
  no policy 1.5%: lending spread 35 (37) bp, output -0.055 (-0.052)%; 2%: 115 (119) bp, -0.103 (-0.100)%.
  rule 1.5%: 0 (0) bp, -0.020 (-0.020)%, buys 1.06%; 2%: 0 (8) bp, -0.024 (-0.027)%, buys 1.42 (1.49)%.
  Incremental crisis effect at 2%: lending spread -115 (-110) bp, output +0.079 (+0.073) pp, sovereign spread -13
  (-11) bp, D bank net worth +2.9 (+2.5) pp, investment +0.84 (+0.73) pp. Cumulative 24-quarter output loss 1.12 ->
  0.64 (exact IRF).
- Thresholds (exact / 2q): no policy 1.31 / 1.30%; rule 2.39 / 1.92% (at 0.25: 1.11 / 1.12% and 1.74 / 1.42%). At 2%
  the no-policy economy binds robustly; the policy economy sits at its threshold, which decides only whether the
  residual spread is 0 or 8 bp. At 2.5% the rule leaves 7 bp, at 3% 37 bp (exact).
- At 0.25 the 2% cut was -120 (exact) / -85 (2q) bp because the rule left 21 / 63 bp; 0.15 is the more robust case.
- Fitted rules overstate every spread (no policy 2%: 201 bp; rule 2%: 75 bp): not quotable. Off-grid (simulated 1200
  quarters, p^d max 2.2%): mean log10 Euler error -3.6, max -1.6, box never left. Fitted next-quarter franchise value
  off by up to 3.4%.
- D consumption rises on impact with and without the rule (+0.13% at 2%, no policy), as in PR #38; pre-existing,
  tied to the fall in the D deposit rate under union deposit clearing; not investigated here.
- Earlier draft table purchase row corrected: 1.1% / 4.7% of the debt (it had reported the bond quantity, 0.0199 /
  0.0868, as a percentage).

## Canonical calibration: Bohn 0.15, headline 2% (2026-11-13, approved)

`calibration/global_projection.py` bohn_annual 0.25 -> 0.15; `run.py` RISK_SHOCK_PD 0.015 -> 0.02 (both uncommitted
on new-rules). Old worktree run kept as `results/GLOBAL_2026-11-12_pre_bohn015_announced_cap`, old solve cache as
`results/GLOBAL_cache_bohn025`. Canonical outputs: `results/GLOBAL_2026-11-13_bohn015/` (main repo): `no_policy/`
(the pipeline run), `rules_bft_k2.pkl`, `headline.json`, `book_sign.json`, `accuracy_exact.json`. Tables and figures:
`bohn015/canonical/viva/` (make_viva.py).

- No-policy pipeline from the 0.25 warm start: 1.2e-9 at the nodes; reproduces the continuation rules of
  tpi_root_cause to 8e-9 (relative). Its step 6 (surprise 200 bp cap at 2%) stops: no purchase up to 60% of the debt
  reaches 200 bp, because the fair (risk-neutral) spread at 2% is 243 bp (207 at 1.5%, 494 at 6%). Economic, not
  numerical: a cap below the fair spread asks the Eurosystem to pay more than the expected payoff.
- Announced k = 2: 1.3e-8, 41 min; reproduces the diagnostic k = 2 to 4e-11. Rule reference ER 15.449 bp.
- Headline numbers identical to the diagnostic section above (115 -> 0 (119 -> 8) bp at 2%, etc.).
- Identities at every reported state (rest, 1.5%, 2%, 6%; one- and two-quarter): net worth 2e-15, D bond market
  1e-11, mu = max(psi,0) 3e-11, lending spread 2e-18, Eurosystem claim 0, F goods 4e-10, book law 1e-20.
- Eurosystem book: >= 0 along every exact IRF (min 0.022% of the debt) and on an exactly cleared 1,200-quarter
  simulation (0.010-1.75%). At the nodes, 3/75 carry a negative book in (the box is symmetric in M); every node with a
  non-negative book in ends non-negative. The FITTED simulation (accuracy.simulate) turns the book negative in its
  first quarter and the rule then switches itself off: 89% of that sample (87% at 0.25) is an economy without the
  backstop, so the solver's own accuracy report for any BFT solution, including the 2026-11-12 0.25 one, is not an
  accuracy report for the policy economy. Replaced by accuracy_exact.py: mean / max log10 Euler error -4.4 / -1.9
  (no policy), -3.8 / -1.7 (rule); worst equations the D franchise value and the purchase rule. The ergodic set binds
  in 0.1% (no policy) and 0% (rule) of quarters, so this measures accuracy mostly away from the crisis region.
- Purchases normalised by the steady-state face value of D debt (B_gov_D_ss 1.857); at 2% 1.42% of it, 1.38% of the
  current debt, market value 1.2% of the s.s. market value (0.25 checks: 1.06 / 1.03 / 0.97 at 1.5%).

## Correction: why the 200 bp surprise cap clears at 1.5% but not at 2% (2026-11-16; arithmetic on saved results, no run)

The 2026-11-13 line "a cap below the fair spread asks the Eurosystem to pay more than the expected payoff" is wrong as
stated. The fair spread (Q_rn = E_payD/(1 + rdep_D)) is evaluated at the allocation WITHOUT the purchase; the purchase
raises next quarter's bond price and so moves its own bound. At 1.5%: fair 207.1 bp, the cap needs a price 0.28% above
Q_rn; the 17.1% purchase lowers next quarter's spread 188 -> 171 bp (Q' +0.66%), E_pay up ~0.6%, post-purchase fair
spread ~191 bp (holding the F yield and rdep_D fixed; approximation). At 2%: fair 243.4 bp, the cap needs 1.69% above
Q_rn (six times the gap). The "no purchase up to 60%" search is reliable only up to the solved book band (~28% of the
debt); above it next quarter's book is outside the solution box, so infeasibility beyond ~28% is untested.

## Capital-price mechanism on the canonical solution, 2% and 1.5% (2026-11-16; clearings on saved rules, no solve)

Scripts `bohn015/canonical/forensic/` (copies of `forensic/` on the canonical files; policy rest = exact zero-book
fixed point from headline.json). Output `results/GLOBAL_2026-11-13_bohn015/forensic/` (main repo). Reproduces the
headline (A 114.79 bp; D 0 bp, buys 1.418%). Common state = no-policy rest point with s shocked. A no policy; B the
rule's purchase (1.43%) as a surprise under no-policy rules; C announced rule, zero purchase today; D rule.

- Q_K (%), common state, 2% [1.5%]: D +0.458 [+0.321] = expected resale value +0.281 [+0.267], constraint premium
  +0.287 [+0.087] (= the no-policy lending spread, 114.8/4e4: the slack constraint removes it), deposit rate -0.112
  [-0.035], dividends 0.000, covariance +0.001. Identity error 5e-15. C +0.443 [+0.308] (97% [96%] of D);
  B +0.015 [+0.012]. 2q: D +0.40, C +0.40, B +0.02 (2%).
- Lending spread (bp), common state 2%: A 115, B 118, C 0, D 0 (2q: 119, 114, 0, 0). 1.5%: A 35, B 37, C 0, D 0.
  Surprise purchase's own effect is a few bp of either sign depending on clearing: not identified as nonzero.
- Net worth, common state 2%: +3.60% of A's n = capital +3.015 (legacy capital revaluation 0.458% x 6.55 = 3.000,
  identity; mpk +0.015) + sovereign capital gain +0.554 + entry +0.042 - F bonds 0.008. 2q: +3.14% (capital 2.63,
  sovereign 0.48). Own rest points: +3.30% of A's shocked n = capital +3.48, sovereign gain +0.54, legacy bond
  quantity -0.16, deposits owed -0.62, Eurosystem claim +0.03, entry +0.05; the state shift (policy rest vs no-policy
  rest, policy rules) is -0.29%: more debt-financed capital and fewer bonds carried in. 3.30% of the shocked n =
  0.0429 = 2.89% of rest n: the headline +2.9 pp. The "+4.3%" in the 2026-11-16 revision plan was the level 0.043,
  mislabelled as a percentage.
- Shapley of the 115 bp cut (common, 2%): capital revaluation -107, sovereign capital gain -12.5, swap into
  non-divertible reserves -4.5, larger divertible base +9.1, banker SDF +3.7, franchise -2.7, deposit rate -2.6.
- Source of the Q_K rise (continuation groups switched to policy, today re-cleared; Shapley, 2% common): next
  quarter's capital and hours +0.233, households / terms of trade / F side +0.200, bond price +0.015, franchise
  +0.002, today's purchase +0.017. Impulse at A's allocation +0.335%, equilibrium +0.457% (amplification 1.36);
  dmu/dlnQ_K = -5.66 at A, so +0.458% x 5.66 = 0.026 > mu_A 0.025: slack.
