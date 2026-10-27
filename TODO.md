# TODO

The single list of open work, in four parts: **Code** (the model and its data), **Paper** (the
viva paper's text), **VIVA** (what we must be able to defend in the examination) and
**Post-VIVA** (the research agenda after it). Within each part, items run in priority order.

Each item is a prompt; detail lives in STATE, DATA, GLOSSARY and the audits in
`Claude files/diagnostics/model_audit/`. Every item carries an owner: [joint] by default;
[Adam] or [Daniel] only where the task is clearly one author's. The tag says who owns the open
task, not who did earlier work.

Completed work never stays here. Ticking a box (`- [x]`) marks an item done but changes
nothing else; at the start of each session Claude deletes ticked items and records each in
`Claude files/docs/PROGRESS.md` with the date and the PR or commit. Last revised 2026-10-27.

## Code

- [ ] [joint] **Decide whether reserves are divertible.** If divertible, Proposition 1 fails
  (mu = mu_hat / (1 - s_R)) and a purchase tightens the constraint; recommended: non-divertible.
- [ ] [Daniel] **Fix the working-capital wedge** to zeta (r_wc - r)/(1 + r), with dividends
  carrying producers' surplus (impact output roughly 2.5x larger; backstop impact sign flips).
- [ ] [Daniel] **Implement national deposit clearing**, including the steady-state external
  balance (p_bar ≠ 1). The paper calibration (now in the code) cannot be solved without it; its
  sovereign stocks are sized at p = 1 and need re-sizing at p_bar (STATE section 1).
- [ ] [joint] **Bring the default branch onto its solved domain** (D-bank holdings 0.238 against
  a floor of 0.690): stock-and-share state space, then choose the post-default fiscal rule.
- [ ] [joint] **Decide the remittance rule and capital-account closure**: mark-to-market vs
  realised income, the 2012 SMP-income return to Greece, the untargeted psi.
- [ ] [joint] **Decide the reserve remuneration rate** (national deposit rate vs a union rate).
- [ ] [joint] **Decide HA vs RA framing**: bring the household distribution into equilibrium or
  call the model representative-agent with an overlay; fix the overlay's r_wc B_L income bias.
- [ ] [joint] Fix the small discrepancies: p_t / p_{t+1} in the cross-border bond conditions,
  fixed F taxes, gamma_tau at the risk-free price, consumer prices in the banker's discount factor.
- [ ] [Daniel] **Solve the paper calibration** once the items above are in; report fitted and
  exactly cleared output side by side (kink risk, STATE section 2).
- [ ] [joint] Run `TPI_CAP = "rest"` (defend the no-backstop rest-point spread, about 56 bp).
- [ ] [Daniel] Record sources for the unscripted calibration values (DATA.md section 4).
- [ ] [joint] Reconcile the data mismatches: lending spreads, overdraft rise, rest-point hazard,
  duration (DATA.md).
- [ ] [joint] Script the validation moments computed only in chat (DATA.md section 4).
- [ ] [joint] Decide which derived outputs in `Empirics/outputs/` to track in git.

## Paper

- [ ] [joint] **Close the backstop economy**: Eurosystem holdings M in D-bond clearing, M and
  the reserve obligation O in the state, complementarity conditions on the purchase rule.
- [ ] [joint] **Fix the government budget** (each treasury gets the whole profit) and state
  kappa (code 0.071; size share 0.081).
- [ ] [joint] **Section 3's 2×2 matrix has no solved counterpart**: implement the warranted
  floor and uncertain execution, or cut to the credible cap.
- [ ] [joint] **Reverse the funding-risk sentence**: sovereign risk relaxes the constraint via
  franchise value; show the three margins.
- [ ] [joint] **Rewrite the premium proposition**: its proof is incomplete and its premise fails
  (Omega(1)/Omega(0) = 1.29 from franchise value alone); call it a leverage-risk premium.
- [ ] [joint] **Say the backstop insures banks against default** (pari passu haircuts), and
  report outcomes including default.
- [ ] [joint] Cut the flight-to-safety claim; note the total-spread cap responds to r_D - r_F.
- [ ] [joint] Fill every [TBC] once Code is done.
- [ ] [joint] Fill the \placeholder items (income process source, HFCS moments, zeta
  sensitivity, payout share f, year of the 11.3 size ratio).

## VIVA

- [ ] [joint] **Why output falls more with the cap than without**, though bank net worth falls
  half as much (STATE section 3); re-establish the capital-key transfer explanation.
- [ ] [joint] **The credit spread is a near-cancellation** of two terms near 0.127 a quarter, so
  E[Omega] accuracy is first order.
- [ ] [joint] **The floor at the corner** depends on the untargeted psi.
- [ ] [joint] **Defend exogenous default** and the absence of fiscal moral hazard.

## Post-VIVA

- [ ] [joint] **Warranted floor**: cap(s) = EL(s) + tau, removing only the balance-sheet premium.
  Prototype on the old model: commit `bc1a83a`.
- [ ] [joint] **Welfare of backstop rules** (level vs change), redone on the paper's model.
  Record: `Claude files/docs/shelved/welfare_and_policy_rules.md`; branch
  `shelved/policy-rules-welfare`; former section 6 in `Claude files/docs/shelved/`.
- [ ] [joint] **Endogenous non-strategic default** via a stochastic fiscal limit (Bi 2012; Bi,
  Foerster and Traum 2025).
- [ ] [joint] **Bayesian estimation** of the risk process and bank block (particle filter on the
  global solution). Design note on commit `c2d2337`.
