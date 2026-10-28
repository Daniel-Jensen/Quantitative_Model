# TODO

The single list of open work, in three parts: **Code** (the model and its data), **Paper** (the
viva paper's text) and **Post-VIVA** (the research agenda after it). Within each part, items run in priority order.

Each item is a prompt; detail lives in STATE, DATA, GLOSSARY and the audits in
`Claude files/diagnostics/model_audit/`. Every item carries an owner: [joint] by default;
[Adam] or [Daniel] only where the task is clearly one author's. The tag says who owns the open
task, not who did earlier work.

Completed work never stays here. Ticking a box (`- [x]`) marks an item done but changes
nothing else; at the start of each session Claude deletes ticked items and records each in
`Claude files/docs/PROGRESS.md` with the date and the PR or commit. Pushing a tick to `main`
notifies the co-author on GitHub (both authors for a [joint] item), so tick and push before the
item is deleted. Last revised 2026-10-28.

## Code

- [ ] [joint] **Decide whether reserves are divertible.** If divertible, Proposition 1 fails
  (mu = mu_hat / (1 - s_R)) and a purchase tightens the constraint; recommended: non-divertible.
- [ ] [Daniel] **Fix the working-capital wedge** to zeta (r_wc - r)/(1 + r), with dividends
  carrying producers' surplus (impact output roughly 2.5x larger; backstop impact sign flips).
- [ ] [Daniel] **Implement national deposit clearing**, including the steady-state external
  balance (p_bar ≠ 1) and D banks' holdings of F bonds (the paper is silent; the six Greek EBA
  banks held Bunds worth 0.11% of their assets). Precondition for the paper calibration: the
  steady state imposes one union deposit market, which clears only a symmetric economy.
- [ ] [Adam] **Put the paper's calibration into the code** (`calibration/global_projection.py`;
  rho_s is hard-coded 0.95 in the grid builder). Values: STATE section 1. Propose the changes
  for approval first (a first attempt, 7c683bb, was reverted). Size the sovereign stocks to the
  own-sovereign shares (in closed form at p = 1, D debt is 52.6% of annual GDP against the 40%
  EBA banks held) and state the tax rule by its Bohn coefficient (gamma_tau = 0.01).
- [ ] [joint] **Bring the default branch onto its solved domain** (D-bank holdings 0.238 against
  a floor of 0.690): stock-and-share state space, then choose the post-default fiscal rule.
- [ ] [joint] **Decide the remittance rule and capital-account closure**: mark-to-market vs
  realised income, the 2012 SMP-income return to Greece, the untargeted psi.
- [ ] [joint] **Decide the reserve remuneration rate** (national deposit rate vs a union rate).
- [ ] [joint] **Decide HA vs RA framing**: bring the household distribution into equilibrium or
  call the model representative-agent with an overlay; fix the overlay's r_wc B_L income bias.
- [ ] [Daniel] **Solve the paper calibration** once the items above are in; report fitted and
  exactly cleared output side by side (kink risk, STATE section 2).
- [ ] [joint] Run `TPI_CAP = "rest"` (defend the no-backstop rest-point spread, about 56 bp).
- [ ] [Adam] Update the constants of `bohn.py` and `risk_process.py` to the paper (committed with
  the old calibration's; at the paper's they reproduce appendix B's table and roots).
- [ ] [joint] Re-solve the reference (Bocola) calibration to measure the three
  equilibrium-condition fixes of 2026-10-27 (terms of trade in the cross-border bond
  conditions, F budget balance, gamma_tau at the traded price; STATE section 4).
- [ ] [Daniel] Record sources for coupon decay (0.056, 0.039) and the 32% import share behind
  home bias 0.68, which no national-accounts ratio reproduces (DATA.md section 4).
- [ ] [joint] Recover or drop the basis-adjusted bond hazard (0.02-0.05%) in appendix B.
- [ ] [Adam] Stop `greek_cds/greek.py` appending duplicate columns to its cached series on
  re-runs (DATA.md section 6).

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
- [ ] [joint] State that the banker's discount factor omits the ratio of consumer price indices
  (a choice; of order the import share times the terms-of-trade movement).
- [ ] [joint] Fill every [TBC] once Code is done.
- [ ] [joint] Fill the \placeholder items (income process source, HFCS moments, zeta
  sensitivity, payout share f; the 11.3 size ratio is German over Greek GDP in 2006).

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
