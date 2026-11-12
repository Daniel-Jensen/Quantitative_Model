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
item is deleted. Last revised 2026-11-12.

## Code

- [ ] [joint] **Reconcile the fiscal rule and re-tune the rest point**: the code's Bohn
  coefficient is 0.25 against the paper's 0.04 (Greek estimates 0.06-0.15); at 0.25 leverage
  (7.96) and D exposure (13.0% against 14.4%) miss their rest-point targets and the
  transmission row overshoots its bound (142 against 108 bp). STATE section 1.
- [ ] [Daniel] **Fix the working-capital wedge** to zeta (r_wc - r)/(1 + r), with dividends
  carrying producers' surplus (impact output roughly 2.5x larger; backstop impact sign flips).
- [ ] [joint] **Decide the reserve remuneration rate** (national deposit rate vs a union rate).
- [ ] [Daniel] **Implement national deposit clearing in the dynamics.** The steady state
  clears nationally since PR #38 (p_bar = 0.9893), but the dynamic deposit market is still one
  union market with deposit-UIP; D consumption's rise on impact runs through its deposit-rate
  fall (STATE section 3).
- [ ] [joint] **Solve the backstop as an announced rule** (the paper's framing): the reference
  cap is unanticipated purchases; the anticipated cap has not solved below about 1100 bp on
  the 12% box.
- [ ] [joint] **Make exact clearing exact**: today's franchise value, loan rate and banker
  consumption still come from the fitted rules, and fitted and exactly cleared impact output
  differ by a factor of two (STATE section 3).
- [ ] [joint] **Bring the default branch onto its solved domain** (D-bank holdings 0.238 against
  a floor of 0.690): stock-and-share state space, then choose the post-default fiscal rule.
- [ ] [joint] **Decide the remittance rule and capital-account closure**: mark-to-market vs
  realised income, the 2012 SMP-income return to Greece, the untargeted psi.
- [ ] [joint] **Decide HA vs RA framing**: bring the household distribution into equilibrium or
  call the model representative-agent with an overlay; fix the overlay's r_wc B_L income bias.
- [ ] [Daniel] **Re-solve the paper calibration** once the items above are in (the first
  solve without them is the 12 Nov reference, PR #38); report fitted and exactly cleared output
  side by side.
- [ ] [Adam] Update the constants of `bohn.py` and `risk_process.py` to the paper (committed with the old calibration's; at the paper's they reproduce appendix B's table and roots).
- [ ] [Daniel] Record sources for coupon decay (0.056, 0.039) and the 32% import share behind
  home bias 0.68, which no national-accounts ratio reproduces (DATA.md section 4).
- [ ] [joint] Recover or drop the basis-adjusted bond hazard (0.02-0.05%) in appendix B.
- [ ] [Adam] Stop `greek_cds/greek.py` appending duplicate columns to its cached series on
  re-runs (DATA.md section 6).

## Paper

- [ ] [joint] **Describe the solved economy**: the D bank is unconstrained at the rest point
  (the 5 bp wedge holds only at the deterministic steady state), the constraint binds from
  p^d of about 1.2-1.5%, and the solver selects the calm root where a self-fulfilling crisis
  root exists (the "no bad equilibrium" sentences). STATE sections 2 and 4.
- [ ] [joint] **Close the backstop economy**: Eurosystem holdings M in D-bond clearing, M and
  the reserve obligation O in the state, complementarity conditions on the purchase rule.
- [ ] [joint] **Fix the government budget** (each treasury gets the whole profit) and state
  kappa (code 0.071; size share 0.081).
- [ ] [joint] **Section 3's 2×2 matrix has no solved counterpart**: implement the warranted
  floor and uncertain execution, or cut to the credible cap.
- [ ] [joint] **Rewrite the premium proposition**: its proof is incomplete and its premise fails
  (Omega(1)/Omega(0) = 1.29 from franchise value alone); call it a leverage-risk premium.
- [ ] [Adam] **Say the backstop insures banks against default** (pari passu haircuts), and
  report outcomes including default.
- [ ] [Adam] Cut the flight-to-safety claim; note the total-spread cap responds to r_D - r_F.
- [ ] [Adam] Fill every [TBC] once Code is done.
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
