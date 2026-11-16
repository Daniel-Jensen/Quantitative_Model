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
item is deleted. Last revised 2026-11-16.

## Code

- [ ] [joint] **Re-tune the rest-point targets at Bohn 0.15**: leverage 7.95 against 8, D
  exposure 13.2% against 14.4%, D wedge 0 against 5 bp (STATE section 1).
- [ ] [joint] **Decide whether reserves are divertible**: the code has lambda_cb_D = 0, the paper's
  incentive constraint counts them in the divertible asset base; the rule's results lean on the
  code's choice (STATE sections 1 and 3).
- [ ] [joint] **Switch off the surprise-cap step of the pipeline** (`run.py` RUN_TPI): at the 2%
  headline shock no purchase within the book band holds 200 bp, so the run stops after saving
  the no-policy results.
- [ ] [joint] **Stronger rules (k >= 4)** fail at the 27.9% book band; widen the band or leave k = 2.
- [ ] [Adam] Regenerate the 2% impulse-response figure with the lending spread measured as in the
  headline table (the figure shows 112 bp on impact, the table 115).
- [ ] [Daniel] **Fix the working-capital wedge** to zeta (r_wc - r)/(1 + r), with dividends
  carrying producers' surplus (impact output roughly 2.5x larger; backstop impact sign flips).
- [ ] [joint] **Decide the reserve remuneration rate** (national deposit rate vs a union rate).
- [ ] [Daniel] **Implement national deposit clearing in the dynamics.** The steady state
  clears nationally since PR #38 (p_bar = 0.9893), but the dynamic deposit market is still one
  union market with deposit-UIP; D consumption's rise on impact runs through its deposit-rate
  fall (STATE section 3).
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
- [ ] [joint] **Review the rule-only results section and appendix D** drafted on 2026-11-16
  (`Claude files/diagnostics/bft_rule/viva_draft/`), then paste them into Overleaf.
- [ ] [joint] **Rewrite section 3 around the proportional rule**: replace the cap problem, the
  price-floor equation and the 2×2 cap matrix (no solved counterpart) with the rule; keep the
  announcement-effect definition, reserves, run-off and profit sharing.
- [ ] [joint] **Align the introduction and calibration with the rule and Bohn 0.15**: the
  introduction's "benchmark spread cap", the Bohn 0.04 and 23-year half-life in section 4 and
  appendix B, and the [TBC] impact spread (62 to 269 bp at 2%).
- [ ] [joint] **Solution appendix**: describe the rule (continuation in k, the never-sell
  constraint, the book band, exactly cleared accuracy) instead of the cap's price-floor kink.
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
