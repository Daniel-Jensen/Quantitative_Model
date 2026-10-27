# TODO

The single list of open work. Sections 1-7 block the viva paper, in order; section 8 is the
long-term agenda. Each item says what is open and why it matters; detail lives where linked.
Every item carries an owner: [joint] by default; [Adam] or [Daniel] only where the task is
clearly one author's (his own calibration values, his own local setup or branch). The tag says
who owns the open task, not who did earlier work.

Completed work never stays here: delete the item and record it in
`Claude files/docs/PROGRESS.md` with the date and the PR or commit. Git records who changed
what; PROGRESS adds the economics. Last revised 2026-10-27.

Symbols: `GLOSSARY.md`. Data provenance: `DATA.md`. Current results:
`Claude files/docs/STATE.md`. The two audits cited below are
`Claude files/diagnostics/model_audit/THEORETICAL_REAUDIT.md` (14 Oct) and
`RESERVES_DEPOSITS_CB_REVIEW.md` (15 Oct); the first audit's reports were deleted, and its only
record is the model-section-audit memory note, summarised here.

## 1. Bring the global code to the paper's model (blocks every model number in the paper)

Every saved global result belongs to an economy the paper no longer describes (union deposit
clearing, leverage 5, a 100 bp wedge). Until this section is done the paper's model values stay
[TBC]. Comparison table: STATE section 1.

- [ ] [joint] **National deposit clearing.** Each bank funded only by its own households
      (dep_X = P^c_X HA_X), replacing the union market with real deposit-rate parity. Open: the
      steady-state external balance is unspecified (p_bar, the terms of trade at rest, is not 1
      under the asymmetric calibration), and D banks' holdings of F bonds were never calibrated.
- [ ] [joint] **Paper calibration values** into `calibration/global_projection.py`: leverage 8, wedge
      5 bp/yr, own-sovereign 14.4% / 6.5% of bank assets, foreign share 12.7%, size 11.3,
      rho_s 0.9 (hard-coded 0.95 in the grid builder), sigma_s 0.41, sigma 1.35, Frisch 0.5,
      home bias 0.68, trade elasticity 1.5, epsilon 11, G 0.204 / 0.191, delta_b 0.056 / 0.039,
      annual Bohn coefficient 0.04.
- [ ] [joint] **The model's open issues (section 2) that change equations** must be settled first:
      reserves, the working-capital wedge, the default-branch domain, the backstop closure.
- [ ] [joint] **Solve the paper calibration (needs an explicit go-ahead; about 4-5 hours with the
      backstop).** Then check the kink risk: report fitted and exactly cleared output responses
      side by side; if mu (the multiplier on the bank constraint) is near zero at the rest
      point, the output response is not numerically identified (STATE section 2).

## 2. Open issues in the model and the paper's algebra

Each is classified as an accounting error (the calculation does not represent the intended
economics), an error of interpretation (equations right, text wrong), or a modelling choice that
needs a decision or a justification.

**Households**
- [ ] [joint] **Heterogeneous-agent versus representative-agent framing (inconsistency).** The solved
      aggregates are representative-agent: the heterogeneous-household block enters no aggregate
      equilibrium condition, so it cannot affect prices or quantities. Either bring it into the
      equilibrium (deposit supply from the household distribution, which under national
      clearing must also price D deposits alone) or describe the model as representative-agent
      with a distributional overlay. Until then, never describe the solved dynamics as HANK.
      The overlay also counts r_wc B_L (working-capital interest) as household income without
      the repayment: a 0.09-point bias in consumption, the order of the distributional effects.

**Banks and reserves**
- [ ] [joint] **Authors' decision: are reserves divertible?** Under the paper's specification
      (divertible at lambda, the divertible fraction, and paid the deposit rate) Proposition 1
      (the closed form for mu) fails once reserves are outstanding: mu = mu_hat / (1 - s_R),
      with s_R the reserve share of assets, so a purchase tightens the constraint on impact.
      If reserves are non-divertible, D banks become unconstrained for about five years after
      buying about 11% of the stock. The composition effect is +1 bp/yr of lending spread per
      point of reserves if divertible, -50 bp/yr if not. Recommended: non-divertible (Gertler
      and Karadi 2011; Sims and Wu 2021). Reserves are also a forced holding while mu > 0, which
      the paper should state.
- [ ] [joint] **Reserves remunerated at the national deposit rate (choice).** A union central bank
      paying a country-specific rate; under national clearing r_D and r_F diverge in a crisis,
      so a union rate would give reserves an excess return.
- [ ] [joint] **Sign of the funding-risk effect (interpretation).** The paper says sovereign risk
      tightens the constraint through E[Omega] (the bank's expected marginal value of wealth);
      the algebra says it relaxes it through the franchise value. Reverse the sentence and
      show the three margins (franchise value, leverage, deposit rate).
- [ ] [joint] **The premium proposition is unproven (interpretation).** The appendix proof is about a
      rise in priced risk, not realised default, cites undefined conditions, and assumes its key
      step. Its premise fails in the solution: consumption is 2.2% higher in default, so
      Omega(1) > Omega(0) (ratio 1.29) comes from the franchise value alone. The "risk premium"
      is a premium for leverage risk on bank balance sheets, not a precautionary premium on
      household consumption.
- [ ] [joint] **The credit spread is a near-cancellation (numerical).** It is the difference of two
      terms near 0.127 a quarter, so the accuracy of E[Omega] is first order for every priced
      number.

**Firms**
- [ ] [joint] **Working-capital wedge (accounting error, priority 1).** The labour condition charges
      w(1 + zeta r_wc), but the loan is inter-period; the consistent wedge is
      zeta (r_wc - r)/(1 + r), with zeta the pre-financed share of the wage bill. Correcting it
      removes the deposit-rate cushion: no-backstop impact output -0.069% to -0.174%, and the
      backstop's impact effect changes sign. Dividends must then carry producers' actual
      surplus, or the goods market leaks about 0.17% of output.

**Default**
- [ ] [joint] **Default branch off its solved domain (accounting error).** After default, D-bank
      holdings are 0.238 against a box floor of 0.690, so every value that prices expected loss
      is extrapolated. Fix with a stock-and-share state space and a friction normalised by the
      current stock, then choose the post-default fiscal rule.

**The backstop and the central bank**
- [ ] [joint] **The backstop economy is not closed as written (omission).** D-bond clearing must
      include Eurosystem holdings (b^D_D + b^D_F + M = B_D); the state must carry the book M and
      the reserve obligation O; the purchase rule needs its complementarity conditions.
- [ ] [joint] **The backstop insures bank balance sheets against default; the text says it does not
      (inconsistency).** Pari passu with private holders, the Eurosystem bears haircuts, which is
      the only channel through which a purchase supports the price. Reporting conditional on no
      default shows the insurer's premium income without its claims.
- [ ] [joint] **Remittance rule and capital-account closure (choices that drive the output results).**
      Profit is remitted mark-to-market each quarter, whereas the Eurosystem remits realised
      income on historical cost. The capital-key split ignores the November 2012 Eurogroup
      commitment to return SMP income to Greece. The capital account is nearly shut by the
      untargeted psi (cross-border holding cost). Together these generate the transfer problem.
- [ ] [joint] **Government budget gives each treasury the whole profit (error); kappa (the capital
      key) is never stated in the paper** (code 0.071; size share 0.081).
- [ ] [joint] **The floor at the corner** is met by pushing F holdings below their anchor, so required
      purchases depend on the untargeted psi.
- [ ] [joint] **Unsupported claims:** flight to safety (nothing in the model produces it); under
      national clearing the total-spread cap now responds to r_D - r_F.

**Implementation discrepancies (small)**
- [ ] [joint] Cross-border bond conditions omit p_t / p_{t+1}; F taxes are held fixed (a budget leak of
      about 2e-4 of F output); gamma_tau (tax response to debt) is evaluated at the risk-free
      price (root 0.9294); the banker's discount factor omits the change in consumer prices.

## 3. The viva paper

- [ ] [joint] Fill every [TBC] once sections 1-2 are done: rest-point mu, rest-point and impact
      sovereign spread, bank-held debt / GDP, bank assets / GDP, deposits / GDP, the
      post-default tax cut, beta_X, chi_X, Z_X, T_X.
- [ ] [joint] Fill the \placeholder items: source for the idiosyncratic income process (0.9, 0.2);
      HFCS distributional moments; sensitivity of output to zeta; the reason for the payout
      share f = 0.08; the year of the relative-size ratio 11.3.
- [ ] [Adam] Pull the local Overleaf clone before editing it; it lags Overleaf.

## 4. Global model, numerical

- [ ] [joint] `TPI_CAP = "rest"` (defend the no-backstop rest-point spread, about 56 bp) never run;
      expect the corner almost everywhere.

## 5. Sequence-space (SSJ) model

- [ ] [joint] Regenerate E1-E4 and the `experiments/paper/` figures; all predate the structural
      refactor (order in STATE section 6).
- [ ] [joint] Delete `psi_lambda_B` (collateral-friction amplification dial) and the risk-sensitive
      branch of `collateral_quality_D/F`; add the name to the banned-names scan in
      `Claude files/linear_ssj/tests/test_nkpc_blocks.py`.

## 6. Data and replication (detail in `DATA.md`)

- [ ] [Adam] **Recover the empirical code from the shelved branch.** The scripts behind the EBA
      exposures, foreign share, leverage 7.98, the 5 bp wedge, the CDS series and the lending
      rates exist only in commit `c2d2337` on `cross-border-empirics`. Restore
      `Empirics/eba_descriptives.py`, `Empirics/bocola_greece/`, `Empirics/cross_border/`, the
      four missing modules of `Empirics/greek_cds/`, `data/bocola2016/` and the Greek bank
      balance sheets.
- [ ] [Adam] Commit `Empirics/greek_fiscal/bohn.py` and `Empirics/greek_cds/risk_process.py`.
- [ ] [joint] Decide which derived outputs to track (`Empirics/outputs/` is git-ignored).
- [ ] [Adam] Update stale model constants: `risk_process.py` uses rho_s 0.95 and sigma_s 0.4455;
      `bohn.py` uses q_bar 0.91057 and debt root 0.93.
- [ ] [Daniel] Record the sources and calculations behind his calibration values that have no
      script: G/Y 0.204 / 0.191, coupon decay 0.056 / 0.039, import share 32%, relative size 11.3
      (DATA.md section 4).
- [ ] [joint] Script the validation moments computed in chat sessions and never saved (deposits
      84% of GDP, bond-implied hazards, bank assets / GDP and the rest of DATA.md section 4).
- [ ] [joint] Reconcile lending spreads 290 / 180 bp (stored data: 280 / 176), the 108 bp overdraft
      rise (124), the 0.092% rest-point hazard (0.102%) and the 4.4-year duration.

## 7. Documentation

- [ ] [joint] Reconcile `Claude files/docs/SPEC.md` with the paper (last revised 2026-10-01).
- [ ] [joint] The commit check (`.claude/hooks/require-docs-before-commit.sh`) inspects the code
      repository even when a commit runs inside the Overleaf clone, so `git commit -a` there is
      blocked by unrelated uncommitted code changes. It should use the repository the command
      runs in.
- [ ] [Adam] The handoff plugin (user-level Claude settings) writes `docs/handoff/HANDOFF.md` before
      compaction; disable it there or delete the file when it reappears.

## 8. Long term

- [ ] [joint] **Welfare of the backstop: level versus change rules (shelved 2026-10-27).** Holdings
      target theta = gamma_L (sp - s_bar) + gamma_D (sp - sp_lag): respond to the spread's
      level, its change, or both, with welfare in consumption equivalents from a Monte Carlo
      cleared exactly each quarter. Findings on the old calibration: the union gains, entirely
      from the announcement; D loses and F gains, because profits flow to F by the capital key;
      all rule classes lie on one frontier and no interior optimum was found. Open: robustness
      to the trade elasticity (0.5 in the code, 1.5 in the paper), why the level rule saturates
      when D keeps the profits, and the household overlay bias (section 2). Must be redone on
      the paper's model after sections 1-2. Code: local branch `shelved/policy-rules-welfare`;
      record: `Claude files/docs/shelved/welfare_and_policy_rules.md`; the paper's former
      section 6 (removed from Overleaf 2026-10-27): `Claude files/docs/shelved/viva_section6_welfare/`.
- [ ] [joint] **Warranted floor.** A cap that moves with fundamentals: cap(s) = EL(s) + tau, where EL(s)
      is the spread an unconstrained risk-neutral investor would charge for expected default
      loss at risk state s. The backstop then removes only the balance-sheet premium and never
      prices default risk below its expected loss, so it does not subsidise insolvency: the
      natural formalisation of a TPI that targets "unwarranted" spreads. Prototype (on the old
      model) in the shelved branch, commit `bc1a83a`.
- [ ] [joint] **Endogenous non-strategic default.** Make the default probability depend on fiscal
      fundamentals without sovereign optimisation, for example a stochastic fiscal limit as in
      Bi (2012) and Bi, Foerster and Traum (2025): default occurs when debt exceeds a random
      limit, so pi^d rises with debt. This lifts the current principle that pi^d is never a
      function of debt, and lets the backstop feed back on default risk through debt dynamics
      (and raises fiscal moral hazard). The SSJ model has an estimated Greek fiscal limit
      (`Empirics/fiscal_limit.py`).
- [ ] [joint] **Bayesian estimation.** Estimate the risk process and the bank block jointly from
      Greek spreads, CDS, bank leverage, output and investment instead of calibrating them. The
      constraint binds only occasionally, so the likelihood is nonlinear (particle filter on the
      global solution). Needs the data work in section 6; observables already assembled are
      listed in DATA.md section 8; a design and data-requirements note is on commit `c2d2337`
      (`diagnostics/bayesian_estimation/`).
