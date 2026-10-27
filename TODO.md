# TODO

The single list of open work. Ordered by what blocks the viva paper. Each item says what is
open and why it matters; the detail lives where linked. When an item is done, delete it here
and record it in `Claude files/docs/PROGRESS.md`. Last revised 2026-10-27.

Symbols: see `GLOSSARY.md`. Data provenance: see `DATA.md`. Current results: see
`Claude files/docs/STATE.md`.

## 1. Bring the global code to the paper's model (blocks every model number in the paper)

Every saved global result belongs to an economy the paper no longer describes (union deposit
clearing, leverage 5, a 100 bp wedge). Until this section is done the paper's model values stay
[TBC]. Comparison table: STATE section 1.

- [ ] **National deposit clearing.** Each bank funded only by its own households
      (dep_X = P^c_X HA_X), replacing the union market with real deposit-rate parity. Open
      questions from the 15 Oct review: the steady-state external balance is
      unspecified (p_bar is not 1 under asymmetric calibration; D banks' F-bond holdings never
      calibrated), and the household block now prices D deposits alone.
- [ ] **Paper calibration values** into `calibration/global_projection.py`: leverage 8, wedge
      5 bp/yr, own-sovereign 14.4% / 6.5% of bank assets, foreign share 12.7%, size 11.3,
      rho_s 0.9 (currently hard-coded 0.95 in the grid builder), sigma_s 0.41, sigma 1.35,
      Frisch 0.5, home bias 0.68, trade elasticity 1.5, epsilon 11, G 0.204 / 0.191,
      delta_b 0.056 / 0.039, annual Bohn coefficient 0.04.
- [ ] **Authors' decision: are reserves divertible?** Recommended non-divertible. Under the
      paper's current specification the closed form for mu (Proposition 1) fails once reserves
      are outstanding; non-divertible reserves make D banks unconstrained for about five years
      after a large purchase, so this decides whether purchases relax the constraint at all.
      `Claude files/diagnostics/model_audit/RESERVES_DEPOSITS_CB_REVIEW.md`.
- [ ] **Fix the accounting errors found in the 14 Oct audit** (the only record is the
      model-section-audit memory note; the reports were deleted):
  - working-capital wedge: the labour condition charges w(1 + zeta r_wc) but the loan is
    inter-period; the consistent wedge is zeta (r_wc - r)/(1 + r). Removing the deposit-rate
    cushion moves no-backstop impact output from -0.069% to -0.174% and flips the sign of the
    backstop's impact effect, so this changes the headline;
  - default branch evaluated off its solved domain (D-bank holdings 0.238 against a box floor
    of 0.690): every value that prices expected loss is extrapolated. Fix with a stock-and-share
    state space and then choose the post-default fiscal rule;
  - dividends must carry producers' actual surplus, or the goods market leaks about 0.17% of
    output;
  - central-bank block: appendix A's D-bond clearing omits Eurosystem holdings M; the government
    budget gives each treasury the whole profit; the bond conditions need complementarity
    pairs; the paper never states the capital key kappa (code 0.071; size share 0.081);
  - cross-border bond conditions omit p_t / p_{t+1} (small); F taxes fixed (budget leak).
- [ ] **Solve the paper calibration (needs an explicit go-ahead; about 4-5 hours with the
      backstop).** Then check the kink risk: report the fitted and exactly cleared output
      responses side by side; if mu is near zero at the rest point the output response is not
      numerically identified (STATE section 2).
- [ ] Re-solve the welfare economies and the section 5 policy economies on the new model.

## 2. The viva paper

- [ ] Fill every [TBC] once section 1 is done: rest-point mu, rest-point and impact sovereign
      spread, bank-held debt / GDP, bank assets / GDP, deposits / GDP, the post-default tax
      cut, beta_X, chi_X, Z_X, T_X.
- [ ] Fill the \placeholder items: source for the idiosyncratic income process (0.9, 0.2);
      HFCS distributional moments; sensitivity of output to zeta (wage-bill share pre-financed);
      the reason for the payout share f = 0.08; the year of the relative-size ratio 11.3.
- [ ] Welfare robustness to the trade elasticity. The incidence result rests on 0.5 (code)
      while the paper calibrates 1.5; the terms-of-trade externality scales with it.
- [ ] Decompose why the level rule's announcement saturates when D keeps the profits (measured,
      not decomposed).
- [ ] Do not describe the solved dynamics as heterogeneous-agent: the household block enters no
      aggregate condition, and the welfare overlay counts r_wc B_L as income without repayment
      (a 0.09-point bias, the order of the reported distributional effects).
- [ ] Pull the local Overleaf clone before editing it; it lags Overleaf.

## 3. Global model, numerical and housekeeping

- [ ] `TPI_CAP = "rest"` (defend the no-backstop rest-point spread, about 56 bp) never run;
      expect the corner almost everywhere.
- [ ] Section 5 policy economies: `results/POLICY/` holds solved economies for the base, the
      warranted floor at tau = 16 and 0 bp, caps 300 and 360 bp, the 200 bp cap with profits
      remitted to D (kappa_D = 1) and kappa_nfa 0.0025 and 0.04, all on the code calibration.
      No saved economy for the tightest non-triggering cap or the 9-node refinement.

## 4. Sequence-space (SSJ) model

- [ ] Regenerate E1-E4 and the `experiments/paper/` figures; all predate the structural
      refactor (order in STATE section 8).
- [ ] Delete `psi_lambda_B` (collateral-friction amplification dial) and the risk-sensitive
      branch of `collateral_quality_D/F`; add the name to the banned-names scan in
      `Claude files/linear_ssj/tests/test_nkpc_blocks.py`.

## 5. Data and replication (detail in `DATA.md`)

- [ ] **Recover the empirical code from the shelved branch.** The scripts behind the EBA
      exposures, foreign share, leverage 7.98, the 5 bp wedge, the CDS series and the lending
      rates exist only in commit `c2d2337` on `cross-border-empirics` (a WIP snapshot); the
      working tree holds only their compiled caches. Restore `Empirics/eba_descriptives.py`,
      `Empirics/bocola_greece/`, `Empirics/cross_border/`, the four missing modules of
      `Empirics/greek_cds/`, and `data/bocola2016/` and the Greek bank balance sheets.
- [ ] Commit `Empirics/greek_fiscal/bohn.py` and `Empirics/greek_cds/risk_process.py`
      (both untracked).
- [ ] Decide which derived outputs to track. `Empirics/outputs/` is git-ignored, so the files
      the paper's numbers come from exist only on this laptop.
- [ ] Update stale model constants in the empirical scripts: `risk_process.py` uses
      rho_s 0.95 and sigma_s 0.4455 for its model-side columns; `bohn.py` uses q_bar 0.91057
      and debt root 0.93, so its implied roots differ from appendix B's (which use 0.947).
- [ ] Script the numbers that were computed in chat sessions and never saved (DATA.md,
      section 4), starting with those the calibration targets.
- [ ] Reconcile the lending spreads: appendix B quotes 290 / 180 bp (2003-2007) and a 108 bp
      overdraft-spread rise; the stored ECB series give 280 / 176 bp and 124 bp on the obvious
      windows.
- [ ] Before any estimation, read `diagnostics/bayesian_estimation/DATA_REQUIREMENTS.md` on
      commit `c2d2337` (also not on the working branch).

## 6. Documentation

- [ ] Reconcile `Claude files/docs/SPEC.md` with the paper (last revised 2026-10-01: predates
      the calibration merge and national clearing).
- [ ] The handoff plugin (enabled in the user-level Claude settings) writes
      `docs/handoff/HANDOFF.md` before compaction. Disable it there, or delete that file when
      it reappears; this project keeps no handoff file.
