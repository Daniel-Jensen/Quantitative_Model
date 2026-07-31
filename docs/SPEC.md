# Project Specification

## Goal

Build a tractable two-country general equilibrium model with heterogeneous households, financial intermediaries, sovereign debt, and cross-border portfolio frictions. Primary application: the 2010–2012 Greek sovereign debt crisis and ECB policy (TPI).

## Functional requirements

- Solve a steady-state equilibrium for two countries (D=Greece, F=Germany), each with:
  - Heterogeneous households optimizing deposits and labour supply (GHH preferences, EIS=0.5, Frisch=0.5)
  - Gertler-Karadi financial intermediaries holding domestic and foreign government bonds and productive capital
  - Production with capital accumulation, price of investment, and capital adjustment costs
  - Government fiscal budget with Bohn fiscal rule and Hatchondo-Martinez geometric-decay perpetuity bonds
- Endogenously determine:
  - Sovereign bond prices, yields, and spreads (default risk feedback via `def_scale`)
  - Deposit returns and banking net worth dynamics (GK incentive constraint and P1 Bellman)
  - Terms of trade and net exports (CES consumption basket over D-goods and F-goods)
- Cross-border portfolio adjustment costs (`psi_bF_D`, `psi_bD_F`) anchoring bilateral sovereign positions.
- Impulse response functions (SSJ Jacobian, T=100–500 periods) for:
  - Sovereign default probability shocks
  - TFP shocks
  - TPI (closed-loop CB bond purchase rule)

## Model code structure

- `code/equations_D.py` — country D household (EGM het block), bank, production, capital, government equations
- `code/equations_F.py` — country F analogues
- `code/equations_global.py` — global goods market, external account, bond clearing, portfolio adjustment costs
- `code/main.py` — production pipeline: calibration → SS solve → Jacobian → IRFs → TPI experiment (orchestrates `calibration.py`, `steady_state.py`, `full_model.py`, `tpi.py`)
- `routines/` — auxiliary: grid construction, income Markov chain, Gini calculation
- ~~`audit_artifacts/`~~ — removed 2026-07-30; findings retained in `docs/audit.md` and `docs/STATE.md`
- `plots/` — output figures (TPI welfare and spread panels)

## Research objectives

1. Study the interaction between sovereign risk and bank portfolio choice in a two-country monetary union.
2. Quantify how investment and deposit returns propagate through global goods markets (terms-of-trade channel).
3. Assess portfolio adjustment costs and cross-border sovereign exposure (bilateral Greece-Germany via EBA data).
4. Analyse the TPI (Transmission Protection Instrument): welfare and spread effects of a closed-loop CB bond purchase rule.
5. Evaluate the doom-loop mechanism (debt → spread → bank net worth → output → debt) and its dependence on bond duration, bank leverage, and fiscal rule strength.

## Key modelling choices (with rationale)

- **Y = F(K_t):** production uses current-period capital (not K(−1)). Capital producer receives `mpk·(K−K(−1))` to close capital income accounting (W-1 fix). Alternative: K(−1) timing eliminates the term; both are internally consistent.
- **Predetermined deposit rate:** `Rgross = (1+rdep(−1))·P(−1)/P`. Deposit contracts signed at t−1, so funding costs are predetermined — standard NK timing for bank liabilities.
- **Hatchondo-Martinez perpetuity:** bond coupon decays geometrically at rate `1−delta_b`. Duration = 1/delta_b quarters. Captures MTM capital losses on bank balance sheets.
- **GK agency problem:** divertable fraction `Delta` drives the IC constraint binding. Multi-asset IC requires separate `Delta` for each asset class.
- **Walras redundancy:** equations `ca_res_D` and `goods_mkt_F` are dropped from targets. Post-fix they hold to machine tolerance; see `docs/walras_forensics.md`.

## Calibration strategy

**Current (2026-07-22), see `docs/eba_calibration.md` for the full parameter →
moment map and `docs/STATE.md` for the live calibration table:**
- Bilateral GR/DE bank exposures from the EBA 2011 stress-test disclosure
  (31 Dec 2010 actual): own-book concentration `phi_bD_D_ss=2.39` (GR),
  `phi_bF_F_ss=2.76` (DE); cross-holdings `phi_bF_D_ss=0.018`, `phi_bD_F_ss=0.069`
  (q·sovereign-book / bank-capital, not asset-normalised — do not conflate with
  the pre-EBA "b_D_D/asset≈24.47%" moments this section previously cited).
- `psi_lambda_B_D/F = 1.1793` (re-tuned 2026-07-22 after resolving `EL_price`
  below), data-disciplined to the 2010 GR-DE spread (~150bp on a 1pp
  default-probability shock) — the amplification dial has no other empirical
  anchor and must be re-verified against this target after any structural
  change to the collateral/IC block, or to `recovery_rate`/`EL_price`
  (`audit_artifacts/psilam_moment_sweep_postC1.py`,
  `audit_artifacts/psilam_verify_postEL1.py`). **Values above ~1.5-2.0
  currently sit in a linear-approximation-breakdown region on this model and
  must not be used without re-checking stability at that value.**
- Bond duration: `delta_b_D/F=0.10` (2.5yr); empirical target `0.036/0.038`
  (Hatchondo-Martinez matching GR/DE 2011 avg maturities ~7yr/6.5yr) not yet
  ported to the committed calibration — see Finding F-1 in `docs/STATE.md`.
- Bohn fiscal coefficient: `phi_lamb_D/F=0.60`; literature 0.025–0.038 quarterly
  for EA periphery (Staehr 2008) is far below what this model needs for
  stability at current amplification.
- `EL_price_D/F≈0.0717`, from `(1-recovery)·delta_b/q_b`. **`recovery_rate_D/F=0.30`
  (resolved 2026-07-22)**, an NPV-recovery estimate for the actual March 2012
  Greek PSI: Zettelmeyer, Trebesch & Gulati ("The Greek Debt Restructuring: An
  Autopsy", PIIE WP13-8) find actual investor NPV losses of 59-65% (below the
  ~75% commonly quoted); contemporaneous bank estimates (Credit Suisse, Morgan
  Stanley) put NPV haircuts at 73-78%. 0.30 recovery (70% haircut) is central
  in this range. The pass-through ("consistency") moment — bank net worth
  response per 100bp of spread, ≈−4.5%/100bp — was checked against
  Acharya-Drechsler-Schnabl (2014 JF)'s bank-equity-return-on-sovereign-CDS
  elasticity and sits within its literature-implied range (−1.8% to
  −8.6%/100bp depending on baseline CDS level); see `docs/STATE.md` issue PT-1.

## Out of scope (current phase)

- Formal welfare analysis with non-linear transition paths (Jacobian-linearised only)
- External habit formation
- Macroprudential capital requirements (conceptual only, not implemented)
- Estimation / Bayesian identification

## Theoretical framing (paper narrative)

*Merged from `docs/FRAMING_HANDOFF.md` (2026-07-22), which is retired — this is
now the single source for the paper's argument. Numbers below were updated
where this session's C-1 fix / EBA calibration / `psi_lambda_B` recalibration
changed them; anything not explicitly flagged as re-verified should be treated
as needing a fresh check before being quoted in the paper.*

### The lane

Not a paper about multiplicity, self-fulfilling crises, strategic default, or
the *decision* to bail out. Those belong to Bocola-Dovis (origin/decomposition),
Gourinchas-Martin-Messer (the bailout decision), Fornaro-Grosse-Steffen (origin
of divergence via sunspots), and to a coauthor's parallel nonlinear global
model.

This paper: **transmission and incidence of a backstop of exogenously-varying
strength**, in the fundamental-risk regime the multiplicity literature
abstracts from. Backstop aggressiveness (`phi_TPI`, `gamma` in `code/tpi.py`)
is a policy stance, not an equilibrium object.

Positioning between the two papers that bracket it: Bi-Foerster-Traum (2025)
have the sovereign-bank nexus, two-country monetary union, asset purchases —
but representative households. Chiang-Zoch have GK-in-HANK with distributional
output — but no sovereign debt, closed economy. Nobody sits at the
intersection. **Caveat:** the fundamental-bailout space is crowded —
Gourinchas-Philippon-Vayanos (2017) is the closest macrofinance ancestor (Greek
crisis, heterogeneous households, fiscal/sudden-stop decomposition) but is
TANK, no GK block, no creditor country, no backstop-aggressiveness object. Cite
it explicitly; a referee will know it.

### Structural constraints that bind the theory

These are not caveats to bury — several framings die on them.

- **First-order certainty equivalence.** No risk aversion, no risk premium
  proper, no covariance term. The model prices the tail's *first moment*
  exactly and cannot price its second moment at all.
- **`def_rate` is exogenous.** No endogenous default, therefore no
  moral-hazard channel.
- **Always-binding IC.** No risk-shifting/gamble-for-resurrection;
  renationalisation cannot be endogenised (though its *quantity signature*
  emerges from market clearing — see "Model facts" in `docs/STATE.md`).
- **`writeoff_enabled = 0`, the risk-premium framing.** `def_rate` is a genuine
  probability; agents price expected loss; the IRF traces the **no-default
  branch**. This is standard (risk-premium-shock device), *not* "default is
  impossible."
- **Consequence — a live implementation hazard:** the excess-return flow
  `EL_price × def_rate_t × b_ss` is a first-order deviation times a
  *steady-state level*, hence **first-order and does not vanish**. Bondholders
  earn the premium along the computed path with no offsetting loss. That is an
  artifact of inserting a premium without the compensating branch. **CB
  expected loss must be hand-computed off-path** (`Σ β^t · EL_price ·
  def_rate_t · q_b · cb_buy_t`, `code/tpi.py`'s `cb_pnl`), never read off the
  realised path, which will mechanically show the CB profiting.

### The Region 1 / Region 2 structure (the organising device)

- **Region 1 (sunspot):** backstop eliminates the bad equilibrium, never
  executes. `def_rate → 0`, exposure → 0, tail → 0. Costless insurance.
  Draghi's claim, coherent *in this region*.
- **Region 2 (fundamental):** no bad equilibrium to eliminate. Intervention
  means holding real risk. Costless-insurance fails.
- **"Unwarranted, disorderly market dynamics" is a jurisdictional claim that
  the ECB is in Region 1.** One phrase does two jobs at once: satisfies
  Article 123 (transmission repair = monetary policy, not fiscal financing)
  *and* guarantees zero Haftungsrisiken, since Region-1 interventions are
  never executed. Load-bearing, not sloppy drafting.
- **The factual predicate is mostly false.** Bocola-Dovis: 12% non-fundamental
  for Italy; Greece almost certainly *more* fundamental (actual restructuring,
  ~10% primary deficit at peak, ~15% CA deficit, statistical misreporting =
  pure fundamental-information shock). Their OMT result is more direct still:
  spreads fell ~100bp *below* the no-rollover counterfactual — a measurement of
  the ECB acting against fundamentals.
- **The model lives entirely in Region 2 by construction.** That is a feature:
  it is the counterfactual world where the ECB's legal claim is false, which is
  exactly why it can price what the Court worried about. Import the 12% for the
  bridge: ~88% of intervention generates tail, ~12% is free.

### The German litigation mapping

The Bundesverfassungsgericht case was **never primarily "you're transferring
money."** It was **Haftungsrisiken / Haushaltsverantwortung** — Bundestag
budgetary autonomy compromised by liability it didn't authorise. A claim about
Region-2 *exposure*. PSPP proportionality is the same demand restated: show us
you weighed the fiscal side-effects, i.e. show us you're in Region 1.

Where a genuine "third thing" survives: the *Treaty's* categories are monetary
policy vs. fiscal financing, and the Court's actual concern — unauthorised
risk-bearing — maps onto neither. Which is why it had to be litigated through
German constitutional law rather than Article 123. That's an observation about
why the litigation was tortured, not a claim to have out-theorised the Court.

### Live claims (survived scrutiny)

1. **Expected P&L favours the CB.** `EL_price·def_rate` is actuarially fair by
   construction — the expected loss is *fully* compensated, not partially.
   `psi_spread·def_rate` sits **on top**. **Current calibration (2026-07-22,
   `psi_lambda_B=1.1793`, `recovery_rate=0.30`): loading (TPI premium PV /
   expected-loss PV) is 3.59/3.03/2.47 at gamma=2/5/10** — over-compensated,
   declining in aggressiveness. (This number moved twice the same day: first
   recalibrating `psi_lambda_B` to the 150bp target gave 2.54/2.14/1.74;
   resolving `recovery_rate` afterward — which shrinks `EL_price`, the
   denominator — raised it to the current 3.59/3.03/2.47. Both supersede the
   pre-fix ~7-7.6x figure at `psi_lambda_B=2.8-3.0`, which is no longer a valid
   calibration on this model — see `docs/eba_calibration.md`. Re-verify this
   number after any further recalibration.) **The monetary-financing objection
   fails on the model's own terms.**
2. **The `psi_spread` ambiguity — preserve it, do not resolve it.** Since the
   model has *no* risk-aversion channel, `psi_spread` is the only place a
   real-world risk premium could hide. Either (a) genuine agency rent
   (constraint shadow price, CB not entitled) or (b) the reduced-form stand-in
   for the risk premium the market demanded (CB fairly paid).
   **Indistinguishable from inside the model.** The *sign* survives either way.
3. **The right question is efficiency, not fairness.** Not "is the premium
   actuarially fair" but **"is the risk moving to a cheaper holder?"** The
   price is set by the marginal holder: a balance-sheet-constrained Greek bank
   at crisis peak, whose IC binds hard — which is *why* the wedge is a large
   share of the baseline spread. A diversified, unleveraged taxpayer base is
   plausibly a cheaper holder. Germany buys at the *constrained seller's*
   price. **This is the only version of "Germany doesn't lose" that survives
   risk aversion** — it's the standard risk-sharing-efficiency argument. The
   model gives you the seller's shadow price (`psi_spread`); it cannot give
   you Germany's valuation.
4. **Consent, not price.** TPI is a compulsory, fairly-priced-or-better
   insurance contract reallocating Greek sovereign tail risk from Greek balance
   sheets to euro-area taxpayers by capital key (Germany: capital-key conduit
   share `kappa_cb_F=0.929` of the two-country renormalised key; ~26.1% of the
   full euro-area key). The German objection is not that the price is unfair.
   **It is that the Bundestag never signed the contract.** The consent-vs-price
   distinction is the genuinely novel bit and the thing the economics
   literature cannot see.
5. **The profit self-extinguishes — a testable prediction.** The wedge exists
   *because* the marginal holder is constrained; TPI relieves the constraint.
   So intervention erodes its own profit source: more credible backstop →
   spreads compress toward fundamentals → `psi_spread` shrinks → less earned
   per unit. **Confirmed post-recalibration**: loading declines monotonically
   in gamma at the current calibration (3.59→3.03→2.47 at gamma=2/5/10).
   **"Germany profits" and "TPI works" are in tension.**
6. **The `EL_price`/`psi_spread` decomposition must not be confused with
   Bocola-Dovis's.** Ours is expected-loss vs collateral-friction; theirs is
   fundamental vs rollover (88%/12%). **Orthogonal decompositions.** State
   this explicitly or the numbers will be read as contradicting the solvency
   argument.

### Retracted — do not reuse (each sounds right and is wrong)

- **R1. "Unpriced tail insurance."** Wrong. The CB *collects* the premium, via
  the discounted `q_b`. It buys cheap precisely because the market demands
  compensation. The insurance is priced.
- **R2. "The tail isn't priced in our model."** Wrong and a referee will catch
  it. The first moment is priced exactly. What's absent is the
  *risk-aversion component*. Correct phrasing: "prices expected loss but has
  no risk premium proper."
- **R3. "ECB says costless / Court says financing / the truth is a third
  thing."** This strawmanned the Court by giving it a cruder position than the
  one it held, in order to have something to dissolve. Honest version: **the
  ECB's Region-1 claim was factually false and the Court's concern was
  well-founded.** Not "both wrong" — one side mischaracterised, the other was
  substantively right. The contribution is *quantifying* what the Court
  worried about, not adjudicating between two errors.
- **R4. "Fairly priced, therefore in Germany's economic interest."** Does not
  follow — it's contradicted by it. A fair premium leaves you indifferent only
  if risk-neutral. Forced fair insurance makes a risk-averse party *worse
  off*. And the model, being certainty-equivalent, **structurally cannot**
  adjudicate: it speaks only to the first moment when the entire dispute is
  about the second.
- **R5. "We prove Germany is better off bailing out Greece."** Three faults:
  "prove" (no CB balance sheet or capital key existed at time of writing —
  since resolved, see the capital-key conduit above); "better off" (see R4);
  "bailing out" (concedes the framing you just disproved — if expected P&L
  favours the CB it isn't a bailout).
- **R6. "Germany profits from buying distressed Greek debt."** Invites a
  vulture-fund reading. The CB isn't picking off panicked sellers; it earns
  the wedge because the marginal holder is *balance-sheet constrained* and
  pays rent to shed risk consuming its IC slack. Risk-sharing efficiency, not
  opportunism.
- **R7. "The 'unwarranted' wording is purely legal."** Too strong.
  Conditionality has a real economic function: it limits moral hazard. If
  markets expect fundamentals to be rescued, `def_rate` rises endogenously and
  the tail fattens. **Concede this** — it costs nothing (exogenous `def_rate`
  means we can't model it anyway; it's the coauthor's territory) and it's the
  difference between a sharp paper and advocacy.
- **R8. "The agency wedge is unrelated to the tail."** Both scale with
  `def_rate` — same driver, different channel. `psi_spread` is a *different
  pricing of the same underlying risk* by a constrained holder, not a windfall
  from an orthogonal source.
- **R9. "The wedge partially compensates the risk-aversion premium."**
  Presumes the answer. Whether the wedge over- or under-compensates a
  risk-averse taxpayer is **indeterminate** — it depends on a parameter the
  model doesn't contain. Preserve the indeterminacy; it's what makes the
  efficiency framing (Live Claim 3) the live question instead of the fairness
  one.

### Where Germany genuinely benefits: the trade channel

Greek stabilisation preserves German export demand, which survives
two-country bond-market clearing (the bond-portfolio spillback mostly nets
out: German banks fleeing *into* Bunds dominate Greek banks dumping them, so
the backstop, by calming the crisis, likely *raises* German yields). But GMM
already established creditor-interest-in-bailouts via portfolio
diversification, so the trade channel is *positioning*, not novelty. Cite them
and differentiate on channel.

**Do not lead with ΔY.** A small headline output number can be *only* small
because two large channels (investment contraction, NX cushion) are netting
out — and they land on different households. That's a strength (a RANK model
can't see it) only if framed as the reallocation it is, not as "nothing
happened." Check the current model's investment/NX decomposition before
asserting this — the specific magnitudes reported in earlier drafts predate
the EBA recalibration and need re-verification.

### The TL;DR as it currently stands

> The "unwarranted" condition does legal work, not economic work: on the
> model's terms intervention is fairly priced or better, so the
> monetary-financing objection fails on its own terms. But the German
> objection was never about price — it was about *consent* to risk-bearing.
> We quantify the tail Germany assumes as a function of `phi_TPI`. Whether
> that constitutes a cost depends on risk preferences our framework cannot
> price; whether Germany nets out ahead turns on the trade channel, not the
> bond P&L. And the CB's expected profit self-extinguishes as the policy
> succeeds.

Survives contact, but every specific number in it needs to be re-quoted from
`docs/STATE.md`'s current calibration rather than from memory of earlier
drafts — this model's numbers moved substantially in July 2026 (EBA
calibration, C-1 fix, `psi_lambda_B` recalibration) and will likely move again
before the paper is final.
