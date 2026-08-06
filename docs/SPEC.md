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

*Added 2026-08-06 with the `add-nkpc` workstream (sticky prices + nominal deposits). Full
numbers in `docs/STATE.md`.*

- **Rotemberg price Phillips curve, subsidy-neutralised.** `pi = beta*pi(+1) + kappa_p*(mu_p*mc − 1)`
  in both countries, with the markup wedge `mu_p*mc` entering labour demand
  (`w = mu_p*mc*(1−alpha)*Y/N`). Wages stay **flexible**. Writing the gap as the *ratio*
  `mu_p*mc − 1` rather than a level difference makes it unit-free: it linearises to exactly
  `mc_hat` for any `mu_p`, so published Calvo slopes map straight onto `kappa_p` with no
  steady-state rescaling, and `mu_p` is a free normalisation to first order. The steady state
  sets `mc_ss = 1/mu_p`, i.e. a production subsidy `tau_s = 1 − 1/mu_p` neutralising the
  markup, so `mu_p*mc = 1`, `profit_ss = 0`, `pi_ss = 0`, and **the entire steady state is
  bit-identical to the flexible-price model**. That is what makes the sticky-price results
  comparable to the earlier ones rather than confounded by a re-solved SS. `kappa_p → ∞`
  recovers flexible prices exactly, which is the standing equivalence gate.
  *Why it is needed at all:* under flexible prices, flexible labour supply plus competitive
  labour demand eliminate `Y` from the labour block entirely and pin `N` on `Z`, `K`, `P_CES`
  alone — there is nothing for aggregate demand to act on, and the crisis response is two
  orders of magnitude too small.
- **Union-inflation normalisation as the nominal anchor; no Taylor rule.** The monetary-union
  identity `p/p(-1) = (1+pi_F)/(1+pi_D)` pins the inflation *differential* off the existing
  unknown `p` (the nominal exchange rate is fixed at 1, so terms-of-trade movement **is** the
  inflation differential). The *level* is pinned by
  `omega_pi_D*pi_D + (1−omega_pi_D)*pi_F = 0` — the `phi_pi → ∞` limit of an ECB rule on
  union-wide PPI inflation. **There is deliberately no modelled policy rate**: no financial
  contract in this model carries one, so no Fisher relation is required to close the nominal
  side, and adding a Taylor rule would introduce a free parameter with no additional
  discipline. State it in the paper as an abstraction (perfectly credible union-inflation
  targeting), not as a modelled reaction function. `omega_pi_D = 0.071` is the renormalised
  two-country ECB capital key, **not** GDP weights: GDP weights would split any
  terms-of-trade move ~50/50 (because the model normalises `Y_D_ss ~ Y_F_ss ~ 1`) and erase
  the 93/7 Greek-deflation / German-inflation pattern that the 2010–12 internal devaluation
  actually took.
- **Nominal deposits against real sovereign bonds — a deliberate asymmetry.** Deposit
  contracts are nominal (`i_dep` is the contracted rate; `rdep_expost` carries the inflation
  surprise into `bank_return` and `capital_fund`), while sovereign bonds remain real. This is
  a choice, not an oversight: it makes banks nominal debtors and real creditors, which
  **maximises their Fisher exposure** and is the configuration under which the deflation
  channel does the most work. **It must be stated as such in the paper**, since the opposite
  convention (nominal sovereign debt) would give the sovereign an inflation-erosion channel
  and flip the sign of the bank's net Fisher position. Nominal sovereign bonds are a
  candidate extension, not a correction.
- **The markup rent is distributed in proportion to productivity `e`, not lump-sum.** Once
  labour is paid `mu_p*mc*(1−alpha)*Y` and capital keeps `alpha*Y`, the residual
  `(1 − mu_p*mc)(1−alpha)Y` must be routed somewhere or it is a Walras leak of the W-1/W-2
  class. Routing it proportional to `e` makes labour-plus-profit income per unit of `e`
  exactly `(1−alpha)*Y*e` — identical to the flexible-price model — so the markup wedge bites
  only on the *firm's hiring decision* and never on household income, and `labor_market_D/F`
  (labour supply) needs no change. A **lump-sum rebate was rejected because it is
  countercyclical**: markup rents rise when `mc` falls, so a lump-sum transfer would hand the
  largest windfall to the poorest households exactly in the downturn, manufacturing a
  progressive incidence result as an artifact of the rebate rule rather than of the
  transmission mechanism this paper is about.

## Calibration strategy

> **Values in this section date from 2026-07-22 and several are superseded.** The live
> calibration table is `docs/STATE.md`. In particular: `psi_lambda_B_D/F = 7.85` (not
> 1.1793 — re-tuned 2026-07-31 to 8.5 for `BANK_SCOPE="broad"`, then 2026-08-06 to 7.85
> once sticky prices and the Fisher channel pushed the spread response to 162bp);
> `EL_price_D/F = 0.056134` (not 0.0717 — that predates the EBA `delta_b=0.0777`,
> `q_b=0.969`); `delta_b_D/F = 0.0777/0.0568`, measured from the sovereign maturity
> ladder; `phi_lamb_D/F = 0.15`. The *reasoning* below is still the reasoning; the
> numbers are not all current. **Re-derive, do not copy.**

**As of 2026-07-22, see `docs/eba_calibration.md` for the full parameter →
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
   `psi_spread·def_rate` sits **on top**. **Current calibration (2026-08-06,
   sticky prices + nominal deposits, `psi_lambda_B=7.85`, `recovery_rate=0.30`):
   loading (TPI premium PV / expected-loss PV) is 3.82 at the medium regime and
   2.90 at the aggressive one, and above 1 at all 59 grid points of the schedule**
   — over-compensated, declining in aggressiveness. **The claim survived the move
   to sticky prices**; the flex-price values were 4.00/3.17. (An earlier number
   moved twice on 2026-07-22: first
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
   per unit. **Confirmed, and it survived the move to sticky prices**: at the
   current calibration (2026-08-06) the loading declines monotonically in gamma
   at all 59 finite grid points, **4.43 → 1.49 over γ ∈ [0.51, 30.00]**, staying
   above 1 throughout. The *decline* is the claim; the level is not.
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

**Report the decomposition, not the headline ΔY — because the channels cancel,
and they land on different households.** A consumption expansion, an investment
contraction and a net-export cushion offset each other; that offsetting is the
economics, and a RANK model cannot see it. Frame it as the reallocation it is,
never as "nothing happened."

*Restated 2026-08-06.* The earlier version of this caution said the headline was
*only* small because it was the residue of channels ~4× its size. **That is no
longer true and must not be repeated.** Under flexible prices `dY[0]` moved
+4.9e−04 passive → aggressive against an investment channel of +2.2e−03 and net
exports of −1.9e−03. Under sticky prices `dY[0]` moves **+1.38e−02** while
investment moves +3.41e−03 and net exports −2.92e−03 — the largest single channel
is now **0.25× the headline, not 4×**. The instruction is unchanged; its
justification has inverted. Arguably this is an improvement: output is no longer
a numerically fragile residue of nearly-cancelling terms, so the headline can be
quoted without the earlier caveat that it is an artifact of near-cancellation —
but the decomposition still carries the distributional content, which is the
reason to lead with it. Current numbers: `docs/experiments_results.md` (E2).

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
