# Project specification

What the viva paper argues and the model it argues with. Written from the paper (Overleaf
project 698b4f88aeef1d0e1d08cc0c, commit `5639c71`, 27 Oct 2026): sections 1-4 and appendices
A-C. When the paper and this file disagree, the paper wins and this file is updated. The
previous specification (sequence-space era, last revised 2026-10-01) is kept verbatim in
`shelved/SPEC_2026-10-01.md`; most of its framing rested on a model without a risk premium and
on results that have since reversed.

Current numbers: `STATE.md`. Open work: `TODO.md` at the repository root. Notation and code
names: `GLOSSARY.md`. Data provenance: `DATA.md`.

## 1. Research question

Should a central bank promise to buy a government's debt when its borrowing costs rise? Part of
a sovereign spread compensates for expected default losses; the rest is the premium the banks
that hold the debt require because their net worth is most valuable in the states where the
sovereign defaults. A backstop can compress the second part, and because the same banks finance
firms, that relieves private credit and output. It cannot remove expected losses without a
fiscal transfer: that risk moves onto the common central-bank balance sheet and, through the
profit split, onto the taxpayers of the other member. The paper asks which risks a backstop
removes, which it only redistributes, and how the design of the rule decides between them. The
spread alone is not a sufficient statistic for intervention, because the same spread can mix
the two parts in any proportion.

## 2. Positioning

- **Sovereign-bank transmission:** the doom loop (Farhi and Tirole; Brunnermeier et al.) and
  Bocola (2016), whose banker's problem and closed-form multiplier the model inherits, with the
  bank balance sheet re-anchored to Greek data. Empirical counterparts (Gennaioli et al.,
  Acharya et al., Altavilla et al.) are single-country and take sovereign risk as given.
- **Structural sovereign pricing:** Bocola and Dovis (fundamental versus rollover risk), the
  Eaton-Gersovitz and Arellano tradition (endogenous repayment). These price the bond without
  an intermediary, so they cannot ask what a purchase does to the balance sheet that holds it.
- **Gertler-Karadi banks in a monetary union:** Bi, Foerster and Traum (2026), the closest in
  architecture (two countries, endogenous default, liquidity risk, representative households).
- **The paper's object:** one constrained intermediary per country holding the sovereign and
  the capital stock on one balance sheet, in both countries, so a purchase reaches activity
  through the selling bank's constraint and its incidence can be traced across both.
- **Not this paper's lane:** multiplicity and self-fulfilling crises, strategic default, the
  decision to bail out. Default risk is exogenous and there is no bad equilibrium to remove.

## 3. The model (paper section 2, appendix A)

- **Two countries,** D (calibrated to Greece) and F (Germany), in a monetary union with one
  central bank; F is ϰ = 11.3 times D's size. Each produces a distinct traded good; p_t is the
  relative price of the F good. Prices are flexible.
- **Households:** a continuum per country with idiosyncratic productivity, GHH preferences and
  real deposits as the only asset. As solved, the aggregates are representative-agent; whether
  the household distribution enters equilibrium is open (TODO).
- **Production:** intermediate producers pre-finance a share ζ_X of the wage bill with bank
  loans at r^wc, the only channel from bank balance sheets to labour demand on impact; capital
  producers with Jermann adjustment costs and Tobin's Q.
- **Banks (the core):** Gertler-Karadi intermediaries financing all productive capital, both
  sovereigns and working capital out of net worth and deposits. One divertible fraction λ_X for
  every asset; the incentive constraint gives a closed-form multiplier μ = max{·, 0} (kinked at
  zero). Pricing conditions for capital, the domestic bond, the foreign bond (with the change in
  the terms of trade and a holding cost ψ on deviations from steady-state holdings) and the
  working-capital loan, whose rate is the deposit rate plus the credit spread λμ / E[Ω].
- **Deposit market:** national clearing; each bank is funded by its own households.
- **Sovereign risk:** long bonds with geometric coupon decay δ^b; the priced default probability
  π^d_t is the logistic transform of an AR(1) risk factor s_t and never a function of debt.
  Default writes the bond down to the recovery ϱ_D (0.45); there is no exogenous default cost.
  Every expectation averages over the default regime and the next risk factor.
- **Premium decomposition:** the bond price splits into expected payoff, a risk premium from the
  covariance of the payoff with the bank's discount factor, and a liquidity discount from the
  binding constraint. Following the 14 Oct audit, the risk premium is a premium for leverage
  risk on bank balance sheets, not a precautionary premium on household consumption.
- **Government:** a lump-sum tax rule with annual Bohn coefficient 4γ_τ = 0.04; the F debt
  stock is fixed.
- **Solution:** global; Chebyshev collocation on a Smolyak grid, Gauss-Hermite quadrature over
  s and the two default regimes, the multiplier from its closed form at every point.

## 4. The policy experiment (paper section 3)

- **Baseline: a credible spread cap.** The central bank buys the least needed to keep the D
  spread at or below a cap in the no-default regime, nothing in default; held to maturity, never
  sold; the cap is a floor on the D bond price that moves with the safe yield.
- **Financing and profits:** purchases are paid with reserves held by D banks at the D deposit
  rate, ranking pari passu; profit is remitted each period, κ^cb_D to D and the rest to F.
- **Two effects:** the announcement effect (the two economies' stochastic rest points, the
  backstop in place, nothing bought) and the intervention effect (their responses to the same
  rise in default risk).
- **Two departures, giving a 2×2 matrix:** an *unwarranted-spread* target (the cap applies to the
  yield over the default-warranted price, the value to a risk-neutral bank with the same funding
  costs, so the bank's balance-sheet premium is removed and the expected-loss part is not) and
  *uncertain execution* (the backstop is carried out with probability 𝔐𝔅 when triggered). The
  four cells map to a credible cap, OMT, the SMP and the TPI.
- Section 3 is marked work in progress, and neither departure is solved in code on `main`
  (TODO, Paper).

## 5. Calibration (paper section 4, appendix B)

Quarterly; D to Greece, F to Germany. Set directly from the literature or data, normalised in
the steady state, and three bank moments matched jointly (leverage 8, a 5 bp intermediation
wedge, own-sovereign holdings of 14.4% and 6.5% of bank assets). The risk process is estimated
on Greek ten-year spreads and one-year CDS (ρ_s = 0.9, σ_s = 0.41, π^d 0.10% a quarter at rest).
The crisis experiment raises π^d to 2% a quarter, the level Greek CDS priced between October
2009 and April 2010. Values and their status against the code: `STATE.md` section 1; sources:
`DATA.md`.

## 6. Framing rules

- The paper studies a **promise**: a conditional commitment priced before it is used, traced
  along two histories, announced and never used, and activated. It is not a pure
  announcement-effect paper, and never write that every effect comes from anticipation; the
  defended claim is that, with no bad equilibrium, the promise works through fundamentals
  (state-contingent bank wealth, hence credit, investment and capital).
- Never write that a backstop **eliminates** the spread: with exogenous default the
  expected-loss part cannot be removed without a transfer.
- Default risk is exogenous and **never a function of debt**; no balance-sheet friction is a
  function of the default probability; no spread equation takes the form spread += c·π.
- Describe **channels**, never "x% fundamental, y% non-fundamental".
- The Greek figures fix the **level of risk** only: Greece in 2010 would not have met the TPI's
  eligibility criteria.
- The fiscal rule is the same with and without the backstop, so the model abstracts from the
  fiscal moral hazard a backstop could create; say so rather than claim it away.
- Name the spread: intermediation wedge, credit (working-capital) spread and sovereign spread
  are different objects.

## 7. Out of scope (for the viva)

Strategic or endogenous default (post-viva: a stochastic fiscal limit), multiplicity and
sunspots, nominal rigidities, fiscal moral hazard, estimation (post-viva), and the welfare
evaluation of backstop rules (shelved 2026-10-27; post-viva).
