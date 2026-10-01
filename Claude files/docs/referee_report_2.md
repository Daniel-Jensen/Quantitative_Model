# Referee Report — Second Referee
## Focused on motivation, framing and contribution

Scope note: findings about the correspondence between Sections 2–3 and the
sequence-space code are parked pending the authors' reconciliation of the two
implementations. Everything below is a framing objection that survives that
reconciliation, because it bears on the paper's argument rather than on either
solver.

---

## Verdict

The introduction is well written and its argument is legible, which is why the
objections are structural rather than cosmetic. As it stands the motivation has
three load-bearing weaknesses: the model rules out the mechanism the policy
exists to address; the chosen laboratory is the one country the instrument
excludes; and the paper's second headline question is answered with the wrong
object. Each is fixable, and none is fixed by more results.

---

## M1. The model assumes away the reason the policy exists

Section 2 states that "We model default risk as *exogenous*", with a latent
factor `s_t` (sovereign-risk state) following an AR(1) and a logistic map to the
priced probability. Two paragraphs earlier the same section says the perpetuities
are "subject to *rollover risk*". Those are incompatible. Rollover risk is a
multiplicity concept — the bad equilibrium exists because lenders expect it to.
An exogenous AR(1) has one equilibrium and no self-fulfilling component at all.
The term must go, or the mechanism must arrive.

The framing problem is larger than the word. The introduction's own account of
OMT and TPI is a multiplicity account: instruments "designed to be effective
without being deployed", removing premia that "impaired the transmission and
singleness of monetary policy", never used and never needing to be. That is the
Calvo/Cole–Kehoe/Corsetti–Dedola story of a backstop that selects the good
equilibrium at zero cost. In a model with exogenous default risk there is no bad
equilibrium to select away from, so the instrument cannot do the thing the
introduction spends two paragraphs saying it was built to do. What the paper can
evaluate is a price-support operation acting on intermediary balance sheets. That
is a real and interesting object, but it is not what the introduction promises.

This also reframes the paper's most graceful passage. The introduction declines
the fundamental/non-fundamental decomposition, positioning that as a considered
methodological stance — "one equilibrium object seen twice, so no share of the
premium can be labelled non-fundamental" — and as complementarity with
Bocola and Dovis (2019). A hostile reading is simpler: the non-fundamental
component has been assumed to be zero, so of course no share of the premium can
be attributed to it. Declining a decomposition you have ruled out by assumption
is not a discipline. Either endogenise the default decision — a fiscal-limit
formulation is the obvious route and is standard in this literature — or state
plainly, in the introduction, that the paper studies the fundamental-risk
transmission channel conditional on an exogenous risk process, and that the
announcement/multiplicity channel is out of scope. The second is an acceptable
paper. The current framing tries to have both.

## M2. Greece is the one sovereign the instrument would have excluded

"Greece is the cleanest available laboratory" is the paper's central design
choice and it is defended in one paragraph, on the grounds that the episode was
large and long. Largeness is not cleanliness. Two objections, and the second is
serious enough that a referee will not let it pass.

**Confounding.** The Greek collapse ran alongside a troika-administered fiscal
consolidation, three official programmes, IMF conditionality, structural
conditionality across labour and product markets, an actual restructuring,
redenomination risk that was openly priced from 2011, and eventually capital
controls. The model contains none of these: no official-sector financing, no
programme, no austerity path distinct from the tax rule, no exit risk. Ascribing
a 25% output collapse to a bank–sovereign feedback in that setting requires an
argument about the other channels, and the paper does not make one. The
introduction currently uses the size of the Greek depression as evidence that the
mechanism is identifiable — "violently enough, to identify the mechanism" — which
inverts the inference. A large collapse with many simultaneous causes identifies
less, not more.

**Eligibility.** This is the harder problem. OMT was available only under an
appropriate EFSF/ESM programme and, for a country under a full macroeconomic
adjustment programme, only once it had regained bond-market access. Greece in
2010–12 had lost market access entirely, which is precisely why it was on official
financing. TPI's published criteria are stricter still: compliance with the EU
fiscal framework, absence of severe macroeconomic imbalances, fiscal
sustainability, and sound and sustainable macroeconomic policies. Greece in
2010–12 fails all four on any reading. The paper's title promises "an evaluation
of the Transmission Protection Instrument" and its research question is "what
would a TPI-style backstop have done had it existed during the Greek episode".
The answer a referee will supply is: nothing, because it would not have been
deployed. Italy 2011–12 and 2018, or Spain 2012, are the market-access,
plausibly-eligible cases, which is precisely why Bocola (2016), Bocola–Dovis and
Bi et al. all work on Italy.

The authors have three routes, and must take one explicitly in the introduction.
Reframe the object as a generic sovereign backstop and drop TPI from the title
and the question. Keep TPI and move the laboratory to Italy. Or keep both and
argue that the eligibility criteria are not binding for the counterfactual — but
that argument has to be made, at length, and it is not currently attempted.

## M3. The opening paragraph promises a loop; the model delivers a chain

Paragraph 1 sets out four legs: repricing reduces bank net worth, tightens the
leverage constraint, contracts credit, and deteriorates "the tax base that
justified the original revision". The last clause closes the circuit and is what
earns the phrase "doom loop", explicitly contrasted with "a metaphor for
contagion in general".

With exogenous default risk, that last leg does not exist. Output and the tax
base can fall as much as they like and the default probability does not move. The
paper's transmission is a one-way chain: sovereign risk → bank net worth →
credit → activity. That is a legitimate object and much of the empirical
literature the paper cites documents exactly that leg. But the opening paragraph
describes a feedback the model cannot produce, and the phrase "doom loop" is then
used throughout on that basis. Either close the loop or rewrite paragraph 1 to
describe the chain and say why the open version is the right object to quantify.

## M4. The motivating statistic is a model artefact, and it is off by roughly a factor of five

"Greek banks entered the crisis holding the overwhelming majority of their own
government's paper — in the calibration below, 87%." This sits in the paper's
first paragraph, in a sentence otherwise made of facts, and reads as one.

It is not. The 87% is an equilibrium share in a model whose only holders of
sovereign debt are two banking systems: there is no official sector, no domestic
non-bank sector, no foreign non-bank investor. It is 87% of the paper the model
issues, not 87% of Greek general government debt. The paper's own EBA input has
Greek banks holding roughly EUR 54bn of Greek sovereign exposure against a
general government debt stock of about EUR 330bn at end-2010 — on the order of
16%, and falling as official financing displaced private holders through 2011–12.
The empirical claim the sentence makes is false by roughly five times.

The underlying point — Greek banks were *concentrated*, holding an outsized
position relative to their own thin equity — is true, is the right point, and is
supportable directly from the EBA disclosures the paper already uses. Make that
claim, with that source, and keep model shares out of the motivating paragraph.
Note also that the asymmetry of the two banking systems, which the results lean
on heavily, is a calibrated feature the introduction never mentions.

## M5. "No free parameter appears anywhere in it" oversells

The claim attaches to the introduction's central identity, which pins the
sovereign premium as a fraction `Δ` (pledgeability of sovereign paper relative to
capital) of the capital premium. The identity has exactly one coefficient, and
that coefficient is a calibrated number. Saying no free parameter appears in an
equation whose only parameter is free is a rhetorical move, and it is the kind a
referee reports. What is true and worth claiming is narrower and better: the
sovereign premium is not an *additive* wedge chosen to hit a spread moment, but a
proportionality whose single coefficient is a collateral primitive with an
independent interpretation. Say that, then defend `Δ` — which the introduction
never states a value for, and which the body never introduces at all.

The accompanying swipe — "a discipline earlier versions of this model, and much
of the applied literature, did not impose" — should be cut or named. Bocola (2016)
derives the sovereign premium from intermediary first-order conditions and is
cited approvingly two paragraphs later, so "much of the applied literature" is
carrying weight it has not earned.

## M6. "No direct lever on the spread at all" is an assumption presented as a restriction

The introduction bills this as the paper's contribution to what backstops can and
cannot do: purchases move the spread only by moving the return on capital, so the
spread effect and the investment effect are one effect. It is a nice result and I
believe it holds in the model. But it holds *because* a single representative
constrained intermediary is the marginal holder of every asset in the economy, and
because its problem is linear in portfolio shares under a constraint that always
binds. Relax any one of those — an occasionally slack constraint, a non-bank or
official marginal holder, segmentation, or any portfolio adjustment cost on the
domestic leg — and the proportionality breaks. The model itself already carries an
adjustment cost on the cross-border legs, which is an admission that the
proportionality cannot be imposed everywhere at once.

Presented as "a restriction on what a sovereign backstop can and cannot do", this
reads as a general theoretical claim. It is a property of a particular market
structure. State the market structure in the same sentence, and the result becomes
defensible rather than overreaching. The phrasing also risks being read as an
efficacy claim — "no direct lever on the spread at all" — when the paper's own
results have the backstop compressing the spread substantially. Direct versus
indirect is a decomposition point, not an efficacy point, and the sentence does
not currently distinguish them.

## M7. "Who would have paid for it?" is answered with exposure, not incidence

This is half of the stated research question, and the introduction answers it with
a capital-key share: 93% of the resulting exposure sits with the German treasury.
Exposure is not payment. On the paper's own baseline, default is priced but never
realised, so along the traced path the central bank buys distressed paper that
subsequently performs and the operation is *profitable*. On that path the German
taxpayer does not pay; the German taxpayer is paid. What a reader wants is the
expected fiscal transfer integrated over the default distribution the model
prices, together with its distribution across the two countries and — given the
paper's fourth advertised contribution — across households within them.

"It is a transfer, it is signed and scaled" is asserted three times in the
introduction and never given a number of the right kind. Either supply the
state-contingent expected cost, or retire the second half of the research question
and the "who would have paid for it" framing with it.

## M8. Heterogeneity is advertised, and then not used

The fourth contribution claim is that embedding this in a heterogeneous-agent
economy "makes the distributional incidence of an unconventional policy
computable". The three findings that follow are: aggregate impulse responses,
bank sovereign concentration, and who sells to the central bank under the capital
key. Not one is distributional. A reader finishes the introduction unable to say
why the model needed heterogeneous households, and the honest answer for these
three results is that it did not — every one of them would survive a
representative-agent household block.

This is the most easily fixed objection in the report, because the incidence
results exist: there is a quintile incidence table and two distributional figures
in the project's own results set, and none of them reach the introduction. A
two-country heterogeneous-agent model with Gertler–Karadi intermediaries and
defaultable perpetuities is an expensive object, and the introduction has to
justify the expense with a result only that object can produce. Lead with the
incidence finding, or demote the claim.

## M9. The marginal contribution over the closest analogue is thin as stated

Against Bi et al. — named as "the nearest published analogue" — the introduction
claims three differences: heterogeneous households, a premium from portfolio
optimality rather than a calibrated wedge, and an output contraction from the
intermediary constraint rather than a loan-in-advance requirement. The third is
contradicted by the paper's own Section 2, which installs a working-capital
requirement and states in bold that it is the only impact channel from spreads to
output. The second is at least arguably true of the antecedent literature as
well. That leaves the first, which per M8 produces none of the reported findings.

The contribution is almost certainly real — the proportionality restriction of M6
and the selling-side asymmetry of the results are both genuinely novel. But they
are not what the paragraph claims. Rewrite it around what the paper actually does
that its antecedents do not.

## M10. "Unanswered by the data" overclaims against literature the paper itself cites

"Both were designed to be effective without being deployed, and both have been.
That leaves the central quantitative question unanswered by the data." Two
paragraphs later the paper cites the announcement-window evaluations of exactly
these programmes. The OMT announcement is one of the most heavily measured policy
events in modern macroeconomics, and those measurements are data about precisely
this instrument. What the announcement studies cannot deliver is the general
equilibrium counterfactual and the distributional incidence — which is the paper's
real claim, and which it makes correctly later in the same paragraph.

The stronger move is to stop positioning the paper as filling a vacuum and start
using those estimates as external validation. A model of a never-used instrument
has no natural target; the OMT announcement effects are the closest thing
available, and matching them would materially raise the credibility of the
counterfactual. At present the introduction advertises that none of the reported
numbers is tuned to a target, which is honest, but leaves the reader with no
independent check on any of them.

---

## Smaller framing points

1. The section titled "Empirical Analysis" contains a calibration discussion and a
   solution method, and no empirical analysis. Retitle.
2. There is no conclusion, and the policy section is a title and a label. The paper
   currently cannot be assessed as an evaluation of anything.
3. "Evaluation" in the title implies welfare. If welfare is secondary in this
   project — and the project's own notes say it is — the title overclaims.
4. The abstract is `Lorem Ipsum`. The abstract is where M1, M2 and M8 will be
   judged, so it should be written after they are resolved, not before.
5. The road map promises households, production, intermediaries, the default event
   and the government. It does not mention the central bank, which is the paper's
   subject.
6. The paper says the recovery rate is the 30% "implied by the 2012 restructuring",
   while the body cites a 53.5% face-value cut and a 59–65% net-present-value
   haircut. 30% recovery is a 70% haircut, outside the range cited as its own
   authority. Fix the number or drop the appeal.

---

## What would make this reviewable

In order of importance:

1. Choose, in the introduction, between the multiplicity framing and the exogenous
   risk process, and make the paper consistent with the choice (M1, M3).
2. Defend Greece against the eligibility objection, or move the laboratory, or
   drop TPI from the title and question (M2).
3. Replace the 87% with a sourced concentration statistic from the EBA
   disclosures (M4).
4. Lead the contribution paragraph with a distributional result, using the
   incidence output that already exists (M8, M9).
5. Answer "who would have paid" with an expected fiscal transfer integrated over
   the priced default distribution, or retire the question (M7).
6. Qualify the proportionality result by the market structure that produces it,
   and state and defend the pledgeability coefficient (M5, M6).
