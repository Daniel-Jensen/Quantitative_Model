# Referee Report — Fifth Referee

## On the introduction, read against a full audit of Section 2

**Scope.** `VIVA/sections/01-introduction.tex` as of 2026-08-28, reviewed from scratch.
This report differs from `docs/referee_report_3.md` in one respect that changes most of
the findings: Section 2 has since been audited equation by equation
(`docs/referee_report_4.md`), so every claim the introduction makes about the model can
now be checked against what the model actually contains rather than against what
Section 2 says about itself.

The draft has moved twice since report 3: that report's fixes were applied, and five
Overleaf web-editor commits on 2026-08-27 trimmed the ECB, OMT-evidence and
exogenous-default paragraphs. Lines 100–117 were rewritten by me at the author's
instruction on the same day; they are reviewed here on the same terms as everything
else.

Line numbers refer to the current file.

---

## Verdict

**Major revision, and the problems are no longer in the prose.**

Report 3 found an introduction that promised things Section 3 did not deliver. Those
promises have been repaired. What the Section 2 audit reveals is worse and more
specific: **the introduction's three stated research questions map one-to-one onto the
three places Section 2 is incomplete.** Not approximately — exactly.

| Question (ll. 151–154) | What answering it requires | Status in Section 2 |
|---|---|---|
| "what the backstop does to the equilibrium" | an equilibrium | never defined; no market clearing of any kind |
| "which of the two banking systems supplies the bonds" | the cross-border Euler conditions | omitted at l. 290; `b^F_D` determined by nothing |
| "how the exposure and the gains fall … across households" | a distribution that moves aggregates | no distribution in the state vector; aggregate SDF asserted |

An introduction is allowed to promise what a later section delivers. It is not allowed
to promise three things the model is structurally unable to produce. That is **I1**,
and it subsumes report 3's C7.

Underneath it sits a quantitative exposure the paper has not confronted (**I2**), a set
of motivating statistics that are not model objects (**I3**), and one inconsistency
created by yesterday's edits (**I4**).

Three findings run the other way, and are recorded as such: on the risk-premium
mechanism, on rollover risk, and on the status of the recessionary-default claim, **the
introduction is right and Section 2 is wrong**. Fix Section 2 to match the
introduction, not the reverse (**I9**).

---

# Findings

## I1. All three stated questions are currently unanswerable

Lines 150–154:

> \Cref{sec:policy} adds an unlimited, country-specific backstop of the TPI kind to
> this environment. It asks what the backstop does to the equilibrium, which of the two
> banking systems supplies the bonds bought under it, and how the resulting exposure
> and the resulting gains fall across the two countries and across households.

Taking them in order.

**"What the backstop does to the equilibrium."** Section 2 defines no equilibrium. It
contains no goods-market clearing, no bond clearing tying the banks' holdings to the
governments' outstanding stocks, no deposit clearing, no labour clearing, and no
definition of a recursive competitive equilibrium. A purchase programme is an
intervention *in* market clearing; there is nothing for it to intervene in.

**"Which of the two banking systems supplies the bonds."** This is a portfolio
question, and it is answered by the cross-border Euler conditions. Section 2 imposes
its asset-pricing condition for every asset class, then writes out only the own-leg
conditions and silently drops both cross-border ones. The consequence is that
`b^F_D` (D-bank holdings of F government paper) and its mirror are determined by
nothing at all. The paper's second question is answered by the one equation Section 2
omits.

**"How the exposure and the gains fall … across households."** This needs the wealth
distribution to affect aggregates. Section 2 never aggregates the household block; the
banker's stochastic discount factor is written on the *aggregate* consumption–labour
composite, which with a binding borrowing constraint is nobody's marginal rate of
substitution; and Section 3's seven-dimensional state vector contains no distribution.
As the model is solved, the aggregate block is representative-agent.

The third question also carries a design problem the paper has not noticed. Taxes are
lump-sum in the household budget, and the fiscal rule responds to the surviving debt
stock. A central bank that changes the *ownership* of `D` paper without changing its
quantity therefore has almost no fiscal consequence for households under the current
rule. Whatever household incidence the paper eventually reports will be a property of
that financing assumption more than of the instrument.

Report 3 recorded Section 4 as "still a heading" and the results paragraph as dropped
at the author's instruction. That decision stands and is not re-litigated. The point
here is different: **the machinery required to answer the questions does not exist
yet**, so the gap is not a writing gap that Section 4 will close.

## I2. The introduction opens on 22 points of output; the model's only channel is worth fractions of a point

Lines 18–23 are the paper's motivating claim, and they are the right motivation for the
model's design:

> By the end of 2011 Greek output had fallen 22 percentage points below its 2007Q4
> level and investment 60 points. … The contraction ran ahead of the event that was
> supposed to cause it.

Lines 123–126 then state the transmission channel, correctly and exclusively:

> That common excess return is the model's credit spread, and it reaches output through
> a single equation: firms pre-finance a fraction of the wage bill, so the spread enters
> labour demand directly.

Section 2 makes the same claim more strongly, calling labour demand "the **only**
channel through which financial spreads reach output on impact". Combining that
equation with Section 3's GHH labour supply and the production function at
predetermined capital gives a closed form for the impact response, and it is small.
At conventional parameters (`\alpha` = 0.35, Frisch elasticity `\nu` = 1, pre-financed
share `\zeta` = 1):

| Rise in the credit spread | Impact response of output |
|---|---|
| +100 bp p.a. | −0.120% |
| +1000 bp p.a. | −1.204% |
| +2500 bp p.a. | −3.009% |

Reproducing the 22-point fall through this channel requires a credit spread of roughly
**18,300 bp per annum**.

The introduction does hedge, at l. 165–168: the Greek depression "had causes the model
omits" and "its aggregate path lies outside what the model is built to reproduce". That
hedge is honest and it should stay. But it is 145 lines after the claim it qualifies,
and it says the paper does not *target* the level. It does not say that the mechanism
is two orders of magnitude away from it. Those are different admissions, and a referee
who does the arithmetic will feel the difference.

The paper should state the elasticity, in Section 2 where the exclusivity claim is
made, and decide whether the motivating fact belongs at the top of the introduction on
those terms. The honest framing is available and is arguably stronger: the model is
built to isolate and price one channel, not to account for the Greek depression.

## I3. The motivating statistics are not model objects, and the paper later says so

Three facts carry the introduction's motivation. None of them is something the model
matches, and in one case the paper explicitly disclaims the correspondence.

**The debt ratio (l. 12).** "Debt rose from 105 per cent of GDP at the end of 2007 to
175 by the end of 2011." The steady-state debt stock `\bar B_X` is uncalibrated —
it appears in no table and no paragraph, Section 3's Government paragraph gives only the
coupon decay `\delta_b`, and `B_F` (F debt stock) is not even a state variable.
More fundamentally, the model's shock is a draw of the innovation to the latent risk
factor `s_t`; the debt stock plays no part in generating it. **The model cannot
represent the event the paper opens with** — a fiscal deterioration causing a
repricing. It represents the repricing alone.

**The EBA concentration (ll. 25–28).** "Greek banks … held €54.4 billion of their own
government's paper against €22.8 billion of Core Tier 1 capital." This is the
introduction's central piece of balance-sheet evidence. Then l. 163–165 says: "no
bank-level disclosure identifies the leverage of an intermediary defined to hold an
entire capital stock; the consolidated portfolio here has no counterpart in observed
balance sheets." Both statements are defensible in isolation. Together they say that
the fact motivating the mechanism is measured on an object the model does not contain.
Nothing in the calibration targets the sovereign-to-capital ratio on the model bank's
balance sheet, so the concentration the EBA number establishes has no counterpart in
the calibrated model.

**The 50 bp pre-crisis spread (l. 8).** This one is a missed opportunity rather than an
error. The model has a computable steady-state sovereign spread: the constraint
component is exactly 8 bp p.a. at the reported calibration, and the expected-loss
component implied by a 0.1% quarterly default probability and a Greek recovery rate is
roughly 21–26 bp p.a., for a total near **30–34 bp p.a.** against the 50 bp the paper
opens with. That is a respectable untargeted match and the paper does not claim it,
because Section 2 never writes an expression for the sovereign spread and Section 3
conflates the total with its 8 bp constraint component.

## I4. The paper now states the ECB's rationale and no longer says the model cannot represent it

Lines 62–64 attribute a purpose to the instruments:

> … exist to remove sovereign risk premia the Bank judges unwarranted by fundamentals.

Until yesterday, ll. 173–176 closed the loop:

> Premia unwarranted by fundamentals cannot arise here, so the model does not speak to
> that part of the ECB's stated rationale for these instruments. Endogenising the
> default decision is left for later work.

Commit `c302a09` deleted both sentences. What remains is "Default risk follows an
exogenous latent factor. The government's repayment decision is not modelled, so there
is no self-fulfilling equilibrium. Amplification comes only from the intermediary
balance sheet, which keeps the channel identified."

The reader is now told what the ECB says the instrument is for, and is not told that the
model cannot evaluate it on those terms. Every premium in this model is warranted by
fundamentals by construction. A paper proposing to evaluate an instrument whose stated
purpose is removing unwarranted premia must say, near the claim, that it does not model
them. The deletion also removed the only signal that endogenous default is a planned
extension.

This is the one finding in this report created by a recent edit rather than surviving
from an earlier draft, and it is a one-sentence repair.

## I5. "Occasionally binding" is not exercised by the calibration

Lines 51–54 justify the solution method:

> Sovereign risk inside a currency union is episodic … We therefore solve the model
> globally. The mechanism is an occasionally binding constraint interacting with a rare
> event, and a local approximation preserves neither.

The constraint is not occasionally binding on the reported path. Section 3 calibrates
the steady-state multiplier to `\bar\mu` = 0.0010, which is strictly positive: the
constraint binds at the rest point. The introduction says so itself at l. 139–141 —
"The steady state is calibrated to a barely-binding constraint, so the premium is
already active at the point the experiment starts from." The experiment then raises
risk, which reduces net worth and tightens the constraint further. The `max{·,0}`
operator in the closed-form multiplier never activates on the path the paper reports.

The global solution is well justified — by the default fork, by the branch-contingent
kernel, and by the curvature of the logistic default probability. It is not justified
by the occasional binding of a constraint that always binds. Say what actually requires
it.

Relatedly, l. 138–139 claims the risk premium "survives when the constraint is slack".
True of the model, and never visited by the experiment.

## I6. The contribution claim over the literature is unearned in both halves

Lines 94–98:

> All three work with a representative household, and none holds the sovereign and the
> productive capital stock on the same constrained balance sheet in both countries at
> once. Putting them on one balance sheet makes the incidence of a sovereign backstop
> computable: across the two countries, across the two banking systems, and across
> households.

Two claims of novelty, each currently undelivered.

**"In both countries at once."** The two balance sheets are linked only by the
cross-border sovereign positions, and those positions are pinned down by nothing
(**I1**). The structure is on the page; the equation that makes it operative is not.

**"All three work with a representative household."** As solved, so does this model
(**I1**). The comparison is only available once the household block does something to
aggregates.

One further item I could not verify and flag rather than assert: the claim that all
three cited papers use a representative household should be checked against
`bi2026asset` specifically before it survives another draft. A referee who works on
that paper will check it, and a wrong characterisation of the closest published
analogue is expensive.

## I7. "The safe sovereign does not covary" rules out a result the paper probably wants

Lines 135–136:

> Domestic sovereign bonds covary with that valuation; the safe sovereign does not.

The intended point is right: only the `D` sovereign carries default risk, so only its
payoff has the survival factor. But "does not covary" is much stronger than "does not
default", and it is false in this model. The `F` bond's payoff includes its own
continuation price, which moves with the `F` intermediary's constraint and its
stochastic discount factor. A `D` default hits the `F` bank, which holds `D` paper,
tightening `F`'s constraint and moving the `F` bond price. From the `D` bank's
perspective the cross-border position additionally carries terms-of-trade risk.

This matters beyond precision. A flight-to-quality result — safe yields falling when
periphery risk rises — is exactly a covariance effect on the safe bond, and it is one of
the more interesting things a two-country model of this kind can produce. The sentence
as written rules it out a priori. Weaken it to the default-risk point.

## I8. "The bond price separates into" three pieces, and it does not

Line 136–137 restates Section 2's three-way decomposition. The underlying algebra is
correct but the result is not a decomposition: the bond price appears on both sides,
and the constraint term is the price itself scaled by a factor. Solving properly puts
the constraint in the denominator, which is the pricing equation rewritten. The two
terms labelled "expected payoff" and "risk premium" are therefore not the corresponding
components of the actual price.

At the calibrated steady state the constraint term is about 2 bp of the price, so a
first-order treatment is numerically harmless and the economics of the sentence
survives. The word "separates" does not.

## I9. Three places where the introduction is right and Section 2 is wrong

Recorded so the repairs go in the correct direction.

**The risk-premium mechanism.** Section 2 states that because the banker values wealth
more in the default branch, the expected augmented discount factor *falls* when risk
rises. It rises: raising the priced default probability shifts weight onto the branch
just declared to have the higher value. The introduction (ll. 132–135) states only the
correct half and draws no sign conclusion. **The introduction is clean here; Section 2
needs the repair.**

**Rollover risk.** Section 2 calls the debt "subject to rollover risk" two lines before
declaring default exogenous. The introduction (ll. 173–174) states the exogeneity
plainly and never uses the phrase. **Delete it from Section 2.**

**Recessionary default.** Section 2 states this as a titled Proposition with an empty
body and a proof deferred to an appendix that does not exist. The introduction
(ll. 145–148) states it as a property of the solved decision rules, which is the honest
version. **Section 2 should adopt the introduction's framing**, not the other way round.

One caveat on the introduction's version: no decision rules are reported anywhere in the
paper, so the claim is currently unverifiable by the reader. That is downstream of
Section 4.

## I10. Not one number in the introduction comes from the model

Every quantity in eighteen paragraphs is a data moment: 50 bp, 105 and 175 per cent, 29
per cent, 59–65 per cent, 22 and 60 points, €54.4bn and €22.8bn. There is no
steady-state spread, no impulse magnitude, no elasticity, no loading, no welfare
number. For a quantitative paper this is unusual enough that a referee will read it as a
signal about the state of the results.

Some of this is the deliberate decision to drop the results paragraph until Section 4
exists, and that decision is not reopened here. But several numbers are available
*today* from Sections 2 and 3 and would cost nothing: the steady-state sovereign spread
(**I3**), the size of the impulse, and the transmission elasticity (**I2**).

On the impulse: Section 3 reports that the experiment raises the priced quarterly
default probability from 0.10% to 1.98% on impact. Inverting the logistic, that is a
move of 3.005 in the latent factor against an innovation standard deviation of 0.63 —
a **4.8-standard-deviation draw**. It sits inside the solution grid, so the numerics are
sound, but the introduction describes the experiment only as the probability being
"positive and rising". The magnitude should be owned in the introduction, in
standard-deviation units, because a reader who finds it in Section 3 will wonder why it
was not mentioned.

Separately: `fig:greece`, a full-page four-panel figure, is **never referenced in the
text**. `fig:decoupling` is referenced twice. Whatever the figure is for, the argument
does not currently reach for it.

## I11. Smaller items

- **"Disciplined throughout by the Greek episode" (l. 157).** Contradicted two sentences
  later, where leverage, payout rate, intermediation wedge and risk process are
  attributed to an estimated model of Italy. "Throughout" is the overclaiming the
  author has been stripping elsewhere; the Greek episode disciplines two parameters.
- **The recovery rate (ll. 158–160).** The introduction leads with it as the parameter
  Greece uniquely identifies. Section 3 says only "We accordingly set the recovery rate
  to match that event" and never gives a number. The one Greek-identified parameter the
  paper advertises is unreported.
- **Duration.** Section 2 gives the perpetuity's duration as the reciprocal of the
  coupon decay (17.9 quarters); Section 3 uses a discounted formula (16.9 quarters).
  Since the coupon decay is calibrated by inverting a measured Greek duration, the two
  formulas imply different parameter values. Not a typo, an identification question.
- **TPI eligibility (l. 168).** "Greece in 2010–12 would also have failed TPI's
  eligibility criteria" is the obvious attack on the whole exercise, and it is
  currently a subordinate clause in the middle of a long paragraph. It deserves to be
  met head on, because a referee will not let it pass as an aside.
- **The road map (ll. 178–180)** stops at Section 3 while the document contains a
  Section 4 heading and no conclusion. Deliberate per commit `17ded8a`, but the
  asymmetry will read as an oversight to anyone who has not seen that decision.

---

# On lines 100–117, which I wrote

Reviewed on the same terms.

The paragraphs no longer track Section 2's wording and no longer assert the
union-wide funding market and real interest parity that Section 2 never writes. Those
were the defects they were written to fix, and they are fixed.

Two criticisms stand against them.

**"Those claims sit on a single constrained intermediary in each country" (l. 107)**
opens on the cross-border positions and builds the paragraph on them. Per **I1** those
positions are determined by nothing. The sentence is not false — the claims are on the
balance sheet — but it gives structural prominence to the least well-defined object in
the model. If **I1** is resolved by fixing the cross-border position exogenously rather
than by adding an adjustment cost, this paragraph will need rewriting.

**"A two-country monetary union" (l. 100)** is doing work the model may not support.
By the paragraph's own account, what links the countries is goods trade and
cross-border bank holdings of sovereigns; there is no common policy rate, no common
funding market, no nominal instrument anywhere, and no monetary authority. That is a
two-country real model with financial integration. Calling it a monetary union is a
labelling choice, and the paper's policy question — the behaviour of a union-wide
central bank instrument — presumes a union with more structure than the model has. This
is the same gap as the absent central bank, seen from the introduction's side, and it
is worth confronting explicitly rather than resolving by vocabulary.

---

# Priority

1. **I4** — restore the sentence saying the model cannot speak to unwarranted premia.
   One sentence, and the paper currently overstates its own scope without it.
2. **I9** — repair Section 2 in the three places where the introduction is already
   right. No change to the introduction.
3. **I2** — decide how to frame the 22-point opening against a channel worth fractions
   of a point, and put the elasticity in Section 2.
4. **I5**, **I7**, **I8** — three sentences that claim more than the model delivers.
5. **I3** — either connect the motivating statistics to calibrated objects or stop
   presenting them as evidence about the model's bank.
6. **I1**, **I6** — not writing tasks. The introduction's questions and its contribution
   claim both become true when Section 2 is completed, and not before. Until then the
   introduction is writing cheques against work that has not been specified.

---

# Not reviewed

Per the author's instruction, proof-level items are excluded: notation and
sub/superscript consistency, citation-command choice, LaTeX mechanics, and the
`references.bib` metadata. Those from Section 2 are recorded in
`docs/referee_report_4.md`; the introduction has few.
