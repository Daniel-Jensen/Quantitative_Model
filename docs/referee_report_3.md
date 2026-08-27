# Referee Report — Third Referee

## On the introduction, read against Sections 2 and 3 of the Overleaf draft

**Scope.** The introduction (`VIVA/sections/01-introduction.tex`) checked for language,
internal coherence, and consistency with the compiled paper — `02-model.tex`,
`03-calibration.tex`, `04-policy experiment.tex`. No comparison to the repository's
sequence-space implementation: Section 2 describes the coauthors' global model and the
two are separate objects. No data verification.

The comparison draft throughout is `01-introduction.tex.bak`, which
`docs/referee_report_2.md` reviewed.

---

## Disposition (2026-08-27)

The introduction was rewritten against these findings and pushed to Overleaf. What
was actioned, and what was left:

**Actioned in `01-introduction.tex`.** C1 — the measurement justification is replaced
by what Section 3 does: periphery/core blocs, the sovereign side disciplined to the
Greek episode, the financial block openly attributed to Bocola. C2 — Bocola credited
in the literature paragraph and on the closed-form multiplier. C3 and C4 — the two
transmission properties reordered, the constraint channel labelled asset-neutral, the
risk premium identified as the sovereign-specific one, and the barely-binding steady
state stated. C6 — a paragraph added on flexible prices, real allocations and the
absence of a policy rate. C9 — the opening claim restated as an average rather than a
bound, and the debt ratio given a base date. C10 — the recessionary default branch
demoted from theorem to a property of the solved decision rules. Plus the language
items: "the paper" the central bank buys, the credit spread named, the episodic-risk
sentence connected to the global solution method, GHH preferences added, the
`fig:decoupling` redundancy removed, and the ECB's unwarranted-premia rationale
brought adjacent to the exogeneity restriction that rules it out.

**Actioned in `04-policy experiment.tex`.** Section title corrected to "Outright
Monetary Transactions and the TPI"; the comment header no longer reads "Conclusion".

**Left open, deliberately.**

- **C8, the results paragraph.** Not written. The findings that would fill it exist
  only in the sequence-space implementation, which is a different solution of a
  different model from the one Sections 2--3 describe; importing those numbers would
  create exactly the kind of inconsistency this report is about. The position is
  marked `[ADD RESULTS]`, matching the draft's existing placeholder convention, with
  a comment block giving the intended structure.
- **C7, Section 4.** Still a heading. The introduction's three questions and the road
  map both still promise it.
- **C5 and the notation ledger.** In `02-model.tex` and `03-calibration.tex`, which
  are the coauthors' global-model sections; not edited here. Outstanding: the
  "rollover risk" phrase that contradicts the exogenous-default framing, the
  `omega_X^f`/`omega_X` and `f_X`/`f` splits, the `lambda_X`/`Lambda` and `mu`/Smolyak
  collisions, the two duration formulas, the unstated recovery rate, and the empty
  recessionary-default proposition.
- **Compile defects in `02-model.tex`.** `\D` and `\F` are used at lines 43, 63 and
  142 but defined nowhere — `main.tex` carries only a commented-out macro block — so
  the document currently builds with five errors and those symbols drop silently from
  the PDF. There is also a stray `\\` after the proof at line 288. Left for the
  coauthors, since defining the macros fixes the build but fixes their notation for
  them. Three `XXX` cross-references remain undefined.
- **`references.bib`.** `bi2026asset` renders as (2026) against FRBSF Working Paper
  2025-10.

The paper compiles to 23 pages with all citations resolved; every remaining error and
undefined reference predates this edit and sits in Section 2.

---

## Verdict

**Major revision.** The introduction is a faithful and often elegant summary of Section
2's *mechanics* — the two-property transmission story, the branch-contingent kernel, the
three-way price decomposition, the priced-but-unrealised default framing all match the
model section closely, in places almost verbatim. The problems are not there.

They are in the two places where the introduction makes a promise the rest of the paper
does not keep: it justifies the choice of Greece by a measurement argument that Section 3
does not carry out, and it describes a policy experiment that exists nowhere in the
document. Underneath both sits a pattern of under-attribution to Bocola (2016), on whose
posterior the paper's entire financial block rests. Add to that a policy section that is a
title, an abstract that is Lorem Ipsum, and an introduction that reports no findings.

Findings C1, C2 and C7 are load-bearing. The rest are fixable in an afternoon.

---

## C1. The stated reason for choosing Greece is contradicted by the calibration section

This is the most serious inconsistency in the paper, and it is entirely internal.

The introduction (ll. 138–141) justifies the laboratory on measurement grounds, and does
so emphatically:

> We calibrate to Greece and Germany over 2010–12 **for reasons of measurement**. The
> 2011 EBA disclosures report sovereign exposures, maturity ladders and bilateral
> cross-holdings bank by bank, and the March 2012 restructuring supplies a realised
> recovery rate. **Elsewhere these objects have to be assumed.**

Section 3 assumes them anyway. The EBA is not cited in the calibration section at all. The
financial block is instead lifted wholesale from an estimated model of *Italy*:

| Object | Section 3's source |
|---|---|
| `theta_bar` (steady-state intermediary leverage) = 5 | "posterior mean of \citet{bocola2016pass}" |
| `f` (banker payout/exit share) = 0.04 | same |
| `varsigma_bar` (steady-state intermediation wedge) = 8 bp p.a. | "following \citet{bocola2016pass}" |
| `rho_s` (persistence of the latent risk factor) = 0.95 | "posterior means from \citet{bocola2016pass}" |
| `sigma_s` (innovation s.d. of the risk factor) = 0.63 | same |
| Smolyak grid half-width on `s` | "matching the coverage in \citet{bocola2016pass}" |

Section 3 is candid about it — "Neither has a clean accounting counterpart", and the wedge
"is the least identified parameter of the model, and we treat it as such" — which is
exactly the right register, and makes the introduction's claim look worse by contrast.

Three of the four EBA objects the introduction advertises never appear. Sovereign
exposures: absent. Bilateral cross-holdings: absent — the cross-border position `b^F_D`
that Section 2 puts on the balance sheet is **never calibrated anywhere in the paper**.
Maturity ladders: `delta_b` = 0.056 is set to a Greek debt duration, but no ladder and no
EBA citation. The recovery rate is the fourth, and Section 3 says only "We accordingly set
the recovery rate to match that event" after quoting *two* different numbers — a 53.5%
face-value reduction and a 59–65% NPV haircut. **The recovery rate is never stated.**

Either write the calibration the introduction promises, or rewrite the paragraph to say
what Section 3 actually does: the sovereign and fiscal side are disciplined to the Greek
episode, the intermediary block is imported from Bocola's Italian posterior, and Greece is
chosen because the default event is observed and the recovery rate realised. That second
version is still a good reason to pick Greece. It is just not the reason currently given.

## C2. The introduction under-credits Bocola (2016) three times

C1 is one instance of a pattern the referee will notice and dislike.

**First**, the literature paragraph (ll. 87–90) introduces Bocola as a limitation:
"\citet{bocola2016pass} measures the pass-through from sovereign risk to lending for
Italy, **in a model with a representative household**." The reader is invited to file it
under superseded. Section 3 then takes five parameters and a grid specification from it.
You cannot distance yourself from a paper's household block and adopt its financial
posterior in the same document without saying so.

**Second**, Section 2 labels its central analytical result
"\begin{proposition}[Closed-form multiplier, **from \cite{bocola2016pass}}]" — an explicit
attribution. The introduction (ll. 110–112) presents the same object with none: "The first
is the multiplier on the incentive constraint, **which has a closed form.**" Standing
alone in a contribution paragraph, that reads as a claim of novelty for a result the
paper's own model section credits elsewhere.

**Third**, `mu_closed` is described in Section 2 as "the analytical heart of the
transmission mechanism", and the introduction gives it top billing. If the heart is
borrowed, say so in the introduction — it costs one clause and removes the referee's best
line of attack.

## C3. The two transmission properties are presented in the order the calibration reverses

The introduction gives the multiplier channel primacy — "**The first** is the multiplier
on the incentive constraint" (l. 111) — the longer treatment, and the strongest verb: "The
working-capital wedge is **the channel** through which financial conditions reach output
on impact." The risk-premium channel comes "second" (l. 119) and gets four lines.

Section 3 then calibrates the steady state to `mu_bar` (steady-state IC multiplier) =
0.0010 and calls it, in its own words, "**a barely-binding constraint**". Section 2's
Endogenous Premiums subsection says of the other channel that "the *risk premium* operates
even when current constraints are slack", and that it is what generates "contractionary
deleveraging pressure prior to any actual default event."

So the paper's own calibration puts the economy at the edge of slackness, where the
channel the introduction ranks first is weakest and the channel it ranks second does the
work. A referee who reads Sections 2 and 3 will reverse the ordering and ask why the
introduction did not. The honest and more interesting version is that the constraint is
barely binding at the rest point and the anticipation channel is what moves the economy —
which is also a sharper contribution claim than the one currently made.

## C4. "Premium" is written as though sovereign-specific; in the model the constraint channel is asset-neutral

The introduction's second-property paragraph (ll. 119–124) is consistent with Section 2 —
`Omega^{(1)} > Omega^{(0)}`, the covariance term, the three-way decomposition, all of it
matches `eq:branch-sdf` and `eq:bond-decomposition`. But the prose reads as if the
banker's discount factor is what makes *sovereigns* special, and Section 2 says something
more specific.

`eq:foc-general` imposes `E[Omega(R_j − R)] = lambda_X mu_X` for **every** asset class —
capital, home sovereign, foreign sovereign, and the working-capital book alike. A single
divertibility parameter `lambda_X` (fraction of assets the banker can abscond with), one
constraint, one multiplier: the liquidity discount is identical across the portfolio.
Nothing in the constraint channel distinguishes sovereign paper from capital. What is
sovereign-specific is the covariance term, because only `Xi^D` (the D-sovereign per-unit
payoff) carries the survival factor `h_{t+1}`, and Section 2 states explicitly that
"$F$-bonds are safe in both states."

Say it that way. The clean statement — the constraint channel is asset-neutral, the risk
premium is sovereign-specific, and only the second can move a *spread* — is stronger than
what the introduction currently says, and it is what licenses Section 3's calibration
step, where a measured sovereign spread is used to bound `varsigma_bar`, a *capital*
excess return. As written, that step in Section 3 arrives unmotivated.

## C5. Section 2 says "rollover risk"; the introduction correctly says the model has none

Section 2 opens the sovereign block with: "Both governments issue [Hatchondo /
Chatterjee] perpetuities that are subject to \emph{rollover risk}." Four lines later:
"We model default risk as \emph{exogenous}", with a latent AR(1) factor `s_t` and a
logistic priced probability.

The introduction (ll. 149–152) gets this right and states the cost plainly: "Default risk
follows an exogenous latent factor. The government's own solvency calculus plays no part
in it, so the model admits no self-fulfilling component."

Rollover risk *is* the self-fulfilling channel. The introduction cites
\citet{bocola2019self} two paragraphs earlier for precisely the fundamental/rollover
decomposition, so the collision is visible to any reader who reaches Section 2. Delete the
phrase from Section 2; the introduction's version is the correct one and should not be
softened to accommodate it.

## C6. "Monetary union" is asserted in the introduction and quietly withdrawn in Section 2

The introduction's model paragraph opens: "We build a two-country monetary union"
(l. 97). It says nothing further about prices, inflation, nominal contracts, or monetary
policy — not one word in 158 lines.

Section 2's third sentence: "Prices are fully flexible, so all equilibrium objects below
are real. This implies that the ``monetary union'' is modelled at the level at which it
binds real allocations, namely a single union-wide funding market for intermediary
liabilities together with real interest parity."

Section 2's own scare quotes are the tell. A paper titled *Uncertain Unconventional
Policy*, whose object is a central bank instrument, cannot leave the reader to discover in
Section 2 that the model has no nominal side. State it in the introduction, as a modelling
choice with a defence — the mechanism is a real balance-sheet mechanism, and the nominal
block would add transmission the paper is not studying. Stated up front it is a
simplification; discovered later it looks like something withheld.

Two smaller casualties of the same silence. Section 2 posits "a single union-wide funding
market", yet `eq:foc-general` and `eq:rwc` carry country-indexed deposit rates `r_{X,t}`,
and Section 3 solves for `r_D` and `r_F` as separate market-clearing unknowns. And the
household budget deflates the predetermined deposit return by
`P^c_{X,t−1}/P^c_{X,t}`, an inflation term, one paragraph after "all equilibrium objects
below are real" — defensible, since the consumption basket price moves with the terms of
trade `p_t` under flexible prices, but it needs the half-sentence that says so.

## C7. The policy experiment described in the introduction exists nowhere in the paper

Lines 133–136:

> Into this environment we introduce a purchase rule of the TPI kind. We ask what an
> unlimited, country-specific backstop does to the equilibrium, which banking system
> supplies the paper the central bank buys, and how the resulting exposure and the
> resulting gains fall across the two countries and across households.

Section 2 has no central bank. It ends at the fiscal rule — no purchase rule, no
Eurosystem, no capital key. It also has no market-clearing conditions and no equilibrium
definition. Section 3 calibrates no policy parameter. Section 4 is a section heading:

```
% 4. Conclusion
\section{Outright Money Transations/TPI}
\label{sec:policy}
```

Six lines, of which two are a comment header that says "Conclusion", and a title
misspelling "Transactions" — and naming the instrument *Outright Money* Transactions,
where the introduction gets it right at l. 58. The road map then promises the reader that
"\cref{sec:policy}" delivers "the backstop counterfactual and its incidence."

The three questions at ll. 133–136 are the paper's reason for existing. Until Section 4
exists, the introduction is writing cheques on it.

## C8. The introduction reports no findings

A reader reaches the road map having been told what the model contains, why Greece, and
what will be asked — and not one thing that was learned. Every sentence in the results
position is interrogative.

`01-introduction.tex.bak` had a full results paragraph. Whatever prompted its deletion,
the replacement is an introduction that cannot be assessed: a referee cannot grade a
contribution that has not been stated, and an examiner will ask why the paper does not
want to say what it found. Restore a results paragraph, with `docs/referee_report_2.md`'s
M5 and M6 qualifications applied to the claims that overreached.

## C9. The 50 bp statistic does opposite work in the introduction and in Section 3

The introduction's opening sentence rests on it: "For the first decade of the euro, Greek
and German ten-year yields traded within half a percentage point of each other", footnoted
"Mean spread of 50 basis points over 1999Q1–2007Q4." It is offered as evidence of
convergence — the thing whose unwinding the paper is about.

Section 3 reaches for the same number for the opposite purpose: as an *upper bound* on
`varsigma_bar`, the intermediation wedge — "The Greek–German ten-year differential
averaged 50 basis points over 1999Q1–2007Q4" — and then discards it, setting
`varsigma_bar` = 8 bp on Bocola's posterior instead.

The same statistic cannot be the paper's motivating fact and a bound the paper declines to
use. Also note the logical form of the introduction's sentence: "traded within half a
percentage point" asserts a *bound*; the footnote supplies a *mean*. A mean of 50 bp is
consistent with wide excursions. Either state the claim as an average, or defend the bound
with a range.

Section 3's footnote on this line is broken mid-sentence and will compile as such:
"\footnote{FRED series ...; the $2010\text{Q}1$--$2012\text{Q}2$}."

## C10. The introduction asserts a proposition Section 2 states without content or proof

Introduction, ll. 129–131: "The recession in the default branch is itself endogenous,
arising from the same balance-sheet mechanism operating on a smaller asset base. The model
imposes no exogenous output cost of default." This tracks Section 2 almost word for word,
so as a summary it is faithful.

But in Section 2 the claim is a formal object:

```
\begin{proposition} [Under standard assumptions, default is recessionary]
\end{proposition}
Proof in the appendix. [ADD]
```

An empty proposition body, a proof that does not exist, and no appendix containing it. The
introduction presents as settled the one result the paper has flagged as owed. Either
prove it, demote it to a numerical finding, or state it in the introduction as a property
of the solved model rather than a theorem.

---

## Notation and cross-section consistency

Purely internal; all of these will be caught on a careful read.

| Item | Problem |
|---|---|
| `omega_X^f` vs `omega_X` | The entrant transfer is `omega_X^f` in `eq:nw-lom` and `omega_X` in `eq:bank-div` — adjacent equations — and `omega_X` again in Section 3. One object, two symbols. |
| `f_X` vs `f` | Subscripted in `eq:branch-sdf` and the banker's problem, bare in `eq:omega`, `eq:nw-lom` and `eq:ss-closedform`. |
| `delta_{X,b}` vs `delta_b` | Subscripted once, in the default subsection; bare in `eq:bond-payoff`, `eq:govt-budget` and throughout Section 3. |
| `lambda_X` vs `Lambda_{X,t,t+1}` | Divertible asset fraction and household stochastic discount factor, differing only by case, and appearing in the same display (`eq:omega` beside `eq:foc-general`). Rename one. |
| `mu_{X,t}` vs Smolyak `mu` | Section 3's solution method sets "an isotropic Smolyak level $\mu=1$" — the same glyph as the rescaled IC multiplier, which Section 2 calls the analytical heart of the model. |
| `a_{it}` vs `a_{j,t}` | Household deposits in `eq:hh-budget`; the market value of the bank's position in asset class *j* in `eq:ng-excess`. |
| `w_{X,t}` units | `eq:hh-budget` deflates it by `P^c_{X,t}`, implying nominal; `eq:labour-demand` delivers it in own-good units, implying real. One of the two is wrong. |
| Bond duration | Section 2: "approximately $1/\delta_{X,b}$ quarters". Section 3: "$(1+\bar r)/(\bar r+\delta_b)$ quarters". Two formulas for one object, never reconciled. |
| Country naming | Greece/Germany (introduction), `D`/`F` (Section 2), "Periphery"/"Core" (Section 3, which announces a grouping strategy the introduction's single-country-pair framing contradicts). |
| Section 3 title | "Empirical Analysis", for a section containing a calibration and a solution method. The introduction's road map describes the contents correctly; the title matches neither. |
| Weak vs strict IC | Introduction: depositors lend "only while the franchise is worth **more than** the divertible proceeds". `eq:IC` is weak (`\ge`). |

## Language and copy

- **The abstract is "Lorem Ipsum."**
- Three unresolved cross-references in Section 2 — `\cref{XXX}` (l. 35),
  `Section~\ref{XXX}` (l. 150), `Section~\ref{sec:XXX}` (l. 158) — two of which point at
  the bond-payoff and predetermined-rate conventions the introduction relies on.
- Section 4's comment header says "4. Conclusion" above a section titled "Outright Money
  Transations/TPI". Two errors and a missing section in six lines.
- Section 3 typos and broken sentences: "Perihpery"; "broadstrokes"; "To study its
  quantitative needs" (garbled); "the wealth distribution, is pinned down" (stray comma);
  "Since \eqref{eq:labour-demand} is the only channel from the credit spread to output on
  impact." (fragment, no main clause); "the spread of the sovereign yield over the
  risk-free rate, and any measured spread is an upper bound on it" (missing verb — the
  sentence never says what the spread *equals*).
- Introduction, l. 135: "which banking system supplies **the paper** the central bank
  buys" — "paper" meaning sovereign bonds, in a document that calls itself the paper four
  times. Use "the bonds".
- The introduction never names the object Section 2 calls "the model's *credit spread*"
  and "the single instantaneous quantity that links the financial block to the production
  block." It has a name; use it.
- Redundancy: ll. 46–51 and ll. 65–72 both make the no-variation-to-estimate-against
  point off `fig:decoupling`, and the figure is discussed at l. 46 but placed at l. 74.
- Missed connection: "Sovereign risk inside a currency union is episodic, close to zero
  for long stretches and very large for short ones" (ll. 50–51) is the argument for
  Section 3's global nonlinear solution — an occasionally binding constraint interacting
  with a rare event. The introduction never cashes it, and never mentions the solution
  method at all, though the road map promises one.
- The introduction omits two substantive choices from its model paragraph: GHH
  preferences, which remove the wealth effect on labour supply and which Section 3 leans
  on for identification, and the monopolistically competitive retail layer, which under
  flexible prices contributes only a constant markup.
- ll. 138–147 and ll. 149–153 remain the strongest writing in the draft — the omitted
  causes of the Greek depression, the concession on TPI eligibility, and the exogeneity
  restriction stated with its cost. `docs/referee_report_2.md`'s M1 and M2, answered about
  as well as they can be without changing the model. Keep every line.

---

## What to do

1. **Rewrite ll. 138–141.** Say what Section 3 calibrates. Greece is the right laboratory
   because the default event is observed and the recovery realised — not because the bank
   block is measured, which it is not. (C1)
2. **Credit Bocola (2016) in the introduction**, once in the literature paragraph and once
   on the closed-form multiplier. (C2)
3. **Swap the order of the two transmission properties**, or add the sentence that says
   the constraint is barely binding at the rest point and the anticipation channel is what
   moves the economy. (C3, C4)
4. **Add two sentences on the nominal side** — flexible prices, real allocations, no
   policy rate, and why that is the right abstraction here. (C6)
5. **Write Section 4, or cut ll. 133–136 back to a single forward-looking sentence.** The
   present text promises a counterfactual, an incidence analysis, and a
   who-sells-to-the-central-bank result, none of which the document contains. (C7)
6. **Restore a results paragraph.** (C8)
7. **Delete "rollover risk" from Section 2**, fix the notation ledger, resolve the three
   `XXX` references, state the recovery rate, and either prove or demote the recessionary-
   default proposition. (C5, C10)
