# Referee Report — Fourth Referee

## On Section 2, "The Model", read against Sections 1 and 3

**Scope.** `VIVA/sections/02-model.tex` checked for internal consistency, algebraic
correctness, completeness as a model section, and agreement with
`01-introduction.tex` and `03-calibration.tex`. Every equation in the section was
re-derived. No comparison to the repository's sequence-space implementation:
`docs/referee_report_3.md` established that Section 2 describes the coauthors' global
model and that the two are separate objects, and that finding is respected here.

Line numbers refer to `02-model.tex` as of 2026-08-27. Equation numbers are the
compiled ones from `main.pdf` (23 pages) where I give them; otherwise I use the
labels.

---

## Verdict

**Major revision.** The algebra that is present is, with one exception, correct — I
re-derived the CES demand system, the Jermann adjustment-cost block, the
excess-return representation, Proposition 1, the bond-pricing formula, the credit
spread and the two government equations, and they all check out. Section 3's
closed-form steady state reproduces Section 2's equations exactly (`\bar\mu` = 0.000996
against the reported 0.0010, the `\omega_X` > 0 bound at 3.8667% against the reported
3.87%). That is more than most model sections survive.

The problems are of a different kind. **Section 2 does not define an equilibrium.**
There is no market clearing of any kind in it — no resource constraint, no bond
clearing, no deposit clearing, no labour clearing, no parity condition, no
productivity process, no household first-order condition, and no statement of how the
household continuum aggregates. What the section contains is a list of agents'
problems. A reader cannot tell what the model *is*, and cannot check that the seven
unknowns Section 3 solves for are the right seven.

Underneath that sit two substantive defects. The bank's portfolio problem as written
is **over-determined**: four Euler equations per bank against two bond prices, with no
adjustment cost and no corner, and the cross-border position `b^F_D` (D-bank holdings
of F-government paper) is never determined by anything. And the sentence that
explains the risk-premium mechanism (l. 411) **contradicts itself and inverts the
sign** of the object it is explaining.

Findings **M1, M2, M3** are load-bearing. **M4–M9** need decisions rather than
repairs. The rest is a long afternoon.

Two questions the section provokes and does not answer, taken up as **M9**: default is
**not** endogenised — the sovereign's repayment decision is not modelled anywhere, and
the section says so (l. 349) while also calling the debt subject to "rollover risk"
(l. 345, see **m3**) — and there is **no central bank in the paper at all**, so the
instrument the title names has nowhere to attach.

---

# Major findings

## M1. Section 2 contains no equilibrium

The section ends at the fiscal rule. It never states what clears. Missing, in full:

| Object | Where it should be | Consequence of its absence |
|---|---|---|
| Goods-market clearing / resource constraint, either country | after §2.2 | `p_t` (the terms of trade) has nothing to clear; the trade balance is undefined |
| Bond-market clearing, `b^{X}_{D,t+1}+b^{X}_{F,t+1}=B_{X,t+1}` | §2.5 | the link between the government's `B_{X,t}` (total stock outstanding) and the banks' `b` holdings is never made — they are separate symbols that never meet |
| Deposit-market clearing, `\int a_{it}\,di=\mathrm{dep}_{X,t}` | §2.1 or §2.3 | `r_{X,t}` (the deposit rate) has nothing to clear; `\mathrm{dep}` is defined twice, once as the bank's residual (l. 335) and once implicitly as household saving (l. 32) |
| Labour-market clearing, `N_{X,t}=\int e_{it}n_{it}\,di` | §2.1 | `N_{X,t}` (aggregate effective hours) in \eqref{eq:production} is never connected to `n_{it}` (household `i`'s hours) in \eqref{eq:hh-objective} |
| Capital-market clearing | §2.3 | asserted in prose at l. 132 ("intermediates *all* of the economy's productive capital"), never as an equation |
| Real interest parity | §2 opening | see **M5** |
| The process for `Z_{X,t}` (country-`X` TFP) | §2.2 | `Z_{X,t}` enters \eqref{eq:production} and is a state variable in Section 3's grid, and its law of motion appears nowhere in the paper |
| The household Euler equation and the KKT conditions on `a_{i,t+1} \ge 0` | §2.1 | the household problem is posed (eqs. 1–3) and never solved |
| A definition of recursive competitive equilibrium | §2 close | — |

Line 290 says "the equations imposed in equilibrium", which is the only occurrence of
the word. Section 3 then asserts (l. 78–80) that at each collocation point "the seven
market-clearing unknowns `[N_D, N_F, K_D', K_F', r_D, r_F, p]` are solved" — against
market-clearing conditions the paper never writes. A referee cannot verify the count,
and neither can a reader trying to replicate.

This is the single largest thing wrong with the section, and it is also the easiest to
fix: the conditions exist in the authors' code, and they need about half a page.

## M2. The bank's portfolio problem is over-determined, and the cross-border position is determined by nothing

\eqref{eq:foc-general} imposes
`\E_t[\Omega_{X,t,t+1}(R_{j,t+1}-R_{X,t})] = \lambda_X\mu_{X,t}` for **every**
`j \in \mathcal{J}_X = \{K_X, b^D, b^F, L_X\}` — four conditions for the D-bank and
four for the F-bank. Both banks hold both sovereigns (l. 132, l. 136). So both banks
price both bonds:

- D-bank, own leg: `\E_t[\Omega_D(\Xi^D_{t+1}/q^D_t - R_{D,t})] = \lambda_D\mu_{D,t}`
- F-bank, cross leg: `\E_t[\Omega_F(\tfrac{1}{p_{t+1}}\Xi^D_{t+1}/(\tfrac{1}{p_t}q^D_t) - R_{F,t})] = \lambda_F\mu_{F,t}`

Two conditions, one price `q^D_t`. Symmetrically for `q^F_t`. `\Omega_D \ne \Omega_F`
(different countries, different consumption paths, and only D's sovereign defaults),
and the cross leg additionally carries the `p_{t+1}/p_t` conversion, so these are not
the same equation. **Generically they cannot both hold.**

Section 2 does not notice. \eqref{eq:bond-pricing} presents `q^D_t` as determined by
D-bank optimality and `q^F_t` by F-bank optimality, and simply never writes the two
cross-border conditions down. Line 290 promises "Specialising \eqref{eq:foc-general}
to each class yields the equations imposed in equilibrium" and then delivers three of
the four classes: Capital, Sovereign bonds (own only), Working capital. **The
cross-border leg is silently dropped.**

The consequence is not cosmetic. Nothing in the paper determines `b^F_{D,t+1}` or
`b^D_{F,t+1}`. They are on the balance sheet \eqref{eq:balance-sheet}, they enter
`\mathcal{A}_{X,t}` (total assets) and therefore the incentive constraint
\eqref{eq:IC} and the multiplier \eqref{eq:mu-closed}, they enter gross wealth
\eqref{eq:asset-payoff}, and they are absent from Section 3's seven unknowns. The
introduction (l. 99–101) rests the paper's claimed contribution on exactly this
position — "none holds the sovereign and the productive capital stock on the same
constrained balance sheet in both countries at once" — and the model section leaves it
free.

Three ways out, and the paper must pick one and say so:

1. **A portfolio adjustment cost on the cross-border leg**, which turns the cross
   conditions into quantity restrictions and leaves the own legs to pin prices. This
   is the standard fix and it is what the authors' own sequence-space implementation
   does.
2. **A corner**: cross-border holdings fixed exogenously, in which case
   \eqref{eq:foc-general} must be restricted to `j \in \{K_X, b^{own}, L_X\}` and
   `\mathcal{J}_X` redefined, and the balance-sheet leg becomes an endowment.
3. **Segmented markets** with the cross leg priced off a separate condition.

Whichever is chosen, `\mathcal{J}_X` at l. 189 is currently wrong.

## M3. The sentence explaining the risk premium contradicts itself, and the sign is backwards

Line 411, immediately after \eqref{eq:branch-sdf}:

> Because the GHH composite is lower in the default branch, `x^{(1)}<x^{(0)}`, and
> hence `\Omega^{(1)}>\Omega^{(0)}`: the bank values wealth more in the bad state.
> This is what makes the bond carry a genuine *risk premium* rather than a pure
> actuarial discount, and it is what makes `\E_t[\Omega]` in \eqref{eq:mu-closed}
> **fall** when risk rises --- tightening the constraint.

The first half is right. The second half does not follow from it; it follows the
opposite way. By the section's own \eqref{eq:expectation},

```
E_t[Omega] = (1 - pi^d_t) Omega^(0) + pi^d_t Omega^(1),
d E_t[Omega] / d pi^d = Omega^(1) - Omega^(0) > 0
```

by the premise stated one clause earlier. Raising the priced default probability
shifts weight onto the branch the sentence has just declared to have the **higher**
`\Omega` (the banker's augmented stochastic discount factor). `\E_t[\Omega]` **rises**.

This matters because \eqref{eq:mu-closed} is
`\mu = \max\{1-\E_t[\Omega]R_{X,t}n_{X,t}/(\lambda_X\mathcal{A}_{X,t}),0\}`, in which
`\partial\mu/\partial\E_t[\Omega] = -R n/(\lambda\mathcal{A}) < 0`. A higher
`\E_t[\Omega]` **loosens** the constraint. So the clause is wrong twice over: the
premise gives the opposite sign, and the conclusion drawn from it is the opposite of
what \eqref{eq:mu-closed} says.

The clause is also unnecessary. Two paragraphs earlier (l. 288) the section already
states the correct mechanism — "A fall in sovereign bond prices reduces `n_{X,t}`
(intermediary net worth) through \eqref{eq:ng}, which *mechanically* raises
`\mu_{X,t}` (the multiplier on the incentive constraint)" — and that is the channel
that does the work. Deleting the clause after the em-dash loses nothing.

There *is* a defensible statement in the neighbourhood, and the paper should make it
rather than assume it: `\Lambda^{(d')}=\beta(x_{X,t}/x^{(d')}_{X,t+1})^{\sigma_X}`
also contains today's composite `x_{X,t}`, which falls on impact, and if today's
contraction outweighs the expected future one then `\E_t[\Omega]` can fall on net.
That is a horse race between a level effect and a composition effect, it is
quantitative, and the section currently asserts the answer while giving the reason
that points the other way. Worth noting for the authors: because the composition
effect pushes `\E_t[\Omega]` up, the risk-premium channel *offsets* the net-worth
channel inside the multiplier. That is an interesting property of the model. It is
currently written as the reverse.

## M4. The household heterogeneity does nothing, and Section 3's state vector confirms it

The paper is titled "A Two-Country **Heterogeneous-Agent** Evaluation". The
introduction (l. 119–121) promises "a continuum differing in individual labour
productivity, with Greenwood--Hercowitz--Huffman preferences, a no-borrowing
constraint, and bank deposits as their only savings vehicle", and Section 2.1
delivers exactly that: household `i`, idiosyncratic `e_{it}`, `a_{i,t+1}\ge0`.

Then the heterogeneity is never used again.

- Section 2 never aggregates. `C_{X,t}` appears in \eqref{eq:ces-aggregator} with no
  `i`; `N_{X,t}` appears in \eqref{eq:production} with no `i`; nothing connects them
  to `c_{it}`, `n_{it}`.
- Line 194 says the banker discounts "with the aggregate household stochastic discount
  factor `\Lambda_{t,t+1}`", and \eqref{eq:branch-sdf} then *defines* that object as
  `\beta_X(x_{X,t}/x^{(d')}_{X,t+1})^{\sigma_X}` on the **aggregate** composite. With
  `a_{i,t+1}\ge0` binding for some households, `\beta(X_t/X_{t+1})^{\sigma}` is not
  the marginal rate of substitution of any agent in the model, constrained or
  unconstrained. It prices nothing. It is an assumption, and a substantive one,
  because the *entire* risk premium in \eqref{eq:bond-decomposition} is the curvature
  of this object.
- Section 3's aggregate state is `[K_D, K_F, P_D, P_F, B_D, s, Z_D]` (l. 78). **There
  is no distribution in it.** A Smolyak grid over seven scalars cannot carry a wealth
  distribution.

Either (a) the aggregate block is representative-agent and the heterogeneity is
layered on ex post to compute the incidence the introduction promises — in which case
Section 2 must say so explicitly, state that the distribution does not feed back into
aggregates, and justify the SDF; or (b) it is not, and the solution method is
misspecified. As written the reader cannot tell which, and the title asserts the
answer that the state vector denies.

If (a), note the cost: the whole point of a HANK block is that the distribution moves
aggregate demand. If it is shut off, the paper should not claim it, and the phrase
"heterogeneous-agent" in the title is doing work the model does not do.

## M5. "A single union-wide funding market together with real interest parity" is asserted and never written, and Section 3 contradicts it

Line 8, the paragraph that defines what a monetary union means here:

> the "monetary union" is modelled at the level at which it binds real allocations,
> namely a single union-wide funding market for intermediary liabilities together with
> real interest parity.

Repeated verbatim in the introduction (l. 104–106). **Neither object appears in the
model.** There is no parity condition anywhere in Section 2, and there are two deposit
rates: `r_{D,t}` and `r_{F,t}` are separate symbols throughout, and Section 3 solves
for **both** as separate unknowns (l. 78). Two rates in two markets is the definition
of a *segmented* funding market. Households, moreover, hold only domestic deposits
(l. 32, "Real deposits `a_{it}` are the only saving vehicle"), with no cross-border
household position.

So the sentence that tells the reader what makes this a monetary union rather than two
open economies is false as stated. Either write the parity condition
(`1+r_{D,t}=(1+r_{F,t})\,p_{t+1}/p_t`, or whatever the intended form is) and drop one
unknown, or delete the claim and describe the union as what it actually is: a common
numéraire with cross-border bank holdings of sovereign paper and no other integration.
The choice changes the count in **M1**.

## M6. Working capital is intra-period in the firm's problem and a one-period asset on the balance sheet

Line 79: the firm "must pre-finance a fraction `\zeta_X` of its wage bill with
**intra-period** bank credit at the gross rate `1+r^{wc}_{X,t}`". The static problem
\eqref{eq:firm-problem} charges `(1+\zeta_X r^{wc}_{X,t})w_{X,t}N_{X,t}` at date `t`,
undiscounted. That is the Neumeyer–Perri convention and it is correct *for
intra-period credit*: borrow at the start of `t`, repay at the end of `t`.

But the same loan then appears:

- as an end-of-period-`t` **asset** in the balance sheet, `+L_{D,t}` inside
  `\mathcal{A}_{D,t}` \eqref{eq:balance-sheet};
- as a `t+1` receivable in the obligation state,
  `P_{X,t+1}=(1+r_{X,t})\mathrm{dep}_{X,t}-(1+r^{wc}_{X,t})L_{X,t}` \eqref{eq:P-state};
- as an asset class `L_X \in \mathcal{J}_X` with gross return `R_{L,t+1}=1+r^{wc}_{X,t}`
  in the excess-return representation (l. 189–191), which is definitionally a
  `t \to t+1` claim;
- and with a `t \to t+1` Euler equation \eqref{eq:rwc}.

These cannot both be true. If the credit is intra-period it is not on the end-of-period
balance sheet, it does not consume incentive-constraint capacity across periods, and
`\mathcal{J}_X` has three members. If it is a `t \to t+1` asset then the firm's
date-`t` cost in \eqref{eq:firm-problem} is wrong — the repayment falls at `t+1` and
must be discounted — and the "intra-period" wording at l. 79 must go.

This is load-bearing because `L_X` sits inside `\mathcal{A}_{X,t}`, which is the
divertible base in \eqref{eq:IC} and the denominator of \eqref{eq:mu-closed}. Its
treatment changes the multiplier, and the multiplier is the whole transmission
mechanism.

Line 165 is also simply broken English: "Because the working-capital loan is extended
at rates the period at the rate locked at $t$".

## M7. Equation \eqref{eq:bond-decomposition} is not a decomposition

The algebra is right — I checked each step from \eqref{eq:euler-bond} through the
covariance identity — but the result is

```
q = E[Xi]/R  +  Cov(Omega,Xi)/(E[Omega] R)  -  [lambda*mu/(E[Omega] R)] * q
```

with `q^D_t` on **both sides**. The third term is not an additive component of the
price; it is the price itself, scaled. Solving properly gives

```
q = ( E[Xi] + Cov(Omega,Xi)/E[Omega] ) / ( R + lambda*mu/E[Omega] )
```

which is \eqref{eq:bond-pricing} rewritten — the constraint enters the *denominator*,
not as a subtracted term. As a consequence the two terms the paper labels "expected
payoff" and "risk premium" are not the corresponding components of `q^D_t`; they are
the components of a different, higher number.

The introduction leans on this hard: "The bond price separates into an expected
payoff, a risk premium and a constraint discount" (l. 140–141). It does not separate.
It can be *bounded* or *approximated* to first order — for the calibrated steady state
the constraint term is 1.99 bp of `q`, so a first-order treatment is numerically
harmless — but the paper should present it as such rather than as an identity.

Two smaller things in the same passage. Line 310: "even a default-free bond trades
below its risk-neutral present value **by the factor** `\lambda\mu`" — the discount
factor is `\E_t[\Omega]R/(\E_t[\Omega]R+\lambda\mu)`, not `\lambda\mu`, which is not
even dimensionless in the right way. And the `(<0)` annotation on the risk-premium
term in \eqref{eq:bond-decomposition} is asserted; it is defensible (in the default
branch both `h_{t+1}` and `q^D_{t+1}` fall while `\Omega` rises, and within the
no-default branch a higher `s_{t+1}` does the same) but it should be shown in a line.

## M8. The section asserts its central quantitative result and never quantifies it

Lines 99–102:

> Notice that equation \eqref{eq:labour-demand} is the **only** channel through which
> financial spreads, affecting `r^{wc}_{X,t}`, reach output on impact. Setting
> `\zeta_X=0` nests the model without it exactly, and doing so collapses the output
> response to a sovereign-risk shock to approximately zero even when bond prices fall
> significantly.

Three problems.

**(a) The claim is a result, stated in the model section, with no number and no
forward reference.** Section 4 does not exist. There is nothing for the reader to
check it against.

**(b) "Nests exactly" is not exact.** At `\zeta_X = 0` there is no loan book, so
`L_X` leaves `\mathcal{J}_X` and leaves `\mathcal{A}_{X,t}`. Leverage `\theta_{X,t}`,
the divertible base, `\lambda_X` and `\omega_X` (all calibrated *to* `\bar\theta = 5`)
all change. It nests the transmission channel, not the model.

**(c) A referee will compute the elasticity, and the authors should do it first.**
Combining \eqref{eq:labour-demand} with Section 3's GHH labour supply
`\chi_X n_{it}^{1/\nu_X}=(w/P^c)e_{it}` and \eqref{eq:production} at predetermined
`K_{X,t}`, holding the terms of trade fixed:

```
d ln N = -[ nu / (1 + alpha*nu) ] * zeta * d r_wc,     d ln Y = (1-alpha) d ln N
```

At `\alpha=0.35`, `\nu=1`, `\zeta=1`:

| Rise in the credit spread `\lambda\mu/\E[\Omega]` | Impact `d\ln Y` |
|---|---|
| +100 bp p.a. | −0.120% |
| +500 bp p.a. | −0.602% |
| +1000 bp p.a. | −1.204% |
| +2500 bp p.a. | −3.009% |

The introduction opens on a 22-point output fall. Reproducing it through this channel
needs a credit spread of roughly **18,300 bp per annum**. The introduction does hedge
— "its aggregate path lies outside what the model is built to reproduce" (l. 171–172)
— and that hedge is honest, but it is in Section 1 and the reader meets the claim in
Section 2. Section 2 should state the elasticity, because it is a closed form, it
takes two lines, and every referee will derive it.

(The Frisch elasticity `\nu_X`, `\alpha_X` and `\delta_X` are never given a value
anywhere in the paper — see **m12** — so I have used conventional ones. If `\nu_X = 2`
the impact response roughly doubles and the conclusion is unchanged.)

## M9. There is no central bank, and the model has no slot to put one in

Section 2's agents are: households, retailers, intermediate producers, capital
producers, one representative bank per country, and two governments. That is the whole
list. The words "central bank", "ECB", "purchase", "backstop", "TPI" and "policy"
appear **nowhere** in the section. The introduction is consistent with this — "There
is no policy rate and no inflation" (l. 106) — and the instrument enters the paper only
as a promise: "\Cref{sec:policy} adds an unlimited, country-specific backstop of the
TPI kind to this environment" (l. 154). `04-policy experiment.tex` is five lines: a
comment banner, a `\section` and a `\label`.

Ordinarily "Section 4 is empty" is report 3's finding (its **C7**), not a model-section
finding. It becomes one here because **the two holes identified above are exactly the
two places a sovereign backstop has to attach**, and neither exists yet:

| What a TPI purchase does | Where it enters | Status in Section 2 |
|---|---|---|
| Absorbs a quantity of `D` paper | bond-market clearing, `b^D_{D}+b^D_{F}+b^D_{CB}=B_D` | **the condition is not in the paper** (**M1**) |
| Takes it from *one of the two banking systems* | the banks' cross-border and own legs | **`b^F_D` and `b^D_F` are determined by nothing** (**M2**) |
| Relieves the incentive constraint | `\mathcal{A}_{X,t}` in \eqref{eq:IC}, hence `\mu_{X,t}` | mechanism present, entry point absent |
| Books a gain or loss, and remits it | the government budget \eqref{eq:govt-budget} | no remittance term, and no rule for who bears CB losses |
| Is funded by *something* | the union funding market | **asserted at l. 8, never written** (**M5**) |

Which banking system supplies the bonds is one of the paper's three stated questions
(l. 154–158, "which of the two banking systems supplies the bonds bought under it").
That question is a *portfolio* question. It is answered by the cross-border Euler
conditions — the ones \eqref{eq:foc-general} promises and l. 290 quietly omits. So
**M2 is not a tidiness complaint: it is the equation that answers the paper's own
research question**, and until it is written the policy section cannot be.

The same goes for the fiscal side. \eqref{eq:hh-budget} taxes households lump-sum, and
\eqref{eq:bohn} makes the tax a function of the surviving debt stock. A central bank
holding `D` paper changes that stock's ownership but not its size, so under the
current rule the backstop has *no* fiscal consequence for households except through
`q^X_t`. Whether that is the intended design or an accident of the rule needs
deciding before Section 4, not during it — and it will drive the household incidence
result (see **m8**).

A related point of framing. Section 2 opens (l. 8) by defining the monetary union as
"a single union-wide funding market … together with real interest parity", which is
the *minimum* structure that lets a union-wide institution exist. The paper then never
writes that structure down (**M5**) and never introduces the institution. As it
stands, Section 2 describes two open economies sharing a numéraire, with cross-border
bank holdings of sovereign paper and nothing else in common. That is a defensible
model. It is not yet one in which "the ECB" is a well-defined agent.

---

# Moderate findings

## m1. The section does not compile, and the failure corrupts the balance sheet

`\D` and `\F` are used at ll. 43, 63 and 142 and defined nowhere; `main.tex`'s macro
block (ll. 7–9) is commented out. The build emits four `Undefined control sequence`
errors and LaTeX drops the symbols silently. This was flagged in
`docs/referee_report_3.md` and left "for the coauthors". It has not been fixed, and it
is worse than a warning, because line 142 is inside \eqref{eq:balance-sheet}. What the
PDF actually prints is

```
Q_{D,t} K_{D,t+1} + q^{D}_{t} b^{D}_{D,t+1} + p_t q_{t} b_{D,t+1} + L_{D,t}
```

The cross-border leg loses its issuer superscript entirely, so the two bond terms
become typographically indistinguishable and the reader has no way to recover the
`b^{issuer}_{holder}` convention — which, note, is never stated in words either. The
most important equation in the section is unreadable in the compiled paper.

Lines 43 and 63 print "The consumption basket of country ␣ aggregates…" and "expressed
in ␣-good units".

Fix: `\newcommand{\D}{D}` `\newcommand{\F}{F}` in `main.tex`, or replace the three
uses. Two minutes, and it should not wait for another round.

## m2. \eqref{eq:bank-lagrangian} does not attach a multiplier to \eqref{eq:IC}

Line 225 says "attach a multiplier `\tilde\mu_{X,t}\ge0` to \eqref{eq:IC}", i.e. to
`\varphi_{X,t}n_{X,t}\ge\lambda_X\mathcal{A}_{X,t}`. The Lagrangian written is

```
L = (1 + mu~) E[Omega n^g] - mu~ lambda A  =  E[Omega n^g] + mu~ ( E[Omega n^g] - lambda A )
```

whose constraint is `\E_t[\Omega n^g_{t+1}]\ge\lambda_X\mathcal{A}_{X,t}`, not
\eqref{eq:IC}. The two coincide only *at the optimum*, via
`\varphi n = \E_t[\Omega n^g]`, which is the value function the banker is in the
middle of solving for. The manipulation is standard and the answer is right — I
verified `\partial L/\partial a_j` gives \eqref{eq:foc-general} with the rescaling
\eqref{eq:mu-rescale}, and Proposition 1 follows — but the sentence describing it is
wrong and a careful reader will stall on it. State the constraint in the
`\E_t[\Omega n^g]` form and note the equivalence.

## m3. "Rollover risk" contradicts the default specification two lines later

Line 345: "Both governments issue \citet{Hatchondo2009}/\citet{Chatterjee2011}
perpetuities that are subject to *rollover risk*." Line 349: "We model default risk as
*exogenous*." The introduction is explicit that this is deliberate: "The government's
repayment decision is not modelled, so there is no self-fulfilling equilibrium"
(l. 177–178). Rollover risk in this literature (Cole–Kehoe, Bocola–Dovis) *is* the
self-fulfilling kind. Neither cited paper is about rollover risk. Delete the phrase.
Flagged in report 3; still present.

## m4. Proposition 2 is an empty environment with no proof and no appendix

Lines 376–378:

```latex
\begin{proposition} [Under standard  assumptions, default is recessionary]
\end{proposition}
Proof in the appendix. [ADD]
```

The environment has no body, so the compiled paper prints "Proposition 2 (Under
standard assumptions, default is recessionary)." followed by nothing. There is no
appendix: `main.tex` has `\startappendix` immediately before `\end{document}`. The
proposition is never cross-referenced.

The introduction has already retreated from it — "in the solved decision rules the
default state is recessionary" (l. 150–152), a numerical claim, not a theorem. Either
state and prove the proposition or delete it and cite the decision rules. Leaving a
titled, numbered, empty theorem in a submitted draft is the kind of thing that decides
a desk reject.

## m5. \eqref{eq:expectation} hides the branch-dependent state — which is where the risk premium comes from

The quadrature is written `\E_t[g] = (1-\pi^d_t)\E^s_t[g(0,s_{t+1})] +
\pi^d_t\E^s_t[g(1,s_{t+1})]`, as though `s_{t+1}` were the only argument. It is not.
The aggregate state is seven-dimensional, and at least three of its components are
**branch-dependent**: on default the debt stock is written down through `h_t` in
\eqref{eq:debt-lom}, and both intermediaries' obligation states `P_D`, `P_F` differ
because \eqref{eq:asset-payoff} pays `\Xi^D_{t+1}=h_{t+1}[\delta_b+(1-\delta_b)q^D_{t+1}]`.
That difference is precisely what makes `x^{(1)}<x^{(0)}` and therefore
`\Omega^{(1)}>\Omega^{(0)}` at l. 411 — the risk premium in
\eqref{eq:bond-decomposition} is *nothing but* that state difference. Writing
`g(d', s_{t+1})` erases it from the notation at the exact point where the reader needs
to see it. Section 3 (l. 80) is better on this than Section 2 is.

Related, and unstated anywhere: **what is the economy after a default?** Line 369
calls it "a single deterministic bond face value write-off", which reads as
once-and-for-all, but then the post-default economy must have a different risk process
(`\pi^d \to 0`? an absorbing state? `s` continues and default can recur?). If `s`
continues unchanged, the model prices repeated defaults on an already-written-down
stock, which is a different object from the Greek 2012 exchange the calibration is
built on. Section 2 must say which.

## m6. \eqref{eq:mpk} is not a first-order condition of \eqref{eq:firm-problem}

Line 87 says the firm's problem "deliver[s] the factor-demand conditions" and lists
both \eqref{eq:labour-demand} and \eqref{eq:mpk}. But \eqref{eq:firm-problem} is a max
over `N_{X,t}` alone — `K_{X,t}` is predetermined and `mpk_{X,t}K_{X,t}` is a
subtracted constant. Differentiating in `N` gives \eqref{eq:labour-demand} only.
\eqref{eq:mpk} is the **zero-profit** condition: imposing
`mc\,Y-(1+\zeta r^{wc})wN-mpk\,K=0` and substituting \eqref{eq:labour-demand} gives
`mpk\,K = mc\,Y - mc(1-\alpha)Y = mc\,\alpha Y`, which is \eqref{eq:mpk}. Correct
result, wrong derivation as described. One sentence to fix.

## m7. \eqref{eq:nw-lom} is simultaneous in `\mathcal{A}_{X,t}`

`n_{X,t}=(1-f)n^{g}_{X,t}+\omega^f_X\mathcal{A}_{X,t}`, where
`\mathcal{A}_{X,t}=n_{X,t}+\mathrm{dep}_{X,t}` is *end-of-period* assets, chosen after
net worth is known. So `n_t` depends on `\mathcal{A}_t` depends on `n_t`. In
Gertler–Karadi the entrant transfer is normally a fraction of the previous period's
assets, or of exiting bankers' assets, precisely to avoid this.

It is not fatal — I verified that Section 3's `\omega_X = \mathcal{D}/\bar\theta -
(1-f)\bar\varsigma` with `\mathcal{D}=1-(1-f)(1+\bar r)` follows exactly from
\eqref{eq:nw-lom} with the contemporaneous timing, so the timing is deliberate and the
steady state is consistent. But out of steady state it is a within-period fixed point
interacting with an occasionally binding constraint, and the section should say so and
say how it is resolved in the solution.

## m8. The government block is written symmetrically for a country that never defaults

\eqref{eq:govt-budget}, \eqref{eq:debt-lom} and \eqref{eq:bohn} are all indexed
`X\in\{D,F\}` and all carry an **unsubscripted** `h_t` — the D-specific survival
factor, defined at l. 150 with `\varrho_D` and stated at l. 374 to be identically one
for F ("The `F`-sovereign never defaults"). As written the F government's budget
constraint applies a Greek haircut to German debt. Subscript it, or state `h^F_t\equiv1`.

Same three equations use `\delta_b` unsubscripted while l. 345 defines
`\delta_{X,b}`; \eqref{eq:bond-payoff} uses `\delta_b` for both countries. Section 3
gives one value, `\delta_b = 0.056`, for both. Pick one convention.

And the whole block is uncalibrated: `G_X` (government spending), `\bar B_X`
(steady-state debt stock), `\gamma_\tau` (the Bohn elasticity) and `\bar T^\tau_X`
appear in no table and no paragraph. Section 3's "Government" paragraph gives only
`\delta_b`. The introduction opens on a debt ratio going from 105 to 175 per cent of
GDP; `\bar B_D` is the parameter that would carry that, and it has no value. Section 3
also lists `B_D` as a state variable but not `B_F`, so German debt is apparently
fixed — never stated.

Finally, note for Section 4: `T^{\tau}_{X,t}` is **lump-sum** in the household budget
\eqref{eq:hh-budget}. With heterogeneous, borrowing-constrained households, a lump-sum
tax is mechanically regressive in consumption terms. Since the paper's stated question
includes "how the resulting exposure and the resulting gains fall … across households"
(l. 155–158), the incidence result is going to be largely a property of the financing
rule. That needs to be confronted rather than discovered.

## m9. Home bias is country-specific in \eqref{eq:ces-aggregator}, common in \eqref{eq:ces-D}–\eqref{eq:ces-F}, and country-specific again in Section 3

\eqref{eq:ces-aggregator} has `\varpi_X` and `\eta_X` — and, in the same equation,
`\varpi_X^{1/\eta}` with an **unsubscripted** `\eta` in the exponent next to
`(1-\varpi_X)^{1/\eta_X}` with a subscripted one. \eqref{eq:ces-D} and
\eqref{eq:ces-F} then use bare `\varpi` and `\eta` for both countries, which is a
*symmetry assumption*: the F price index as written,
`P^c_F=[\varpi+(1-\varpi)p_t^{\eta-1}]^{1/(1-\eta)}`, is correct only if F's home bias
equals D's. Section 3 (l. 14) then says "a home bias `\varpi_X` calibrated to match
the import penetration", country by country. Greek and German import penetration are
not equal, so the two are inconsistent.

The demand system itself is right: I verified `IM_{D,t}=(1-\varpi)(P^c_D/p_t)^\eta C_D`
and `IM_{F,t}=(1-\varpi)(P^c_F p_t)^\eta C_F` are the correct CES demands given the
`1/\eta` weight convention and the `1/p_t` price of the D-good in F. Only the
subscripting is wrong.

Also: \eqref{eq:ces-aggregator} is written for `C_{D,t}` over `c_{DD,t}`, `c_{FD,t}` —
aggregate on the left, lower-case on the right, no `i`. Homotheticity makes this
harmless but it should be stated once.

## m10. Two different duration formulas

Section 2, l. 345: "Duration is therefore approximately `1/\delta_{X,b}` quarters" —
at `\delta_b = 0.056` that is **17.9 quarters** (4.46 years). Section 3, l. 70: "its
expected duration is `(1+\bar r)/(\bar r+\delta_b)` quarters. Setting
`\delta_b = 0.056` delivers a duration of 16.9 quarters, or 4.2 years" — I get 17.0,
so the number is right modulo rounding, but it is not Section 2's formula. Since
`\delta_b` is calibrated *to* a measured Greek duration, the two formulas imply
different `\delta_b`. Use the exact one in both places. Flagged in report 3; still
present.

## m11. The steady state is not deterministic, and Section 2 offers no expression for the sovereign spread

Section 3 (l. 63) refers to "the deterministic steady state used for calibration". It
is not deterministic: \eqref{eq:pd} sets `\pi^d(\bar s)=0.1\%` per quarter at the rest
point, and by \eqref{eq:expectation} that probability is priced in
\eqref{eq:bond-pricing}. The steady state is a risky one.

This has a consequence Section 3 then trips over. It calls `\bar\varsigma` "the spread
of the sovereign yield over the risk-free rate" and calibrates it to 8 bp p.a. But by
\eqref{eq:bond-pricing} the steady-state sovereign spread has **two** components: the
constraint term `\lambda\bar\mu/\bar\Omega`, which I confirm is exactly 8.00 bp p.a.
at the reported calibration, **plus** the expected-loss term
`\approx\pi^d(1-\varrho_D)`, which at `\pi^d=0.1\%` per quarter is 21–26 bp p.a. for a
recovery rate in the 0.35–0.465 range. Total: roughly **30–34 bp p.a.**, against the
50 bp measured 1999–2007. That is a defensible number and arguably a point in the
paper's favour — but Section 3 does not report it, because Section 2 never writes down
an expression for the sovereign spread. It should: two lines after
\eqref{eq:bond-pricing}, separating the expected-loss and constraint components, and
naming the object the calibration targets.

(`\bar\varsigma` itself is never defined in Section 2. It appears for the first time
in Section 3. It *is* consistent with \eqref{eq:capital-euler}, which at a
deterministic point gives `\bar r^k-\bar r=\lambda\bar\mu/\bar\Omega` — but Section 2
should name it.)

## m12. The reported impulse is a 4.8-standard-deviation innovation, unremarked

Section 3 (l. 67): "The impulse we report raises `\pi^d` on impact from 0.10% to 1.98%
per quarter." Inverting \eqref{eq:pd}, that is a move in the latent factor from
`\bar s=-6.9068` to `-3.9021`, a jump of **3.005**, against `\sigma_s=0.63`. That is a
**4.77-standard-deviation innovation** — a one-in-a-million draw under the process
\eqref{eq:s-process} the paper has just calibrated. (It is 1.49 unconditional standard
deviations, and it sits inside Section 3's `\pm4.35` Smolyak box, so the solution is
not being extrapolated. The issue is the characterisation, not the numerics.)

Either the impulse is a sequence of innovations rather than one, or it is a large
deliberate scenario, or `\sigma_s` is too small. The paper should say which, and
report the shock in standard-deviation units as well as in probability units.

Related and unfixable by hand-waving: **no parameter of the production or household
block is given a value anywhere in the paper.** `\alpha_X`, `\delta_X`, `\xi_X`,
`\nu_X`, `\sigma_X`, `\epsilon_X`, `\chi_X`, `\beta_X`, `n_e`, `\varpi_X`,
`\varrho_D`, `\gamma_\tau`, `G_X`, `\bar B_X` — all absent. Section 3 l. 12 says "All
specific calibration values, as well as empirical targets, can be found in the
appendix [ADD APPENDIX SECTION]", and there is no appendix. `\varrho_D` in particular
is the recovery rate, it is the one parameter the introduction says Greece uniquely
identifies (l. 162–164), and Section 3 says only "We accordingly set the recovery rate
to match that event" without giving a number.

---

# Notation ledger

Report 3 recorded `\mu`/Smolyak-`\mu`, `\lambda_X`/`\Lambda`, `\omega^f_X`/`\omega_X`
and `f_X`/`f`. All four are still live. Adding the collisions internal to Section 2,
worst first:

| Symbol | Use A | Use B | Severity |
|---|---|---|---|
| `n_{X,t}` | **bank net worth**, \eqref{eq:ng-excess}–\eqref{eq:nw-lom} | **aggregate hours** inside the GHH composite `x_{X,t}\equiv c_{X,t}-v(n_{X,t})`, l. 401 | **severe** — same subscripts, same section, both aggregate |
| `a_{it}` | **household deposits**, \eqref{eq:hh-budget} | `a_{j,t}` = **bank asset position** `j`, \eqref{eq:ng-excess} | severe |
| `Q_{X,t}` / `q^{X}_t` | Tobin's q | bond price | severe — both appear in \eqref{eq:balance-sheet} |
| `\mu_{X,t}` | KKT multiplier | `\mu=1`, Smolyak level (§3) | high |
| `\lambda_X` | divertible fraction | `\Lambda_{X,t,t+1}`, the SDF | high |
| `\mathcal{D}` | `1-(1-f)(1+\bar r)` (§3) | `D`, the country index | high |
| `\Xi^{X}_{t+1}` | bond payoff | `\xi_X`, adjustment-cost elasticity | moderate |
| `\varrho_D` | recovery rate | `\rho_s`, risk persistence | moderate |
| `\varpi` | home bias | `\omega_X`, `\omega^f_X`, entrant transfer | moderate |
| `\sigma_X` | risk aversion | `\sigma_s`, innovation s.d. | moderate |
| `\epsilon_X` | retail elasticity | `\varepsilon_{t+1}` shock, `\varepsilon_m` quadrature node | low |
| `\delta_X` | depreciation | `\delta_b` / `\delta_{X,b}`, coupon decay | low |
| `b^{X}_{Y,t+1}` / `B_{X,t}` | bank holdings / total stock | convention never stated in words | **and destroyed by m1** |

Two more: `V_t(n)=\varphi_{X,t}n_{X,t}` (l. 204) writes the argument as `n` on the left
and `n_{X,t}` on the right; and `\mathrm{Cov}_t` is used at ll. 423, 425, 431 where
`econpaper.sty` already defines `\Cov`.

---

# Language, references and mechanics

**Register.** "let me describe the behaviour" (l. 323) and "let us think about"
(l. 401) in a two-author paper, alongside "we" elsewhere. Pick one; in a submitted
draft it should be "we".

**Broken or garbled sentences.**
- l. 79: "Thus the static problem is of the firm:"
- l. 165: "Because the working-capital loan is extended at rates the period at the rate locked at $t$" — unparseable.
- l. 204: "Following a standard result from \cite{gertler2011model}, we obtain that \eqref{eq:ng-excess} is linear in $n$ and the constraint below is linear in assets, the value function is linear" — needs "because … , the value function is linear".
- l. 288: "requirment"; and the sentence attributes a fall in `n_{X,t}` to \eqref{eq:ng}, which delivers `n^g_{X,t+1}` — the mechanism runs \eqref{eq:ng} → \eqref{eq:nw-lom}.
- l. 310: "This is because even a default-free bond trades below…" — the "because" does not connect to the preceding sentence.
- l. 459: "Taxes follow a modified \cite{bohn1998behavior} as a constant elasticity" — missing a noun ("rule"), and the citation is being used as one.

**Citation commands.** `\cite` is used where `\citet` or `\citep` is needed at
ll. 105, 204, 256, 345 (twice) and 459. Line 345 currently renders as
"issue (Hatchondo and Martinez, 2009)/(Chatterjee and Eyigungor, 2011) perpetuities".
Only l. 79's `\citep{neumeyer2005business}` is right. All seven keys resolve in
`references.bib`; the "undefined citation" warnings in `main.log` are first-pass
artefacts, not real.

**Cross-references.** Three `XXX` placeholders survive: `\cref{XXX}` (l. 35, with the
author's own reminder comment "% <-- point at your Section 1.2"), `Section~\ref{XXX}`
(l. 150), `Section~\ref{sec:XXX}` (l. 158). All three point at things the reader needs:
the predetermined-rate convention, the survival factor, and the bond payoff.

**LaTeX mechanics.**
- Stray `\\` after `\end{proof}` (l. 287). Also `\\`+`\\` used as paragraph breaks at
  ll. 178–179, 194–195, 362–363, 411–412 — use a blank line.
- l. 69: `\label{eq:mc}` sits outside any numbered environment, so it silently attaches
  to the section counter. Nothing references it yet; it will misfire when something
  does.
- l. 12: `\label{sec:Households}` capitalised against `sec:wc`, `sec:model`,
  `sec:policy`.
- `\subparagraph` (ll. 291, 299, 311) inside `\paragraph` inside `\subsection` is four
  levels deep and renders as an unnumbered run-in; consider a list.
- l. 38, l. 178: `taxes.\\` and `separately.\\` end paragraphs with `\\`, which
  produces an underfull hbox warning each time.

---

# What is right, and should not be touched

Worth recording, because the list of findings is long and the section is not bad:

- The CES demand system \eqref{eq:ces-D}–\eqref{eq:ces-F} is correct, including the
  `p_t^{\eta-1}` in the F price index and the `(P^c_F p_t)^\eta` in F's import demand.
- The Jermann block is exactly right: `Q=\iota^{\xi}/(\gamma_0(1-\xi))`, the inverted
  investment rate, and `\Pi^k_{X,t}` as the maximised objective. Section 3's
  `\gamma_{0,X}=\delta^{\xi}/(1-\xi)`, `\gamma_{1,X}=-\delta\xi/(1-\xi)` are exactly
  the coefficients that deliver `\bar Q=1` and `\bar\iota=\delta` from
  \eqref{eq:capital-lom}. I verified both.
- The excess-return representation \eqref{eq:ng-excess} follows correctly from
  \eqref{eq:balance-sheet} and \eqref{eq:ng}, and the `P_{X,t+1}` device genuinely is
  a sufficient statistic for the liability side. It is a nice piece of construction and
  the section is right to say so.
- Proposition 1 and its proof are correct, including the slack case and the
  `\max\{\cdot,0\}` construction.
- \eqref{eq:bond-pricing} and \eqref{eq:rwc} follow correctly from
  \eqref{eq:foc-general}.
- The logistic calibration is right: `\bar s=\log(0.001/0.999)` gives
  `\pi^d(\bar s)` = 0.1000% to machine precision.
- \eqref{eq:govt-budget}, \eqref{eq:debt-lom} and the anchor
  `\bar T^\tau_X=G_X+\delta_b\bar B_X(1-\bar q^X)` are mutually consistent.
- Section 3's entire GK steady state reproduces Section 2's equations. I get
  `\bar\mu=0.000996` (reported 0.0010), `\bar\varsigma\bar\theta=0.1000\%` per quarter
  (reported 0.10%), the `\omega_X>0` bound at 3.8667% (reported 3.87%),
  `\bar\varphi=1.0255`, `\lambda_X=0.2051`. The triangularity argument at l. 46–48 is
  correct and the identification point — that `\bar\mu` depends only on the product
  `\bar\varsigma\bar\theta`, so leverage and the wedge are not separately identified —
  is exactly the right thing to have noticed.

---

# Priority

1. **m1** — define `\D` and `\F`. Two minutes, and \eqref{eq:balance-sheet} is
   currently corrupt in the compiled paper.
2. **M3** — delete or repair the clause at l. 411. One sentence, and it currently
   inverts the paper's own mechanism.
3. **M1** — write the equilibrium. Half a page, from the code.
4. **M2** — decide how the cross-border position is pinned down and say so. This one
   is a modelling decision, not an edit, and per **M9** it is the equation that answers
   the paper's own second research question.
5. **M5**, **M6** — reconcile the union/parity claim and the working-capital timing.
   **M5** is also a precondition for **M9**.
6. **M4** — decide what the household heterogeneity is for, and make the title and
   Section 2 agree with Section 3's state vector.
7. **M7**, **M8** — present the decomposition as an approximation, and state the
   transmission elasticity.
8. **m4** — remove or complete Proposition 2 before anyone external sees the draft.
9. Everything else.
