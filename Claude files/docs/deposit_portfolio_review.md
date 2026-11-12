# Deposit funding and cross-border sovereign portfolios

**Diagnosis, tractable alternatives, and a stress test for the global model**

*Working note, 12 October 2026, branch `analysis` (from `main` at `55e174d`). Markdown companion to
[deposit_portfolio_review.tex](deposit_portfolio_review.tex) and its PDF. Every stress-test number below is
from the global-projection solver on the quick grid (25 points) with the union deposit market, unless a
line says otherwise.*

---

## Contents

1. [Summary](#1-summary)
2. [Issue 1: the deposit funding rate](#2-issue-1-the-deposit-funding-rate)
3. [Issue 2: cross-border holdings](#3-issue-2-cross-border-holdings)
4. [Literature](#4-literature)
5. [Alternatives for the funding rate](#5-alternatives-for-the-funding-rate)
6. [Alternatives for cross-border holdings](#6-alternatives-for-cross-border-holdings)
7. [Stress test: design](#7-stress-test-design)
8. [Stress test: results](#8-stress-test-results)
9. [Next steps](#9-next-steps)
10. [Files](#10-files)
11. [References](#11-references)

---

## 1. Summary

| | Finding | Recommendation |
|---|---|---|
| **Issue 1: deposit funding rate** | The deposit rate is a market-clearing price that **falls** when sovereign risk rises. On impact it cancels **51–61%** of the credit-spread rise. The deposit-market choice flips the sign of 25-quarter cumulative output (−3.9 national against +0.5 union). | **D3:** fund the marginal unit of bank liabilities through a Eurosystem refinancing facility at a policy rate, against sovereign collateral, with a penalty schedule. |
| **Issue 2: cross-border holdings** | The portfolio cost ψ (the cross-border adjustment-cost coefficient) does four jobs at once, its anchor is a level rather than a share, and the D bank's calibrated foreign share is 20% against 0.75% in the data. | **P1 now:** anchor on shares with a measured return elasticity. **P3 (fixed shares) later**, once D3 carries the external margin. |
| **Stress test** | Under union funding, per log unit of the dial, κ_nfa (the union closure's debt-elastic premium) moves impact output **4×** as much as ψ. ψ mainly moves *who holds* the bonds. No solution exists at ψ ≤ 0.5. | Run the same sweep under national funding, where ψ is predicted to matter far more; confirm on the full grid. |

---

## 2. Issue 1: the deposit funding rate

### 2.1 What the code does now

Two closures exist. R_X (gross deposit rate of country X, set at t, paid at t+1) is pinned by:

- **National market** (default on the unmerged branch `domestic-funding-fix`): each bank is funded by its own households only,

$$\mathrm{dep}_{X,t} = P^c_{X,t}\,HA_{X,t}, \qquad X\in\{D,F\},$$

  where dep is bank deposits, P^c the consumption price index and HA household deposits in basket units.

- **Union market** (`main`, and the TPI's economy on every branch): one union-wide clearing plus real interest
  parity with a debt-elastic premium,

$$R_{D,t} = R_{F,t}\,\frac{\mathbb{E}_t[p_{t+1}]}{p_t} - \kappa_{\rm nfa}\,\frac{V^{\rm flow}_{D,t}}{\bar Y_D},$$

  where p is the price of the F good in D goods, V^flow (D's net cross-border deposit position) is carried as a
  state, and κ_nfa = 0.01 (`point_map.py:451`, `calibration/global_projection.py:46`).

Output responds to the financial block only through the loan rate r^wc (the rate on working-capital loans):

$$r^{wc}_{X,t} = r_{X,t} + \underbrace{\frac{\lambda_X\,\mu_{X,t}}{\mathbb{E}_t[\Omega_{X,t,t+1}]}}_{\text{credit spread}},$$

where λ is the divertible fraction of bank assets, μ the multiplier on the incentive constraint and Ω the
banker's discount factor. A fall in r offsets the spread one for one.

### 2.2 Fault 1: the funding rate is an automatic stabiliser with the wrong sign

Headline shock: the quarterly default probability p^d rises to 2%. Impact response, **full grid**:

| | union, `main` | union, current | national, current |
|---|---:|---:|---:|
| credit spread (bp/yr) | +79.4 | +76.6 | +65.0 |
| deposit rate r_D (bp/yr) | −47.9 | −46.6 | −33.1 |
| loan rate r^wc_D (bp/yr) | +34.5 | +32.9 | +34.1 |
| **share of the spread cancelled** | **60%** | **61%** | **51%** |
| output: credit-spread channel | −0.179% | −0.173% | −0.146% |
| output: deposit-rate channel | +0.104% | +0.101% | +0.072% |
| output: total impact | −0.073% | −0.061% | −0.068% |

*"Current" means `domestic-funding-fix`, with the household-budget and cross-border p′/p fixes. Sources: the run logs and
output decompositions in `results/GLOBAL`, `results/UNION_current` and `results/NATIONAL`.*

- **The spread channel alone matches Bocola.** On its own it gives −0.15% to −0.18%, which is Bocola's
  like-for-like open-economy target (−0.157%). In his economy the funding rate is a world rate,
  R = 1/β + 0.01·B_for/gdp, which hardly moves on impact. So, to first order, **the gap between our response and
  his is the endogenous funding rate.** (This is an accounting decomposition, not a counterfactual.)
- **The sign is counterfactual.** In 2010–12, periphery bank funding costs *rose*: Greek deposit rates rose
  relative to the euro area, deposits ran off, and the Eurosystem became the marginal funder. In the model,
  deposits are riskless in both regimes, so nothing ties the deposit rate to the sovereign.

### 2.3 Fault 2: each closure buys solvability with a device

| National | Union |
|---|---|
| Needs three solver devices: no-short-sale pairs on all three bank sovereign books, a collapsed TPI band, and a least-squares Newton fallback. | Needs an extra state V whose only anchor is κ_nfa, a stationarity device borrowed from the small-open-economy literature and never calibrated here. |
| A TPI purchase barely moves the D bond price, so the cap cannot be held below about 490 bp. **The backstop is therefore evaluated on a different economy (union) from the baseline (national).** | Puts an expectation, E[p′], inside a residual. |

### 2.4 Fault 3: the closure, not the friction, decides the medium run

| Current code, headline shock | national | union |
|---|---:|---:|
| impact output | −0.068% | −0.061% |
| trough output | −0.211% (q7) | −0.061% (q0) |
| cumulative output, 25 quarters | **−3.89** | **+0.51** |
| cumulative investment, 25 quarters | −6.00 | +4.07 |

The impact responses are almost the same, but the medium run has opposite signs. A reader cannot tell from the paper which one is the
model's answer.

---

## 3. Issue 2: cross-border holdings

### 3.1 What the code does now

The F bank's first-order condition for the D bond carries a non-pecuniary cost (`point_map.py:356–357, 465`):

$$\mathbb{E}_t\!\left[\Omega_{F}\,\Xi^D_{t+1}\,\frac{p_t}{p_{t+1}}\right]
= q^D_t\left(\mathbb{E}_t[\Omega_F]R_{F,t} + \lambda_F\mu_{F,t}\right)
\left[1 + \psi_F\,\frac{b^D_{F,t+1} - \bar b^D_F}{\bar B_D}\right].$$

Here Ξ^D is the bond's payoff (coupon plus continuation value, times recovery in default), q^D its price, b^D_F the
F bank's holding, b̄^D_F = 0.196 its anchor (per D capita), B̄_D the steady-state debt stock, and ψ_F = ψ_D = 2.
The p_t/p_{t+1} factor is in the paper draft and on `domestic-funding-fix`, **not on `main`**.

### 3.2 What ψ is

Define e (the F bank's expected excess return on the D bond over its marginal funding cost, constraint
included):

$$e^D_{F,t} \equiv \frac{\mathbb{E}_t[\Omega_F\,\Xi^D_{t+1}\,p_t/p_{t+1}]}{q^D_t(\mathbb{E}_t[\Omega_F]R_{F,t} + \lambda_F\mu_{F,t})} - 1
\qquad\Longrightarrow\qquad
\frac{b^D_{F,t+1} - \bar b^D_F}{\bar B_D} = \frac{e^D_{F,t}}{\psi_F}\ \ \text{(exactly)}.$$

**1/ψ is the semi-elasticity of cross-border holdings with respect to the quarterly excess return.** At ψ = 2,
an excess return of 100 bp/yr moves holdings by 0.125% of the stock.

### 3.3 Four jobs for one parameter

| Job | What ψ does |
|---|---|
| (a) Selection | Single-λ bankers' value is linear in their portfolio, so the split of each bond between the two banks is indeterminate at the steady state and a corner off it. ψ > 0 makes it determinate. |
| (b) Elasticity | By the identity above, ψ alone scales the retrenchment or flight-home response. |
| (c) External margin | Under national funding, D's external position is the net cross-border sovereign book, p·q^F·b^F_D − q^D·b^D_F, where q^F is the F bond price and b^F_D the D bank's F-bond holding. ψ is then the inverse elasticity of D's access to foreign saving: the national "S = I trap". |
| (d) Solvability | In the default branch, at the 55% haircut, F banks want to sell their whole D book. A finite ψ keeps that branch interior: **no solution at ψ ≤ 0.5** (Section 8). |

### 3.4 Two calibration faults

**(e) The anchor is a level.** When the risk shock raises D's debt, ψ makes F banks resist holding more than
b̄, so new issuance is absorbed at home. That is home bias produced by the units, not by a mechanism. The sweep
confirms it: a *larger* ψ gives a *larger* fall in the F share.

**(f) The steady-state position contradicts the data on the D side.** Symmetric gross positions, needed for a
zero net position at p = 1, give:

| | model SS | EBA 2011 stress test (Dec 2010) |
|---|---:|---:|
| F-bank share of bank-held D debt | 20.0% | 12.7% |
| foreign share of the **D** bank's sovereign book | **20.0%** | **0.75%** |
| foreign share of the F bank's sovereign book | 2.5% | 2.45% |

*EBA, gross direct long: Greek banks hold €54.4bn of Greek and €0.41bn of German sovereigns; German banks hold €7.93bn
of Greek and €315.3bn of German (`data/eba_moments.json`).*

The D bank therefore carries a buffer of foreign safe bonds about **27 times** the one Greek banks actually had.
Under national funding, that buffer is D's only way to finance a current-account gap.

---

## 4. Literature

**Funding rates in intermediary models**

- *Riskless deposits from a household Euler equation:* Gertler & Karadi (2011) and Gertler & Kiyotaki (2010).
  The latter add an interbank market whose friction spans segmented and frictionless funding.
- *Bocola (2016):* deposits are riskless; his open-economy section uses a world rate plus a Schmitt-Grohé &
  Uribe (2003) debt-elastic premium.
- *Domestic plus foreign funding:* Aoki, Benigno & Kiyotaki (2016) price the funding mix with a convex
  divertibility of foreign borrowing. Two-country intermediary models with cross-border funding:
  Akinci & Queralto (2018), Dedola, Karadi & Lombardo (2013), Kollmann, Enders & Müller (2011).
- *Country-spread closures:* Neumeyer & Perri (2005) and Uribe & Yue (2006). In Neumeyer & Perri the spread
  works through working capital, which is this model's own output channel.
- *Sovereign-to-private pass-through:* Corsetti, Kuester, Meier & Müller (2013).
- *Sticky deposit rates:* Drechsler, Savov & Schnabl (2017) and Begenau (2020).
- *Euro-area plumbing:* Drechsler, Drechsel, Marques-Ibanez & Schnabl (2016) on lender-of-last-resort
  borrowing against risky collateral; Crosignani, Faria-e-Castro & Fonseca (2020) on LTRO funding used to buy
  domestic sovereigns.
- *Retail versus wholesale funding, with runs:* Gertler, Kiyotaki & Prestipino (2016).

**Cross-border portfolios**

- *Portfolios from higher-order terms:* Devereux & Sutherland (2011) and Tille & van Wincoop (2010). Not needed
  in a global solution, but they show why a device is needed at first order.
- *Finite-elastic holdings tied to risk capacity:* Blanchard, Giavazzi & Sa (2005), Gabaix & Maggiori (2015),
  Itskhoki & Mukhin (2021). Survey of home bias: Coeurdacier & Rey (2013).
- *Measuring the elasticity directly:* demand-system estimates (Koijen & Yogo 2019; Koijen, Koulischer,
  Nguyen & Yogo 2021) estimate exactly the object 1/ψ for euro-area sovereigns, by holder sector.
- *Euro-crisis repatriation:* documented by Broner, Erce, Martin & Ventura (2014) and Brutti & Sauré (2016).
  Proposed explanations:
  - creditor discrimination (Broner et al. 2014);
  - risk-shifting and carry trades (Acharya & Steffen 2015);
  - moral suasion (Ongena, Popov & Van Horen 2019);
  - financial repression (Becker & Ivashina 2018);
  - the doom loop (Farhi & Tirole 2018; Gennaioli, Martin & Rossi 2014).
- *Exposures and transmission:* Altavilla, Pagano & Simonelli (2017); Battistini, Pagano & Simonelli (2014).
- *A structural basis for risk-sensitive exposure limits:* value-at-risk constraints (Adrian & Shin 2014).

---

## 5. Alternatives for the funding rate

### D1. Two-tier funding with costly cross-border wholesale funding

The D bank funds itself with net worth, national retail deposits at R_D, and union wholesale funding w at R^w
(the union money-market rate). Wholesale funding raises divertibility, where x = w/𝒜 is the wholesale share of
assets 𝒜:

$$\alpha_t n_t \;\ge\; \lambda\,\Theta(x_t)\,\mathcal{A}_t, \qquad \Theta(x) = 1 + \tfrac{\gamma}{2}x^2 .$$

Here α is the franchise value per unit of net worth n. Because Θ(x)𝒜 is homogeneous of degree one, Bocola's
closed form survives with λ → λΘ(x):

$$\mu_t = \max\!\left\{1 - \frac{\mathbb{E}_t[\Omega]\,R_{D,t}\,n_t}{\lambda\,\Theta(x_t)\,\mathcal{A}_t},\,0\right\},
\qquad
R_{D,t} - R^w_t = \frac{\lambda\gamma\,\mu_t\,x_t}{\mathbb{E}_t[\Omega]} + \kappa_w\,(x_t - \bar x).$$

Every asset first-order condition carries λμ(1 − γx²/2). The small κ_w term (a determinacy device) matters only
where μ = 0.

- **Nests:** γ → ∞ gives national; γ = 0 gives union.
- **Sign:** when μ_D rises, the domestic deposit rate is pushed *up*, the deposit competition of 2010–12.
- **Closed-form calibration:** at the steady state, γ ≈ (R_D − R^w) / [x̄ · (r^wc − r_D)], from the retail
  deposit spread over Euribor, the lending spread and the wholesale share x̄.
- **Cost:** one more unknown and residual (15 in total); V is reused as the wholesale position.

### D2. Defaultable deposits: sovereign risk in the deposit rate

Deposits pay a fraction h^dep (the deposit recovery rate) in the default regime: bail-in, as in Cyprus in 2013,
or redenomination. With x_D the GHH consumption composite, σ the curvature of utility and d_{t+1} the default
indicator, the household Euler equation weights the default branch:

$$x_{D,t}^{-\sigma} = \beta\,\mathbb{E}_t\!\left[\big(1-(1-h^{\rm dep})\,d_{t+1}\big)\,R_{D,t}\,\frac{P^c_{D,t}}{P^c_{D,t+1}}\,x_{D,t+1}^{-\sigma}\right].$$

The bank's deposit leg is written down by the same factor in that branch, so the depositor's loss is the bank's
gain and the budget does not leak. To first order:

$$r_{D,t} - r^{\,h=1}_{D,t} \;\approx\; \pi_t\,(1-h^{\rm dep})\,\frac{\mathbb{E}_t[x^{-\sigma}_{t+1}\mid d'=1]}{\mathbb{E}_t[x^{-\sigma}_{t+1}]}.$$

- **Calibration:** matching the Corsetti et al. reduced form (deposit spread = ξ × sovereign spread) gives
  h^dep ≈ 1 − 0.55·ξ.
- **Cost:** no new state or unknown. The default branch is already in the quadrature; only the d′ = 1
  continuation state becomes regime-specific.
- **Caveat:** the bail-in transfer cushions n_D in default, which lowers the bank-side risk premium.
- Works with any deposit market.

### D3. Eurosystem refinancing at a policy rate, against sovereign collateral (recommended)

The D bank's marginal funding is a Eurosystem facility ℓ (in D goods). The policy rate R^ℓ rises with usage
(κ_ℓ is the slope of the penalty schedule, as ELA was priced above the MRO rate), and borrowing is limited by
haircut-adjusted collateral (h^c is the haircut):

$$\mathcal{A}_{D,t} = n_{D,t} + \mathrm{dep}_{D,t} + \ell_t, \qquad
R^\ell_t = \bar R + \kappa_\ell\,\frac{\ell_t - \bar\ell}{\bar Y_D}, \qquad
\ell_t \le (1-h^c)\,q^D_t\,b^D_{D,t+1}.$$

$$R_{D,t} = R^\ell_t + \frac{\eta_t}{\mathbb{E}_t[\Omega_D]},
\qquad \eta_t \ge 0,\ \ (1-h^c)q^D_t b^D_{D,t+1} - \ell_t \ge 0,\ \ \text{complementary}.$$

$$q^D_t\big(\mathbb{E}_t[\Omega_D]R_{D,t} + \lambda\mu_{D,t} - (1-h^c)\,\eta_t\big) = \mathbb{E}_t[\Omega_D\,\Xi^D_{t+1}].$$

Here η is the collateral multiplier. The Eurosystem funds ℓ with reserves held by F banks and remits the
facility's profit and loss, including the terms-of-trade valuation, by the capital key. It uses the accounts the TPI already has, so
a rate fixed in D goods is consistent in a two-country model.

- **What it does:** while collateral is slack, R_D = R^ℓ is set by policy. The deposit rate stops absorbing the
  shock, and D draws on union saving through the facility (TARGET2 in all but name). When the sovereign is marked
  down, collateral binds, η > 0 and R_D **rises**: the pass-through has the right sign. Binding collateral also
  raises the D bank's valuation of its own bond, so repatriation becomes an outcome.
- **Nests:** with slack collateral and R^ℓ = R̄ + κ_ℓ·ℓ/Ȳ it is Bocola's open-economy closure (κ_ℓ = 0.01);
  κ_ℓ → ∞ with ℓ̄ = 0 is the national market.
- **Calibration:** h^c from the published Eurosystem haircut schedule; ℓ̄ from pre-crisis Eurosystem borrowing
  by Greek banks; κ_ℓ from the ELA–MRO spread.
- **Cost:** one state (replacing the inert V under national funding), one unknown, and one complementarity
  pair. The national solver devices, the collapsed TPI band and the separate union baseline should no longer be
  needed.

### Comparison

| | Sign of Δr_D when risk rises | Nests | New parameters (data source) | States / unknowns |
|---|---|---|---|---|
| national (now) | falls | — | — | 12 / 14, plus 3 devices |
| union (now) | falls | — | κ_nfa (none) | 12 / 14 |
| D1 wholesale | rises if μ rises | national, union | γ, x̄ (ECB MIR, Bank of Greece) | 12 / 15 |
| D2 defaultable | rises with π | any closure | h^dep (deposit spread) | unchanged |
| **D3 Eurosystem** | rises if collateral binds | national, Bocola §V.C | h^c, ℓ̄, κ_ℓ (ECB) | 12 / 15, plus 1 complementarity pair |

**Recommendation: D3.** It is the only option in which the marginal funding rate is a policy rate rather than
an equilibrium price, and it puts the baseline and the TPI on one funding structure. D2 can be layered on top
for deposit flight; D1 is the cleanest robustness dial.

---

## 6. Alternatives for cross-border holdings

### P1. Share anchor with a measured return elasticity (recommended now)

Replace the level anchor with s_F (the F bank's share of the outstanding D stock):

$$s_{F,t} \equiv \frac{b^D_{F,t+1}}{B^D_{t+1}}, \qquad
1 + \tilde\psi_F\,(s_{F,t} - \bar s_F)\ \text{replaces the level term}
\quad\Longleftrightarrow\quad
s_{F,t} - \bar s_F = \frac{e^D_{F,t}}{\tilde\psi_F}.$$

- New issuance is shared at s̄_F unless returns say otherwise.
- **Calibration:** 1/ψ̃ comes from demand-system estimates for foreign banks, or from the 2010–11 change in
  German banks' Greek exposure (EBA 2011 stress test against the 2011 capital exercise; BIS consolidated banking
  statistics) divided by the model-implied change in e. Set s̄_F = 12.7% (EBA).
- **Note:** matching the D side (0.75%) breaks the zero net position at p = 1; under D3 the facility can carry
  D's steady-state net debt instead.
- **Cost:** none; same residuals and states.

### P2. Limited risk-bearing capacity (a risk budget)

GK kernels treat D's default as nearly idiosyncratic for an F banker, so nothing makes foreign banks *shed*
when risk rises. A risk budget adds that force. Let χ_F be the exposure as a share of F net worth, and 𝒱 its
conditional return variance, computed on the existing quadrature:

$$\chi_{F,t} \equiv \frac{q^D_t\,b^D_{F,t+1}}{\varkappa\,p_t\,n_{F,t}}, \qquad
\chi_{F,t}\,\big(\bar\sigma^2 + \mathcal{V}_t\big) = \chi_{\rm ss}\,\bar\sigma^2 + \frac{e^D_{F,t}}{\gamma_F},$$

with ϰ = 8 the size of F relative to D. It is determinate at the deterministic steady state (𝒱 = 0, e = 0).
When π rises, 𝒱 rises roughly as π·0.55², and the exposure falls by the factor σ̄²/(σ̄² + 𝒱). When F net worth
falls, the level falls too, which carries contagion. σ̄² is the variance at which the budget halves the
exposure; γ_F is the return sensitivity. **Cost:** the same residual slots and a few quadrature sums.

### P3. Fixed shares (segmentation, the ψ → ∞ limit)

b^D_F = s̄_F·B^D and b^F_D = s̄_D·B^F (B^F is fixed in the code).

- **Gain:** residuals 9 and 13 and their unknowns go, and the states b_DF and b_FD become functions of b_DD
  and M^cb (the central bank's D-bond book). The grid drops from 12 states to 10, and the full grid becomes
  about **30% cheaper** (115 → 95 points; the solve scales as n²).
- **What is kept and lost:** F banks still take marked-to-market losses, so contagion survives. Cross-bank
  rebalancing is lost, and with it the split of the TPI's purchases by seller.
- **Only admissible with an external funding margin** (union, D1 or D3). Under national funding it closes D
  to foreign saving altogether.

**The general point:** the model needs exactly one external adjustment margin with an elasticity that can be
calibrated. Today that margin is either a frictionless deposit position (union, governed by κ_nfa) or the
ψ-governed bond books (national). D3 puts it on the funding side, where 2010–12 actually put it, and leaves ψ
with one job instead of four.

---

## 7. Stress test: design

Implemented in [`psi_sweep.py`](../diagnostics/portfolio_stress/psi_sweep.py). Each point is a full global
re-solve, never a perturbation of fixed rules.

| Step | What |
|---|---|
| 1. Sweep | ψ ∈ {0.25, 0.5, 1, 2, 4, 8, 32}, with ψ_F and ψ_D moved together. The calibrated point is solved cold, then the sweep walks outward with each solve warm-started from its neighbour, retrying in 2 and then 4 log-steps before a cold restart. |
| 2. Read | Bocola's differenced IRF from **each solution's own** stochastic rest point; shock p^d = 2%/qtr, 25 quarters. |
| 3. Metrics | Impact, trough and cumulative output; cumulative investment; impact bank net worth n_D, n_F and bond price q^D; impact credit spread and deposit rate; F-bank share of D debt at q4; D's net external position at q8; the model's retrenchment elasticity. |
| 4. Sensitivity | Log-elasticity ∂ln\|m\|/∂ln ψ at the calibrated point, and the sweep range (max − min)/\|m(ψ_ref)\|. |
| 5. Isolation | (a) Portfolio channel: m(ψ) − m(ψ → ∞). (b) Within-period versus expectations: re-clear the period map at impact under each ψ with the ψ = 2 rules held fixed. |
| 6. Accuracy guards | Solver verdict and worst residual; largest box clip along the IRF (the bond band is ±12%); impact read both fitted and exactly cleared, since an effect smaller than that gap is not identified. |
| 7. Ranking | The same sweep for other dials (`--param`), ranked with `--compare`. A Morris screen over eight dials is about 90 quick solves (overnight). |
| 8. Funding interaction | The whole sweep in both deposit markets. Prediction: ψ governs the medium run under national funding and only the composition under union. |

```bash
S="Claude files/diagnostics/portfolio_stress/psi_sweep.py"
python3 "$S"                                   # psi, quick grid
python3 "$S" --values 0.5 8 --s-refine 5       # full-grid confirmation
python3 "$S" --param kappa_nfa --values 0.0025 0.04
python3 "$S" --set deposit_market=national     # domestic-funding-fix only
python3 "$S" --compare                         # rank all saved sweeps
```

---

## 8. Stress test: results

*Quick grid (25 points), union deposits, `main`. A preview, not a converged result: at ψ = 2 the full grid gives
an impact output response of −0.073%, against −0.064% here.*

### 8.1 The ψ sweep

| ψ | start | Y₀ fitted | Y₀ exact | Σ Y (25q) | Σ I (25q) | n_D,0 | spread₀ (bp/yr) | r_D,0 (bp/yr) | Δs_F q4 (pp) | nfa q8 (% Y) | Y₀ exact at ψ=2 rules |
|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | — | **no solution** | | | | | | | | | |
| 0.5 | — | **no solution** | | | | | | | | | |
| 1 | warm | −0.047% | −0.060% | +1.64 | +8.24 | −5.09% | +94.8 | −47.8 | −0.342 | +0.701 | −0.069% |
| **2** | cold | **−0.064%** | **−0.076%** | **+1.63** | **+8.49** | **−5.20%** | **+100.1** | **−50.6** | **−0.445** | **+0.756** | **−0.076%** |
| 4 | warm | −0.051% | −0.062% | +1.90 | +8.98 | −5.20% | +102.3 | −52.3 | −0.505 | +0.745 | −0.079% |
| 8 | warm | −0.047% | −0.056% | +2.03 | +9.19 | −5.22% | +103.7 | −53.2 | −0.537 | +0.744 | −0.081% |
| 32 | warm | −0.043% | −0.052% | +2.13 | +9.36 | −5.23% | +104.9 | −53.9 | −0.563 | +0.743 | −0.083% |
| **range / \|value at ψ=2\|** | | **33%** | | **30%** | 13% | 3% | 10% | 12% | **50%** | 7% | |

Column key:
- Y₀ is impact output of country D.
- Σ is the 25-quarter sum of percentage deviations.
- n_D,0 is impact D-bank net worth.
- Δs_F is the change in the F-bank share of D debt.
- nfa is D's net external position.
- The last column re-clears the period map under each ψ with the ψ = 2 continuation rules held fixed.

Every solved point passed the acceptance test (worst residual ≤ 5.2e-9), with no box clipping.

![ψ sweep: six IRF panels](../../results/STRESS/psi_bF_D_default_s0.png)

*(The figure is at `results/STRESS/psi_bF_D_default_s0.png`; `results/` is git-ignored, so the link works
only locally.)*

### 8.2 ψ against κ_nfa

| κ_nfa | Y₀ fitted | Y₀ exact | Σ Y | Σ I | n_D,0 | spread₀ | r_D,0 | Δs_F q4 | nfa q8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0025 | −0.085% | −0.096% | +2.02 | +8.98 | −5.32% | +102.5 | −52.1 | −0.441 | +0.783 |
| **0.01** | **−0.064%** | **−0.076%** | **+1.63** | **+8.49** | **−5.20%** | **+100.1** | **−50.6** | **−0.445** | **+0.756** |
| 0.04 | −0.014% | −0.025% | +1.10 | +7.62 | −4.90% | +95.0 | −47.4 | −0.453 | +0.712 |

| Dial (values swept) | ln range | Y₀ range | Σ Y range | Δs_F range | nfa range | Y₀ per ln unit | Σ Y per ln unit |
|---|---:|---:|---:|---:|---:|---:|---:|
| ψ (1 → 32) | 3.47 | 33% | 30% | 50% | 7% | 9.5% | 8.7% |
| κ_nfa (0.0025 → 0.04) | 2.77 | **110%** | **56%** | 3% | 9% | **40%** | **20%** |

### 8.3 Reading

1. **Solvability bounds ψ from below.** At ψ ≤ 0.5 the warm walk from ψ = 1 stalls at max|F| ≈ 3e-2, and the
   cold ladder hits a non-finite iterate in the default-branch haircut homotopy. In that branch the F bank wants
   to sell its whole D book, and the 1e-4 floors (used instead of a complementarity pair) make the Jacobian
   degenerate. This is the national branch's failure, reached through the portfolio cost.
2. **Under union funding, ψ governs composition, not the current account.** The F share moves by half its value
   across the sweep; D's net external position moves by only 7%. The deposit position V absorbs the external
   adjustment whatever ψ is.
3. **Output moves by a third over a 32-fold range, with a bump at the calibrated point.** A stiffer friction
   gives a smaller impact fall and a larger later boom. Above ψ = 2 both are monotone: the one-sided elasticity
   from 2 to 32 is −0.15 for |Y₀| and +0.10 for cumulative output. The ψ = 1 point breaks monotonicity on
   impact, and this is **not** multiplicity. Cold-solving ψ = 4 reproduces the warm solution exactly
   (−0.0514%), and walking 4 → 2 warm reproduces the cold ψ = 2 (−0.0644%). The full grid at {1, 2, 4} must
   decide whether the bump survives refinement.
4. **The mechanism runs through expectations.** With the rules fixed, a larger ψ *deepens* the impact fall
   (−0.069% → −0.083%): today, a stiffer friction stops F banks buying at the dip, so the D bank absorbs more.
   The fully re-solved effect has the opposite sign, so ψ acts mostly through expected future rebalancing.
5. **The level anchor shows up as predicted.** The F share falls *more* as ψ rises (−0.34 → −0.56 pp at q4),
   because the friction holds the level while debt grows. At ψ = 1 the F share rises on impact (F banks buy at
   the dip) before falling.
6. **The funding closure dominates.**
   - **Size of the gap:** per log unit, κ_nfa moves impact output **4×** as much as ψ and cumulative output
     **2×** as much. The most powerful dial of the two is a stationarity device that was never calibrated.
   - **Where it acts:** the effect runs through the terms of trade and labour supply, not the loan rate. Impact
     r^wc is +50 bp/yr at κ_nfa = 0.0025 and +48 at 0.04, and κ_nfa has almost no effect at fixed rules
     (−0.0756% → −0.0758%).
   - **Against the decomposition:** the deposit-rate channel there (+0.104 pp) is about five times the whole ψ
     range of impact output (0.021 pp). The two numbers come from different grids, so read this as an order of
     magnitude.
   - **Answer, for the union economy:** ψ strongly drives the cross-border composition and moderately drives
     the output path; Issue 1 is first-order. Under national funding the prediction is the reverse, and that
     sweep is the next run.

---

## 9. Next steps

1. **Port the cross-border p′/p fix to `main`** before quoting any ψ result; this sweep uses `main`'s condition.
2. **Run the sweep under national funding** on `domestic-funding-fix`
   (`--set deposit_market=national`), and the **full-grid confirmation** at ψ ∈ {1, 2, 4} (`--s-refine 5`).
3. **Pull four data moments:**
   - the Greece–euro area deposit spread (ECB MIR);
   - the non-resident wholesale share and Eurosystem borrowing of Greek banks (Bank of Greece, ECB);
   - German banks' Greek exposure in 2011 (EBA capital exercise, BIS consolidated banking statistics);
   - the ELA penalty spread.
4. **Implement D3** behind `deposit_market="eurosystem"`. The N1–N5 nesting tests extend directly
   (κ_ℓ → ∞ must reproduce national). Then implement P1, re-run the sweep, and decide on P3.

---

## 10. Files

| File | Contents |
|---|---|
| `Claude files/docs/deposit_portfolio_review.md` | this note |
| `Claude files/docs/deposit_portfolio_review.tex` / `.pdf` | the LaTeX version (13 pages, full bibliography) |
| `Claude files/docs/deposit_portfolio_results.tex` | the results section, `\input` by the .tex |
| `Claude files/diagnostics/portfolio_stress/psi_sweep.py` | the stress-test script |
| `results/STRESS/psi_bF_D_default_s0.{json,png}` | ψ sweep data and figure (git-ignored) |
| `results/STRESS/kappa_nfa_default_s0.{json,png}` | κ_nfa sweep data and figure (git-ignored) |
| `results/STRESS/psi_quick_sweep.log`, `kappa_quick.log` | console logs of both sweeps |

---

## 11. References

- Acharya, Drechsler & Schnabl (2014), "A pyrrhic victory? Bank bailouts and sovereign credit risk", *JF*.
- Acharya & Steffen (2015), "The 'greatest' carry trade ever?", *JFE*.
- Adrian & Shin (2014), "Procyclical leverage and value-at-risk", *RFS*.
- Akinci & Queralto (2018), "Exchange rate dynamics and monetary spillovers with imperfect financial markets", FRBNY Staff Report.
- Altavilla, Pagano & Simonelli (2017), "Bank exposures and sovereign stress transmission", *Review of Finance*.
- Aoki, Benigno & Kiyotaki (2016), "Monetary and financial policies in emerging markets", working paper.
- Battistini, Pagano & Simonelli (2014), "Systemic risk, sovereign yields and bank exposures in the euro crisis", *Economic Policy*.
- Becker & Ivashina (2018), "Financial repression in the European sovereign debt crisis", *Review of Finance*.
- Begenau (2020), "Capital requirements, risk choice, and liquidity provision in a business-cycle model", *JFE*.
- Blanchard, Giavazzi & Sa (2005), "International investors, the U.S. current account, and the dollar", *BPEA*.
- Bocola (2016), "The pass-through of sovereign risk", *JPE*.
- Broner, Erce, Martin & Ventura (2014), "Sovereign debt markets in turbulent times", *JME*.
- Brutti & Sauré (2016), "Repatriation of debt in the euro crisis", *JEEA*.
- Coeurdacier & Rey (2013), "Home bias in open economy financial macroeconomics", *JEL*.
- Corsetti, Kuester, Meier & Müller (2013), "Sovereign risk, fiscal policy, and macroeconomic stability", *EJ*.
- Crosignani, Faria-e-Castro & Fonseca (2020), "The (unintended?) consequences of the largest liquidity injection ever", *JME*.
- Dedola, Karadi & Lombardo (2013), "Global implications of national unconventional policies", *JME*.
- Devereux & Sutherland (2011), "Country portfolios in open economy macro-models", *JEEA*.
- Drechsler, Drechsel, Marques-Ibanez & Schnabl (2016), "Who borrows from the lender of last resort?", *JF*.
- Drechsler, Savov & Schnabl (2017), "The deposits channel of monetary policy", *QJE*.
- Farhi & Tirole (2018), "Deadly embrace", *REStud*.
- Gabaix & Maggiori (2015), "International liquidity and exchange rate dynamics", *QJE*.
- Gennaioli, Martin & Rossi (2014), "Sovereign default, domestic banks, and financial institutions", *JF*.
- Gertler & Karadi (2011), "A model of unconventional monetary policy", *JME*.
- Gertler & Kiyotaki (2010), "Financial intermediation and credit policy in business cycle analysis", *Handbook of Monetary Economics*.
- Gertler, Kiyotaki & Prestipino (2016), "Wholesale banking and bank runs in macroeconomic modeling of financial crises", *Handbook of Macroeconomics*.
- Itskhoki & Mukhin (2021), "Exchange rate disconnect in general equilibrium", *JPE*.
- Koijen, Koulischer, Nguyen & Yogo (2021), "Inspecting the mechanism of quantitative easing in the euro area", *JFE*.
- Koijen & Yogo (2019), "A demand system approach to asset pricing", *JPE*.
- Kollmann, Enders & Müller (2011), "Global banking and international business cycles", *EER*.
- Krishnamurthy, Nagel & Vissing-Jorgensen (2018), "ECB policies involving government bond purchases", *Review of Finance*.
- Morris (1991), "Factorial sampling plans for preliminary computational experiments", *Technometrics*.
- Neumeyer & Perri (2005), "Business cycles in emerging economies: the role of interest rates", *JME*.
- Ongena, Popov & Van Horen (2019), "The invisible hand of the government", *AEJ: Macro*.
- Schmitt-Grohé & Uribe (2003), "Closing small open economy models", *JIE*.
- Tille & van Wincoop (2010), "International capital flows", *JIE*.
- Uribe & Yue (2006), "Country spreads and emerging countries: who drives whom?", *JIE*.
