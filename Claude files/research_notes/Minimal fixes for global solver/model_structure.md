# Structural changes that make a Bocola/GK sovereign-bank model well-behaved at high leverage and exposure

Scope note. The primary PDFs of Bocola (2016), its appendix, Gertler-Karadi (2011) and Schmitt-Grohe-Uribe (2003) returned 403/502 to the fetch tool. The model-specific details of those papers below come from (i) abstracts and landing pages that were reachable, (ii) this repository's own records of Bocola's code and calibration (CLAUDE.md, `calibration/global_projection.py`), and (iii) the author's recall. Recall-only items are marked **[recall, unverified]** and are repeated under Gaps. Everything else marked **[inference]** is reasoning about this model, not something a paper says.

Problem labels used throughout: (a) a 63% net-worth loss in default, a binding floor, post-default states outside the box; (b) a large calm-time premium from a high default-state alpha (~7.5 against 1.3), which leaves mu ~1e-4 at rest; (c) a jumpy, non-unique rest point near 135-150 bp; (d) no deterministic SS above ~240 bp at f = 0.08; (e) slow eigenvalues (0.987-0.994) in the deposit position and the bank balance sheet.

## Q1. Recapitalisation and bailout rules in default (HFSF-like), and their effect on default-state valuation and calm-time premia

### Takeaway
An HFSF-style rule recapitalises banks back to a target capital ratio, with public money, after a sovereign haircut. Greece did exactly this in 2012-13: a 9% Core Tier 1 target, about EUR 37.7bn of PSI-related losses, and a EUR 50bn HFSF envelope. Replacing the model's 15% floor with "restore net worth to a share chi of its pre-default level, financed by new debt" lowers the default-state alpha toward the calm-state alpha. That cuts the calm-time premium in (b) at its source and pulls post-default states back inside the box (a). It is a one-equation change to the default branch, and it is the best payoff per effort in these notes.

### Cited Findings
- Bank of Greece assessed recapitalisation needs against a 9% Core Tier 1 minimum in the baseline and 7% in the adverse scenario, to end-2014 — [Eurobank H1 2012 financial report](https://athexgroup.gr/en/documents/10180/43442/Financial%20Report%20EUROBANK%20ERGASIAS%20S.A.%20%282012%2CSix-Month%20Statement%2CBoth%29.pdf/c5ef1103-5b7c-49f1-9354-cdb8f2be7a6b); the baseline used 9% for 2012 and 10% for 2013-14 — [Eurobank Greece Macro Focus, Jan 2013](https://www.eurobank.gr/Uploads/Reports/GREECE%20MACRFOCUS%20%20January%2025%202013.pdf)
- Bank of Greece put the losses from PSI plus selected state-related loans at about EUR 37.7bn, including EUR 5.8bn of provisions booked before PSI — [Eurobank Greece Macro Focus, Jan 2013](https://www.eurobank.gr/Uploads/Reports/GREECE%20MACRFOCUS%20%20January%2025%202013.pdf)
- Recapitalisation framework (Nov 2012 cabinet act): 9% EBA Core Tier 1, of which 6% in common shares and the rest possibly in HFSF-subscribed CoCos or preference shares. Private investors taking at least 10% of the raise got warrants. Alpha Bank's need was EUR 4.6bn, EUR 2.9bn of it already advanced by the HFSF — [Alpha Bank Q3 2012 results](https://www.athexgroup.gr/documents/10180/43440/64_26366_2012_English_2.pdf/1f379013-8de7-45f0-a82e-1b4587cf53e9); Eurobank's need was EUR 5.8bn, "commensurate with PSI losses" — [Eurobank filing](https://www.athexgroup.gr/documents/10180/43440/278_26554_2012_English_2.pdf/8f5312e2-765d-4d2e-b1f2-57864feb6228)
- The HFSF's total envelope for the banking system was EUR 50bn — [Eurobank 2012 investor presentation, via search summary](https://www.athexgroup.gr/documents/10180/43440/64_12407_2012_English_2.pdf/d5fb4620-a17f-4b55-9326-7ebe66e14185)
- Farhi-Tirole (REStud 2018), "double-decker bailout": the domestic government bails out its banks and foreign creditors or other countries may forgive or support the sovereign. They draw implications for the re-nationalisation of sovereign debt and for banking union — [NBER w21843](https://www.nber.org/system/files/working_papers/w21843/revisions/w21843.rev0.pdf); [TSE](https://www.tse-fr.eu/sites/default/files/TSE/documents/doc/by/tirole/doom_loops_122115.pdf)
- Hur, Sosa-Padilla and Yom, "Optimal Bailouts in Banking and Sovereign Crises": larger bailouts relax financial frictions and raise output but raise fiscal needs and default risk (a diabolic loop). Asset guarantees are the most common form of intervention in the data — [NBER w28412](https://nber.org/system/files/working_papers/w28412/w28412.pdf)
- Sosa-Padilla (2018): bankers hold government debt, so default causes a banking crisis, a corporate credit collapse and output drops of about 7 pp; calibrated to Argentina 2001-02 — [UAB record](https://ddd.uab.cat/record/196768?ln=en)
- Perez: the balance-sheet effect of default on banks dominates. Post-default bailouts are desirable ex post but reduce welfare ex ante, because they weaken the government's commitment not to default — [SIEPR WP](https://siepr.stanford.edu/publications/working-paper/sovereign-debt-domestic-banks-and-provision-public-liquidity); [IDEAS](https://ideas.repec.org/p/sip/dpaper/15-016.html)
- Faria-e-Castro (REStat 2024): bank recapitalisations and transfers carry the largest fiscal multipliers at the height of a financial crisis — [IDEAS / St. Louis Fed WP 2018-023](https://ideas.repec.org/p/fip/fedlwp/2018-023.html); [ECB conference paper](https://ECB.eu/pub/conferences/shared/pdf/20171120_fiscal_conference/3a_paper_Faria-e-Castro.pdf)
- Erce (2020): bank-sovereign pass-through is stronger when banks hold more domestic sovereign debt and when governments bail banks out — [Bocola search summary citing Erce](https://international.vlex.com/vid/bank-and-sovereign-risk-855624978)
- Bocola (2016) models two channels: a liquidity channel (default news tightens funding) and a risk channel (precautionary deleveraging). The risk channel was sizable in Italy — [Minneapolis Fed WP 722 landing page](https://minneapolisfed.org/research/working-papers/the-pass-through-of-sovereign-risk)
- **[recall, unverified]** Bocola's default event is a haircut on the bonds held by banks, with no separate bank bailout; the post-default economy is the same economy at a lower bond value. This repo's CLAUDE.md records the same design ("ONE pure-haircut feared event", recovery 0.45 here against his D = 0.55) — repository CLAUDE.md (local).

### Inferences
- [inference] Why recap fixes (b). The calm-time sovereign premium is roughly E[Omega' * (R^nd - R^d)] in the default branch, where Omega' = beta[f + (1-f) alpha']. With alpha^d ~7.5 against 1.3, the default branch's SDF weight is ~5x the no-default weight. Recapitalising to the pre-default ratio sets n^d close to n^nd, so alpha^d falls toward alpha^nd and the premium collapses toward the pure expected-loss spread (pi * (1 - recovery) * duration-adjusted, a few bp per quarter at pi = 0.1%). The calm-time Greek-German spread then moves toward the 25 bp in the 2002-07 data. The constraint is no longer relaxed by sovereign carry, so mu at rest rises toward its SS value. This plausibly removes the need to push the SS spread to 135-140 bp, and with it the jumpy region in (c).
- [inference] Why recap fixes (a). Post-default n_D is pinned at chi * n_pre rather than 0.37 * n_pre, so the default-branch continuation states lie in or near the box, and the floor stops binding at most nodes (it becomes the rule itself).
- [inference] What it costs the story. The pass-through survives through the mark-to-market loss on the no-default path, which is Bocola's object, because the IRF conditions on no default. It is weakened only through the risk channel: banks fear default less if they are made whole. Two ways to cap that cost. (i) Partial recap, chi < 1. The data supports partial: HFSF recap restored 9% CT1 but wiped out pre-existing shareholders, i.e. the bankers' own franchise was not restored. (ii) Recapitalise with equity that dilutes incumbent bankers: the government injects n, but only a share (1 - dilution) accrues to surviving bankers' value. In Omega-space this is alpha^d multiplied by a dilution factor. That matches the HFSF's actual terms (common shares, warrants to private investors) and keeps a precautionary motive.
- [inference] Implementation size. One equation in the default branch of `point_map` (the net-worth law of motion when d' = 1), plus the government budget term already in place for the floor top-up. Roughly the same as the existing 15% floor, so it is a substitution, not a new ingredient. The fiscal cost (new debt) is already wired.
- [inference] Calibration anchor. EUR 37.7bn of losses against EUR 50bn of HFSF capacity says the Greek recap was essentially full on PSI losses (Eurobank: "commensurate with PSI losses"). chi in 0.6-1.0 of pre-default net worth is defensible. Pre-PSI Greek bank equity (gaps) would let one compute chi exactly.
- [inference] Guarantee alternative (Hur-Sosa-Padilla-Yom: asset guarantees are the most common). A government guarantee on the banks' sovereign book is inconsistent here, since the guarantor is the defaulting sovereign. Recapitalisation through new debt (EFSF-funded in Greece) is the internally consistent version.

### Gaps
- Could not fetch Bocola (2016) or its appendix (403). I could not confirm from the source whether he has any post-default transfer to bankers or how he solves the post-default regime.
- Bocola-Dovis (2019) is a sovereign-maturity model with risk-averse lenders and no banks, so it says nothing about bank bailouts. Its relevance is the pricing kernel (Q2).
- I did not find the primary Bank of Greece "Report on the Recapitalisation and Restructuring of the Greek Banking Sector" (Dec 2012). The figures above are from bank filings and a bank research note that quote it.
- Pre-PSI aggregate Greek bank equity (to pin chi) was not found.

## Q2. Who prices sovereign bonds: foreign or diversified investors against domestic banks

### Takeaway
In much of the literature the marginal pricer of sovereign debt is a deep-pocketed foreign or diversified investor, often risk-neutral (Arellano-type models) or risk-averse with an exogenous kernel (Bocola-Dovis). Domestic banks then hold the bonds at that price, with crowding-out and creditor discrimination as the bank-side forces (Broner-Erce-Martin-Ventura). Letting an outside investor (or the F bank, at a fixed risk-neutral kernel) set Q_bD caps the premium at the expected loss. That kills the calm-time over-pricing (b) without touching the default regime. But it detaches the price from Greek bank net worth, which is part of Bocola's liquidity channel. Medium effort: one pricing equation is replaced, and D-bank holdings become a quantity choice at a given price.

### Cited Findings
- Bocola-Dovis (AER 2019): endogenous maturity, risk-averse lenders, Cole-Kehoe self-fulfilling crises. Lenders price via no-arbitrage. Non-fundamental (rollover) risk played a limited role in Italy 2008-12; drafts reported 20% and then 12% — [NBER w22694](https://www.nber.org/papers/w22694); [AER](https://www.aeaweb.org/doi/10.1257/aer.20161471)
- Broner-Erce-Martin-Ventura (JME 2014): in the periphery spreads rose and the domestic share of debt rose with them. In turbulent times domestic creditors earn higher expected returns than foreign ones (creditor discrimination). Domestic purchases crowd out private investment, and the crowding-out can spread across the euro area — [ESM/JME version](https://www.esm.europa.eu/publications/sovereign-debt-markets-turbulent-times-creditor-discrimination-and-crowding-out-effects); [BdE](https://www.bde.es/wbe/en/publicaciones/analisis-economico-investigacion/documentos-trabajo/sovereign-debt-markets-turbulent-times-creditor-discrimination-and-crowdingout-effects.html)
- Perez: the government issues both external and domestic debt; banks are heterogeneous — [SIEPR](https://siepr.stanford.edu/publications/working-paper/sovereign-debt-domestic-banks-and-provision-public-liquidity)
- Arellano-Bai-Bocola, "Sovereign (Default) Risk and Firm Heterogeneity" (JEEA final version, Oct 2025 per Bocola's site) — [Minneapolis Fed SR 547 listing](https://ideas.repec.org/p/fip/fedmsr/547.html); [Bocola homepage](https://lbocola.people.stanford.edu/). **[recall, unverified]** In that paper bonds are priced by risk-neutral (or exogenous-kernel) foreign lenders, and the bank-side transmission works through sovereign spreads raising firms' borrowing costs.
- Kirchner-van Wijnbergen (JME 2016): leverage-constrained banks holding government bonds crowd out private credit; deficit-financed stimulus weakens when banks hold a large share of the debt — [Tinbergen](https://tinbergen.nl/publication/159452/fiscal-deficits-financial-fragility-and-the-effectiveness-of-government-policies); van der Kwaak-van Wijnbergen: sovereign discounts appear near the limit of bank balance-sheet capacity — [TI 2012-044](https://papers.tinbergen.nl/12044.pdf)

### Inferences
- [inference] There is a minimal version that needs no new agent. The model already has a safe, diversified F bank whose default-state alpha is near 1 (F loses only its b_DF/sz exposure). Restricting the D-bond's marginal pricer to F (D banks hold b_DD at the F-priced Q_bD, with b_DD set by a portfolio cost or a quantity rule) is the segmented-market version. The cheaper option is to keep the D-bank FOC and add a risk-neutral foreign investor with a holding cost (a Broner-type demand curve). With a deep foreign demand the spread is bounded by pi*(1-R)/duration plus the holding cost.
- [inference] What it costs. Bocola's mechanism requires the bond price to move with an exogenous pi, and that survives. What is lost is the feedback from bank net worth to Q_bD (bank-pricing liquidity effects). It also makes the TPI analysis about a different object: the CB would then be competing with a deep market, not relieving a balance-sheet-constrained one. That conflicts with the repo's current TPI finding that price impact runs mostly through the continuation. So this change touches the paper's application; recap (Q1) does not.
- [inference] Broner et al.'s creditor discrimination justifies a model where D banks earn MORE than foreigners in turbulent times. Paired with Q1 (recap), it also gives a rationale for D-bank holdings rising with spreads without the calm-time over-pricing.

### Gaps
- Full texts of Arellano-Bai-Bocola, Perez and Broner et al. were not retrieved, so their exact pricing-kernel specifications are unverified.

## Q3. Constant intermediation wedges against constraint-only spreads; joint calibration of exit rate, leverage and spread; SS existence

### Takeaway
In GK/GKi models the SS spread is the constraint's shadow value alone. Its size is tied mechanically to the exit/survival rate, leverage and the divertable share, and for a given target pair (leverage, spread) an SS exists only for a range of survival rates, which is problem (d). Two small fixes. (i) Add a constant cost-of-intermediation wedge (a resource cost per unit of loans, or a fixed loan-management cost) so the target loan spread splits into wedge plus shadow value. The constraint then needs only a modest mu at rest to bind, which removes the need for the 135-140 bp region (c-d). (ii) Recalibrate f jointly with leverage 8, rather than taking f = 0.08 as given.

### Cited Findings
- GK (2011) SS targets as reported by replicators: leverage 4, spread 100 bp (pre-2007 BAA-government and mortgage-government gaps) — [Dynare forum, user report, not the paper](https://forum.dynare.org/t/other-important-questions-about-gertler-and-karadi-2011/17492). **[recall, unverified]** GK (2011) set banker survival theta = 0.972 (expected horizon ~36 quarters), divertable share lambda = 0.381, and start-up transfer omega = 0.002.
- This repo records Bocola's exit share as 0.0354 (survival ~0.965) against the 0.08 used here, described as "the bank block's one departure" — `calibration/global_projection.py`, line 30 (local).
- Some replications of GK (2011) fail to find a steady state (NaN/Inf), so SS existence is a known practical issue in this class — [Dynare thread](https://forum.dynare.org/t/gertler-karadi-2011/27491)

### Inferences
- [inference] SS existence (d). In GK the SS franchise value solves a quadratic in leverage phi from nu-eta conditions. A real root requires that the survival-weighted excess return times leverage not exceed what the exit rate allows: beta(1-f)(R_k - R)phi + R < 1/(beta(1-f))-type bounds. The repo's own note, "raise f until leverage is the least root", is this condition. Higher spread targets at fixed leverage push the required survival up, until no root exists at f = 0.08. So the existence frontier is a curve in (f, spread, leverage) space. Mapping it numerically before choosing targets (a 2-D grid of SS solves, minutes) is the cheapest diagnostic of all.
- [inference] Constant wedge. Write R_loan - R_dep = tau + lambda*mu/E[Omega], with tau a resource cost (paid to households or lost). Calibrate tau to (data spread - target shadow spread), e.g. tau = 60-80 bp and a shadow part of 20-40 bp. At rest mu is then pinned by a smaller shadow target, bound strictly away from 0 (the cure the repo found at 100 bp), with no need to move to 135-140 bp. Cost: the wedge is not state-contingent, so the pass-through comes only from the shadow part. Impact output responses shrink unless tau is small. It is one line in the loan-rate / working-capital equation and one SS calibration line. Note: the repo has already found that the shadow part must be ~100 bp/yr for mu to bind robustly at rest (CLAUDE.md "CURED 2026-08-29"), so the wedge only helps if Q1 removes the sovereign carry that is relaxing the constraint. Done alone it will push mu back toward the kink.
- [inference] f for stationarity (e). A higher exit share speeds net-worth mean reversion: the net-worth eigenvalue is roughly (1-f)*R_on-equity. The slow eigenvalues 0.987-0.994 are consistent with (1-f)*(1+roe) ~ 0.99. Raising f to ~0.10-0.12 (horizon 8-10 quarters) moves the root, but it lowers alpha and the franchise channel. A lower f (Bocola's 0.035) does the reverse.

### Gaps
- GK (2011) and Gertler-Kiyotaki (2010, 2015) primary texts not retrieved; the survival / divertable / start-up values above are recall.
- No source found that gives an explicit closed-form SS existence condition for the multi-asset GK bank. The bound above is the author's derivation sketch.

## Q4. Priced and physical default probability; lower priced hazard at rest

### Takeaway
Rest-point pi = 0.10%/qtr with recovery 0.45 implies an expected-loss spread of roughly 0.10% x 0.55 x 4 ~ 22 bp/yr before any risk premium, already near the whole 2002-07 Greek-German spread. Every bp of risk premium the banks add (problem b) pushes the model above the data. Lowering the rest hazard to 0.02-0.05%/qtr is a calibration change, not a structural one. It shrinks the default branch's weight and the calm premium proportionally, so it is the cheapest partial fix for (b). It leaves (a), the size of the default-state loss, untouched.

### Cited Findings
- Bocola-Dovis use risk-averse lenders, so priced default exceeds physical default; they find limited non-fundamental risk in Italy — [AER 2019](https://www.aeaweb.org/doi/10.1257/aer.20161471)
- This repo's Bocola design prices default (def_price) separately from realised default (def_real); the IRF prices risk but never realises it — repository CLAUDE.md (local)

### Inferences
- [inference] In this model the D bank is the pricing kernel, so a lower priced pi at rest is the same object as a lower physical pi. The distinction only matters if a separate foreign kernel is added (Q2). A basis-adjusted 0.02-0.05% puts the expected-loss part at 4-11 bp/yr, leaving room for a bank premium inside 25 bp.
- [inference] Cost. A lower rest pi with an unchanged headline shock means a larger proportional shock; the logistic p^d(s) curvature moves, and the s-box and the S_REFINE accuracy should be rechecked.

### Gaps
- No primary source retrieved for basis-adjusted (CDS-bond) Greek default intensities 2002-07; the 0.02-0.05% range in the brief is the caller's and is not verified here.

## Q5. Stationarity devices for cross-border positions and bank net worth (Schmitt-Grohe-Uribe 2003)

### Takeaway
SGU (2003) show that the closing devices (endogenous discounting, a debt-elastic premium, portfolio adjustment costs, complete markets) give near-identical business-cycle dynamics. They are chosen for convenience. kappa = 0.05 not helping suggests the slow root in (e) is not the deposit position's unit root but the bank net-worth / franchise recursion (0.99 per period, which the repo already measured in time iteration). The fix is then on the bank side (exit rate f, recap, a dividend-payout rule tied to excess net worth), not the external side.

### Cited Findings
- SGU (JIE 2003, 61(1): 163-185) compare an endogenous discount factor (Uzawa), a debt-elastic interest rate, portfolio adjustment costs, complete markets, and a no-stationarity-device benchmark. The only noticeable difference is smoother consumption under complete markets — [NBER w9270](https://www.nber.org/papers/w9270); [IDEAS](https://ideas.repec.org/a/eee/inecon/v61y2003i1p163-185.html)
- Users report that the debt-elastic and portfolio-cost versions are easier to replicate than the endogenous-discount ones — [Dynare archive](https://archives.dynare.org/phpBB3/viewtopic.php%3Fp=9548.html)
- In this repo, time iteration's binding mode is "the franchise-value recursion at 0.990 per sweep" — repository CLAUDE.md (local)

### Inferences
- [inference] Diagnose before changing. Compute the eigenvector of the 0.987-0.994 root. If it loads on n_D/V_dep jointly, the deposit position is just the counterpart of bank net worth (banks fund through cross-border deposits), and a stronger kappa cannot help, as observed.
- [inference] Bank-side stationarity device, one equation: make the payout share state-dependent, f_t = f + phi_f*(n_t/n_ss - 1), i.e. banks pay out more when net worth is above target. Gertler-Kiyotaki-Prestipino-type models and many DSGE-bank papers use payout or exit-rate variations, **[recall, unverified]**. It shortens the net-worth half-life without changing the SS. Cost: it dampens the persistence of the bank-capital loss, which is part of the story, so keep phi_f small and report it.
- [inference] Bocola-faithful alternative for the external side: his R = 1/beta + 0.01*(B_for/gdp) is a world rate with a debt-elastic premium (repo CLAUDE.md). That is the SGU option already used. Endogenous discounting would add a state; skip it.

### Gaps
- Full SGU text not retrieved; the eigenvalue magnitudes for each device are not reported here.

## Q6. Keeping the default regime tractable

### Takeaway
The standard tractable designs are: (i) solve the post-default economy as its own regime on its own state domain (an absorbing, or exclusion-then-re-entry, economy with its own rules), and (ii) make default an event whose aftermath is pinned by policy (recap, Q1) so that the post-default state lands on a known manifold. Cole-Kehoe and Arellano-type models use exclusion with a re-entry probability. Bocola's design reads the same rules at the post-default point. That is what fails here when n falls 63%.

### Cited Findings
- Sosa-Padilla and Perez model default with domestic banks holding the debt, with sharp post-default credit collapses (Sosa-Padilla: ~7 pp output drop) — [UAB](https://ddd.uab.cat/record/196768?ln=en); [SIEPR](https://siepr.stanford.edu/publications/working-paper/sovereign-debt-domestic-banks-and-provision-public-liquidity)
- Bocola-Dovis embed Cole-Kehoe self-fulfilling crises — [NBER w22694](https://www.nber.org/papers/w22694)
- **[recall, unverified]** Arellano (2008)-type models use exclusion from markets after default with an exogenous re-entry probability (~0.28/qtr in several calibrations) and an output cost during exclusion.

### Inferences
- [inference] Option A, separate box for the default regime: a second Smolyak grid over the post-default state range (n_D at 30-50% of SS), with d'=1 expectations evaluated there. The economy re-enters the main box via a re-entry probability or by net-worth recovery. This is a large change: a second set of rules, the ladder doubled, and the n^2 Jacobian cost.
- [inference] Option B, absorbing default with a closed-form value: in default, the bank's continuation value is approximated by the steady state of a post-default economy (alpha^d computed once, offline, as a constant). That is one scalar instead of a branch of rules. It is cheap, but it replaces Bocola's "same rules at a reachable point" with a frozen stand-in. The repo's history flags frozen stand-ins as the source of a past sign error, so it needs a check against the full branch at nodes inside the box.
- [inference] Option C, recap (Q1): post-default n is pinned near the box. This is the smallest change and dominates A and B whenever the economics allow it.

### Gaps
- No source retrieved for how Bocola (2016) solves the post-default continuation; his appendix was not accessible.

---

## Ranked list (payoff per effort) — [inference throughout, based on the findings above]

| Rank | Change | Fixes | Size | Main cost to story |
|---|---|---|---|---|
| 1 | **HFSF-style recap in default**: n^d = max(n_haircut, chi*n_pre), chi ~0.6-1.0, debt-financed (replaces the 15% floor); optionally dilute incumbents' share | (a) directly; (b) by lowering alpha^d; likely (c) and (d), since the SS spread need not rise | 1 equation in the default branch, plus the existing fiscal term | Weaker risk channel (smaller fear of default); keep partial or dilutive. Data-backed (9% CT1, EUR 37.7bn losses, EUR 50bn HFSF) |
| 2 | **Lower rest-point priced hazard** to 0.02-0.05%/qtr | (b) proportionally | calibration only | Shock scale and s-box change; recheck accuracy |
| 3 | **Map the SS existence frontier in (f, leverage, spread)** and recalibrate f jointly (Bocola's own 0.035, or higher) | (d); informs (c) and (e) | an offline grid of SS solves | None; diagnostic |
| 4 | **Diagnose the slow root's eigenvector**; if it is bank-side, add a state-dependent payout f_t = f + phi_f*(n/n_ss - 1) | (e) | 1 equation | Less persistent bank-capital loss |
| 5 | **Constant intermediation wedge tau** in the loan spread | (c) and (d), by letting the shadow part be smaller | 1 equation plus 1 SS line | Pass-through only via the shadow part; only useful after #1, otherwise mu returns to the kink |
| 6 | **Outside (foreign/risk-neutral) marginal pricer of D bonds** / segmented market | (b) strongly | replace 1 FOC; holdings become quantities | Detaches Q_bD from bank net worth; changes the TPI's object |
| 7 | **Separate default-regime domain / absorbing default value** | (a) | large (second grid) or medium (frozen alpha^d) | Frozen stand-ins have a bad history here |

The caller's memory note asks that every modelling choice be put to the user with measured numbers. Each row above is a candidate to propose, not a decision.
