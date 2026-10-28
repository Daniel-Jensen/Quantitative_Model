# Numerical methods for converging a Chebyshev-Smolyak collocation solution with an occasionally binding constraint (minimal-change catalogue)

Scope note. Failure codes used throughout: (1) no collocation root once cross terms are added where the constraint switches (sigma_min/sigma_max ~1e-6, LM stalls at 3e-3); (2) wider boxes fail from cold start; (3) time iteration (TI) and Newton select different solutions, TI contracts at 0.99; (4) non-unique no-shock rest point, slack starts drift to box walls; (5) fitted rules read slack where the cleared period map binds; (6) wider eps breaks coarse solve; always-binding gives negative multipliers at a third of nodes.
Source limits: several primary PDFs were fetched and text-extracted (Maliar-Maliar-Winant 2021; Fernandez-Villaverde/Rubio-Ramirez/Schorfheide handbook draft 2015; Judd-Maliar-Maliar-Valero 2014). Christiano-Fisher (2000), Brumm-Scheidegger (2017), Maliar-Maliar (2015), Rendahl (2015), Guerrieri-Iacoviello (2015), Azinovic et al. (2022) were seen at abstract level only. Bocola's replication code and the Gertler-Kiyotaki-Prestipino numerical appendix were NOT found online. Anything not traced to a fetched text is in Inferences or Gaps.

## Q1. Least-squares (overdetermined) projection vs square collocation for kinked policies

### Takeaway
Least squares (more nodes than coefficients, Gauss-Newton/LM on sum of squared residuals) always has a minimiser and removes the "no root exists" failure of square collocation, but it does not remove ill-conditioning by itself; the literature recommends it mainly for stability with oversampling, and the handbook explicitly warns it can be ill-conditioned. It is the cheapest change that directly targets failure (1).

### Cited Findings
- The handbook frames least-squares projection as minimising the integral of R^2 over the domain, "mathematically equivalent to a standard regression problem"; it "always generates symmetric matrices" but "may lead to ill-conditioning and systems of equations complicated to solve numerically." — [Fernandez-Villaverde, Rubio-Ramirez, Schorfheide, handbook draft, s.5.12](https://www.hoover.org/sites/default/files/fernandez-rubio-shorfheide-hbofmacro_solveestimate_v1.pdf)
- Same source: collocation is "simple to compute ... specially attractive when dealing with strong nonlinearities"; orthogonal (Chebyshev-root) collocation has "surprisingly good performance"; Galerkin is "highly accurate and robust but difficult to code." — [handbook draft, s.5.14-5.15](https://www.hoover.org/sites/default/files/fernandez-rubio-shorfheide-hbofmacro_solveestimate_v1.pdf)
- In multivariate discrete least squares with polynomial bases, stability and convergence are proven when the number of sample points scales roughly quadratically in the dimension of the polynomial space (random/unstructured samples). — [Multivariate discrete least-squares approximations with a new type of collocation grid (arXiv 1401.0894)](https://arxiv.org/pdf/1401.0894); weighted-sampling variants in [Christoffel-weighted least squares (arXiv 1412.4305)](https://arxiv.org/pdf/1412.4305)
- Maliar-Maliar-Winant solve their kinked consumption-savings model by minimising squared (Fischer-Burmeister-transformed) Euler residuals over simulated points, i.e. a regression/least-squares formulation, reaching residuals of order 10^-2.75 (0.18%) "given that we solve the model with a kink in the decision rule." — [Maliar, Maliar, Winant (2021, JME), s.4](https://web.stanford.edu/~maliars/Files/JME2021.pdf)

### Inferences
- Existence: sum-of-squares is continuous and bounded below in theta, so a minimiser exists on any compact coefficient set; square collocation needs F(theta)=0 to have a root, and a smallest singular value ~1e-6 at the best LM point with max|F|~3e-3 is the classic signature of a fold / near-tangency where the root has disappeared. Least squares turns "no root" into "nonzero minimum", which you can measure (the residual at the minimum is the approximation error of the kink, not a solver failure).
- It does not by itself avoid folds: the Gauss-Newton Hessian J'J inherits J's near-singularity. What helps is that with oversampling the extra rows typically add rank in exactly the directions where the square system is singular (nodes on both sides of the switching surface constrain the cross-term coefficient that the square grid leaves under-determined).
- Practical sizing (inference from the 1401.0894 result and the usual econ practice): 1.5-3x as many points as coefficients. Minimal implementation in this code: keep the same Smolyak basis (e.g. the 59- or 115-coefficient refined basis) but evaluate residuals on a denser point set (the next Smolyak level's nodes, or a Chebyshev tensor in s x Smolyak rest) and switch the existing LM fallback to minimise ||F||^2 rather than treat it as a root-finder; scipy.optimize.least_squares(method='trf' or 'lm', x_scale='jac') with the existing FD Jacobian is a drop-in. Acceptance test must change from max|F|<=1e-9 to "gradient J'F ~ 0 and max|F| reported as accuracy".
- Weight the rows: residuals in levels of different units (FB residuals vs Euler errors vs identity log-residuals) should be scaled so no equation class dominates; MMW explicitly note the FB arguments may need a relative weight v (Psi_FB(a, v h)). — [MMW 2021, eq. 25 discussion](https://web.stanford.edu/~maliars/Files/JME2021.pdf)
- Fixes: (1) directly; partly (2) (minimiser always exists, so cold starts fail softly); not (3)-(4). Effort: minimal (hours to a day).

### Gaps
- Did not retrieve Judd (1992, JET) or Judd (1998, book) text comparing least squares and collocation; no economics paper found quantifying "extra points needed" for kinked policies specifically.

## Q2. Approximate smooth objects instead of kinked ones (PEA, end-of-period expectations, closed-form multiplier)

### Takeaway
The strongest literature-backed fix for kinks: approximate the conditional expectation (smooth) rather than the policy or multiplier (kinked), and recover the multiplier/policy exactly from the expectation via max/min at each point. Christiano-Fisher (2000) found this parameterised-expectations family dominates the alternatives they tested for occasionally binding constraints. In this solver it means storing E[Omega'R'] (or the post-decision expectation) as the rule and computing mu and alpha in closed form, instead of storing alpha = E[Omega]R/(1-mu), which inherits mu's kink.

### Cited Findings
- Christiano and Fisher (2000, JEDC 24(8):1179-1232) compare several algorithms on the growth model with irreversible investment; "the parameterized expectations approach performs best", their PEA variants "appear to dominate the other algorithms" on speed, accuracy and convenience. — [NBER t0218](https://www.nber.org/papers/t0218); [IDEAS/JEDC](https://ideas.repec.org:443/a/eee/dyncon/v24y2000i8p1179-1232.html)
- The Chebyshev-PEA of Christiano-Fisher is still used for occasionally binding borrowing constraints in recent OLG work (2024). — [IDEAS nzt/nztwps/24-03](https://ideas.repec.org/p/nzt/nztwps/24-03.html)
- Guerrieri-Iacoviello (2015 JME) use a Christiano-Fisher-type global solution as the "virtually exact" benchmark against which OccBin is compared. — [JME 2015 paper](https://www.matteoiacoviello.com/research_files/JME_2015.pdf); [FEDS page](https://www.federalreserve.gov/econres/feds/occbin-a-toolkit-for-solving-dynamic-models-with-occasionally-binding-constraints-easily.htm)
- OccBin-type perturbation decision rules "abstract from precautionary savings and disregard the risk effects of the occasionally binding constraint" — so OccBin is not a substitute here, where the risk effect is the object. — [Mendoza-Villalvazo / PIER WP 20-028](https://www.econ.upenn.edu/system/files/working-papers/20-028%20PIER%20Paper%20Submission.pdf)
- Project context (not literature): the current solver already recovers mu in Bocola's closed form at each point, but stores alpha (which contains 1/(1-mu)) as an approximated rule with an identity residual log(guess/implied). — `/Users/Huawei/Desktop/Quantitative_Model/CLAUDE.md`

### Inferences
- Failure (5) (fitted rules read slack where exact clearing binds, 22 vs 115 bp) is exactly what PEA avoids: if the stored object is smooth and mu is always recomputed by max(.) from it, the interpolant can never "read" a multiplier — mu at any off-grid point is exact given the interpolated expectation. Current design stores something kinked (alpha via 1/(1-mu)) and therefore interpolates across the kink (Gibbs). Changing the stored rule from alpha to the smooth expectation E[Omega' R^k'] (or E[Omega'] and the franchise term separately) and computing mu = max(z,0) then alpha in closed form at every node and every continuation evaluation is a local change to the packing and to point_map's continuation read.
- This also attacks (1): the cross-term coefficients in a smooth expectation are well identified; the near-singular directions likely come from kinked stored rules whose coefficients are fitting the switching surface.
- It also attacks (6) and the always-binding failure: with exact max(.) recovery there is no eps to tune and no sign violation; a small eps can be kept only for Newton differentiability.
- Post-decision (end-of-period) state version (Carroll-style EGM logic): approximate the expectation as a function of end-of-period choices (K', b', deposits) rather than beginning-of-period states; this makes the expectation smooth even if the policy has a kink in states. More invasive here (the unknowns change), so moderate effort.
- Effort: minimal-to-moderate (change which 1-3 of the 20 stored rules are approximated and how alpha/mu are read; days). Fixes (5) strongly, (1) and (6) plausibly; does not address (3)-(4).

### Gaps
- Could not retrieve the Christiano-Fisher full text (404 on the Cleveland Fed PDF) to quote accuracy tables or their exact recommended PEA variant (they compare several PEA variants with Chebyshev collocation).
- No source found on how Bocola's own replication code stores the franchise value; Bocola's replication files were not located online.

## Q3. Ergodic-set grids vs boxes; adaptive/rotated domains

### Takeaway
A grid on the high-probability region (ergodic set) or a box rotated to principal components is consistently more accurate for a given number of coefficients than an axis-aligned hypercube, and avoids fitting the solution far from where it matters. For this solver the cheapest version is JMMV's adaptive parallelotope (PCA rotation of a simulated cloud) feeding the existing Smolyak code; caveat: rotation can create cross terms that mu=1 Smolyak cannot represent.

### Cited Findings
- JMMV (2014, JEDC 44:92-123): build "a parallelotope that encloses a high-probability area of the state space ... reorienting it with a principle-component transformation"; focusing on "a relevant domain" is more accurate there than fitting a larger partially irrelevant domain; anisotropic grid + adaptive domain cut maximum Euler residuals by about an order of magnitude even with two states, "two orders of magnitude in some examples" in multi-dimensional problems, holding the number of coefficients roughly fixed. — [Judd, Maliar, Maliar, Valero (2014)](https://lmaliar.ws.gc.cuny.edu/files/2019/01/JEDC2014_JMMV.pdf)
- Same paper, pitfall: "a Smolyak method with the level of approximation mu=1 has no cross terms"; interactions small in original coordinates "may become important when we rotate ... the quality of approximation can decrease rather than increase. Some experimentation may help." — [JMMV 2014, s.5.3](https://lmaliar.ws.gc.cuny.edu/files/2019/01/JEDC2014_JMMV.pdf)
- JMMV also find TI "may dominate fixed-point iteration in small-scale applications but the situation reverses when the dimensionality increases." — [JMMV 2014](https://lmaliar.ws.gc.cuny.edu/files/2019/01/JEDC2014_JMMV.pdf)
- Maliar-Maliar (2015, QE 6(1):1-47): simulate to estimate the ergodic set, then place a fixed grid of representative points over it (clustering, epsilon-distinguishable set, locally adaptive EDS) and solve by projection; handles hundreds of states, applied to a large NK model with a ZLB. — [QE page](https://www.econometricsociety.org/publications/quantitative-economics/2015/03/01/Merging-simulation-and-projection-approaches-to-solve-high-dimensional-problems-with-an-application-to-a-new-Keynesian-model); [author PDF](https://lmaliar.ws.gc.cuny.edu/files/2019/01/QEMM2015.pdf); [NBER w18501](https://www.nber.org/papers/w18501.pdf)

### Inferences
- EDS/cluster grids are unstructured, so they pair naturally with least-squares regression on a fixed polynomial basis (Q1) rather than with Smolyak interpolation. Combining "polynomial basis + more points drawn from the ergodic cloud + least squares" is a single coherent minimal change.
- Failure (2) (wider boxes fail) suggests the opposite direction to the instinct: make the domain smaller and better shaped, not wider. Failure (4) (slack starts drift to box walls) is then partly an artefact of extrapolating/saturating outside the relevant set; on an ergodic grid those states are by construction not where the solution is fitted, which may remove spurious attractors at the walls (inference; untested).
- Chicken-and-egg: you need a solution to simulate the ergodic set. Iterate: coarse box solve -> simulate -> PCA/ergodic grid -> re-solve -> resimulate (JMMV and MM both describe this loop).
- Note the project already uses a shear in (b_DD, M_cb), i.e. a hand-built parallelotope — so the code supports non-axis domains.
- Effort: PCA box minimal-moderate (one linear map before the existing Smolyak code); EDS grid + regression moderate. Fixes (2), possibly (4), improves (5) accuracy per coefficient.

### Gaps
- No source specifically on ergodic grids for models with two regimes and an absorbing-ish default fork.

## Q4. Adaptive sparse grids with local bases; deep-learning solvers

### Takeaway
Local piecewise-linear adaptive sparse grids (Brumm-Scheidegger) are the textbook answer to kinks in high-dimensional TI, but they replace the global-Chebyshev Newton with TI on a hierarchical local basis — effectively a rewrite of the solver core. Deep-learning solvers (MMW 2021, DEQN 2022) handle kinks via a Fischer-Burmeister loss on simulated paths, but are also a rewrite and give lower accuracy (~1e-3) than the 1e-9 target used here.

### Cited Findings
- Brumm-Scheidegger (2017, Econometrica 85(5):1575-1612): adaptive sparse grids inside time iteration / VFI; grid points are refined where the function has steep slopes or non-differentiabilities; tested on high-dimensional IRBC with irreversible investment (an occasionally binding constraint). — [Econometrica page](https://www.econometricsociety.org/publications/econometrica/2017/09/01/using-adaptive-sparse-grids-solve-high-dimensional-dynamic); [Dynare 2015 WP](https://archives.dynare.org/DynareConference2015/papers/Brumm_Scheidegger.pdf)
- FV-RR-S: finite elements let you "easily handle issues like kinks or constraints", with small elements where the economy spends time. — [handbook draft s.5.10](https://www.hoover.org/sites/default/files/fernandez-rubio-shorfheide-hbofmacro_solveestimate_v1.pdf)
- MMW (2021): rewrite Kuhn-Tucker A>=0, H>=0, AH=0 as the Fischer-Burmeister equation a+h-sqrt(a^2+h^2)=0 with unit-free a = 1-c/w, h = 1-beta r E[u'(c')]/u'(c); "differentiable", same solution as KT; may need a relative weight v; Euler residuals ~10^-2.75 with a kink; TensorFlow code provided. — [MMW 2021 JME](https://web.stanford.edu/~maliars/Files/JME2021.pdf)
- Azinovic-Gaegauf-Scheidegger (2022, IER 63(4):1471-1525): neural nets trained to satisfy equilibrium conditions on simulated paths, directly approximating all equilibrium functions, applied to OLG models with occasionally binding constraints. — [SSRN](https://papers.ssrn.com/abstract=3393482)

### Inferences
- For this code, a low-cost borrowing from Brumm-Scheidegger without a rewrite: keep the global Chebyshev basis for smooth rules and only switch the multiplier-bearing object to exact recovery (Q2). A hybrid "piecewise basis in s only" (e.g. splitting the s interval at the switching point into two Chebyshev elements) is a moderate change that keeps everything else.
- DL: rewrite (weeks), accuracy targets incompatible with the current acceptance thresholds. Not recommended for this deadline.

### Gaps
- Did not confirm whether DEQN uses FB (abstract-level only).

## Q5. Time iteration acceleration, EGM with OBC, equilibrium selection

### Takeaway
Rendahl (2015) proves TI on the Euler equations converges to the true solution under broad conditions even with occasionally binding constraints, which is the theoretical reason to trust TI's limit when it and Newton disagree. Anderson acceleration is a generic, small-code way to speed a linearly contracting fixed point (0.99 per sweep) but there is no economics-specific evidence found.

### Cited Findings
- Rendahl (2015, Economic Journal 125(585):1110-1135): time iteration "delivers a sequence of approximate policy functions that converges to the true solution under a wide range of circumstances", including many endogenous/exogenous states and "a very broad spectrum of occasionally binding constraints"; proof via an envelope theorem establishing differentiability of each iterate's value function. — [Cambridge repository](https://www.repository.cam.ac.uk/handle/1810/255270?show=full); [EUI precursor](https://cadmus.eui.eu/entities/publication/9c10546a-d361-56a4-b4dc-75e01728d7da/full)
- Anderson acceleration (Anderson/Pulay mixing/DIIS) reuses several past iterates; proven to improve the linear convergence rate of linearly convergent fixed-point iterations but not of quadratically convergent ones. — [Wikipedia overview](https://en.wikipedia.org/wiki/Anderson_acceleration); [Linear asymptotic convergence of AA (arXiv 2109.14176)](https://arxiv.org/pdf/2109.14176); [Acceleration methods for fixed point iterations (arXiv 2507.11746)](https://arxiv.org/pdf/2507.11746)
- Mendoza-Villalvazo's FiPIt (fixed-point iteration on a grid with two endogenous states and OBCs) is reported significantly faster than standard TI and hybrids. — [Review of Economic Dynamics / IDEAS](https://ideas.repec.org/a/red/issued/19-424.html)

### Inferences
- Equilibrium selection (3): Rendahl's result applies to single-agent-style Euler systems (or planner problems); this model is a decentralised two-country GE with banks, so uniqueness is not guaranteed and his theorem need not transfer. Still, the standard practical rule is: the TI limit is the limit of finite-horizon economies (the "Markov equilibrium as limit of finite horizon"), so it is the defensible selection; a Newton root that TI does not reach from it is suspect (could be a spurious collocation root created by the approximation, especially with kinked stored rules). Recommended test: start Newton from the TI limit; if Newton walks away, the two disagree on the approximated problem, not the economics — then compare both on off-grid Euler errors.
- A 0.99 contraction means ~2,000 sweeps to 1e-9; Anderson with memory 5-10 on the vector of stored coefficients is ~30 lines in numpy and can plausibly cut this by an order of magnitude (inference). It also offers a cleaner selection device: accelerated TI converges to the TI fixed point, not to whichever Newton root is nearest.
- The slow mode at 0.99 is identified in the project as the franchise-value recursion (alpha'); applying PEA (Q2) to that object, or solving the linear franchise-value recursion directly (it is linear in alpha' given mu), would remove the slow mode at source (inference).
- Effort: Anderson minimal; using TI-as-selection minimal. Fixes (3), speeds the warm start.

### Gaps
- No economics paper found that applies Anderson acceleration to time iteration; evidence is from numerical analysis/RL.
- No source on GE uniqueness for GK-bank two-country models.

## Q6. Homotopy / pseudo-arclength continuation for collocation systems with folds

### Takeaway
No economics source was retrieved in this pass. From numerical-continuation practice (inference): when a parameter or resolution homotopy hits sigma_min -> 0, a natural-parameter ladder (what the solver does) fails at the fold, while pseudo-arclength continuation can pass it; but a fold in the *discretisation* homotopy (coarse -> refined grid) usually signals that the refined collocation problem has no nearby root, which continuation cannot fix — least squares or a smoother stored object can.

### Cited Findings
- None retrieved. (The handbook draft lists homotopy among solver tools but the passage was not located in the extracted text.) — [handbook draft](https://www.hoover.org/sites/default/files/fernandez-rubio-shorfheide-hbofmacro_solveestimate_v1.pdf)

### Inferences
- Diagnostic before implementing: track sigma_min along the existing homotopy (leverage 5 -> 8, exposure, eps). If sigma_min -> 0 at a finite step and the residual min stays >0 after the fold, it is a turning point (arclength can follow it, possibly to a different branch = equilibrium multiplicity, related to (4)). If the refined grid simply has no root near the coarse one, continuation is the wrong tool.
- Effort: moderate (augment F with the arclength equation, tangent predictor, existing Newton as corrector; ~100 lines). Lower payoff than Q1/Q2 here.

### Gaps
- No economics-specific literature on pseudo-arclength for projection solutions located; Eaves-Schmedders / Judd homotopy references not retrieved.

## Q7. Smoothing KKT (FB, penalty, Garcia-Zangwill) and how Bocola / GKP handle the constraint

### Takeaway
FB is the standard differentiable reformulation (used e.g. by MMW 2021) and is scale-invariant in its zero set, but any smoothing replaces a kink by a high-curvature region of width ~eps that a global low-order polynomial still has to fit; smoothing helps Newton, not approximation. Failure (6) (wider eps breaks the coarse solve) is consistent with eps shifting the economics (a smoothed mu>0 where truth is 0 changes the rest point), not with numerics.

### Cited Findings
- FB: same solution set as KT, differentiable; MMW use unit-free arguments and note a possible relative weight; reference Jiang (1996) smoothed FB. — [MMW 2021](https://web.stanford.edu/~maliars/Files/JME2021.pdf)
- GKP (Review of Economic Studies 2020) solve their banking-panic model by time iteration (forum report, not the paper); replication code reportedly linked from Prestipino's site but unverified. — [Dynare forum](https://forum.dynare.org/t/help-for-a-macroeconomic-model-with-financial-panics-code/13378); [IFDP 1219](https://www.federalreserve.gov/econres/ifdp/files/ifdp1219.pdf)
- Bocola (2016, JPE 124(4):879-926): published record found; no numerical appendix or replication code located online. — [JPE](https://www.journals.uchicago.edu/doi/10.1086/686734); [IDEAS](https://ideas.repec.org/a/ucp/jpolec/doi10.1086-686734.html)
- Project context: the solver ports Bocola's residual_model.m + parsolve.m (damped FD Newton) with closed-form mu, and the project's own measurement found Bocola's fitted-mu policy gives a 28.4 bp liquidity premium where the exact multiplier gives 2.1 bp — the same failure (5). — `/Users/Huawei/Desktop/Quantitative_Model/CLAUDE.md`

### Inferences
- The project's own finding that Bocola's published number is the fitted-mu read suggests Bocola approximated a kinked object globally and accepted the error; i.e. his code is not a model for curing (5).
- Garcia-Zangwill (mu = max(x,0)^k, slack = max(-x,0)^k) gives one smooth unknown x across regimes; the project already uses it for the TPI floor. Using it for the leverage constraint with x as the stored rule makes the stored object linear through the switch (x crosses zero smoothly) — this is a cheap alternative to PEA if changing the stored rules is hard. Effort minimal.
- Always-binding gives negative multipliers at a third of nodes: this is direct evidence the constraint is genuinely slack over much of the box, so any method assuming binding is ruled out; eps-smoothing at 2e-2 likely pushes mu enough above zero to move the rest point, which interacts with (4).

### Gaps
- Bocola's and GKP's actual numerical handling of the constraint could not be verified from primary sources.

---

Ranked by payoff per effort (inference, synthesised from above):
1. PEA-style stored object (approximate the smooth expectation in alpha, recover mu and alpha in closed form) or Garcia-Zangwill single variable: minimal-moderate; targets (5), (1), (6).
2. Overdetermined least squares on the same basis (denser evaluation points, existing LM/Gauss-Newton via scipy least_squares): minimal; targets (1), softens (2).
3. Anderson-accelerated TI, used as the equilibrium selector and as Newton's start: minimal; targets (3), speed.
4. PCA/ergodic domain (JMMV parallelotope; MM 2015 EDS grid with regression): minimal-moderate; targets (2), (4), accuracy.
5. Pseudo-arclength continuation with sigma_min tracking: moderate; diagnostic for folds/multiplicity (1),(4).
6. Piecewise (two-element) Chebyshev in s split at the switch: moderate; (1),(5).
7. Adaptive local sparse grids (Brumm-Scheidegger) or deep learning (MMW, DEQN): rewrite.
