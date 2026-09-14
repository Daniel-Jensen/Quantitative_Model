# Bocola (2016) global-solution replication — derivation and conventions

Standalone replication of Luigi Bocola, "The Pass-Through of Sovereign Risk"
(JPE 2016), Section II.C global projection method. Package:
`code/bocola2016/`. This document records the balanced-growth-path (BGP)
derivation, the pinned conventions where the (unavailable) online appendix
would otherwise settle a detail, and the accuracy/replication results as
each stage lands.

Purpose within the wider project: the two-country model in `code/global/`
cannot reproduce Bocola's *risk channel* with its finite-horizon
representative-branch approximation (documented at length in
`~/.claude/plans/…`; every branch variant makes the crisis-state bank lever
*into* capital, giving an expansionary artifact). Bocola's own model has no
"branch" — it is one recursive competitive equilibrium solved globally, in
which the bad state is the same decision rules evaluated at a different
point. This package builds that reference solution.

## Status by stage

- **Stage 1 — Smolyak library (`smolyak.py`): PASSED.** Nested
  Chebyshev-extrema sparse grid (Krueger-Kubler 2004). Gate tests
  (`tests/test_smolyak.py`): exact point counts (85 at d=6/μ=2, 389 at
  μ=3); on-grid exactness 1e-10; complete-quadratic exactness incl. cross
  terms 1e-9; sustained geometric error decay 1e-1→1e-2→5e-4 across
  μ=2,3,4 on a smooth 6-D exponential; anisotropic level caps; and a
  full-loop Brock-Mirman time-iteration smoke test recovering the analytic
  policy `C=(1−αβ)Y` to 1e-6 on the same infrastructure the real solve
  uses.
- **Stage 2 — Calibration/BGP (`calibration.py`): PASSED (with a
  documented tension, below).** Gate tests (`tests/test_calibration.py`):
  every stated Table-1/2 target reproduced to 1e-12 under the faithful
  reading; av-recursion, closed-form spread, Jermann normalization, and
  labor FOC all internally consistent. **Primary calibration is now
  `delta_mode="standard"`** (user decision, 2026-07-23): sensible
  K/Y=2.66yr and a well-conditioned global solve; the faithful `euler`
  reading is retained as a robustness variant (its near-frictionless
  K/Y=10yr makes the global solve markedly harder — see tension below).
- **Stage 3 — Restricted (no-sovereign-risk) 5-state model: EQUATIONS
  VALIDATED; global solve accurate on the ergodic set.**
  - `period_map.py` (detrended statics + transitions), `expectations.py`
    (Gauss-Hermite quadrature over the two shocks), `time_iteration.py`
    (pointwise least-squares on unit-free ratio residuals; μ closed form;
    outer damping).
  - **BGP rest-point gate (`tests/test_restricted.py`): PASSED** — with
    shocks off the period map reproduces every BGP stock exactly and all
    four Eulers vanish to 1e-9; μ, av recover their BGP values.
  - **Deterministic dynamics validated**: with shocks off the nonlinear
    model sits stably at the BGP (μ~0.001, leverage 5.0, no drift) —
    confirming the equations, not just the static BGP.
  - **First-order perturbation (`perturbation.py`, Klein 2000 / QZ): PASSED
    gate (`tests/test_perturbation.py`)** — exact BGP linearization point,
    Blanchard-Kahn satisfied (5 stable roots; the persistent net-worth root
    is 0.98), stable transition, sign-sensible IRFs. Serves two roles: the
    independent IRF cross-check AND the smooth initial guess for the global
    solve.
  - **Perturbation-seeded global solve (`restricted_main.py`)**: converges
    (conv~3e-6), the deterministic path stays pinned at the BGP (no drift),
    and the **ergodic μ mean≈0.0045, median 0.0005, p99 0.031** — matching
    Bocola's "μ near zero most of the time (his μ^bg≈0.001), spikes in
    stress." A constant-BGP seed instead left ~1/3 of collocation points
    infeasible and the simulation drifted into those wrong-rule regions
    (ergodic μ exploded to ~0.3); the perturbation seed fixes this.
  - **Accuracy (log10 max-abs Euler residual), closed by the
    never-visited-corners argument**: along the ERGODIC PATH — the states
    the model actually visits — mean −3.1, median −3.0, p99 −2.3, max −1.75
    (typical error ~0.1%, worst ~1.8%). The **ergodic path is 99.3%
    strictly interior to the box**, so the box corners — where ~50% of grid
    points are infeasible (μ→1) because the near-unit-root capital and the
    constraint make the pointwise solve blow up there — are essentially
    never visited and their inaccuracy is irrelevant to every simulated
    result. `restricted_main.py` reports and proves this (interior
    fraction). Solution saved to `output/restricted_solution.npz`.
  - **Corner-accuracy investigation (why the box is deliberately WIDE)**:
    for this near-unit-root model (perturbation eigenvalue 0.98, so capital
    wanders widely) you cannot simultaneously avoid infeasible corners AND
    avoid clipping the wide ergodic set with a fixed box. Tried and
    measured: (a) tightening the box to the ergodic cloud → clips the
    wandering capital and makes ergodic behavior WORSE; (b) a PCA-rotated
    grid decorrelating the endogenous states (`rotated_grid.py`,
    `solve_restricted_rotated`) → cuts the infeasible-corner fraction
    48%→16% and gives MACHINE-PRECISION on-grid accuracy (median −6.5), but
    still mildly clips capital so ergodic accuracy is slightly worse than
    the wide box. The **wide axis-aligned box is therefore the primary**
    (best ergodic accuracy, 99.3% interior path); the rotated variant is
    kept available for uniform on-grid accuracy. Uniform machine precision
    everywhere is a genuinely hard property for a near-unit-root
    occasionally-binding model — the accepted standard (Bocola's own
    included) is accuracy on the ergodic set, which is met.
  - **μ matches Bocola**: ergodic median 0.003, p99 0.031 — near zero most
    of the time (his estimated μ^bg ≈ 0.001), spiking to ~0.03 in stress,
    exactly his Figure-1 characterization.
- **Stage 4 — Full 6-state, two-regime model with priced default: MACHINERY
  BUILT AND NESTING-VALIDATED; default solve in progress.**
  - `expectations_full.py`: the genuine two-branch quadrature — three shocks
    (ε_z, ε_g, ε_s) × the two default branches d′∈{0,1} weighted by
    `p^d(s)=logistic(s)` (Bocola eq. 11). The priced default risk enters
    through the d′=1 branch's haircut bond payoff (`surv=1−D`). **This is
    the "distribution, not a point mass" structure whose absence made the
    two-country model's risk channel expansionary.**
  - `time_iteration_full.py`: two-regime time iteration; rules carry a
    coefficient vector per regime `coef[rule]=(coef_d0,coef_d1)`; the
    current regime enters only the current statics (haircut on bonds held
    into a default period); the two-branch continuation is shared
    (re-default allowed).
  - **Nesting gate (`tests/test_full.py`): PASSED** — with `p^d≡0` and
    shocks off, the full machinery reduces EXACTLY to the restricted BGP
    rest-point (residuals <1e-8, stocks/μ/av recovered), and the d=0
    residual equals the restricted residual off-BGP. `p^d(s*)≈0.0009`
    matches Bocola's 0.09% BGP default probability. This proves the 6-state
    grid, 3-shock quadrature, two-branch expectations, and default plumbing
    are all correct.
  - **Default solve: solved (conv 6e-3), and the pass-through is CORRECTLY
    SIGNED.** The d=1 regime post-haircut has crushed net worth so the
    constraint binds hard (μ~0.15–0.4; acceptance threshold raised to 0.85
    to admit these real solutions). Along a d=0 slice through the BGP,
    raising the risk factor s (priced default probability):

    | p^d | Q_b | net worth | leverage |
    |---|---|---|---|
    | 0.09% (BGP) | 0.977 | −0.9% | 5.04 |
    | 1.7% | 0.895 | −3.6% | 5.15 |
    | **2.5% (Bocola expt.)** | **0.876 (−12.4%)** | **−4.3%** | 5.18 |

    **At Bocola's 2.5% experiment the bond price falls −12.4%** (his
    Figure 3: ~−15%), net worth falls, leverage rises — the OPPOSITE of the
    two-country model's expansionary artifact, and the reason for building
    the global solution. This is the central validation that Bocola's
    mechanism, faithfully replicated, transmits sovereign risk
    contractionarily. Caveats: this is a static policy slice at fixed BGP
    capital, so μ≈0 (constraint not yet binding — as in Bocola until losses
    are large) and investment/output are flat; the −2% investment response
    needs the Stage-5 dynamic IRF (capital + net worth evolving, the
    risk-premium covariance term) and a tighter solve (the near-unit-root
    capital margin is the most convergence-sensitive).
- **Stage 5 — Dynamic pass-through IRF (`full_irf.py`): financial block
  replicates; real block needs a tighter solve + the working-capital wedge.**
  A one-time priced s-shock (raising the default probability to 2.5%, never
  realized — Bocola's experiment), traced forward on the no-default path:
  - **Financial pass-through correct**: Q_b −10.3% on impact (Bocola ~−15%),
    net worth −3.4% (his ~−9%), reverting over ~30 quarters as the shock
    decays. Correct sign and reasonable magnitude — the OPPOSITE of the
    two-country artifact.
  - **Real block flat** (I +0.04%, Y +0.008%; his −2%, −0.3%), because
    **μ stays 0 on the path** — the constraint never binds. Two causes,
    both understood: (1) μ~0.001–0.03 is below this solve's resolution
    (conv 6e-3; the near-unit-root capital margin is the least accurate, so
    the barely-binding constraint reads as slack) — needs ~1e-4 accuracy in
    the constraint region; (2) the CLOSED-economy real transmission is weak
    by construction — **Bocola himself (§V.C) flags the closed-economy
    comovement problem and adds a Neumeyer-Perri working-capital wedge in
    his open-economy version (Fig. 7) to get clean investment/output
    declines.** This replication is his closed-economy model, so the same
    applies (and my C response comes out slightly negative vs his +0.3%).
  - **Load-bearing result established**: Bocola's global solution transmits
    sovereign risk contractionarily through the financial block, which is
    exactly what the two-country representative-branch approximation got
    backwards. Getting the real-side IRF magnitudes to match his Figure 5
    is the remaining work: (a) tighten the full solve so the constraint
    binds (level-3 grid / more iterations / the decorrelation machinery
    from Stage 3d applied to the 6-state grid), then (b) add the
    working-capital wedge for the open-economy real transmission.
  - **Solve-tightening attempt (decorrelated 6-state grid, `full_rotated.py`,
    `simulate_full.py`)**: the PCA-rotated grid + best-iterate tracking cut
    the infeasible-corner fraction from 53% to **13%** and strengthened the
    pass-through (Q_b −19% at p^d 2.5%, closer to Bocola's ~−15%). But μ
    still reads 0 where leverage rises above the 5.0 ceiling — internally
    inconsistent, which pinpoints the residual difficulty: the two-regime
    near-unit-root iteration is not globally contractive (it reaches a good
    near-fixed-point ~iter 68 then oscillates away; the solver now returns
    the best iterate), and its update-norm floor (~0.09) is far coarser
    than the ~1e-4 needed to resolve μ~0.001. The HIGH-RISK region — exactly
    where the pass-through matters — is the hardest to converge; the bond
    price is robust (pinned by the bond Euler) but μ and the real-side
    response (capital Euler) are not yet reliably resolved there.
- **Honest state of the replication**: the financial pass-through is
  replicated and correctly signed (the load-bearing result — sovereign risk
  lowers bond prices and net worth, the opposite of the two-country
  artifact). Matching Bocola's Figure-5 real-side magnitudes (investment
  −2%) requires two things, both understood and neither a conceptual gap:
  (1) a globally-convergent solver for the barely-binding constraint in the
  coupled two-regime near-unit-root model (candidates: a proper endogenous
  grid / policy-iteration hybrid, Judd-Maliar-style adaptive Smolyak, or
  Bocola's own Smolyak-with-nonlinear-filter tuning) to conv ~1e-4 so μ
  resolves; and (2) the Neumeyer-Perri working-capital wedge that Bocola
  adds in his open-economy version (Fig. 7) for clean real transmission.
- **Solver work toward global convergence (`newton_collocation.py`)**: built
  a Newton solver on the STACKED collocation residual (quadratic convergence,
  indifferent to the near-unit-root eigenvalue that cripples time iteration),
  with (a) masked anchoring of infeasible corners, and (b) **Fischer-Burmeister
  complementarity smoothing** — μ made an explicit unknown and the `max()`
  kink replaced by `μ + slack − √(μ²+slack²)=0` (the technique the
  two-country model uses). This precisely characterized why the problem is
  hard: there are THREE non-smoothness sources, not one — (1) the
  μ-constraint kink (FB fixes this), (2) an **anchoring floor** (the wide box
  the near-unit-root ergodic set requires has ~40–64% infeasible corners that
  must be anchored at frozen values, which propagate through the global fit
  and floor the achievable residual at ~0.03), and (3) **feasibility guards
  inside the expectation** (`np.where(I'>0,…)` at quadrature nodes) that
  give a noisy FD Jacobian. All three were addressed: FB removes (1); the
  DECORRELATED grid removes (2)'s infeasibility (100% feasible vs 87%
  axis-aligned); a softplus floor on next-period investment inside
  `expect_pieces` removes (3) (BGP still a rest point to 1e-12). With all
  three, **FB-Newton UNSTALLS** — it takes a real step reducing the residual
  7× (0.037→0.0053) where the axis-aligned solve was fully stuck. But it then
  **floors at ~5e-3**, and the smoothed guard did NOT lower this floor, which
  localizes the true residual obstacle: a handful (~16%) of **hard anchored
  points** — transitional-μ / near-insolvency states with a large warm-start
  residual — whose frozen values propagate through the GLOBAL collocation fit
  and floor the residual of every polished point. This is intrinsic to a
  global collocation on a wide box (the corners are hard AND coupled to the
  interior through the fit), and 5e-3 does not beat plain time iteration
  (which resolves the restricted μ to median 0.003 on the ergodic set).
- **Honest conclusion on the solver**: Newton/FB/decorrelation/guard-smoothing
  — the full standard toolkit, matching this repo's own two-country solver —
  were implemented and each removes one obstacle, but the combination floors
  at ~5e-3 on the coupled hard corners, short of the ~1e-6 that resolves
  μ~0.001. Beating this needs an approach that decouples the hard corners
  from the interior fit — e.g. a local/finite-element basis instead of global
  Chebyshev (so corner errors stay local), an equilibrium-selection scheme
  that never places collocation nodes in near-infeasible states, or Bocola's
  own tuned Smolyak+nonlinear-filter pipeline. This is genuine research-grade
  numerical work; the machinery here (Newton + FB + decorrelation + smoothed
  guards, all in `newton_collocation.py`) is the right foundation for it.
- Stage 6 (price/quantity-of-risk decomposition per his Table 4; full
  write-up): pending, and requires the resolved solve above.

## BGP back-out (footnote-14 reparameterization)

Structural parameters `[β, λ, ω, δ, χ, ι, τ*, a1, a2]` are backed out from
BGP targets. With `G = e^γ` the gross tech-growth factor:

- **β** from the household deposit Euler `1 = R·E[Λ′]`, `Λ′ = β·e^{−γ}` on
  BGP ⇒ `β = G/R^bg`.
- **Marginal value of wealth** `av` (his α(S)) from eq. (1) on BGP:
  `av = [(1−ψ)+ψ·av]/(1−μ)` ⇒ `av = (1−ψ)/(1−μ^bg−ψ)`.
- **λ** from the binding leverage identity `lev^bg = av/λ`.
- **Credit spread** (capital/bond Euler on BGP, single λ):
  `R_K−R = λμ/E[Λ̂′] = μR/(lev(1−μ))`. At Bocola's estimates this is
  **8bp annualized** — matching his own text (average liquidity premium
  ≈ the ~8bp interbank spread).
- **δ, K/Y** jointly from `R_K = (1−δ)+αY/K`, the Jermann BGP replacement
  `I/K = G−1+δ` (with `Φ(x*)=x*`), and the `i/y=0.213` target.
- **a1, a2** from the Jermann normalization `Φ′(x*)=1` (Q_K=1) and
  `Φ(x*)=x*` (adj^bg=0).
- **χ** from the intratemporal labor FOC at `l^bg`.
- **ι** from `q_b^bg=1` and the BGP bond Euler (`R_B=R_K` under single λ).
- **B/K** from the exposure target `exp^bg = Q_B B/(Q_K K+Q_B B)`.
- **ω, τ*** (entrant endowment, fiscal intercept) from the net-worth and
  government-budget BGP fixed points — wired in Stage 3 with the full
  period map.

## ⚠ Documented tension: implausibly high capital-output ratio

The faithful ("euler") reading hits **every stated target exactly**, but
implies **δ = 1.65%/yr and K/Y = 10.2 years** — a very capital-intensive
BGP. This is *forced*, not a bug: Bocola's tiny estimated agency friction
(μ^bg=0.001) makes the BGP capital premium only 8bp, and an 8bp user-cost
wedge with a 21% investment share and 30% capital share mechanically
requires a large capital stock. The identities reproduce his numbers; the
implication is what it is.

A `delta_mode="standard"` variant (δ=0.025) gives a sensible K/Y=2.66yr but
then *misses* i/y (0.277 vs 0.213). The two readings bracket the target.

**Resolution taken:** default to the faithful `euler` reading (hits stated
targets; the user's overriding requirement is fidelity to Bocola). The high
K/Y makes capital dynamics sluggish, which could dampen the investment IRF;
**this is the first thing to check when the restricted-model and full-model
IRFs are validated in Stages 3–5.** If capital dynamics come out
implausibly slow relative to his Figure 5 (−2% investment on impact), that
is the signal that the online appendix likely fixes δ at a conventional
value and treats i/y as an approximate moment — a fork to raise with the
user at that point, with concrete IRF evidence rather than a priori.

## Pinned conventions (documented; each sensitivity-checkable)

- **Approximated controls** = the paper's four: `{C̃, R, α(S), Q_B}`,
  separate coefficient vectors for d∈{0,1}. Labor, portfolios, Q_K, μ, N,
  K′, B′ recovered exactly from static conditions — not approximated.
- **Occasionally binding constraint** via the paper's own closed form
  `μ = max{1 − E[Λ̂′]RN/(λ(Q_K K′+Q_B B′)), 0}`; no smoothing.
- **Entrant endowment**, post-haircut GK-standard: with
  `X ≡ [(1−δ)Q_K+Z]K + (1−dD)[π+(1−π)(ι+Q_B)]B`,
  `N = ψ(X−P) + ωX`.
- **Detrending**: period-t flows/prices in lagged-tech units; end-of-period
  stocks scaled by `e^{−Δz_t}` once (K̃′, B̃′, P̃′); SDF factor `e^{−Δz_t}`
  dated t.
- **Default process**: `d′ ~ Bernoulli(logistic(s_t))`, independent of
  current d (re-default allowed).
