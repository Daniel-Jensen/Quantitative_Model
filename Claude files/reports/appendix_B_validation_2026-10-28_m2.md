# Appendix-B global model, rest-point calibration (M2): why the D bank cannot bind at rest (2026-10-28)

Everything here runs on the clean, converged coarse solve (`runs/m1_calm_hard_mu_eps_book_eps_floor_eps.pkl`:
genuine IC max, hard books, hard floor, 2.5e-10, no crisis node). The drivers are in `diagnostics/tex_calibration/`:
`m2_warm.py` (warm-started evaluations), `probe_f.py`, `check_calm_port.py`, `m3_fine.py`.

## The targets at the rest point, before calibration

| target | rest point | target |
|---|---|---|
| leverage D / F | 7.81 / 8.04 | 8 / 8 |
| loan spread D / F | 0 bp (IC slack) / 13.6 bp | 5 / 5 bp |
| own-sovereign share D / F | 15.6% / 6.50% | 14.4% / 6.5% |

## The D spread instrument cannot bind the D bank at rest

The bind phase of `m2_warm.py` doubled the SS spread instrument, each step a warm hard solve that converged:

| SS spread instrument (bp/yr) | 5 | 20 | 40 | 80 | 160 | 320 |
|---|---|---|---|---|---|---|
| SS mu_D | 0.0010 | 0.0040 | 0.0079 | 0.0157 | 0.0309 | SS infeasible ("fold-blocked: raise f") |
| rest slack ratio E[Ω]R·n/assets (1 = binding) | 1.123 | 1.125 | 1.136 | 1.176 | 1.234 | — |
| rest E[Ω] | 1.104 | 1.148 | 1.220 | 1.409 | 1.953 | — |

- A tighter SS constraint raises the rest-point franchise value faster than it tightens the constraint.
- The rest point's `E[Ω]` (banks' value of a euro of equity) sits 10–20% above the SS franchise value. Today's default
  probability at rest is only 0.1%, but the logit risk process is volatile (ρ_s 0.9, σ_s 0.41), so high-risk quarters
  are likely ahead, and equity is scarce and valuable then.
- The franchise value is forward-looking with weight ~(1 − f), and its sensitivity to the discounted future is ~1/f.

## The banker exit rate f_D works, but only at implausible values

`probe_f.py`, at the 5 bp instrument, all solves converged (exact reads):

| f_D | 0.08 | 0.12 | 0.16 | 0.24 |
|---|---|---|---|---|
| rest slack ratio | 1.123 | 1.099 | 1.080 | 1.054 |
| slack ratio at p^d 2% | 1.112 | 1.084 | 1.060 | 1.027 |
| E[Ω] at rest | 1.104 | 1.084 | 1.066 | 1.048 |

- Binding at rest needs f_D ≈ 0.35–0.4, a banker horizon under a year (Bocola: 0.035).
- Even at f_D 0.24, the 2% shock leaves the bank slack, so the loan spread r_wc − r stays at 0.

## Bocola's lending spread moves, but by little

Bocola's lending spread is the expected excess return on capital, E[r^k′] − r_dep. Even with μ = 0, it carries the
risk premium −Cov(Ω, r^k)/E[Ω]:

| p^d | rest (0.1%) | 0.5% | 1% | 2% |
|---|---|---|---|---|
| E[r^k] − r, f_D 0.08 (bp/yr) | 0.4 | 1.4 | 2.3 | 3.7 |
| E[r^k] − r, f_D 0.24 (bp/yr) | 0.4 | 1.7 | 2.9 | 5.1 |
| r_wc − r (the model's reported spread) | 0 | 0 | 0 | 0 |

The sign is right; the size is a few bp at the 2% shock.

## Status

The calibration loop cannot reach its six targets at appendix B: no instrument value makes the D bank bind at rest.
This is the plan's Gate 5, and the choice is the author's. The fine check (M3) runs meanwhile, on the uncalibrated clean
solution.

## The fine grid (M3) overturns the coarse conclusion

`m3_fine.py` re-solved the clean coarse solution on Bocola's fine grid (9 nodes in s, cross terms, 75 points, 3000
unknowns), through `solve_recursive(base=...)` with the calm guard. It converged at every rung, including the hard polish
(2.9e-10). Fitted and exact reads now agree.

| exact read, p^d 2% | coarse (3 s nodes) | fine (9 s nodes) |
|---|---|---|
| Q_bD | −9.28% | −6.78% |
| n_D | −11.3% | −5.40% |
| n_F | −0.14% | −0.25% |
| Y_D | −0.011% | +0.031% |
| C_D | +0.082% | −0.005% |
| E[Ω] | +9.8% | −2.9% |
| lending spread r_wc − r | 0 | 0 |

**The D constraint on the fine grid, at the appendix's f_D = 0.08 (exact reads):**

| p^d | rest (0.1%) | 0.5% | 1% | 2% | 3% |
|---|---|---|---|---|---|
| slack ratio | 1.116 | 1.080 | 1.059 | 1.033 | 1.015 |
| Bocola's spread E[r^k] − r (bp/yr) | 0.4 | 2.0 | 4.1 | 8.3 | 12.6 |
| rdep_D (bp/yr) | 120.0 | 120.7 | 119.6 | 116.4 | 112.8 |

- The coarse grid's rise of E[Ω] with risk, and with it the "f_D ≈ 0.35 needed" conclusion, were artefacts of the
  3-node s axis.
- On the fine grid the D bank moves toward binding as risk rises, and by extrapolation binds at about **p^d 4%**: the
  occasionally-binding design (decided 2026-10-28) is met at the appendix's own f_D = 0.08.
- The coarse M2 table above (spread instrument and f_D) should be re-read on the fine grid before any use.
- **Y_D rises slightly** at p^d ≤ 3% while the IC is slack. With μ = 0 the working-capital rate equals the deposit rate,
  which falls with risk (rdep_D 120 → 113 bp at 3%), so labour demand rises. The output channel needs the IC to bind.

**Running:** `m2_ob.py --fine` extends the s box to a 12% p^d top, locates the threshold, and reads impact rows at
p^d 2/4/6/9.1%.

## The extended fine box: the D constraint binds from p^d 1.3% (2026-10-29)

`m2_ob.py --fine` raised the s box's top on the fine grid from 3% to 6%, then to 12% p^d (the 2010–11 average is 9.1%).
Both rungs converged (2.3e-12; 1.6e-9) through the hard conditions. Runs: `runs/m2_ob_fine_box6.pkl`, `m2_ob_fine_box12.pkl`, `m2_ob_fine_f.pkl`.

- **Guard bug (fixed):** the crisis guard first broke this. At high p^d, low D equity is the genuine solution. The old
  guard polished anyway, sent converged nodes to another root (alpha_D 40) and crashed on NaNs. A first fix still
  switched a p^d-12% node to a centre-seeded root at 138% of n_ss, which cascaded.
- **The guard now:**
  - flags only normal-times nodes, p^d ≤ 1% (`CRISIS_PD_MAX`);
  - polishes only after a real switch;
  - refuses polish moves larger than 0.05 in fit space;
  - turns non-finite node residuals into a failed solve.
- It is in both `recursive_experiment` and `branch_select.py`; `check_calm_port.py` and the fast tests pass.

**Rest point (f_D 0.08, deterministic-SS instruments).** Leverage D/F 7.94/8.03, F spread 12.5 bp, exposures D/F
15.5%/6.50%, sovereign spread 39 bp; IC slack (μ_D 0). Inside the box (b_DD at u +0.78).

**The D IC binds from p^d = 1.3%** (other states at the rest point). Impact rows, exact reads (fitted within ~0.2 pp):

| p^d | 2% | 4% | 6% | 9.1% |
|---|---|---|---|---|
| Q_bD | −10.7% | −17.0% | −21.7% | −27.5% |
| n_D | −20.9% | −32.4% | −39.7% | −48.2% |
| n_F | +0.29% | +0.39% | +0.37% | +0.34% |
| Y_D | −0.30% | −0.58% | −0.71% | −0.82% |
| C_D | −0.10% | −0.33% | −0.36% | −0.37% |
| lending spread r_wc − r | +452 bp | +972 bp | +1204 bp | +1410 bp |
| μ_D | 0.099 | 0.215 | 0.273 | 0.336 |

- Bocola's mechanism, end to end: bond prices fall, D bank equity falls, the IC binds, the lending spread rises, and
  output and consumption fall.
- Every acceptance sign is right **except n_F**, which rises 0.3–0.4%.
- **Magnitudes are large.** Appendix B's transmission row puts the GR-DE overdraft-spread rise at 108 bp at p^d 2% (an
  upper bound); the model gives 452 bp.
- **Two crisis targets conflict:**
  - "binding from p^d ≈ 4%" (the 2026-10-28 choice) means slack at 2%, so 0 bp there;
  - matching ≤ 108 bp at 2% needs a threshold just under 2%.
  - The threshold rises as f_D falls (toward Bocola's 0.035), since a smaller f makes the franchise value more forward-looking.
- **Author decision pending:** which crisis-state target pins f_D.

## f_D against the 108 bp transmission row (2026-10-29 decision), fine grid, 12% box top

`m2_ob.py --target spread2 --search`. The moment is the lending-spread rise r_wc − r at p^d 2% from the rest point (exact
reads). Each f step is a warm fine solve, capped at 0.02 per step: one jump of 0.09 (0.06 → 0.152) stalled.

| f_D | 0.06 | 0.08 (appendix) | 0.10 | 0.12 |
|---|---|---|---|---|
| IC binds from p^d | 1.27% | 1.30% | 1.36% | — |
| spread rise at p^d 2% | 547 bp | 452 bp | 380 bp | hard solve stalls at 9.5e-5 (line search shrinking) |

- **Direction:** on this grid a higher f lowers the response (the coarse probe suggested the opposite); a lower f, toward
  Bocola's 0.035, raises it.
- **Size:** each +0.02 removes less (−95, then −71 bp). Extrapolated, f_D alone bottoms out near 150–200 bp, above the 108 bp
  bound. The threshold barely moves (1.27–1.36%), so the 2% shock always lands well inside the binding region.
- **Why the response is large:** the D bank loses ~21% of its equity at 2% (bond price −10.7% on a 15.5% book at leverage
  ~8). That is data-driven (exposure 14.4%, duration 4.4 years), and it pushes μ to ~0.1, so λ_K·μ/E[Ω] ≈ 450 bp/yr.
- **Status:** stopped after f_D 0.12 stalled. The f 0.10 solution is saved (`runs/m2_ob_fine_f.pkl`). Author decision pending.
