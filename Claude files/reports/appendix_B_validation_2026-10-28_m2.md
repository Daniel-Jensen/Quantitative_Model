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
