# The 2x2 sovereign-holdings matrix, before and after TPI

Generated 2026-08-19 10:09:49 by `diagnostics/cb_audit/probe_portfolio.py`. `size_F` = 11.696651.

All entries are **aggregate market value in D goods**, `q_b * quantity`, with per-F-capita legs (`b_D_F`, `b_F_F`, `b_gov_F`) scaled by `size_F`. `q_b_D` and `q_b_F` are both D-good prices, so no terms-of-trade conversion enters.

## Clearing identities (must be ~0)

| point | D paper residual | F paper residual |
|---|---|---|
| steady state | +1.110e-16 | -3.227e-16 |
| gamma=0, t=0 | +5.551e-17 | -5.022e-16 |
| gamma=0, t=4 | +8.327e-17 | -5.291e-16 |
| gamma=0, t=20 | +2.220e-16 | -5.317e-16 |
| gamma=2, t=0 | +1.527e-16 | -5.499e-16 |
| gamma=2, t=4 | -8.500e-17 | +7.286e-17 |
| gamma=2, t=20 | -4.207e-17 | -6.791e-16 |
| gamma=5, t=0 | +4.857e-17 | -8.864e-16 |
| gamma=5, t=4 | +1.180e-16 | +8.674e-17 |
| gamma=5, t=20 | +1.145e-16 | -1.130e-15 |
| gamma=10, t=0 | +1.874e-16 | -1.180e-16 |
| gamma=10, t=4 | +1.249e-16 | -3.695e-16 |
| gamma=10, t=20 | +2.429e-17 | -7.702e-16 |

### Steady state (TPI dormant, `cb_buy_ss = 0`)

| holder | D paper | F paper | total | D paper, % of D issue |
|---|---|---|---|---|
| D banks | +0.967383 | +0.007282 | +0.974666 |  87.27% |
| F banks | +0.141118 | +5.596796 | +5.737913 |  12.73% |
| **CB** | +0.000000 | 0.000000 | +0.000000 |   0.00% |
| **total held** | +1.108501 | +5.604078 | +6.712579 | |
| **issued** | +1.108501 | +5.604078 | +6.712579 | |

### Impact of the 1pp default shock, t=0, gamma=0

| holder | D paper | F paper | total | D paper, % of D issue |
|---|---|---|---|---|
| D banks | +0.930292 | +0.007584 | +0.937875 |  87.76% |
| F banks | +0.129756 | +5.625571 | +5.755327 |  12.24% |
| **CB** | +0.000000 | 0.000000 | +0.000000 |   0.00% |
| **total held** | +1.060048 | +5.633155 | +6.693202 | |
| **issued** | +1.060048 | +5.633155 | +6.693202 | |

### Impact of the 1pp default shock, t=0, gamma=2

| holder | D paper | F paper | total | D paper, % of D issue |
|---|---|---|---|---|
| D banks | +0.929302 | +0.008184 | +0.937486 |  87.54% |
| F banks | +0.122907 | +5.626666 | +5.749574 |  11.58% |
| **CB** | +0.009420 | 0.000000 | +0.009420 |   0.89% |
| **total held** | +1.061630 | +5.634850 | +6.696480 | |
| **issued** | +1.061630 | +5.634850 | +6.696480 | |

### Impact of the 1pp default shock, t=0, gamma=5

| holder | D paper | F paper | total | D paper, % of D issue |
|---|---|---|---|---|
| D banks | +0.928030 | +0.008954 | +0.936985 |  87.25% |
| F banks | +0.114123 | +5.627939 | +5.742062 |  10.73% |
| **CB** | +0.021539 | 0.000000 | +0.021539 |   2.02% |
| **total held** | +1.063692 | +5.636893 | +6.700586 | |
| **issued** | +1.063692 | +5.636893 | +6.700586 | |

### Impact of the 1pp default shock, t=0, gamma=10

| holder | D paper | F paper | total | D paper, % of D issue |
|---|---|---|---|---|
| D banks | +0.926341 | +0.009973 | +0.936314 |  86.86% |
| F banks | +0.102520 | +5.629347 | +5.731867 |   9.61% |
| **CB** | +0.037622 | 0.000000 | +0.037622 |   3.53% |
| **total held** | +1.066483 | +5.639320 | +6.705803 | |
| **issued** | +1.066483 | +5.639320 | +6.705803 | |

## What TPI moves: holdings at t=0 relative to gamma=0 (same shock)

Aggregate market value in D goods, and as a % of the D-paper stock at SS.

| leg | SS level | g=0 | g=2 | g=5 | g=10 | d(g=10 - g=0) | as % of SS D issue |
|---|---|---|---|---|---|---|---|
| D banks / D paper | +0.967383 | +0.930292 | +0.929302 | +0.928030 | +0.926341 | -0.003951 | -0.356% |
| F banks / D paper | +0.141118 | +0.129756 | +0.122907 | +0.114123 | +0.102520 | -0.027236 | -2.457% |
| CB / D paper | +0.000000 | +0.000000 | +0.009420 | +0.021539 | +0.037622 | +0.037622 | +3.394% |
| D banks / F paper | +0.007282 | +0.007584 | +0.008184 | +0.008954 | +0.009973 | +0.002389 | +0.216% |
| F banks / F paper | +5.596796 | +5.625571 | +5.626666 | +5.627939 | +5.629347 | +0.003777 | +0.341% |
| D govt issue | +1.108501 | +1.060048 | +1.061630 | +1.063692 | +1.066483 | +0.006435 | +0.581% |
| F govt issue | +5.604078 | +5.633155 | +5.634850 | +5.636893 | +5.639320 | +0.006166 | +0.556% |

## Sovereign concentration `phi = q_b*b / n_inter` — what the IC sees

This is the object `intermediation_IC_D/F` reads, so it is where the portfolio shift becomes a constraint effect.

| ratio | SS | g=0 t0 | g=2 t0 | g=5 t0 | g=10 t0 |
|---|---|---|---|---|---|
| phi_bD_D (D bank, own paper) | 0.452489 | 0.486757 | 0.479736 | 0.470614 | 0.458345 |
| phi_bF_D (D bank, F paper) | 0.003406 | 0.003936 | 0.004167 | 0.004463 | 0.004853 |
| phi_bD_F (F bank, D paper) | 0.007415 | 0.006774 | 0.006428 | 0.005983 | 0.005397 |
| phi_bF_F (F bank, own paper) | 0.294078 | 0.293864 | 0.294447 | 0.295200 | 0.296203 |

