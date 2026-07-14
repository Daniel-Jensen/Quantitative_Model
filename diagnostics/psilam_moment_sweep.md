# psi_lambda_B moment sweep (EL_price anchored) — macro-pru-fix

Generated 2026-07-14 09:07:36. `EL_price_D = 0.102491` held FIXED (the empirically-anchored
fundamental expected-loss loading, from recovery/duration). `psi_spread` scales
linearly with `psi_lambda_B`. Moments are the response to a **1pp sovereign default
shock**. Choose `psi_lambda_B` to match an external moment, NOT the old IRFs.

| psi_lambda_B | psi_spread | FOC load (EL+ψs) | peak spread (bp, ann) | peak Δn_inter (%SS) | peak ΔY (%SS) | Δn per 100bp | spread amp vs ψλ=0 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 ← fundamental floor / Case-3 check | 0.000 | 0.102 | 13.1 | -0.089 | -0.0018 | -0.674 | 1.00× |
| 0.50 | 0.130 | 0.232 | 30.8 | -0.247 | -0.0048 | -0.801 | 2.35× |
| 1.00 | 0.259 | 0.362 | 50.2 | -0.421 | -0.0082 | -0.838 | 3.82× |
| 1.50 | 0.389 | 0.491 | 71.8 | -0.613 | -0.0122 | -0.854 | 5.47× |
| 2.00 | 0.518 | 0.621 | 96.5 | -0.831 | -0.0172 | -0.861 | 7.35× |
| 2.60 ← loading-matched to old | 0.674 | 0.777 | 132.7 | -1.141 | -0.0271 | -0.860 | 10.10× |
| 3.00 ← current baseline | 0.778 | 0.880 | 163.6 | -1.390 | -0.0389 | -0.850 | 12.45× |
| 4.00 | 1.037 | 1.139 | 310.6 | -2.139 | -0.1335 | -0.689 | 23.64× |
| 5.00 | 1.296 | 1.399 | 1705.3 | -8.523 | -0.7901 | -0.500 | 129.78× |

## How to read this
- **Disciplining moment = the spread level / amplification** (highly `psi_lambda_B`-sensitive):
  pick `psi_lambda_B` so the peak spread matches the observed Greek response to a
  comparable default-probability move (2010 GR–DE spread ≈ 150bp is the paper's target).
- **Consistency moment = the bank-networth-to-spread pass-through** `Δn per 100bp`
  (roughly `psi_lambda_B`-robust): compare to bank-equity/sovereign-spread event studies
  (e.g. Altavilla–Pagano–Simonelli; Acharya–Drechsler–Schnabl). If it sits in the
  empirical range across the sweep, the doom-loop transmission is the right size.
- **`psi_lambda_B = 0` row is the fundamental floor**: nonzero, correctly-signed response
  from `EL_price` alone — the Case-3 null is resolved independent of the dial.

NOTE: exact empirical target values are a literature-retrieval task (flagged in the
handoff); this table supplies the model side of the mapping so the dial is set by data.
