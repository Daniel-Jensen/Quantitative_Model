# STATE — as of 2026-06-04

## Branch
`walras-plumbing`, clean. Main = `main`.

## What works
- SS solves; **linear Walras leak at machine precision** (`51d9d09`). Remaining ~1e-5 at large shock = quadratic IRF noise, not a leak.
- Auclert intermediary cap-adj `Phi` live; `Phi_D/F=0` SS override removed (`35d21fd`).
- `cap_profit` routed to HH (not bank); production uses `K(−1)` (`51d9d09`).
- Macroprudential tax `tau_mp = T0+T1·def_rate` replaces `risk_weight × mp_wedge` (inactive, `T0=T1=0`).
- `rho_theta` bank-releveraging friction removed; full re-leveraging hardwired.
- Market-value balance sheets throughout.
- `q_b` primary, `rb` derived interpretation-only.

## Walras-plumbing trail (Jun)
| Commit | What |
|---|---|
| `a8129cc` | 8-bug sweep: SS leak 2.7e-6 → 1.2e-7. |
| `f810776` | Added `cap_profit` to resource constraint — later reverted (masked Phi-override bug). |
| `35d21fd` | Removed `Phi_D/F=0` SS override contradicting `smart_steady`'s real Φ≈1.7e-4. |
| `09ccac0` | Cleanup; `def_scale_D` reverted to 0. |
| `51d9d09` | `cap_profit`→HH; production→`K(−1)`. Linear leak at machine precision. |
| `8b46fd7` | Formalised README, added files. |

## Known issues / open
- **IRF amplification**: 1% TFP shock → 6% Y / 369% rb_D moves. Cause: Daniel's `8f190f1`
  removed `rho_theta` partial-adj + `23f35f4` stripped Taylor rule + NKPCs + bank PAC.
  Restoring `rho_theta=0.9 + ksi=0.9` (Adam's `5f34a51`) drops dynamic residuals to 1e-7.
  **Calibration call, not a Walras fix.**
- **SS anchor check**: `excess_return_*_ss` subtracts only `T0`, not `T0+T1·def_rate_ss`.
  OK iff `def_rate_ss=0`; else residual `T1·def_rate_ss` mismatch. Verify (see DEBUG_BRIEF §6).
- **Cross-border share inconsistency**: `phi_bF_D_ss`/`phi_bD_F_ss` mix `ss['q_b']` (solved)
  with `calibration_start['b']`/`['n_inter']` (guesses). Check if `b_F_D`/`b_D_F` solved by
  `smart_steady`; if so switch to `ss[...]` (DEBUG_BRIEF §7).

## Planned extensions
- ECB/ESM backstop: TPI block exists; needs threshold non-linear activation + counterfactuals.
- Bank capital requirements (Basel III): scaffolded, inactive (`T0=T1=0`).
- External habit formation: disabled (`habit_D=0`).
- Yields: `rb` = implied YTM; interpretation work pending.

## Side task this session
VSCode auto-loads Julia env on restart. No `settings.json` flag stops it (extension
`julialang.language-julia` activates on startup). Fix = disable per-workspace via
Extensions UI → gear → Disable (Workspace). Not yet applied.
