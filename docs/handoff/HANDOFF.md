# HANDOFF

---

## Handoff: 2026-06-04 13:19 BST

### Current Task State
Session was light: (1) created project docs — `docs/SPEC.md`, `docs/STATE.md`, this handoff;
(2) diagnosed VSCode auto-loading Julia env on restart. No active model coding this session.
Branch `walras-plumbing`, clean. Walras plumbing complete: linear leak at machine precision (`51d9d09`).

### Key Decisions
- Split project context into SPEC (model definition / guardrails) + STATE (current status / open issues) + HANDOFF (resume). SPEC ≈ stable, STATE ≈ volatile.
- VSCode Julia: no settings.json fix exists; only per-workspace extension disable works.

### Modified Files
- `docs/SPEC.md` — new: full model spec, blocks, files, guardrails.
- `docs/STATE.md` — new: current status, Walras trail, open issues, planned extensions.
- `docs/handoff/HANDOFF.md` — new: this file.
- (no `.py` / notebook changes this session)

### Blockers / Open Questions
- SS anchor: does `excess_return_*_ss` need `T0+T1·def_rate_ss` instead of just `T0`? Depends on whether `def_rate_ss=0`. Unverified.
- Cross-border share `phi_bF_D_ss`/`phi_bD_F_ss` mix solved `ss['q_b']` with guess `calibration_start['b']`/`['n_inter']` — verify if `b_F_D`/`b_D_F` are solved by `smart_steady`.
- IRF amplification (6% Y / 369% rb on 1% TFP) — calibration, awaiting decision to restore `rho_theta=0.9 + ksi=0.9`.

### Next Steps
1. Verify SS anchor `T0` vs `T0+T1·def_rate_ss` (DEBUG_BRIEF §6) — print `def_rate_ss`.
2. Fix cross-border SS share inconsistency if `b_F_D`/`b_D_F` solved (DEBUG_BRIEF §7).
3. Decide on IRF amplification: restore `rho_theta`/`ksi` tuning or keep stripped.
4. ECB TPI: add threshold non-linear activation + counterfactual experiments.
5. (Optional) Disable Julia ext per-workspace in VSCode.

### Critical Context
- Python env: `/opt/anaconda3/envs/ssj/bin/python` (base broken liblapack).
- Token budget: edit `.py`, NEVER Read full `model_v12.ipynb`; batch experiments.
- `_F` = strict mirror of `_D`. No commits unless asked. No comments by default.
- Do NOT pull wage-NKPC from `add-nkwpc` (PR #16) unless asked.
- `Phi`/`T` enter `intermediation_P2` not `bank_return`/`rn`.
- `q_b` primary; `rb` derived interpretation-only (not in any residual). `rb_actual` IS used.
- `zeta_writeoff=1` → default = bank→gov transfer.

### Model Summary
- Goal: 2-country monetary-union HANK-DSGE w/ sovereign default; study doom loop + ECB TPI breaker.
- Framework: sequence-jacobian 1.0.0; 23 unknowns/23 targets; main = `model_v12.ipynb`.
- Walras plumbing DONE — linear leak machine precision; residual 1e-5 = quadratic IRF noise.
- Auclert cap-adj `Phi` live; `Phi=0` SS override removed; `cap_profit`→HH; production `K(−1)`.
- Macroprudential tax `tau_mp=T0+T1·def_rate` replaces old friction; inactive (`T0=T1=0`).
- `rho_theta` removed → full re-leveraging. Market-value balance sheets throughout.
- Open: SS anchor T0 vs T0+T1·def_rate_ss; cross-border share mix; IRF amplification (calibration).
- Planned: ECB/ESM backstop activation rule; Basel III cap reqs; habit; yield interpretation.
- Guardrails: terse, no fallbacks, no comments, `_F` mirrors `_D`, no commit unless asked.

### Handoff Context (paste into next session)
Repo: `/Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL`, branch `walras-plumbing` (clean).
Read first: `docs/SPEC.md`, `docs/STATE.md`, `DEBUG_BRIEF.md`, `README.md`. Do NOT Read full notebook.
Python: `/opt/anaconda3/envs/ssj/bin/python`. Equation files: `equations_D.py`, `equations_F.py`, `equations_global.py`.
Walras work is finished (machine precision). Remaining work is verification + calibration, not leak hunting.
Priority probes: (1) print `def_rate_ss` to settle SS-anchor `T0` vs `T0+T1·def_rate_ss` question;
(2) check whether `smart_steady` solves `b_F_D`/`b_D_F` — if yes, switch cross-border SS shares to `ss[...]`.
IRF amplification is a known calibration tradeoff: restoring `rho_theta=0.9 + ksi=0.9` (`5f34a51`) tames it
but that is a modelling call — ask Adam before changing. `_F` strict mirror of `_D`. No commits unless asked.

---
