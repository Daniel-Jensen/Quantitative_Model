---

## Handoff: 2026-06-08T17:05:17Z (auto-saved before compaction)

### Compaction Metadata
- Trigger: auto
- Transcript: /Users/Adam/.claude/projects/-Users-Adam-Documents-uni-phd-research-QUANTITATIVE-MODEL/34671378-4e1e-49f3-8e25-bc81b0f04612.jsonl
- CWD: /Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL/code

### Git Snapshot
- Branch: bank-cal
- Most recent commit: 26baed2 Recalibrate phi_lamb=0.03 (Bohn=0.12) and update def_scale diagnostic cell

### Model Summary

- Two-country HANK model with GK bank intermediation: D = Greece, F = Germany.
- Core files: `code/model_v12.ipynb` (active), `code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py`.
- SSJ framework (sequence_jacobian): steady state solved via `ha_full.solve_jacobian`, IRFs via `G.apply`.
- EBA bilateral sovereign calibration: b_D_D=24.47%, b_F_F=25.79%, b_F_D=0.18%, b_D_F=0.65% of respective bank assets.
- Portfolio frictions psi_bF_D = psi_bD_F = 0.5 (reduced from 1.5 to match small cross-border positions).
- Steady state validated: IC/P1/labour residuals < 1e-15; goods-market < 1e-7.
- **phi_lamb = 0.03** (annual Bohn coefficient = 0.12, EA periphery range; was 0.40 → Bohn=1.6).
- Endogenous default: `def_scale_D` controls spread-debt feedback; slope ≈ def_scale × 430 bps per pp annual B/Y.
- Bifurcation at def_scale ≈ 0.13 with phi_lamb=0.03 (safe paper scenarios: 0.00, 0.05, 0.10).
- Diagnostic cell (cell `7d7d6cbe`) sweeps three scenarios and checks SS residual, spectral radius (<0.99), IRF smoothness, and range constraints.

### Handoff Context (paste into next session)

Active notebook: `code/model_v12.ipynb` on branch `bank-cal`.

Key calibration state (cell `f042652d`):
- `phi_lamb_D = phi_lamb_F = 0.03` — annual Bohn = 0.12 (EA periphery)
- `def_scale_D = 0.0` at SS (changed per-scenario in diagnostic cell only)
- Bank leverage `theta = 4` (conservative; historical 2011 was 10–25×)
- EBA bond shares applied in cell `96c6bd50`

Diagnostic sweep in cell `7d7d6cbe` (after cell `41ff1316`):
- Scenarios: def_scale ∈ {0.00, 0.05, 0.10}
- Thresholds: spectral radius < 0.99, spread[0] ∈ [0,600] bps, n_inter_D[5] ∈ [-40,+20]%, b_gov peak ∈ [0,15]%
- nw_D[0] > 0 at def_scale=0.05-0.10 is expected (GK forward-looking Bellman, franchise value rises on impact)

Known issues / next steps:
1. theta=4 is 2.5–6× below historical 2011 leverage; consider raising to 10–15 for crisis scenarios.
2. Policy extensions not yet implemented: ECB backstop, macroprudential bond tax.
3. goods_mkt_F not in targets_tp (Walras law; linearisation residual ~7e-4 at impact is acceptable).

Run environment: `/opt/anaconda3/envs/ssj/bin/python`; all scripts must run from `code/` directory.

---
---

## Handoff: 2026-06-09T08:48:09Z (auto-saved before compaction)

### Compaction Metadata
- Trigger: auto
- Custom instructions: (none)
- Transcript: /Users/Adam/.claude/projects/-Users-Adam-Documents-uni-phd-research-QUANTITATIVE-MODEL/34671378-4e1e-49f3-8e25-bc81b0f04612.jsonl
- CWD: /Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL/code

### Last User Message (transcript tail)
(unavailable)

### Last Assistant Message (transcript tail)
(unavailable)

### Git Snapshot
- Branch: bank-cal
- Status:
?? .vscode/
?? docs/handoff/
- Recent commits:
ff3afe1 Calibrate bond duration to historical 2011 avg maturities; add slide figures
26baed2 Recalibrate phi_lamb=0.03 (Bohn=0.12) and update def_scale diagnostic cell
61cf2e8 Add empirical calibration targets for def_scale and phi_lamb
7529353 Document def_scale bifurcation: three regimes in endogenous default feedback
c281346 Validate SS and IRFs; update STATE.md with calibration audit results

### Model Summary
(TODO: fill after compaction — 8–12 bullets)

### Handoff Context (paste into next session)
(TODO: fill after compaction — 10–20 lines of concrete resume instructions)

---
---

## Handoff: 2026-06-09T08:48:12Z (auto-saved before compaction)

### Compaction Metadata
- Trigger: auto
- Custom instructions: (none)
- Transcript: /Users/Adam/.claude/projects/-Users-Adam-Documents-uni-phd-research-QUANTITATIVE-MODEL/34671378-4e1e-49f3-8e25-bc81b0f04612.jsonl
- CWD: /Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL/code

### Git Snapshot
- Branch: bank-cal
- Most recent commit: 17b5ff3 Fix SS convergence for long-duration bonds; update safe def_scale scenarios

### Model Summary

- Two-country HANK model with GK bank intermediation: D = Greece, F = Germany.
- Core files: `code/model_v12.ipynb` (active), `code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py`.
- SSJ framework: SS via `ha.solve_steady_state`; Jacobian via `ha_full.solve_jacobian`; IRFs via `G.apply`.
- EBA bilateral sovereign calibration: b_D_D=24.47%, b_F_F=25.79%, b_F_D=0.18%, b_D_F=0.65% of respective bank assets.
- Bond duration: delta_b_D=0.036, delta_b_F=0.038 (Hatchondo-Martinez perpetuity matching GR/DE 2011 avg maturities ~7yr/6.5yr). q_b_D≈0.529, q_b_F≈0.536.
- SS solved with nDep_F=1500 (raised from 500 to eliminate HA grid quantisation noise in deposit_mkt_F). Beta guesses: beta_D=0.96898, beta_F=0.96817.
- phi_lamb=0.03 (annual Bohn=0.12, EA periphery). Paper scenarios: def_scale ∈ {0.00, 0.02, 0.07}.
- Spectral radius landscape is non-monotonic with new delta_b: fragility zone at def_scale≈0.03–0.06 (rho briefly >1). Safe windows: 0.00 (ρ=0.97), 0.02 (ρ=0.98), 0.07 (ρ=0.97).
- IRF (1pp default shock to D, ρ=0.8): spread[0]=+193 bps, n_inter_D[0]=−4.0%, n_inter_F[0]=−0.2%.
- Slide figures saved: `docs/figures/slide_default_shock_main.pdf`, `slide_default_shock_balance.pdf`.

### Handoff Context (paste into next session)

Active notebook: `code/model_v12.ipynb` on branch `bank-cal`.

Key calibration state (cell `f042652d`):
- `delta_b_D=0.036`, `delta_b_F=0.038` — duration 26.8Q/25.3Q matching 2011 GR/DE avg maturities
- `nDep_D=500`, `nDep_F=1500` — F grid raised to eliminate deposit_mkt_F quantisation noise
- `phi_lamb_D = phi_lamb_F = 0.03` — annual Bohn = 0.12 (EA periphery)
- `def_scale_D = 0.0` at SS; varied per-scenario in diagnostic cell only
- Bank leverage `theta = 4` (conservative vs 2011 historical 10–25×)

SS solve cell `e28719a3`:
- Solver: `broyden_custom`
- Beta initial guesses: `beta_D=0.9689841123`, `beta_F=0.9681712547`, `p=0.997170`
- Converges to machine precision with nDep_F=1500

Diagnostic sweep in cell `7d7d6cbe`:
- Scenarios: def_scale ∈ {0.00, 0.02, 0.07} — all ✓ (rho < 0.99)
- Avoid def_scale ∈ [0.03, 0.06] — fragility zone, rho briefly >1 with long-duration bonds

Known issues / next steps:
1. theta=4 is 2.5–6× below historical 2011 leverage; consider raising to 10–15 for crisis scenarios.
2. Policy extensions not yet implemented: ECB backstop, macroprudential bond tax.
3. goods_mkt_F not in targets_tp (Walras law; linearisation residual ~7e-4 at impact is acceptable).

Run environment: `/opt/anaconda3/envs/ssj/bin/python`; all scripts must run from `code/` directory.

---
