# Handoff Notes

## Where to start

- Working branch: `main` (all audit fixes merged). Modular entry points in `code/` (`main.py`, `full_model.py`, `steady_state.py`, `tpi.py`); legacy notebook `code/model_v12.ipynb` retained.
- Core equations: `code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py`.
- Read `docs/STATE.md` for current model status and open issues.
- Read `docs/bank_cal_review.md` for the calibration porting roadmap (next major task).

## Quick start

```bash
conda activate ssj   # or use /opt/anaconda3/envs/ssj/bin/python explicitly
cd /path/to/QUANTITATIVE_MODEL
jupyter notebook code/model_v12.ipynb
```

Install deps if needed:
```bash
pip install sequence-jacobian numpy scipy matplotlib nbstripout nbdime
nbstripout --install && nbdime config-git --enable
```

## Current model state (2026-06-22)

`main` is the working baseline. The six structural bugs fixed in the 2026-06-11 audit were merged into `main` via PR #27 (`AB-audit`, 2026-06-11); `main` was then reorganised into modular Python files via PR #28 (2026-06-22). Use `main` for new work; the `audit` / `AB-audit` branches are historical.

**Calibration (main):**
- `phi_lamb_D = phi_lamb_F = 0.15` — min stable at current amplification (Bohn=0.60/yr; over-strong vs literature, see STATE.md)
- `def_scale_D = 0.25` — strong, above GR crisis peak
- `delta_b_D = delta_b_F = 0.10` — 2.5yr duration; needs porting to 0.036/0.038 (7yr/6.5yr)
- `theta_D = theta_F = 4.0`, `psi_lambda_B = 3.0`, `f = 0.12`
- `writeoff_enabled = 0`, `Delta_cross = 1.45` (both open issues, see below)

**Verified post-fix:**
- `ca_res_D ≤ 5.8e−8`, `goods_mkt_F ≤ 8e−10` (Walras holds to solver tolerance)
- Default shock: `n_inter_D[0] = −3.5%` (negative, correct sign)
- TPI: `ca_res_D` γ-invariant ≈5e−8 (CB budget closed)

## Open issues (author decisions required)

**C-1 — degenerate IC constraint:** `Delta_cross = 1.45 > 1`. The multi-asset incentive constraint back-solve (`_ic_delta`, ratio=2.0) produces a divertable fraction above 1. Resolution: hardcode `Delta_D = 0.2, Delta_F = 0.4` per `bank-cal` branch (avoids the back-solve, accepts small SS IC slack absorbed by lambda_gk).

**S-1 — no realized losses:** `writeoff_enabled = 0`. Default shock affects bond prices but produces no balance-sheet write-off. `recovery_rate = 0.40` and `zeta_writeoff = 1.0` are already set on bank-cal; porting them + flipping `writeoff_enabled = 1` gives the doom-loop balance-sheet channel. Alternative: keep writeoff=0 and reframe the paper as a pure risk-premium story.

## Next session priorities

In order from `docs/bank_cal_review.md` §Recommendation:

1. **Port calibration values from bank-cal onto `main`** (not a merge — port values only):
   - `delta_b_D = 0.036`, `delta_b_F = 0.038`
   - `f_D = f_F = 0.03`
   - `Delta_D = 0.2`, `Delta_F = 0.4` (hardcoded, resolves C-1)
   - `recovery_rate = 0.40`, `zeta_writeoff = 1.0`
   - EBA bilateral shares: `b_F_D/asset = 0.0018`, `b_D_F/asset = 0.0065`, `b_D_D/asset = 0.2447`, `b_F_F/asset = 0.2579`
   - Frisch: decide 0.5 (current) vs 1.0 (bank-cal)
   - Income process: consider bank-cal's `rho_z_F=0.97/sigma=0.8` but keep `nZ=15` (bank-cal's 7 is too coarse)

2. **Decide S-1** (writeoff_enabled).

3. **Re-map (phi_lamb, def_scale) bifurcation** on the fixed model at the ported amplification (lower psi_lambda_B, longer duration). The bank-cal bifurcation diagram is invalid post-T-2-fix (accidental stabilizer removed). Target: lowest phi_lamb consistent with stability at a non-trivial doom loop strength. Use `audit_artifacts/philamb_test.py` as the template.

4. **Re-generate all figures** from `main` once calibration is settled.

## Important file locations

| File | Purpose |
|------|---------|
| `code/model_v12.ipynb` | Active notebook |
| `code/equations_D.py` | Country D blocks |
| `code/equations_F.py` | Country F blocks |
| `code/equations_global.py` | Global clearing + portfolio costs |
| `docs/STATE.md` | Current model status, calibration table, open issues |
| `docs/SPEC.md` | Research goals and modelling choices |
| `docs/PROCESS.md` | Workflow, debugging, regression test |
| `docs/audit.md` | Master audit log (ranked findings, fix history) |
| `docs/verification_report.md` | Post-fix verification with numerical evidence |
| `docs/bank_cal_review.md` | bank-cal branch analysis; calibration porting roadmap |
| `docs/walras_forensics.md` | Analytical Walras derivation; all leaks proven |
| `audit_artifacts/run_audit.py` | Full Walras regression pipeline |
| `code/tpi_plots.py`, `code/irf_plots.py` | Figure-generation scripts (regenerate from `main`) |
| Overleaf | https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c |

## Run environment

```
/opt/anaconda3/envs/ssj/bin/python   ← always use this
```

Base env has a broken `liblapack` symlink. Each Jacobian solve at current calibration (T=500) takes ~3 min.
