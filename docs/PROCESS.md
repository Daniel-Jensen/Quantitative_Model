# Development Process

## Workflow

1. Use `code/model_v12.ipynb` as the active development notebook.
2. Core model equations live in `code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py`. Edit these first; the notebook imports them.
3. Figures output to `plots/`. Audit reproduction scripts in `audit_artifacts/`.
4. Update `docs/STATE.md` after any calibration or structural change.

## Python environment

Always use `/opt/anaconda3/envs/ssj/bin/python`. The base environment has a broken `liblapack` symlink that causes silent numerical failures.

```bash
/opt/anaconda3/envs/ssj/bin/python script.py
# or in notebook kernel settings: select the 'ssj' conda environment
```

## Notebook hygiene

`.gitattributes` registers `nbstripout` and `nbdime`. Run once per clone:

```bash
pip install nbstripout nbdime
nbstripout --install
nbdime config-git --enable
```

Committed notebooks will have outputs and execution counts stripped automatically.

## Typical iteration

1. Edit calibration parameters in the calibration cell (`f042652d` or equivalent).
2. Re-run the steady-state solve section.
3. Inspect diagnostic residuals: `goods_mkt_D`, `goods_mkt_F`, `ca_res_D`, `deposit_mkt_D/F` — all should be ≤1e-7 at any shock magnitude, or you have a structural problem.
4. Run IRFs; check that `n_inter_D[0]` and `Y_D[0]` fall on a default shock (positive = timing bug).
5. Commit cleaned notebook. Update `docs/STATE.md`.

## Structural regression test

To verify Walras accounting after any equation change, run:

```bash
/opt/anaconda3/envs/ssj/bin/python audit_artifacts/run_audit.py
```

This replicates the full steady-state → Jacobian → IRF pipeline (TFP + default shocks) and prints max residuals for `ca_res_D`, `goods_mkt_F`, `goods_mkt_D`, `deposit_mkt_*`. Reference outputs in `audit_artifacts/run4_kt_log.txt`. Each Jacobian solve is ~3 min on the current calibration.

Acceptance thresholds (from `docs/verification_report.md`):
- `goods_mkt_D` ≤ 1e−14
- `goods_mkt_F` ≤ 1e−7
- `ca_res_D` ≤ 1e−7
- `deposit_mkt_D/F` ≤ 1e−13

## EBA calibration verification

After any change to bank balance-sheet parameters, verify bilateral bond shares in the calibration cell:

```python
# Verify EBA shares
assert b_D_D + b_F_D <= asset_D
assert b_F_F + b_D_F <= asset_F
assert abs(b_D_D + b_D_F - B_supply_D) < 1e-6
assert abs(b_F_F + b_F_D - B_supply_F) < 1e-6
# Cross-border: b_F_D/asset_D ≈ 0.0018, b_D_F/asset_F ≈ 0.0065
```

## Branch convention

- `audit`: working branch with all structural fixes. Use this for new work.
- `main`: pre-fix state preserved for PR diff. Do not commit new model work to main.
- `bank-cal`: calibration branch predating all structural fixes. Do not merge; port calibration VALUES only (see `docs/bank_cal_review.md`).

## Key debugging resources

- `audit_artifacts/run_audit.py` — full regression pipeline
- `audit_artifacts/fix_test.py` — isolated W-1/W-2 Walras repair test
- `audit_artifacts/tpi_test.py` — TPI CB accounting verification
- `audit_artifacts/philamb_test.py` — phi_lamb stability sweep
- `audit_artifacts/bankcal_stability_test.py` — low-amplification probe at bank-cal params
- `docs/walras_forensics.md` — analytical derivation of all Walras leaks
- `docs/audit.md` — master audit log with ranked finding list and fix history

## Version history

- `model_v11.ipynb` (in `OLD models/`): predecessor with free bond trade between intermediaries
- `model_v12`: current version; adds NK labour, portfolio adjustment costs, TPI extension
- Structural fixes applied on `audit` branch (2026-06-11): W-1, W-2, W-3, T-2, A-2, TPI-1 — see `docs/audit.md`
