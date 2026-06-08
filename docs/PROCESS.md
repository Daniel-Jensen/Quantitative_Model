# Development Process

## Workflow

1. Use `model_v12.ipynb` as the active development notebook.
2. Keep equations in `equations_D.py`, `equations_F.py`, and `equations_global.py` as the primary model code.
3. Keep `model_notes` as a legacy research log; use the new documentation files for current status and handover.
4. Use notebook hygiene with `nbstripout` and `nbdime` to avoid noisy notebook diffs.

## Notebook hygiene setup

```bash
pip install nbstripout nbdime
nbstripout --install
nbdime config-git --enable
```

## Typical iteration

- edit calibration and equation blocks in the notebook
- re-run the steady-state solve sections
- inspect diagnostic residuals and impulse responses
- verify the corrected EBA bilateral sovereign calibration in `code/model_v12.ipynb`:
  - b_D_D = 0.2447 × assets (Greek sovereigns in Greek banks)
  - b_F_F = 0.2579 × assets (German sovereigns in German banks)
  - b_F_D = 0.0018 × assets (German sovereigns in Greek banks)
  - b_D_F = 0.0065 × assets (Greek sovereigns in German banks)
  - Verify: `b_D_D + b_F_D <= asset_D` and `b_F_F + b_D_F <= asset_F`
  - Verify: `b_D_D + b_D_F == B_supply_D` and `b_F_F + b_F_D == B_supply_F`
- confirm that Greek sovereign shocks hit Greek bank net worth more than German bank net worth, and vice versa
- note that corrected calibration reduces bond supplies 47-63%, changing debt-to-output ratios
- adjust `psi_bF_D`, `psi_bD_F`, and bank leverage if needed given smaller cross-border positions
- commit the cleaned notebook and update documentation

## Key debugging resources

- `diagnose_default_shock.py`: use this for testing default shock responses and solver behavior
- `Discretisation/Outputs`: source data for Markov chain approximations used by household grids
- `model_v11.ipynb`: reference implementation of free bond trade between intermediaries

## Code structure

- Root notebooks: main experiment and calibration environment
- `code/model_v12.ipynb`: active notebook with the bilateral sovereign calibration
- `code/equations_D.py` / `code/equations_F.py`: functions decorated with `sequence_jacobian.simple` and `sequence_jacobian.het`
- `code/equations_global.py`: cross-country clearing and portfolio friction blocks

## Version management

- Keep `model_v12.ipynb` as the leading version.
- Preserve previous versions in `model_v11.ipynb` and `OLD models /` for historical comparison.
- When adding major new features, document the change in `STATE.md` and `PROCESS.md`.
