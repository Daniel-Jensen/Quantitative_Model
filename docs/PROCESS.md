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
- adjust `psi_bF_D`, `psi_bD_F`, and related parameters if calibration is unstable
- commit the cleaned notebook and update documentation

## Key debugging resources

- `diagnose_default_shock.py`: use this for testing default shock responses and solver behavior
- `Discretisation/Outputs`: source data for Markov chain approximations used by household grids
- `model_v11.ipynb`: reference implementation of free bond trade between intermediaries

## Code structure

- Root notebooks: main experiment and calibration environment
- `equations_D.py` / `equations_F.py`: functions decorated with `sequence_jacobian.simple` and `sequence_jacobian.het`
- `equations_global.py`: cross-country clearing and portfolio friction blocks

## Version management

- Keep `model_v12.ipynb` as the leading version.
- Preserve previous versions in `model_v11.ipynb` and `OLD models /` for historical comparison.
- When adding major new features, document the change in `STATE.md` and `PROCESS.md`.
