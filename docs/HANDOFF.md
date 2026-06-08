# Handoff Notes

## Where to start

- Open `code/model_v12.ipynb` first — this is the current active version.
- Review `code/equations_D.py`, `code/equations_F.py`, and `code/equations_global.py` for the core model structure.
- Use `OLD models /model_v11.ipynb` as a secondary reference for the free-bond-trade setup.

## Quick start

1. Install dependencies:

```bash
pip install sequence-jacobian numpy scipy matplotlib
```

2. Open `model_v12.ipynb`.
3. Run the calibration and steady-state sections.
4. Inspect IRF sections to validate shock responses.

## Important files

- `code/model_v12.ipynb` — active development notebook
- `OLD models /model_v11.ipynb` — previous version with free trade in bank bond holdings
- `code/equations_D.py` / `code/equations_F.py` — country-specific household and bank equations
- `code/equations_global.py` — global clearing, external account, and portfolio adjustment costs
- `code/diagnose_default_shock.py` — diagnostic tool for shocks and default timing
- `Discretisation/Outputs` — auxiliary data for household grid construction

## Core concepts

- `psi_bF_D`, `psi_bD_F`: portfolio adjustment cost intensities for foreign bond holdings
- `b_F_D`, `b_D_F`: cross-border bond stock positions
- `b_D_D`, `b_F_F`: domestic sovereign holdings for Greek and German banks
- `rb_actual_D`, `rb_actual_F`: actual bond yields used in return equations
- `rdep_D`, `rdep_F`: deposit returns for D and F banks
- `p`: relative price / terms of trade between country goods
- Calibration note: the current model uses EBA asset-share targets for greek/german sovereign holdings, with D=Greek banks and F=German banks.

## Known risks and issues

- Solver fragility around portfolio adjustment calibration.
- Not all policy extensions are implemented; ECB backstop and macroprudential rules are currently conceptual.
- `model_notes` is stale and should not be relied upon as the current project plan.

## Handoff advice

- Treat `STATE.md`, `SPEC.md`, `PROCESS.md`, and `HANDOFF.md` as the new project documentation.
- Keep `model_v12.ipynb` as the primary working notebook and record major decisions there.
- If you need the Overleaf link, it is referenced in `model_notes`, but focus on code and markdown docs as the active handover.
