# HHBANK 2-Country Macro-Financial Model

A two-country heterogeneous-agent / financial-intermediary model built with `sequence_jacobian`.

## Project overview

This repository implements a macro-financial model with:
- domestic and foreign households optimizing deposits and labour supply
- domestic and foreign financial intermediaries holding government bonds, capital, and deposits
- sovereign debt issuance, default risk, and endogenous bond spreads
- cross-border bond positions and portfolio adjustment costs
- global goods markets and trade balance dynamics
- a modular equation structure in `equations_D.py`, `equations_F.py`, and `equations_global.py`

## Current development

- `model_v12.ipynb` — latest active version, with:
  - New Keynesian labour-side features
  - foreign-bond portfolio adjustment cost
  - two-good consumption basket and trade price dynamics
- `model_v11.ipynb` — previous working version, with free trade in bonds between intermediaries and portfolio friction analysis

## Repository structure

- `code/model_v12.ipynb` — current active notebook
- `OLD models /model_v11.ipynb` — archived previous version
- `code/equations_D.py`, `code/equations_F.py`, `code/equations_global.py` — core model equations
- `code/diagnose_default_shock.py` — diagnostic tools for default and bond shock analysis
- `docs/` — project documentation, including `STATE.md`, `SPEC.md`, `PROCESS.md`, `HANDOFF.md`, and legacy `model_notes`
- `Discretisation/Outputs` — Markov approximation data used by household grids
- `requirements.txt` — reproducible Python dependency list

## Quick start

1. Install required Python packages from the dependency file:

```bash
pip install -r requirements.txt
```

2. Set up notebook hygiene once per clone:

```bash
pip install nbstripout nbdime
nbstripout --install
nbdime config-git --enable
```

3. Open `code/model_v12.ipynb` and run the calibration + steady-state sections first.

## Notes

- `model_v12.ipynb` is the current master notebook for development.
- `model_notes` contains older research notes and a stale log; the new documentation files in the repository should be treated as the current reference.
- The Overleaf project link previously used is documented in `model_notes` but is not the active source of truth for code structure.

