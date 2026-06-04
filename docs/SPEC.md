# SPEC — Two-Country Monetary-Union HANK-DSGE with Sovereign Default

## Purpose
Study **doom-loop dynamics** (bank ↔ sovereign feedback) in a monetary union and the
**ECB TPI** (Transmission Protection Instrument) as a doom-loop breaker. PhD project 1.
Overleaf: https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c

## Framework
- `sequence-jacobian` 1.0.0 (shade-econ). Linearised transition: **23 unknowns / 23 targets** per `model_v12.ipynb`.
- Python env: `/opt/anaconda3/envs/ssj/bin/python` (base env has broken liblapack symlink).

## Blocks (per country D / F + global)

### Households (HET)
- Heterogeneous, Rouwenhorst / GMAR income process (matrices in `Discretisation/`).
- Hold **deposits only** with domestic intermediary. GHH preferences.
- Progressive income tax `z = λ·y^(1−τ)`.
- `income_D`: `y_pre = (w·N·e + div)/P_CES`. Default hits HH only via labour + dividend.
- Habit formation disabled (`habit_D = 0`).

### Gertler–Karadi intermediary
- Multi-asset incentive constraint over capital, domestic LT sovereign bonds, cross-border bonds.
- Three ν's per country: `ν_K`, `ν_bD`, `ν_bF`.
- **Auclert (2019) capital adjustment cost `Phi`**: `arg = (K − (1+r^k)·K(-1))/(K(-1)+chi0)`;
  `Phi = (chi1/chi2)·(arg²)^(chi2/2)·(K(-1)+chi0)`. `chi0` regularises div-zero at SS; `chi2=2`.
  Enters net-worth LoM in `intermediation_P2` directly (NOT via `rn`).
- **Macroprudential bond tax** `tau_mp = T0 + T1·def_rate` (replaced old `risk_weight × mp_wedge`).
  `T = (T0+T1·def_rate)·(b_D_D+b_F_D)`. Marginal rate → bond Euler FOC wedge; total `T` → net-worth drain.
  Currently inactive (`T0=T1=0`).
- **Market-value balance sheets** throughout: `kappa=Q·K/n`, `phi_b=q_b·b/n`; IC: `Q·K+q_b·b_D+q_b·b_F=θ·n`.
- Bank releveraging `rho_theta` removed → full re-leveraging `ic_res = theta − theta_tgt`.

### Production
- Cobb-Douglas in begin-of-period capital `K(−1)` and labour `N`. TFP shocks `Z`.
- Capital producers with `Φ(I/K)` convex adj cost pinning Tobin's `Q`.
- `firm_profit_D = (1−mc_D)·Y_D` (no Rotemberg on main).

### Government
- Calvo perpetual bond, decay δ_b ≈ 0.05 (~5y avg maturity).
- Endogenous default rate: smooth power function of debt-gap. Tax rule on outstanding debt.
- `zeta_writeoff=1`: default = transfer bank→gov; `zeta=0`: gov pays full face value.

### Trade / external
- CES home/foreign bundle, terms-of-trade `p` (pinned by `goods_mkt_D`, NOT a target on this branch).
- Cross-border bonds via external-account / NFA block.

### Bond pricing
- `q_b` is the **primary** variable (forward-looking PV). `rb = 1/q_b − 1` derived as
  interpretation-only (`bond_yield`, NOT consumed by any residual). `rb_actual` = realised return after default haircut (used).

### ECB TPI (headline policy)
- Buys distressed sovereign bonds when spreads breach trigger. Block exists in `equations_global.py`.

## Key files
| File | Role |
|---|---|
| `model_v12.ipynb` | Main: calibration → SS solve → Jacobian → IRFs → Walras probe. |
| `equations_D.py` / `equations_F.py` | Country-block formulas (`_F` strict mirror of `_D`). |
| `equations_global.py` | Trade, NFA, global bond clearing, portfolio adj cost, ECB TPI. |
| `Discretisation/` | Pre-computed GMAR/Rouwenhorst income matrices. |
| `probe_walras.py` | Walras-leak diagnostic. |
| `diagnose_default_shock.py` | Default-shock IRF diagnostics. |

## Guardrails
- Personal research repo, single user (Adam, PhD). Terse. No trailing summaries.
- No error handling / fallbacks / speculative abstractions. Trust the framework.
- Default no comments; preserve existing equation-block comments.
- `_F` edits = strict mirror of `_D` unless told otherwise.
- Do NOT commit unless asked.
- Do NOT pull wage-NKPC from `add-nkwpc` (PR #16) unless asked.
- Token budget: edit `.py`, never Read full notebook; batch experiments.
