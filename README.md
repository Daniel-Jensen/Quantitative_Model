Code for first project 

Link to overleaf: 


https://www.overleaf.com/project/698b4f88aeef1d0e1d08cc0c

## Model description

A **two-country monetary union HANK-DSGE** with sovereign default risk, used to study
doom-loop dynamics and the ECB's TPI as a "doom-loop breaker." Built on the
[`sequence-jacobian`](https://github.com/shade-econ/sequence-jacobian) framework.

Headline blocks (per country D/F, plus global):

- **Heterogeneous households** (Rouwenhorst / GMAR income process) holding deposits
  with the domestic intermediary; GHH preferences; progressive income tax
  `z = λ·y^(1−τ)`.
- **Gertler–Karadi intermediary** with a multi-asset incentive constraint over
  capital, domestic long-term sovereign bonds, and cross-border sovereign bonds.
  Auclert-style intermediary capital adjustment cost `Φ(·)`. Optional partial
  adjustment on leverage `θ`.
- **Production**: Cobb-Douglas in begin-of-period capital `K(−1)` and labour `N`;
  TFP shocks `Z`; capital producers with a `Φ(I/K)` convex adjustment cost
  pinning Tobin's `Q`.
- **Government**: Calvo-style perpetual bond (decay δ_b ≈ 0.05 ⇒ ~5y avg
  maturity); endogenous sovereign default rate driven by a smooth power
  function of the debt-gap; tax rule on outstanding debt.
- **Trade and external position**: CES bundle of home/foreign goods, terms-of-trade
  `p`, cross-border bond holdings tracked via an external-account / NFA block.
- **Aggregates**: 23 unknowns / 23 targets in the linearised transition system
  (per `model_v12.ipynb`).

The headline policy block is the **ECB TPI** (Transmission Protection
Instrument) — added in `bdfff9f` — which buys distressed sovereign bonds when
spreads breach a trigger, breaking the bank ↔ sovereign feedback loop.

### Key files

| File | Role |
|---|---|
| `model_v12.ipynb` | Current main notebook: calibration → SS solve → Jacobian → IRFs → Walras probe. |
| `equations_D.py` / `equations_F.py` | All country-block formulas (HET, GK bank, gov, production, market clearing). |
| `equations_global.py` | Trade balance, NFA / external account, global bond clearing, portfolio adj. cost, ECB TPI. |
| `Discretisation/` | Pre-computed GMAR/Rouwenhorst income transition matrices. |
| `solve_chi.py` | Optional: re-estimates portfolio adjustment costs after calibration changes. |

## Timeline

### Single-country phase

| Date / Version | What changed |
|---|---|
| 27 May (v3) | Aggregate SDF added to FI problem; monetary shock; gradual tax adjustment. |
| 28 May | Structural overhaul: HHs only hold deposits; FI holds gov bonds + capital. |
| 30 May (v6) | Progressive taxation; wage rigidity. |
| 1 Jun | Default risk + endogenous labour supply (timing convention unstable). |
| 4 Jun | Reworked FI block; AR(1) risk; wage rigidities stripped (uncertain CB response). |
| v9 | Stable: progressive taxation + FI + purely exogenous default. |

### Two-country phase

| Version | What changed |
|---|---|
| v10 | Two countries trading goods only — first contagion element. |
| v11 | Free trade in bonds between intermediaries; portfolio adjustment cost as main friction. |
| v12 | NK features on the labour side; portfolio adj-cost on FI's foreign-bond deviation from SS; CES consumption bundle. ECB TPI added (`bdfff9f`). |

### Recent Walras-plumbing work (Jun)

| Commit | What |
|---|---|
| `a8129cc` | 8-bug sweep: SS Walras leak 2.7e-6 → 1.2e-7. |
| `f810776` | Added `cap_profit` to resource constraint (later reverted — was masking the Phi-override bug). |
| `35d21fd` | Removed `Phi_D/F = 0` SS override that contradicted `smart_steady`'s real Φ ≈ 1.7e-4. |
| `09ccac0` | Cleanup; `def_scale_D` reverted to 0 after Jacobian crash root cause resolved. |
| `51d9d09` | `cap_profit` routed to HH (not bank); production reverted to `K(−1)` (Adam's earlier convention). Linear Walras leak now at machine precision; remaining ~1e-5 at big shock is quadratic IRF noise. |

## Open questions / planned extensions

- **IRF amplification**: 1% TFP shock currently generates 6% Y / 369% rb_D moves
  because Daniel's `8f190f1` removed the `rho_theta` partial-adjustment on bank
  leverage and `23f35f4` stripped Taylor rule + NKPCs + bank PAC. Restoring
  `rho_theta=0.9 + ksi=0.9` (Adam's `5f34a51` tuning) drops every dynamic
  residual to 1e-7 range — calibration call, not a Walras fix.
- **ECB / ESM backstop**: TPI block exists; needs threshold-based non-linear
  activation rule and counterfactual experiments.
- **Bank capital requirements** (Basel III / macroprudential): scaffolded but
  not active (`T0=T1=0`).
- **External habit formation**: in earlier drafts, currently disabled
  (`habit_D = 0`).
- **Yields**: rb is currently the implied YTM; interpretation work pending.


## Notebook hygiene setup

`.gitattributes` declares `nbstripout` as the clean filter and `nbdime` as the
diff/merge driver for `*.ipynb` files. The filters only take effect after each
collaborator runs the setup commands once per local clone:

```bash
pip install nbstripout nbdime
nbstripout --install        # registers clean filter in .git/config
nbdime config-git --enable  # registers diff/merge drivers in .git/config
```

After setup, committed notebooks will have outputs, execution counts, and
volatile metadata stripped automatically, and `git diff` / `git merge` on
notebooks will use nbdime's cell-aware rendering.
