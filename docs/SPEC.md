# Project Specification

## Goal

Build a tractable two-country general equilibrium model with heterogeneous households, financial intermediaries, sovereign debt, and cross-border portfolio frictions. Primary application: the 2010–2012 Greek sovereign debt crisis and ECB policy (TPI).

## Functional requirements

- Solve a steady-state equilibrium for two countries (D=Greece, F=Germany), each with:
  - Heterogeneous households optimizing deposits and labour supply (GHH preferences, EIS=0.5, Frisch=0.5)
  - Gertler-Karadi financial intermediaries holding domestic and foreign government bonds and productive capital
  - Production with capital accumulation, price of investment, and capital adjustment costs
  - Government fiscal budget with Bohn fiscal rule and Hatchondo-Martinez geometric-decay perpetuity bonds
- Endogenously determine:
  - Sovereign bond prices, yields, and spreads (default risk feedback via `def_scale`)
  - Deposit returns and banking net worth dynamics (GK incentive constraint and P1 Bellman)
  - Terms of trade and net exports (CES consumption basket over D-goods and F-goods)
- Cross-border portfolio adjustment costs (`psi_bF_D`, `psi_bD_F`) anchoring bilateral sovereign positions.
- Impulse response functions (SSJ Jacobian, T=100–500 periods) for:
  - Sovereign default probability shocks
  - TFP shocks
  - TPI (closed-loop CB bond purchase rule)

## Model code structure

- `code/equations_D.py` — country D household (EGM het block), bank, production, capital, government equations
- `code/equations_F.py` — country F analogues
- `code/equations_global.py` — global goods market, external account, bond clearing, portfolio adjustment costs
- `code/model_v12.ipynb` — active development notebook; calibration, SS solve, Jacobian, IRFs, TPI experiments
- `routines/` — auxiliary: grid construction, income Markov chain, Gini calculation
- `audit_artifacts/` — forensic audit reproduction scripts and result logs (see `docs/audit.md`)
- `plots/` — output figures (TPI welfare and spread panels)

## Research objectives

1. Study the interaction between sovereign risk and bank portfolio choice in a two-country monetary union.
2. Quantify how investment and deposit returns propagate through global goods markets (terms-of-trade channel).
3. Assess portfolio adjustment costs and cross-border sovereign exposure (bilateral Greece-Germany via EBA data).
4. Analyse the TPI (Transmission Protection Instrument): welfare and spread effects of a closed-loop CB bond purchase rule.
5. Evaluate the doom-loop mechanism (debt → spread → bank net worth → output → debt) and its dependence on bond duration, bank leverage, and fiscal rule strength.

## Key modelling choices (with rationale)

- **Y = F(K_t):** production uses current-period capital (not K(−1)). Capital producer receives `mpk·(K−K(−1))` to close capital income accounting (W-1 fix). Alternative: K(−1) timing eliminates the term; both are internally consistent.
- **Predetermined deposit rate:** `Rgross = (1+rdep(−1))·P(−1)/P`. Deposit contracts signed at t−1, so funding costs are predetermined — standard NK timing for bank liabilities.
- **Hatchondo-Martinez perpetuity:** bond coupon decays geometrically at rate `1−delta_b`. Duration = 1/delta_b quarters. Captures MTM capital losses on bank balance sheets.
- **GK agency problem:** divertable fraction `Delta` drives the IC constraint binding. Multi-asset IC requires separate `Delta` for each asset class.
- **Walras redundancy:** equations `ca_res_D` and `goods_mkt_F` are dropped from targets. Post-fix they hold to machine tolerance; see `docs/walras_forensics.md`.

## Calibration strategy

Bilateral Greece-Germany calibration targeting:
- EBA (2012) bilateral sovereign exposures: b_F_D/asset≈0.18%, b_D_F/asset≈0.65% (cross-border), b_D_D≈24.47%, b_F_F≈25.79%
- Bond duration: `delta_b_D=0.036`, `delta_b_F=0.038` (Hatchondo-Martinez matching GR/DE 2011 avg maturities ~7yr/6.5yr) — target for next calibration phase; current value 0.10
- Bohn fiscal coefficient: literature 0.025–0.038 quarterly for EA periphery (Staehr 2008); stability requirement on the fixed model must be re-mapped (see `docs/bank_cal_review.md`)
- Spread-debt slope: `def_scale × 430 bps per pp annual B/Y`; crisis peak ~0.12–0.23

## Out of scope (current phase)

- Formal welfare analysis with non-linear transition paths (Jacobian-linearised only)
- External habit formation
- Macroprudential capital requirements (conceptual only, not implemented)
- ECB/ESM backstop (separate from TPI; not implemented)
- Estimation / Bayesian identification
