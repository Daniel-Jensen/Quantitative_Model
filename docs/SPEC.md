# Project Specification

## Goal

Build a tractable two-country general equilibrium model with financial intermediaries, sovereign debt, and cross-border portfolio frictions.

## Functional requirements

- Solve a steady-state equilibrium for two countries, each with:
  - households optimizing deposits and labour
  - financial intermediaries holding domestic and foreign government bonds
  - production, capital, and government fiscal budgets
- Endogenously determine:
  - sovereign bond prices and spreads
  - deposit returns and banking profits
  - trade prices and net exports
- Include frictional portfolio preferences for cross-border bond holdings.
- Produce impulse response functions for shocks to sovereign spreads, portfolio costs, and policy variables.

## Model components

- `equations_D.py`: domestic household, bank, production, and government equations
- `equations_F.py`: foreign-country household, bank, production, and government equations
- `equations_global.py`: global goods market, external account, bond clearing, and portfolio adjustment costs
- `model_v12.ipynb`: active notebook tracking the latest model version
- `model_v11.ipynb`: earlier version with free bond trade across intermediaries

## Research objectives

- Study the interaction between sovereign risk and bank portfolio choice in a two-country setting.
- Evaluate how investment and deposit returns propagate through global goods markets.
- Assess the role of portfolio adjustment costs for cross-border bond holdings.
- Explore policy extensions such as ECB-style backstops and macroprudential capital requirements.

## Out of scope for the current phase

- Full welfare analysis and estimation
- Fully calibrated long-term bond maturity structure
- External habit preferences unless integrated with stable calibration
- Formal policy experiment library beyond IRF-driven exploration
