# Glossary

Symbols, names and labels used in the paper, the code and the project documents. Paper
notation follows the viva paper (Overleaf project 698b4f88aeef1d0e1d08cc0c); code names are
the global model's (`global_projection/`, `calibration/global_projection.py`) unless marked
SSJ. X ∈ {D, F} indexes countries throughout.

## 1. The economy and the two solvers

| Term | Meaning |
|---|---|
| D | The distressed periphery, calibrated to Greece. Only D's sovereign bond is risky. |
| F | The core, calibrated to Germany (in older SSJ documents, sometimes the rest of the euro area). |
| Union, monetary union | D and F with one central bank (the Eurosystem). |
| Global model, global projection | The nonlinear model the viva paper describes: occasionally binding bank constraint, solved by Chebyshev collocation on a Smolyak grid. Code in `global_projection/`. |
| SSJ model, sequence space | The first-order model with sticky prices and nominal deposits, solved with sequence-space Jacobians; used for the linear policy experiments E1-E4. Code in `linear_ssj/`. Calibrated differently; its impulse magnitudes are not comparable with the global model's. |
| Code calibration | The parameter values in `calibration/global_projection.py`, under which every saved global result was produced. |
| Paper calibration | The values in the paper's section 4 and appendix B. Ahead of the code as of 2026-10-27 (see STATE section 1). |

## 2. Rates and spreads (three different objects; always say which)

| Term | Definition |
|---|---|
| Deposit rate r_X (`rdep_X`) | Real rate on household deposits; r̄ = 0.003 a quarter (1.2% a year). |
| Intermediation wedge ς̄ = r̄^k - r̄ | Expected return on bank assets over the deposit rate at the deterministic steady state. Paper 5 bp/yr (one-week Euribor-Eurepo); code 100 bp/yr. Prices the agency friction at rest. |
| Lending (working-capital) spread r^wc - r (`r_wc_X`) | Rate on firms' wage-bill loans over the deposit rate. The only channel from bank balance sheets to labour demand on impact. Moves with mu. |
| Sovereign spread | D bond yield over the safe F bond yield (bp/yr). Prices expected default loss plus the bank-balance-sheet premium. |
| Excess return on bank equity ς̄θ̄ | Wedge times leverage; 0.0010 a quarter in the paper. |
| Credit spread (older documents) | Usually the lending spread; check context. |

## 3. Paper notation and code names (global model)

| Paper | Code | Meaning |
|---|---|---|
| s_t | `s` | Latent sovereign-risk factor, AR(1): s_t = c + ρ_s s_{t-1} + σ_s ε_t. |
| π^d_t = π^d(s_t) | `default_prob(s)` | Priced quarterly default probability, logistic in s; never a function of debt. 0.10% at rest. |
| ρ_s, σ_s | `rho_s` (in the grid builder), `sigma_s` | Persistence and innovation s.d. of s. Paper 0.9, 0.41; code 0.95, 0.4455. |
| d_t, d' | `d` | Default indicator; regimes d' ∈ {0, 1} in every expectation. |
| ϱ_D | `recovery_rate_D` | Recovery rate in default (0.45: the 55% PSI haircut). |
| h_t = 1 - d_t(1 - ϱ_D) | | Survival factor on D-sovereign claims. |
| q^D_t, q^F_t | `Q_bD`, `Q_bF` | Sovereign bond prices. |
| δ^b_X | `delta_b_X` | Coupon decay of the long bond; duration ≈ 1/δ^b quarters. |
| B̄_X | `B_gov_X_ss` | Steady-state bank-held sovereign stock (banks are the only holders). |
| b^D_F | `b_FD` / `b_F_D_ss` | F banks' holding of the D bond (the foreign share of D debt). |
| θ_X = 𝒜_X / n_X | `leverage_target_X` | Bank leverage: assets over net worth. |
| n_X | `n_X` | Bank net worth. |
| μ_X | `mu_X` | Multiplier on the bank's incentive constraint; μ = max{·, 0} is kinked at zero. μ̄ is its deterministic steady-state value. |
| λ_X | `lambda_*_X` | Divertible fraction of bank assets (one value for every asset in the paper). |
| φ_X | `alpha_X` | Marginal value of a unit of bank net worth (franchise value); V = φ n. |
| Ω_{X,t+1} | `Omega_X` | Bank's stochastic discount factor weighted by the marginal value of wealth next period. |
| f_X | `f_X` | Payout share of bankers (exit probability); horizon 1/f quarters. |
| ω_X | `omega_ent_X` | Transfer to entering bankers. |
| ψ | `psi_bF_D`, `psi_bD_F` | Cost on cross-border bond holdings; pins how each bond is split between the two banks. |
| ϰ = 𝓜_F / 𝓜_D | `size_F` | Size of F relative to D. Paper 11.3, code 8. |
| p_t | `p` | Relative price of the F good (terms of trade); the D good is the numeraire. |
| ϖ_X | `omega_home_X` | Home bias in consumption. |
| η | `epsilon_trade` | Trade elasticity between D and F goods. |
| ζ_X | `zeta_wc_X` | Share of the wage bill firms pre-finance with bank loans. |
| ξ_X | `ksi_X` | Elasticity of the Jermann capital-adjustment cost. |
| γ_τ, 4γ_τ | `debt_root_X` (stated as the root) | Tax response to debt; 4γ_τ is the annual Bohn coefficient. |
| κ^cb_D | `tpi_key_D` | D's share of Eurosystem profit and loss (capital key 0.071). |
| m^cb_t | `m_cb` | Eurosystem purchases of the D bond. |
| M, O | `M_cb`, `O_cb` | Eurosystem holdings of D bonds; its obligation to D banks (reserves). |
| | `x_cb` | One variable that carries the purchase complementarity: its positive part is purchases, its negative part the price premium over the floor. |
| | `K_X`, `Kp_X`, `N_X`, `A_X`, `P_X`, `V_dep`, `Z_D` | Capital, next-period capital, hours, household savings, banks' gross deposit obligations, the cross-border deposit position, D productivity. |

## 4. The backstop and policy terms

| Term | Meaning |
|---|---|
| TPI | Transmission Protection Instrument (ECB, July 2022). The paper's modelled instrument: a promise to buy D bonds that caps the spread. |
| OMT, SMP, LTRO, PSI | Outright Monetary Transactions (2012); Securities Markets Programme (2010); Longer-Term Refinancing Operations (the global model's former backstop, deleted 2026-10-01); Private Sector Involvement, the March 2012 Greek debt exchange. |
| Spread cap | The backstop buys whatever holds the D spread at or below the cap (200 bp/yr in the reference run). |
| Warranted floor | A price floor at the default-warranted value of the D bond (its value to a risk-neutral bank with the same funding costs), plus tau bp. |
| Level rule, change rule | Holdings targets θ = γ_L (sp - s̄) and γ_Δ (sp - sp_lag): respond to the spread's level or to its change. Combined rules use both. |
| Capital key | Eurosystem profit shared by ECB capital shares (D receives 7.1%). |
| kappa_D = 1 | Counterfactual in which D receives all profit and loss on its own bonds. |
| Announcement effect | Difference between the economies with and without the backstop at their stochastic rest points, where the backstop exists but buys nothing. |
| Intervention effect | Difference in their responses to the same rise in default risk. |
| Headline shock | Default probability raised from 0.10% to 2% a quarter, then mean reversion at ρ_s. |
| Transfer problem | The backstop's profit on the no-default path flows to F by the capital key, so D must run a trade surplus through a real depreciation. |

## 5. Solution and numerical terms

| Term | Meaning |
|---|---|
| Deterministic steady state | The model with no shocks and no priced risk. Not where impulse responses start. |
| Stochastic rest point | The point where the economy settles when shocks are expected but none are realised. Global impulse responses start here and are differenced against an unshocked path. |
| Kink | The point where μ hits zero. Polynomial rules approximate a kinked function poorly, so results near it are numerically fragile. |
| Fitted rules vs exact clearing | Fitted: every variable read off the polynomial decision rules. Exact clearing: the current quarter's equilibrium solved exactly, with only expectations taken from the fitted rules. Backstop and welfare numbers are quoted on exact clearing. |
| Euler error | log10 of the relative error in an equilibrium condition between grid points; -4 means one part in 10,000. |
| Smolyak grid, s-refinement | Sparse grid over the states; the risk factor s gets extra nodes (5 by default, 9 as a check). |
| Corner | The backstop buys the D bank's whole book; the cap is met only by exhausting purchases. |
| Fischer-Burmeister / Garcia-Zangwill form | Two smooth ways of writing the purchase complementarity; the solver roots in the first and stores rules in the second. |
| Household-fix run | The 6 Oct 2026 global solution (households receive the working-capital flow); every earlier impulse-response level is void. |

## 6. SSJ model names

| Name | Meaning |
|---|---|
| `psi_lambda_B` | Collateral-friction amplification dial: made divertibility depend on default risk. Zero and scheduled for deletion. |
| `Delta` | Divertibility of a sovereign bond relative to capital (0.20: sovereign paper is a fifth as divertible, hence better collateral). |
| `zeta_writeoff`, `writeoff_enabled` | Share of principal written off in default (1); whether the impulse response realises the default (0: it traces the no-default branch). |
| `EL_load_D` | Endogenous expected-loss loading (0.7014). |
| `BANK_SCOPE="broad"` | Banks intermediate the whole capital stock. |
| `phi_lamb`, `mv_rule` | Bohn coefficient of the SSJ tax rule; whether the rule targets market (1) or par (0) value of debt. |
| TPI loading | Eurosystem compensation per unit of expected loss; below one means no over-compensation. |
| E1-E4 | SSJ policy experiments: backstop schedule; output decomposition; write-off variant; distributional incidence. |
| Regimes (dove, base, hawk) | Markov-switching states of the ECB's intervention stance, dated from peripheral spreads. |
| Live Claims | Numbered claims in `SPEC.md` that survived scrutiny. Reported, never targeted. |

## 7. Issue codes in the archive

| Code | Meaning (all resolved unless stated) |
|---|---|
| C-1 | Back-solved divertible fraction above one; the bank block was silently degenerate. |
| S-1 | Whether default is realised along the impulse response (framing: it is not). |
| GK-1 | Well-posedness of the bank's collateral condition at the measured EBA moments. |
| GK-2 | Scope of bank net worth; explosive dynamics fixed by `BANK_SCOPE="broad"`. |
| EL-1 | Recovery rate in the SSJ model (set to 0.30). |
| PT-1 | Pass-through of bank net worth to spreads, checked against event studies. |
| F-1 | Market-value fiscal rule: a hard break at `phi_lamb` 0.15 (open as a caveat). |

## 8. Data and institutions

| Term | Meaning |
|---|---|
| EBA 2011 | European Banking Authority EU-wide stress test, published July 2011; bank-level sovereign exposures at 31 December 2010. The source for 14.4%, 6.5% and 12.7%. |
| CT1 | Core Tier 1 capital, the stress test's narrow capital measure. |
| Euribor-Eurepo | Unsecured minus secured one-week interbank rate; the paper's measure of the wedge. |
| MIR | ECB MFI interest-rate statistics (bank lending and deposit rates). |
| LSEG | London Stock Exchange Group data (formerly Refinitiv); source of the CDS series. Licensed; never committed. |
| Bocola (2016) | Bocola, "The pass-through of sovereign risk", JPE 124(4). The global model's template. |
| BFT | Bi, Foerster and Traum (FRBSF WP 2025-10), the closest published analogue. |
| GK | Gertler and Karadi (2011), the bank block. |
| GHH | Greenwood-Hercowitz-Huffman preferences (no wealth effect on hours). |
| HFCS | ECB Household Finance and Consumption Survey. |

## 9. Branches

| Branch | Role |
|---|---|
| `main` | Merge target. |
| `shelved/policy-rules-welfare` | Local only. The welfare section, the policy-rule economies and the warranted floor, shelved on 2026-10-27 (`Claude files/docs/shelved/welfare_and_policy_rules.md`). |
| `cross-border-empirics` | Shelved; commit `c2d2337` holds empirical code not yet on the working branch (see `DATA.md`). |
| `bocola-rewrite`, `OMT-fix`, `hh-budget-fix` | Earlier global-model branches, merged into `main` (PRs #32, #35, #36). |
| `add-nkpc`, `gk-structural-foc` | SSJ branches (sticky prices; the structural refactor), merged into `main` (PR #33 for the refactor). |
| `bank-cal` | Old calibration branch. Do not merge; port values only. |
| `audit`, `AB-audit`, `global` | Historical. Do not reuse. |
