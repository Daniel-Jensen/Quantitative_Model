# Global-projection model: current calibration inventory (branch `calibration`, working tree, 2026-10-23)

Scope: `global_projection/` + `calibration/global_projection.py` only, read from the UNCOMMITTED working
tree. Sources are repository files cited as `path:line`. Line numbers are those of the working tree on
2026-10-23. Nothing was run except `git diff`/`git log`, JSON printing, and a hand-written recomputation
of the closed-form SS formulas (marked "recomputed"). No model solve or test was run. Note that
`Claude files/diagnostics/calibration_walk/runs/x1.log` was being written by another session at 15:28
on 2026-10-23, so a calibration walk is in progress.

## Q1. What is in `calibration/global_projection.py`, key by key, and what changed versus `main`

### Takeaway
The working-tree file has about 90 keys. The Oct-2026 changes are these: leverage goes from 5 to 8. Six
stochastic-rest-point (SRP) targets are added (`srp_*`), with six deterministic-SS instruments set
equal to them. `sov_exposure_D/F` (14.4% / 6.5%) and `foreign_share_D/F` (20% / 0.5%) replace the
hard-coded `B_gov_*_ss = 0.98` and `b_*_ss = 0.196`. Three keys are new: `mu_eps = 1e-2`,
`pd_rest = 0.001` and `book_eps = 1e-4`. The six instruments still equal the targets, so the Broyden
SRP calibration result has NOT been pasted in. Measured rest-point moments at these values miss the
D-spread target badly (see Q6).

### Cited Findings

**Diff vs main** ([git diff main -- calibration/global_projection.py](calibration/global_projection.py)):
- Removed on main: `leverage_target_D/F=5.0` with the comment "leverage 5 (Bocola)". Also removed:
  `b_F_D_ss=b_D_F_ss=0.196` ("~20% of each stock held abroad"), `excess_return_*_ss=0.0`, and
  `B_gov_D_ss=B_gov_F_ss=0.98` ("debt sized so sovereigns are 7.6% of D-bank assets (Bocola)").
- Added: `srp_leverage_D/F=8.0`, `srp_spread_D/F=0.0025`, `srp_exposure_D=0.144`,
  `srp_exposure_F=0.065` (lines 40-42). `leverage_target_D/F=8.0`,
  `credit_spread_target_D/F=0.0025` and `sov_exposure_D/F=0.144/0.065` (lines 44-46).
  `foreign_share_D=0.20`, `foreign_share_F=0.005` (line 60). `book_eps=1e-4` (line 65).
  `mu_eps=1e-2` (line 82). `pd_rest=0.001` (line 85).
- Main computed the B and b stocks by hand. In the working tree they are SOLVED in the SS
  (`steady_state._size_books`, comment lines 61-62).
- `git diff --stat main -- calibration/ global_projection/ run.py`: 17 files, +1350/-450.
  `steady_state.py` +362/-..., `calibrate_stochastic.py` and `recursive_experiment.py` about 336 each.

**Parameter-by-parameter table** (values from `calibration/global_projection.py`. "How set" is as
documented in code or docs. D = Greece, F = Germany/rest of union.)

| Key (line) | Value | Meaning | How set / documented source | Greek data? |
|---|---|---|---|---|
| `sigma_D/F` (8) | 1.0 | CRRA, log utility | "log utility" comment, no source | n/a (standard) |
| `frisch_D/F` (9) | 2.0 | Frisch elasticity, GHH | "Frisch elasticity 2 (Bocola)" | Borrowed, Bocola (Italy) |
| `chi_D/F` (10) | 0.5417 (placeholder) | labour disutility | Overwritten in the SS: `chi = w_ss/P_CES` so N_ss = 1 (`steady_state.py:299-301`) | internal normalisation |
| `n_e` (13) | 2 | income states | Rouwenhorst, no source | none |
| `rho_e` (14-15) | 0.9 | income persistence | no source | **no data source** |
| `sigma_e` (14-15) | 0.2 | innovation sd of log e | no source | **no data source** |
| `a_min` (18-19) | 0.0 | borrowing limit (no borrowing) | no source | none |
| `a_max` (18-19) | 87.2 | asset-grid top | numerical, no source | none |
| `n_a`, `a_curve` | 250, 2.0 | grid size, curvature | numerical | n/a |
| `epsilon_D/F` (22) | 6.0 | goods demand elasticity (mc = 5/6, markup 20%) | no source | **no data source** |
| `Z_ss_D/F` (23) | 0.45 (placeholder) | TFP | Overwritten so Y_ss = 1 (`steady_state.py:88-92`) | normalisation |
| `alpha_D/F` (26) | 0.30 | capital share | "(Bocola)" | Borrowed, Bocola |
| `delta_D/F` (27) | 0.025 | depreciation per quarter | no source | **no data source** (standard 10%/yr) |
| `ksi_D/F` (28) | 0.50 | Jermann adjustment-cost curvature | "inside Bocola's posterior" | Borrowed, Bocola |
| `f_D/F` (31) | 0.08 | banker exit/payout share | "vs Bocola's 0.0354 (the bank block's one departure)". The reason (SS fold) is in `steady_state.py:193-195` | **no data source**; departure from Bocola |
| `r_dep_*_target` (33) | 0.003/qtr (1.2%/yr) | real deposit (risk-free) rate | "risk-free rate 1.003 a quarter (Bocola)" | Borrowed, Bocola. The table shows data 1.11% "Real deposit rate 1999-2007" (`calibrate_stochastic.py:202-204`) |
| `beta_inter_D/F` (35) | 0.997 | banker discount factor | "beta*R = 1 at the SS" (1/1.003 = 0.99701) | derived from r_dep |
| `srp_leverage_D/F` (40) | 8.0 | TARGET: assets/net worth at the SRP | "GR financial accounts 2002-07: 7.98; not Bocola's 5" (36-37) | Greek. **F uses the Greek number too** |
| `srp_spread_D/F` (41) | 0.0025/qtr = 100 bp/yr | TARGET: loan wedge r_wc - r_dep at the SRP | "not his 8, so the constraint binds at rest" (38). The table cites Gertler-Kiyotaki and 2003-07 lending spreads GR 290 / DE 180 bp as "upper bound" (`calibrate_stochastic.py:197-201`) | modelling choice, **no direct data match** |
| `srp_exposure_D` (42) | 0.144 | TARGET: own-sovereign / own-bank assets | "EBA 2010: GR 14.4%" | Greek EBA (end-2010 data from the EBA 2011 exercise) |
| `srp_exposure_F` (42) | 0.065 | same, F | "DE 6.5%" | German EBA |
| `leverage_target_*`, `credit_spread_target_*`, `sov_exposure_*` (44-46) | 8.0 / 0.0025 / 0.144, 0.065 | deterministic-SS INSTRUMENTS | Meant to be solved by Broyden (`calibrate_stochastic.py:28-35, 95-146`). Currently **equal to the targets**, i.e. not yet calibrated | — |
| `lambda_*`, `omega_ent_*` (48-51) | 0.22, 0.002 (warm starts) | divertable share (single λ), entrant transfer | Overwritten by `calibrate_bank_targets` (see Q2) | internal |
| `kappa_nfa` (54) | 0.01 | debt-elastic premium | "Bocola's debt-elastic premium on the cross-border position". In code it multiplies the cross-border DEPOSIT position over Y in deposit-UIP (`point_map.py:488-489`) | Borrowed, Bocola |
| `psi_bF_D`, `psi_bD_F` (56) | 2.0 | cross-border bond adjustment cost | "they pin each sovereign's split between the banks" | **no data source** |
| `foreign_share_D` (60) | 0.20 | share of the D stock held by F banks | EBA 2010 bracket: "German banks 12.7%, all foreign banks 39.6%" | Greek, but a value chosen inside a bracket |
| `foreign_share_F` (60) | 0.005 | share of the F stock held by D banks | "~0 in the data; positive so the book stays clear of the solver's floors" | numerical choice. EBA data gives 0.13% (Q3) |
| `book_eps` (65) | 1e-4 | Fischer-Burmeister smoothing of the no-short-sale pairs | numerical | n/a |
| `delta_b_D/F` (69) | 0.056 | HM perpetuity decay per quarter | "Bocola's pi: Macaulay duration 4.1 years at the 2.2% SS yield (Greek banks' Greek book, EBA 2010: 4.4 years)" | **Borrowed, Bocola.** The EBA ladder in the repo implies 0.0777 for the Greek book (Q3) |
| `recovery_rate_D` (72) | 0.45 | recovery on default (55% haircut) | "Greek PSI; Bocola" | Greek PSI. The table brackets face 46.5% and NPV 35-41% (`calibrate_stochastic.py:214-215`) |
| `sigma_s` (74) | 0.4455 | innovation sd of the logit-default factor s | "Bocola's 0.63 under his physicists' Gauss-Hermite rule, i.e. 0.63/sqrt(2)" | **Borrowed, Bocola (Italy)** |
| `rho_s` (not in cal) | 0.95 | persistence of s | HARD-CODED at `state_grid.py:407` "Bocola's estimate" | **Borrowed, Bocola (Italy)**; not a calibration key |
| `mu_eps` (82) | 1e-2 | smoothing of the KKT switch mu = (z+sqrt(z^2+eps^2))/2 | numerical, user decision 2026-10-20. "raises a binding multiplier ~6%" | n/a (affects economics) |
| `pd_rest` (85) | 0.001/qtr | priced default probability at rest (the AR(1) mean of s) | "the table's 0.10%". The table's data column says "Default probability, GR 1y CDS, 2008Q3: 0.092%", "Set directly" (`calibrate_stochastic.py:205-206`) | Greek CDS, rounded |
| `tpi_key_D` (87) | 0.071 | D share of Eurosystem P&L | "capital key, as in calibration/ssj.py", = BoG 2.0 / (2.0 + BuBa 26.1) (`calibration/ssj.py:435, 502`) | ECB capital key |
| `tpi_on`, `tpi_cap_bp`, `tpi_eps`, `tpi_gz` (89-93) | False, 200, 1e-4, True | TPI switch, spread cap (bp/yr), smoothing, rule form | policy setting. `run.py:52-53` (`RUN_TPI=True`, `TPI_CAP=200.0`) overrides at run time | policy choice |
| `zeta_wc_D/F` (96) | 1.0 | share of the wage bill pre-financed by bank loans | no source (Bocola's phi is cited in `bank.py:101-104` for the mechanism only) | **no data source** |
| `debt_root_D/F` (99) | 0.93 | strength of the Bohn tax rule, stated as the debt root dB'/dx | no source | **no data source** |
| `G_D/F` (100) | 0.0 | government consumption | no source | none |
| `size_D`, `size_F` (103) | 1.0, 8.0 | country masses | "F is eight times D's size" | **no documented data source.** DE/GR 2010 GDP in `data/eba_moments.json` = 11.7 |
| `omega_home_D` (106) | 0.85 | D home-goods weight (D imports 15% of its basket) | no source | **no data source** |
| `epsilon_trade` (106) | 0.5 | trade elasticity | no source. `steady_state.py:7-8` calls it "on the Marshall-Lerner boundary" | **no data source** |
| `omega_home_F` (118) | 0.98125 (derived) | F home-goods weight | `1-(1-0.85)*1/8` so trade balances at p = 1 | derived |
| `T`, `tol_*`, `n_jobs`, `use_numba` (109-115) | 300, 1e-12/1e-10, 0, True | solver settings | numerical | n/a |

**Hidden or hard-coded calibration outside the dict:**
- `rho_s = 0.95`, `rho_z = 0.9`, `z_star = Z_ss_D` at `state_grid.py:403-411`.
- `S_BOX_PD = 0.001`, `S_COVER_SD = 2.75` at `state_grid.py:265-266`. The s box is always centred on
  p^d 0.1%, "whatever the rest point's p^d".
- `NW_FLOOR = 0.15`, "Bocola's net-worth floor, as a share of n_ss" (`run.py:43`). It is passed as
  `cal["nw_floor_frac"]` (`calibrate_stochastic.py:72`).
- `RISK_SHOCK_PD = 0.02` (headline shock, `run.py:51`) and `TFP_SHOCK = 0.01` (`run.py:50`).
- Consumption floors `X_FLOOR = 0.05*C_ss` (`point_map.py:110`), `alpha_cap = 40` (`point_map.py:356`)
  and `_MU_CAP = 0.95` (`point_map.py:25`).
- Box bands `BOX_KW` (`recursive_experiment.py:33-34`). Ladders: `HAIRCUT_LADDER`, `WC_LADDER`,
  `FX_LADDER`, `MU_EPS_START=1e-2` (`recursive_experiment.py:174-190`).
- `WALK_FROM` = the six instruments at their target values (`recursive_experiment.py:196-198`), with
  `WALK_STEP` 0.1 leverage, 10 bp spread and 0.5 pp exposure per rung (200-202).

### Inferences
- The six `*_target`/`sov_exposure_*` instruments are byte-identical to the six `srp_*` targets and to
  `WALK_FROM`. The calibration loop's converged instruments have therefore not been written back.
  The file currently encodes "targets = instruments". The docs say this is wrong at the rest point:
  the D spread instrument must rise to about 135-140 bp.
- Leverage 8 is a Greek financial-accounts number, but it is applied to the F (German) bank too. The
  EBA alternatives for F are 6.94 (credit book/CT1) and 13.3 (total book/CT1) (Q3).

### Gaps
- No saved run at the current calibration exists. `results/GLOBAL/data/calibration.json` is dated
  2026-10-12 and holds the OLD main calibration (leverage 5, `B_gov 0.98`, `p_ss = 1.0000`). Solved
  values of lambda, omega_ent, B_gov_D/F_ss, beta_D/F, Z_ss, chi and p_ss at the new calibration are
  therefore not readable from a file. The worklog quotes p_ss = 1.151 (`calibration_worklog_2026-10.md`, 2026-10-20 entry).
- `results/CALIBRATION/broyden_state.json` (the Broyden checkpoint, `calibrate_stochastic.py:244-247`)
  does not exist. Either the loop has not run, or it ran under a different results root.

## Q2. How targets become parameters (SS solver, bank, government, firms, trade, SRP loop)

### Takeaway
The deterministic SS is solved in sequence:
1. λ and ω_ent per country from (leverage, spread) in closed form.
2. r_k = r + spread.
3. Z_ss so that Y = 1.
4. A 2x2 linear solve for the debt stocks that put own-sovereign/assets at the exposure instruments,
   inside a bracketed root for p from D's current account.
5. ω_ent is re-solved with the foreign-bond leg.
6. χ = w/P_CES.
7. Cross-border friction anchors.
8. β_D, β_F (households) from each country's deposit market.
9. γ_τ from the target debt root at the SS bond price.

On top of this, an outer Broyden loop (`calibrate_stochastic.py`) moves the six instruments so the
STOCHASTIC rest point hits the six `srp_*` targets.

### Cited Findings
- **λ and ω_ent** (`bank.py:60-93`). The franchise value α solves α = Ω(1+r)/(1-μ), with
  Ω = β_inter[f + (1-f)α] and the hard-max multiplier m0 = Ω·s·θ/α, smoothed by `_ic_multiplier`.
  The least root is taken (`_least_root`, `bank.py:11-23`). Then λ = α/(θ(1+δ)), with the smoothing
  slack δ = ε²/(4μ(1-μ)) (`bank.py:42-48, 84`). Then
  ω_ent = D_val/θ - (1-f)(s + (s_for - s)·share_for), with D_val = 1 - (1-f)(1+r) (`bank.py:86-87`).
  This runs twice: first without the foreign leg, then with each bank's foreign-bond share of assets
  (`steady_state.py:36-53, 280, 288`).
  - Recomputed at θ = 8, s = 0.0025, f = 0.08, r = 0.003, β = 0.997, ε = 1e-2: α ≈ 1.3238,
    Ω ≈ 1.2940, μ ≈ 0.0196, λ ≈ 0.1653. Without smoothing λ = 0.1654.
  - ω_ent before the foreign leg = 0.07724/8 - 0.92·0.0025 = 0.007355.
  - The old main run (θ = 5) saved λ = 0.23634 and ω_ent = 0.013148 (`results/GLOBAL/data/calibration.json`).
    0.013148 is reproduced by the formula.
- **Capital FOC** (`steady_state.py:281-282`): r_k = r_dep + credit_spread_target, because a single λ
  makes every asset class earn the same IC excess return. Every SS bond therefore yields r + s. At
  equal D/F instruments the deterministic-SS sovereign spread y_D - y_F is zero.
  - Recomputed: Q_ss = δ_b/(r+δ_b+s) = 0.9106, flow yield 2.2%/yr, Macaulay duration 4.09 years.
- **Z_ss** (`steady_state.py:88-92`): Z = ((r_k + δ)/(mc·α))^α, so Y_ss = 1 exactly.
- **Firm** (`firms.py:9-24`): N_ss = 1. Working capital: w_ss = mc(1-α)Y/(1+ζ·r_wc), with
  r_wc = r_dep + spread. mc = (ε-1)/ε.
- **Debt stocks** (`steady_state.py:95-110`): own/assets = e is linear in (B_D, B_F), so the two
  exposure instruments give one 2x2 solve at a trial p. Then b_D_F = φ_D·B_D and b_F_D = sz·φ_F·B_F.
  `_check_targets` asserts leverage, spread, exposure (1e-12), Y = 1 and the foreign shares
  (`steady_state.py:184-205`).
- **p_ss** (`steady_state.py:129-155`): brentq on D's current account over Y on [0.5, 2.0]. NFI uses
  coupons r + issuer spread (113-119). Consumption comes from P_CES·C = Y - I - G + NFI (122-126).
  The header says p is weakly identified: "trade balance moves ~0.012 of output per unit of log p"
  (lines 7-8).
- **Friction anchors** (`steady_state.py:168-181`): b_ref is set so the calibrated books are optimal
  when spreads differ. With equal spreads, b_ref = b.
- **Household β** (`steady_state.py:208-255`): brentq for β_c in [0.5, 1/(1+r) - 1e-4] so that
  household assets equal bank deposit supply (θ-1)n/P_CES, nation by nation. The SS cross-border
  deposit position is zero. Household income is y_e = (w/P_CES)·e + (Div - Tax)/P_CES, with
  Div = (1-mc)Y + bank dividends (`steady_state.py:227-240`).
- **Bohn rule** (`government.py:7-21`): Tax_ss = G + δ_b·B(1-Q). γ_τ = δ_b - (debt_root - (1-δ_b))·Q_B_ss,
  evaluated at the SS bond price. The comment says it was the risk-free price until 2026-10-19
  ("root 0.9294"). Recomputed at Q = 0.9106: γ_τ ≈ 0.0687.
- **Trade** (`trade.py`): a CES basket with elasticity `epsilon_trade` and per-country home bias.
  `omega_home_F` is derived for size consistency (`calibration/global_projection.py:118`).
- **Capital** (`capital.py:5-11`): Jermann γ0 = δ^ξ/(1-ξ) and γ1 = -δξ/(1-ξ), pinned by δ and ξ.
- **SRP loop** (`calibrate_stochastic.py:28-35, 51-56, 95-146`). It pairs 6 instruments with 6
  moments one to one. Tolerances: leverage 0.005, spread 0.1 bp, exposure 0.01 pp. Bounds: leverage
  [6, 11], spread [10, 150] bp/qtr-equivalent (0.0010-0.0150), exposure_D [8, 25]%, exposure_F
  [3, 12]%. The method is Broyden in tolerance-scaled units, at most 12 evaluations, each a full
  fine-grid solve. The moments are leverage = (deposits + n)/n, spread = r_wc - r_dep, and
  exposure = Q_bD·b_DD/assets, all read at the rest point.
- **Rest-point convention**: the rest point reached from the deterministic SS by a zero-shock
  simulation (worklog, 2026-10-21/22 "User decision").

### Inferences
- In the model the exposure numerator is a MARKET value (Q_bD·b_DD), but the EBA data numerator is
  a FACE value (gross direct long exposure). At end-2010 the Greek banks' own book had
  MV/face = 41,856/54,447 = 0.77 (`data/eba_moments.json` ladder.b_D_D). The target is therefore not
  measured on the same basis as the model moment.
- The model's "assets" are capital + working-capital loans + bonds (deposits + n). The EBA
  denominator is total assets (TA), which includes interbank, reserves and retail. The leverage
  target, by contrast, is not TA/CT1 (Q3).

### Gaps
- β_D/β_F (household discount factors) at the new calibration are not recorded in any saved file.

## Q3. `data/eba_moments.json`: contents, and what the global model uses or ignores

### Takeaway
The file holds the EBA 2011 stress-test base-date data (31-Dec-2010) for 6 Greek and 12 German
banks: raw EUR-million balances, maturity ladders, MTM losses, leverage alternatives, and two sets of
"model_targets" built for the OLD SSJ calibration. The global model reads nothing from this file at
run time. Its exposure targets (14.4% / 6.5%) and the 12.7% foreign-share bracket are reproduced
from `raw_EURm` (face/TA). Its leverage (8), duration (δ_b 0.056) and foreign-share values (20%,
0.5%) are NOT the json's numbers.

### Cited Findings (all from `data/eba_moments.json`)
- **meta**: "EBA 2011 EU-wide stress test", base date "2010-12-31 (actual; scenario blank)",
  sovereign gross direct long exposures (34010) in 7 maturity buckets.
  - Country map: D = Greece (GR030-035), F = Germany (12 banks DE017-DE029).
  - GDP from Eurostat nama_10_gdp 2010; yields from FRED 10y via `Empirics/outputs/spreads_fred.csv`.
  - `K_over_Y_annual` 2.7, `beta_inter` 0.9975155088, and coupon assumptions GGB 4.7% / Bund 3.5%.
  - `NOT_used`: adverse-scenario CT1 depletion.
- **raw_EURm**:
  - CT1: D 22,778; F 114,317.
  - RWA: D 222,466; F 1,222,402.
  - TA: D 377,200; F 4,872,189.
  - Sovereign books: b_D_D 54,447; b_D_F (DE banks' GGBs) 7,934; b_F_D (GR banks' Bunds) 410.7;
    b_F_F 315,313.
  - Corporate/CRE own and total; K_bank D 70,665 / F 470,272; GK_assets D 125,523 / F 793,519.
  - GDP_ann: D 223,590.5; F 2,615,260.
- **ladder**:
  - b_D_D: weighted residual maturity 5.13 y; yield 12.01% (Dec-2010) and 21.14% (Dec-2011); MV
    41,856; **Macaulay 3.49 y**; modified 3.12; 2011 MTM loss 21.6% of MV (EUR 9,054m);
    **delta_b_implied 0.0777**.
  - b_D_F: Macaulay 4.54 y; δ_b 0.0593.
  - b_F_F: Macaulay 4.34 y; δ_b 0.0568.
  - b_F_D: Macaulay 6.64 y; δ_b 0.0363.
- **mtm**:
  - dNW per 100 bp: D -0.0573, F -0.119.
  - 2011 MTM loss / CT1: D -0.397, F +0.123.
  - Fair-valued share of own book: D 11.1%, F 45.4%.
- **leverage_alternatives**: TA/CT1 D 16.56, F 42.62; total book/CT1 D 7.21, F 13.30.
- **model_targets** (narrow; SSJ-era):
  - n_inter D 0.407 / F 0.175; φ_bD_D 2.39; φ_bF_F 2.76; φ_bD_F 0.069; φ_bF_D 0.018.
  - θ_D 5.51, θ_F 6.94; ω_K D 0.117 / F 0.067; δ_b D 0.0777 / F 0.0568.
  - B_supply_D_qgdp 1.116, F 0.483.
- **model_targets_broad**: same θ and δ_b; n_inter D 2.14; ω_K = 1.
- **identification** block (SSJ-era; stale for the global model):
  - "f_D/F: banker exit rate; 0.12, from the GK literature".
  - "recovery_rate_D/F: 0.30 from the PSI NPV literature".
  - Δ divertability weights "Hardcoded 0.2/0.4"; `psi_lambda_B`, `def_scale` 0.25, `phi_lamb`.

**Use by the global model:**
- `srp_exposure_D = 0.144` = b_D_D/TA_D = 54,447/377,200 = 14.43% (recomputed). Used.
- `srp_exposure_F = 0.065` = 315,313/4,872,189 = 6.47%. Used.
- 12.7% bracket = b_D_F/(b_D_D + b_D_F) = 7,934/62,381 = 12.72%. Used as the lower end of the
  `foreign_share_D` bracket. The model value is 20%. The 39.6% "all foreign banks" figure is not in
  the json.
- `foreign_share_F` data = 410.7/(315,313 + 410.7) = 0.13%. The model uses 0.5%.
- "Leverage, GR credit book / CT1, 2010: 5.5" (`calibrate_stochastic.py:195`) = model_targets.theta_D
  5.51. Reported as a "narrower measure", not targeted.
- "GR sovereign / GR Core Tier 1: 2.39" (`calibrate_stochastic.py:218`) = φ_bD_D_ss. Labelled
  "Limitation".
- The repricing-loss row "40%" matches mtm_loss_2011_over_CT1 D = 0.397. The model revalues at
  9.1% (`calibrate_stochastic.py:40, 185`).
- **Ignored**: δ_b implied (0.0777 GR), Macaulay 3.49 y, TA/CT1 leverage (16.6), n_inter, ω_K,
  B_supply, the fair-valued share, F leverage alternatives, the beta_inter in meta, and the
  identification block.

### Inferences
- The code comment "Greek banks' Greek book, EBA 2010: 4.4 years" (`calibration/global_projection.py:68`)
  does not match the json's Greek-book Macaulay duration of 3.49 y (at the Dec-2010 12.01% yield).
  It is close to the German banks' Bund book (4.34 y). The table also has a row "Duration at 4.5%:
  4.4 years" (`calibrate_stochastic.py:226`), so 4.4 y may be the ladder's duration at a 4.5% yield.
  No file reproduces that number. The model's own duration at 4.5% is 3.76 y (recomputed:
  (1+y)/(y+δ_b)/4).
- "EBA 2010" in the code comments refers to end-2010 data from the EBA 2011 exercise.
- Leverage 8 (financial accounts) is about the json's total-book θ_D of 7.21. It is not the
  credit-book 5.51 or the TA-based 16.56, and the json does not contain it.

### Gaps
- No file or script in `Empirics/`, `data/` or `Claude files/` derives the "39.6%", "7.98",
  "0.092%", "25.3 bp", "84%", "108 bp", "290/180 bp" or "40% debt/GDP" data values. The search
  returned only coincidental numeric hits in `Empirics/outputs/ms_regime_*_summary.csv`. They exist
  only as hard-coded strings in `calibrate_stochastic.py:193-236` (and `reporting/export.py`).

## Q4. What the docs say about the new calibration, and doc/code inconsistencies

### Takeaway
The only up-to-date record of the Oct-2026 calibration is `Claude files/docs/calibration_worklog_2026-10.md`
plus the user's memory file. `STATE.md` Part II, `guide_global.md`, `SPEC.md`, `HANDOFF.md`, `PROGRESS.md`
and the paper's `03-calibration.tex` still describe the leverage-5 / 7.6% / p = 1 calibration, or
older ones still.

### Cited Findings
- **Worklog** (`calibration_worklog_2026-10.md`, 2026-10-20 entry):
  - Leverage 8 "in BOTH banks"; wedge 100 bp in both; exposures 14.4% / 6.5% "(EBA 2010)"; F banks
    hold 20% of the D stock, D banks 0.5% of the F stock.
  - Six targets at the SRP with tolerances 0.005 / 0.1 bp / 0.01 pp.
  - "Kept: r_dep 1.2%, p^d 0.1%/qtr at rest, recovery 45%, delta_b 0.056 (4.1-year duration at the
    2.2% yield), sigma_s 0.4455, rho_s 0.95". p_ss = 1.151.
  - γ_τ fixed at the SS bond price (the realised root was 0.9294).
- **Memory** (`~/.claude/.../memory/calibration-branch-2026-10.md`):
  - The repricing loss revalues at 9.1% p.a. = 2010-average Greek 10y (FRED 9.09%).
  - "D banks are INSOLVENT in default (63% direct bond loss + more); the 15% NW floor binds at
    ~22/25 nodes". The user chose to FUND the floor.
  - Final mu_eps = 1e-2 (3.7e-3 sat at a fold).
- **PROGRESS.md**: the newest entry is 2026-10-06 (household budget, `hh-budget-fix`). There is no
  entry for the Oct-19+ calibration. The worklog says it is the "Working draft of the PROGRESS.md entry".
- **HANDOFF.md** working-tree diff: only the branch convention and history text moved from CLAUDE.md
  on 2026-10-22. No calibration content.
- **STATE.md Part II** (line 2470 onward): the newest note is 2026-10-06. There are no
  leverage-8 / exposure / SRP numbers.
- **guide_global.md** (stale):
  - Line 100: "leverage 5, exposure 7.6%, recovery 0.45 are his [Bocola's]".
  - Line 184: "per-capita SS is UNCHANGED (p_ss = 1, identical n_ss, leverage 5, μ_ss = 0.001...)".
  - Lines 198-202: "Symmetric steady state".
  - Line 211: "beta_inter ≈ β_hh ≈ 0.99" (code 0.997).
  - Lines 337-345: "UNIQUE GLOBAL ATTRACTOR" (now non-unique, see Q6).
  - Lines 362-375: mu_rest 0.0098 at 100 bp (leverage-5 era).
- **03-calibration.tex** (paper; header "CORRECTED against code/global/ (2026-08-22)"):
  - Lines 134, 191: "f_X=0.04", "θ̄=5 and f=0.04, both posterior means of Bocola". Code f = 0.08,
    and the code comment says Bocola's f is 0.0354.
  - Line 207: "ς̄=8 basis points" (code 100 bp).
  - Line 280: exposure "8.9% [tbc]", citing Bocola (code 14.4%).
  - Line 297: "γ_τ=1, the value in Bocola" (code: γ_τ solved from debt_root 0.93, ≈0.069).
- **SPEC.md** "Calibration strategy" (line 103 onward) points to STATE.md and eba_calibration.md
  (SSJ, 2026-07-22): δ_b 0.10, `psi_lambda_B` and so on. These are SSJ objects, not the global
  calibration.
- **Code-internal inconsistencies:**
  - (a) The `steady_state.py:228-231` comment says the working-capital financing income stays with
    the BANK, "NOT the household". `recursive_main.py:48-51` and PROGRESS 2026-10-06 say "The
    working-capital flow is now in the budget". The worklog records "welfare_quintiles adds wc_ss to
    SS dividends, which the SS household block does not" (recorded, not fixed).
  - (b) The `WALK_FROM` comment says d = 0 rungs stall "above a D SS-spread instrument of ~125 bp"
    (`recursive_experiment.py:192-195`). The worklog says the SS does not exist above ~240 bp and the
    SS cold start fails at 200 bp.
  - (c) The `S_REFINE` comment in run.py:46 says 27 points (axis mode). guide_global line 88 and
    the STATE note say 115 (tensor mode, now unused).
  - (d) Banker exit f is cited three ways: code 0.08 vs "Bocola's 0.0354"; paper 0.04 "posterior
    mean of Bocola"; eba_moments identification 0.12 "GK literature".
  - (e) Recovery is cited three ways: code 0.45 "Greek PSI; Bocola"; guide "Bocola's D = 0.55"
    (consistent, as a haircut); eba_moments identification 0.30 "PSI NPV literature". The calibration
    table brackets face 46.5% and NPV 35-41%, so 0.45 lies above the NPV range.
  - (f) pd_rest is cited as "the table's 0.10%", but the table data column is 0.092% (GR 1y CDS
    2008Q3). The unit (per quarter or per year) of the 0.092% is not stated in the code.
  - (g) The real deposit rate is 1.2%/yr in the model vs 1.11% "data" in the table row labelled
    "Set (Bocola)" (`calibrate_stochastic.py:202-204`).

### Inferences
- Per CLAUDE.md "Committing", a calibration change requires STATE.md, PROGRESS.md and HANDOFF.md
  updates. None has been made yet for this calibration.
- Values borrowed from Bocola (estimated on Italy) rather than from Greek data:
  - frisch 2;
  - alpha 0.30;
  - ksi 0.50;
  - r_dep 0.003;
  - kappa_nfa 0.01;
  - delta_b 0.056 (Bocola's pi);
  - sigma_s 0.4455 (0.63/sqrt 2);
  - rho_s 0.95;
  - the net-worth floor 0.15;
  - the working-capital mechanism.
- Recovery 0.45 is attributed both to the Greek PSI and to Bocola. The F bank's leverage 8 is a
  Greek number applied to Germany.

### Gaps
- I did not find the referenced "user's calibration table" source document (data column with sources)
  outside the hard-coded `TABLE_ROWS`.

## Q5. Heterogeneous-household block calibration

### Takeaway
Both countries get the same block: a 2-state Rouwenhorst income chain (rho 0.9, innovation sd 0.2),
no borrowing, a 250-point curved grid to 87.2, log utility, GHH with Frisch 2, and χ and β solved
in the SS. No Greek micro data disciplines any of these. The only household "validation" moment is
liquid deposits/GDP of 84% (`calibrate_stochastic.py:231`), and the docs do not report whether it is
matched.

### Cited Findings
- Income: `rouwenhorst.py:33-42`. The log-e half-range is psi = sigma·sqrt(n-1)/sqrt(1-rho²), so
  sigma_e is the INNOVATION sd. The grid is normalised to mean 1.
  - Recomputed: psi = 0.4588 (unconditional sd of log e). e = {0.571, 1.429}. Stay probability
    (1+rho)/2 = 0.95.
- Grid: `household.py:5-11`. a = a_min + (a_max-a_min)·linspace(0,1,250)^2, with a_min 0 and a_max 87.2.
- EGM with the GHH composite c - v(N). v(N_ss) = χ/(1+1/frisch) (`steady_state.py:238`). The
  borrowing constraint is a_min (`household.py:14-38`).
- χ_c = w_ss/P_CES at N = 1 (`steady_state.py:299-301`). β_c is solved by deposit-market clearing
  (`steady_state.py:243-255`). Households do not price default (guide_global.md:335-338).
- `calibrate_household_anchors` sets `hh_T_D/F` so the period-map budget reproduces C_ss. It is now
  a ~0 residual (`recursive_main.py:47-65`).

### Inferences
- With N fixed at 1 and GHH, the HA block affects the deposit-supply/β mapping and consumption
  responses. Its income-risk parameters have no documented data source (Greek or otherwise).

### Gaps
- The solved β_D and β_F, the MPC, and wealth-distribution moments at the new calibration are not in
  any saved file. Whether deposits/GDP = 84% is matched is not recorded.

## Q6. Known open calibration issues

### Takeaway
At the current (target-valued) instruments the six SRP targets are NOT hit: the D bank is
essentially slack at rest (D wedge about 2 bp vs 100 bp). The fix is to raise the D spread
instrument to about 135-140 bp. That walk is fragile: the rest point is non-unique and
non-monotone in the instrument. D banks are insolvent in default (net worth at the 15% floor), and
the deterministic SS does not exist above about 240 bp.

### Cited Findings
- **Rest-point moments at the current instruments** (`Claude files/diagnostics/calibration_walk/runs/floor_sweep.json`,
  floor 0.15, ok = true, worst 1.8e-8):
  - lev_D 7.62 and lev_F 7.98 (target 8).
  - D wedge 2.6 fitted / 2.0 exact bp, and F wedge 95.3 bp (target 100).
  - exp_D 14.31% and exp_F 6.50% (targets 14.4 / 6.5).
  - mu_D 0.0005 / 0.0004.
  - Sovereign spread at rest 184 / 181 bp. The table validation value is 25.3 bp
    (`calibrate_stochastic.py:212-213`).
  - n_D in default 0.150 (on the floor); α_D in default 7.48.
  - Rest point differs across starts by 2.07 half-widths.
  - With floor 0.60 the solve fails (worst 3.1e-5).
- **D slack at rest**: "At the target-valued instruments the D bank is slack at the stochastic rest
  point (mu_D ~1e-4, D loan spread ~2 bp): it earns the sovereign risk premium (net worth 5-9% above
  the SS)" (worklog, 2026-10-21/22). See also `calibrate_stochastic.py:30-32`.
- **SS non-existence**: "The deterministic SS does not exist above a D SS-spread instrument of ~240 bp
  at f_D = 0.08 (fold-blocked)... SS cold start fails at 200 bp" (worklog). The assertion is at
  `steady_state.py:191-195`.
- **Walk**: rest-point D wedge, fitted/exact, by D instrument:
  - 100 bp: 2.6/2.0
  - 110 bp: 7.5/4.2
  - 120 bp: 29.7/9.3
  - 130 bp: 77.5/31.3
  - 140 bp: 101.7/163.7 (lev_D 7.69, exp_D 14.04%, F 97.3 bp)
  - 150 bp: 74.9/59.9

  The worklog calls this "steep and non-monotone".
- **Non-unique rest point**: at 150 bp the D spread at rest is 75, 159, 53, 4 or 0.1 bp depending
  on the start, and D leverage ranges 7.7-11.5. The slowest eigenvalue is 0.987-0.994. kappa_nfa 0.05
  does not fix it. User decision: "binding at rest", with the rest point taken from the deterministic
  SS and the D spread instrument found by a 1-D search near 135-140 bp (worklog).
- **D insolvent in default**: "63% direct bond loss"; the 15% NW floor binds at ~22/25 nodes; the
  user chose to FUND the floor (memory file). Post-default states lie 3-8 half-widths outside the box
  and are read at the V5 box edge (`TWO_BOXES = False`, `recursive_experiment.py:36-41`).
- **pd_rest sensitivity** (worklog 2026-10-23):
  - The s box is centred on p^d 0.1% whatever `pd_rest` is (`state_grid.py:263-266, 379`).
  - At 0.02% the coarse grid misprices the D bond (0.980 vs 0.909 fundamental), and the axis-refined
    grid is about 0.03 below the fixed point.
  - The coarse "binding at 0.02%" came from the mispriced bond.
- **Smoothing**: mu_eps = 1e-2 raises a binding multiplier about 6%. Smaller values sit at a fold
  (`calibration/global_projection.py:76-82`).
- **Validation misses noted in the table design**:
  - Bank-held GR debt/GDP (40%) is "Not matched".
  - Euribor-Eurepo 8 bp is "Not matched".
  - Sov/CT1 2.39 is a "Limitation".
  - Lending spreads 290/180 bp and the overdraft spread rise of 108 bp are "Upper bound"
    (`calibrate_stochastic.py:193-236`).

### Inferences
- Even after the D instrument is raised, the deterministic-SS sovereign spread is zero, because both
  bonds yield r + their bank's spread and a single λ applies. Any GR-DE spread at rest comes from
  priced risk alone. It currently overshoots the 25.3 bp validation value about 7x (181 bp at
  p^d 0.1%).
- The F wedge at rest (95 bp) is close to its target, but D's is not. The six-instrument Broyden is
  therefore effectively a 1-D problem in `credit_spread_target_D`, in a region where the SS is near
  its existence limit (~240 bp) and the moment is non-monotone.

### Gaps
- The final instrument values from the in-progress walk/search (`x1.log`, still being written on
  2026-10-23) are unknown.
- No converged six-target result exists in the repository.
