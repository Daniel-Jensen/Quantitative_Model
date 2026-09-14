"""Estimate the Greek fiscal limit and back out def_scale_D / def_curvature_D.

Mirrors Bi-Foerster-Traum (FRBSF WP 2025-10) eq (3.1): they fit
    log(P/(1-P)) = eta0 + eta_s * s_{t-1} + u
on PRE-CRISIS Italian CDS, then let the fiscal-limit shifter carry the crisis.
Our shock_def_D is the shifter analogue, so the same split applies.

Data
  debt   Eurostat gov_10q_ggdebt, geo=EL, Maastricht GG consolidated gross debt,
         % of GDP, quarterly (the same concept BFT use).
  spread repo's Empirics/outputs/spreads_fred.csv, GRC 10y minus DEU 10y, monthly,
         averaged to quarterly.

Spread -> default probability uses the credit triangle at the MODEL's own
recovery_rate_D = 0.30, so the calibration is internally consistent:
    lambda_annual = spread / (1 - recovery);   P_quarterly = 1 - exp(-lambda/4)
These are RISK-NEUTRAL probabilities, as are BFT's CDS-implied ones — the
risk-premium wedge is inherited from their design, not introduced here.
"""
import json, urllib.request
import numpy as np
import pandas as pd

RECOVERY = 0.30                 # model recovery_rate_D
DEF_OFFSET = 0.05               # model def_offset_D, held fixed (2 moments, 2 params)
DR_SS = 1.152                   # model b_gov_ss_D, quarterly-GDP units (bank-held)
EBA_DATE = "2010-Q4"            # EBA moment base date, 2010-12-31

# ── data ─────────────────────────────────────────────────────────────────────
url = ("https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/gov_10q_ggdebt"
       "?format=JSON&lang=EN&geo=EL&unit=PC_GDP&sector=S13&na_item=GD")
with urllib.request.urlopen(url, timeout=60) as r:
    d = json.load(r)
idx = d["dimension"]["time"]["category"]["index"]
inv = {v: k for k, v in idx.items()}
debt = pd.Series({inv[int(i)]: v for i, v in d["value"].items()}).sort_index()
debt.index = pd.PeriodIndex(debt.index.str.replace("-Q", "Q"), freq="Q")
debt = debt / 100.0             # -> share of annual GDP

sp = pd.read_csv("/Users/Adam/Documents/uni/phd/research/QUANTITATIVE_MODEL/"
                 "Empirics/outputs/spreads_fred.csv", parse_dates=["date"])
sp = sp.set_index("date")["GRC_spread"].dropna() / 100.0     # pp -> decimal
spq = sp.resample("QE").mean()
spq.index = pd.PeriodIndex(spq.index, freq="Q")

df = pd.DataFrame({"spread": spq, "debt": debt}).dropna()
df["debt_lag"] = df["debt"].shift(1)
df = df.dropna()

lam = df["spread"] / (1.0 - RECOVERY)          # annual hazard
df["P_q"] = 1.0 - np.exp(-lam / 4.0)           # quarterly default probability
df = df[(df["P_q"] > 1e-9) & (df["P_q"] < 1 - 1e-9)]
df["logodds"] = np.log(df["P_q"] / (1 - df["P_q"]))

print(f"merged sample: {df.index.min()}..{df.index.max()}  n={len(df)}")
print(f"debt range {df['debt_lag'].min():.3f}..{df['debt_lag'].max():.3f} of annual GDP")
print(f"P_q  range {df['P_q'].min():.6f}..{df['P_q'].max():.4f}\n")


def ols(y, x):
    X = np.column_stack([np.ones_like(x), x])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    n, k = len(y), 2
    s2 = res @ res / (n - k)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    r2 = 1 - (res @ res) / ((y - y.mean()) @ (y - y.mean()))
    return b, se, r2, n


SAMPLES = {
    "pre-crisis 2000Q1-2009Q3 (BFT design)": (None, "2009Q3"),
    "pre-OMT 2000Q1-2012Q3 (preferred)": (None, "2012Q3"),
    "full 2000Q1-2026Q1": (None, None),
    "crisis 2009Q4-2013Q4": ("2009Q4", "2013Q4"),
}
fits = {}
for name, (lo, hi) in SAMPLES.items():
    s = df.copy()
    if lo: s = s[s.index >= pd.Period(lo, "Q")]
    if hi: s = s[s.index <= pd.Period(hi, "Q")]
    if len(s) < 6:
        print(f"{name}: n={len(s)}, skipped"); continue
    b, se, r2, n = ols(s["logodds"].values, s["debt_lag"].values)
    fits[name] = b
    print(f"{name}\n  eta0 = {b[0]:8.3f} ({se[0]:.3f})   eta_s = {b[1]:7.3f} ({se[1]:.3f})"
          f"   R2 = {r2:.3f}  n = {n}")
    print(f"  [BFT Italian pre-crisis: eta0 = -10.70, eta_s = 5.25]\n")


# ── map the fitted logistic onto the model's power function ──────────────────
# Model:  def_rate_D = def_scale * [ (dr + off)^curv - (dr_ss + off)^curv ],
#         dr = b_gov_D(-1)/Y_ss_D  in quarterly-GDP units of BANK-HELD debt.
# Scope:  bank-held = sigma * total, so dr = 4 * sigma * s.  sigma is pinned by
#         requiring the model SS debt to be the EBA base-date debt ratio.
s_eba = float(debt.loc[pd.Period(EBA_DATE.replace("-Q", "Q"), "Q")])
sigma = (DR_SS / 4.0) / s_eba
print(f"scope factor: model SS bank-held debt = {DR_SS/4:.4f} of annual GDP; "
      f"total debt at {EBA_DATE} = {s_eba:.4f}  ->  sigma = {sigma:.4f}")
print(f"  (model dr = 4*sigma*s, so dr_ss = {4*sigma*s_eba:.4f})\n")


def logistic_slope(b, s):
    z = b[0] + b[1] * s
    P = 1.0 / (1.0 + np.exp(-z))
    return b[1] * P * (1 - P)


for name, b in fits.items():
    # Two moments: the data slope dP/ds at two debt levels -> def_scale, curv.
    s1, s2 = 1.05, s_eba                       # pre-crisis norm and EBA base date
    m1, m2 = logistic_slope(b, s1), logistic_slope(b, s2)
    if m1 <= 0 or m2 <= 0:
        print(f"{name}: degenerate slopes, skipped"); continue
    # data slope in MODEL units: dP/d(dr) = (dP/ds) / (4*sigma)
    g1, g2 = m1 / (4 * sigma), m2 / (4 * sigma)
    dr1, dr2 = 4 * sigma * s1, 4 * sigma * s2
    # model slope: def_scale*curv*(dr+off)^(curv-1); ratio kills def_scale
    curv = 1.0 + np.log(g2 / g1) / np.log((dr2 + DEF_OFFSET) / (dr1 + DEF_OFFSET))
    scale = g2 / (curv * (dr2 + DEF_OFFSET) ** (curv - 1.0))
    print(f"{name}\n  implied def_scale_D = {scale:.4f}   def_curvature_D = {curv:.4f}"
          f"   [current 0.25 / 0.50]")
    # what the model's CURRENT calibration implies for the same slope
    cur = 0.25 * 0.5 * (dr2 + DEF_OFFSET) ** (0.5 - 1.0)
    print(f"  slope at EBA debt: data {g2:.5f} vs current calibration {cur:.5f} "
          f"({g2/cur:.2f}x)\n")
