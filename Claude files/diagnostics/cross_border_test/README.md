# Cross-border bond adjustment: theory prediction against data (2026-10-30)

**The term.** `psi_bD_F` (= 2) enters the F bank's D-bond FOC as `adj_F = 1 + psi*(b_DF - b_ref)/B_D`, so
`(b_DF - b_ref)/B_D = (F banks' perceived excess return on D bonds)/psi`. The excess return is the gap between the F bank's
valuation and the price, which the D bank sets.

**Prediction.** The gap is positive in a Greek crisis: D banks value D bonds less, because they are hit in exactly the
default states (risk premium) and become constrained (liquidity premium). F banks therefore BUY as risk rises; psi only
scales the inflow.

| model (fine grid, exact reads; `model_side.py`) | rest | p^d 2% | 5% | 12% |
|---|---|---|---|---|
| F banks' D book vs rest | 0 | +2.3% | +7.1% | +18% |
| F excess return (bp/yr) | −12 | +219 | +705 | +1812 |

**Data** (`compare.py`; BIS CBS claims on the Greek official sector, ultimate risk, EUR; ECB BSI Greek MFIs' euro-area
government securities; FRED 10y spread):

| | 2010Q1 | 2010Q4 | 2011Q4 |
|---|---|---|---|
| GR 10y spread (bp) | 314 | 910 | 1921 |
| German banks (EUR bn) | 17.1 | 11.0 | 5.2 |
| all foreign banks (EUR bn) | 68.6 | 34.6 | 17.6 |
| Greek banks (EUR bn) | 40.7 | 46.2 | 46.0 |
| German share of German + Greek holdings | 29.6% | 19.3% | 10.2% |

**Model over the same spread rise** (314 → 910 bp, p^d 2.6% → 10.9%):

| | model | data |
|---|---|---|
| German holdings, value | −7% | −36% |
| German holdings, quantity | +13.5% | about −21% (the value deflated by the model's −18% price change) |
| Greek banks' holdings, value | −18% | +14% |
| German share | 11.6% → 13.0% | 29.6% → 19.3% |

**Verdict.** The sign is wrong: the model has foreign banks buying, the data show them selling and Greek banks absorbing (home
bias rising). A finite psi only scales the wrong-signed inflow and cannot flip it.

**Caveats.**
- The ECB's SMP bought Greek bonds from May 2010, largely from foreign banks, which explains part of the 2010 fall; German
  banks kept selling through 2011.
- The BIS official sector is broader than government bonds.
- The 2011 spreads lie beyond the model's box.
