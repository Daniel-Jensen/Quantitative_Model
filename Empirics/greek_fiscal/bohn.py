# THE GREEK FISCAL REACTION FUNCTION (BOHN 1998) AND THE DEBT ROOT IT IMPLIES IN THE MODEL'S TAX RULE.
# pb_t = a + b*debt_{t-1} + c*gap_t, annual, primary balance and debt in % of GDP, gap from HP(100) on
# log real GDP. The model's rule T = Tbar + gamma*(B - Bbar) on quarterly output has an annual Bohn
# coefficient 4*gamma, and debt root rho_B = (1 - delta_b) + (delta_b - gamma)/qbar, where qbar
# prices the bond at the deposit rate plus the intermediation wedge (the paper's deterministic SS).
# Writes Empirics/outputs/bohn_estimates.json (the calibration appendix's fiscal table).
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.filters.hp_filter import hpfilter

ROOT = Path(__file__).resolve().parents[1]
DELTA_B, R_BAR, WEDGE = 0.056, 0.003, 0.000125      # D coupon decay, deposit rate, wedge (quarterly)
QBAR = DELTA_B / (R_BAR + WEDGE + DELTA_B)          # deterministic-SS bond price, 0.947
BOHN_MODEL = 0.04                                    # the paper's annual Bohn coefficient (Mendoza-Ostry)
SAMPLES = [(1995, 2008), (1995, 2019), (2001, 2019), (2010, 2019)]


def series(raw, code):
    v = raw["series"][code]["values"]
    return pd.Series({int(k): float(x) for k, x in v.items()})


def root_from_bohn(b):
    # quarterly tax response gamma = b/4; the root of B' in (B h) at the steady-state price
    return (1.0 - DELTA_B) + (DELTA_B - b / 4.0) / QBAR


def main():
    raw = json.loads((ROOT / "greek_fiscal/data/imf_greece.json").read_text())
    pb, debt, g = series(raw, "GGXONLB_G01_GDP_PT"), series(raw, "GGXWDG_NGDP"), series(raw, "NGDP_RPCH")
    y = np.log((1.0 + g.loc[1980:2025] / 100.0).cumprod())
    cyc, _ = hpfilter(y, lamb=100)
    df = pd.DataFrame({"pb": pb, "debt_lag": debt.shift(1), "gap": 100.0 * cyc}).dropna()

    rho_model = root_from_bohn(BOHN_MODEL)
    out = dict(model=dict(rho_B=rho_model, gamma_quarterly=BOHN_MODEL / 4.0, bohn_annual=BOHN_MODEL, qbar=QBAR,
                          half_life_quarters=float(np.log(0.5) / np.log(rho_model))),
               stationarity_threshold_bohn=4.0 * DELTA_B * (1.0 - QBAR), estimates={})
    for a, b in SAMPLES:
        s = df.loc[a:b]
        r = sm.OLS(s["pb"], sm.add_constant(s[["debt_lag", "gap"]])).fit(cov_type="HAC", cov_kwds={"maxlags": 2})
        bohn = float(r.params["debt_lag"])
        rho = root_from_bohn(bohn)
        out["estimates"][f"{a}-{b}"] = dict(T=int(len(s)), bohn=bohn, se_hac=float(r.bse["debt_lag"]),
                                            gap_coef=float(r.params["gap"]), implied_root=rho,
                                            half_life_quarters=float(np.log(0.5) / np.log(rho)) if rho < 1 else None)
    path = ROOT / "outputs/bohn_estimates.json"
    path.write_text(json.dumps(out, indent=2))

    m = out["model"]
    print(f"model: root {m['rho_B']:.4f}, gamma {m['gamma_quarterly']:.4f}/qtr, Bohn {m['bohn_annual']:.3f}, "
          f"half-life {m['half_life_quarters']:.1f} qtrs; stationary iff Bohn > {out['stationarity_threshold_bohn']:.3f}")
    for k, e in out["estimates"].items():
        hl = f"{e['half_life_quarters']:.0f} qtrs" if e["half_life_quarters"] else "explosive"
        print(f"{k}: T={e['T']} Bohn {e['bohn']:.3f} (HAC se {e['se_hac']:.3f})  root {e['implied_root']:.3f}  {hl}")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
