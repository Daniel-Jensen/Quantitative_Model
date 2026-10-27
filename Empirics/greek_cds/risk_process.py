# THE SOVEREIGN-RISK PROCESS ON GREEK AND ITALIAN DATA: AR(1) IN THE LOGIT OF THE QUARTERLY
# DEFAULT PROBABILITY, s_t = c + rho*s_{t-1} + sigma*e_t, AS IN THE MODEL'S pi^d = 1/(1+exp(-s)).
# Writes Empirics/outputs/risk_process_estimates.json (the calibration appendix's risk-process table).
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

ROOT = Path(__file__).resolve().parents[1]
LGD = 0.55            # loss given default (Bocola; the model's 1 - recovery)
RHO_MODEL, SIGMA_MODEL, PI_REST = 0.9, 0.41, 0.001   # the paper's process (appendix B)


def logit(p):
    return np.log(p / (1.0 - p))


def quarterly_pd(spread_bp):
    # quarterly default probability from an annualised spread; equals Bocola's x/(4 LGD) for small x
    return 1.0 - np.exp(-spread_bp / 1e4 / 4.0 / LGD)


def ar1(s, rho_fixed=RHO_MODEL):
    s = s.dropna()
    y, x = s.values[1:], s.values[:-1]
    T = len(y)
    X = np.column_stack([np.ones(T), x])
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ b
    df = T - 2
    sig = float(np.sqrt(e @ e / df))
    se_rho = float(np.sqrt(sig**2 * np.linalg.inv(X.T @ X)[1, 1]))
    return dict(first=str(s.index[0].to_period("Q")), last=str(s.index[-1].to_period("Q")), T=T,
                rho=float(b[1]), se_rho=se_rho,
                sigma=sig,
                sigma_ci90=[sig * float(np.sqrt(df / chi2.ppf(0.95, df))),
                            sig * float(np.sqrt(df / chi2.ppf(0.05, df)))],
                sigma_given_rho_model=float((y - rho_fixed * x).std(ddof=1)))


def main():
    cds = pd.read_csv(ROOT / "greek_cds/data/cds_panel.csv", parse_dates=["Date"]).set_index("Date")
    bond = pd.read_csv(ROOT / "outputs/spreads_fred.csv", index_col=0, parse_dates=True)["GRC_spread"] * 100.0

    s_cds_gr = logit(quarterly_pd(cds["GR_1Y_PAR_MID1"].resample("QE").mean()))
    s_cds_it = logit(quarterly_pd(cds["IT_1Y_PAR_MID1"].resample("QE").mean()))
    # euro entry (2001) on: earlier spreads carry drachma exchange-rate risk
    s_bond = logit(quarterly_pd(bond.resample("QE").mean().loc["2001":].clip(lower=1.0)))

    est = {
        "GR bond 10y, full": ar1(s_bond),
        "GR bond 10y, tranquil": ar1(s_bond.loc[:"2007"]),
        "GR bond 10y, post-PSI": ar1(s_bond.loc["2013":]),
        "GR bond 10y, CDS overlap": ar1(s_bond.loc["2008-07":"2011-06"]),
        "GR CDS 1y, to 2011Q2": ar1(s_cds_gr.loc[:"2011-06"]),
        "GR CDS 1y, to 2011Q4": ar1(s_cds_gr.loc[:"2011-12"]),
        "IT CDS 1y": ar1(s_cds_it),
    }
    # a ten-year spread averages the hazard over the bond's life: rescale its volatility by the
    # CDS/bond ratio over the quarters both are observed
    ratio = est["GR CDS 1y, to 2011Q2"]["sigma"] / est["GR bond 10y, CDS overlap"]["sigma"]
    rescaled = est["GR bond 10y, full"]["sigma"] * ratio

    # how often Greece sat above the crisis state, against the model's ergodic frequency
    p_bond = quarterly_pd(bond.resample("QE").mean().loc["2001":])
    sd_model = SIGMA_MODEL / np.sqrt(1.0 - RHO_MODEL**2)
    s_bar = logit(PI_REST)
    freq = {f"{100*x:g}%": dict(data=float((p_bond > x).mean()),
                                model=float(norm.sf((logit(x) - s_bar) / sd_model)))
            for x in (0.02, 0.05, 0.091)}

    # the 2008Q3-2011Q2 run-up as draws from the model's process (mean imposed at the rest point)
    x = s_cds_gr.loc[:"2011-06"].dropna().values
    v = x[1:] - RHO_MODEL * x[:-1] - (1.0 - RHO_MODEL) * s_bar
    runup = dict(innovations=len(v), positive=int((v > 0).sum()), cumulative=float(v.sum()),
                 z_at_model_sigma=float(v.sum() / (SIGMA_MODEL * np.sqrt(len(v)))))

    out = dict(lgd=LGD, estimates=est, cds_over_bond_sigma_ratio=ratio,
               rescaled_full_sample_sigma=rescaled, share_of_quarters_above=freq, greek_runup=runup)
    path = ROOT / "outputs/risk_process_estimates.json"
    path.write_text(json.dumps(out, indent=2))

    for k, r in est.items():
        lo, hi = r["sigma_ci90"]
        print(f"{k:<26} {r['first']}-{r['last']}  T={r['T']:>3}  rho {r['rho']:.3f} ({r['se_rho']:.3f})"
              f"  sigma {r['sigma']:.3f} [{lo:.2f}, {hi:.2f}]  sigma|rho={RHO_MODEL:g} {r['sigma_given_rho_model']:.3f}")
    print(f"CDS/bond sigma ratio {ratio:.2f}; rescaled full-sample sigma {rescaled:.3f}")
    print("share of quarters above:", {k: (round(v['data'], 3), round(v['model'], 4)) for k, v in freq.items()})
    print("run-up:", runup)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
