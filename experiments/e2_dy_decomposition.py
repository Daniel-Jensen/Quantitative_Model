"""E2 — decompose the output response into its goods-market components.

Identity, from market_clearing_D (code/equations_D.py:139):

    Y_D = P_CES_D*C_D + I_D + G_D + Phi_D + T_D + NX_D

linearised, with the product rule on the consumption term:

    dY = P_ss*dC + C_ss*dP_CES + dI + dG + dPhi + dT + dNX

goods_mkt_D is a TARGETED residual held to <=1e-14, so this closes to solver
tolerance. That makes the decomposition self-verifying: a non-closing residual
means a missing term, not a small error to tolerate.

Gates two flagged claims:
  * SPEC "Where Germany genuinely benefits" — do not assert the trade channel
    without checking the investment/NX split.
  * STATE.md watch item — Y_D[0] is POSITIVE under both intervening regimes and
    the A5 dY_D trough never goes negative. This shows which component flips.
"""
import numpy as np

from common import load_cache, provenance, regime_irfs, write_results

CLOSURE_TOL = 1e-7

# Display order = economic reading order, not dict order.
COMPONENTS = ["consumption_quantity", "consumption_price", "investment",
              "government", "portfolio_cost", "macropru_tax", "net_exports"]


def decompose_dY(irf, ss):
    """Return ({component: path}, residual). Residual must be ~0 by the identity."""
    T = len(np.asarray(irf["Y_D"]))
    zero = np.zeros(T)

    def get(name):
        # An absent term is structurally zero (G_D is constant and absent from the
        # Jacobian; T_D is identically zero at T0=T1=0). Carried explicitly so the
        # identity stays complete if either is ever switched on.
        return np.asarray(irf[name]) if name in irf else zero

    components = {
        "consumption_quantity": float(ss["P_CES_D_ss"]) * get("C_D"),
        "consumption_price": float(ss["C_D_ss"]) * get("P_CES_D"),
        "investment": get("I_D"),
        "government": get("G_D"),
        "portfolio_cost": get("Phi_D"),
        "macropru_tax": get("T_D"),
        "net_exports": get("NX_D"),
    }
    residual = np.asarray(irf["Y_D"]) - sum(components.values())
    return components, residual


def run(horizon=40):
    cache = load_cache()
    ss = {"P_CES_D_ss": float(cache["P_CES_D_ss"]), "C_D_ss": float(cache["C_D_ss"])}
    Y_ss = float(cache["Y_D_ss"])
    beta = float(cache["beta_D"])
    disc = beta ** np.arange(horizon)

    payload = {"provenance": provenance(), "horizon": horizon,
               "closure_tol": CLOSURE_TOL, "regimes": {}}

    for name, (gamma, irf) in regime_irfs(cache).items():
        components, residual = decompose_dY(irf, ss)
        max_resid = float(np.max(np.abs(residual[:horizon])))
        assert max_resid < CLOSURE_TOL, (
            f"E2 identity does not close for regime '{name}': max|residual| = "
            f"{max_resid:.3e} > {CLOSURE_TOL:.0e}. A term is missing from the "
            f"decomposition — do not report these numbers.")

        dY = np.asarray(irf["Y_D"])[:horizon]
        payload["regimes"][name] = {
            "gamma": gamma,
            "max_abs_residual": max_resid,
            "dY_impact_pct_ss": float(dY[0] * 100.0 / Y_ss),
            "dY_trough_pct_ss": float(dY.min() * 100.0 / Y_ss),
            "dY_pv": float((disc * dY).sum()),
            "components_pv": {k: float((disc * components[k][:horizon]).sum())
                              for k in COMPONENTS},
            "components_impact": {k: float(components[k][0]) for k in COMPONENTS},
            "paths": {k: components[k][:horizon].tolist() for k in COMPONENTS},
            "dY_path": dY.tolist(),
        }

    write_results("e2_dy_decomposition", payload)
    return payload


if __name__ == "__main__":
    res = run()
    print(f"{'regime':>12} {'gamma':>9} {'dY[0] %SS':>11} {'dI PV':>12} "
          f"{'dNX PV':>12} {'dC_q PV':>12} {'max|resid|':>11}")
    print("-" * 84)
    for name, r in res["regimes"].items():
        c = r["components_pv"]
        print(f"{name:>12} {r['gamma']:>9.4f} {r['dY_impact_pct_ss']:>+11.4f} "
              f"{c['investment']:>+12.3e} {c['net_exports']:>+12.3e} "
              f"{c['consumption_quantity']:>+12.3e} {r['max_abs_residual']:>11.2e}")
    print("-" * 84)
    print("Identity: dY = P_ss*dC + C_ss*dP + dI + dG + dPhi + dT + dNX "
          "(market_clearing_D, targeted to <=1e-14)")
