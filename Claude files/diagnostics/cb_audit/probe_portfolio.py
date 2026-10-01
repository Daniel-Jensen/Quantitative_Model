"""CB-block audit — the decisive portfolio diagnostic.

Reports the FULL 2x2 sovereign-holdings matrix (plus the CB row) before and after
TPI, in consistent aggregate market-value units.

WHY THIS NEEDS ITS OWN SOLVE. probe_pipeline.py dumped b_D_D and b_D_F only. The
matrix needs all four legs plus both government stocks, and the per-capita
conventions differ across them (see UNITS below), so partial data cannot be
rescaled after the fact.

UNITS -- the whole point of this diagnostic, and the thing that has bitten this
repo before (CLAUDE.md's "percentages must divide by their own SS level"):

  b_D_D   D banks holding D paper      D aggregate      (per D capita, D size = 1)
  b_F_D   D banks holding F paper      per D capita
  b_D_F   F banks holding D paper      PER F CAPITA  -> x size_F for an aggregate
  b_F_F   F banks holding own paper    PER F CAPITA  -> x size_F for an aggregate
  b_gov_D D government stock           D aggregate
  b_gov_F F government stock           PER F CAPITA  -> x size_F for an aggregate
  cb_buy_D CB book                     D aggregate (already; takes no weight)

  q_b_D, q_b_F are BOTH D-good prices (intermediation_P3_D adds q_b_F*b_F_D
  straight onto Q_D*K_D; external_account_D books the F-bond receipt with no p),
  so market values in D goods need no terms-of-trade conversion.

The two clearing identities the matrix must satisfy:
  D paper:  b_gov_D       = b_D_D + size_F*b_D_F + cb_buy_D
  F paper:  size_F*b_gov_F = size_F*b_F_F + b_F_D

Audit-only. Writes portfolio_matrix.md / .json / .npz. No model source touched.
"""
import os, sys, json, copy, datetime
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "code"))

LEGS = ["b_D_D", "b_F_D", "b_D_F", "b_F_F", "b_gov_D", "b_gov_F",
        "q_b_D", "q_b_F", "n_inter_D", "n_inter_F", "cb_buy_D",
        "spread_rb", "Y_D", "K_D", "theta_D", "theta_F"]
GAMMAS = [0, 2, 5, 10]


def ts():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(m):
    print(f"[{ts()}] {m}", flush=True)


def main():
    from calibration import get_calibration
    from steady_state import solve_steady_state
    from ic_delta_calibration import calibrate_ic_delta
    from depreciation_calibration import calibrate_depreciation
    from full_model import build_and_solve
    import tpi as tpi_mod

    cal = get_calibration()
    log("solving steady state...")
    ssr = calibrate_depreciation(calibrate_ic_delta(solve_steady_state(cal)))
    ss = ssr["ss_final"] if "ss_final" in ssr else ssr
    log("steady state solved; building dynamic model...")
    mr = build_and_solve(ssr)
    log("running TPI...")
    tr = tpi_mod.run_tpi(mr)
    log("TPI done.")

    G = tr["G_tpi"]
    T = tr["T"]
    sF = float(ss["size_F"])

    def g(k):
        try:
            return float(ss[k])
        except Exception:
            return float(ss.toplevel[k])

    SS = {k: (0.0 if k == "cb_buy_D" else g(k)) for k in LEGS}
    irfs = {gam: {k: np.asarray(tr["irfs_tpi"][gam][k]) if k in tr["irfs_tpi"][gam]
                  else np.zeros(T) for k in LEGS} for gam in GAMMAS}
    for gam in GAMMAS:
        irfs[gam]["cb_buy_D"] = np.asarray(tr["irfs_tpi"][gam]["cb_buy_D"])

    # ── the matrix, in AGGREGATE D-GOOD MARKET VALUE ─────────────────────────
    # weight[leg] converts the model variable to a D-capita aggregate quantity.
    W = {"b_D_D": 1.0, "b_F_D": 1.0, "b_D_F": sF, "b_F_F": sF,
         "b_gov_D": 1.0, "b_gov_F": sF, "cb_buy_D": 1.0}
    PRICE = {"b_D_D": "q_b_D", "b_F_D": "q_b_F", "b_D_F": "q_b_D",
             "b_F_F": "q_b_F", "b_gov_D": "q_b_D", "b_gov_F": "q_b_F",
             "cb_buy_D": "q_b_D"}

    def mv_level(leg, gam=None, t=None):
        """Market value in D goods. Level if gam is None, else level + deviation.

        Everything is a linearised deviation, so the market value moves as
        q*db + b*dq  (the dq*db cross term is second order and is dropped, in
        line with cb_pnl's convention in code/tpi.py)."""
        q, b = SS[PRICE[leg]], SS[leg]
        base = W[leg] * q * b
        if gam is None:
            return base
        dq, db = irfs[gam][PRICE[leg]][t], irfs[gam][leg][t]
        return base + W[leg] * (q * db + b * dq)

    def qty_level(leg, gam=None, t=None):
        base = W[leg] * SS[leg]
        if gam is None:
            return base
        return base + W[leg] * irfs[gam][leg][t]

    # ── clearing identities ──────────────────────────────────────────────────
    checks = {}
    checks["SS D paper"] = (qty_level("b_gov_D")
                            - qty_level("b_D_D") - qty_level("b_D_F")
                            - qty_level("cb_buy_D"))
    checks["SS F paper"] = (qty_level("b_gov_F")
                            - qty_level("b_F_F") - qty_level("b_F_D"))
    for gam in GAMMAS:
        for t in (0, 4, 20):
            checks[f"g{gam} t{t} D paper"] = (
                qty_level("b_gov_D", gam, t) - qty_level("b_D_D", gam, t)
                - qty_level("b_D_F", gam, t) - qty_level("cb_buy_D", gam, t))
            checks[f"g{gam} t{t} F paper"] = (
                qty_level("b_gov_F", gam, t) - qty_level("b_F_F", gam, t)
                - qty_level("b_F_D", gam, t))

    L = []
    def p(s=""):
        print(s); L.append(s)

    p("# The 2x2 sovereign-holdings matrix, before and after TPI\n")
    p(f"Generated {ts()} by `diagnostics/cb_audit/probe_portfolio.py`. "
      f"`size_F` = {sF:.6f}.\n")
    p("All entries are **aggregate market value in D goods**, `q_b * quantity`, with "
      "per-F-capita legs (`b_D_F`, `b_F_F`, `b_gov_F`) scaled by `size_F`. "
      "`q_b_D` and `q_b_F` are both D-good prices, so no terms-of-trade conversion "
      "enters.\n")

    p("## Clearing identities (must be ~0)\n")
    p("| point | D paper residual | F paper residual |")
    p("|---|---|---|")
    p(f"| steady state | {checks['SS D paper']:+.3e} | {checks['SS F paper']:+.3e} |")
    for gam in GAMMAS:
        for t in (0, 4, 20):
            p(f"| gamma={gam}, t={t} | {checks[f'g{gam} t{t} D paper']:+.3e} | "
              f"{checks[f'g{gam} t{t} F paper']:+.3e} |")
    p("")

    # ── the matrix itself ────────────────────────────────────────────────────
    ROWS = [("D banks", "b_D_D", "b_F_D"),
            ("F banks", "b_D_F", "b_F_F")]

    def render(title, gam, t):
        p(f"### {title}\n")
        p("| holder | D paper | F paper | total | D paper, % of D issue |")
        p("|---|---|---|---|---|")
        tot_D = qty_level("b_gov_D", gam, t) * (SS["q_b_D"] if False else 1.0)
        mvD_issue = mv_level("b_gov_D", gam, t)
        colD = colF = 0.0
        for name, lD, lF in ROWS:
            vD, vF = mv_level(lD, gam, t), mv_level(lF, gam, t)
            colD += vD; colF += vF
            p(f"| {name} | {vD:+.6f} | {vF:+.6f} | {vD+vF:+.6f} | "
              f"{100*vD/mvD_issue:6.2f}% |")
        vcb = mv_level("cb_buy_D", gam, t)
        colD += vcb
        p(f"| **CB** | {vcb:+.6f} | 0.000000 | {vcb:+.6f} | "
          f"{100*vcb/mvD_issue:6.2f}% |")
        p(f"| **total held** | {colD:+.6f} | {colF:+.6f} | {colD+colF:+.6f} | |")
        p(f"| **issued** | {mvD_issue:+.6f} | {mv_level('b_gov_F', gam, t):+.6f} | "
          f"{mvD_issue+mv_level('b_gov_F', gam, t):+.6f} | |")
        p("")

    render("Steady state (TPI dormant, `cb_buy_ss = 0`)", None, None)
    for gam in GAMMAS:
        for t in (0,):
            render(f"Impact of the 1pp default shock, t={t}, gamma={gam}", gam, t)

    # ── deltas vs the no-TPI counterfactual, the decisive table ──────────────
    p("## What TPI moves: holdings at t=0 relative to gamma=0 (same shock)\n")
    p("Aggregate market value in D goods, and as a % of the D-paper stock at SS.\n")
    mvD_ss = mv_level("b_gov_D")
    p("| leg | SS level | g=0 | g=2 | g=5 | g=10 | d(g=10 - g=0) | as % of SS D issue |")
    p("|---|---|---|---|---|---|---|---|")
    for label, leg in (("D banks / D paper", "b_D_D"),
                       ("F banks / D paper", "b_D_F"),
                       ("CB / D paper", "cb_buy_D"),
                       ("D banks / F paper", "b_F_D"),
                       ("F banks / F paper", "b_F_F"),
                       ("D govt issue", "b_gov_D"),
                       ("F govt issue", "b_gov_F")):
        v = {gam: mv_level(leg, gam, 0) for gam in GAMMAS}
        d = v[10] - v[0]
        p(f"| {label} | {mv_level(leg):+.6f} | {v[0]:+.6f} | {v[2]:+.6f} | "
          f"{v[5]:+.6f} | {v[10]:+.6f} | {d:+.6f} | {100*d/mvD_ss:+.3f}% |")
    p("")

    # ── concentration ratios: what the IC actually sees ──────────────────────
    p("## Sovereign concentration `phi = q_b*b / n_inter` — what the IC sees\n")
    p("This is the object `intermediation_IC_D/F` reads, so it is where the "
      "portfolio shift becomes a constraint effect.\n")
    p("| ratio | SS | g=0 t0 | g=2 t0 | g=5 t0 | g=10 t0 |")
    p("|---|---|---|---|---|---|")
    for label, leg, n in (("phi_bD_D (D bank, own paper)", "b_D_D", "n_inter_D"),
                          ("phi_bF_D (D bank, F paper)", "b_F_D", "n_inter_D"),
                          ("phi_bD_F (F bank, D paper)", "b_D_F", "n_inter_F"),
                          ("phi_bF_F (F bank, own paper)", "b_F_F", "n_inter_F")):
        # phi is per-capita on both numerator and denominator -> no size_F weight
        q, b, nn = SS[PRICE[leg]], SS[leg], SS[n]
        base = q * b / nn
        row = [f"{base:.6f}"]
        for gam in GAMMAS:
            dq, db, dn = (irfs[gam][PRICE[leg]][0], irfs[gam][leg][0],
                          irfs[gam][n][0])
            row.append(f"{base + (q*db + b*dq)/nn - q*b*dn/nn**2:.6f}")
        p(f"| {label} | " + " | ".join(row) + " |")
    p("")

    payload = {"size_F": sF, "SS": SS, "checks": checks,
               "matrix_mv": {("ss" if gam is None else f"g{gam}"):
                             {leg: mv_level(leg, gam, None if gam is None else 0)
                              for leg in W}
                             for gam in [None] + GAMMAS}}
    json.dump(payload, open(os.path.join(HERE, "portfolio_matrix.json"), "w"),
              indent=2, default=str)
    np.savez_compressed(os.path.join(HERE, "portfolio_matrix.npz"),
                        **{f"irf_g{gam}_{k}": irfs[gam][k]
                           for gam in GAMMAS for k in LEGS})
    open(os.path.join(HERE, "portfolio_matrix.md"), "w").write("\n".join(L) + "\n")
    log("WROTE portfolio_matrix.md / .json / .npz")


if __name__ == "__main__":
    main()
