# COMPARISON EXPORT: THE GLOBAL IRFS IN THE SCHEMA THAT compare.py READS.
import json
import os

import numpy as np

from global_projection.solver_recursive.state_grid import default_prob

# comparison name -> key in dynamic_irf / tfp_irf paths (already % deviations)
SERIES = {"Y_D": "Y", "C_D": "C", "I_D": "I", "n_D": "n", "Q_bD": "dQ_bD",
          "Y_F": "Y_F", "C_F": "C_F", "I_F": "I_F", "n_F": "n_F"}


def _pick(path):
    # SHARED SERIES ONLY, AS PLAIN LISTS.
    return {name: np.asarray(path[k]).tolist() for name, k in SERIES.items() if k in path}


def _pct(v, digits=1):
    # A PERCENTAGE CELL IN LaTeX MATH.
    return f"${v:.{digits}f}\\%$"


def export_calibration_tex(v, out_path):
    # THE PAPER'S CALIBRATION TABLE WITH THE MODEL COLUMN FILLED FROM THE REST POINT.
    small = (f"${v['crisis_share']:.3f}\\%$" if v["crisis_share"] >= 1e-3
             else f"${v['crisis_share']:.1e}\\%$")
    rows = [
        ("Banks", None),
        ("Leverage, GR fin.\\ accounts, 2002--07", "$7.98$", f"${v['lev_D']:.2f}$", "Target"),
        ("Leverage, GR credit book / CT1, 2010", "$5.5$", f"${v['lev_D']:.2f}$",
         "Narrower measure"),
        ("Intermediation wedge", "$100$bp", f"${v['loan_D']:.0f}$bp", "Target"),
        ("Lending spread, 2003--07, GR; DE", "$290$; $180$bp",
         f"${v['loan_D']:.0f}$; ${v['loan_F']:.0f}$bp", "Upper bound"),
        ("Euribor--Eurepo, 2002--2012", "$8$bp", f"${v['loan_D']:.0f}$bp", "Not matched"),
        ("Real deposit rate, 1999--2007", "$1.11\\%$", _pct(v["rdep_par"], 2), "Set (Bocola)"),
        ("Sovereign risk and recovery", None),
        ("Default probability, Greek 1y CDS, 2008Q3", "$0.092\\%$", _pct(v["pd_rest"], 2),
         "Set directly"),
        ("Default probability, 2010--2011 average", "$9.1\\%$", small, "Validation"),
        ("Bond-implied hazard, 1999--2007", "$0.23\\%$", _pct(v["pd_rest"], 2), "Upper bound"),
        ("Bond-implied hazard, basis-adjusted", "$0.02$--$0.05\\%$", _pct(v["pd_rest"], 2),
         "Lower bound"),
        ("Greek--German 10y spread, 2002--2007", "$25.3$bp", f"${v['sov_spread']:.1f}$bp",
         "Validation"),
        ("PSI recovery, face value; NPV", "$46.5\\%$; $35$--$41\\%$", _pct(v["recovery"], 0),
         "Bracket"),
        ("Sovereign exposure", None),
        ("GR sovereign / GR bank assets", "$14.4\\%$", _pct(v["exp_D"]), "Target"),
        ("DE sovereign / DE bank assets", "$6.5\\%$", _pct(v["exp_F"]), "Target"),
        ("GR sovereign / GR Core Tier 1", "$2.39$", f"${v['sov_ct1']:.2f}$", "Limitation"),
        ("Foreign share, German banks; all foreign", "$12.7\\%$; $39.6\\%$",
         _pct(v["foreign_D"], 0), "Bracket"),
        ("Duration of Greek banks' Greek book", "$4.4$ years", f"${v['dur_rest']:.1f}$ years",
         "Validation"),
        ("EBA-held Greek debt / GDP, 2010", "$40\\%$", _pct(v["debt_mkt"], 0), "Not matched"),
        ("End-2011 repricing loss / GR CT1", "$40\\%$", _pct(v["reprice"], 0), "Validation"),
        ("Households", None),
        ("Liquid deposits / GDP, Greece", "$84\\%$", _pct(v["deposits"], 0), "Validation"),
        ("Transmission", None),
        ("GR$-$DE rise in overdraft spread", "$108$bp", f"${v['overdraft_q0']:.0f}$bp",
         "Upper bound"),
    ]
    body, first = [], True
    for r in rows:
        if r[1] is None:
            body.append(("    \\addlinespace\n" if not first else "")
                        + f"    \\multicolumn{{4}}{{l}}{{\\emph{{{r[0]}}}}} \\\\")
            first = False
        else:
            body.append(f"    {r[0]} & {r[1]} & {r[2]} & {r[3]} \\\\")
    note = (
        "Note: rates and spreads are annualised; default probabilities are\n"
        "  quarterly. Bank stocks are gross direct long positions at 31 December 2010\n"
        "  \\citep{eba2011}. Bond-implied hazards assume a flat hazard at $\\varrho_D=0.45$;\n"
        "  the basis adjustment deflates the ten-year spread by the ten-year/one-year CDS\n"
        "  ratio. CT1 is Core Tier 1 capital. The German foreign share is of Greek and\n"
        "  German banks' combined holdings; the all-foreign share is of the EBA total.\n"
        "  Durations are Macaulay, discounted at $4.5\\%$ in the data and at the model's\n"
        f"  rest-point yield (${v['y_D'] / 100:.2f}\\%$; ${v['dur_45']:.1f}$ years at $4.5\\%$). "
        "The repricing loss revalues\n"
        "  the 2010 book at a ten-year yield of $9.1\\%$, the 2010 average. Model moments are\n"
        "  read at the stochastic rest point, where leverage, the intermediation wedge\n"
        "  $r^{wc}-r$ and both exposure shares are targeted; model debt is held by banks,\n"
        "  at market value. Against the $2010$--$2011$ average, the model entry is the\n"
        "  ergodic share of quarters with $\\pi^d_t$ above it; the overdraft-spread rise is\n"
        "  the model's $D-F$ lending-spread rise on impact when the default probability\n"
        "  jumps to $2\\%$ a quarter.")
    tex = ("\\begin{table}[H]\n  \\centering\n"
           "  \\caption{Empirical targets and validation benchmarks}\n"
           "  \\label{tab:targets}\n  \\footnotesize\\setstretch{1}\n"
           "  \\begin{tabular}{llll}\n    \\toprule\n    Moment & Data & Model & Role \\\\\n"
           "    \\midrule\n" + "\n".join(body) + "\n    \\bottomrule\n  \\end{tabular}\n"
           f"  \\tabnote{{{note}}}\n\\end{{table}}\n")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(tex)
    return out_path


def export_irfs(risk_path, tfp_path, sproc, pd_shock, tfp_dz, out_path, grid_note=""):
    # WRITE BOTH EXPERIMENTS TO ONE JSON FILE.
    risk = _pick(risk_path)
    # pd is recorded as a level (%/qtr); the schema wants the deviation from the rest point
    pd_rest = 100.0 * float(default_prob(sproc["s_star"]))
    risk["pd"] = (np.asarray(risk_path["pd"]) - pd_rest).tolist()
    out = {
        "method": "global",
        "label": "Global projection (nonlinear, occasionally binding)" + grid_note,
        "experiments": {
            "risk": {"shock": f"p^d {pd_rest:.2f}% -> {100 * pd_shock:.2f}%/qtr, "
                              f"decays through s at rho_s = {sproc['rho_s']}",
                     "series": risk},
            "tfp": {"shock": f"+{100 * tfp_dz:g}% Z_D, decays at rho_z = {sproc['rho_z']}",
                    "series": _pick(tfp_path)},
        },
    }
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1)
    return out_path
