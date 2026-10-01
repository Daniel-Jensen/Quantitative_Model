"""Write the SSJ IRFs in the comparison schema that compare.py reads.

Every series is a percent deviation from its own steady-state level (SSJ IRFs are
LEVEL deviations, so dividing by the SS level is required -- n_inter_D_ss = 2.138 and
K_D_ss = 10.8 are not ~1). `pd` is the quarterly default probability in pp.
"""
import json
from pathlib import Path

import numpy as np

# comparison name -> SSJ variable
SERIES = {"Y_D": "Y_D", "C_D": "C_D", "I_D": "I_D", "n_D": "n_inter_D", "Q_bD": "q_b_D",
          "Y_F": "Y_F", "C_F": "C_F", "I_F": "I_F", "n_F": "n_inter_F"}
HORIZON = 40


def _ss_level(ss, k):
    # I_D/I_F have no SS entry of their own; SS investment is delta*K (as in dump_irfs.py)
    if k in ("I_D", "I_F"):
        return float(ss[f"delta_{k[-1]}"]) * float(ss[f"K_{k[-1]}"])
    return float(ss[k])


def _pct(irfs, ss):
    return {name: (np.asarray(irfs[k][:HORIZON]) * 100.0 / _ss_level(ss, k)).tolist()
            for name, k in SERIES.items() if k in irfs}


def export_irfs(model_results, path):
    ss = model_results["ss_final"]
    irf_def = model_results["irfs_def_D"]
    risk = _pct(irf_def, ss)
    # the equilibrium default rate if the model reports it, else the exogenous shock path
    pd = irf_def["def_rate_D"] if "def_rate_D" in irf_def else model_results["dShock_def_D"]
    risk["pd"] = (100.0 * np.asarray(pd[:HORIZON])).tolist()
    out = {
        "method": "ssj",
        "label": "Sequence-space (linear, sticky prices)",
        "experiments": {
            "risk": {"shock": "+1pp quarterly default rate (shock_def_D), AR(1) decay",
                     "series": risk},
            "tfp": {"shock": "+1% Z_D, AR(1) decay",
                    "series": _pct(model_results["irfs_Z_D"], ss)},
        },
    }
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(out, indent=1))
    return path
