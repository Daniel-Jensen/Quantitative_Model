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
