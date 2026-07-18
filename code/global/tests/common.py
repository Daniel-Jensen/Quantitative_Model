"""Shared test fixtures: solve the steady state once per process."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402
from calibration import get_calibration  # noqa: E402
from steady_state import solve_steady_state  # noqa: E402

_CACHE = {}


def get_ss():
    """Calibration + steady state, solved once and cached for the process."""
    if "ss" not in _CACHE:
        cal = get_calibration()
        _CACHE["cal"] = cal
        _CACHE["ss"] = solve_steady_state(cal, verbose=False)
    return _CACHE["cal"], _CACHE["ss"]


def ss_input_paths(cal, ss):
    """Constant steady-state input paths for the bank block (length T)."""
    T = cal["T"]
    return dict(
        Kap_D=np.full(T, ss["Kap_D_ss"]), Kap_F=np.full(T, ss["Kap_F_ss"]),
        Q_D=np.ones(T), Q_F=np.ones(T),
        rk_D=np.full(T, ss["rk_D_ss"]), rk_F=np.full(T, ss["rk_F_ss"]),
        rdep_D=np.full(T, cal["r_dep_D_target"]),
        rdep_F=np.full(T, cal["r_dep_F_target"]),
        p_path=np.full(T, ss["p_ss"]),
    )


def transition_residuals(out, cal):
    """Max absolute residuals of the imposed and Walras-redundant markets.
    The IC is occasionally binding: cap_* is the complementarity product
    μ·slack (0 at an exact solution), not the old always-binding gap."""
    goods_D = np.max(np.abs(out["Y_D"] - out["P_CES_D"] * out["C_D"] - out["I_D"]
                            - out["NX_D"] - cal["G_D"]))
    goods_F = np.max(np.abs(out["Y_F"] - out["P_CES_F"] * out["C_F"] - out["I_F"]
                            - out["NX_F"] - cal["G_F"]))
    # union deposit market: one clearing (D-good units) + real-rate parity
    dep_union = np.max(np.abs((out["P_CES_D"] * out["A_D"] - out["Dep_supply_D"])
                              + out["p"] * (out["P_CES_F"] * out["A_F"]
                                            - out["Dep_supply_F"])))
    p_next = np.append(np.asarray(out["p"])[1:], np.asarray(out["p"])[-1])
    uip = np.max(np.abs((1.0 + np.asarray(out["rdep_D"]))
                        - (1.0 + np.asarray(out["rdep_F"])) * p_next / out["p"]))
    slack_D = np.asarray(out["alpha_D"]) * (np.asarray(out["n_D"]) - np.asarray(out["n_IC_D"]))
    slack_F = np.asarray(out["alpha_F"]) * (np.asarray(out["n_F"]) - np.asarray(out["n_IC_F"]))
    cap_D = np.max(np.abs(np.asarray(out["mu_D"]) * slack_D))
    cap_F = np.max(np.abs(np.asarray(out["mu_F"]) * slack_F))
    return dict(goods_D=goods_D, goods_F=goods_F, dep_union=dep_union, uip=uip,
                cap_D=cap_D, cap_F=cap_F,
                slack_min_D=float(np.min(slack_D)), slack_min_F=float(np.min(slack_F)),
                mu_min_D=float(np.min(out["mu_D"])), mu_min_F=float(np.min(out["mu_F"])))
