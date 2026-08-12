# DIAGNOSTIC HELPERS FOR THE RECURSIVE SOLUTION (DELIVERABLE C, TIME-ITERATION).
# The collocation residual is assembled POINTWISE in point_map.point_residuals
# (the per-period image of transition.py's stacked system) and driven to a fixed
# point by recursive_main.time_iteration -- a Coleman/time-iteration operator, not
# one monolithic Newton over all grid points (which is far stiffer at the
# occasionally-binding kink; the earlier make_collocation_residual approach was
# superseded once the pointwise structure landed). This module keeps the
# grid-wide diagnostics the gates and the E-experiment read off the solved rules.
import numpy as np

from solver_recursive.decision_rules import SOLVE7
from solver_recursive.point_map import point_residuals


def grid_diagnostics(rules, cal, ss, sproc, regime=0, no_default=False, n_gh=7):
    # PER-POINT mu/slack/Y_D/C_D/residual MAPS FOR GATES N3 AND THE E DIAGNOSTIC.
    n = rules.grid.n
    keys = ("mu_D", "mu_F", "slack_D", "slack_F", "n_D", "n_F",
            "Y_D", "C_D", "I_D", "Q_bD", "alpha_D")
    out = {k: np.empty(n) for k in keys}
    resid = np.empty(n)
    for i in range(n):
        S = rules.grid.points[i]
        x = np.array([rules.vals[k][regime][i] for k in SOLVE7])
        r, o = point_residuals(S, regime, x, rules, cal, ss, sproc,
                               n_gh=n_gh, no_default=no_default)
        resid[i] = np.max(np.abs(r))
        for k in keys:
            out[k][i] = o[k]
    out["resid"] = resid
    return out
