# RULE ACCURACY ALONG THE 9% CRISIS PATHS ON THE 12% BOX (2026-11-06, READ-ONLY): AT EACH QUARTER'S STATE, THE FITTED
# RULES (WHAT THE PREVIOUS QUARTER'S EXPECTATIONS READ) AGAINST THE EXACT SOLVE THERE, IN BOTH REGIMES, WITH AND WITHOUT
# THE ONE-OFF PURCHASE THAT SWITCHES THE D CONSTRAINT OFF (15.6% OF THE DEBT).
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

KEYS = ("Q_bD", "C_D", "alpha_D", "N_D", "Kp_D", "p")
T = 9

if __name__ == "__main__":
    from crisis_shock import setup12
    from one_off_purchase import read_forced
    from regime_anatomy import solve_regime, exact_val
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS, IM
    cal, ss, sproc, rules = setup12()
    S0 = np.load(HERE / "runs" / "crisis_rest_12pct.npy")
    for lab, m0 in (("no purchase", None), ("relief purchase 15.6%", 0.156)):
        S = S0.copy(); S[IS] = RE.s_from_pd(0.09)
        print(f"\n{lab}: fitted minus exact at each quarter's state (%), regime 0 | regime 1")
        print("   t    M%  " + " ".join(f"{k:>8s}" for k in KEYS) + "  | " + " ".join(f"{k:>8s}" for k in KEYS))
        for t in range(T):
            o = read_forced(rules, cal, ss, sproc, S, m0) if (t == 0 and m0) else RE.read_exact(rules, cal, ss, sproc, S)
            cells = []
            for d in (0, 1):
                ox = solve_regime(rules, cal, ss, sproc, S, d) if t > 0 or d == 1 else o
                cells.append(" ".join(f"{100*(float(rules.eval(k, d, np.atleast_2d(S))[0])/exact_val(ox, k)-1):+8.3f}"
                                      for k in KEYS))
            print(f"  {t:2d} {100*S[IM]/cal['B_gov_D_ss']:5.1f}  " + cells[0] + "  | " + cells[1])
            S = rules.grid.clip(RE.advance(o, S, sproc))[0]
