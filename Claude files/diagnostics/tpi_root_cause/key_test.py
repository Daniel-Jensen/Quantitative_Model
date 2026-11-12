# DOES THE REMITTANCE KEY CAUSE THE POST-IMPACT DAMAGE OF A CRISIS PURCHASE? (2026-11-06) THE NO-TPI ECONOMY ON THE
# APPENDIX-B 12% BOX RE-SOLVED WITH THE EUROSYSTEM'S P&L KEPT IN D (tpi_key_D = 1), WARM FROM data/global_warm_start.pkl,
# THEN THE ONE-OFF PURCHASE THAT SWITCHES THE D CONSTRAINT OFF, SIDE BY SIDE WITH THE CAPITAL-KEY RULES (0.071).
#   python3 -u key_test.py        (a collocation solve: run detached; workers re-import this file)
import copy
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

OUT = HERE / "runs" / "rules_12pct_key1.pkl"
SHOW = (0, 1, 2, 4, 8, 12)


def max_resid(rules, cal, ss, sproc):
    # max|F| OF THE COLLOCATION SYSTEM AT THE RULES' OWN NODE VALUES.
    from global_projection.solver_recursive.collocation import make_residual, pack
    F = make_residual(copy.deepcopy(rules), cal, ss, sproc)
    return float(np.max(np.abs(F(pack(rules)))))


def m_switch(rules, cal, ss, sproc, S):
    # THE ONE-OFF PURCHASE (SHARE OF THE DEBT) AT WHICH THE D CONSTRAINT STOPS BINDING, RULES HELD FIXED.
    from scipy.optimize import brentq
    from one_off_purchase import read_forced
    g = lambda m: read_forced(rules, cal, ss, sproc, S, m)["mu_D"] - 1e-6
    return brentq(g, 0.0, 0.4, xtol=1e-5) if g(0.0) > 0 else 0.0


if __name__ == "__main__":
    from crisis_shock import setup12, walk, dev
    from exact_rest import exact_rest
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    cal, ss, sproc, r071 = setup12()
    print(f"warm rules at key_D 0.071: max|F| {max_resid(r071, cal, ss, sproc):.2e}", flush=True)
    if OUT.exists():
        r1 = pickle.load(open(OUT, "rb"))
    else:
        r1 = copy.deepcopy(r071)
        c1 = dict(cal, tpi_key_D=1.0)
        print(f"same rules at key_D 1: max|F| {max_resid(r1, c1, ss, sproc):.2e}", flush=True)
        ok, its, w = RE._stage(r1, c1, ss, sproc, (0, 1), False, "key_D 1", True, warm=0, maxit=20)
        print(f"key_D 1 solve: ok {ok}, {its} iterations, max|F| {w:.2e}", flush=True)
        if not ok:
            sys.exit(1)
        with open(OUT, "wb") as fh:
            pickle.dump(r1, fh)
    for key, rules in ((0.071, r071), (1.0, r1)):
        c = dict(cal, tpi_key_D=key)
        S0, err = exact_rest(rules, c, ss, sproc, RE.stochastic_rest_point(rules, c, ss, sproc, verbose=False))
        ref = walk(rules, c, ss, sproc, S0)
        for pd in (0.09, 0.06):
            S = S0.copy(); S[IS] = RE.s_from_pd(pd)
            m = m_switch(rules, c, ss, sproc, S)
            a = dev(walk(rules, c, ss, sproc, S0, RE.s_from_pd(pd)), ref)
            b = dev(walk(rules, c, ss, sproc, S0, RE.s_from_pd(pd), m), ref)
            print(f"\nkey_D {key}: p^d {100*pd:.0f}%, one-off purchase {100*m:.1f}% (rest sov {ref['sov'][0]:.1f} bp, "
                  f"n_D {ref['n'][0]:.4f}):   no purchase | purchase, quarters {SHOW}", flush=True)
            for k in ("sov", "n", "Y", "C", "I", "N", "mu", "lend", "QK", "Pi"):
                print(f"   {k:5s} " + " ".join(f"{a[k][t]:+8.3f}" for t in SHOW) + "  |  "
                      + " ".join(f"{b[k][t]:+8.3f}" for t in SHOW), flush=True)
