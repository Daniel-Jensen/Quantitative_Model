# WHY THE CRISIS IRFS HAVE NOT RETURNED BY QUARTER 24 (2026-11-06, READ-ONLY; run.py's 12% box rules): THE EIGENVALUES
# OF THE EXACT LAW OF MOTION AT THE REST POINT, THE SLOW DIRECTIONS, AND THE 6% SHOCK'S IRF OVER 160 QUARTERS.
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

FREE = (0, 1, 2, 3, 4, 5, 6, 7, 10, 11)
T_LONG = 160

if __name__ == "__main__":
    from crisis_shock import setup12
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import STATE_NAMES, IS
    cal, ss, sproc, rules = setup12()
    S0 = RE.stochastic_rest_point(rules, cal, ss, sproc, verbose=False)
    idx = np.array(FREE)

    def g(S):
        return RE.advance(RE.read_exact(rules, cal, ss, sproc, S.copy()), S, sproc)[idx]

    # the exact law of motion's Jacobian at the rest point, in relative units (central differences)
    scale = np.maximum(np.abs(S0[idx]), 1e-2)
    J = np.empty((idx.size, idx.size))
    for j, i in enumerate(idx):
        h = 1e-5 * scale[j]
        Sp, Sm = S0.copy(), S0.copy()
        Sp[i] += h; Sm[i] -= h
        J[:, j] = (g(Sp) - g(Sm)) / (2 * h) * scale[j] / scale
    lam, V = np.linalg.eig(J)
    order = np.argsort(-np.abs(lam))
    print("EIGENVALUES OF THE EXACT LAW OF MOTION AT THE REST POINT (s at its mean)")
    for k in order:
        hl = np.log(0.5) / np.log(np.abs(lam[k])) if 0 < np.abs(lam[k]) < 1 else np.inf
        v = np.real(V[:, k]); v = v / np.max(np.abs(v))
        top = ", ".join(f"{STATE_NAMES[idx[i]]} {v[i]:+.2f}" for i in np.argsort(-np.abs(v))[:4])
        print(f"   |lambda| {np.abs(lam[k]):.5f} ({lam[k]:.4f})  half-life {hl:6.1f} q   direction: {top}")

    # the 6% shock over 160 quarters, exact, against the unshocked path (which stays at the rest point)
    o0 = RE.read_exact(rules, cal, ss, sproc, S0.copy())
    S = S0.copy(); S[IS] = RE.s_from_pd(0.06)
    rows = []
    for t in range(T_LONG):
        o = RE.read_exact(rules, cal, ss, sproc, S)
        rows.append((t, 100 * (o["Y_D"] / o0["Y_D"] - 1), 100 * (S[0] / S0[0] - 1), 100 * (o["I_D"] / o0["I_D"] - 1),
                     100 * (o["n_D"] / o0["n_D"] - 1), 100 * (o["C_D"] / o0["C_D"] - 1),
                     4e4 * (o["rdep_D"] - o0["rdep_D"]), 100 * (o["Bp_D"] / o0["Bp_D"] - 1),
                     100 * (o["nfa_dep_D"] - o0["nfa_dep_D"]) / ss["ss_firm_D"]["Y_ss"],
                     100 * (S[4] / S0[4] - 1)))
        S = rules.grid.clip(RE.advance(o, S, sproc))[0]
    print("\nTHE 6% SHOCK OVER 160 QUARTERS (exact; % vs the rest point, rates bp/yr, NFA in % of quarterly GDP)")
    print("   q      GDP       K       I     n_D     C_D    rdep    debt     NFA    b_DD")
    for r in rows:
        if r[0] in (0, 4, 8, 12, 16, 20, 24, 32, 40, 60, 80, 100, 120, 140, 159):
            print(f"  {r[0]:3d} " + " ".join(f"{x:+7.3f}" for x in r[1:]))
