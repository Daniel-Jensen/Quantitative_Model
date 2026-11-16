# THE CAP AND THE SHOCK ON THE FASTER-FISCAL-RULE RULES (2026-11-06, READ-ONLY; runs/rules_12pct_bohn0.250.pkl):
# THE SLOWEST ROOTS AND A 60-QUARTER NO-TPI RESPONSE (DOES IT NOW RETURN?), THEN, FOR EACH SHOCK AND CAP, THE NO-TPI
# RECESSION, THE PURCHASE THAT HOLDS THE CAP ON IMPACT AND WHETHER IT FITS THE EUROSYSTEM-BOOK BOX.
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

SHOCKS = (0.05, 0.055, 0.06, 0.07)
CAPS = (450.0, 475.0, 500.0, 525.0)
FREE = (0, 1, 2, 3, 4, 5, 6, 7, 10, 11)

if __name__ == "__main__":
    from crisis_shock import setup12
    from bohn_homotopy import economy_at
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS, IM
    base, _, _, _ = setup12()
    cal, ss, sproc = economy_at(base, 0.25)
    rules = pickle.load(open(HERE / "runs" / "rules_12pct_bohn0.250.pkl", "rb"))
    S0 = RE.stochastic_rest_point(rules, cal, ss, sproc, verbose=True)
    sov = lambda o: 4e4 * (cal["delta_b_D"] * (1 / o["Q_bD"] - 1) - cal["delta_b_F"] * (1 / o["Q_bF"] - 1))

    # 1. the slowest roots of the exact law of motion at the rest point
    idx = np.array(FREE)
    scale = np.maximum(np.abs(S0[idx]), 1e-2)
    g = lambda S: RE.advance(RE.read_exact(rules, cal, ss, sproc, S.copy()), S, sproc)[idx]
    J = np.empty((idx.size, idx.size))
    for j, i in enumerate(idx):
        h = 1e-5 * scale[j]
        Sp, Sm = S0.copy(), S0.copy()
        Sp[i] += h; Sm[i] -= h
        J[:, j] = (g(Sp) - g(Sm)) / (2 * h) * scale[j] / scale
    lam = np.sort(np.abs(np.linalg.eigvals(J)))[::-1]
    print("slowest roots:", ", ".join(f"{l:.4f} (half-life {np.log(0.5) / np.log(l):.1f} q)" for l in lam[:4]))

    # 2. the 6% no-TPI response over 60 quarters
    o0 = RE.read_exact(rules, cal, ss, sproc, S0.copy())
    S = S0.copy(); S[IS] = RE.s_from_pd(0.06)
    print("\n6% shock, no TPI (exact; % vs rest):   q   GDP      K      I    debt    n_D   spread")
    for t in range(61):
        o = RE.read_exact(rules, cal, ss, sproc, S)
        if t in (0, 1, 2, 4, 8, 12, 16, 20, 24, 32, 40, 60):
            print(f"   {t:3d} {100*(o['Y_D']/o0['Y_D']-1):+7.3f} {100*(S[0]/S0[0]-1):+7.3f} {100*(o['I_D']/o0['I_D']-1):+7.3f}"
                  f" {100*(o['Bp_D']/o0['Bp_D']-1):+7.3f} {100*(o['n_D']/o0['n_D']-1):+7.2f} {sov(o):7.1f}")
        S = rules.grid.clip(RE.advance(o, S, sproc))[0]

    # 3. shock x cap: the no-TPI impact and the purchase that holds the cap
    edge = rules.grid.hi[IM] / cal["B_gov_D_ss"]
    print(f"\nrest spread {sov(o0):.1f} bp; the book box ends at {100*edge:.1f}% of the debt")
    for pd in SHOCKS:
        Sx = S0.copy(); Sx[IS] = RE.s_from_pd(pd)
        ox = RE.read_exact(rules, cal, ss, sproc, Sx)
        dY = 100 * (ox["Y_D"] / o0["Y_D"] - 1)
        cells = []
        for cap in CAPS:
            try:
                m = RE.cap_purchase(rules, cal, ss, sproc, Sx, cap)
            except RuntimeError:
                cells.append(f"{cap:.0f}: unreachable")
                continue
            if m <= 0:
                cells.append(f"{cap:.0f}: not binding")
                continue
            oc = RE.read_forced(rules, cal, ss, sproc, Sx, m)
            cells.append(f"{cap:.0f}: buy {100*m:4.1f}%{'*' if m > edge else ''} GDP {100*(oc['Y_D']/o0['Y_D']-1):+.3f}")
        print(f"  p^d {100*pd:4.1f}%: no TPI spread {sov(ox):6.1f} bp, GDP {dY:+.3f}%, mu_D {ox['mu_D']:.3f} | " + " | ".join(cells))
    print("  (* past the book box edge)")
