# WHY THE D CAPITAL PRICE FALLS AFTER A CRISIS PURCHASE, AND WHAT THAT DOES TO D BANK NET WORTH (2026-11-06, READ-ONLY).
# Q_t (1 + r_t + psi_t) = mpk_{t+1} + (1 - delta) Q_{t+1} on the no-default path, psi the realised capital excess return;
# linearised against the unshocked path, dQ_0 splits into the deposit-rate path, the mpk path, the psi path and q24 on.
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

T = 25


def walk(rules, cal, ss, sproc, S0, s_shock=None, m0=None):
    # THE PATH WITH Q_K, mpk, rdep AND THE IMPACT BALANCE-SHEET LEGS.
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive import point_map as PM
    from global_projection.solver_recursive.state_grid import IS, IK_D, IBDD, IBFD, IO, IP_D
    from one_off_purchase import read_forced
    S = S0.copy()
    if s_shock is not None:
        S[IS] = s_shock
    out = {k: [] for k in ("Q", "mpk", "r", "legs")}
    for t in range(T):
        o = (read_forced(rules, cal, ss, sproc, S, m0) if (t == 0 and m0 is not None)
             else RE.read_exact(rules, cal, ss, sproc, S))
        v = PM._read_state_and_regime(S, 0, o["_x"], rules, cal, ss)
        PM._production(v, cal)
        PM._government(v, cal, ss)
        PM._balance_sheets(v, rules, cal)
        out["Q"].append(v.Q_D); out["mpk"].append(v.mpk_D); out["r"].append(v.rdep_D)
        if t == 0:
            payF = cal["delta_b_F"] + (1 - cal["delta_b_F"]) * v.Q_bF
            out["legs"] = dict(capital=(v.mpk_D + (1 - cal["delta_D"]) * v.Q_D) * S[IK_D],
                               D_bonds=v.Xi_D * S[IBDD], F_bonds=v.p * payF * S[IBFD], O=S[IO], P=S[IP_D],
                               ent=cal["omega_ent_D"] * v.assets_D, n=v.n_D)
        S = rules.grid.clip(RE.advance(o, S, sproc))[0]
    return {k: (np.array(v) if k != "legs" else v) for k, v in out.items()}


def decompose(a, ref, delta):
    # dQ_0 FROM THE LINEARISED NO-DEFAULT CAPITAL EULER, AGAINST THE UNSHOCKED PATH.
    def psi(p):
        return (p["mpk"][1:] + (1 - delta) * p["Q"][1:]) / p["Q"][:-1] - 1 - p["r"][:-1]
    Qr, R = ref["Q"], ref["r"][:-1] + psi(ref)
    dr = a["r"][:-1] - ref["r"][:-1]
    dpsi = psi(a) - psi(ref)
    dmpk = a["mpk"][1:] - ref["mpk"][1:]
    D = np.concatenate([[1.0], np.cumprod((1 - delta) / (1 + R))[:-1]])
    parts = dict(deposit_rate=-np.sum(D * Qr[:-1] * dr / (1 + R)),
                 risk_premium=-np.sum(D * Qr[:-1] * dpsi / (1 + R)),
                 mpk=np.sum(D * dmpk / (1 + R)),
                 beyond_q24=np.prod((1 - delta) / (1 + R)) * (a["Q"][-1] - Qr[-1]))
    parts = {k: 100 * v / Qr[0] for k, v in parts.items()}
    parts["sum"] = sum(parts.values())
    parts["actual"] = 100 * (a["Q"][0] / Qr[0] - 1)
    return parts, 4e4 * dr, 4e4 * dpsi


if __name__ == "__main__":
    from post_impact import setup
    from global_projection.solver_recursive import recursive_experiment as RE
    cal, ss, sproc, r_off, r_on, c_off, c_on, run = setup()
    rests = np.load(HERE / "runs" / "exact_rest_points.npy", allow_pickle=True).item()
    s_sh = RE.s_from_pd(run.RISK_SHOCK_PD)
    ref_off = walk(r_off, c_off, ss, sproc, rests["off"])
    ref_on = walk(r_on, c_on, ss, sproc, rests["on"])
    runs = {"A no TPI": (walk(r_off, c_off, ss, sproc, rests["off"], s_sh), ref_off),
            "B TPI": (walk(r_on, c_on, ss, sproc, rests["on"], s_sh), ref_on),
            "C no TPI + 41% buy": (walk(r_off, c_off, ss, sproc, rests["off"], s_sh, 0.4071), ref_off)}
    print("THE IMPACT CAPITAL PRICE, % vs the unshocked path, split by the no-default capital Euler")
    for lab, (a, ref) in runs.items():
        parts, dr, dpsi = decompose(a, ref, cal["delta_D"])
        print(f"  {lab:20s} " + "  ".join(f"{k} {v:+.4f}" for k, v in parts.items()))
        print(f"  {'':20s} deposit rate bp/yr q0..q8: {np.round(dr[:9], 2)}")
        print(f"  {'':20s} capital excess return bp/yr q0..q8: {np.round(dpsi[:9], 1)}")
    print("\nTHE IMPACT D BANK NET WORTH: n = (1 - f)(capital + D bonds + F bonds + O - P) + omega * assets")
    f = cal["f_D"]
    for lab, (a, ref) in runs.items():
        L, Lr = a["legs"], ref["legs"]
        n0 = Lr["n"]
        line = {k: 100 * (1 - f) * (L[k] - Lr[k]) / n0 for k in ("capital", "D_bonds", "F_bonds", "O")}
        line["deposit bill"] = -100 * (1 - f) * (L["P"] - Lr["P"]) / n0
        line["new bankers"] = 100 * (L["ent"] - Lr["ent"]) / n0
        line["total"] = 100 * (L["n"] / n0 - 1)
        print(f"  {lab:20s} " + "  ".join(f"{k} {v:+.3f}pp" for k, v in line.items()))
