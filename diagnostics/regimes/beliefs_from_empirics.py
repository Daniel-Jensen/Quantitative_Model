"""
beliefs_from_empirics.py — map the monthly 3-state MS estimation
(Empirics/outputs/ms_regime_COMPOSITE.npz: states ordered by mean spread,
0=dove/intervention, 1=base, 2=hawk/stress) to quarterly beliefs over POLICY
types {aggressive, medium, passive}.

Mapping (spec §10.4, judgment call, swappable): dove->aggressive, base->medium,
hawk-persistence->passive. NAMING: aggressive/medium/passive are the policy
regimes; dove/base/hawk are empirical market regimes — never mix (spec header).
"""
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(os.path.dirname(os.path.dirname(HERE)),
                   "Empirics", "outputs", "ms_regime_COMPOSITE.npz")

def main():
    d = np.load(NPZ)
    P_m, means, erg = d["P"], d["means"], d["ergodic"]
    assert np.all(np.diff(means) > 0), "states must be mean-ordered dove<base<hawk"
    assert np.allclose(P_m.sum(axis=1), 1.0, atol=1e-10), "P_m rows must sum to 1"
    P_q = np.linalg.matrix_power(P_m, 3)          # monthly -> quarterly

    # Onset beliefs: conditional distribution of next quarter's regime given the
    # economy sits in the stress (hawk) regime today — mapped to policy types.
    hawk_row_q = P_q[2]                            # (dove, base, hawk)
    beliefs = {
        "pi_onset_1q":  {"aggressive": float(hawk_row_q[0]), "medium": float(hawk_row_q[1]),
                         "passive": float(hawk_row_q[2])},
        "pi_ergodic":   {"aggressive": float(erg[0]), "medium": float(erg[1]),
                         "passive": float(erg[2])},
        "pi_onset_k":   {str(k): {n: float(v) for n, v in
                         zip(("aggressive", "medium", "passive"),
                             np.linalg.matrix_power(P_m, 3 * k)[2])} for k in (1, 2, 4)},
        "P_quarterly":  P_q.tolist(),
        "provenance":   "ms_regime_COMPOSITE.npz (FRED 10y peripheral-Bund composite, "
                        "1995-2026 monthly); P_q = P_m^3; mapping dove->aggressive, "
                        "base->medium, hawk->passive (judgment call, spec §10.4)",
    }
    for name, pi in [("pi_onset_1q", beliefs["pi_onset_1q"]), ("pi_ergodic", beliefs["pi_ergodic"])]:
        s = sum(pi.values())
        assert abs(s - 1.0) < 1e-8, f"{name} does not sum to 1: {s}"
        print(f"{name}: " + ", ".join(f"{k}={v:.4f}" for k, v in pi.items()))
    out = os.path.join(HERE, "beliefs.json")
    json.dump(beliefs, open(out, "w"), indent=2)
    print(f"written: {out}")

if __name__ == "__main__":
    main()
