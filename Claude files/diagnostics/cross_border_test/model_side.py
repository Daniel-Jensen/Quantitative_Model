# THE MODEL SIDE OF THE CROSS-BORDER TEST (2026-10-30, READ-ONLY): F BANKS' HOLDINGS OF D BONDS AND THE EXCESS RETURN THEY SEE,
# BY RISK LEVEL; THEN THE MODEL AT THE DEFAULT PROBABILITIES WHOSE SOVEREIGN SPREAD MATCHES THE DATA (2010Q1 314 bp, 2010Q4
# 910 bp), FOR THE HEAD-TO-HEAD WITH compare.py.
#   python3 -u model_side.py
import os, pickle, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT); sys.path.insert(0, str(ROOT / "Claude files/diagnostics/tex_calibration"))
import numpy as np
from scipy.optimize import brentq


def main():
    from tpi_purchase_tests import setup
    from global_projection.solver_recursive import recursive_experiment as RE
    from global_projection.solver_recursive.state_grid import IS
    cal, ss, sp, r = setup()
    S0 = RE.stochastic_rest_point(r, cal, ss, sp, verbose=False)
    B, psi, ref = cal["B_gov_D_ss"], cal["psi_bD_F"], cal["b_D_F_ref"]

    def read(pd):
        # THE EXACT READ AT p^d (None: THE REST POINT).
        S = S0.copy()
        if pd:
            S[IS] = RE.s_from_pd(pd)
        return RE.read_exact(r, cal, ss, sp, S)
    sov = lambda o: 4e4 * (cal["delta_b_D"] * (1 / o["Q_bD"] - 1) - cal["delta_b_F"] * (1 / o["Q_bF"] - 1))
    rest = read(None)
    print(f"psi_bD_F {psi}: (b_DF - b_ref)/B = (F banks' excess return on D bonds)/psi")
    for pd in (None, 0.02, 0.05, 0.12):
        o = read(pd)
        print(f"   p^d {'rest' if pd is None else f'{100 * pd:.0f}%':>5s}: F book {100 * (o['b_D_F_new'] / rest['b_D_F_new'] - 1):+6.2f}% "
              f"vs rest, foreign share of D stock {100 * o['b_D_F_new'] / o['Bp_D']:5.2f}%, F excess return "
              f"{4e4 * psi * (o['b_D_F_new'] - ref) / B:+7.0f} bp/yr, sovereign spread {sov(o):6.0f} bp")
    rows = {}
    for tag, target in (("2010Q1", 314.0), ("2010Q4", 910.0)):
        pd = brentq(lambda p: sov(read(p)) - target, 0.002, 0.12, xtol=1e-5)
        rows[tag] = o = read(pd)
        print(f"{tag}: model p^d {100 * pd:.2f}%, German share of German + Greek bank holdings "
              f"{100 * o['b_D_F_new'] / (o['b_D_F_new'] + o['b_D_D_new']):.1f}%, price {o['Q_bD']:.3f}")
    a, b = rows["2010Q1"], rows["2010Q4"]
    val = lambda k: 100 * (b[k] * b["Q_bD"] / (a[k] * a["Q_bD"]) - 1)
    print(f"2010Q1 -> 2010Q4: German book {100 * (b['b_D_F_new'] / a['b_D_F_new'] - 1):+.1f}% (quantity), {val('b_D_F_new'):+.1f}% "
          f"(value); Greek banks' book {val('b_D_D_new'):+.1f}% (value); price {100 * (b['Q_bD'] / a['Q_bD'] - 1):+.1f}%")


if __name__ == "__main__":
    main()
