"""Verify the switchable tax_rule promotion on the DEFAULT calibration.

(1) mv_rule=0 must reproduce current behavior: default-shock n_inter_D[0]<0 and
    Y_D[0]<0, and Walras residuals (goods_mkt_F, ca_res_D) at solver tolerance.
(2) mv_rule=1 must build/solve and produce finite IRFs with mv_gov_ss set.
"""
import sys, copy
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))

from calibration import get_calibration
from steady_state import solve_steady_state
from ic_delta_calibration import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model import build_and_solve

ss = get_calibration()
ss = solve_steady_state(ss); ss = calibrate_ic_delta(ss); ss = calibrate_depreciation(ss)
base = build_and_solve(ss)                       # mv_rule=0 (default)
G0 = base['G']; ha = base['ha_full']; ssf = base['ss_final']
un = base['unknowns_tp']; tg = base['targets_tp']; T = base['T']
exo = ['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F']
dS = 0.01 * 0.8 ** np.arange(T); z = np.zeros(T)

print(f"\n[mv_gov_ss check] mv_gov_ss_D={ssf.toplevel['mv_gov_ss_D']:.4f} "
      f"(= q_b_D_ss*b_gov_ss_D = {float(ssf['q_b_D'])*float(ssf['b_gov_ss_D']):.4f})")

d0 = G0 @ {'Z_D': z, 'Z_F': z, 'shock_def_D': dS, 'shock_def_F': z}
print("\n=== mv_rule=0 (par value, default) — regression check ===")
print(f"  n_inter_D[0] = {d0['n_inter_D'][0]*100:+.4f}%   (must be < 0)")
print(f"  Y_D[0]       = {d0['Y_D'][0]*100:+.4f}%   (must be < 0)")
print(f"  max|goods_mkt_F| = {np.max(np.abs(d0['goods_mkt_F'])):.2e}   (target <=1e-7)")
print(f"  max|ca_res_D|    = {np.max(np.abs(d0['ca_res_D'])):.2e}   (target <=1e-7)")
ok0 = (d0['n_inter_D'][0] < 0 and d0['Y_D'][0] < 0
       and np.max(np.abs(d0['goods_mkt_F'])) < 1e-7 and np.max(np.abs(d0['ca_res_D'])) < 1e-7)
print(f"  REGRESSION {'PASS' if ok0 else 'FAIL'}")

# mv_rule=1: re-solve (mv_gov_ss already set by build_and_solve)
st = copy.deepcopy(ssf)
st.toplevel['mv_rule_D'] = 1.0; st.toplevel['mv_rule_F'] = 1.0
G1 = ha.solve_jacobian(st, unknowns=un, targets=tg, inputs=exo, T=T)
d1 = G1 @ {'Z_D': z, 'Z_F': z, 'shock_def_D': dS, 'shock_def_F': z}
finite = np.all(np.isfinite(d1['n_inter_D'])) and np.all(np.isfinite(d1['b_gov_D']))
differs = np.max(np.abs(np.asarray(d1['T_ls_D']) - np.asarray(d0['T_ls_D']))) > 1e-9
print("\n=== mv_rule=1 (market value) — switch works ===")
print(f"  finite IRFs: {finite}   T_ls path differs from par: {differs}")
print(f"  n_inter_D[0] = {d1['n_inter_D'][0]*100:+.4f}%   max|goods_mkt_F|={np.max(np.abs(d1['goods_mkt_F'])):.2e}")
print(f"  SWITCH {'PASS' if (finite and differs) else 'FAIL'}")
print("\nPROMOTION VERIFY COMPLETE")
