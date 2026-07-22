"""Order-selection diagnostic: the fine-grid Prony(p=8) estimate produced an
implausible modulus (~3.2) that contradicts the tiny observed b_gov_D[499]
levels (~1e-6) -- a classic overfitting symptom (excess model order absorbs
numerical noise into spurious high-modulus "ghost" poles with negligible
actual contribution to the signal). This script re-solves ONE representative
cell (phi_lamb=0.10, mv_rule=1, same setup as philamb_sweep_postC1_fine.py),
saves the raw b_gov_D IRF, and checks fit quality (R^2) vs. AR order p=1..12
to find the order actually supported by the data, instead of guessing p=4/p=8.
"""
import sys, copy
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))

from calibration              import get_calibration
from steady_state             import solve_steady_state
from ic_delta_calibration     import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model                import build_and_solve

cali = get_calibration()
cali['delta_b_D'] = 0.036; cali['delta_b_F'] = 0.038
cali['q_b_D'] = 0.93; cali['q_b_F'] = 0.93

ss = solve_steady_state(cali)
ss = calibrate_ic_delta(ss)
ss = calibrate_depreciation(ss)
base = build_and_solve(ss)
ha_full, ss_final = base['ha_full'], base['ss_final']
unknowns_tp, targets_tp, T = base['unknowns_tp'], base['targets_tp'], base['T']
exogenous = ['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F']

ss_try = copy.deepcopy(ss_final)
for c in ['D', 'F']:
    ss_try.toplevel[f'phi_lamb_{c}'] = 0.10
    ss_try.toplevel[f'mv_rule_{c}']  = 1.0
G = ha_full.solve_jacobian(ss_try, unknowns=unknowns_tp, targets=targets_tp, inputs=exogenous, T=T)

rho = 0.8
dZ, zero = 0.01 * rho ** np.arange(T), np.zeros(T)
irf_d = G @ {'Z_D': zero, 'Z_F': zero, 'shock_def_D': dZ, 'shock_def_F': zero}
bgov_d = np.asarray(irf_d['b_gov_D'])
np.savez(ROOT / 'audit_artifacts' / 'philamb_order_check_bgov.npz', bgov_d=bgov_d)
print(f"b_gov_D[499] = {bgov_d[499]:.6e}   b_gov_D[150] = {bgov_d[150]:.6e}   b_gov_D[300] = {bgov_d[300]:.6e}")


def fit_order(x, p, t0=150, t1=500):
    x = np.asarray(x, dtype=float)
    seg = x[t0:t1]
    scale = np.max(np.abs(seg))
    seg_n = seg / scale
    n = len(seg_n)
    X = np.column_stack([seg_n[p - 1 - k: n - 1 - k] for k in range(p)])
    y = seg_n[p:n]
    a, res, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ a
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    C = np.zeros((p, p)); C[0, :] = a
    for i in range(1, p):
        C[i, i - 1] = 1.0
    moduli = np.sort(np.abs(np.linalg.eigvals(C)))[::-1]
    n_obs = len(y)
    aic = n_obs * np.log(ss_res / n_obs + 1e-300) + 2 * p
    return r2, aic, moduli, rank


print(f"\n{'p':>3} {'R^2':>12} {'AIC':>10} {'rank':>5} {'max_modulus':>12}  top-4 moduli")
print("-" * 90)
for p in range(1, 13):
    r2, aic, moduli, rank = fit_order(bgov_d, p)
    print(f"{p:>3} {r2:>12.9f} {aic:>10.2f} {rank:>5} {moduli[0]:>12.4f}  {np.round(moduli[:4], 4)}")
