"""Corrected re-run of philamb_sweep_postC1_fine.py.

That script's Prony(p=8) estimate gave an implausible modulus (~3.2) that
contradicted the tiny observed b_gov_D[499] levels. philamb_order_selection_check.py
diagnosed why: R^2 already saturates to machine precision at order p=2-3, and
every order beyond that fits pure numerical noise into spurious high-modulus
"ghost" poles (contribute negligibly to the actual signal, but have huge
eigenvalues in the fitted companion matrix) -- classic AR/Prony overfitting.
Order p=4 (used as the "conservative" cross-check in the first fine sweep) was
ALREADY past the onset of this overfitting for at least one grid point.

This version fits orders p=1..4 for every cell, reports R^2 alongside modulus
so the saturation point is visible and auditable per-cell (not assumed from a
single reference point), and takes the modulus at the SMALLEST order whose R^2
exceeds 0.99999 as the trusted estimate (falling back to p=3 if none clears
the bar, which would itself be a red flag worth reporting).
"""
import sys, json, copy
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))

from calibration              import get_calibration
from steady_state             import solve_steady_state
from ic_delta_calibration     import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model                import build_and_solve

PHI_GRID = [0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25]
MV_RULE  = 1.0
R2_BAR   = 0.99999
rho = 0.8


def dom_modulus_energy(x):
    x = np.asarray(x)
    e1 = np.sum(x[300:400] ** 2)
    e2 = np.sum(x[400:500] ** 2)
    return float((e2 / e1) ** (1 / 200)) if e1 > 1e-30 else 0.0


def fit_order(x, p, t0=150, t1=500):
    x = np.asarray(x, dtype=float)
    seg = x[t0:t1]
    scale = np.max(np.abs(seg))
    if scale < 1e-14:
        return 0.0, 1.0, np.array([0.0])
    seg_n = seg / scale
    n = len(seg_n)
    X = np.column_stack([seg_n[p - 1 - k: n - 1 - k] for k in range(p)])
    y = seg_n[p:n]
    a, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ a
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    C = np.zeros((p, p)); C[0, :] = a
    for i in range(1, p):
        C[i, i - 1] = 1.0
    moduli = np.sort(np.abs(np.linalg.eigvals(C)))[::-1]
    return moduli[0], r2, moduli


def trusted_modulus(x, orders=(1, 2, 3, 4), r2_bar=R2_BAR):
    """Smallest order clearing r2_bar; falls back to p=3 (flagged) if none does."""
    fits = {p: fit_order(x, p) for p in orders}
    for p in orders:
        mod, r2, moduli = fits[p]
        if r2 >= r2_bar:
            return p, mod, r2, fits
    # none cleared the bar -- fall back, flag it
    p_fb = 3
    mod, r2, moduli = fits[p_fb]
    return p_fb, mod, r2, fits


print("=== SS solve: today's committed calibration, delta_b -> empirical long "
      "duration (0.036/0.038) [identical setup to philamb_sweep_postC1.py] ===", flush=True)
cali = get_calibration()
cali['delta_b_D'] = 0.036
cali['delta_b_F'] = 0.038
cali['q_b_D'] = 0.93
cali['q_b_F'] = 0.93

ss = solve_steady_state(cali)
ss = calibrate_ic_delta(ss)
ss = calibrate_depreciation(ss)

base        = build_and_solve(ss)
ha_full     = base['ha_full']
ss_final    = base['ss_final']
unknowns_tp = base['unknowns_tp']
targets_tp  = base['targets_tp']
T           = base['T']
exogenous   = ['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F']

dZ   = 0.01 * rho ** np.arange(T)
zero = np.zeros(T)

RESULTS = {}
for phi in PHI_GRID:
    key = f"phi{phi}"
    print(f"--- long_empirical duration | mv_rule=1  phi_lamb={phi} ---", flush=True)
    try:
        ss_try = copy.deepcopy(ss_final)
        for c in ['D', 'F']:
            ss_try.toplevel[f'phi_lamb_{c}'] = phi
            ss_try.toplevel[f'mv_rule_{c}']  = MV_RULE
        G = ha_full.solve_jacobian(ss_try, unknowns=unknowns_tp, targets=targets_tp,
                                    inputs=exogenous, T=T)
        irf_d = G @ {'Z_D': zero, 'Z_F': zero, 'shock_def_D': dZ, 'shock_def_F': zero}
        irf_z = G @ {'Z_D': dZ, 'Z_F': zero, 'shock_def_D': zero, 'shock_def_F': zero}
        bgov_d = np.asarray(irf_d['b_gov_D'])
        bgov_z = np.asarray(irf_z['b_gov_D'])

        e_mod_d = dom_modulus_energy(bgov_d)
        e_mod_z = dom_modulus_energy(bgov_z)
        p_d, mod_d, r2_d, fits_d = trusted_modulus(bgov_d)
        p_z, mod_z, r2_z, fits_z = trusted_modulus(bgov_z)
        mod = max(mod_d, mod_z)

        m = {
            'energy_mod_def': e_mod_d, 'energy_mod_Z': e_mod_z,
            'trusted_order_def': p_d, 'trusted_mod_def': float(mod_d), 'trusted_r2_def': r2_d,
            'trusted_order_Z':   p_z, 'trusted_mod_Z':   float(mod_z), 'trusted_r2_Z':   r2_z,
            'r2_by_order_def': {str(p): fits_d[p][1] for p in (1, 2, 3, 4)},
            'mod_by_order_def': {str(p): float(fits_d[p][0]) for p in (1, 2, 3, 4)},
            'bgov_D_499_def': float(bgov_d[499]), 'bgov_D_499_Z': float(bgov_z[499]),
            'n_inter_D_0_def': float(irf_d['n_inter_D'][0]),
        }
        RESULTS[key] = m
        v_energy = 'STABLE' if max(e_mod_d, e_mod_z) < 0.999 else 'EXPLOSIVE'
        v_trust  = 'STABLE' if mod < 0.999 else 'EXPLOSIVE'
        print(f"  energy: {v_energy} (mod={max(e_mod_d,e_mod_z):.4f})   "
              f"trusted (order={max(p_d,p_z)}): {v_trust} (mod={mod:.4f}, R2_def={r2_d:.9f})", flush=True)
        print(f"  by-order mod (def shock): "
              + ", ".join(f"p{p}={fits_d[p][0]:.4f}(R2={fits_d[p][1]:.6f})" for p in (1, 2, 3, 4)), flush=True)
    except Exception as e:
        RESULTS[key] = {'error': repr(e)}
        print(f"  FAILED: {e!r}", flush=True)

out = ROOT / 'audit_artifacts' / 'philamb_sweep_postC1_fine_v2_results.json'
with open(out, 'w') as fh:
    json.dump(RESULTS, fh, indent=2, default=str)

print("\n" + "=" * 110)
print("FINE phi_lamb GRID v2 (mv_rule=1, empirical duration): energy-ratio proxy vs. order-selected Prony estimate")
print("=" * 110)
print(f"{'phi_lamb':>9} {'energy_v':>10} {'e_mod':>8} {'trust_v':>10} {'t_mod':>8} {'order':>6} {'R2':>12} {'bgov[499]':>11}")
for phi in PHI_GRID:
    m = RESULTS.get(f"phi{phi}", {})
    if 'error' in m:
        print(f"{phi:>9} ERROR: {m['error'][:80]}")
        continue
    e_mod = max(m['energy_mod_def'], m['energy_mod_Z'])
    t_mod = max(m['trusted_mod_def'], m['trusted_mod_Z'])
    order = max(m['trusted_order_def'], m['trusted_order_Z'])
    v_e = 'STABLE' if e_mod < 0.999 else 'EXPLOSIVE'
    v_t = 'STABLE' if t_mod < 0.999 else 'EXPLOSIVE'
    print(f"{phi:>9} {v_e:>10} {e_mod:>8.4f} {v_t:>10} {t_mod:>8.4f} {order:>6} "
          f"{m['trusted_r2_def']:>12.8f} {m['bgov_D_499_def']:>11.4g}")
print("=" * 110)
print(f"Results -> {out}")
print("FINE SWEEP v2 COMPLETE")
