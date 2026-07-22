"""Finer re-test of the mv_rule=1 "explosive island" found in philamb_sweep_postC1.py.

That sweep used a crude two-window energy-ratio proxy for the dominant modulus
(compare sum-of-squares over t=[300,400) vs t=[400,500), take the 1/200-th
power). That proxy is known to be fooled by oscillatory near-unit-circle modes
(complex conjugate pole pairs can make windowed energy fluctuate independent of
the true decay/growth rate). This script replaces it with a proper estimate of
the dominant pole modulus via linear-prediction / Prony analysis:

  1. Take the tail of the b_gov_D IRF (t in [150, 500), well past the initial
     shock transient).
  2. Fit an order-p AR recursion x[t] = sum_k a_k x[t-k] by least squares.
  3. Build the companion matrix from the fitted a_k and take its eigenvalues --
     these ARE the (estimated) poles of the underlying linear system restricted
     to this observable, and max(|eigenvalue|) is the proper dominant modulus,
     correctly separating oscillatory modes from a simple energy trend.

Fit at two orders (p=4, p=8) as an internal robustness check -- if they
disagree sharply, the estimate itself is suspect and needs a longer horizon or
different observable, not just a bigger p.

Also refines the phi_lamb grid to step 0.025 over [0.05, 0.25] (mv_rule=1 only;
mv_rule=0 was uniformly explosive on the coarse grid, so refining it adds
little). Reuses one SS solve (fiscal-rule params are SS-invariant) and one
ha_full build; only the Jacobian is re-solved per phi_lamb cell.
"""
import os, sys, json, copy
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
rho = 0.8


def dom_modulus_energy(x):
    """Original coarse two-window energy-ratio proxy (kept for comparison)."""
    x = np.asarray(x)
    e1 = np.sum(x[300:400] ** 2)
    e2 = np.sum(x[400:500] ** 2)
    return float((e2 / e1) ** (1 / 200)) if e1 > 1e-30 else 0.0


def dominant_pole_modulus(x, p=8, t0=150, t1=500):
    """Order-p linear-prediction (Prony-style) estimate of the dominant pole
    modulus from the tail of a real IRF. Returns (max modulus, all moduli desc).
    """
    x = np.asarray(x, dtype=float)
    seg = x[t0:t1]
    scale = np.max(np.abs(seg))
    if scale < 1e-14:
        return 0.0, np.array([])
    seg = seg / scale
    n = len(seg)
    if n <= p + 5:
        return float('nan'), np.array([])
    X = np.column_stack([seg[p - 1 - k: n - 1 - k] for k in range(p)])
    y = seg[p:n]
    a, *_ = np.linalg.lstsq(X, y, rcond=None)
    C = np.zeros((p, p))
    C[0, :] = a
    for i in range(1, p):
        C[i, i - 1] = 1.0
    eigvals = np.linalg.eigvals(C)
    moduli = np.sort(np.abs(eigvals))[::-1]
    return float(moduli[0]), moduli


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
        p8_d, moduli8_d = dominant_pole_modulus(bgov_d, p=8)
        p8_z, moduli8_z = dominant_pole_modulus(bgov_z, p=8)
        p4_d, moduli4_d = dominant_pole_modulus(bgov_d, p=4)
        p4_z, moduli4_z = dominant_pole_modulus(bgov_z, p=4)

        m = {
            'energy_mod_def': e_mod_d, 'energy_mod_Z': e_mod_z,
            'prony8_mod_def': p8_d, 'prony8_mod_Z': p8_z,
            'prony4_mod_def': p4_d, 'prony4_mod_Z': p4_z,
            'prony8_top3_def': [float(v) for v in moduli8_d[:3]],
            'prony8_top3_Z':   [float(v) for v in moduli8_z[:3]],
            'bgov_D_499_def': float(bgov_d[499]), 'bgov_D_499_Z': float(bgov_z[499]),
            'n_inter_D_0_def': float(irf_d['n_inter_D'][0]),
        }
        RESULTS[key] = m
        v_energy = 'STABLE' if max(e_mod_d, e_mod_z) < 0.999 else 'EXPLOSIVE'
        v_prony8 = 'STABLE' if max(p8_d, p8_z) < 0.999 else 'EXPLOSIVE'
        v_prony4 = 'STABLE' if max(p4_d, p4_z) < 0.999 else 'EXPLOSIVE'
        print(f"  energy: {v_energy} (mod={max(e_mod_d,e_mod_z):.4f})   "
              f"prony(p=8): {v_prony8} (mod={max(p8_d,p8_z):.4f})   "
              f"prony(p=4): {v_prony4} (mod={max(p4_d,p4_z):.4f})", flush=True)
        print(f"  prony8 top-3 moduli (def shock): {m['prony8_top3_def']}", flush=True)
    except Exception as e:
        RESULTS[key] = {'error': repr(e)}
        print(f"  FAILED: {e!r}", flush=True)

out = ROOT / 'audit_artifacts' / 'philamb_sweep_postC1_fine_results.json'
with open(out, 'w') as fh:
    json.dump(RESULTS, fh, indent=2, default=str)

print("\n" + "=" * 110)
print("FINE phi_lamb GRID (mv_rule=1, empirical duration): energy-ratio proxy vs. Prony eigenvalue estimate")
print("=" * 110)
print(f"{'phi_lamb':>9} {'energy_v':>10} {'e_mod':>8} {'prony8_v':>10} {'p8_mod':>8} "
      f"{'prony4_v':>10} {'p4_mod':>8} {'bgov[499]':>11}")
for phi in PHI_GRID:
    m = RESULTS.get(f"phi{phi}", {})
    if 'error' in m:
        print(f"{phi:>9} ERROR: {m['error'][:80]}")
        continue
    e_mod = max(m['energy_mod_def'], m['energy_mod_Z'])
    p8_mod = max(m['prony8_mod_def'], m['prony8_mod_Z'])
    p4_mod = max(m['prony4_mod_def'], m['prony4_mod_Z'])
    v_e  = 'STABLE' if e_mod  < 0.999 else 'EXPLOSIVE'
    v_p8 = 'STABLE' if p8_mod < 0.999 else 'EXPLOSIVE'
    v_p4 = 'STABLE' if p4_mod < 0.999 else 'EXPLOSIVE'
    print(f"{phi:>9} {v_e:>10} {e_mod:>8.4f} {v_p8:>10} {p8_mod:>8.4f} "
          f"{v_p4:>10} {p4_mod:>8.4f} {m['bgov_D_499_def']:>11.4g}")
print("=" * 110)
print(f"Results -> {out}")
print("FINE SWEEP COMPLETE")
