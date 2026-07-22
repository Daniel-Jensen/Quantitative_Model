"""Re-test Finding F-1 (duration <-> fiscal-rule stability tension) now that C-1
is fixed (2026-07-22, see docs/eba_calibration.md / docs/STATE.md).

Finding F-1 was diagnosed BEFORE the C-1 fix: empirical long-duration bonds
(delta_b_D/F=0.036/0.038, ~7yr/6.5yr) were reported explosive at every
phi_lamb in [0.02, 0.50] under the par-value Bohn rule (mv_rule=0), with the
market-value rule (mv_rule=1) said to restore a stable plateau at phi_lamb in
[0.07, 0.12] -- but only in a risk-premium parameterization that was never the
production calibration. Separately, docs/eba_calibration.md's EBA-calibration
work found (again pre-C-1-fix) that even the SHORT current duration (0.10) was
mildly explosive under EBA concentration. The C-1 fix (multi-asset lambda_gk,
see equations_D.py/F.py) resolved the latter finding outright with no other
calibration change. This script asks whether it also changes the duration
finding: is the model still explosive at every phi_lamb with empirical long
duration, using TODAY'S actual committed calibration (EBA-anchored,
psi_lambda_B=0.31, def_scale=0.25, writeoff=0) unchanged except delta_b (the
axis under test) and phi_lamb/mv_rule (the two dimensions of the original
finding)?

Methodology follows audit_artifacts/philamb_sweep_mktval.py: fiscal-rule
parameters (phi_lamb, mv_rule) are SS-invariant (def_rate=0, T_ls=0 at SS), so
ONE steady-state solve serves the whole phi_lamb x mv_rule grid; only the
Jacobian is re-solved per cell (ha_full is built once and reused).
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
from full_model               import build_and_solve

PHI_GRID = [0.05, 0.10, 0.15, 0.25, 0.40, 0.60]
MV_RULES = [0.0, 1.0]
rho = 0.8


def dom_modulus(x):
    """Per-quarter modulus of the dominant mode via energy ratio over two late
    100q windows (phase-independent; robust for oscillatory series)."""
    x = np.asarray(x)
    e1 = np.sum(x[300:400] ** 2)
    e2 = np.sum(x[400:500] ** 2)
    return float((e2 / e1) ** (1 / 200)) if e1 > 1e-30 else 0.0


def peak(x, n=80):
    seg = np.asarray(x[:n])
    return float(seg[np.argmax(np.abs(seg))])


print("=== SS solve: today's committed calibration, delta_b -> empirical long "
      "duration (0.036/0.038) ===", flush=True)
cali = get_calibration()
cali['delta_b_D'] = 0.036
cali['delta_b_F'] = 0.038
cali['q_b_D'] = 0.93   # warm-start guess for the solver under long duration;
cali['q_b_F'] = 0.93   # q_b is solved endogenously, not a fixed calibration constant.

ss = solve_steady_state(cali)
ss = calibrate_ic_delta(ss)      # consistency check only post-C-1-fix; no overwrite
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
for mv in MV_RULES:
    for phi in PHI_GRID:
        key = f"mv{int(mv)}_phi{phi}"
        print(f"--- long_empirical duration | mv_rule={mv:.0f}  phi_lamb={phi} ---", flush=True)
        try:
            ss_try = copy.deepcopy(ss_final)
            for c in ['D', 'F']:
                ss_try.toplevel[f'phi_lamb_{c}'] = phi
                ss_try.toplevel[f'mv_rule_{c}']  = mv
            G = ha_full.solve_jacobian(ss_try, unknowns=unknowns_tp, targets=targets_tp,
                                        inputs=exogenous, T=T)
            irf_d = G @ {'Z_D': zero, 'Z_F': zero, 'shock_def_D': dZ, 'shock_def_F': zero}
            irf_z = G @ {'Z_D': dZ, 'Z_F': zero, 'shock_def_D': zero, 'shock_def_F': zero}
            bgov_d = np.asarray(irf_d['b_gov_D'])
            bgov_z = np.asarray(irf_z['b_gov_D'])
            md = dom_modulus(bgov_d)
            mz = dom_modulus(bgov_z)
            mm = max(md, mz)
            finite = bool(np.all(np.isfinite(bgov_d)) and np.all(np.isfinite(bgov_z)))
            m = {
                'stable':   bool(finite and mm < 0.999),
                'marginal': bool(finite and 0.999 <= mm <= 1.001),
                'mod_def':  md, 'mod_Z': mz,
                'bgov_D_499_def': float(bgov_d[499]), 'bgov_D_499_Z': float(bgov_z[499]),
                'n_inter_D_0_def': float(irf_d['n_inter_D'][0]),
                'Y_D_0_def':       float(irf_d['Y_D'][0]),
                'spread_peak_def': peak(irf_d['spread_rb']),
                'max_ca_res_def':  float(np.max(np.abs(irf_d['ca_res_D']))),
            }
            RESULTS[key] = m
            v = 'STABLE' if m['stable'] else ('marginal' if m['marginal'] else 'EXPLOSIVE')
            print(f"  {v}  mod_def={md:.4f} mod_Z={mz:.4f}  bgov[499]def={m['bgov_D_499_def']:.4g}  "
                  f"n_inter[0]={m['n_inter_D_0_def']:.4e}  Y_D[0]={m['Y_D_0_def']:.4e}", flush=True)
        except Exception as e:
            RESULTS[key] = {'error': repr(e)}
            print(f"  FAILED: {e!r}", flush=True)

out = ROOT / 'audit_artifacts' / 'philamb_sweep_postC1_results.json'
with open(out, 'w') as fh:
    json.dump(RESULTS, fh, indent=2, default=str)

print("\n" + "=" * 100)
print("POST-C1-FIX DURATION SWEEP: delta_b=0.036/0.038 (empirical), today's committed calibration otherwise")
print("=" * 100)
print(f"{'mv_rule':>8} {'phi_lamb':>9} {'verdict':>10} {'mod_def':>8} {'mod_Z':>8} "
      f"{'bgov[499]':>11} {'n_inter[0]':>11} {'Y_D[0]':>11}")
for mv in MV_RULES:
    for phi in PHI_GRID:
        key = f"mv{int(mv)}_phi{phi}"
        m = RESULTS.get(key, {})
        if 'error' in m:
            print(f"{mv:>8.0f} {phi:>9} {'ERROR':>10}  {m['error'][:60]}")
            continue
        v = 'STABLE' if m['stable'] else ('marginal' if m['marginal'] else 'EXPLOSIVE')
        print(f"{mv:>8.0f} {phi:>9} {v:>10} {m['mod_def']:>8.4f} {m['mod_Z']:>8.4f} "
              f"{m['bgov_D_499_def']:>11.4g} {m['n_inter_D_0_def']:>11.4e} {m['Y_D_0_def']:>11.4e}")
print("=" * 100)
print(f"Results -> {out}")
print("SWEEP COMPLETE")
