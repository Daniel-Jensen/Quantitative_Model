"""phi_lamb sweep on the PORTED bank-cal calibration (lower amplification).

Path 2: reduce amplification so the model is stable at a literature-plausible
Bohn coefficient (phi_lamb ~ 0.02-0.10) WITH a live doom loop. Maps two
transmission channels (author chose "map both"):

  Config A  balance-sheet  : psi_lambda_B=0, write-off ON (recovery=0.40,
            zeta=1.0), def_scale=0.05. Default -> realized MTM losses on bank
            net worth -> Gertler-Karadi doom loop. (Resolves S-1.)
  Config B  risk-premium   : psi_lambda_B=1.0 (reduced from 3.0), write-off OFF,
            def_scale=0.10. Doom loop via divertability spread-loading only.

Shared structural port (SS-affecting bank-cal values):
  delta_b_D=0.036, delta_b_F=0.038 (7yr/6.5yr duration); f=0.03 (GK exit);
  Delta hardcoded 0.2/0.4 (skip calibrate_ic_delta -> resolves C-1).

The amplification params (psi_lambda_B, def_scale, writeoff/recovery/zeta) are
SS-INVARIANT (they vanish at SS where def_rate=0), so ONE steady state serves
both configs and every phi_lamb; only the Jacobian is re-solved per cell.

DEFERRED (not part of this experiment, noted in report): EBA bilateral exposure
re-targeting, income process (rho_z/sigma), frisch, chi1.

Run SS-only de-risk:  SWEEP_SS_ONLY=1 /opt/anaconda3/envs/ssj/bin/python audit_artifacts/philamb_sweep_bankcal.py
Run full sweep:                       /opt/anaconda3/envs/ssj/bin/python audit_artifacts/philamb_sweep_bankcal.py
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

SS_ONLY        = os.environ.get('SWEEP_SS_ONLY') == '1'
# Default: back-solve Delta for a consistent IC-binding SS (valid linearization).
# Set SWEEP_HARDCODE_DELTA=1 to instead keep Delta=0.2/0.4 (resolves C-1 but
# leaves IC slack at SS — only valid if the slack proves numerically benign).
BACKSOLVE_DELTA = os.environ.get('SWEEP_HARDCODE_DELTA') != '1'

# ── Shared structural port (affects SS) ──────────────────────────────────────
cali = get_calibration()
cali['delta_b_D'] = 0.036;  cali['delta_b_F'] = 0.038   # 7yr / 6.5yr duration
cali['f_D']       = 0.03;   cali['f_F']       = 0.03    # GK exit rate
cali['q_b_D']     = 0.93;   cali['q_b_F']     = 0.93    # better SS guess for long duration
# Delta stays hardcoded 0.2/0.4 from calibration.py (calibrate_ic_delta SKIPPED) -> C-1 resolved

mode = "Delta back-solved (consistent IC)" if BACKSOLVE_DELTA else "Delta hardcoded 0.2/0.4 (IC slack)"
print(f"=== SS solve: ported structural calibration (delta_b=0.036/0.038, f=0.03) | {mode} ===", flush=True)
ss = solve_steady_state(cali)
if BACKSOLVE_DELTA:
    ss = calibrate_ic_delta(ss)                          # consistent IC-binding SS
ss = calibrate_depreciation(ss)
ss_final = ss['ss_final']
print(f"\n[SS] theta_D={float(ss_final['theta_D']):.4f}  theta_F={float(ss_final['theta_F']):.4f}")
print(f"[SS] q_b_D={float(ss_final['q_b_D']):.6f}  rb_D={float(ss_final['rb_D']):.6f}  "
      f"delta_b_D={float(ss_final['delta_b_D']):.4f}")
_dmax = max(float(ss_final['Delta_bD_D']), float(ss_final['Delta_bF_D']),
            float(ss_final['Delta_bF_F']), float(ss_final['Delta_bD_F']))
print(f"[SS] Delta_bD_D={float(ss_final['Delta_bD_D']):.3f}  Delta_bF_D={float(ss_final['Delta_bF_D']):.3f}  "
      f"(max Delta={_dmax:.3f} => C-1 {'resolved' if _dmax <= 1 else 'PERSISTS (needs lambda_gk re-solve)'})")

if SS_ONLY:
    print("\nSWEEP_SS_ONLY=1 -> stopping before Jacobian sweep.")
    sys.exit(0)

# ── Build model once (reuses build_and_solve's construction; baseline G unused) ──
base        = build_and_solve(ss)
ha_full     = base['ha_full']
ss_final    = base['ss_final']
unknowns_tp = base['unknowns_tp']
targets_tp  = base['targets_tp']
T           = base['T']
exogenous   = ['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F']

rho  = 0.8
dZ   = 0.01 * rho ** np.arange(T)
zero = np.zeros(T)

CONFIGS = {
    'A_balance_sheet': dict(psi_lambda_B=0.0, writeoff=1.0, recovery=0.40, zeta=1.0, def_scale=0.05),
    'B_risk_premium':  dict(psi_lambda_B=1.0, writeoff=0.0, recovery=0.40, zeta=1.0, def_scale=0.10),
}
PHI_GRID = [0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.35, 0.50]

# Smoke mode: one solve (config A, phi=0.15) to verify the IC-slack SS gives a
# sane Jacobian before committing to the full sweep.
if os.environ.get('SWEEP_SMOKE') == '1':
    CONFIGS  = {'A_balance_sheet': CONFIGS['A_balance_sheet']}
    PHI_GRID = [0.15]


def set_params(ss_obj, cfg, phi):
    for c in ['D', 'F']:
        ss_obj.toplevel[f'psi_lambda_B_{c}']    = cfg['psi_lambda_B']
        ss_obj.toplevel[f'writeoff_enabled_{c}'] = cfg['writeoff']
        ss_obj.toplevel[f'recovery_rate_{c}']    = cfg['recovery']
        ss_obj.toplevel[f'zeta_writeoff_{c}']    = cfg['zeta']
        ss_obj.toplevel[f'def_scale_{c}']        = cfg['def_scale']
        ss_obj.toplevel[f'phi_lamb_{c}']         = phi


def peak(x, n=80):
    seg = x[:n]
    return float(seg[np.argmax(np.abs(seg))])


def tail_rho(x, t0=440, t1=498):
    """Dominant-mode geometric decay rate over the late horizon (robust to
    oscillatory crossings via median). rho<1 stationary, rho>1 explosive."""
    a = np.abs(np.asarray(x[t0:t1])); b = np.abs(np.asarray(x[t0 + 1:t1 + 1]))
    m = a > 1e-13
    return float(np.median(b[m] / a[m])) if np.any(m) else 0.0


PATHS = {}  # (config, phi) -> dict of downsampled b_gov paths, for offline reclassification


def compute_metrics(G, key):
    irf_d = G @ {'Z_D': zero, 'Z_F': zero, 'shock_def_D': dZ,   'shock_def_F': zero}
    irf_z = G @ {'Z_D': dZ,   'Z_F': zero, 'shock_def_D': zero, 'shock_def_F': zero}
    bgov_d = np.asarray(irf_d['b_gov_D']); bgov_z = np.asarray(irf_z['b_gov_D'])
    PATHS[f'{key}__bgov_def'] = bgov_d.copy()
    PATHS[f'{key}__bgov_Z']   = bgov_z.copy()

    rho_d = tail_rho(bgov_d); rho_z = tail_rho(bgov_z)
    finite = np.all(np.isfinite(bgov_d)) and np.all(np.isfinite(bgov_z))
    rho_max = max(rho_d, rho_z)
    # stationary if the dominant root is below 1 (small tolerance); marginal if ~1.
    stable   = bool(finite and rho_max < 0.9999)
    marginal = bool(finite and 0.9999 <= rho_max <= 1.0002)
    n_pk      = peak(irf_d['n_inter_D'])
    loop_live = abs(n_pk) > 1e-5
    return {
        'stable':          stable,
        'marginal':        marginal,
        'loop_live':       bool(loop_live),
        'rho_tail_def':    rho_d,
        'rho_tail_Z':      rho_z,
        'b_gov_D_499_def': float(bgov_d[499]),
        'b_gov_D_max_def': float(np.max(np.abs(bgov_d))),
        'spread_peak_def': peak(irf_d['spread_rb']),
        'n_inter_D_pk_def':n_pk,
        'Y_D_peak_def':    peak(irf_d['Y_D']),
        'T_ls_D_peak_def': peak(irf_d['T_ls_D']),
        'U_D_0_def':       float(irf_d['U_D'][0]),
        'max_ca_res_def':  float(np.max(np.abs(irf_d['ca_res_D']))),
    }


RESULTS = {}
for name, cfg in CONFIGS.items():
    RESULTS[name] = {'config': cfg, 'sweep': {}}
    for phi in PHI_GRID:
        print(f"\n--- {name}  phi_lamb={phi}  (psi_lambda_B={cfg['psi_lambda_B']}, "
              f"writeoff={cfg['writeoff']}, def_scale={cfg['def_scale']}) ---", flush=True)
        try:
            ss_try = copy.deepcopy(ss_final)
            set_params(ss_try, cfg, phi)
            G = ha_full.solve_jacobian(ss_try, unknowns=unknowns_tp, targets=targets_tp,
                                       inputs=exogenous, T=T)
            m = compute_metrics(G, f'{name}_{phi}')
            RESULTS[name]['sweep'][f'{phi}'] = m
            print(f"  stable={m['stable']} marginal={m['marginal']} live={m['loop_live']}  "
                  f"rho_def={m['rho_tail_def']:.4f} rho_Z={m['rho_tail_Z']:.4f}  "
                  f"spread_pk={m['spread_peak_def']:.2e}  n_inter_pk={m['n_inter_D_pk_def']:.2e}", flush=True)
        except Exception as e:
            RESULTS[name]['sweep'][f'{phi}'] = {'stable': False, 'error': repr(e)}
            print(f"  FAILED: {e!r}", flush=True)

# minimal stationary(or marginal)-AND-live phi per config
for name in CONFIGS:
    good = sorted(float(k) for k, v in RESULTS[name]['sweep'].items()
                  if (v.get('stable') or v.get('marginal')) and v.get('loop_live'))
    RESULTS[name]['min_stable_live_phi'] = good[0] if good else None

out_path = ROOT / 'audit_artifacts' / 'philamb_sweep_bankcal_results.json'
with open(out_path, 'w') as fh:
    json.dump(RESULTS, fh, indent=2, default=str)
np.savez_compressed(ROOT / 'audit_artifacts' / 'philamb_sweep_bankcal_paths.npz', **PATHS)

# ── Summary ──────────────────────────────────────────────────────────────────
for name, cfg in CONFIGS.items():
    print("\n" + "=" * 104)
    print(f"CONFIG {name}: {cfg}")
    print(f"{'phi':>5} {'verdict':>9} {'live':>5} {'rho_def':>8} {'rho_Z':>8} {'spread_pk':>11} "
          f"{'n_inter_pk':>11} {'Y_D_pk':>10} {'T_ls_pk':>10} {'U_D[0]':>10}")
    print("-" * 104)
    for phi in PHI_GRID:
        m = RESULTS[name]['sweep'].get(f'{phi}', {})
        if 'error' in m:
            print(f"{phi:>5} {'ERR':>9}  {m['error'][:74]}");  continue
        verdict = 'STABLE' if m['stable'] else ('marginal' if m['marginal'] else 'explosive')
        print(f"{phi:>5} {verdict:>9} {str(m['loop_live']):>5} {m['rho_tail_def']:>8.4f} "
              f"{m['rho_tail_Z']:>8.4f} {m['spread_peak_def']:>11.2e} {m['n_inter_D_pk_def']:>11.2e} "
              f"{m['Y_D_peak_def']:>10.2e} {m['T_ls_D_peak_def']:>10.2e} {m['U_D_0_def']:>10.2e}")
    print(f"min stationary/marginal + live phi_lamb: {RESULTS[name]['min_stable_live_phi']}")
print("=" * 104)
print(f"Results -> {out_path}")
print("BANKCAL PHILAMB SWEEP COMPLETE")
