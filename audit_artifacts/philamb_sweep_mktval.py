"""Experiment: fiscal rule on debt-level at MARKET value (mark-to-market).

Replaces the face/par-value Bohn rule
    T_ls = phi * (b_gov(-1) - b_gov_ss)                       [face value, current model]
with the mark-to-market value of the inherited debt stock at the CURRENT price
    T_ls = phi * (q_b * b_gov(-1) - mv_gov_ss)                [market value, this experiment]
where mv_gov_ss = q_b_ss * b_gov_ss.

Rationale: the long-duration instability is driven by MTM swings in q_b that the
face-value rule cannot see until par debt accumulates. A market-value rule reacts
to the current spread. Sign is a priori ambiguous (q_b falls when spreads widen,
so market-value debt FALLS -> rule loosens), so we measure it.

Same ported structural calibration as philamb_sweep_bankcal.py (delta_b=0.036/0.038,
f=0.03, back-solved Delta for a consistent IC-binding SS). Compared head-to-head
against the face-value benchmark (philamb_sweep_bankcal_results.json).

Implementation: monkeypatch full_model's tax-rule blocks (no edits to production
equations). SS is independent of the tax rule (T_ls=0 at SS), so one SS serves all.

SMOKE: SWEEP_SMOKE=1 runs one cell (config A, phi=0.10) to validate the swap.
"""
import os, sys, json, copy
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))

import sequence_jacobian as sj
from sequence_jacobian import simple
import full_model
from calibration              import get_calibration
from steady_state             import solve_steady_state
from ic_delta_calibration     import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation

SS_ONLY = os.environ.get('SWEEP_SS_ONLY') == '1'


# ── Market-value tax rules (swap into the dynamic model) ─────────────────────
@simple
def tax_rule_mv_D(q_b_D, b_gov_D, mv_gov_ss_D, phi_lamb_D):
    T_ls_D = phi_lamb_D * (q_b_D * b_gov_D(-1) - mv_gov_ss_D)
    return T_ls_D


@simple
def tax_rule_mv_F(q_b_F, b_gov_F, mv_gov_ss_F, phi_lamb_F):
    T_ls_F = phi_lamb_F * (q_b_F * b_gov_F(-1) - mv_gov_ss_F)
    return T_ls_F


# ── Shared structural port (same as bank-cal sweep) ──────────────────────────
cali = get_calibration()
cali['delta_b_D'] = 0.036;  cali['delta_b_F'] = 0.038
cali['f_D']       = 0.03;   cali['f_F']       = 0.03
cali['q_b_D']     = 0.93;   cali['q_b_F']     = 0.93

print("=== SS solve (ported structural cal, face-value rule; SS is rule-independent) ===", flush=True)
ss = solve_steady_state(cali)
ss = calibrate_ic_delta(ss)            # consistent IC-binding SS (Delta back-solved)
ss = calibrate_depreciation(ss)
ss_final = ss['ss_final']

# SS market value of debt -> reference level for the MV rule
mv_ss_D = float(ss_final['q_b_D']) * float(ss_final['b_gov_ss_D'])
mv_ss_F = float(ss_final['q_b_F']) * float(ss_final['b_gov_ss_F'])
ss_final.toplevel['mv_gov_ss_D'] = mv_ss_D
ss_final.toplevel['mv_gov_ss_F'] = mv_ss_F
print(f"[SS] q_b_D={float(ss_final['q_b_D']):.4f}  b_gov_ss_D={float(ss_final['b_gov_ss_D']):.4f}  "
      f"mv_gov_ss_D={mv_ss_D:.4f}")

if SS_ONLY:
    print("SWEEP_SS_ONLY=1 -> stop."); sys.exit(0)

# ── Monkeypatch the dynamic model's tax rules to market value ────────────────
full_model.tax_rule_D = tax_rule_mv_D
full_model.tax_rule_F = tax_rule_mv_F
print(">>> Patched full_model tax rules to MARKET VALUE (q_b * b_gov(-1)).", flush=True)

base        = full_model.build_and_solve(ss)     # builds model w/ patched rule; baseline G unused
ha_full     = base['ha_full']
ss_final    = base['ss_final']
ss_final.toplevel['mv_gov_ss_D'] = mv_ss_D        # ensure present post-deepcopy
ss_final.toplevel['mv_gov_ss_F'] = mv_ss_F
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
PHI_GRID = [0.03, 0.05, 0.10, 0.15, 0.25]

# Env overrides for focused re-runs (don't clobber prior outputs)
_sel = os.environ.get('SWEEP_CONFIGS')
if _sel:
    CONFIGS = {k: CONFIGS[k] for k in _sel.split(',') if k in CONFIGS}
_grid = os.environ.get('SWEEP_PHI_GRID')
if _grid:
    PHI_GRID = [float(x) for x in _grid.split(',')]
TAG = os.environ.get('SWEEP_TAG', '')

if os.environ.get('SWEEP_SMOKE') == '1':
    CONFIGS = {'A_balance_sheet': CONFIGS['A_balance_sheet']}; PHI_GRID = [0.10]


def set_params(ss_obj, cfg, phi):
    for c in ['D', 'F']:
        ss_obj.toplevel[f'psi_lambda_B_{c}']    = cfg['psi_lambda_B']
        ss_obj.toplevel[f'writeoff_enabled_{c}'] = cfg['writeoff']
        ss_obj.toplevel[f'recovery_rate_{c}']    = cfg['recovery']
        ss_obj.toplevel[f'zeta_writeoff_{c}']    = cfg['zeta']
        ss_obj.toplevel[f'def_scale_{c}']        = cfg['def_scale']
        ss_obj.toplevel[f'phi_lamb_{c}']         = phi
    ss_obj.toplevel['mv_gov_ss_D'] = mv_ss_D
    ss_obj.toplevel['mv_gov_ss_F'] = mv_ss_F


def peak(x, n=80):
    seg = np.asarray(x[:n]); return float(seg[np.argmax(np.abs(seg))])


def dom_modulus(x):
    """Per-quarter modulus of the dominant mode via energy ratio over two late
    100q windows (phase-independent; robust for oscillatory series)."""
    x = np.asarray(x); e1 = np.sum(x[300:400] ** 2); e2 = np.sum(x[400:500] ** 2)
    return float((e2 / e1) ** (1 / 200)) if e1 > 1e-30 else 0.0


PATHS = {}


def compute_metrics(G, key):
    irf_d = G @ {'Z_D': zero, 'Z_F': zero, 'shock_def_D': dZ,   'shock_def_F': zero}
    irf_z = G @ {'Z_D': dZ,   'Z_F': zero, 'shock_def_D': zero, 'shock_def_F': zero}
    bgov_d = np.asarray(irf_d['b_gov_D']); bgov_z = np.asarray(irf_z['b_gov_D'])
    PATHS[f'{key}__bgov_def'] = bgov_d.copy(); PATHS[f'{key}__bgov_Z'] = bgov_z.copy()
    md = dom_modulus(bgov_d); mz = dom_modulus(bgov_z); mm = max(md, mz)
    finite = np.all(np.isfinite(bgov_d)) and np.all(np.isfinite(bgov_z))
    n_pk = peak(irf_d['n_inter_D'])
    return {
        'stable':   bool(finite and mm < 0.999),
        'marginal': bool(finite and 0.999 <= mm <= 1.001),
        'loop_live':bool(abs(n_pk) > 1e-5),
        'mod_def':  md, 'mod_Z': mz,
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
        print(f"\n--- MV rule | {name}  phi={phi} ---", flush=True)
        try:
            ss_try = copy.deepcopy(ss_final); set_params(ss_try, cfg, phi)
            G = ha_full.solve_jacobian(ss_try, unknowns=unknowns_tp, targets=targets_tp,
                                       inputs=exogenous, T=T)
            m = compute_metrics(G, f'{name}_{phi}')
            RESULTS[name]['sweep'][f'{phi}'] = m
            v = 'STABLE' if m['stable'] else ('marginal' if m['marginal'] else 'explosive')
            print(f"  {v} live={m['loop_live']}  mod_def={m['mod_def']:.4f} mod_Z={m['mod_Z']:.4f}  "
                  f"spread_pk={m['spread_peak_def']:.2e}  n_inter_pk={m['n_inter_D_pk_def']:.2e}", flush=True)
        except Exception as e:
            RESULTS[name]['sweep'][f'{phi}'] = {'error': repr(e)}
            print(f"  FAILED: {e!r}", flush=True)

out = ROOT / 'audit_artifacts' / f'philamb_sweep_mktval{TAG}_results.json'
with open(out, 'w') as fh:
    json.dump(RESULTS, fh, indent=2, default=str)
np.savez_compressed(ROOT / 'audit_artifacts' / f'philamb_sweep_mktval{TAG}_paths.npz', **PATHS)

for name, cfg in CONFIGS.items():
    print("\n" + "=" * 96)
    print(f"MARKET-VALUE RULE | CONFIG {name}: {cfg}")
    print(f"{'phi':>5} {'verdict':>9} {'live':>5} {'mod_def':>8} {'mod_Z':>8} {'spread_pk':>11} "
          f"{'n_inter_pk':>11} {'T_ls_pk':>10} {'U_D[0]':>10}")
    print("-" * 96)
    for phi in PHI_GRID:
        m = RESULTS[name]['sweep'].get(f'{phi}', {})
        if 'error' in m:
            print(f"{phi:>5} {'ERR':>9}  {m['error'][:70]}");  continue
        v = 'STABLE' if m['stable'] else ('marginal' if m['marginal'] else 'explosive')
        print(f"{phi:>5} {v:>9} {str(m['loop_live']):>5} {m['mod_def']:>8.4f} {m['mod_Z']:>8.4f} "
              f"{m['spread_peak_def']:>11.2e} {m['n_inter_D_pk_def']:>11.2e} "
              f"{m['T_ls_D_peak_def']:>10.2e} {m['U_D_0_def']:>10.2e}")
print("=" * 96)
print(f"Results -> {out}")
print("MARKET-VALUE RULE SWEEP COMPLETE")
