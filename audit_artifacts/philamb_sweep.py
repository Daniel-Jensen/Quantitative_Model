"""phi_lamb sweep: impact + stability of the Bohn fiscal-rule coefficient.

phi_lamb is the response of the lump-sum tax to the lagged debt gap
(tax_rule_D: T_ls = phi_lamb * (b_gov(-1) - b_gov_ss)). It is the *only*
channel by which a sovereign default-risk shock passes through to the
primary balance (T_ls enters household income -> aggregate TAX), and it
is the debt-stabilising root of the doom loop.

phi_lamb does NOT affect the steady state (at SS the debt gap is zero, so
T_ls=0 for any phi_lamb). So we run the calibration + SS pipeline ONCE,
then re-solve only the Jacobian (~3 min each) for each phi_lamb value.

Run:  /opt/anaconda3/envs/ssj/bin/python audit_artifacts/philamb_sweep.py
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
from full_model               import build_and_solve

# Sweep grid: brackets the literature-plausible low end (0.03-0.10),
# the claimed stability frontier (~0.12), the current calibration (0.15),
# and comfortably-stable high values (0.20, 0.30).
PHI_GRID = [0.03, 0.06, 0.10, 0.12, 0.15, 0.20, 0.30]

# ── Pipeline (once) ──────────────────────────────────────────────────────────
print("=== phi_lamb sweep: running calibration + SS pipeline once ===", flush=True)
cali    = get_calibration()
ss      = solve_steady_state(cali)
ss      = calibrate_ic_delta(ss)
ss      = calibrate_depreciation(ss)
print("=== building model + baseline Jacobian (phi_lamb=0.15) ===", flush=True)
base    = build_and_solve(ss)            # baseline solve at calibration phi_lamb=0.15

ha_full     = base['ha_full']
ss_final    = base['ss_final']
unknowns_tp = base['unknowns_tp']
targets_tp  = base['targets_tp']
T           = base['T']
exogenous   = ['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F']

rho  = 0.8
dZ   = 0.01 * rho ** np.arange(T)
zero = np.zeros(T)


def peak(x, n=80):
    """Signed peak (max abs) of x over the first n quarters."""
    seg = x[:n]
    return float(seg[np.argmax(np.abs(seg))])


def compute_metrics(G):
    irf_d = G @ {'Z_D': zero, 'Z_F': zero, 'shock_def_D': dZ,   'shock_def_F': zero}
    irf_z = G @ {'Z_D': dZ,   'Z_F': zero, 'shock_def_D': zero, 'shock_def_F': zero}

    bgov_d = irf_d['b_gov_D']
    bgov_z = irf_z['b_gov_D']
    # Stability: debt must decay back to SS at long horizon under BOTH shocks.
    bgov_end_d   = float(bgov_d[499])
    bgov_end_z   = float(bgov_z[499])
    bgov_max_d   = float(np.max(np.abs(bgov_d)))
    decayed      = (abs(bgov_end_d) < 1e-3) and (abs(bgov_end_z) < 1e-3)
    exploded     = (not np.all(np.isfinite(bgov_d))) or (bgov_max_d > 10.0)
    stable       = bool(decayed and not exploded)

    return {
        # --- stability ---
        'stable':          stable,
        'b_gov_D_499_def': bgov_end_d,
        'b_gov_D_499_Z':   bgov_end_z,
        'b_gov_D_max_def': bgov_max_d,
        # --- impact: default-risk shock (1pp, rho=0.8) ---
        'spread_0_def':    float(irf_d['spread_rb'][0]),
        'spread_peak_def': peak(irf_d['spread_rb']),
        'n_inter_D_0_def': float(irf_d['n_inter_D'][0]),
        'n_inter_D_pk_def':peak(irf_d['n_inter_D']),
        'Y_D_0_def':       float(irf_d['Y_D'][0]),
        'Y_D_peak_def':    peak(irf_d['Y_D']),
        'T_ls_D_peak_def': peak(irf_d['T_ls_D']),     # primary-balance tightening (passthrough)
        'U_D_0_def':       float(irf_d['U_D'][0]),
        # --- residual sanity ---
        'max_ca_res_def':  float(np.max(np.abs(irf_d['ca_res_D']))),
        'max_goodsF_def':  float(np.max(np.abs(irf_d['goods_mkt_F']))),
    }


RESULTS = {}
for phi in PHI_GRID:
    print(f"\n--- phi_lamb = {phi} : solving Jacobian ---", flush=True)
    try:
        if abs(phi - float(ss_final.toplevel['phi_lamb_D'])) < 1e-12:
            G = base['G']                      # reuse baseline solve
        else:
            ss_try = copy.deepcopy(ss_final)
            ss_try.toplevel['phi_lamb_D'] = phi
            ss_try.toplevel['phi_lamb_F'] = phi
            G = ha_full.solve_jacobian(
                ss_try, unknowns=unknowns_tp, targets=targets_tp,
                inputs=exogenous, T=T)
        m = compute_metrics(G)
        RESULTS[f'{phi}'] = m
        print(f"  stable={m['stable']}  b_gov_D[499]_def={m['b_gov_D_499_def']:.3e}  "
              f"spread_pk={m['spread_peak_def']:.3e}  n_inter_pk={m['n_inter_D_pk_def']:.3e}  "
              f"T_ls_pk={m['T_ls_D_peak_def']:.3e}", flush=True)
    except Exception as e:
        RESULTS[f'{phi}'] = {'stable': False, 'error': repr(e)}
        print(f"  FAILED: {e!r}", flush=True)

# Minimal stable phi on the grid
stable_phis = sorted(float(k) for k, v in RESULTS.items() if v.get('stable'))
RESULTS['minimal_stable_phi_grid'] = stable_phis[0] if stable_phis else None

out_path = ROOT / 'audit_artifacts' / 'philamb_sweep_results.json'
with open(out_path, 'w') as fh:
    json.dump(RESULTS, fh, indent=2, default=str)

# ── Summary table ────────────────────────────────────────────────────────────
print("\n" + "=" * 96)
print(f"{'phi':>5} {'stable':>7} {'bgov[499]def':>13} {'spread_pk':>11} "
      f"{'n_inter_pk':>11} {'Y_D_pk':>10} {'T_ls_pk':>10} {'U_D[0]':>10}")
print("-" * 96)
for phi in PHI_GRID:
    m = RESULTS.get(f'{phi}', {})
    if 'error' in m:
        print(f"{phi:>5} {'ERR':>7}  {m['error'][:70]}")
        continue
    print(f"{phi:>5} {str(m['stable']):>7} {m['b_gov_D_499_def']:>13.3e} "
          f"{m['spread_peak_def']:>11.3e} {m['n_inter_D_pk_def']:>11.3e} "
          f"{m['Y_D_peak_def']:>10.3e} {m['T_ls_D_peak_def']:>10.3e} {m['U_D_0_def']:>10.3e}")
print("=" * 96)
print(f"Minimal stable phi_lamb on grid: {RESULTS['minimal_stable_phi_grid']}")
print(f"Results -> {out_path}")
print("PHILAMB SWEEP COMPLETE")
