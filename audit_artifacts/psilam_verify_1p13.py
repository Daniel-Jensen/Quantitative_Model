"""Direct verification of the interpolated psi_lambda_B=1.1284 (from the smooth,
pre-breakdown segment of psilam_moment_sweep_postC1.py's grid) against the 150bp
target, on today's C-1-fixed EBA-anchored model.
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

cal = get_calibration()
ssr = solve_steady_state(cal)
ssr = calibrate_ic_delta(ssr)
ssr = calibrate_depreciation(ssr)
res = build_and_solve(ssr)
ha, ss = res['ha_full'], res['ss_final']
unk, tgt, T, dshock = res['unknowns_tp'], res['targets_tp'], res['T'], res['dShock_def_D']

psi_spread_base = float(ss['psi_spread_D'])
psilam_base     = float(ss['psi_lambda_B_D'])
n_ss, Y_ss = float(ss['n_inter_D']), float(ss['Y_D'])

for g in [1.10, 1.1284, 1.15]:
    ssg = copy.deepcopy(ss)
    ssg.toplevel['psi_lambda_B_D'] = g
    ssg.toplevel['psi_lambda_B_F'] = g
    ssg.toplevel['psi_spread_D']   = psi_spread_base * g / psilam_base
    ssg.toplevel['psi_spread_F']   = psi_spread_base * g / psilam_base
    Gg = ha.solve_jacobian(ssg, unknowns=unk, targets=tgt,
                           inputs=['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F'], T=T)
    irf = Gg @ {'Z_D': np.zeros(T), 'Z_F': np.zeros(T), 'shock_def_D': dshock, 'shock_def_F': np.zeros(T)}
    sp, nD, YD = np.asarray(irf['spread_rb'])[:100], np.asarray(irf['n_inter_D'])[:100], np.asarray(irf['Y_D'])[:100]
    spread_bp = float(sp.max()) * 4.0 * 1e4
    n_pct = float(nD.min()) / n_ss * 100.0
    Y_pct = float(YD.min()) / Y_ss * 100.0
    passthru = n_pct / (spread_bp / 100.0)
    print(f"psi_lambda_B={g:.4f}  spread={spread_bp:.2f}bp  dn={n_pct:+.3f}%SS  "
          f"dY={Y_pct:+.4f}%SS  passthru={passthru:+.3f}%/100bp  "
          f"n_inter_D[0]={float(irf['n_inter_D'][0])*100:+.4f}%  Y_D[0]={float(irf['Y_D'][0])*100:+.4f}%")
