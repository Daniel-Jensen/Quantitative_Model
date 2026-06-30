"""Verify economic sensibility of the market-value-rule stable cells (phi=0.10).

Checks impact (t=0) and early-path signs of the default-shock IRF:
  spread should WIDEN (+), bank net worth n_inter and output Y should FALL (-),
  debt b_gov should RISE then be stabilised. A 'stable' model with perverse
  impact signs would be a false victory.
"""
import sys, copy
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))

from sequence_jacobian import simple
import full_model
from calibration import get_calibration
from steady_state import solve_steady_state
from ic_delta_calibration import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation


@simple
def tax_rule_mv_D(q_b_D, b_gov_D, mv_gov_ss_D, phi_lamb_D):
    T_ls_D = phi_lamb_D * (q_b_D * b_gov_D(-1) - mv_gov_ss_D)
    return T_ls_D


@simple
def tax_rule_mv_F(q_b_F, b_gov_F, mv_gov_ss_F, phi_lamb_F):
    T_ls_F = phi_lamb_F * (q_b_F * b_gov_F(-1) - mv_gov_ss_F)
    return T_ls_F


cali = get_calibration()
cali['delta_b_D'] = 0.036; cali['delta_b_F'] = 0.038
cali['f_D'] = 0.03; cali['f_F'] = 0.03
cali['q_b_D'] = 0.93; cali['q_b_F'] = 0.93
ss = solve_steady_state(cali); ss = calibrate_ic_delta(ss); ss = calibrate_depreciation(ss)
ss_final = ss['ss_final']
mv_ss_D = float(ss_final['q_b_D']) * float(ss_final['b_gov_ss_D'])
mv_ss_F = float(ss_final['q_b_F']) * float(ss_final['b_gov_ss_F'])
ss_final.toplevel['mv_gov_ss_D'] = mv_ss_D   # inject before build_and_solve's baseline solve
ss_final.toplevel['mv_gov_ss_F'] = mv_ss_F

full_model.tax_rule_D = tax_rule_mv_D
full_model.tax_rule_F = tax_rule_mv_F
base = full_model.build_and_solve(ss)
ha_full = base['ha_full']; ss_final = base['ss_final']
unknowns_tp = base['unknowns_tp']; targets_tp = base['targets_tp']; T = base['T']
exo = ['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F']
dZ = 0.01 * 0.8 ** np.arange(T); zero = np.zeros(T)

CFG = {
    'A_balance_sheet': dict(psi_lambda_B=0.0, writeoff=1.0, recovery=0.40, zeta=1.0, def_scale=0.05),
    'B_risk_premium':  dict(psi_lambda_B=1.0, writeoff=0.0, recovery=0.40, zeta=1.0, def_scale=0.10),
}
PHI = 0.10
VARS = ['def_rate_D', 'spread_rb', 'q_b_D', 'n_inter_D', 'Y_D', 'b_gov_D', 'T_ls_D', 'U_D']

for name, cfg in CFG.items():
    st = copy.deepcopy(ss_final)
    for c in ['D', 'F']:
        st.toplevel[f'psi_lambda_B_{c}'] = cfg['psi_lambda_B']
        st.toplevel[f'writeoff_enabled_{c}'] = cfg['writeoff']
        st.toplevel[f'recovery_rate_{c}'] = cfg['recovery']
        st.toplevel[f'zeta_writeoff_{c}'] = cfg['zeta']
        st.toplevel[f'def_scale_{c}'] = cfg['def_scale']
        st.toplevel[f'phi_lamb_{c}'] = PHI
    st.toplevel['mv_gov_ss_D'] = mv_ss_D; st.toplevel['mv_gov_ss_F'] = mv_ss_F
    G = ha_full.solve_jacobian(st, unknowns=unknowns_tp, targets=targets_tp, inputs=exo, T=T)
    irf = G @ {'Z_D': zero, 'Z_F': zero, 'shock_def_D': dZ, 'shock_def_F': zero}
    print(f"\n=== {name}  phi={PHI}  (MARKET-VALUE rule) — default-shock IRF ===")
    print(f"{'var':>11} " + ' '.join(f't={t:>3}' and f'{t:>9}' for t in [0, 1, 2, 4, 8, 20, 60]))
    for v in VARS:
        x = np.asarray(irf[v])
        print(f"{v:>11} " + ' '.join(f'{x[t]:>+9.2e}' for t in [0, 1, 2, 4, 8, 20, 60]))
print("\nVERIFY COMPLETE")
