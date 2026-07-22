"""psi_lambda_B moment sweep, redone on today's C-1-fixed, EBA-anchored model.

Direct descendant of diagnostics/psilam_moment_sweep.py (2026-07-14), which
produced the "psi_lambda_B=2.8, data-disciplined to ~150bp 2010 GR spread"
number cited in docs/FRAMING_HANDOFF.md. That sweep ran on the OLD model:
placeholder portfolio shares (phi_bD_D=0.25, not EBA's 2.39), and Delta
back-solved to a DEGENERATE 1.4545 (the C-1 bug, present and un-fixed at the
time -- visible in its own stdout: "WARNING (C-1): back-solved Delta > 1").
Neither that 2.8 nor the pre-existing round-number 3.0 was ever calibrated
against the current model: EBA portfolio concentration, the multi-asset-
consistent lambda_gk (C-1 fix), or the omega_K capital-fund split. This script
reruns the identical moment-matching procedure on the current
code/calibration.py base to find whatever psi_lambda_B now hits the same
external target (2010 GR-DE spread ~150bp on a 1pp default-probability shock).

Methodology (unchanged from the original): EL_price_D is anchored (SS-
determined, independent of psi_lambda_B -- it only enters def_rate(+1)=0 terms
at SS). psi_spread is EXACTLY linear in psi_lambda_B holding lambda_gk/Omega/
beta_inter fixed (all SS objects, unaffected by the counterfactual psi_lambda_B
used only in the dynamic Delta_eff = Delta + psi_lambda_B*def_rate(+1) term).
So one SS solve + one ha_full build serves the whole grid; only the Jacobian is
re-solved per psi_lambda_B value.
"""
import sys, copy, datetime
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))

from calibration              import get_calibration
from steady_state             import solve_steady_state
from ic_delta_calibration     import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model                import build_and_solve

GRID = [0.0, 0.31, 0.5, 1.0, 1.5, 2.0, 2.6, 2.8, 3.0, 4.0, 5.0]
ANNUAL = 4.0
BP = 1e4
TARGET_BP = 150.0  # 2010 GR-DE spread, the paper's external target (unchanged)


def ts():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


print(f"[{ts()}] pipeline: today's calibration.py (EBA-anchored, C-1-fixed) -> SS -> "
      f"ic_delta (consistency check) -> depreciation", flush=True)
cal = get_calibration()
ssr = solve_steady_state(cal)
ssr = calibrate_ic_delta(ssr)
ssr = calibrate_depreciation(ssr)

print(f"[{ts()}] build_and_solve (today's committed psi_lambda_B={cal['psi_lambda_B_D']})", flush=True)
res = build_and_solve(ssr)
ha, ss = res['ha_full'], res['ss_final']
unk, tgt, T, dshock = res['unknowns_tp'], res['targets_tp'], res['T'], res['dShock_def_D']

psi_spread_base = float(ss['psi_spread_D'])
psilam_base     = float(ss['psi_lambda_B_D'])
EL_D            = float(ss['EL_price_D'])
n_ss            = float(ss['n_inter_D'])
Y_ss            = float(ss['Y_D'])
print(f"[{ts()}] anchors: EL_price_D={EL_D:.6f} (FIXED), "
      f"psi_spread_D(base)={psi_spread_base:.6f} at psi_lambda_B={psilam_base}", flush=True)

zero = np.zeros(T)
rows = []
for g in GRID:
    ssg = copy.deepcopy(ss)
    ssg.toplevel['psi_lambda_B_D'] = g
    ssg.toplevel['psi_lambda_B_F'] = g
    ssg.toplevel['psi_spread_D']   = psi_spread_base * g / psilam_base
    ssg.toplevel['psi_spread_F']   = psi_spread_base * g / psilam_base
    Gg = ha.solve_jacobian(ssg, unknowns=unk, targets=tgt,
                           inputs=['Z_D', 'shock_def_D', 'Z_F', 'shock_def_F'], T=T)
    irf = Gg @ {'Z_D': zero, 'Z_F': zero, 'shock_def_D': dshock, 'shock_def_F': zero}
    sp = np.asarray(irf['spread_rb'])[:100]
    nD = np.asarray(irf['n_inter_D'])[:100]
    YD = np.asarray(irf['Y_D'])[:100]
    spread_bp = float(sp.max()) * ANNUAL * BP
    n_pct     = float(nD.min()) / n_ss * 100.0
    Y_pct     = float(YD.min()) / Y_ss * 100.0
    passthru  = n_pct / (spread_bp / 100.0) if abs(spread_bp) > 1e-9 else float('nan')
    rows.append(dict(psilam=g, psi_spread=psi_spread_base * g / psilam_base,
                     foc_load=EL_D + psi_spread_base * g / psilam_base,
                     spread_bp=spread_bp, n_pct=n_pct, Y_pct=Y_pct, passthru=passthru))
    print(f"[{ts()}]  psi_lambda_B={g:4.2f}  spread={spread_bp:8.2f}bp  "
          f"dn={n_pct:+.3f}%SS  dY={Y_pct:+.4f}%SS  passthru={passthru:+.3f}%/100bp", flush=True)

floor_bp = rows[0]['spread_bp']
for r in rows:
    r['amp_vs_floor'] = r['spread_bp'] / floor_bp if abs(floor_bp) > 1e-12 else float('nan')

# ---- find the psi_lambda_B that hits the 150bp target by linear interpolation ----
g_arr = np.array([r['psilam'] for r in rows])
sp_arr = np.array([r['spread_bp'] for r in rows])
order = np.argsort(g_arr)
g_arr, sp_arr = g_arr[order], sp_arr[order]
if TARGET_BP <= sp_arr.max() and TARGET_BP >= sp_arr.min():
    target_psilam = float(np.interp(TARGET_BP, sp_arr, g_arr))
else:
    target_psilam = float('nan')

print("\n" + "=" * 100)
print(f"psi_lambda_B moment sweep (today's model) -- target: 2010 GR-DE spread ~{TARGET_BP:.0f}bp")
print("=" * 100)
print(f"{'psi_lambda_B':>12} {'psi_spread':>10} {'FOC load':>9} {'spread(bp)':>11} "
      f"{'dn(%SS)':>9} {'dY(%SS)':>9} {'passthru':>9} {'amp':>7}")
for r in rows:
    print(f"{r['psilam']:>12.2f} {r['psi_spread']:>10.3f} {r['foc_load']:>9.3f} {r['spread_bp']:>11.2f} "
          f"{r['n_pct']:>+9.3f} {r['Y_pct']:>+9.4f} {r['passthru']:>+9.3f} {r['amp_vs_floor']:>6.2f}x")
print("=" * 100)
print(f"Interpolated psi_lambda_B hitting {TARGET_BP:.0f}bp target: {target_psilam:.4f}")
print("=" * 100)
