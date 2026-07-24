"""
PAC (cross-border portfolio-adjustment cost, psi_bF_D / psi_bD_F) sweep.

Purpose: quantify the financial-accelerator COMPLEX-PAIR eigenvalue — the ~25q
"ring" visible in n_inter_D / q_b_D / q_b_F in the baseline IRFs — as a function
of the PAC level. This is the mode phi_lamb does NOT govern (phi_lamb=0.60 leaves
the debt/fiscal mode well-damped; the ring lives in the bank-networth <-> bond-price
block, whose stickiness is set by the PAC).

The PAC cost term enters the portfolio FOCs as `- psi * (b - b_ss)`, which vanishes
at the steady state, so the SS is PAC-invariant. We solve the SS pipeline ONCE and
only re-solve the baseline Jacobian per PAC value.

Committed value: psi_bF_D = psi_bD_F = 0.5.

Eigenvalue extraction: Prony / linear-prediction. Fit an AR(p) to a ring-dominated
IRF over a window, take the companion-matrix eigenvalues (the system poles), and
select the complex-conjugate pair whose period 2*pi/|arg(lambda)| falls in a
business-cycle band [10, 60]q and has the largest modulus (= the least-damped =
the visible ring). Done at p=6 and p=8; agreement across order is the trust check
(F-1's overfitting lesson). Empirical peak-spacing period is reported as an
assumption-light cross-check.
"""
import sys, copy, json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'code'))

from calibration import get_calibration
from steady_state import solve_steady_state
from ic_delta_calibration import calibrate_ic_delta
from depreciation_calibration import calibrate_depreciation
from full_model import build_and_solve

PAC_GRID = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0]   # committed = 0.5
RING_VARS = ['n_inter_D', 'q_b_F', 'q_b_D']
PERIOD_BAND = (10.0, 60.0)
WINDOW = (1, 140)   # skip t=0 impact; ~5-6 ring cycles, ring amplitude strong


def prony_poles(y, t0, t1, p):
    """AR(p) least-squares (Prony). Returns (poles, R^2)."""
    seg = np.asarray(y[t0:t1], float)
    n = len(seg)
    Y = seg[p:]
    X = np.column_stack([seg[p - k:n - k] for k in range(1, p + 1)])
    a, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ a
    ss_tot = np.sum((Y - Y.mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan
    C = np.zeros((p, p))
    C[0, :] = a
    C[1:, :-1] = np.eye(p - 1)
    return np.linalg.eigvals(C), r2


def select_ring(poles, band):
    """Complex-conjugate pole in period band with largest modulus."""
    best = None
    for lam in poles:
        if abs(lam.imag) < 1e-9:
            continue
        ang = abs(np.angle(lam))
        if ang < 1e-9:
            continue
        period = 2 * np.pi / ang
        if band[0] <= period <= band[1]:
            mod = abs(lam)
            if best is None or mod > best[0]:
                best = (mod, period, lam)
    return best


def empirical_period(y, t0, t1):
    """Peak-to-peak spacing of the (crudely detrended) signal — sanity anchor."""
    seg = np.asarray(y[t0:t1], float)
    t = np.arange(len(seg))
    detr = seg - np.polyval(np.polyfit(t, seg, 4), t)  # remove slow envelope
    # local maxima
    peaks = [i for i in range(1, len(detr) - 1)
             if detr[i] > detr[i - 1] and detr[i] > detr[i + 1]]
    if len(peaks) >= 2:
        return float(np.mean(np.diff(peaks)))
    return np.nan


def ring_estimate(irf, var, band=PERIOD_BAND, window=WINDOW):
    y = irf[var]
    out = {}
    for p in (6, 8):
        poles, r2 = prony_poles(y, window[0], window[1], p)
        sel = select_ring(poles, band)
        out[f'p{p}'] = None if sel is None else {
            'modulus': float(sel[0]), 'period': float(sel[1]), 'r2': float(r2)}
    out['emp_period'] = empirical_period(y, window[0], window[1])
    return out


def main():
    print("=" * 70)
    print("PAC sweep — financial-accelerator complex-pair eigenvalue")
    print("=" * 70)
    print("Solving SS pipeline once (PAC is SS-neutral)...")
    cali = get_calibration()
    ss = solve_steady_state(cali)
    ss = calibrate_ic_delta(ss)
    ss = calibrate_depreciation(ss)
    print("SS ready.\n")

    results = {}
    for pac in PAC_GRID:
        tag = f"psi_bF_D=psi_bD_F={pac}"
        print("\n" + "-" * 70)
        print(f"PAC = {pac}   ({tag}, committed = 0.5)")
        print("-" * 70)
        try:
            ss_mod = copy.deepcopy(ss)
            ss_mod['ss_final'].toplevel['psi_bF_D'] = pac
            ss_mod['ss_final'].toplevel['psi_bD_F'] = pac
            ss_mod['calibration_start']['psi_bF_D'] = pac
            ss_mod['calibration_start']['psi_bD_F'] = pac
            mr = build_and_solve(ss_mod)
            irf = mr['irfs_def_D']            # eigenvalues are shock-invariant
            bgov499 = float(irf['b_gov_D'][499])
            est = {v: ring_estimate(irf, v) for v in RING_VARS}
            results[str(pac)] = {'bgov_D_499': bgov499, 'rings': est}
            print(f"\n  b_gov_D[499] = {bgov499:.3e}")
            for v in RING_VARS:
                e = est[v]
                s6, s8 = e['p6'], e['p8']
                fmt = lambda s: ("  --  " if s is None
                                 else f"|λ|={s['modulus']:.4f}  T={s['period']:.1f}q  R²={s['r2']:.5f}")
                print(f"  {v:>10}:  p6[{fmt(s6)}]  p8[{fmt(s8)}]  emp_T={e['emp_period']:.1f}q")
        except Exception as ex:
            results[str(pac)] = {'error': repr(ex)}
            print(f"  FAILED: {ex!r}")

    # Summary table (use n_inter_D, p8 as headline)
    print("\n" + "=" * 70)
    print("SUMMARY — ring (complex pair) from n_inter_D, p=8")
    print("=" * 70)
    print(f"{'PAC':>6} {'|λ|':>9} {'period(q)':>10} {'R²':>9} {'half-life(q)':>12} {'b_gov[499]':>12}")
    print("-" * 70)
    for pac in PAC_GRID:
        r = results.get(str(pac), {})
        if 'error' in r:
            print(f"{pac:>6}   FAILED")
            continue
        s = r['rings']['n_inter_D']['p8']
        if s is None:
            print(f"{pac:>6}   (no ring pole in band)")
            continue
        mod = s['modulus']
        hl = np.log(0.5) / np.log(mod) if 0 < mod < 1 else float('inf')
        mark = "  <-- committed" if pac == 0.5 else ""
        print(f"{pac:>6} {mod:>9.4f} {s['period']:>10.1f} {s['r2']:>9.5f} "
              f"{hl:>12.1f} {r['bgov_D_499']:>12.2e}{mark}")
    print("-" * 70)

    out_path = Path(__file__).parent / 'pac_sweep_results.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == '__main__':
    main()
