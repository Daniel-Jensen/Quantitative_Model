"""Order-selected Prony / linear-prediction eigenvalue estimator.

The estimator docs/STATE.md refers to lived in the retired `audit_artifacts/`
harness and is not in the working tree (checked 2026-08-19: no file in the repo
matches /prony/i). Reimplemented here, self-contained, with the two properties
the original was validated on:

  * order SELECTION rather than a fixed order — textbook AR/Prony overfits and
    manufactures spurious near-unit-circle roots once the model order exceeds
    what the decay actually supports (docs/STATE.md round 2);
  * a synthetic self-test (`_selftest`) that must recover known moduli.

Used only to read the dominant modulus off an IRF tail. Never feeds a model.
"""
import numpy as np


def prony_modulus(x, orders=range(1, 9), t0=None, t1=None, tol_r2=0.999):
    """Dominant |eigenvalue| of the linear recursion that generates x[t0:t1].

    Fits x_t = sum_{k=1..p} a_k x_{t-k} by least squares for each p in `orders`,
    keeps the smallest p whose in-sample R^2 clears `tol_r2` (falling back to the
    best R^2 if none does), and returns the largest companion-matrix eigenvalue
    modulus. Returns (modulus, order, r2).
    """
    x = np.asarray(x, dtype=float)
    if t0 is None:
        t0 = int(np.argmax(np.abs(x)))          # start after the impact spike
    if t1 is None:
        t1 = len(x)
    y = x[t0:t1]
    # Drop a numerically dead tail: rows of ~0 make the design matrix singular
    # and the fit reports a meaningless modulus.
    scale = np.max(np.abs(y)) if len(y) else 0.0
    if scale == 0.0 or not np.all(np.isfinite(y)):
        return float("nan"), 0, float("nan")
    keep = np.abs(y) > 1e-14 * scale
    if keep.any():
        y = y[: int(np.max(np.nonzero(keep))) + 1]
    best = (float("nan"), 0, -np.inf)
    for p in orders:
        if len(y) < 3 * p + 5:
            continue
        X = np.column_stack([y[p - k - 1: len(y) - k - 1] for k in range(p)])
        z = y[p:]
        if X.shape[0] <= p:
            continue
        a, *_ = np.linalg.lstsq(X, z, rcond=None)
        resid = z - X @ a
        ss = float(np.sum((z - z.mean()) ** 2))
        r2 = 1.0 - float(np.sum(resid ** 2)) / ss if ss > 0 else 1.0
        C = np.zeros((p, p))
        C[0, :] = a
        if p > 1:
            C[1:, :-1] = np.eye(p - 1)
        mod = float(np.max(np.abs(np.linalg.eigvals(C))))
        if r2 > best[2]:
            best = (mod, p, r2)
        if r2 >= tol_r2:
            return mod, p, r2
    return best


def _selftest():
    t = np.arange(300)
    cases = {
        "single real 0.90": 0.90 ** t,
        "two real 0.95/0.60": 0.95 ** t + 0.5 * 0.60 ** t,
        "damped osc r=0.97": 0.97 ** t * np.cos(0.35 * t),
    }
    truth = {"single real 0.90": 0.90, "two real 0.95/0.60": 0.95,
             "damped osc r=0.97": 0.97}
    ok = True
    for k, v in cases.items():
        m, p, r2 = prony_modulus(v)
        err = abs(m - truth[k])
        ok &= err < 1e-6
        print(f"  {k:24s} -> |lam|={m:.9f} (truth {truth[k]:.2f}, "
              f"err {err:.2e}, order {p}, R2 {r2:.9f})")
    print(f"  selftest {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    _selftest()
