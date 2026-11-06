# THE ACTUARIALLY FAIR D-F SPREAD AGAINST THE TPI CAP: RISK-NEUTRAL PRICING OF THE HM PERPETUITY ON THE s PROCESS ALONE.
# no rules are read: the price is the fixed point Q(s) = (1 - p(s) h) E[delta + (1 - delta) Q(s')] / (1 + r), on a fine s grid
import numpy as np

rho, sig, pd_rest, top = 0.9, 0.41, 0.001, 0.062
dD, dF, rec, r = 0.056, 0.039, 0.45, 0.003
h = 1.0 - rec
s_star = np.log(pd_rest / (1 - pd_rest))
s_lo, s_hi = np.log(3.2e-5 / (1 - 3.2e-5)), np.log(top / (1 - top))
sg = np.linspace(s_star - 6 * sig / np.sqrt(1 - rho ** 2), np.log(0.3 / 0.7), 4001)
p = 1 / (1 + np.exp(-sg))
x, w = np.polynomial.hermite_e.hermegauss(31); w = w / w.sum()


def solve(Qfloor=None, box=None):
    # FIXED POINT OF THE RISK-NEUTRAL PRICE; Qfloor: TOMORROW'S NO-DEFAULT PRICE IS max(Q', Qfloor) (THE TPI PUT)
    Q = np.full_like(sg, dD / (dD + r))
    for _ in range(5000):
        sn = (1 - rho) * s_star + rho * sg[:, None] + sig * x[None, :]
        if box is not None:
            sn = np.clip(sn, *box)
        Qn = np.interp(sn, sg, Q)
        Qnd = Qn if Qfloor is None else np.maximum(Qn, Qfloor)
        new = ((1 - p) * ((dD + (1 - dD) * Qnd) @ w) + p * rec * ((dD + (1 - dD) * Qn) @ w)) / (1 + r)
        if np.max(np.abs(new - Q)) < 1e-13:
            break
        Q = new
    return new


yF = r  # the F bond is safe and priced at the deposit rate
spread = lambda Q: 4e4 * (dD * (1 - Q) / Q - yF)
Qrn = solve()
Qrn_box = solve(box=(s_lo, s_hi))
cap = 200.0
Qfl = dD / (dD + yF + cap / 4e4)
print(f"cap {cap:.0f} bp -> price floor Q_floor = {Qfl:.4f}")
print(f"{'p^d/qtr':>8} {'fair spread (bp)':>17} {'(s clipped to box)':>19} {'1-qtr E[R]-r at the floor, bp/yr':>34}")
for pdq in (0.0003, 0.0015, 0.006, 0.009, 0.0125, 0.021, 0.03, 0.047, 0.062, 0.09, 0.12):
    i = np.argmin(np.abs(p - pdq))
    # a risk-neutral holder buying at the floor today, the floor defended tomorrow in no default, the fair price after default
    sn = (1 - rho) * s_star + rho * sg[i] + sig * x
    Qn = np.interp(sn, sg, Qrn)
    ER = ((1 - p[i]) * (dD + (1 - dD) * np.maximum(Qn, Qfl)) + p[i] * rec * (dD + (1 - dD) * Qn)) @ w / Qfl
    print(f"{100 * p[i]:7.2f}% {spread(Qrn[i]):17.0f} {spread(Qrn_box[i]):19.0f} {4e4 * (ER - 1 - r):34.0f}")
# the p^d above which a risk-neutral bank will not hold the bond at the 200 bp floor
ER_all = []
for i in range(sg.size):
    sn = (1 - rho) * s_star + rho * sg[i] + sig * x
    Qn = np.interp(sn, sg, Qrn)
    ER_all.append(((1 - p[i]) * (dD + (1 - dD) * np.maximum(Qn, Qfl)) + p[i] * rec * (dD + (1 - dD) * Qn)) @ w / Qfl - 1 - r)
ER_all = np.array(ER_all)
j = np.argmax(ER_all < 0)
print(f"\na risk-neutral bank holds no D bonds at the 200 bp floor once p^d > {100 * p[j]:.2f}%/qtr")
j2 = np.argmax(spread(Qrn) > cap)
print(f"the fair (expected-loss-only) spread passes 200 bp at p^d {100 * p[j2]:.2f}%/qtr")
