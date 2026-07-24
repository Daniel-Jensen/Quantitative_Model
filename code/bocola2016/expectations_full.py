# FULL-MODEL EXPECTATIONS: THREE SHOCKS (eps_z, eps_g, eps_s) x TWO DEFAULT
# BRANCHES d' in {0,1}. This is the genuine two-branch quadrature of the
# banker's Euler -- the "distribution, not a point mass" structure. The
# priced default risk enters through the d'=1 branch, where the bond payoff
# is haircut (surv=1-D): a rise in the current risk factor s raises the
# weight p^d(s)=logistic(s) on that low-payoff branch, so the bond price
# falls and (through bank net worth) the constraint tightens.
# Rules are per-regime: coef[rule] = (coef_d0, coef_d1); next-period rules are
# evaluated in whichever regime d' realizes. Re-default is allowed (d'
# independent of current d), so the two current regimes share this continuation.
import numpy as np

from period_map import statics


def gauss_hermite_3d(n=3):
    # TENSOR GAUSS-HERMITE NODES/WEIGHTS FOR THREE INDEP. STANDARD NORMALS.
    x, w = np.polynomial.hermite_e.hermegauss(n)
    w = w / w.sum()
    xz, xg, xs = np.meshgrid(x, x, x, indexing="ij")
    wz, wg, ws = np.meshgrid(w, w, w, indexing="ij")
    return (xz.ravel(), xg.ravel(), xs.ravel(), (wz * wg * ws).ravel())


def default_prob(s):
    # PROB OF SOVEREIGN DEFAULT NEXT PERIOD (Bocola eq. 11): LOGISTIC IN s.
    return 1.0 / (1.0 + np.exp(-s))


def next_exog_full(Dz, g, s, ez, eg, es, cal):
    # AR(1) TRANSITIONS FOR THE THREE EXOGENOUS STATES AT THE QUADRATURE NODES.
    Dz_n = (1 - cal["rho_z"]) * cal["gamma"] + cal["rho_z"] * Dz + cal["sig_z"] * ez
    logg_n = (1 - cal["rho_g"]) * np.log(cal["g_star"]) + cal["rho_g"] * np.log(g) \
        + cal["sig_g"] * eg
    s_n = (1 - cal["rho_s"]) * cal["s_star"] + cal["rho_s"] * s + cal["sig_s"] * es
    return Dz_n, np.exp(logg_n), s_n


def _branch_pieces(S_next, d_next, grid, coef, cal):
    # PER-BRANCH INTEGRAND FACTORS AT THE NEXT STATES FOR REGIME d_next.
    S_eval = grid.clip(S_next)
    Cn = np.maximum(grid.eval(coef["C"][d_next], S_eval), 1e-8)
    avn = grid.eval(coef["av"][d_next], S_eval)
    QBn = np.maximum(grid.eval(coef["QB"][d_next], S_eval), 1e-4)
    st = statics(S_eval[:, 0], S_eval[:, 1], S_eval[:, 2], S_next[:, 3],
                 S_next[:, 4], Cn, np.ones(len(Cn)), QBn, cal, d=float(d_next))
    Qk = np.where(st["I"] > 0.0, st["Q_K"], 1.0)
    wprime = (1.0 - cal["psi"]) + cal["psi"] * avn
    surv = 1.0 - d_next * cal["haircut"]
    RKnum = (1.0 - cal["delta"]) * Qk + st["Zrent"]
    RBpay = surv * (cal["pi_mat"] + (1.0 - cal["pi_mat"]) * (cal["iota"] + QBn))
    base = wprime / Cn
    return 1.0 / Cn, base, base * RKnum, base * RBpay


def expect_pieces_full(Kp, Bp, Pp, Dz, g, s, grid, coef, cal, quad):
    # FOUR INTEGRALS AT ONE POINT OVER 3 SHOCKS x 2 DEFAULT BRANCHES.
    # q_inv=E[1/C'], q_lh=E[w'/C'], q_lh_RKnum=E[w'/C'((1-d)Q_K'+Z')],
    # q_lh_RBnum=E[w'/C' * surv(pi+(1-pi)(iota+Q_B'))]  (haircut in d'=1).
    ez, eg, es, wq = quad
    n = ez.size
    Dz_n, g_n, s_n = next_exog_full(Dz, g, s, ez, eg, es, cal)
    S_next = np.column_stack([np.full(n, Kp), np.full(n, Bp), np.full(n, Pp),
                              Dz_n, g_n, s_n])
    pd = default_prob(s)                        # weight on the default branch
    inv0, lh0, rk0, rb0 = _branch_pieces(S_next, 0, grid, coef, cal)
    inv1, lh1, rk1, rb1 = _branch_pieces(S_next, 1, grid, coef, cal)
    w0, w1 = (1.0 - pd) * wq, pd * wq
    return (np.dot(w0, inv0) + np.dot(w1, inv1),
            np.dot(w0, lh0) + np.dot(w1, lh1),
            np.dot(w0, rk0) + np.dot(w1, rk1),
            np.dot(w0, rb0) + np.dot(w1, rb1))
