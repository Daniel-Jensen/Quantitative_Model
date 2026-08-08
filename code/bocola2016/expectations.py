# CONDITIONAL EXPECTATIONS FOR THE BANKER/HOUSEHOLD EULERS (RESTRICTED MODEL).
# Genuine quadrature over the two Gaussian innovations (eps_z, eps_g); the
# default branch is added in the full model. Returns the four shock-integrated
# pieces that the residuals combine with current-period Q_K, Q_B, C, R. The
# SDF growth wedge e^{-Dz_t} is DATED t (known, factored out of E_t) -- the
# classic bug is to date it t+1.
import numpy as np

from period_map import statics


def gauss_hermite_2d(n=5):
    # TENSOR GAUSS-HERMITE NODES/WEIGHTS FOR TWO INDEP. STANDARD NORMALS.
    x, w = np.polynomial.hermite_e.hermegauss(n)      # E[f], eps ~ N(0,1)
    w = w / w.sum()
    xz, xg = np.meshgrid(x, x, indexing="ij")
    wz, wg = np.meshgrid(w, w, indexing="ij")
    return xz.ravel(), xg.ravel(), (wz * wg).ravel()


def next_exog(Dz, g, eps_z, eps_g, cal):
    # AR(1) TRANSITIONS FOR THE EXOGENOUS STATES AT THE QUADRATURE NODES.
    Dz_next = (1.0 - cal["rho_z"]) * cal["gamma"] + cal["rho_z"] * Dz \
        + cal["sig_z"] * eps_z
    logg_next = (1.0 - cal["rho_g"]) * np.log(cal["g_star"]) \
        + cal["rho_g"] * np.log(g) + cal["sig_g"] * eps_g
    return Dz_next, np.exp(logg_next)


def expect_pieces(Kp, Bp, Pp, Dz, g, grid, coef, cal, quad):
    # FOUR INTEGRALS AT ONE POINT GIVEN NEXT-PERIOD RULES (RESTRICTED, d'=0).
    # q_inv       = E[1/C']
    # q_lh        = E[w'/C']                       w' = (1-psi)+psi*av'
    # q_lh_RKnum  = E[w'/C' * ((1-delta)Q_K' + Z')]
    # q_lh_RBnum  = E[w'/C' * (pi+(1-pi)(iota+Q_B'))]
    eps_z, eps_g, wq = quad
    Dz_n, g_n = next_exog(Dz, g, eps_z, eps_g, cal)
    n = eps_z.size
    S_next = np.column_stack([np.full(n, Kp), np.full(n, Bp), np.full(n, Pp),
                              Dz_n, g_n])
    # clip next states into the grid box before evaluating rules: extrapolation
    # of Chebyshev rules outside the box can produce negative C'/Q_B' and NaNs
    S_eval = grid.clip(S_next)
    Cn = grid.eval(coef["C"], S_eval)
    avn = grid.eval(coef["av"], S_eval)
    QBn = grid.eval(coef["QB"], S_eval)
    # guard the next-period real block against negative consumption/investment
    Cn = np.maximum(Cn, 1e-8)
    QBn = np.maximum(QBn, 1e-4)
    st = statics(S_eval[:, 0], S_eval[:, 1], S_eval[:, 2], Dz_n, g_n,
                 Cn, np.ones(n), QBn, cal, d=0.0)
    # SMOOTH floor on the next-period investment rate so Q_K' is differentiable
    # everywhere (a hard np.where(I'>0,...) kink makes the FD Jacobian noisy and
    # stalls Newton). Softplus floor: EXACTLY ik for feasible ik (BGP stays a
    # rest point to 1e-12), smoothly -> eps only near/below I'=0.
    ik = st["ik"]
    eps, w = 1e-6, 1e-4
    ik_safe = eps + w * np.logaddexp(0.0, (ik - eps) / w)
    Qk_n = ik_safe ** cal["xi"] / (cal["a1"] * (1.0 - cal["xi"]))
    Zrent_n = st["Zrent"]

    wprime = (1.0 - cal["psi"]) + cal["psi"] * avn
    base = wprime / Cn
    RKnum = (1.0 - cal["delta"]) * Qk_n + Zrent_n
    RBpay = cal["pi_mat"] + (1.0 - cal["pi_mat"]) * (cal["iota"] + QBn)

    q_inv = np.dot(wq, 1.0 / Cn)
    q_lh = np.dot(wq, base)
    q_lh_RKnum = np.dot(wq, base * RKnum)
    q_lh_RBnum = np.dot(wq, base * RBpay)
    return q_inv, q_lh, q_lh_RKnum, q_lh_RBnum
