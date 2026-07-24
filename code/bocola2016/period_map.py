# DETRENDED STATIC PERIOD MAP: (STATE, RULE VALUES) -> ALL DERIVED QUANTITIES.
# This is the ONLY module that applies the growth relabeling e^{-Dz}: exactly
# three end-of-period stocks (Kp, Bp, Pp) carry the factor, and the SDF growth
# wedge (also e^{-Dz}, dated t) is applied in expectations.py. Everything here
# is a deterministic function of the current state and the current rule values
# (C, R, Q_B); the multiplier mu additionally needs E[Lambda_hat'] and is a
# separate closed-form helper. Vectorized over points (arrays) or scalar.
import numpy as np


def statics(K, B, P, Dz, g, C, R, Q_B, cal, d=0.0):
    # REAL BLOCK: LABOR, OUTPUT, INVESTMENT, Q_K, NEXT STOCKS, NET WORTH.
    alpha = cal["alpha"]; nu = cal["nu"]; xi = cal["xi"]
    delta = cal["delta"]; pi = cal["pi_mat"]; iota = cal["iota"]
    a1 = cal["a1"]; a2 = cal["a2"]; chi = cal["chi"]; Zbar = cal["Zbar"]
    psi = cal["psi"]; omega = cal["omega"]; D = cal["haircut"]
    tau_star = cal["tau_star"]; gamma_tau = cal["gamma_tau"]
    eDz = np.exp(Dz)
    surv = 1.0 - d * D                        # bond survival share (d=0 -> 1)

    A = Zbar * K ** alpha * eDz ** (1.0 - alpha)     # TFP-scaled capital block
    l = ((1.0 - alpha) * A / (chi * C)) ** (1.0 / (nu + alpha))
    Y = A * l ** (1.0 - alpha)
    w = (1.0 - alpha) * Y / l
    I = Y * (1.0 - g) - C
    ik = I / K                                # investment rate
    Q_K = 1.0 / (a1 * (1.0 - xi) * ik ** (-xi))
    Phi = a1 * ik ** (1.0 - xi) + a2
    Kp = ((1.0 - delta) + Phi) * K / eDz      # relabel to A_t units

    Zrent = alpha * Y / K                     # profit per unit capital
    tau = tau_star * eDz + gamma_tau * B
    # government budget (eq. 9), then relabel B' to A_t units
    Bp = ((1.0 - pi) * surv * B
          + ((pi + (1.0 - pi) * iota) * surv * B + g * Y - tau) / Q_B) / eDz

    # realized asset payoffs at current prices -> beginning-of-period net worth
    X = ((1.0 - delta) * Q_K + Zrent) * K + surv * (pi + (1.0 - pi) * (iota + Q_B)) * B
    N = psi * (X - P) + omega * X
    assets = Q_K * Kp + Q_B * Bp              # book carried into next period
    dep = assets - N
    Pp = R * dep                              # gross payment owed next period
    return dict(l=l, Y=Y, w=w, I=I, ik=ik, Q_K=Q_K, Phi=Phi, Kp=Kp,
                Zrent=Zrent, tau=tau, Bp=Bp, X=X, N=N, assets=assets,
                dep=dep, Pp=Pp)


def mu_closed_form(E_Lamhat, R, N, assets, cal):
    # OCCASIONALLY BINDING MULTIPLIER (BOCOLA PROP. 1, eq. 2), NO SMOOTHING.
    # Capped just below 1: mu>=1 means the banker value function does not exist
    # (guess unverified), so it flags an infeasible trial rather than blowing up.
    mu = 1.0 - E_Lamhat * R * N / (cal["lam"] * assets)
    return np.minimum(np.maximum(mu, 0.0), 0.999999)
