# FIRST-ORDER PERTURBATION (KLEIN 2000) OF THE RESTRICTED MODEL.
# Two uses: (i) an INDEPENDENT cross-check of the global solution's small-shock
# IRFs, and (ii) a smooth, globally-roughly-correct INITIAL GUESS for the global
# time iteration (a constant-BGP start leaves ~1/3 of collocation points
# infeasible; the linear rules seed every point near the truth).
# Near the BGP the leverage constraint binds smoothly (mu_bg>0), so the max()
# is on its positive branch and a perturbation is valid.
import numpy as np
from scipy.linalg import ordqz

from period_map import statics, mu_closed_form

# w = [K, B, P, Dz, g | C, R, QB, av] : 5 predetermined states + 4 controls.
NX = 5
NW = 9


def w_bgp(cal):
    # BGP VALUE OF THE STACKED VECTOR.
    return np.array([cal["K_bg"], cal["B_bg"], cal["P_bg"], cal["gamma"],
                     cal["g_star"], cal["C_bg"], cal["R_bg"], cal["qb_bg"],
                     cal["av_bg"]])


def _static(w, cal):
    # STATIC BLOCK FROM A STACKED w (av = w[8] not needed for the real block).
    return statics(w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7], cal, d=0.0)


def F_det(wp, w, cal):
    # DETERMINISTIC EQUILIBRIUM RESIDUAL F(w_{t+1}, w_t) (CERTAINTY EQUIVALENT).
    st = _static(w, cal)                    # current statics
    stp = _static(wp, cal)                  # next-period statics (for Q_K', Z')
    beta, psi, lam = cal["beta"], cal["psi"], cal["lam"]
    delta, pi, iota = cal["delta"], cal["pi_mat"], cal["iota"]

    Lam = beta * (w[5] / wp[5]) * np.exp(-w[3])         # SDF t->t+1
    Lamhat = Lam * ((1.0 - psi) + psi * wp[8])          # augmented SDF
    R_K = ((1.0 - delta) * stp["Q_K"] + stp["Zrent"]) / st["Q_K"]
    R_B = (pi + (1.0 - pi) * (iota + wp[7])) / w[7]
    mu = mu_closed_form(Lamhat, w[6], st["N"], st["assets"], cal)

    res = np.empty(NW)
    res[0] = wp[0] - st["Kp"]                            # capital transition
    res[1] = wp[1] - st["Bp"]                            # bond transition
    res[2] = wp[2] - st["Pp"]                            # deposit transition
    res[3] = wp[3] - ((1.0 - cal["rho_z"]) * cal["gamma"] + cal["rho_z"] * w[3])
    res[4] = np.log(wp[4]) - ((1.0 - cal["rho_g"]) * np.log(cal["g_star"])
                              + cal["rho_g"] * np.log(w[4]))
    res[5] = w[6] * Lam - 1.0                            # deposit Euler
    res[6] = w[8] * (1.0 - mu) - Lamhat * w[6]           # av recursion
    res[7] = Lamhat * R_K - Lamhat * w[6] - lam * mu     # capital Euler
    res[8] = Lamhat * R_B - Lamhat * w[6] - lam * mu     # bond Euler
    return res


def _jac(fun, x, eps=1e-6):
    # CENTRAL-DIFFERENCE JACOBIAN OF fun AT x.
    n = x.size
    f0 = fun(x)
    J = np.empty((f0.size, n))
    for j in range(n):
        dx = np.zeros(n); h = eps * max(1.0, abs(x[j])); dx[j] = h
        J[:, j] = (fun(x + dx) - fun(x - dx)) / (2.0 * h)
    return J


def klein(A, B, nk):
    # SOLVE  A E_t[w_{t+1}] = B w_t  FOR POLICY u=F k AND TRANSITION k'=P k.
    # Dynamic eigenvalues lambda solve B v = lambda A v (pencil (B,A)); stable
    # |lambda|<1 ordered top-left. With B=Q AA Z^H, A=Q BB Z^H the recursion in
    # z=Z^H w is BB E[z']=AA z, and the stable block gives F=Z21 Z11^{-1},
    # P=Z11 BB11^{-1} AA11 Z11^{-1} (Klein 2000).
    def stable(alpha, beta):
        # SELECT |lambda|=|alpha/beta|<1 (beta~0 => infinite => excluded).
        return np.abs(alpha) < np.abs(beta) * (1.0 - 1e-9)

    AA, BB, a, b, Q, Z = ordqz(B, A, sort=stable, output="complex")
    n_stable = int(np.sum(np.abs(a) < np.abs(b) * (1.0 - 1e-9)))
    if n_stable != nk:
        raise RuntimeError(f"Blanchard-Kahn fails: {n_stable} stable roots, "
                           f"{nk} predetermined (need equal).")
    Z11 = Z[:nk, :nk]; Z21 = Z[nk:, :nk]
    AA11 = AA[:nk, :nk]; BB11 = BB[:nk, :nk]
    if np.linalg.matrix_rank(Z11) < nk:
        raise RuntimeError("Z11 singular; no stable solution.")
    Z11i = np.linalg.inv(Z11)
    F = np.real(Z21 @ Z11i)                                  # policy u = F k
    P = np.real(Z11 @ np.linalg.solve(BB11, AA11) @ Z11i)    # transition k'=P k
    return F, P


def solve_perturbation(cal):
    # LINEARIZE AT THE BGP AND RETURN {w_bg, F, P, A, B} FOR THE RESTRICTED MODEL.
    w0 = w_bgp(cal)
    A = _jac(lambda wp: F_det(wp, w0, cal), w0)          # dF/dw_{t+1}
    B = -_jac(lambda w: F_det(w0, w, cal), w0)           # -dF/dw_t
    F, P = klein(A, B, NX)
    return dict(w_bg=w0, F=F, P=P, A=A, B=B, cal=cal)


def policy(sol, x):
    # LINEAR RULES (C, R, QB, av) AT STATE x=[K,B,P,Dz,g] (LEVELS).
    dx = np.atleast_2d(x) - sol["w_bg"][:NX]
    du = dx @ sol["F"].T
    return sol["w_bg"][NX:] + du


def irf(sol, shock_idx, size, T=40):
    # IMPULSE RESPONSE OF ALL 9 VARS TO A ONE-TIME STATE-INNOVATION SHOCK.
    F, P, w0 = sol["F"], sol["P"], sol["w_bg"]
    dk = np.zeros(NX); dk[shock_idx] = size
    path = np.empty((T, NW))
    for t in range(T):
        du = sol["F"] @ dk
        path[t, :NX] = dk; path[t, NX:] = du
        dk = P @ dk
    return path + w0
