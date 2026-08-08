# FULL-MODEL TIME ITERATION: 6 STATES [K,B,P,Dz,g,s], TWO DEFAULT REGIMES.
# Rules {C,R,QB,av} carry a coefficient vector per regime d in {0,1}:
# coef[rule] = (coef_d0, coef_d1). The current regime enters only the current
# statics (a haircut on the bond held into a default period); the two-branch
# continuation is shared (re-default allowed). pi^d(s)=0 (cal["no_default"])
# nests the restricted model exactly -- the Stage-4 gate.
import numpy as np
from scipy.optimize import least_squares

from period_map import statics, mu_closed_form
from expectations_full import expect_pieces_full, gauss_hermite_3d, default_prob

# In the default regime (d=1) the 55% bond haircut crushes net worth and the
# constraint genuinely binds hard (leverage ~6 => mu ~ 0.15-0.4), so the
# acceptance threshold must admit real high-mu default-state solutions; only
# near-insolvent mu -> 1 is non-physical. (Restricted model used 0.10 because
# default never happens there.)
MU_ACCEPT = 0.85


def point_block_full(x, state, d_cur, grid, coef, cal, quad):
    # UNIT-FREE RESIDUALS + mu + av AT ONE POINT IN CURRENT REGIME d_cur.
    C, R, Q_B = x
    K, B, P, Dz, g, s = state
    st = statics(K, B, P, Dz, g, C, R, Q_B, cal, d=float(d_cur))
    if cal.get("no_default"):
        s_eff = -50.0                          # forces pi^d ~ 0 in the quadrature
    else:
        s_eff = s
    q_inv, q_lh, q_lh_RKnum, q_lh_RBnum = expect_pieces_full(
        st["Kp"], st["Bp"], st["Pp"], Dz, g, s_eff, grid, coef, cal, quad)

    disc = cal["beta"] * C * np.exp(-Dz)
    E_Lam = disc * q_inv
    E_Lamhat = disc * q_lh
    E_Lamhat_RK = disc * q_lh_RKnum / st["Q_K"]
    E_Lamhat_RB = disc * q_lh_RBnum / Q_B
    mu = mu_closed_form(E_Lamhat, R, st["N"], st["assets"], cal)
    required = E_Lamhat * R + cal["lam"] * mu
    res = np.array([R * E_Lam - 1.0,
                    E_Lamhat_RK / required - 1.0,
                    E_Lamhat_RB / required - 1.0])
    av = min(E_Lamhat * R / (1.0 - mu), 8.0 * cal["av_bg"])
    return res, mu, av, st


def _feasible(C, R, Q_B, state, d_cur, cal):
    # DOMAIN GUARD FOR REGIME d_cur.
    if not (C > 0 and R > 0 and Q_B > 0):
        return False
    st = statics(state[0], state[1], state[2], state[3], state[4],
                 C, R, Q_B, cal, d=float(d_cur))
    return st["I"] > 0.0 and st["N"] > 0.0 and st["assets"] > 0.0


def solve_point_full(x0, state, d_cur, grid, coef, cal, quad):
    # BOUNDED LEAST-SQUARES SOLVE AT ONE POINT IN REGIME d_cur.
    K, B, P, Dz, g, s = state
    Ymax = statics(K, B, P, Dz, g, 1e-6, x0[1], x0[2], cal, d=float(d_cur))["Y"]
    lb = np.array([1e-6, 0.90, 0.20]); ub = np.array([0.999 * Ymax, 1.10, 2.50])
    x0c = np.clip(x0, lb + 1e-9, ub - 1e-9)

    def resid(x):
        if not _feasible(x[0], x[1], x[2], state, d_cur, cal):
            return np.array([1e2, 1e2, 1e2])
        r, _, _, _ = point_block_full(x, state, d_cur, grid, coef, cal, quad)
        return r if np.all(np.isfinite(r)) else np.array([1e2, 1e2, 1e2])

    sol = least_squares(resid, x0c, bounds=(lb, ub), xtol=1e-13, ftol=1e-13,
                        gtol=1e-13, method="trf")
    ok = sol.success and np.max(np.abs(sol.fun)) < 1e-6 \
        and _feasible(sol.x[0], sol.x[1], sol.x[2], state, d_cur, cal)
    if ok:
        _, mu, _, _ = point_block_full(sol.x, state, d_cur, grid, coef, cal, quad)
        ok = mu < MU_ACCEPT
    return sol.x, ok


def initial_coef_full(grid, cal, restricted_pert):
    # SEED BOTH REGIMES FROM THE RESTRICTED PERTURBATION (FLAT IN s).
    from perturbation import policy
    u = policy(restricted_pert, grid.points[:, :5])    # ignore the s-column
    base = dict(C=u[:, 0], R=u[:, 1], QB=u[:, 2], av=u[:, 3])
    rv = {k: (v.copy(), v.copy()) for k, v in base.items()}     # (d0, d1)
    coef = {k: (grid.fit(v[0]), grid.fit(v[1])) for k, v in rv.items()}
    return coef, rv


def initial_coef_full_from_saved(grid, coef_saved):
    # WARM START FROM SAVED COEFFICIENTS: RE-DERIVE POINT VALUES AT THE GRID.
    rv = {k: (grid.eval(coef_saved[k][0], grid.points),
              grid.eval(coef_saved[k][1], grid.points))
          for k in ("C", "R", "QB", "av")}
    coef = {k: (grid.fit(rv[k][0]), grid.fit(rv[k][1])) for k in rv}
    rv = {k: (rv[k][0].copy(), rv[k][1].copy()) for k in rv}
    return coef, rv


def time_iterate_full(grid, cal, restricted_pert=None, quad=None, damp=0.3,
                      tol=1e-6, maxit=300, verbose=True, init=None):
    # TWO-REGIME OUTER TIME ITERATION.
    if quad is None:
        quad = gauss_hermite_3d(3)
    if init is not None:
        coef, rv = init
    else:
        coef, rv = initial_coef_full(grid, cal, restricted_pert)
    n = grid.n
    conv = np.inf
    # keep the BEST iterate: the rotated two-regime iteration is not globally
    # contractive (near-unit-root capital + regime coupling), so it can reach a
    # good near-fixed-point then oscillate away -- return the best, not the last
    best = dict(score=np.inf)
    for it in range(1, maxit + 1):
        new = {k: (rv[k][0].copy(), rv[k][1].copy()) for k in rv}
        mu_r = (np.zeros(n), np.zeros(n))
        n_bad = 0
        for p in range(n):
            for d_cur in (0, 1):
                x0 = np.array([rv["C"][d_cur][p], rv["R"][d_cur][p],
                               rv["QB"][d_cur][p]])
                x, ok = solve_point_full(x0, grid.points[p], d_cur, grid,
                                         coef, cal, quad)
                blk = point_block_full(x if ok else x0, grid.points[p], d_cur,
                                       grid, coef, cal, quad)
                mu_r[d_cur][p] = blk[1]
                if not ok:
                    n_bad += 1
                    continue
                new["C"][d_cur][p], new["R"][d_cur][p], new["QB"][d_cur][p] = x
                new["av"][d_cur][p] = blk[2]
        conv = max(np.max(np.abs(new[k][d] - rv[k][d]) / (np.abs(rv[k][d]) + 1e-8))
                   for k in ("C", "R", "QB", "av") for d in (0, 1))
        for k in ("C", "R", "QB", "av"):
            for d in (0, 1):
                rv[k][d][:] = (1 - damp) * rv[k][d] + damp * new[k][d]
            coef[k] = (grid.fit(rv[k][0]), grid.fit(rv[k][1]))
        # score favours few bad points then low update norm (after a warm-up)
        score = n_bad + min(conv, 1.0)
        if it >= 5 and score < best["score"]:
            best = dict(score=score, it=it, n_bad=n_bad, conv=conv,
                        coef={k: (coef[k][0].copy(), coef[k][1].copy()) for k in coef},
                        rv={k: (rv[k][0].copy(), rv[k][1].copy()) for k in rv},
                        mu=mu_r[0].copy(), mu_d1=mu_r[1].copy())
        if verbose and (it <= 3 or it % 25 == 0 or conv < tol):
            print(f"  iter {it:4d}  conv={conv:.2e}  bad={n_bad:3d}  "
                  f"mu_d0[max]={mu_r[0].max():.4f}  best@{best.get('it',0)}"
                  f"(bad={best.get('n_bad','-')})", flush=True)
        if conv < tol and it >= 3:
            break
    if best["score"] < np.inf:                  # return the best iterate
        return dict(coef=best["coef"], rv=best["rv"], quad=quad,
                    iters=best["it"], conv=best["conv"], mu=best["mu"],
                    mu_d1=best["mu_d1"], grid=grid)
    return dict(coef=coef, rv=rv, quad=quad, iters=it, conv=conv,
                mu=mu_r[0], mu_d1=mu_r[1], grid=grid)
