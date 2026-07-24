# TIME ITERATION FOR THE RESTRICTED (NO-SOVEREIGN-RISK) 5-STATE MODEL.
# Rules {C, R, Q_B, av} on a Smolyak grid; av is a pure read-off from the
# marginal-value recursion. Pointwise 3x3 damped Newton for (C, R, Q_B) with a
# finite-difference Jacobian; outer damping on the fitted coefficients.
import numpy as np
from scipy.optimize import least_squares

from period_map import statics, mu_closed_form
from expectations import expect_pieces, gauss_hermite_2d


def point_block(x, state, grid, coef, cal, quad):
    # UNIT-FREE RESIDUALS (household/capital/bond Euler) + mu + av AT ONE POINT.
    # The two asset Eulers are written as RATIO residuals so all three are O(1)
    # (the raw asset Eulers are O(lam*mu) ~ 2e-4; scaling is essential for the
    # pointwise solve to be well-conditioned near the tiny estimated spread).
    C, R, Q_B = x
    K, B, P, Dz, g = state
    st = statics(K, B, P, Dz, g, C, R, Q_B, cal, d=0.0)
    q_inv, q_lh, q_lh_RKnum, q_lh_RBnum = expect_pieces(
        st["Kp"], st["Bp"], st["Pp"], Dz, g, grid, coef, cal, quad)

    disc = cal["beta"] * C * np.exp(-Dz)      # beta*C*e^{-Dz}: E[Lam']=disc*q_inv
    E_Lam = disc * q_inv
    E_Lamhat = disc * q_lh
    E_Lamhat_RK = disc * q_lh_RKnum / st["Q_K"]
    E_Lamhat_RB = disc * q_lh_RBnum / Q_B
    mu = mu_closed_form(E_Lamhat, R, st["N"], st["assets"], cal)

    required = E_Lamhat * R + cal["lam"] * mu   # RHS both asset Eulers share
    res = np.array([
        R * E_Lam - 1.0,                        # deposit Euler (already O(1))
        E_Lamhat_RK / required - 1.0,           # capital Euler (ratio form)
        E_Lamhat_RB / required - 1.0,           # bond Euler (ratio form)
    ])
    # marginal-value read-off, capped: av = E[Lambda_hat']R/(1-mu) blows up as
    # mu->1, but mu that large is non-physical here (Bocola's mu ~ 0.001), so
    # cap keeps a spurious pointwise trial from poisoning the fit
    av = min(E_Lamhat * R / (1.0 - mu), 8.0 * cal["av_bg"])
    return res, mu, av, st


# ANY mu above this is non-physical for Bocola's tiny estimated friction; a
# pointwise trial reaching it is spurious and rejected (keep the prior value).
MU_ACCEPT = 0.10


def _feasible(C, R, Q_B, state, cal):
    # DOMAIN GUARD: POSITIVE INVESTMENT, NET WORTH, ASSETS, AND mu < 1.
    if not (C > 0 and R > 0 and Q_B > 0):
        return False
    st = statics(state[0], state[1], state[2], state[3], state[4],
                 C, R, Q_B, cal, d=0.0)
    return st["I"] > 0.0 and st["N"] > 0.0 and st["assets"] > 0.0


def _residual(x, state, grid, coef, cal, quad):
    # RESIDUAL VECTOR FOR least_squares; LARGE PENALTY OFF THE FEASIBLE DOMAIN.
    C, R, Q_B = x
    if not _feasible(C, R, Q_B, state, cal):
        return np.array([1e2, 1e2, 1e2])
    res, _, _, _ = point_block(x, state, grid, coef, cal, quad)
    if not np.all(np.isfinite(res)):
        return np.array([1e2, 1e2, 1e2])
    return res


def solve_point(x0, state, grid, coef, cal, quad):
    # SOLVE THE 3-EQ SYSTEM AT ONE POINT (BOUNDED LEAST-SQUARES, ROBUST).
    # Bounds keep the trial physical; least_squares tolerates the near-flat
    # asset Eulers far better than a bare Newton at the tiny BGP spread.
    K, B, P, Dz, g = state
    Ymax = statics(K, B, P, Dz, g, 1e-6, x0[1], x0[2], cal, d=0.0)["Y"]
    lb = np.array([1e-6, 0.90, 0.20])
    ub = np.array([0.999 * Ymax, 1.10, 2.50])
    x0c = np.clip(x0, lb + 1e-9, ub - 1e-9)
    sol = least_squares(_residual, x0c, args=(state, grid, coef, cal, quad),
                        bounds=(lb, ub), xtol=1e-13, ftol=1e-13, gtol=1e-13,
                        method="trf")
    ok = sol.success and np.max(np.abs(sol.fun)) < 1e-6 \
        and _feasible(sol.x[0], sol.x[1], sol.x[2], state, cal)
    if ok:
        _, mu, _, _ = point_block(sol.x, state, grid, coef, cal, quad)
        ok = mu < MU_ACCEPT              # reject non-physical maximal binding
    return sol.x, ok


def initial_coef(grid, cal):
    # CONSTANT-BGP STARTING RULES (C, R, Q_B, av) FITTED ON THE GRID.
    n = grid.n
    vals = dict(C=np.full(n, cal["C_bg"]), R=np.full(n, cal["R_bg"]),
                QB=np.full(n, cal["qb_bg"]), av=np.full(n, cal["av_bg"]))
    return {k: grid.fit(v) for k, v in vals.items()}, \
        {k: v.copy() for k, v in vals.items()}


def initial_coef_perturbation(grid, cal, pert_sol):
    # STARTING RULES FROM THE FIRST-ORDER PERTURBATION POLICY AT THE GRID POINTS.
    # A smooth, globally-roughly-correct seed: unlike constant-BGP it keeps the
    # continuation rules right, so far fewer pointwise solves land infeasible.
    from perturbation import policy
    u = policy(pert_sol, grid.points)           # (n, 4) = C, R, QB, av
    vals = dict(C=u[:, 0], R=u[:, 1], QB=u[:, 2], av=u[:, 3])
    return {k: grid.fit(v) for k, v in vals.items()}, \
        {k: v.copy() for k, v in vals.items()}


def initial_coef_from_coef(grid, old_coef, old_grid):
    # WARM START ON A NEW (RE-BOXED) GRID FROM A PREVIOUS SOLUTION'S RULES.
    # Evaluate the old interpolant at the new grid points (clipped into the old
    # box so no extrapolation), then re-fit -- far faster than re-seeding from
    # the perturbation each adaptation round.
    pts = old_grid.clip(grid.points)
    vals = {k: old_grid.eval(old_coef[k], pts) for k in ("C", "R", "QB", "av")}
    return {k: grid.fit(v) for k, v in vals.items()}, \
        {k: v.copy() for k, v in vals.items()}


def time_iterate(grid, cal, quad=None, damp=0.5, tol=1e-7, maxit=800,
                 verbose=True, init=None):
    # OUTER TIME ITERATION TO A FIXED POINT OF THE RULE VALUES.
    # init=(coef, rv) supplies a warm start (e.g. from the perturbation).
    if quad is None:
        quad = gauss_hermite_2d(5)
    if init is None:
        coef, rv = initial_coef(grid, cal)      # coefficients + point values
    else:
        coef, rv = {k: v.copy() for k, v in init[0].items()}, \
            {k: v.copy() for k, v in init[1].items()}
    n = grid.n
    conv = np.inf
    for it in range(1, maxit + 1):
        new = dict(C=rv["C"].copy(), R=rv["R"].copy(), QB=rv["QB"].copy(),
                   av=rv["av"].copy(), mu=np.zeros(n))
        n_bad = 0
        for p in range(n):
            x0 = np.array([rv["C"][p], rv["R"][p], rv["QB"][p]])
            x, ok = solve_point(x0, grid.points[p], grid, coef, cal, quad)
            if not ok:            # keep previous value; don't poison the fit
                n_bad += 1
                _, mu, av, _ = point_block(x0, grid.points[p], grid, coef, cal, quad)
                new["mu"][p] = mu
                continue
            _, mu, av, _ = point_block(x, grid.points[p], grid, coef, cal, quad)
            new["C"][p], new["R"][p], new["QB"][p] = x
            new["av"][p] = av
            new["mu"][p] = mu
        # convergence on rule VALUES (kink makes coefficients cycle)
        conv = max(np.max(np.abs(new[k] - rv[k]) / (np.abs(rv[k]) + 1e-8))
                   for k in ("C", "R", "QB", "av"))
        for k in ("C", "R", "QB", "av"):
            rv[k] = (1.0 - damp) * rv[k] + damp * new[k]
            coef[k] = grid.fit(rv[k])
        if verbose and (it <= 5 or it % 25 == 0 or conv < tol):
            print(f"  iter {it:4d}  conv={conv:.3e}  bad={n_bad:2d}  "
                  f"mu[min,max]=[{new['mu'].min():.4f},{new['mu'].max():.4f}]")
        if conv < tol and it >= 3:
            break
    return dict(coef=coef, rv=rv, quad=quad, iters=it, conv=conv,
                mu=new["mu"])
