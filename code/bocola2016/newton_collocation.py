# NEWTON ON THE STACKED COLLOCATION SYSTEM (replaces slow time iteration).
# Time iteration is successive approximation on the policy: its contraction rate
# is the model's persistence (~0.98 here, near-unit-root capital), so it crawls
# and oscillates and cannot resolve the barely-binding mu~0.001. Newton solves
# the collocation fixed point F(V)=0 DIRECTLY -- quadratic convergence,
# indifferent to the eigenvalue. Unknowns V = {C,R,Q_B,av} at every grid point
# (av is now a genuine unknown pinned by its recursion residual, not a read-off);
# residuals = the four Euler equations at every point. Damped Newton with a
# reusable LU-factored FD Jacobian + Broyden updates (the code/global/solvers.py
# pattern, serial, self-contained).
import numpy as np
from scipy.linalg import lu_factor, lu_solve

from period_map import statics, mu_closed_form
from expectations import expect_pieces, gauss_hermite_2d

_SQRT_EPS = np.sqrt(np.finfo(float).eps)
_WALL = 50.0


def _point_residual(C, R, QB, av, state, grid, coef, cal, quad):
    # FOUR UNIT-FREE RESIDUALS AT ONE POINT (three Eulers + the av recursion).
    # Infeasible corners get a SMOOTH bounded penalty (not a hard wall) that
    # pushes C/QB toward feasibility, so the FD Jacobian stays well-defined.
    K, B, P, Dz, g = state
    st = statics(K, B, P, Dz, g, C, R, QB, cal, d=0.0)
    Yc = max(st["Y"], 1e-6)
    if st["I"] <= 1e-6 * Yc or st["N"] <= 1e-6 or st["assets"] <= 1e-6:
        pen = 2.0 + 10.0 * (max(-st["I"], 0.0) / Yc + max(-st["N"], 0.0))
        return np.full(4, np.tanh(pen) * 5.0)      # bounded in (-5, 5), smooth-ish
    q_inv, q_lh, q_lh_RK, q_lh_RB = expect_pieces(
        st["Kp"], st["Bp"], st["Pp"], Dz, g, grid, coef, cal, quad)
    disc = cal["beta"] * C * np.exp(-Dz)
    E_Lam = disc * q_inv
    E_Lamhat = disc * q_lh
    E_Lamhat_RK = disc * q_lh_RK / st["Q_K"]
    E_Lamhat_RB = disc * q_lh_RB / QB
    mu = mu_closed_form(E_Lamhat, R, st["N"], st["assets"], cal)
    required = E_Lamhat * R + cal["lam"] * mu
    r = np.array([R * E_Lam - 1.0,
                  E_Lamhat_RK / required - 1.0,
                  E_Lamhat_RB / required - 1.0,
                  av * (1.0 - mu) / (E_Lamhat * R) - 1.0])
    return np.clip(r, -5.0, 5.0)                   # keep the FD Jacobian sane


def _fb(a, b):
    # FISCHER-BURMEISTER: phi(a,b)=a+b-sqrt(a^2+b^2)=0 iff a>=0,b>=0,ab=0.
    # Smooth everywhere except (0,0), so it turns the occasionally-binding
    # max() kink into a Newton-friendly equation.
    return a + b - np.sqrt(a * a + b * b)


def _point_residual_fb(C, R, QB, av, mu, state, grid, coef, cal, quad):
    # FIVE UNIT-FREE RESIDUALS: 3 Eulers + av recursion + FB complementarity,
    # with mu an EXPLICIT unknown (no max() kink).
    K, B, P, Dz, g = state
    st = statics(K, B, P, Dz, g, C, R, QB, cal, d=0.0)
    Yc = max(st["Y"], 1e-6)
    if st["I"] <= 1e-6 * Yc or st["N"] <= 1e-6 or st["assets"] <= 1e-6:
        return np.full(5, np.tanh(2.0) * 5.0)
    q_inv, q_lh, q_lh_RK, q_lh_RB = expect_pieces(
        st["Kp"], st["Bp"], st["Pp"], Dz, g, grid, coef, cal, quad)
    disc = cal["beta"] * C * np.exp(-Dz)
    E_Lam = disc * q_inv
    E_Lamhat = disc * q_lh
    E_Lamhat_RK = disc * q_lh_RK / st["Q_K"]
    E_Lamhat_RB = disc * q_lh_RB / QB
    required = E_Lamhat * R + cal["lam"] * mu
    # constraint slack (>=0 when leverage <= av/lam): av*N - lam*assets, scaled
    slack = (av * st["N"] - cal["lam"] * st["assets"]) / (cal["lam"] * st["assets"])
    r = np.array([R * E_Lam - 1.0,
                  E_Lamhat_RK / required - 1.0,
                  E_Lamhat_RB / required - 1.0,
                  av * (1.0 - mu) / (E_Lamhat * R) - 1.0,
                  _fb(mu, slack)])
    return np.clip(r, -5.0, 5.0)


def make_residual_fb(grid, cal, quad, mask, target):
    # STACKED FB RESIDUAL, V = [C, R, QB, av, mu] each length n.
    n = grid.n
    lo = np.concatenate([np.full(n, 1e-4), np.full(n, 0.90), np.full(n, 0.2),
                         np.full(n, 0.2), np.full(n, 0.0)])
    hi = np.concatenate([np.full(n, 5.0), np.full(n, 1.10), np.full(n, 2.5),
                         np.full(n, 5.0), np.full(n, 0.95)])

    def residual(V):
        if not np.all(np.isfinite(V)):
            return np.full(5 * n, 5.0)
        V = np.clip(V, lo, hi)
        C, R, QB = V[:n], V[n:2 * n], V[2 * n:3 * n]
        av, mu = V[3 * n:4 * n], V[4 * n:]
        coef = {"C": grid.fit(C), "R": grid.fit(R), "QB": grid.fit(QB),
                "av": grid.fit(av)}
        F = np.empty(5 * n)
        for p in range(n):
            if mask[p]:
                r5 = _point_residual_fb(C[p], R[p], QB[p], av[p], mu[p],
                                        grid.points[p], grid, coef, cal, quad)
            else:
                r5 = np.array([C[p] - target[p], R[p] - target[n + p],
                               QB[p] - target[2 * n + p], av[p] - target[3 * n + p],
                               mu[p] - target[4 * n + p]])
            for j in range(5):
                F[j * n + p] = r5[j]
        return F

    return residual


def solve_restricted_newton_fb(grid, cal, init_rv, quad=None, accept_tol=1e-9,
                               resid_thresh=0.05, verbose=True):
    # FB-SMOOTHED NEWTON POLISH: mu explicit, no kink. init_rv needs C,R,QB,av;
    # mu seeded from the closed form at the warm start.
    if quad is None:
        quad = gauss_hermite_2d(5)
    n = grid.n
    from time_iteration import point_block
    mu0 = np.array([point_block([init_rv["C"][p], init_rv["R"][p],
                    init_rv["QB"][p]], grid.points[p], grid,
                    {k: grid.fit(init_rv[k]) for k in ("C", "R", "QB", "av")},
                    cal, quad)[1] for p in range(n)])
    y0 = np.concatenate([init_rv["C"], init_rv["R"], init_rv["QB"],
                         init_rv["av"], mu0])
    lo = np.concatenate([np.full(n, 1e-4), np.full(n, 0.90), np.full(n, 0.2),
                         np.full(n, 0.2), np.full(n, 0.0)])
    hi = np.concatenate([np.full(n, 5.0), np.full(n, 1.10), np.full(n, 2.5),
                         np.full(n, 5.0), np.full(n, 0.95)])
    y0 = np.clip(y0, lo, hi)                    # anchors start at exactly zero
    r0 = np.abs(make_residual_fb(grid, cal, quad, np.ones(n, bool), y0)(y0))
    r0 = r0.reshape(5, n).max(axis=0)
    mask = feasible_mask(grid, cal, init_rv) & (r0 < resid_thresh)
    if verbose:
        print(f"    [newton-fb] polishing {mask.sum()}/{n} points "
              f"({mask.mean():.0%}); anchoring {n - mask.sum()}", flush=True)
    residual = make_residual_fb(grid, cal, quad, mask, y0)
    y, F, ok = newton_solve(residual, y0, accept_tol=accept_tol, verbose=verbose)
    coef = {"C": grid.fit(y[:n]), "R": grid.fit(y[n:2 * n]),
            "QB": grid.fit(y[2 * n:3 * n]), "av": grid.fit(y[3 * n:4 * n])}
    Ff = F.reshape(5, n)[:, mask]
    return dict(coef=coef, y=y, mu=y[4 * n:], resid=np.max(np.abs(Ff)), ok=ok,
                quad=quad, grid=grid, mask=mask)


def feasible_mask(grid, cal, rv):
    # POINTS WHOSE STATICS ARE FEASIBLE AT THE WARM-START RULES (Newton solves
    # these; infeasible never-visited corners are anchored to the warm start).
    ok = np.ones(grid.n, dtype=bool)
    for p in range(grid.n):
        s = grid.points[p]
        st = statics(s[0], s[1], s[2], s[3], s[4], rv["C"][p], rv["R"][p],
                     rv["QB"][p], cal, d=0.0)
        ok[p] = (st["I"] > 1e-6 * max(st["Y"], 1e-6) and st["N"] > 1e-6
                 and st["assets"] > 1e-6)
    return ok


def make_residual(grid, cal, quad, mask=None, target=None):
    # STACKED RESIDUAL F(V), V = [C(n), R(n), QB(n), av(n)]. At feasible points
    # (mask True) the four Euler residuals; at infeasible corners the anchor
    # residual (V - target), so the square system stays well-conditioned.
    n = grid.n
    if mask is None:
        mask = np.ones(n, dtype=bool)
    lo = np.concatenate([np.full(n, 1e-4), np.full(n, 0.90), np.full(n, 0.2),
                         np.full(n, 0.2)])
    hi = np.concatenate([np.full(n, 5.0), np.full(n, 1.10), np.full(n, 2.5),
                         np.full(n, 5.0)])

    def residual(V):
        if not np.all(np.isfinite(V)):
            return np.full(4 * n, 5.0)
        V = np.clip(V, lo, hi)
        C, R, QB, av = V[:n], V[n:2 * n], V[2 * n:3 * n], V[3 * n:]
        coef = {"C": grid.fit(C), "R": grid.fit(R), "QB": grid.fit(QB),
                "av": grid.fit(av)}
        F = np.empty(4 * n)
        for p in range(n):
            if mask[p]:
                r4 = _point_residual(C[p], R[p], QB[p], av[p], grid.points[p],
                                     grid, coef, cal, quad)
            else:                                   # anchor the corner to target
                r4 = np.array([C[p] - target[p], R[p] - target[n + p],
                               QB[p] - target[2 * n + p], av[p] - target[3 * n + p]])
            F[p] = r4[0]; F[n + p] = r4[1]; F[2 * n + p] = r4[2]; F[3 * n + p] = r4[3]
        return F

    return residual


def _fd_jacobian(residual, y0, F0):
    # SERIAL FORWARD-DIFFERENCE JACOBIAN WITH PENALTY-WALL BACK-OFF.
    n = y0.size
    steps = _SQRT_EPS * np.maximum(np.abs(y0), 1.0)
    J = np.empty((F0.size, n))
    for j in range(n):
        h = steps[j]; y = y0.copy(); y[j] += h
        Fj = residual(y)
        if np.max(np.abs(Fj - F0)) >= _WALL:
            y = y0.copy(); y[j] += h / 8.0
            Fs = residual(y)
            if np.max(np.abs(Fs - F0)) < _WALL:
                Fj, h = Fs, h / 8.0
        J[:, j] = (Fj - F0) / h
    return J


def newton_solve(residual, y0, accept_tol=1e-9, max_iter=60, max_rebuilds=6,
                 verbose=True):
    # DAMPED NEWTON ON max|F| WITH LU-FACTORED JACOBIAN + BROYDEN UPDATES.
    y = np.array(y0, float)
    F = residual(y); fnorm = np.max(np.abs(F))
    if verbose:
        print(f"    [newton] start max|F|={fnorm:.3e}", flush=True)
    if fnorm <= accept_tol:
        return y, F, True
    J = _fd_jacobian(residual, y, F); lu = lu_factor(J)
    rebuilds = 0; just_rebuilt = False
    for it in range(max_iter):
        dy = -lu_solve(lu, F)
        accepted = False
        for lam in (1.0, 0.5, 0.25, 0.1, 0.03):
            y_new = y + lam * dy; F_new = residual(y_new)
            fn = np.max(np.abs(F_new))
            if np.isfinite(fn) and fn < fnorm:
                accepted = True; break
        if accepted:
            just_rebuilt = False
            s = y_new - y; dF = F_new - F
            J += np.outer((dF - J @ s) / (s @ s), s)   # Broyden good update
            try:
                lu = lu_factor(J)
            except np.linalg.LinAlgError:
                J = _fd_jacobian(residual, y_new, F_new); lu = lu_factor(J)
            y, F, fnorm = y_new, F_new, fn
            if verbose:
                print(f"    [newton] it {it+1:2d}: max|F|={fnorm:.3e} (step {lam})",
                      flush=True)
            if fnorm <= accept_tol:
                return y, F, True
        else:
            if just_rebuilt or rebuilds >= max_rebuilds:
                return y, F, False
            rebuilds += 1
            if verbose:
                print(f"    [newton] rebuild {rebuilds}/{max_rebuilds} "
                      f"at max|F|={fnorm:.3e}", flush=True)
            J = _fd_jacobian(residual, y, F); lu = lu_factor(J)
            just_rebuilt = True
    return y, F, fnorm <= accept_tol


def solve_restricted_newton(grid, cal, init_rv, quad=None, accept_tol=1e-9,
                            resid_thresh=0.05, verbose=True):
    # NEWTON POLISH OF THE RESTRICTED COLLOCATION SYSTEM FROM A WARM START.
    # Newton solves only the points already near their Euler solution at the
    # warm start (residual < resid_thresh); the never-visited corners a warm
    # start leaves with a large residual (infeasible or high-mu) are anchored,
    # so Newton polishes the ergodic interior to machine precision.
    if quad is None:
        quad = gauss_hermite_2d(5)
    n = grid.n
    y0 = np.concatenate([init_rv["C"], init_rv["R"], init_rv["QB"], init_rv["av"]])
    # per-point warm-start residual (all points feasible => full Euler residual)
    full = make_residual(grid, cal, quad, mask=np.ones(n, bool), target=y0)
    r0 = np.abs(full(y0)).reshape(4, n).max(axis=0)
    mask = feasible_mask(grid, cal, init_rv) & (r0 < resid_thresh)
    if verbose:
        print(f"    [newton] polishing {mask.sum()}/{n} well-behaved points "
              f"({mask.mean():.0%}); anchoring {n - mask.sum()} corners", flush=True)
    residual = make_residual(grid, cal, quad, mask=mask, target=y0)
    y, F, ok = newton_solve(residual, y0, accept_tol=accept_tol, verbose=verbose)
    coef = {"C": grid.fit(y[:n]), "R": grid.fit(y[n:2 * n]),
            "QB": grid.fit(y[2 * n:3 * n]), "av": grid.fit(y[3 * n:])}
    # residual on the FEASIBLE points only (the meaningful accuracy)
    Ff = F.reshape(4, n)[:, mask]
    return dict(coef=coef, y=y, resid=np.max(np.abs(Ff)), ok=ok, quad=quad,
                grid=grid, mask=mask)
