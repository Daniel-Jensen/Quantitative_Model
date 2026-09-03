# GLOBAL CHEBYSHEV COLLOCATION SOLVE -- BOCOLA'S OWN DESIGN (residual_model.m + parsolve.m).
# The POLICY VALUES AT THE COLLOCATION POINTS ARE THE UNKNOWNS. There is no inner root
# find: one residual evaluation fits the Chebyshev coefficients to the guess, walks the
# period map at every grid point with THAT interpolant as the continuation, and returns
# the equilibrium conditions. The whole coefficient vector then goes to a single Newton.
#
# WHY THIS REPLACES TIME ITERATION. recursive_main.time_iteration freezes the previous
# iterate as the continuation, solves 13 unknowns pointwise, and damps. Its binding mode
# is the franchise-value recursion alpha = E[Om]R/(1-mu) with Om = beta[f + (1-f)alpha'],
# whose slope is beta(1-f)R/(1-mu) ~ 0.96; at damp = 0.25 that is 0.990 per sweep, i.e.
# 235 sweeps per decade of rule change. Runs reported max|F| = 1e-14 and "rule-change tol
# not reached" at the same time -- every point cleared against a continuation that was
# still moving. A fixed point approached at rate rho carries a level error ~ change/(1-rho),
# so the 1e-6 exit test admits ~1e-4 in alpha. Newton has no contraction rate to leak
# through: it drives the residual itself to machine zero, which is what parsolve.m does.
#
# THE UNKNOWN SET IS BOCOLA'S, NOT OURS. He carries four unknowns (c, R, alp, q) whose
# residuals are LOG RATIOS of guess to Euler-implied value; alpha is an unknown with
# residual log(alp/alp_impl), not an object read off a frozen recursion. Here EVERY
# stored rule is an unknown: the 13 market-clearing/Euler unknowns keep their existing
# residuals from point_map, and the 6 objects point_map used to READ OFF the recursions
# (alpha_D/F, C_D/F, r_wc_D/F) get Bocola's identity residual log(guess/implied). That
# closes the knot-breaking freezes -- point_map's `cont.eval(...)` at a collocation node
# returns the guess exactly, because square collocation interpolates its nodes exactly --
# so the solved object is a genuine recursive equilibrium rather than a fixed point of a
# damped map.
#
# Nothing inside point_map.py changes. Residuals per (point, regime): 13 + 6 = 19.
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import newton_krylov
from scipy.optimize.nonlin import NoConvergence

from solver_recursive.decision_rules import (SOLVE, DERIVED, STORE_RULES,
                                             to_fit, from_fit)
from solver_recursive.point_map import point_residuals

# THE STACKING CONVENTION, in one place. theta is ordered rule-major, then regime,
# then grid point -- the image of Bocola's [cons; R; alp; q] per regime block.
N_RES_POINT = len(SOLVE)                 # residuals point_map itself returns
N_RES = N_RES_POINT + len(DERIVED)       # + the identity residuals for the read-offs
RES_NAMES = ("cap_D", "cap_F", "lab_D", "lab_F", "euler_D", "uip", "goods_D",
             "bondD_D", "bondD_F", "euler_F", "dep_clear", "bondF_F", "bondF_D"
             ) + tuple(f"id_{k}" for k in DERIVED)
assert len(RES_NAMES) == N_RES, (len(RES_NAMES), N_RES)
_BIG = 1e3                               # sentinel for an unevaluable point


def pack(rules, regimes=(0, 1)):
    # RULE VALUES -> THE FLAT UNKNOWN VECTOR, in the FITTED (log) transform.
    return np.concatenate([to_fit(k, rules.vals[k][d])
                           for k in STORE_RULES for d in regimes])


def unpack(theta, rules, regimes=(0, 1)):
    # THE FLAT VECTOR -> {rule: {regime: level values}}, the inverse of pack.
    n = rules.grid.n
    out, i = {}, 0
    for k in STORE_RULES:
        out[k] = {}
        for d in regimes:
            out[k][d] = from_fit(k, theta[i:i + n])
            i += n
    return out


def write_back(theta, rules, regimes=(0, 1)):
    # INSTALL A SOLUTION INTO A RuleSet (values in levels, coefficients refitted).
    vals = unpack(theta, rules, regimes)
    for k in STORE_RULES:
        for d in regimes:
            rules.set_values(k, d, vals[k][d])
    return rules


def make_residual(rules, cal, ss, sproc, regimes=(0, 1), no_default=False, n_gh=5,
                  scale=None, no_cb=False):
    # BUILD THE GLOBAL RESIDUAL F(theta) -- the image of residual_model.m.
    # `rules` is used as the working RuleSet: each call overwrites its values and
    # coefficients from theta, so the continuation is ALWAYS the current guess.
    # `scale` (optional, per-equation) row-scales the residuals; the zero set is
    # unchanged, only the conditioning of the linear solve and the norm are.
    n = rules.grid.n
    pts = rules.grid.points
    sw = np.ones(N_RES) if scale is None else np.asarray(scale, dtype=float)

    def F(theta):
        vals = unpack(theta, rules, regimes)
        # EXACT SQUARE FIT ONLY. A masked/ridge fit would stop the interpolant from
        # reproducing its own nodes, and point_map reads alpha/C/r_wc at the CURRENT
        # state through cont.eval -- exactness there is what makes those the unknowns.
        for k in STORE_RULES:
            for d in regimes:
                rules.set_values(k, d, vals[k][d])
        res = np.empty((len(regimes), n, N_RES))
        for jd, d in enumerate(regimes):
            for i in range(n):
                x = np.array([vals[k][d][i] for k in SOLVE])
                try:
                    r, out = point_residuals(pts[i], d, x, rules, cal, ss, sproc,
                                             n_gh=n_gh, no_default=no_default,
                                             no_cb=no_cb)
                except (ValueError, RuntimeError, ArithmeticError):
                    # ArithmeticError covers ZeroDivisionError/OverflowError too: a
                    # trial Newton step can drive an expectation to zero, and one
                    # unevaluable point must cost the STEP, not the whole solve.
                    res[jd, i, :] = _BIG
                    continue
                res[jd, i, :N_RES_POINT] = r
                # Bocola's identity residual log(guess/implied), in the same transform
                # the rule is fitted in (log level, or log gross rate).
                for q, k in enumerate(DERIVED):
                    res[jd, i, N_RES_POINT + q] = (to_fit(k, vals[k][d][i])
                                                   - to_fit(k, out[k]))
        res = np.where(np.isfinite(res), res, _BIG)
        return (res * sw).ravel()

    return F


def residual_table(theta, F):
    # PER-EQUATION max|R| FROM A STACKED RESIDUAL VECTOR (diagnostics).
    r = np.abs(F(theta)).reshape(-1, N_RES)
    return dict(zip(RES_NAMES, r.max(axis=0)))


def parsolve(F, x0, cc=1.0, tol=1e-20, maxcount=60, eps=1e-6, verbose=True,
             blowup=100.0, backtrack=True, label="", stall_step=1e-3, jac_every=1,
             floor_tol=1e-8):
    # DAMPED NEWTON WITH A FORWARD-DIFFERENCE JACOBIAN -- the port of parsolve.m.
    # Ryan Decker's routine as Bocola ships it: build J column by column at step eps,
    # take x <- x - cc*(J\f), stop when sum(f^2) <= tol, abort if it exceeds `blowup`.
    # TWO ADDITIONS over the original, both about stopping rather than about the step.
    # (1) A backtrack on the damping when a full step does not reduce the norm. Bocola
    # does this by hand instead -- he calls parsolve with cc = 1, 0.8 and 0.5 in
    # different stages -- so this automates his own practice. backtrack=False is his
    # code exactly.
    # (2) A STALL EXIT. His test is sum(f^2) <= 1e-20, which for a system this size is
    # ~1e-12 per residual and is below the arithmetic floor of the period map: the
    # capital and bond Eulers are O(1) expectations differenced to O(1e-4), so ~1e-10
    # is as small as max|F| goes and further Newton steps cannot improve it. Without
    # this the solve burns its whole iteration budget backtracking at step ~1e-5 on a
    # converged answer, which is exactly what it did on first run (stalled at
    # max|F| = 4.9e-10 and kept going). Reaching the floor IS convergence -- but ONLY
    # if the residual is actually at the floor: a stall at max|F| = 1e-3 is a failure
    # and is reported as one, which is what floor_tol decides. The achieved norm is
    # returned either way.
    # jac_every > 1 REUSES the factorised Jacobian for that many steps (a chord /
    # Shamanskii iteration). Bocola refreshes every step, which is jac_every = 1 and the
    # default; the option exists because a dense FD Jacobian costs m+1 residual
    # evaluations, and on the s-refined grid m is 6498 rather than 798.
    x = np.asarray(x0, dtype=float).copy()
    m = x.size
    f = F(x)
    gap = float(f @ f)
    ok = False
    lu = None
    for count in range(1, maxcount + 1):
        if gap <= tol:
            ok = True
            break
        if lu is None or (count - 1) % jac_every == 0:
            J = np.empty((f.size, m))
            for i in range(m):
                xp = x.copy()
                xp[i] += eps
                J[:, i] = (F(xp) - f) / eps
            try:
                lu = lu_factor(J)
            except (np.linalg.LinAlgError, ValueError):
                lu = None
        try:
            dx = lu_solve(lu, f) if lu is not None else np.linalg.lstsq(J, f, rcond=None)[0]
        except (np.linalg.LinAlgError, ValueError):
            dx = np.linalg.lstsq(J, f, rcond=None)[0]
        step = cc
        while True:
            xn = x - step * dx
            fn = F(xn)
            gn = float(fn @ fn)
            if (not backtrack) or gn < gap or step < 1e-5:
                break
            step *= 0.5
        improved = gn < gap
        if improved:
            x, f, gap = xn, fn, gn
        if verbose:
            print(f"    [parsolve{label} {count:2d}] sum|F|^2 = {gap:.3e}   "
                  f"max|F| = {np.max(np.abs(f)):.3e}   step = {step:.3g}"
                  f"{'' if improved else '   (no improvement)'}")
        if not np.isfinite(gap) or gap > blowup:
            return x, False, count, np.max(np.abs(f))
        if not improved or step < stall_step:
            # no admissible step reduces the norm: converged if that is because the
            # residual is already at the period map's arithmetic floor, failed if not
            worst = float(np.max(np.abs(f)))
            return x, worst <= floor_tol, count, worst
    if gap <= tol:
        ok = True
    return x, ok, count, float(np.max(np.abs(f)))


def krylov_solve(F, x0, f_tol=1e-11, maxiter=60, verbose=True, label=""):
    # JACOBIAN-FREE NEWTON-KRYLOV. Same Newton, solving the linear system by GMRES
    # instead of forming J -- attractive because a dense finite-difference Jacobian
    # costs m+1 residual evaluations, which is 0.7 min at 21 points and 48 min at 171.
    #
    # IT DOES NOT WORK ON THIS SYSTEM. Measured on the s-refined grid, from a warm start
    # at max|F| = 1.9e-3: 25 outer iterations move the norm to 1.88e-3, i.e. nowhere,
    # and the same happens from a cold seed at 2.4e-2. The residual rows span four
    # orders (the capital and bond Eulers are in return units ~1e-4; the identity
    # residuals are log ratios ~1e-2) and the unknowns are log-policies of very
    # different scales, so unpreconditioned GMRES returns a direction the Armijo search
    # then cuts to nothing. parsolve -- Bocola's dense factorisation of the same
    # Jacobian -- converges on the identical system. Kept, tested and NOT the default;
    # see BACKEND_DENSE_MAX.
    it = {"n": 0}

    def cb(x, fx):
        it["n"] += 1
        if verbose:
            print(f"    [krylov{label} {it['n']:2d}] max|F| = {np.max(np.abs(fx)):.3e}")
    try:
        x = newton_krylov(F, x0, f_tol=f_tol, maxiter=maxiter, callback=cb,
                          method="lgmres", line_search="armijo", verbose=False)
        f = F(x)
        return x, bool(np.max(np.abs(f)) <= f_tol * 10), it["n"], float(np.max(np.abs(f)))
    except NoConvergence as e:
        x = np.asarray(e.args[0], dtype=float)
        f = F(x)
        return x, False, it["n"], float(np.max(np.abs(f)))


# ACCEPTANCE ON max|F|. 1e-9 is four orders below any economic signal in this model
# (the headline risk shock moves mu by 7.5e-3 and the labour residual by ~7e-4) and it
# is reachable: the period map's own arithmetic floor is ~1e-10, because the capital and
# bond Eulers difference O(1) expectations down to O(1e-4).
TOL_MAXF = 1e-9
# ABOVE THIS MANY UNKNOWNS "auto" would switch to the Jacobian-free backend. It is set
# beyond any grid this model uses, because Newton-Krylov measured non-convergent here
# (see krylov_solve) -- the dense Newton is the default at every size. Lower it only to
# re-test the Krylov path.
BACKEND_DENSE_MAX = 10 ** 9


def solve_collocation(rules, cal, ss, sproc, regimes=(0, 1), no_default=False,
                      n_gh=5, backend="auto", tol=TOL_MAXF, maxit=60, cc=1.0,
                      verbose=True, label="", jac_every=1, no_cb=False):
    # SOLVE THE WHOLE COEFFICIENT VECTOR AT ONCE AND WRITE THE ANSWER BACK INTO `rules`.
    # backend "parsolve" is Bocola's dense-Jacobian Newton (the default whenever the
    # system is small enough to afford it); "krylov" is the Jacobian-free variant for
    # refined grids; "auto" picks on the unknown count.
    F = make_residual(rules, cal, ss, sproc, regimes=regimes,
                      no_default=no_default, n_gh=n_gh, no_cb=no_cb)
    x0 = pack(rules, regimes)
    m = x0.size
    if backend == "auto":
        backend = "parsolve" if m <= BACKEND_DENSE_MAX else "krylov"
    # jac_every STAYS AT BOCOLA'S 1. Reusing the factorisation is cheaper per step but
    # only linearly convergent, and the Jacobian dominates the cost: measured on the
    # s = 5 grid, jac_every = 2 went 1.9e-3 -> 1.8e-4 -> 1.1e-4 (the reused step barely
    # moved), where a fresh Jacobian each step converges quadratically. Fewer, better
    # steps beat more, cheaper ones here. jac_every > 1 is left available for a grid
    # where the Jacobian genuinely cannot be afforded.
    if verbose:
        print(f"  collocation solve{label}: {m} unknowns "
              f"({len(STORE_RULES)} rules x {len(regimes)} regime(s) x {rules.grid.n} "
              f"points), backend = {backend}")
    if backend == "parsolve":
        # tol is on max|F|; parsolve's own test is on the SUM OF SQUARES, as in the
        # original, so convert. m*tol^2 is the sum that corresponds to a uniform max.
        x, ok, its, worst = parsolve(F, x0, cc=cc, tol=m * tol ** 2, maxcount=maxit,
                                     verbose=verbose, label=label, jac_every=jac_every,
                                     floor_tol=10.0 * tol)
    else:
        x, ok, its, worst = krylov_solve(F, x0, f_tol=tol, maxiter=maxit,
                                         verbose=verbose, label=label)
        if not ok:
            # FALLBACK, deliberately Bocola's method: a dense-Jacobian Newton with the
            # Jacobian reused for a few steps. Newton-Krylov is only a cheaper way of
            # solving the same linear system; when it cannot make progress the dense
            # factorisation still can, at m+1 residual evaluations per refresh.
            if verbose:
                print(f"    krylov stopped at max|F| = {worst:.2e}; "
                      f"falling back to the dense-Jacobian Newton (jac_every=3)")
            x, ok, its2, worst = parsolve(F, x, cc=cc, tol=m * tol ** 2,
                                          maxcount=max(6, maxit // 4), verbose=verbose,
                                          label=label + " dense", jac_every=3)
            its += its2
    write_back(x, rules, regimes)
    if verbose:
        tab = residual_table(x, F)
        worst_eq = max(tab, key=tab.get)
        print(f"    -> {'converged' if ok else 'STOPPED'} in {its} iterations, "
              f"max|F| = {worst:.2e}  (worst equation: {worst_eq} {tab[worst_eq]:.1e})")
    return ok, its, worst
